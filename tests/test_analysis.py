"""Statistical invariants and CLI behavior, using small independent fixtures."""
import csv
import gzip
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import uuid
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from go_enrichment import GeneMapper, bh_adjust, download, enrich, load_annotations, read_list, read_loading, safe_name
from ontology import read_go


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(os.environ.get('GO_TEST_TMP', tempfile.gettempdir())) / ('go-test-' + uuid.uuid4().hex)
        self.tmp.mkdir(parents=True)
        self.obo = self.tmp / 'go.obo'
        self.obo.write_text('''format-version: 1.2
data-version: test-fixture

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0000001
name: parent
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000002
name: child
namespace: biological_process
alt_id: GO:9000002
relationship: part_of GO:0000001 ! parent

[Term]
id: GO:0000003
name: other
namespace: biological_process
relationship: regulates GO:0000001 ! parent

[Term]
id: GO:0000004
name: obsolete
namespace: biological_process
is_obsolete: true
''', encoding='utf-8')
        self.info = self.tmp / 'gene_info.gz'
        rows = []
        for i in range(1, 11):
            aliases = 'SHARED' if i in (1, 2) else ('ALT3' if i == 3 else '-')
            rows.append(f'9606\t{i}\tG{i}\t-\t{aliases}\tEnsembl:ENSG{i:011d}\t1\t-\tgene\tprotein-coding\tG{i}\t-\tO\t-\t20260101\t-\n')
        self.info.write_bytes(gzip.compress(''.join(rows).encode()))
        self.g2g = self.tmp / 'gene2go.gz'
        self.g2g.write_bytes(gzip.compress(('''#tax_id\tGeneID\tGO_ID\tEvidence\tQualifier\tGO_term\tPubMed\tCategory
9606\t1\tGO:9000002\tIDA\t-\tchild\t-\tProcess
9606\t2\tGO:0000002\tIDA\tNOT\tchild\t-\tProcess
9606\t3\tGO:0000002\tND\t-\tchild\t-\tProcess
9606\t4\tGO:0000002\tIEA\t-\tchild\t-\tProcess
9606\t5\tGO:0000003\tIDA\t-\tother\t-\tProcess
10090\t6\tGO:0000002\tIDA\t-\tchild\t-\tProcess
9606\t7\tGO:0000004\tIDA\t-\tobsolete\t-\tProcess
''').encode()))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_hypergeometric_matches_combinatorial_tail(self):
        # N=10, M=4, n=3, k=2: (C(4,2)C(6,1)+C(4,3)C(6,0))/C(10,3)=1/3.
        terms = {'GO:1': {'name': 'term'}}
        df = enrich('q', 'BP', {1, 2, 9}, set(range(1, 11)), {'GO:1': {1, 2, 3, 4}},
            terms, {1: ['a'], 2: ['b']}, 0.05)
        self.assertAlmostEqual(df.iloc[0].p_value, 1/3)
        self.assertAlmostEqual(df.iloc[0].GeneRatio, 2/3)
        self.assertAlmostEqual(df.iloc[0].FoldEnrichment, (2/3)/(4/10))

    def test_bh_known_values_in_original_order(self):
        np.testing.assert_allclose(bh_adjust([0.03, 0.001, 1.0, 0.01]), [0.04, 0.004, 1.0, 0.02])

    def test_zero_hit_terms_are_in_bh_family(self):
        terms = {'GO:1': {'name': 'hit'}, 'GO:2': {'name': 'zero'}}
        df = enrich('q', 'BP', {1}, set(range(1, 101)), {'GO:1': {1}, 'GO:2': {2}}, terms, {1:['1']}, .05)
        self.assertAlmostEqual(df.iloc[0].p_adjust_BH, .02)
        self.assertEqual(df.iloc[1].p_value, 1)

    def test_mapping_precedence_alias_ambiguity_and_ensembl_versions(self):
        audit, mapped = GeneMapper(self.info, 9606).map(['G1', '1', 'ENSG00000000001.7', 'ALT3', 'SHARED', 'missing'])
        self.assertEqual(set(mapped.values()), {1, 3})
        self.assertEqual(audit.set_index('gene').loc['SHARED', 'status'], 'ambiguous')
        self.assertEqual(audit.set_index('gene').loc['missing', 'status'], 'unmapped')

    def test_species_mismatch_is_an_error(self):
        with self.assertRaises(ValueError):
            GeneMapper(self.info, 10090)

    def test_propagation_and_evidence_filters(self):
        terms, members, version, unknown = load_annotations(self.obo, self.g2g, 9606, set(range(1, 11)))
        self.assertEqual(members['GO:0000001'], {1, 4})
        self.assertEqual(members['GO:0000002'], {1, 4})
        self.assertNotIn(5, members['GO:0000001'])  # regulates is not propagated
        self.assertEqual(unknown, ['GO:0000004'])

    def test_optional_iea_exclusion(self):
        _, members, _, _ = load_annotations(self.obo, self.g2g, 9606, set(range(1,11)), ['ND','IEA'])
        self.assertEqual(members['GO:0000002'], {1})

    def test_query_outside_background_rejected_by_statistical_core(self):
        with self.assertRaises(ValueError):
            enrich('q', 'BP', {2}, {1}, {}, {}, {}, .05)

    def test_empty_query_returns_unit_pvalues(self):
        df = enrich('q', 'BP', set(), {1,2}, {'GO:1': {1}}, {'GO:1': {'name':'term'}}, {}, .05)
        self.assertEqual(df.iloc[0].p_value, 1)
        self.assertEqual(df.iloc[0].GeneRatio, 0)
        self.assertFalse(df.iloc[0].significant)

    def test_loading_sample_sd_constant_rows_and_strict_threshold(self):
        path = self.tmp/'loading.csv'
        path.write_text('gene,A,B,C\nG1,0,0,9\nG2,4,4,4\nG3,0,0,0\n')
        bg, queries, selection = read_loading(path, 1.0)
        self.assertEqual(queries['C'], ['G1'])
        self.assertAlmostEqual(selection['C'].iloc[0].zscore, 2/math.sqrt(3))
        _, queries, _ = read_loading(path, 2/math.sqrt(3))
        self.assertEqual(queries['C'], [])

    def test_loading_invalid_and_duplicate_identifiers(self):
        path = self.tmp/'bad.csv'
        for text in ['gene,A,A\nG1,1,2\n', 'gene,A,B\nG1,1,2\nG1,1,2\n', 'gene,A,B\nG1,-1,2\n']:
            path.write_text(text)
            with self.assertRaises(ValueError):
                read_loading(path, 3)

    def test_gene_lists_preserve_na_and_leading_zero_strings(self):
        path = self.tmp/'genes.csv'
        path.write_text('gene\nNA\n001\nNA\n')
        self.assertEqual(read_list(path, 'gene'), ['NA', '001'])

    def test_filename_safety(self):
        for name in ['../../x', 'CON', 'NUL.txt', '中文']:
            value = safe_name(name)
            self.assertNotIn('/', value)
            self.assertNotIn('\\', value)
            self.assertNotIn(value.upper(), {'CON', 'NUL.TXT'})

    def cli(self, *extra):
        bg = self.tmp/'bg.txt'; bg.write_text('\n'.join(f'G{i}' for i in range(1, 11)))
        query = self.tmp/'query.txt'; query.write_text('G1\n1\nmissing\n')
        return subprocess.run([sys.executable, str(SCRIPTS/'go_enrichment.py'), 'analyze',
            '--background', str(bg), '--query', f'query={query}', '--obo', str(self.obo),
            '--gene-info', str(self.info), '--gene2go', str(self.g2g), '--min-term-size', '1',
            '--out', str(self.tmp/'results'), *extra], capture_output=True, text=True)

    def test_cli_outputs_deduplication_unannotated_background_and_empty_aspects(self):
        run = self.cli()
        self.assertEqual(run.returncode, 0, run.stderr)
        summary = pd.read_csv(self.tmp/'results/summary.csv')
        self.assertTrue((summary.mapped_query == 1).all())
        self.assertTrue((summary.background_size == 10).all())
        self.assertEqual(summary.set_index('source').loc['MF','status'], 'no_testable_terms')
        self.assertTrue((self.tmp/'results/plots/query_BP.png').is_file())
        self.assertEqual(json.loads((self.tmp/'results/settings.json').read_text())['status'], 'completed')

    def test_cli_never_overwrites_results(self):
        self.assertEqual(self.cli('--no-plots').returncode, 0)
        second = self.cli('--no-plots')
        self.assertNotEqual(second.returncode, 0)
        self.assertIn('must be empty', second.stderr)

    def test_download_resumes_truncated_response_and_reuses_verified_cache(self):
        class Response(io.BytesIO):
            def __init__(self, payload, status=200, headers=None):
                super().__init__(payload)
                self.status = status
                self.headers = headers or {'Content-Length': str(len(payload))}
                self.url = 'https://example.org/annotation'
        args = SimpleNamespace(out=self.tmp/'cache', taxid=9606, gene_info_url=None,
            obo_url='https://example.org/go-basic.obo')
        replies = [Response(b'ab', headers={'Content-Length':'4', 'Last-Modified':'Wed, 01 Jan 2025 00:00:00 GMT'}),
            Response(b'cd', 206, {'Content-Range':'bytes 2-3/4', 'Content-Length':'2'}),
            Response(b'gene_info'), Response(b'gene2go')]
        with patch('go_enrichment.urllib.request.urlopen', side_effect=replies) as request, patch('go_enrichment.time.sleep'):
            download(args)
            self.assertEqual((args.out/'go-basic.obo').read_bytes(), b'abcd')
            self.assertEqual(request.call_args_list[1].args[0].get_header('Range'), 'bytes=2-')
            download(args)
            self.assertEqual(request.call_count, 4)
            (args.out/'go-basic.obo').write_bytes(b'changed')
            with self.assertRaises(ValueError):
                download(args)


if __name__ == '__main__':
    unittest.main()
