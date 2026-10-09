#!/usr/bin/env python3
"""Local Gene Ontology over-representation analysis and reproducible demos."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import platform
import re
import sys
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import hypergeom

from ontology import read_go
from plotting import dotplot

VERSION = '0.1.0'
ASPECTS = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}
ROOTS = {'GO:0008150', 'GO:0003674', 'GO:0005575'}
SPECIES = {9606: 'Homo_sapiens', 10090: 'Mus_musculus', 10116: 'Rattus_norvegicus'}
RESULT_COLUMNS = ['list', 'source', 'GO_ID', 'description', 'p_value', 'p_adjust_BH',
    'significant', 'Count', 'query_size', 'background_term_size', 'effective_background_size',
    'GeneRatio', 'BgRatio', 'FoldEnrichment', 'overlap_entrez_ids', 'overlap_input_ids']


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8-sig') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8-sig')


def normalize_id(value):
    # Strip only Ensembl gene-version suffixes; retain case and symbol punctuation.
    return re.sub(r'^(ENS[A-Z]*G\d+)\.\d+$', r'\1', str(value).strip())


def read_list(path, column=None):
    """TXT is headerless, one ID per line; CSV/TSV require a named column."""
    if column:
        with open_text(path) as f:
            delimiter = '\t' if '.tsv' in Path(path).suffixes else ','
            reader = csv.DictReader(f, delimiter=delimiter)
            if column not in (reader.fieldnames or []):
                raise ValueError(f'{path}: missing column {column!r}')
            values = [row[column].strip() for row in reader if row.get(column, '').strip()]
    else:
        with open_text(path) as f:
            values = [line.strip() for line in f if line.strip() and not line.lstrip().startswith('#')]
        if any('\t' in v or ',' in v for v in values):
            raise ValueError(f'{path}: use --column for CSV/TSV, or one gene per line for TXT')
    return list(dict.fromkeys(values))


class GeneMapper:
    def __init__(self, path, taxid):
        self.exact, self.aliases, self.ensembl = defaultdict(set), defaultdict(set), defaultdict(set)
        self.symbols = {}
        with open_text(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if fields[0] != str(taxid):
                    continue
                gid, symbol = int(fields[1]), fields[2]
                self.symbols[gid] = symbol
                self.exact[symbol].add(gid)
                if fields[10] != '-':
                    self.exact[fields[10]].add(gid)
                for alias in fields[4].split('|'):
                    if alias != '-':
                        self.aliases[alias].add(gid)
                for xref in fields[5].split('|'):
                    if xref.startswith('Ensembl:'):
                        self.ensembl[normalize_id(xref.split(':', 1)[1])].add(gid)
        if not self.symbols:
            raise ValueError(f'No gene_info records for taxid={taxid}')

    def map(self, genes):
        rows, mapping = [], {}
        for gene in genes:
            key = normalize_id(gene)
            if key.isdigit() and int(key) in self.symbols:
                found, mode = {int(key)}, 'Entrez'
            elif re.fullmatch(r'ENS[A-Z]*G\d+', key):
                found, mode = self.ensembl.get(key, set()), 'Ensembl'
            elif key in self.exact:
                found, mode = self.exact[key], 'symbol'
            else:
                found, mode = self.aliases.get(key, set()), 'alias'
            status = 'mapped' if len(found) == 1 else ('ambiguous' if found else 'unmapped')
            gid = next(iter(found)) if status == 'mapped' else None
            if gid is not None:
                mapping[gene] = gid
            rows.append(dict(gene=gene, lookup_id=key, status=status, mapping_type=mode,
                entrez_id=gid if gid is not None else '', symbol=self.symbols.get(gid, ''),
                candidate_entrez_ids=';'.join(map(str, sorted(found)))))
        return pd.DataFrame(rows, columns=['gene', 'lookup_id', 'status', 'mapping_type',
            'entrez_id', 'symbol', 'candidate_entrez_ids']), mapping


def load_annotations(obo, gene2go, taxid, background, excluded=('ND',)):
    terms, alt, ancestors, version = read_go(Path(obo))
    members, unknown = defaultdict(set), set()
    direct = defaultdict(set)
    prefix = str(taxid) + '\t'
    with open_text(gene2go) as f:
        for line in f:
            if not line.startswith(prefix):
                continue
            fields = line.rstrip('\n').split('\t')
            gid = int(fields[1])
            if gid not in background or 'NOT' in fields[4].split('|') or fields[3] in excluded:
                continue
            go = alt.get(fields[2], fields[2])
            if go not in terms:
                unknown.add(go)
                continue
            direct[go].add(gid)
    for go, ids in direct.items():
        for parent in ancestors(go):
            members[parent].update(ids)
    return terms, members, version, sorted(unknown)


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def enrich(name, aspect, query, background, candidates, terms, reverse, fdr):
    if not query <= background:
        raise ValueError('Mapped query must be a subset of the background')
    n, total = len(query), len(background)
    rows = []
    for go, ids in sorted(candidates.items()):
        overlap = query & ids
        k, m = len(overlap), len(ids)
        ratio, bg_ratio = (k / n if n else 0.0), m / total
        rows.append([name, aspect, go, terms[go]['name'],
            float(hypergeom.sf(k - 1, total, m, n)) if k else 1.0, 1.0, False,
            k, n, m, total, ratio, bg_ratio, ratio / bg_ratio,
            ';'.join(map(str, sorted(overlap))),
            ';'.join(sorted({g for gid in overlap for g in reverse[gid]}))])
    df = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    df['p_adjust_BH'] = bh_adjust(df.p_value)
    df['significant'] = (df.p_adjust_BH < fdr) & (df.FoldEnrichment > 1)
    return df.sort_values(['p_adjust_BH', 'GO_ID'])


def read_loading(path, threshold, selected=None):
    with open_text(path) as f:
        header = next(csv.reader(f))
    if len(header) != len(set(x.strip() for x in header)):
        raise ValueError('Duplicate matrix column names')
    frame = pd.read_csv(path, index_col=0, keep_default_na=False, dtype=str)
    frame.index = frame.index.map(str.strip)
    frame.columns = frame.columns.map(str.strip)
    if not frame.index.is_unique or '' in frame.index or '' in frame.columns:
        raise ValueError('Empty or duplicate matrix identifiers')
    h = frame.to_numpy(dtype=float)
    if h.shape[0] == 0 or h.shape[1] < 2 or not np.isfinite(h).all() or (h < 0).any():
        raise ValueError('Loading matrix needs genes, >=2 programs, and finite nonnegative values')
    sd = h.std(axis=1, ddof=1, keepdims=True)
    z = np.divide(h - h.mean(axis=1, keepdims=True), sd, out=np.zeros_like(h), where=sd > 0)
    names = selected or list(frame.columns)
    if not set(names) <= set(frame.columns):
        raise ValueError('Selected program is missing from loading matrix')
    queries, selections = {}, {}
    for name in names:
        j = frame.columns.get_loc(name)
        mask = z[:, j] > threshold
        queries[name] = frame.index[mask].tolist()
        selections[name] = pd.DataFrame({'gene': queries[name], 'loading': h[mask, j], 'zscore': z[mask, j]})
    return frame.index.tolist(), queries, selections


def safe_name(name):
    value = re.sub(r'[^A-Za-z0-9_.-]', '_', name).strip(' .')
    if not value or value in {'.', '..'} or re.fullmatch(r'(?i)(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\..*)?', value):
        value = 'list_' + hashlib.sha256(name.encode()).hexdigest()[:8]
    return value


def analyze(args):
    selections = {}
    if args.loading:
        if args.background:
            raise ValueError('--loading uses all matrix genes as background; omit --background')
        background_genes, queries, selections = read_loading(args.loading, args.z_threshold, args.program)
        inputs = [args.loading]
    else:
        if not args.background:
            raise ValueError('Supply --background: all genes eligible for selection in your experiment')
        background_genes = read_list(args.background, args.background_column)
        queries, inputs = {}, [args.background]
        for spec in args.query:
            if '=' not in spec:
                raise ValueError('--query requires NAME=PATH')
            name, path = spec.split('=', 1)
            if not name or name in queries:
                raise ValueError('Query names must be nonempty and unique')
            queries[name] = read_list(path, args.column)
            inputs.append(Path(path))
    names = [safe_name(n) for n in queries]
    if len(set(n.lower() for n in names)) != len(names):
        raise ValueError('Query names collide after filename normalization')
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output directory must be empty: {out}')
    mapper = GeneMapper(args.gene_info, args.taxid)
    bg_audit, bg_map = mapper.map(background_genes)
    background = set(bg_map.values())
    if not background:
        raise ValueError('No background genes mapped')
    terms, members, version, unknown = load_annotations(args.obo, args.gene2go, args.taxid,
        background, args.exclude_evidence)
    candidates = {aspect: {go: ids for go, ids in members.items()
        if ASPECTS[terms[go]['namespace']] == aspect and go not in ROOTS
        and args.min_term_size <= len(ids) <= args.max_term_size} for aspect in args.aspects}
    if not any(candidates.values()):
        raise ValueError('No testable GO terms; check species, annotation coverage and term-size limits')
    out.mkdir(parents=True, exist_ok=True)
    for sub in ['tables', 'plots', 'gene_lists', 'mapping']:
        (out / sub).mkdir(exist_ok=True)
    bg_audit.to_csv(out / 'mapping/background.csv', index=False)
    pd.DataFrame({'GO_ID': unknown}).to_csv(out / 'mapping/unknown_go_ids.csv', index=False)
    coverage = []
    for aspect in args.aspects:
        annotated = set().union(*(ids for go, ids in members.items() if ASPECTS[terms[go]['namespace']] == aspect))
        coverage.append(dict(source=aspect, background_size=len(background), annotated_genes=len(annotated),
            tested_terms=len(candidates[aspect])))
    pd.DataFrame(coverage).to_csv(out / 'annotation_coverage.csv', index=False)
    metadata = dict(tool_version=VERSION, created_utc=datetime.now(timezone.utc).isoformat(),
        taxid=args.taxid, ontology_version=version, ontology_doi=args.ontology_doi,
        parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        inputs=[dict(file=Path(p).name, sha256=sha256(p)) for p in inputs],
        annotations={key: dict(file=Path(p).name, sha256=sha256(p))
            for key, p in [('obo', args.obo), ('gene_info', args.gene_info), ('gene2go', args.gene2go)]},
        background_policy='unique mapped genes, including genes without GO annotations',
        multiple_testing='BH separately per query and GO aspect, including zero-hit terms',
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        status='running')
    save_json(out / 'settings.json', metadata)
    summary, failures = [], []
    for name, genes in queries.items():
        stem = safe_name(name)
        audit, mapping = mapper.map(genes)
        outside = set(mapping.values()) - background
        audit['in_background'] = [mapping.get(g) in background for g in audit.gene]
        audit.to_csv(out / 'mapping' / f'{stem}.csv', index=False)
        if outside and args.outside_background == 'error':
            metadata['status'] = 'failed: query outside background'
            save_json(out / 'settings.json', metadata)
            raise ValueError(f'{name}: {len(outside)} mapped genes outside background; inspect mapping or explicitly use --outside-background drop')
        query = set(mapping.values()) & background
        reverse = defaultdict(list)
        for gene, gid in mapping.items():
            if gid in query:
                reverse[gid].append(gene)
        selected = selections.get(name, pd.DataFrame({'gene': genes})).copy()
        selected['entrez_id'] = [mapping.get(g, '') for g in selected.gene]
        selected['included'] = [mapping.get(g) in query for g in selected.gene]
        selected.to_csv(out / 'gene_lists' / f'{stem}.csv', index=False)
        for aspect in args.aspects:
            result = enrich(name, aspect, query, background, candidates[aspect], terms, reverse, args.fdr)
            result.to_csv(out / 'tables' / f'{stem}_{aspect}_all.csv', index=False)
            result.loc[result.significant].to_csv(out / 'tables' / f'{stem}_{aspect}_significant.csv', index=False)
            if not args.no_plots:
                try:
                    dotplot(result, name, aspect, out / 'plots' / f'{stem}_{aspect}', len(query),
                        args.fdr, args.top_terms, args.dpi, input_n=len(genes))
                except Exception as exc:
                    failures.append(dict(list=name, source=aspect, error=str(exc)))
            status = 'no_testable_terms' if result.empty else ('empty_query' if not genes else ('no_mapped_query' if not query else 'ok'))
            summary.append(dict(list=name, source=aspect, input_genes=len(genes), mapped_query=len(query),
                unmapped_or_ambiguous=int((audit.status != 'mapped').sum()), outside_background=len(outside),
                background_size=len(background), tested_terms=len(result), significant_terms=int(result.significant.sum()), status=status))
            print(f'{name} / {aspect}: {len(query)} genes, {int(result.significant.sum())} enriched terms', flush=True)
    pd.DataFrame(summary).to_csv(out / 'summary.csv', index=False)
    pd.DataFrame(failures, columns=['list', 'source', 'error']).to_csv(out / 'plot_failures.csv', index=False)
    metadata['status'] = 'completed_with_plot_errors' if failures else 'completed'
    save_json(out / 'settings.json', metadata)
    (out / 'README.md').write_text(
        '# GO analysis results\n\nGO ontology: ' + version + '\n\n'
        'See summary.csv, settings.json, mapping/, tables/ and plots/.\n'
        'Hypergeometric over-representation; BH correction per query and GO aspect.\n'
        'Background includes mapped genes without annotation. GeneRatio = overlap / mapped query.\n'
        'No enriched terms is a valid result. Enrichment does not establish pathway activation.\n'
        'GO data: https://geneontology.org/ ; CC BY 4.0: https://creativecommons.org/licenses/by/4.0/\n', encoding='utf-8')
    if failures:
        raise RuntimeError('Statistics saved, but some plots failed; see plot_failures.csv')


def download(args):
    directory = Path(args.out)
    directory.mkdir(parents=True, exist_ok=True)
    species = SPECIES.get(args.taxid)
    gene_url = args.gene_info_url or (f'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/{species}.gene_info.gz' if species else None)
    if gene_url is None:
        raise ValueError('For this taxid, supply --gene-info-url or use manually downloaded annotations')
    urls = {'go-basic.obo': args.obo_url, 'gene_info.gz': gene_url,
        'gene2go.gz': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'}
    manifest_path = directory / 'downloads.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    for name, url in urls.items():
        path = directory / name
        previous = manifest.get(name)
        if path.exists():
            if previous and previous['url'] == url and previous['sha256'] == sha256(path):
                print(f'Reusing {name}', flush=True)
                continue
            raise ValueError(f'{path} exists without matching provenance; use a new cache directory')
        print(f'Downloading {url}', flush=True)
        temporary = path.with_name(name + '.part')
        partial_meta = path.with_name(name + '.part.json')
        old_partial = json.loads(partial_meta.read_text()) if partial_meta.exists() else {}
        if temporary.exists() and old_partial.get('url') != url:
            raise ValueError(f'{temporary} has unknown provenance; use a new cache directory')
        save_json(partial_meta, dict(url=url, last_modified=old_partial.get('last_modified')))
        for attempt in range(3):
            offset = temporary.stat().st_size if temporary.exists() else 0
            request_headers = {'User-Agent': 'Mozilla/5.0'}
            if offset:
                request_headers['Range'] = f'bytes={offset}-'
                if old_partial.get('last_modified'):
                    request_headers['If-Range'] = old_partial['last_modified']
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers=request_headers), timeout=120) as response:
                    if response.status == 206:
                        span = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                        if not span or int(span[1]) != offset:
                            raise IOError('Unexpected download byte range')
                        expected_total = int(span[3])
                    else:
                        offset = 0  # Server ignored Range or the source changed: restart this partial.
                        expected_total = int(response.headers.get('Content-Length', 0))
                    headers = {'last_modified': response.headers.get('Last-Modified'), 'resolved_url': response.url}
                    old_partial = dict(url=url, **headers)
                    save_json(partial_meta, old_partial)
                    n, last_progress = offset, time.monotonic()
                    with temporary.open('ab' if offset else 'wb') as f:
                        while chunk := response.read(1024 * 1024):
                            f.write(chunk)
                            n += len(chunk)
                            if time.monotonic() - last_progress > 15:
                                print(f'  {name}: {n / 1024**2:.0f} MiB', flush=True)
                                last_progress = time.monotonic()
                    if n == 0 or (expected_total and n != expected_total):
                        raise IOError('Incomplete download')
                temporary.replace(path)
                partial_meta.unlink()
                break
            except (OSError, http.client.IncompleteRead) as exc:
                if attempt == 2:
                    raise OSError(f'{name}: download failed; partial retained for retry: {exc}') from exc
                print(f'Retrying {name}: {exc}', flush=True)
                time.sleep(attempt + 1)
        manifest[name] = dict(url=url, sha256=sha256(path), downloaded_utc=datetime.now(timezone.utc).isoformat(), **headers)
        save_json(manifest_path, manifest)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=VERSION)
    sub = p.add_subparsers(dest='command', required=True)
    d = sub.add_parser('download', help='Download public GO/NCBI annotations; gene2go may exceed 1 GB')
    d.add_argument('--taxid', type=int, default=9606)
    d.add_argument('--out', type=Path, required=True)
    d.add_argument('--obo-url', default='https://current.geneontology.org/ontology/go-basic.obo')
    d.add_argument('--gene-info-url')
    a = sub.add_parser('analyze', help='Analyze gene lists or a gene x program loading matrix')
    mode = a.add_mutually_exclusive_group(required=True)
    mode.add_argument('--query', action='append', metavar='NAME=PATH')
    mode.add_argument('--loading', type=Path)
    a.add_argument('--background', type=Path)
    a.add_argument('--column', help='Gene ID column in query CSV/TSV; otherwise headerless TXT')
    a.add_argument('--background-column')
    a.add_argument('--program', action='append')
    a.add_argument('--z-threshold', type=float, default=3.0)
    a.add_argument('--taxid', type=int, default=9606)
    for name in ['obo', 'gene-info', 'gene2go']:
        a.add_argument('--' + name, type=Path, required=True)
    a.add_argument('--ontology-doi', help='Optional DOI of the exact ontology release')
    a.add_argument('--aspects', nargs='+', choices=['BP', 'MF', 'CC'], default=['BP', 'MF', 'CC'])
    a.add_argument('--exclude-evidence', nargs='*', default=['ND'], help='NOT is always excluded; default excludes ND and retains IEA')
    a.add_argument('--min-term-size', type=int, default=10)
    a.add_argument('--max-term-size', type=int, default=2000)
    a.add_argument('--fdr', type=float, default=0.05)
    a.add_argument('--top-terms', type=int, default=15)
    a.add_argument('--dpi', type=int, default=180)
    a.add_argument('--outside-background', choices=['error', 'drop'], default='error')
    a.add_argument('--out', type=Path, required=True)
    a.add_argument('--no-plots', action='store_true')
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    try:
        if args.taxid <= 0:
            raise ValueError('taxid must be positive')
        if args.command == 'download':
            download(args)
        else:
            if not (0 < args.fdr < 1 and 1 <= args.min_term_size <= args.max_term_size and args.top_terms >= 1 and args.dpi >= 72 and np.isfinite(args.z_threshold)):
                raise ValueError('Invalid FDR, term-size, plot or Z-threshold settings')
            if len(args.aspects) != len(set(args.aspects)):
                raise ValueError('Duplicate GO aspects')
            analyze(args)
    except (ValueError, OSError, RuntimeError) as exc:
        p.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
