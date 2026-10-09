"""Reproduce all published numerical tables from the bundled public-data subset."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import uuid

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from go_enrichment import sha256


class PublishedDemoTests(unittest.TestCase):
    def test_annotation_snapshot_hashes(self):
        demo = ROOT/'examples'
        metadata = json.loads((demo/'provenance.json').read_text(encoding='utf-8'))
        for name, expected in metadata['snapshot_files'].items():
            self.assertEqual(sha256(demo/'annotations'/name), expected, name)

    def test_offline_run_matches_published_tables(self):
        destination = Path(os.environ.get('GO_TEST_TMP', tempfile.gettempdir())) / ('go-demo-' + uuid.uuid4().hex)
        try:
            run = subprocess.run([sys.executable, str(ROOT/'scripts/run_demo.py'), '--out', str(destination)],
                capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            expected = ROOT/'examples/results'
            for path in [expected/'summary.csv', expected/'annotation_coverage.csv', *sorted((expected/'tables').glob('*.csv'))]:
                actual = destination/path.relative_to(expected)
                pd.testing.assert_frame_equal(pd.read_csv(actual), pd.read_csv(path),
                    check_dtype=False, check_exact=False, rtol=1e-10, atol=1e-300)
            self.assertEqual(len(list((destination/'plots').glob('*.png'))), 12)
            self.assertEqual(len(list((destination/'plots').glob('*.pdf'))), 12)
        finally:
            if destination.exists():
                shutil.rmtree(destination)


if __name__ == '__main__':
    unittest.main()
