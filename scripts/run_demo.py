#!/usr/bin/env python3
"""Run the checked-in examples offline. Execute from any working directory."""
import argparse
import json
from pathlib import Path
from go_enrichment import main as analyze, sha256


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', default='demo-results', help='New/empty results directory')
    args = p.parse_args()
    demo = Path(__file__).resolve().parents[1] / 'examples'
    # Prefer portable relative paths in metadata when invoked from the repo root.
    try:
        demo = demo.relative_to(Path.cwd())
    except ValueError:
        pass
    provenance = json.loads((demo/'provenance.json').read_text(encoding='utf-8'))
    for name, expected in provenance['snapshot_files'].items():
        if sha256(demo/'annotations'/name) != expected:
            p.error(f'Example annotation checksum mismatch: {name}')
    argv = ['analyze', '--background', str(demo/'inputs/background.txt'),
        '--obo', str(demo/'annotations/go-basic.obo'), '--gene-info', str(demo/'annotations/gene_info.gz'),
        '--gene2go', str(demo/'annotations/gene2go.gz'), '--out', args.out]
    for name in ['random_control_1', 'random_control_2', 'immune_biased', 'dna_repair_biased']:
        argv += ['--query', f'{name}={demo / "inputs" / (name + ".txt")}']
    analyze(argv)


if __name__ == '__main__':
    main()
