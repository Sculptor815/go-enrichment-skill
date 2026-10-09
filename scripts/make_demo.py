#!/usr/bin/env python3
"""Create seeded synthetic human gene lists and a compact public-data snapshot."""
from __future__ import annotations
import argparse
import gzip
from pathlib import Path
import random
import sys

from go_enrichment import GeneMapper, load_annotations, open_text, save_json, sha256
from ontology import read_go


def write_gzip(path, text):
    Path(path).write_bytes(gzip.compress(text.encode('utf-8'), mtime=0))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--obo', type=Path, required=True)
    p.add_argument('--gene-info', type=Path, required=True)
    p.add_argument('--gene2go', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--seed', type=int, default=815)
    p.add_argument('--background-size', type=int, default=6000)
    p.add_argument('--query-size', type=int, default=120)
    args = p.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        p.error('--out must be empty')
    rng = random.Random(args.seed)
    records = {}
    with open_text(args.gene_info) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[0] == '9606' and fields[9] == 'protein-coding':
                records[int(fields[1])] = line
    if not 1 <= args.query_size <= args.background_size <= len(records):
        p.error('Require 1 <= query size <= background size <= available protein-coding genes')
    background = set(rng.sample(sorted(records), args.background_size))
    for sub in ['inputs', 'annotations']:
        (args.out/sub).mkdir(parents=True, exist_ok=True)
    print('Filtering NCBI gene2go for the synthetic background...', flush=True)
    kept_rows = []
    with open_text(args.gene2go) as f:
        for line in f:
            if not line.startswith('9606\t'):
                continue
            fields = line.rstrip('\n').split('\t')
            if int(fields[1]) in background:
                kept_rows.append(line)
    write_gzip(args.out/'annotations/gene2go.gz', ''.join(kept_rows))
    terms, members, version, unknown = load_annotations(args.obo, args.out/'annotations/gene2go.gz', 9606, background)
    queries, specification = {}, []
    for name in ['random_control_1', 'random_control_2']:
        queries[name] = sorted(rng.sample(sorted(background), args.query_size))
        specification.append(dict(name=name, design='uniform without replacement', size=args.query_size))
    for name, go in [('immune_biased', 'GO:0002376'), ('dna_repair_biased', 'GO:0006281')]:
        intended = round(args.query_size * 0.75)
        available = members.get(go, set())
        if len(available) < intended or len(background - available) < args.query_size - intended:
            p.error(f'Not enough genes for {go}; increase background or decrease query size')
        genes = rng.sample(sorted(available), intended) + rng.sample(sorted(background - available), args.query_size - intended)
        queries[name] = sorted(genes)
        specification.append(dict(name=name, design='GO-biased synthetic sampling', target_go=go,
            size=args.query_size, target_members=intended, other_genes=args.query_size-intended))
    def save_list(name, values):
        (args.out/'inputs'/f'{name}.txt').write_text('\n'.join(map(str, values))+'\n', encoding='utf-8')
    save_list('background', sorted(background))
    for name, ids in queries.items():
        save_list(name, ids)
    # Preserve all source rows for the chosen background, including excluded evidence.
    write_gzip(args.out/'annotations/gene_info.gz', ''.join(records[g] for g in sorted(background)))
    referenced = set()
    all_terms, alt, ancestors, _ = read_go(args.obo)
    for line in kept_rows:
        fields = line.rstrip('\n').split('\t')
        go = alt.get(fields[2], fields[2])
        if go in all_terms:
            referenced.update(ancestors(go))
    # Keep original term blocks, including relation types which the engine ignores.
    content = args.obo.read_text(encoding='utf-8')
    blocks = content.split('\n[Term]\n')
    selected = [blocks[0]]
    for block in blocks[1:]:
        term_id = block.splitlines()[0].removeprefix('id: ').strip()
        if term_id in referenced:
            selected.append(block.split('\n[Typedef]\n', 1)[0])
    (args.out/'annotations/go-basic.obo').write_text('\n[Term]\n'.join(selected), encoding='utf-8')
    urls = dict(obo='https://current.geneontology.org/ontology/go-basic.obo',
        gene_info='https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
        gene2go='https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz')
    save_json(args.out/'provenance.json', dict(seed=args.seed, taxid=9606, ontology_version=version,
        background_design='uniform sample of human NCBI protein-coding GeneIDs',
        background_size=len(background), queries=specification,
        synthetic=True, notice='Annotation-biased demos illustrate workflow, not independent biological discoveries. Random controls are never resampled based on significance.',
        snapshot_scope='Only the demo background and its GO terms / is_a and part_of ancestors; not suitable for arbitrary user gene lists.',
        source_files={key: dict(url=urls[key], sha256=sha256(path)) for key, path in
            [('obo', args.obo), ('gene_info', args.gene_info), ('gene2go', args.gene2go)]},
        snapshot_files={path.name: sha256(path) for path in sorted((args.out/'annotations').iterdir())},
        unknown_go_ids=unknown))
    print(f'Created {len(queries)} synthetic lists and an offline snapshot in {args.out}', flush=True)


if __name__ == '__main__':
    main()
