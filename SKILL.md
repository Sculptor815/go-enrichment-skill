---
name: go-enrichment-skill
description: Run local Gene Ontology over-representation analysis from a user-provided gene list. Produce mapping audits, corrected statistics and GO dot plots, with an automatic annotation-based background or an optional custom background. Use for GO enrichment requests; this is not ranked GSEA or differential expression analysis.
---

# GO enrichment from a gene list

Use the Python tools in this skill directory. Resolve script and reference paths relative to this file; resolve user inputs and outputs relative to the user's working directory. Read [analysis semantics](references/analysis.md) when interpreting background choices, evidence filters or results.

## Accept the list

- The user only needs to supply the genes, as a pasted list or a file. Save pasted IDs to a headerless TXT file, one ID per line. Accept TXT or a CSV/TSV with a gene column. Symbols are case-sensitive; Entrez and Ensembl gene IDs are supported.
- Default to human (`--taxid 9606`) unless the user specifies another organism or the context identifies one. State the organism used. If the data conflict with that assumption, resolve the mismatch rather than interpreting failed mapping as a biological result.
- A background file is optional. By default, use the organism's genes with usable GO annotations in the selected aspects, before term-size filtering. Do not require the user to provide a background. If the user supplies an experiment-specific universe, use `--background` and report that policy.
- Check Python >=3.10 and `requirements.txt` dependencies. Use the selected Python environment or a local virtual environment. Reuse full, species-compatible local annotations when available; the bundled demo snapshot is restricted to its own synthetic background.

## Background reminder

Before analysis, briefly remind the user in their language that background choice affects enrichment results. For experimental data, recommend all genes eligible for selection (for example, all genes tested in the same differential-expression comparison), not only significant genes. If no background is supplied, state that the default GO-annotated universe will be used and continue without requiring extra confirmation. Respect an already specified background and mention any relevant limitation when interpreting results.

## Run

If full annotations are missing, manage their download as part of the analysis:

```sh
python scripts/go_enrichment.py download --taxid 9606 --out annotation-cache
```

Only public annotations are downloaded; gene lists stay local. NCBI gene2go can exceed 1 GB. Human, mouse and rat have built-in gene_info URLs; other supported species need `--gene-info-url` or compatible local files.

Then run:

```sh
python scripts/go_enrichment.py analyze --genes genes.txt
```

This uses `annotation-cache/`, analyzes BP/MF/CC, and writes a new timestamped directory under `results/`. Supply `--obo`, `--gene-info` and `--gene2go` for annotations elsewhere. Use `--column gene` for tabular input. An optional `--out` must identify a new or empty directory. For batches, use repeated `--query NAME=PATH` instead of `--genes`.

## Check and explain

1. Read `settings.json`, `summary.csv`, `annotation_coverage.csv` and mapping audits. A nonzero exit, incomplete status or plot error is not a complete successful analysis.
2. Report the organism, background policy, effective query/background sizes, excluded/unmapped/ambiguous counts, annotation version and significant results. State whether the background reflects the experiment's eligible genes or is an annotation-based fallback; retain the fallback limitation in the final interpretation. In default mode, mapped genes without usable GO annotations are excluded and audited. With a custom background, mapped query genes outside it cause an error unless exclusion was intentionally selected with `--outside-background drop`.
3. Explain that BH correction is separate for each list × GO aspect and includes zero-hit candidate terms. GeneRatio is overlap/eligible mapped query size. Enrichment alone establishes neither activation nor repression.
4. Show relevant plots and link full tables. Nonsignificant and empty results are valid; never add nonsignificant terms to a significant-only plot.

For an offline demonstration, run `python scripts/run_demo.py --out demo-results`. It uses seeded synthetic gene lists, including two uniform controls and two deliberately GO-biased lists, with an explicit synthetic background. Describe them as demonstrations, not biological discoveries. Generate fresh examples with `scripts/make_demo.py --help`; retain provenance and do not select seeds based on significance.
