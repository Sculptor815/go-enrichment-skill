---
name: go-enrichment-skill
description: Run local Gene Ontology over-representation analysis for one or more gene lists against an explicit background, or extract lists from gene-by-program loading matrices. Produce mapping audits, corrected statistics, and GO dot plots. Use for GO enrichment requests; this is not ranked GSEA or differential expression analysis.
---

# GO enrichment

Use the Python tools in this skill directory. Resolve script and reference paths relative to this file; resolve user data and output paths relative to the user's working directory. Read [analysis semantics](references/analysis.md) when choosing background, evidence filters, or interpreting results.

## Choose inputs

- For ordinary gene lists, establish the organism, query file(s), and the background of genes eligible for selection in the experiment. If organism or background is unknown, ask for it before claiming interpretable enrichment. Do not substitute the query itself, a protein-coding universe, or the demo background without a justified, explicit choice.
- Accept headerless TXT (one gene per line) or CSV/TSV with a specified gene column. Multiple lists are named separately. Symbols are case-sensitive; Entrez and Ensembl gene IDs are supported. Preserve input identifiers; the tool audits ambiguous and unmapped values.
- For nonnegative gene × program loading matrices, `--loading` uses all input genes as background and selects strictly Z > threshold across all program columns, with sample SD (`ddof=1`). `--program` selects outputs, not the columns used to calculate Z. Z is descriptive, not a normal-test statistic. With few program columns, Z > 3 may be impossible; do not silently lower it.
- Check Python >=3.10 and dependencies in `requirements.txt`. Install into the user's selected environment or a local virtual environment when needed.

## Run

Existing local GO/NCBI annotation files can be supplied directly. To obtain public annotations, run:

```sh
python scripts/go_enrichment.py download --taxid 9606 --out annotation-cache
```

This downloads public files only; the scripts do not submit gene lists. The all-species `gene2go.gz` can exceed 1 GB. Human, mouse and rat have built-in gene_info download URLs; other NCBI-supported organisms need `--gene-info-url` or local species-compatible files. The included demo snapshot is restricted to its synthetic background and must not be reused as a general annotation database.

```sh
python scripts/go_enrichment.py analyze --query treated=genes.txt --background background.txt --taxid 9606 --obo annotation-cache/go-basic.obo --gene-info annotation-cache/gene_info.gz --gene2go annotation-cache/gene2go.gz --out results/treated
```

Repeat `--query NAME=PATH` for batches. Add `--column gene` and/or `--background-column gene` for tabular inputs. Use `--loading loadings.csv --z-threshold 3` instead of `--query` and `--background` for program matrices. Use a new or empty output directory. Consult `--help` for evidence filters, term-size limits and plot settings.

## Check and explain

1. Read `settings.json`, `summary.csv`, and `annotation_coverage.csv`. A nonzero exit, incomplete status, or plot failure is not a successful complete analysis.
2. Inspect mapping audits for unmapped/ambiguous inputs and genes outside background. The default rejects mapped query genes outside background. Use `--outside-background drop` only when the exclusion is intentional and report its count.
3. Report the organism, effective query/background sizes, annotation version, filters, and BH correction scope (each list × GO aspect separately). Zero-hit candidate terms participate in correction. Mapped genes without annotation remain in the denominator.
4. Show relevant dot plots and link complete result tables. An empty or nonsignificant result is valid; never fill a significant-only plot with nonsignificant terms. GeneRatio is overlap/query size, not fold enrichment. Enrichment alone establishes neither activation nor repression.

For an offline demonstration, run `python scripts/run_demo.py --out demo-results`. These are seeded synthetic lists, including two uniform random controls and two deliberately GO-biased lists. Describe them as demonstrations, not biological discoveries. To generate a fresh seeded example from full annotations, use `scripts/make_demo.py --help`; keep the provenance and do not select seeds by significance.
