# Analysis semantics

## Universe and identifiers

The experimental background should include every gene that could have been selected into a query. For expression-based selection, this usually means the genes eligible after the experiment's detection and testing filters. The software requires this background in gene-list mode; it does not silently pick all genes in the organism. See the [GO Consortium enrichment guide](https://geneontology.org/docs/go-enrichment-analysis/).

NCBI gene_info supplies current Entrez GeneIDs, official symbols, unique aliases and Ensembl cross-references for the chosen taxonomy ID. Exact official symbol matches take priority over aliases. Ambiguous aliases and unmapped identifiers are audited and excluded. Ensembl gene-version suffixes are removed; symbol case and other suffixes are preserved. Only unique mapped Entrez IDs are counted. Numeric inputs are interpreted as Entrez IDs; symbols consisting entirely of numbers are not disambiguated automatically. Historical IDs absent from gene_info are not rescued using gene_history.

The effective background contains all unique mapped background genes, including those without usable GO annotations. Query size is the number of unique mapped query genes in that background. Mapping different aliases to one ID does not increase counts. Out-of-background mapped query genes cause an error by default; explicit `--outside-background drop` records and excludes them. Input file duplicates are removed before reporting `input_genes`.

## Ontology and annotations

Input files are [go-basic.obo](https://geneontology.org/docs/download-ontology/), [NCBI gene_info](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/), and [NCBI gene2go](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz). Full or species-filtered gene2go files are supported. The taxonomy ID must match the data. Availability of gene_info does not guarantee useful gene2go coverage for every organism; inspect annotation_coverage.csv.

Alternate GO IDs resolve to active terms. Obsolete or unknown terms are excluded and audited. Annotations propagate through `is_a` and `part_of`, including the annotated term itself. `regulates`, `has_part`, and other relationships are not propagated. Three GO roots are excluded. `NOT` annotations are always excluded; `ND` is excluded by default, while `IEA` is retained. `--exclude-evidence ND IEA` additionally excludes electronic annotations. Full ontology files with cycles are not supported; use go-basic.

Candidates are chosen by GO aspect (BP/MF/CC) and annotated gene count in the background, before examining the query. Defaults: 10–2000 genes per term. Every query uses the same candidate universe for the same background and aspect.

## Statistics

For background size N, term members M, query size n, and overlap k, the one-sided over-representation p-value is P[X >= k] for X ~ Hypergeometric(N, M, n). `scipy.stats.hypergeom.sf(k - 1, N, M, n)` evaluates this tail. Zero overlaps and empty queries have p=1.

Benjamini–Hochberg correction is applied separately to **all candidate terms within each query × GO aspect**, including terms with zero overlaps. It does not control a single combined family across every query/aspect. GO terms share genes and are dependent; related significant terms do not constitute independent discoveries. Significant means adjusted p < the selected threshold (default 0.05) and fold enrichment > 1. No ranked GSEA, under-representation test, semantic redundancy reduction, or batch-wide FDR is implemented.

- GeneRatio = k/n (0 for an empty query).
- BgRatio = M/N.
- FoldEnrichment = (k/n)/(M/N).

Mapped but unannotated genes remain in N and n. This preserves the original workflow's convention and can differ from tools that restrict the universe to annotated genes. Database versions, evidence filters and universe conventions can also explain cross-tool differences.

## Loading matrices

CSV rows are genes, columns are programs, and the first column is the gene ID. Values must be finite and nonnegative; at least two program columns are required. For each gene, Z is `(loading - mean across all programs) / sample SD across all programs`, with `ddof=1`. Constant rows get Z=0. Selection is strictly greater than the chosen threshold. All matrix genes form the background; neither highly variable genes nor protein-coding genes are imposed as an additional filter. When aliases map to one GeneID, selection of any alias includes the ID once.

For K program columns, the greatest possible sample-standardized Z is (K-1)/sqrt(K). Thus the default Z > 3 cannot select anything when K <= 10. This is a mathematical consequence of the selected scoring rule, not a software error.

## Outputs and reproducibility

`settings.json` stores parameters, input/annotation SHA-256 hashes, ontology version and software versions. Supply `--ontology-doi` if you know the exact release DOI. Keep the exact annotation files: re-downloading a changing URL does not reproduce old results. NCBI gene2go is an independently updated annotation source and need not be synchronized with the GO ontology release.

`summary.csv` reports one row per query/aspect; `annotation_coverage.csv` distinguishes background size from annotated size. `mapping/` contains gene mapping and unknown GO audits; `gene_lists/` contains the submitted/selected query IDs. `tables/*_all.csv` contains every tested term, while `*_significant.csv` is the significant subset. Each dot plot has a corresponding `*_plotted_terms.csv`, PNG and vector PDF. `plot_failures.csv` records failures, which also cause a nonzero exit.

Plots select up to `--top-terms` significant terms by adjusted p, then sort them by GeneRatio. Dot area encodes overlap count; color is -log10(adjusted p). Color and size scales are local to each plot; do not infer cross-panel differences from color/area alone. Empty plots explicitly state that no enrichment passed the threshold.

## Synthetic examples

The published example generator samples 6,000 human protein-coding genes uniformly as an artificial experimental universe. Two 120-gene controls are uniform samples without replacement. Two further lists contain 90 genes sampled from a selected GO category plus 30 genes sampled outside it. Targets are immune system process (GO:0002376) and DNA repair (GO:0006281). Seed: 815. There is no significance-based seed selection.

The GO-biased examples are intentionally circular positive demonstrations of the pipeline. They are not independent validation of biological hypotheses. Random controls may occasionally contain significant terms; the generator never retries to force significance or nonsignificance. The bundled annotation subset is only sufficient for this demo universe. Re-run `scripts/run_demo.py` with the snapshot to reproduce the numerical results; generating with fresh full annotations can change outputs even with the same seed.
