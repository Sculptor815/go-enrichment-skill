# Analysis semantics

## Universe and identifiers

The experimental background should include every gene that could have been selected into a query. For expression-based selection, this usually means the genes eligible after the experiment's detection and testing filters. A custom experimental background is optional. Without one, the software uses every gene from the selected organism that has at least one usable GO annotation in any selected aspect. This union is computed before term-size filtering, after the evidence, qualifier and active-term rules below; it is shared across selected aspects. The default organism is human (9606). A genome-wide annotated background can differ from an experiment-specific tested universe, so the software records the policy and size explicitly. See the [GO Consortium enrichment guide](https://geneontology.org/docs/go-enrichment-analysis/).

NCBI gene_info supplies current Entrez GeneIDs, official symbols, unique aliases and Ensembl cross-references for the chosen taxonomy ID. Exact official symbol matches take priority over aliases. Ambiguous aliases and unmapped identifiers are audited and excluded. Ensembl gene-version suffixes are removed; symbol case and other suffixes are preserved. Only unique mapped Entrez IDs are counted. Numeric inputs are interpreted as Entrez IDs; symbols consisting entirely of numbers are not disambiguated automatically. Historical IDs absent from gene_info are not rescued using gene_history.

When a custom background is provided, the effective background contains all unique mapped background genes, including those without usable GO annotations. With the default background, mapped query genes without usable annotations in the selected aspects are excluded and recorded as `no_usable_go_annotation`. Query size is the number of unique mapped query genes in that background. Mapping different aliases to one ID does not increase counts. For a custom background, out-of-background mapped query genes cause an error by default; explicit `--outside-background drop` records and excludes them. Input file duplicates are removed before reporting `input_genes`.

### Choosing a background

Choose the universe from the experiment's selection process before inspecting enrichment significance. Background composition changes term frequencies, fold enrichment, the candidate-term set and adjusted P-values; a poorly matched universe can create apparent enrichment or hide genuine associations.

| Query origin | Recommended experimental background |
|---|---|
| Differential-expression results | All genes eligible after the relevant detection/filtering steps and tested in that same comparison, including nonsignificant genes |
| Targeted assay or genetic screen | All assayed genes that passed the eligibility/QC criteria and could have been selected |
| Exploratory list without an available experimental universe | The default GO-annotated universe, with its limitation explicitly reported |

The background should contain the query plus eligible nonselected genes, with compatible species and identifiers. Neither the query alone nor only significant genes represents the selection universe. Do not choose or adjust a universe to maximize significance. For multiple comparisons with different eligible gene sets, use separate runs with the appropriate background for each; a shared universe is appropriate only when the selection opportunities are shared.

State the chosen policy before running and in the final interpretation. When the experiment-specific universe is unavailable, proceeding with the default is allowed, but the result is relative to the annotation-based reference rather than a verified experimental sampling universe. Acknowledge this limitation instead of requiring a background file or an extra approval step. This follows the [GO Consortium's recommendation to use the genes from which the query was selected](https://geneontology.org/docs/go-enrichment-analysis/).

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

With a custom background, mapped but unannotated genes remain in N and n. The default background restricts N and n to genes with usable GO annotation in the selected aspects. Database versions, evidence filters and universe conventions can also explain cross-tool differences.

## Outputs and reproducibility

`settings.json` stores parameters, input/annotation SHA-256 hashes, ontology version and software versions. Supply `--ontology-doi` if you know the exact release DOI. Keep the exact annotation files: re-downloading a changing URL does not reproduce old results. NCBI gene2go is an independently updated annotation source and need not be synchronized with the GO ontology release.

`summary.csv` reports one row per query/aspect; `annotation_coverage.csv` distinguishes background size from annotated size. `mapping/` contains gene mapping and unknown GO audits; `gene_lists/` contains the submitted query IDs. `tables/*_all.csv` contains every tested term, while `*_significant.csv` is the significant subset. Each dot plot has a corresponding `*_plotted_terms.csv`, PNG and vector PDF. `plot_failures.csv` records failures, which also cause a nonzero exit.

Plots select up to `--top-terms` significant terms by adjusted p, then sort them by GeneRatio. Dot area encodes overlap count; color is -log10(adjusted p). Color and size scales are local to each plot; do not infer cross-panel differences from color/area alone. Empty plots explicitly state that no enrichment passed the threshold.

## Synthetic examples

The published example generator samples 6,000 human protein-coding genes uniformly as an artificial experimental universe. Two 120-gene controls are uniform samples without replacement. Two further lists contain 90 genes sampled from a selected GO category plus 30 genes sampled outside it. Targets are immune system process (GO:0002376) and DNA repair (GO:0006281). Seed: 815. There is no significance-based seed selection.

The GO-biased examples are intentionally circular positive demonstrations of the pipeline. They are not independent validation of biological hypotheses. Random controls may occasionally contain significant terms; the generator never retries to force significance or nonsignificance. The bundled annotation subset is only sufficient for this demo universe. Re-run `scripts/run_demo.py` with the snapshot to reproduce the numerical results; generating with fresh full annotations can change outputs even with the same seed.
