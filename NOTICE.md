# Data attribution

The source code is licensed under the MIT License. Third-party biological data retain their own terms; they are not relicensed as MIT.

GO term names and ontology content, including the filtered example ontology and derived result tables, originate from the [Gene Ontology Consortium](https://geneontology.org/). GO data are licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/) and are provided without warranties. The demo preserves a subset of the source ontology, selected for its synthetic background and the relevant `is_a` / `part_of` ancestors. The exact source version and file hashes are recorded in [examples/provenance.json](examples/provenance.json). See the [GO citation and license policy](https://geneontology.org/docs/go-citation-policy/).

Gene identifiers and annotations come from [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) public [FTP downloads](https://ftp.ncbi.nlm.nih.gov/gene/DATA/). The bundled gene_info and gene2go files are filtered to the example background; gene2go retains its original evidence/qualifier fields. See [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).

When publishing analyses, cite the GO resource and relevant annotation sources, state their versions, and describe the test and background. This project is independent and is not endorsed by the GO Consortium or NCBI.
