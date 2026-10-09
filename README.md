# GO Enrichment Skill

A Codex skill for local Gene Ontology enrichment analysis. Turn a gene list into BP, MF and CC dot plots, statistical tables and a concise interpretation.

## Install

Requires Python 3.10+.

```sh
git clone https://github.com/Sculptor815/go-enrichment-skill.git ~/.codex/skills/go-enrichment-skill
python -m pip install -r ~/.codex/skills/go-enrichment-skill/requirements.txt
```

For a custom `CODEX_HOME`, use its `skills` directory. Start a new conversation after installation.

## Use

Attach a gene-list file or paste the gene names, then ask:

> Use $go-enrichment-skill to analyze this gene list, generate BP, MF and CC dot plots, and briefly explain the results.

Human is the default; specify other organisms in your prompt. The skill prepares annotations and runs the analysis locally. Initial annotation downloads may exceed 1 GB.

> **Background matters.** Prefer all genes eligible for selection in your experiment, such as all genes tested for differential expression. If none is supplied, the skill uses GO-annotated genes as the background; this choice affects enrichment results.

## Example results

Actual outputs from an immune-biased synthetic gene list, shown for demonstration.

### BP · Biological process

![GO BP enrichment dot plot](examples/results/plots/immune_biased_BP.png)

### MF · Molecular function

![GO MF enrichment dot plot](examples/results/plots/immune_biased_MF.png)

### CC · Cellular component

![GO CC enrichment dot plot](examples/results/plots/immune_biased_CC.png)

Also available as a Python CLI: `python scripts/go_enrichment.py --help`.

[MIT license](LICENSE) · [Data sources and licenses](NOTICE.md)
