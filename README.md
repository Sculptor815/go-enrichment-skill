# GO Enrichment Skill

A Codex skill for GO enrichment analysis. Supply a gene list to generate BP, MF and CC dot plots, full statistical tables and a brief interpretation.

[中文说明](README.zh-CN.md)

## Install

Requires Python 3.10 or newer. Clone the complete repository into your Codex skills directory and install the dependencies:

```sh
git clone https://github.com/Sculptor815/go-enrichment-skill.git ~/.codex/skills/go-enrichment-skill
python -m pip install -r ~/.codex/skills/go-enrichment-skill/requirements.txt
```

If you use a custom `CODEX_HOME`, use its `skills` directory. Start a new conversation after installation.

## Use

Attach a gene-list file or paste the gene names, then ask:

> Use $go-enrichment-skill to run GO enrichment on the provided gene list, generate BP, MF and CC dot plots, and briefly explain the results.

Human is the default organism; specify another organism in your prompt when needed. The skill handles gene mapping, annotation preparation, enrichment and plotting. The first run may download large public annotation files; your gene list is analyzed locally.

> **Choose the background carefully:** for experimental data, prefer all genes that could have entered the query, such as all genes tested in the same differential-expression comparison, rather than significant genes alone. Without a supplied background, the skill uses the default GO-annotated universe. Background choice affects enrichment results.

## Example GO results

These plots were generated from an immune-biased synthetic gene list. They demonstrate the workflow and are not independent biological discoveries.

### BP · Biological process

![GO BP enrichment dot plot](examples/results/plots/immune_biased_BP.png)

### MF · Molecular function

![GO MF enrichment dot plot](examples/results/plots/immune_biased_MF.png)

### CC · Cellular component

![GO CC enrichment dot plot](examples/results/plots/immune_biased_CC.png)

The project also works as a Python CLI; see `python scripts/go_enrichment.py --help` for commands.

[MIT license](LICENSE) · [Data sources and licenses](NOTICE.md)
