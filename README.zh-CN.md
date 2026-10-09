# GO Enrichment Skill

一个用于 GO 富集分析的 Codex Skill。提供 gene list，即可生成 BP、MF、CC 气泡图、完整统计表和简要结果解释。

[English](README.md)

## 安装

需要 Python 3.10 或更新版本。将完整仓库放入 Codex 的 skills 目录，并安装依赖：

```sh
git clone https://github.com/Sculptor815/go-enrichment-skill.git ~/.codex/skills/go-enrichment-skill
python -m pip install -r ~/.codex/skills/go-enrichment-skill/requirements.txt
```

如果设置了自定义 `CODEX_HOME`，请使用其中的 `skills` 目录。安装后在新会话中使用。

## 使用

上传基因列表文件，或直接粘贴基因名，然后输入：

> 使用 $go-enrichment-skill 对提供的 gene list 进行 GO 富集分析，生成 BP、MF、CC 气泡图，并简要解释结果。

默认分析人类基因；其他物种请在提示词中说明。Skill 会处理基因映射、注释准备、富集分析和绘图。首次使用可能需要下载较大的公开注释文件，基因列表在本地分析。

> **请谨慎选择 background：** 真实实验建议使用所有有机会入选目标列表的基因，例如同一次差异分析中实际参与检验的全部基因，而非仅显著基因；未提供时使用默认 GO 注释背景，背景选择会影响富集结果。

## GO 分析示例

以下为免疫功能偏向的模拟基因列表实际运行得到的结果，用于展示分析效果，不代表独立的生物学发现。

### BP · 生物过程

![GO BP 富集气泡图](examples/results/plots/immune_biased_BP.png)

### MF · 分子功能

![GO MF 富集气泡图](examples/results/plots/immune_biased_MF.png)

### CC · 细胞组分

![GO CC 富集气泡图](examples/results/plots/immune_biased_CC.png)

本项目也可作为 Python CLI 使用，可通过 `python scripts/go_enrichment.py --help` 查看命令。

[MIT 许可](LICENSE) · [数据来源与许可](NOTICE.md)
