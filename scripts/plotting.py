"""Headless GO dot plots adapted from the original GEP workflow."""
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np

SOURCE_NAMES = {"BP": "Biological process", "MF": "Molecular function", "CC": "Cellular component"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white"})

def csv_save(frame, path):
    frame.to_csv(path, index=False, encoding="utf-8")

def dotplot(df, gep, source, output, query_n, fdr=0.05, top_terms=15, dpi=180, input_n=None):
    data = df.loc[df["significant"].astype(bool)].sort_values(
        ["p_adjust_BH", "Count", "GO_ID"], ascending=[True, False, True]).head(top_terms)
    data = data.sort_values(["GeneRatio", "p_adjust_BH"], ascending=[False, True])
    csv_save(data, output.with_name(output.name + "_plotted_terms.csv"))
    title = f"{gep}  |  GO {SOURCE_NAMES[source]}"
    if data.empty:
        fig, ax = plt.subplots(figsize=(9, 3.4))
        ax.axis("off")
        msg = "No mapped query genes" if query_n == 0 else f"No enriched terms at BH-adjusted P < {fdr:g}"
        ax.text(0.5, 0.65, title, ha="center", fontsize=14, transform=ax.transAxes)
        ax.text(0.5, 0.4, msg, ha="center", color="#666666", transform=ax.transAxes)
    else:
        labels = [textwrap.fill(f"{name} ({go})", width=48, break_long_words=False)
                  for name, go in zip(data.description, data.GO_ID)]
        # 为多行名称分配真实高度，避免相邻标签重叠。
        heights = np.array([max(1, x.count("\n") + 1) for x in labels], float)
        y = np.cumsum(heights + 0.8) - (heights + 0.8) / 2
        fig_h = max(5.2, float((heights + 0.8).sum()) * 0.23 + 1.6)
        fig = plt.figure(figsize=(12.8, fig_h), layout="constrained")
        grid = fig.add_gridspec(1, 2, width_ratios=[1, 0.22], wspace=0.08)
        ax = fig.add_subplot(grid[0, 0])
        side = grid[0, 1].subgridspec(2, 2, width_ratios=[0.18, 0.82], height_ratios=[1, 1])
        cax = fig.add_subplot(side[0, 0])
        lax = fig.add_subplot(side[1, :])
        lax.axis("off")
        color = -np.log10(np.maximum(data.p_adjust_BH.to_numpy(float), 1e-300))
        lo, hi = float(color.min()), float(color.max())
        if hi - lo < 1e-8:
            lo, hi = lo - 0.15, hi + 0.15
        size_scale = 430.0 / float(data.Count.max())
        points = ax.scatter(data.GeneRatio, y, s=data.Count * size_scale,
            c=color, cmap="viridis", norm=Normalize(lo, hi),
            edgecolors="#36454F", linewidths=0.65, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=10)
        ax.tick_params(axis="y", length=0, pad=10)
        ax.set_ylim(float(y[-1] + heights[-1] / 2 + 0.8), float(y[0] - heights[0] / 2 - 0.8))
        ax.set_xlim(0, max(0.05, float(data.GeneRatio.max()) * 1.16))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
        ax.grid(axis="x", color="#E7EBEF", linewidth=0.8, zorder=0)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color("#9AA4AD")
        ax.set_xlabel("Gene ratio (overlap / mapped query genes)", labelpad=10)
        ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=14)
        cb = fig.colorbar(points, cax=cax)
        cb.set_label(r"$-\log_{10}$(BH-adjusted P)", fontsize=10)
        cb.ax.tick_params(labelsize=9)
        counts = np.unique(np.rint(np.quantile(data.Count, [0, 0.5, 1])).astype(int))
        handles = [ax.scatter([], [], s=v * size_scale, facecolors="white",
                   edgecolors="#36454F", linewidths=0.8) for v in counts]
        lax.legend(handles, [str(v) for v in counts], title="Overlap genes",
                   loc="center left", frameon=False, labelspacing=1.2)
        fig.supxlabel(f"Top {len(data)} enriched terms; input IDs={input_n if input_n is not None else query_n:,}; "
            f"mapped query={int(data.query_size.iloc[0]):,}; "
            f"effective background={int(data.effective_background_size.iloc[0]):,}",
            fontsize=9, color="#59636E")
    try:
        fig.savefig(str(output) + ".png", dpi=dpi, bbox_inches="tight", pad_inches=0.2)
        fig.savefig(str(output) + ".pdf", bbox_inches="tight", pad_inches=0.2)
    finally:
        plt.close(fig)
