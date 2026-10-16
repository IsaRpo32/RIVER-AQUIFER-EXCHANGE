"""
Fig. 9: the Nechí corridor over the 1,058 simulations, each metric against pi3, in one row of five panels.
Data in grey; the simulated Nechí value as a filled accent marker and the value predicted by the pi1-pi4
model as a hollow marker. No pi1 on the axes, so all 1,058 simulations are shown.
Nechí values: data/output/tables/nechi_position.csv (analysis/nechi_position.py).
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import style as S

S.apply()
d = S.load()
nq = pd.read_csv(S.TABLES / "nechi_position.csv").set_index("variable")
ACCENT = "#D55E00"   # vermilion (Okabe-Ito): colour-blind safe and distinct from the data grey

fig, axs = plt.subplots(1, 5, figsize=(S.DOUBLE_WIDTH, 2.2), gridspec_kw=dict(wspace=0.55))
fig.subplots_adjust(bottom=0.31, top=0.9)
for k, (ax, y) in enumerate(zip(axs, S.METRICS)):
    ax.scatter(d.pi3, d[y], s=4, facecolor=S.DATA_GREY, edgecolor="white", linewidths=0.15, alpha=S.DATA_ALPHA,
               zorder=2, rasterized=True)
    ax.scatter(nq.loc["pi3", "nechi"], nq.loc[y, "nechi"], s=26, marker="o", facecolor=ACCENT, edgecolor="white",
               linewidths=0.5, zorder=5)
    ax.scatter(nq.loc["pi3", "nechi"], nq.loc[y, "prediction"], s=22, marker="o", facecolor="none",
               edgecolor=ACCENT, linewidths=0.9, zorder=4)
    ax.set_xscale("log"); ax.set_yscale("log"); S.log_axes(ax)
    ax.set_xticks([1e4, 1e6, 1e8])
    if y in S.METRIC_TICKS:
        S.short_decade_ticks(ax.yaxis, S.METRIC_TICKS[y])
    ax.set_xlabel(S.PI_LABEL["pi3"])
    ax.set_ylabel(S.METRIC_LABEL[y], labelpad=1)
    ax.text(0.0, 1.03, f"({chr(97 + k)})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
            fontweight="bold", color=S.INK)
h = [Line2D([], [], ls="", marker="o", ms=4, mfc=S.DATA_GREY, mec="white", alpha=0.7, label="Simulations (n = 1,058)"),
     Line2D([], [], ls="", marker="o", ms=5.5, mfc=ACCENT, mec="white", label="Nechí, simulated"),
     Line2D([], [], ls="", marker="o", ms=5, mfc="none", mec=ACCENT, mew=0.9, label="Nechí, predicted (Eq. 22)")]
fig.legend(handles=h, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), labelcolor=S.INK_2,
           handletextpad=0.3, columnspacing=1.8)
S.save(fig, "fig09_nechi_pi_space")
