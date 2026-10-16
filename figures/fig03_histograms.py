"""
Fig. 3: distributions of log10 t_c and log10 pi1-pi4 (1,058 simulations; pi1 from the main ensemble only).
Also writes the table of statistics by drainage-area class (t_c, t_L and pi1-pi4):
data/output/tables/pi_statistics.csv.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import style as S

S.apply()
d = S.load()
VARS = ["tc", "pi1", "pi2", "pi3", "pi4"]
LABEL = {"tc": r"$\log_{10}\,t_c$ (s)", **{p: rf"$\log_{{10}}\,\pi_{p[-1]}$" for p in S.PI}}

H = 2.05
fig = plt.figure(figsize=(S.DOUBLE_WIDTH, H))
gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1, 0.1, 1], wspace=0.3, top=1 - 0.5 / H, bottom=0.24)
cols = [0, 1, 2, 3, 5]
axs = []
for k, v in enumerate(VARS):
    ax = fig.add_subplot(gs[0, cols[k]])
    axs.append(ax)
    x = np.log10(S.data_for(d, [v])[v])
    bins = np.linspace(x.min(), x.max(), 26)   # 25 bins in every panel
    ax.hist(x, bins=bins, color=S.DATA_GREY, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=2)
    ax.axvline(np.median(x), color=S.INK, lw=0.8, ls=(0, (3, 2)), zorder=3)
    ax.grid(axis="x", visible=False)
    ax.set_xlabel(LABEL[v])
    if k == 0:
        ax.set_ylabel("Number of simulations")
    ax.text(0.0, 1.02, f"({chr(97 + k)})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
            fontweight="bold", color=S.INK)

y_txt, y_line = 1 - 0.05 / H, 1 - 0.22 / H
for a, b, txt in ((axs[1], axs[3], r"System groups ($\pi_1$–$\pi_3$)"), (axs[4], axs[4], r"Forcing group ($\pi_4$)")):
    xa, xb = a.get_position().x0, b.get_position().x1
    fig.text((xa + xb) / 2, y_txt, txt, ha="center", va="top", fontsize=7, color=S.INK_2)
    fig.add_artist(plt.Line2D([xa, xb], [y_line, y_line], transform=fig.transFigure, color=S.FAINT, lw=0.6))
S.save(fig, "fig03_pi_histograms")

# ---- statistics by drainage-area class and for the whole ensemble
rows = []
for area_class, g in [("10^3 km2", d[d.area_class == "1e3"]), ("10^4 km2", d[d.area_class == "1e4"]),
                      ("10^5 km2", d[d.area_class == "1e5"]), ("All", d)]:
    for v in ["tc", "tL"] + S.PI:      # t_L (scale of delta_tp* and t_r*) only in the table, not in the figure
        x = g[v]
        rows.append(dict(area_class=area_class, n=len(g), variable=v, min=x.min(), median=x.median(), max=x.max(),
                         geometric_mean=10 ** np.log10(x).mean(), mean=x.mean()))
t = pd.DataFrame(rows)
S.TABLES.mkdir(parents=True, exist_ok=True)
t.to_csv(S.TABLES / "pi_statistics.csv", index=False)
print(t.assign(**{c: t[c].map(lambda z: f"{z:.2e}") for c in ["min", "median", "max", "geometric_mean", "mean"]})
      .to_string(index=False))
