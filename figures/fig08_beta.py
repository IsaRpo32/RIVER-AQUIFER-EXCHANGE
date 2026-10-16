"""
Fig. 8: standardized coefficients (beta) of the pi1-pi3 model for each metric, with percentile bootstrap 95 %
CIs (3,000 replicates). Fixed order pi1, pi2, pi3 in every panel.
Data: data/output/tables/hierarchy_beta_pi1_pi3.csv (analysis/regressions.py).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import style as S

S.apply()
tb = pd.read_csv(S.TABLES / "hierarchy_beta_pi1_pi3.csv")
PI3 = ["pi1", "pi2", "pi3"]

fig, axs = plt.subplots(1, len(S.METRICS), figsize=(S.DOUBLE_WIDTH, 1.95), sharey=True, sharex=True,
                        gridspec_kw=dict(wspace=0.16))
y = np.arange(len(PI3))[::-1]
for k, (ax, met) in enumerate(zip(axs, S.METRICS)):
    s = tb[tb.metric == met].set_index("pi").loc[PI3]
    ax.axvline(0, color=S.FAINT, lw=0.7, zorder=1)
    ax.barh(y, s.beta, height=0.52, color=S.DATA_GREY, alpha=0.85, zorder=2)
    ax.errorbar(s.beta, y, xerr=[s.beta - s.ci_boot_low, s.ci_boot_high - s.beta], fmt="none", ecolor=S.INK,
                elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=3)
    for yi, b in zip(y, s.beta):
        # value label above the bar, centred on it, clear of the whiskers and of the zero line
        ax.text(b / 2, yi + 0.36, f"{b:+.2f}".replace("-", "−"), va="bottom", ha="center", fontsize=6,
                color=S.INK_2)
    ax.set_xlim(-1.15, 0.6)
    ax.set_xticks([-1, -0.5, 0, 0.5], ["−1", "−0.5", "0", "0.5"])
    ax.set_ylim(-0.55, 2.75)
    ax.grid(axis="y", visible=False)
    ax.set_title(S.METRIC_LABEL[met] + "\n" + rf"$R^2_{{\mathrm{{adj}}}}$ = {s.R2adj_pi123.iloc[0]:.2f}", fontsize=7.2,
                 color=S.INK, pad=4, linespacing=1.4)
    ax.text(0.0, 1.21, f"({chr(97 + k)})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
            fontweight="bold", color=S.INK)
    ax.set_xlabel(r"$\beta$")
axs[0].set_yticks(y, [S.PI_LABEL[p] for p in PI3])
axs[0].tick_params(axis="y", length=0, labelsize=8)
S.save(fig, "fig08_beta_pi1_pi3")
