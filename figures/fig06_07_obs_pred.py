"""
Fig. 6 (delta_tp*, delta_H*) and Fig. 7 (t_r*, delta_leakage*, delta_storage*): simulated vs predicted by the
log-log model with pi1-pi4. Coefficients and R2 from data/output/tables/predictive_models_pi1_pi4.csv
(analysis/regressions.py).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import style as S

S.apply()
d = S.load()
tp = pd.read_csv(S.TABLES / "predictive_models_pi1_pi4.csv")


def panel(ax, met, letter):
    c = tp[tp.metric == met].set_index("term")
    pred = 10 ** (c.loc["a0", "coef"] + sum(c.loc[p, "coef"] * np.log10(d[p]) for p in S.PI))
    obs = d[met]
    margin = 1.12 if met == "delta_tp*" else 1.6     # delta_tp* spans half a decade: narrow margin
    lo, hi = min(obs.min(), pred.min()) / margin, max(obs.max(), pred.max()) * margin
    ax.plot([lo, hi], [lo, hi], color=S.INK, lw=0.8, zorder=3)
    ax.scatter(obs, pred, s=6, facecolor=S.DATA_GREY, edgecolor="white", linewidths=0.25, alpha=S.DATA_ALPHA,
               zorder=2, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    if met == "delta_tp*":       # y axis up to the tick 2 (predictions do not exceed 1.9)
        ax.set_ylim(lo, 2.08); ax.set_aspect("auto"); ax.set_box_aspect(1)
    S.log_axes(ax)
    if met in S.METRIC_TICKS:
        ticks = [0.8] + S.METRIC_TICKS[met] if met == "delta_tp*" else S.METRIC_TICKS[met]
        for axis in (ax.xaxis, ax.yaxis):
            S.short_decade_ticks(axis, ticks)
    r = c.iloc[0]
    ax.text(0.96, 0.04, rf"$R^2_{{\mathrm{{adj}}}}$ = {r.R2adj:.3f}" + "\n" + rf"$R^2_{{\mathrm{{cv}}}}$ = {r.R2cv:.3f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color=S.INK_2, linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.15", fc="white", ec="none", alpha=0.85), zorder=5)
    ax.set_xlabel("Simulated " + S.METRIC_LABEL[met]); ax.set_ylabel("Predicted " + S.METRIC_LABEL[met])
    ax.text(0.0, 1.02, f"({letter})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
            fontweight="bold", color=S.INK)


for name, mets, width in (("fig06_obs_pred_regulation", ["delta_tp*", "delta_H*"], S.DOUBLE_WIDTH * 0.66),
                          ("fig07_obs_pred_processes", ["t_r*", "delta_leakage*", "delta_storage*"], S.DOUBLE_WIDTH)):
    fig, axs = plt.subplots(1, len(mets), figsize=(width, 2.55), gridspec_kw=dict(wspace=0.42))
    for k, (ax, m) in enumerate(zip(np.atleast_1d(axs), mets)):
        panel(ax, m, chr(97 + k))
    S.save(fig, name)
