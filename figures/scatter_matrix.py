"""
Log-log scatter matrix, metric vs pi1-pi4 (Figs. 4 and 5): grey points, Spearman rho in each panel.
The columns pi1-pi3 (system groups) and pi4 (forcing) are separated by a gap and a label.
"""
import matplotlib.pyplot as plt
import style as S

TICKS_PI2 = [0.02, 0.05, 0.1, 0.2, 0.5]


def matrix(d_all, rows, name, row_height=1.72):
    H = row_height * len(rows) + 0.7
    fig = plt.figure(figsize=(S.DOUBLE_WIDTH, H))
    gs = fig.add_gridspec(len(rows), 5, width_ratios=[1, 1, 1, 0.12, 1], wspace=0.08, hspace=0.2,
                          top=1 - 0.5 / H)
    cols = [0, 1, 2, 4]
    axs = [[None] * 4 for _ in rows]
    for i, met in enumerate(rows):
        for j, pi in enumerate(S.PI):
            share_y = axs[i][0] if j else None
            share_x = axs[0][j] if i else None
            ax = fig.add_subplot(gs[i, cols[j]], sharey=share_y, sharex=share_x)
            axs[i][j] = ax
            d = S.data_for(d_all, [pi])
            ax.scatter(d[pi], d[met], s=6, marker="o", facecolor=S.DATA_GREY, edgecolor="white", linewidths=0.25,
                       alpha=S.DATA_ALPHA, zorder=2, rasterized=True)
            ax.set_xscale("log"); ax.set_yscale("log"); S.log_axes(ax)
            if pi == "pi2":
                S.short_decade_ticks(ax.xaxis, TICKS_PI2)
            if met in S.METRIC_TICKS:
                S.short_decade_ticks(ax.yaxis, S.METRIC_TICKS[met])
            ax.text(0.97, 0.96, rf"$\rho$ = {S.spearman(d[pi], d[met]):+.2f}", transform=ax.transAxes, ha="right",
                    va="top", fontsize=6.5, color=S.INK_2,
                    bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.15", fc="white", ec="none", alpha=0.85),
                    zorder=5)
            ax.text(0.0, 1.015, f"({chr(97 + 4 * i + j)})", transform=ax.transAxes, ha="left", va="bottom",
                    fontsize=8.5, fontweight="bold", color=S.INK)
            if i < len(rows) - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(S.PI_LABEL[pi])
            if j:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(S.METRIC_LABEL[met])
    # block labels above the first row
    y_txt, y_line = 1 - 0.05 / H, 1 - 0.24 / H
    x0 = axs[0][0].get_position().x0; x1 = axs[0][2].get_position().x1; x4 = axs[0][3].get_position()
    for xa, xb, txt in ((x0, x1, r"System groups ($\pi_1$–$\pi_3$)"), (x4.x0, x4.x1, r"Forcing group ($\pi_4$)")):
        fig.text((xa + xb) / 2, y_txt, txt, ha="center", va="top", fontsize=7, color=S.INK_2)
        fig.add_artist(plt.Line2D([xa, xb], [y_line, y_line], transform=fig.transFigure, color=S.FAINT, lw=0.6))
    S.save(fig, name)
    return fig
