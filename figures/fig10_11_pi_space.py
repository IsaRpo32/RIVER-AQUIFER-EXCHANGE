"""
Fig. 10: delta_H* in the space of the system groups (pi1, pi2, pi3).
   - 3D version (paper): axes log10 pi_x, log10 pi_y and delta_H*; delta_H* >= 0.10 highlighted.
   - 2D alternative: pi_x-pi_y planes with delta_H* as a continuous grey scale (FIG_COLOR=1 for a blue scale).
Fig. 11a: pi1-pi2-pi3 space (3D) coloured by log10 pi4 (continuous).
Panels with pi1 use the main ensemble only (758); the others all 1,058 simulations.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.ticker import MaxNLocator, FixedLocator, FuncFormatter
import style as S

S.apply()
d = S.load()
L = {p: np.log10(d[p]).to_numpy() for p in ["pi1", "pi2", "pi3", "pi4", "delta_H*"]}
ALL = np.ones(len(d), bool)
MAIN = (d.group == "main").to_numpy()          # if pi1 is plotted, only the 758 of the main ensemble
mask = lambda *vs: MAIN if "pi1" in vs else ALL
COLOR = os.environ.get("FIG_COLOR") == "1"
SUFFIX = "_color" if COLOR else ""
AXIS = {"pi1": r"$\log_{10}\,\pi_1$", "pi2": r"$\log_{10}\,\pi_2$", "pi3": r"$\log_{10}\,\pi_3$"}


def axes_3d(ax, elev=18, azim=-58):
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.fill = False
        a.pane.set_edgecolor(S.GRID)
        a._axinfo["grid"].update(color=S.GRID, linewidth=0.45)
        a._axinfo["tick"].update(inward_factor=0.0, outward_factor=0.25)
        a.line.set_color(S.FAINT); a.line.set_linewidth(0.6)
    ax.tick_params(labelsize=6, pad=-1, colors=S.INK_2)
    ax.view_init(elev=elev, azim=azim)


# ------------------------------------------------------------------ Fig. 10 (3D)
# delta_H* vertical; in each panel one pi in front (x axis) and another in depth (y axis). Nearly frontal view of
# the x-z plane, to show how the cloud organizes when pi3 is on the front axis (panels b and c).
# Two classes: delta_H* >= 0.10 (peak reduction >= 10 % of d_bkf) in light blue; the rest in grey.
BLUE_10, GREY_REST = "#6da7ec", "#8e8d88"
cl_dh = (d["delta_H*"].to_numpy() >= 0.10).astype(int)
panels = [("pi2", "pi1"), ("pi3", "pi1"), ("pi3", "pi2")]    # (front, depth)
fig = plt.figure(figsize=(S.DOUBLE_WIDTH, 2.45))
z = L["delta_H*"]
for k, (px, py) in enumerate(panels):
    ax = fig.add_subplot(1, 3, k + 1, projection="3d", computed_zorder=False)
    axes_3d(ax, elev=14, azim=-64)
    m = mask(px, py)
    x, y, zz, cl = L[px][m], L[py][m], z[m], cl_dh[m]
    lim = [(v.min(), v.max()) for v in (x, y, zz)]
    pad = [0.04 * (b - a) for a, b in lim]
    lo = [a - p for (a, _), p in zip(lim, pad)]; hi = [b + p for (_, b), p in zip(lim, pad)]
    for kk, col in ((0, GREY_REST), (1, BLUE_10)):          # blue on top
        s = cl == kk
        ax.scatter(x[s], y[s], zz[s], s=6, c=col, edgecolors="white", linewidths=0.15, alpha=0.95,
                   depthshade=False, zorder=2 + kk)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.yaxis.set_major_locator(MaxNLocator(3)); ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.zaxis.set_major_locator(FixedLocator([-3, -2, -1]))
    ax.zaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{10.0 ** v:g}"))   # delta_H* values (log spacing)
    ax.set_xlabel(AXIS[px], labelpad=-3, fontsize=8); ax.set_ylabel(AXIS[py], labelpad=-1, fontsize=8)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r"$\Delta H^{*}$", labelpad=0, fontsize=8.5, rotation=0)
    ax.set_box_aspect((1.6, 0.9, 0.9), zoom=0.8)
    ax.set_facecolor("none")          # a panel background must not hide the delta_H* label of its neighbour
    ax.text2D(0.08, 0.86, f"({chr(97 + k)})", transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=S.INK)
fig.subplots_adjust(left=0.0, right=0.95, top=1.1, bottom=0.1, wspace=0.14)
fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, mfc=BLUE_10, mec="white", mew=0.3,
                           label=r"$\Delta H^{*} \geq 0.10$"),
                    Line2D([], [], ls="", marker="o", ms=4.5, mfc=GREY_REST, mec="white", mew=0.3,
                           label=r"$\Delta H^{*} < 0.10$")],
           loc="lower center", bbox_to_anchor=(0.5, 0.14), ncol=2, labelcolor=S.INK_2, handletextpad=0.3,
           columnspacing=2.0)
S.save(fig, "fig10_deltaH_pi_space_3D", crop=False)

# ------------------------------------------------------------------ Fig. 10 (2D alternative)
cmap = LinearSegmentedColormap.from_list(
    "seq", ["#dce9f9", "#5598e7", "#1c5cab", "#0d2c55"] if COLOR else ["#e4e3df", "#8e8d88", "#1f1e1c"])
fig, axs = plt.subplots(1, 3, figsize=(S.DOUBLE_WIDTH, 2.45), gridspec_kw=dict(wspace=0.32))
pairs = [("pi1", "pi2"), ("pi1", "pi3"), ("pi2", "pi3")]
norm = LogNorm(d["delta_H*"].min(), d["delta_H*"].max())     # same colour scale in the three panels
for k, (ax, (px, py)) in enumerate(zip(axs, pairs)):
    s = d[mask(px, py)].sort_values("delta_H*")              # high values on top
    sc = ax.scatter(s[px], s[py], c=s["delta_H*"], cmap=cmap,
                    norm=norm, s=7, edgecolors="white", linewidths=0.2, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log"); S.log_axes(ax)
    if px == "pi2":
        S.short_decade_ticks(ax.xaxis, [0.02, 0.05, 0.1, 0.2, 0.5])
    if py == "pi2":
        S.short_decade_ticks(ax.yaxis, [0.02, 0.05, 0.1, 0.2, 0.5])
    ax.set_xlabel(S.PI_LABEL[px]); ax.set_ylabel(S.PI_LABEL[py])
    ax.text(0.0, 1.02, f"({chr(97 + k)})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
            fontweight="bold", color=S.INK)
cb = fig.colorbar(sc, ax=axs, fraction=0.025, pad=0.02)
cb.set_label(r"$\Delta H^{*}$"); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6.5)
S.save(fig, "fig10_deltaH_pi_space_2D" + SUFFIX)

# ------------------------------------------------------------------ Fig. 11a (3D, continuous pi4)
m = mask("pi1")
x1, x2, x3, c4 = L["pi1"][m], L["pi2"][m], L["pi3"][m], L["pi4"][m]
order = np.argsort(c4)
cmap4 = LinearSegmentedColormap.from_list("pi4", ["#dce9f9", "#86b6ef", "#2a78d6", "#104281", "#0a2a52"])
fig = plt.figure(figsize=(S.SINGLE_WIDTH * 1.35, 3.3))
ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
axes_3d(ax, elev=20, azim=-55)
sc = ax.scatter(x1[order], x2[order], x3[order], c=c4[order], cmap=cmap4, norm=Normalize(c4.min(), c4.max()), s=6,
                edgecolors="white", linewidths=0.15, depthshade=False)
ax.set_xlim(x1.min(), x1.max()); ax.set_ylim(x2.min(), x2.max()); ax.set_zlim(x3.min(), x3.max())
ax.set_xlabel(AXIS["pi1"], labelpad=-5, fontsize=7.5); ax.set_ylabel(AXIS["pi2"], labelpad=-5, fontsize=7.5)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(AXIS["pi3"], labelpad=-3, fontsize=7.5, rotation=90)
ax.set_box_aspect((1, 1, 0.9), zoom=0.8)
cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02, shrink=0.6)
cb.set_label(r"$\log_{10}\,\pi_4$", fontsize=7.5); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6.5)
fig.subplots_adjust(left=0.0, right=0.92, top=1.0, bottom=0.04)
S.save(fig, "fig11a_pi_space_pi4_continuous", crop=False)
