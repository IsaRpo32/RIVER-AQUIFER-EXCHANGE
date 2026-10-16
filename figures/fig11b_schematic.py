"""
Fig. 11b: conceptual schematic of flood regulation in pi-space (rebuilt from the author's original schematic,
whose native file was lost). pi3 vertical, pi2 to the right, pi1 diagonal; a wedge along each axis (wider
towards the high value, and wider for the higher-ranked group: pi3 > pi1 > pi2); pi4 box as the (unranked)
forcing; box with the direction of increasing regulation (lower pi3, higher pi1, lower pi2).
No "Mixed" classes (continuous gradient, no thresholds). All labels are edited here.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
import style as S

S.apply()
WEDGE, TXT, TXT2 = "#a7ad9c", S.INK, S.INK_2

fig = plt.figure(figsize=(S.DOUBLE_WIDTH, 4.3))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 15.5); ax.set_ylim(-0.6, 10.4); ax.axis("off")

O = np.array([5.0, 3.6])
AXES = {"pi3": np.array([5.0, 9.6]), "pi2": np.array([11.6, 3.6]), "pi1": np.array([1.4, -0.2])}


def axis_arrow(end, label, off):
    ax.add_patch(FancyArrowPatch(O, end, arrowstyle="-|>", mutation_scale=11, lw=1.0, color=TXT, zorder=4))
    ax.text(*(end + off), label, fontsize=11, color=TXT, ha="center", va="center")


def wedge(p0, p1, max_width, normal, text, rot, shift):
    """Triangle along the axis, from width 0 at p0 to max_width at p1, shifted along the normal."""
    n = normal / np.linalg.norm(normal)
    a, b = p0 + n * 0.12, p1 + n * 0.12
    ax.add_patch(Polygon([a, b, b + n * max_width], closed=True, fc=WEDGE, ec="none", alpha=0.85, zorder=2))
    c = (a + b) / 2 + n * (max_width / 2 + shift)
    ax.text(*c, text, fontsize=7.2, color=TXT2, style="italic", rotation=rot, ha="center", va="center",
            rotation_mode="anchor")


def dot(p):
    ax.plot(*p, "o", ms=5, color=TXT, zorder=5)


axis_arrow(AXES["pi3"], r"$\pi_3$", np.array([0, 0.45]))
axis_arrow(AXES["pi2"], r"$\pi_2$", np.array([0.45, 0]))
axis_arrow(AXES["pi1"], r"$\pi_1$", np.array([-0.3, -0.35]))

# --- pi3: primary control (widest wedge)
wedge(O + np.array([0, 0.5]), AXES["pi3"] - np.array([0, 0.5]), 0.42, np.array([1, 0]),
     "kinematic contrast (primary control)", 90, 0.32)
p_high, p_low = O + np.array([0, 5.0]), O + np.array([0, 1.05])
dot(p_high); dot(p_low)
ax.text(p_high[0] - 0.25, p_high[1], "Flood-transmitting", fontsize=8.5, fontweight="bold", color=TXT, ha="right",
        va="bottom")
ax.text(p_high[0] - 0.25, p_high[1] - 0.08, "wave crosses the reach\nbefore the aquifer responds", fontsize=7.5,
        color=TXT2, ha="right", va="top", linespacing=1.25)
ax.text(p_low[0] - 0.25, p_low[1] + 0.05, "Flood-regulating", fontsize=8.5, fontweight="bold", color=TXT,
        ha="right", va="bottom")
ax.text(p_low[0] - 0.25, p_low[1] - 0.03, "aquifer responds while\nthe wave is in the reach", fontsize=7.5,
        color=TXT2, ha="right", va="top", linespacing=1.25)

# --- pi2: third control (narrowest wedge)
wedge(O + np.array([0.6, 0]), AXES["pi2"] - np.array([0.6, 0]), 0.2, np.array([0, 1]),
     "storage partitioning (third control)", 0, 0.25)
q_low, q_high = O + np.array([1.3, 0]), O + np.array([5.4, 0])
dot(q_low); dot(q_high)
for q, t in ((q_low, "Aquifer-buffered\nstorage"), (q_high, "Channel-controlled\nstorage")):
    ax.text(q[0], q[1] - 0.25, t, fontsize=7.5, color=TXT2, ha="center", va="top", linespacing=1.25)

# --- pi1: second control
d1 = (AXES["pi1"] - O) / np.linalg.norm(AXES["pi1"] - O)
n1 = np.array([d1[1], -d1[0]])               # normal towards the upper left of the axis
ang = np.degrees(np.arctan2(d1[1], d1[0])) + 180
wedge(O + d1 * 0.6, AXES["pi1"] - d1 * 0.5, 0.3, n1, "vertical exchange efficiency (second control)", ang, 0.28)
r_low, r_high = O + d1 * 1.6, O + d1 * 4.3
dot(r_low); dot(r_high)
for r, t in ((r_low, "Weak streambed\nexchange"), (r_high, "Strong streambed\nexchange")):
    ax.text(r[0] + 0.3, r[1] - 0.05, t, fontsize=7.5, color=TXT2, ha="left", va="top", linespacing=1.25)

# --- increasing-regulation box: direction for each group, in the order of the hierarchy
BLUE = "#3e5f78"
xk, yk = 11.75, 2.55
ax.plot([xk, xk], [yk - 1.95, yk], color=BLUE, lw=0.9)
for yy in (yk, yk - 1.95):
    ax.plot([xk, xk + 0.18], [yy, yy], color=BLUE, lw=0.9)
ax.text(xk + 0.3, yk - 0.05, "Increasing flood regulation", fontsize=8.6, color=BLUE, fontweight="bold",
        ha="left", va="top")
for k, (arrow, t) in enumerate(((r"$\downarrow$", r"lower $\pi_3$ (primary)"), (r"$\uparrow$", r"higher $\pi_1$ (second)"),
                                 (r"$\downarrow$", r"lower $\pi_2$ (third)"))):
    yy = yk - 0.75 - 0.5 * k
    ax.text(xk + 0.42, yy, arrow, fontsize=9, color=BLUE, ha="center", va="center")
    ax.text(xk + 0.7, yy, t, fontsize=7.8, color=TXT2, ha="left", va="center")

# --- pi4 box (forcing)
x0, y0 = 11.75, 9.65
ax.plot([x0, x0], [y0 - 3.75, y0], color=TXT, lw=0.9)
for yy in (y0, y0 - 3.75):
    ax.plot([x0, x0 + 0.18], [yy, yy], color=TXT, lw=0.9)
ax.text(x0 + 0.3, y0 - 0.05, r"$\pi_4$" + " forcing (not ranked)", fontsize=8.8, color=TXT, ha="left", va="top")
ax.text(x0 + 0.3, y0 - 0.75, "High " + r"$\pi_4$" + ": event short relative\nto aquifer response", fontsize=7.6,
        color=TXT2, ha="left", va="top", linespacing=1.3)
ax.text(x0 + 0.3, y0 - 1.75, "Low " + r"$\pi_4$" + ": event long enough for\nexchange, storage and\ndelayed return flow",
        fontsize=7.6, color=TXT2, ha="left", va="top", linespacing=1.3)
ax.text(x0 + 0.3, y0 - 3.0, "Controls exchanged and stored\nvolumes; little effect on peak\nattenuation", fontsize=7.6,
        color=TXT, style="italic", ha="left", va="top", linespacing=1.3)

S.save(fig, "fig11b_regulation_schematic")
fig.savefig(S.OUTPUT / "fig11b_regulation_schematic.svg", bbox_inches="tight", pad_inches=0.02)   # editable text (svg.fonttype = none)
