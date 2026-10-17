"""
Common style of the paper figures. Every figure script imports this module.

Double column 174 mm / single column 84 mm; sans-serif font (Arial, falling back to DejaVu Sans); thin
ticks; faint grid; data in a single grey; text always in neutral ink; vector output (PDF) plus PNG at 600 dpi.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from results_table import load_results  # noqa: E402

TABLES = ROOT / "data" / "output" / "tables"
OUTPUT = Path(__file__).resolve().parent / "output"

MM = 1 / 25.4
DOUBLE_WIDTH, SINGLE_WIDTH = 174 * MM, 84 * MM

# Inks and surfaces
INK, INK_2, FAINT, GRID = "#1b1b1a", "#52514e", "#898781", "#e7e6e2"

# Data: a single grey (no colouring by basin size in any figure)
DATA_GREY = "#5f5e5a"
DATA_ALPHA = 0.5

PI = ["pi1", "pi2", "pi3", "pi4"]
PI_LABEL = {"pi1": r"$\pi_1$", "pi2": r"$\pi_2$", "pi3": r"$\pi_3$", "pi4": r"$\pi_4$"}
METRICS = ["delta_tp*", "delta_H*", "t_r*", "delta_leakage*", "delta_storage*"]
METRIC_LABEL = {"delta_tp*": r"$\Delta t_p^{*}$", "delta_H*": r"$\Delta H^{*}$", "t_r*": r"$t_r^{*}$",
                "delta_leakage*": r"$\Delta\mathrm{Leakage}^{*}$", "delta_storage*": r"$\Delta\mathrm{Storage}^{*}$"}
# Metrics spanning less than a decade (scaled with t_L): explicit ticks on the log axis
METRIC_TICKS = {"delta_tp*": [1, 1.5, 2, 3], "t_r*": [15, 20, 30, 50, 80]}


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5, "legend.fontsize": 7,
        "axes.edgecolor": FAINT, "axes.linewidth": 0.6, "axes.labelcolor": INK,
        "xtick.color": FAINT, "ytick.color": FAINT, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.4, "ytick.minor.size": 1.4,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "legend.frameon": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })


def load():
    return load_results()


def data_for(d, variables):
    """Any figure that plots pi1 uses only the 758 simulations of the main ensemble."""
    return d[d.group == "main"] if "pi1" in variables else d


def log_axes(ax, x=True, y=True):
    for on, axis in ((x, ax.xaxis), (y, ax.yaxis)):
        if on:
            axis.set_major_locator(LogLocator(base=10, numticks=12))
            axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=12))
            axis.set_minor_formatter(NullFormatter())


def short_decade_ticks(axis, values):
    """For log axes spanning about 1-2 decades: explicit ticks with decimal labels."""
    from matplotlib.ticker import FixedLocator, FuncFormatter
    axis.set_major_locator(FixedLocator(values))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


def spearman(x, y):
    a, b = pd.Series(np.log10(x)).rank(), pd.Series(np.log10(y)).rank()
    return float(np.corrcoef(a, b)[0, 1])


def save(fig, name, crop=True):
    """crop=False for 3D figures: the automatic crop ignores the 3D axis labels."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if crop:
        extra = {"bbox_inches": "tight", "pad_inches": 0.02}
    else:   # crop to the actual content: test render, find the non-white area and use that box
        import io
        from matplotlib.transforms import Bbox
        dpi = 150
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=dpi); buf.seek(0)
        img = plt.imread(buf)[..., :3]
        rows, cols = np.where(img.min(axis=2) < 0.97)
        H, W = img.shape[:2]; pad = 0.04
        x0, x1 = cols.min() / dpi - pad, (cols.max() + 1) / dpi + pad
        y0, y1 = (H - rows.max() - 1) / dpi - pad, (H - rows.min()) / dpi + pad
        extra = {"bbox_inches": Bbox([[max(x0, 0), max(y0, 0)], [min(x1, W / dpi), min(y1, H / dpi)]])}
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
        fig.savefig(OUTPUT / f"{name}.{ext}", **extra, **kw)
    print("->", OUTPUT / f"{name}.png")
