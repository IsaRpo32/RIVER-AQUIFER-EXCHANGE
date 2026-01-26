"""
Visualization Module (Legacy)

NOTE: This is legacy code preserved for reference.
See src/modflow_utils.py for active plotting functions.

Comprehensive Visualization Module for MODFLOW Groundwater-Surface Water Interaction Models.
All measurements are in SI system (m, s, m3/s, m2, etc)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

from .settings import ModflowConfig


class Visualization(ModflowConfig):
    """
    Unified visualization class for MODFLOW simulation results.

    Provides comprehensive plotting capabilities for:
    - Model configuration and geometry
    - Hydrological data and time series
    - MODFLOW results (heads, flows, budgets)
    - SFR stream-aquifer interaction analysis
    """

    def __init__(self, figsize=(15, 8), style='default', dpi=150, verbose=True):
        """
        Initialize visualization with default settings.

        Parameters
        ----------
        figsize : tuple, optional
            Default figure size (default: (15, 8))
        style : str, optional
            Matplotlib style (default: 'default')
        dpi : int, optional
            Figure resolution (default: 150)
        verbose : bool, optional
            Control output verbosity (default: True)
        """
        super().__init__(verbose=verbose)

        self.figsize = figsize
        self.style = style
        self.dpi = dpi

        self.model = None
        self.results = None
        self.topography = None
        self.river_network = None
        self.i_riv = None
        self.j_riv = None
        self.nc_riv = None

        plt.style.use(self.style)
        plt.rcParams.update({
            'figure.dpi': self.dpi,
            'savefig.dpi': self.dpi,
            'font.size': 12,
        })

        if self.verbose:
            print("Visualization module initialized")


def plot_exchange_flow(t, y_raw, metrics, t_delay):
    """
    Plot raw vs smoothed exchange flow, shaded integrals, and mark peaks.
    """
    y_smooth = metrics['y_smooth']
    si, ni, zi = metrics['start_idx'], metrics['neg_idx'], metrics['zero_idx']
    ap, an = metrics['Qin_exc'], metrics['Qout_exc']
    ppi, ppi_val = metrics['peak_pos_idx'], metrics['peak_pos_val']
    pni, pni_val = metrics['peak_neg_idx'], metrics['peak_neg_val']
    ppt, pnt = metrics['peak_pos_time'], metrics['peak_neg_time']

    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'font.size': 12,
        'axes.linewidth': 1.2,
        'legend.fontsize': 11,
        'legend.frameon': False,
        'lines.linewidth': 1.5,
        'figure.dpi': 300,
    })

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(t, y_raw, color='gray', alpha=0.4, label='Raw signal')
    ax.plot(t, y_smooth, color='black', label='Smoothed signal')

    ax.fill_between(t[si:ni], y_smooth[si:ni],
                    where=y_smooth[si:ni]>=0,
                    color='orange', alpha=0.3,
                    label=f"Positive area\n({ap:.2e})")
    ax.fill_between(t[ni:zi], y_smooth[ni:zi],
                    where=y_smooth[ni:zi]<=0,
                    color='purple', alpha=0.3,
                    label=f"Negative area\n({an:.2e})")

    ax.axvline(t_delay, color='green', linestyle='--', label='t_delay')

    ax.scatter([ppt], [ppi_val],
               color='orange', edgecolor='black', zorder=5,
               label=f"Max +: {ppi_val:.2e} @ {ppt:.0f}")
    ax.scatter([pnt], [pni_val],
               color='purple', edgecolor='black', zorder=5,
               label=f"Max -: {pni_val:.2e} @ {pnt:.0f}")

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Exchange flow (Q)')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
