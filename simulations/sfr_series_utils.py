"""
Helpers to rebuild the time vector of a finished run and read the depth / flow series of the
two reaches used for delta_tp, H_ratio and delta_H (2nd and 2nd-to-last reach).
The time discretization is recomputed with the same code path as run_single_simulation.
"""
import os
import numpy as np
import pandas as pd
from flopy.utils.sfroutputfile import SfrFile

from parallel_modflow import (DEFAULT_TIME_CFG, build_adaptive_periods, duration_by_basin,
                              step_end_times)
from src import dimensionless_numbers as dn
from src.modflow_utils import hydrograph_cooper_rorabaugh


def recompute_times(params: pd.Series, cfg: dict) -> np.ndarray:
    """End time [s] of every transient step, measured from the start of the transient period."""
    cfg = {**DEFAULT_TIME_CFG, **cfg}
    d_bkf, w_bkf, L, B, D = params["d_bkf"], params["w_bkf"], params["L"], params["B"], params["D"]
    Ka, Sy, Slope, n = params["Ka"], params["Sy"], params["S"], params["n"]
    H_0, t_d = params["H_0"], params["t_d_days"] * 86400
    nrow, ncol = int(params["nrow"]), 65
    delc, delr = L / nrow, B / ncol

    t = np.arange(0, t_d, 600)
    H_b = hydrograph_cooper_rorabaugh(t, H_0, (d_bkf - H_0) / H_0, 0.5, t_d)
    Q = (1 / n) * H_b**(5 / 3) * Slope**0.5 * w_bkf
    v_max = Q.max() / (w_bkf * d_bkf)
    reach_len = L / nrow
    duration = (cfg['duration_days'] or duration_by_basin(params["da_km2"])) * 86400

    if cfg['mode'] == 'courant':
        tc = dn.get_tc(Sy=Sy, Ka=Ka, B=B, d_bkf=d_bkf)
        tc3 = Sy * (delr * delc) / (4 * Ka * D)
        rates = 1 / np.array([L / v_max, tc, t[np.argmax(Q)], tc3])
        delta_t_TC = 0.25 * len(rates) / np.sum(rates)
        cte = {1000: .5, 10000: .35}.get(params["da_km2"], .3)
        delta_t = min(delta_t_TC, cfg['courant_factor'] * cte * reach_len / v_max)
        perlen, nstp, tsmult = [duration], [int(np.ceil(duration / delta_t))], [1.0]
    else:
        perlen, nstp, tsmult = build_adaptive_periods(duration, 86400, t_d, v_max, reach_len, 50, cfg)
    return step_end_times(perlen, nstp, tsmult)


def read_reach_series(run_dir: str, nrow: int) -> dict:
    """Depth and outflow of reach 2 (in) and reach nrow-1 (out), transient steps only."""
    fname = [f for f in os.listdir(run_dir) if f.endswith(".sfr.out")][0]
    df = SfrFile(os.path.join(run_dir, fname)).df
    out = {}
    for key, reach in (("in", 2), ("out", nrow - 1)):
        r = df[df["reach"] == reach]
        out[f"H_{key}"] = r["depth"].to_numpy()[1:]
        out[f"Q_{key}"] = r["Qout"].to_numpy()[1:]
        out[f"Qin_{key}"] = r["Qin"].to_numpy()[1:]
    return out
