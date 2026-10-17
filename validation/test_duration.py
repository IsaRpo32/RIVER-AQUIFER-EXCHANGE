"""
How long must the transient run so that the net river-aquifer exchange approaches zero?
Long runs (730 d) of rows spanning low / median / high tc per basin, 'fine' scheme, corrected n.
For each run: time at which |Qexc(t)| (aggregated over reaches) stays below 5 %, 2 % and 1 % of its
positive peak, and the fraction of the final returned volume already returned at 90/120/150/365 d.

    python validation/test_duration.py            # run (~10-15 min)
    python validation/test_duration.py --analyze
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "simulations"))
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from flopy.utils.sfroutputfile import SfrFile

from parallel_modflow import DATA_INPUT, DEFAULT_EXE, run_single_simulation
from sfr_series_utils import recompute_times

INPUT_CSV = DATA_INPUT / "inputs_main.csv"

OUT = ROOT / "data" / "output" / "validation"
RUNS = str(OUT / "runs")
OUT.mkdir(parents=True, exist_ok=True)
IDS = [49, 331, 381, 81, 459, 512, 130, 658, 713, 41, 96, 612]
CFG = {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0, 'save_cbc': False, 'duration_days': 730}
SUBDIR = "duration_test"


def exchange_series(sim_id, p):
    T = recompute_times(p, CFG) / 86400
    df = SfrFile(os.path.join(RUNS, SUBDIR, f"run_{sim_id:03d}", f"water_sim_bagre_{sim_id}.sfr.out")).df
    q = df["Qaquifer"].to_numpy().reshape(-1, df["reach"].max()).sum(axis=1)[1:]
    return T, q


def settle_time(T, q, frac):
    """First time after which |q| stays below frac * max(q) until the end (NaN if never)."""
    above = np.flatnonzero(np.abs(q) > frac * q.max())
    return np.nan if above[-1] == len(q) - 1 else T[above[-1] + 1]


def analyze():
    dt = pd.read_csv(INPUT_CSV).set_index("IDMERGE", drop=False)
    res = pd.read_csv(OUT / "test_duration_results.csv").set_index("simulation_id")
    rows = []
    for sim_id in IDS:
        p = dt.loc[sim_id]
        T, q = exchange_series(sim_id, p)
        dtau = np.diff(np.r_[0, T]) * 86400
        ret = np.cumsum(np.where(q < 0, -q, 0) * dtau)          # returned volume (aquifer -> river)
        row = {"ID": sim_id, "da_km2": p["da_km2"], "t_d_d": p["t_d_days"],
               "tc_years": p["Sy"] * p["B"]**2 / (p["Ka"] * p["d_bkf"]) / 86400 / 365,
               "t_r_d": res.loc[sim_id, "t_r"] / 86400,
               "q_end_%peak": 100 * abs(q[-1]) / q.max()}
        for f in (0.05, 0.02, 0.01):
            row[f"t_{int(f*100)}%_d"] = settle_time(T, q, f)
        for d in (90, 120, 150, 365):
            row[f"ret@{d}d_%of730"] = 100 * ret[np.searchsorted(T, d) - 1] / ret[-1]
        # depth spike check at the end of the hydrograph (tabfile fix)
        rows.append(row)
    t = pd.DataFrame(rows).sort_values(["da_km2", "tc_years"])
    t.to_csv(OUT / "test_duration_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(2).to_string(index=False))


def run():
    dt = pd.read_csv(INPUT_CSV).set_index("IDMERGE", drop=False)
    jobs = [(i, dt.loc[i], RUNS, SUBDIR, DEFAULT_EXE, CFG) for i in IDS]
    with ProcessPoolExecutor(max_workers=len(jobs)) as ex:
        results = [f.result() for f in as_completed([ex.submit(run_single_simulation, j) for j in jobs])]
    pd.DataFrame(results).to_csv(OUT / "test_duration_results.csv", index=False)
    print(pd.DataFrame(results)[["simulation_id", "success", "nstp_total", "execution_time",
                                 "PERCENT_DISCREPANCY"]].to_string(index=False))


if __name__ == "__main__":
    analyze() if "--analyze" in sys.argv else run()
