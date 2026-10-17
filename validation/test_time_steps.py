"""
Time-step convergence test for simulations/parallel_modflow.py.

Runs a few representative rows (3 per basin: low / median / high tc) with several
time-step schemes and compares every metric against a very fine reference.
Choose the cheapest scheme whose errors stay below TOLERANCE.

    python validation/test_time_steps.py
"""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "simulations"))
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from parallel_modflow import DEFAULT_EXE, load_inputs, run_single_simulation

MODFLOW_EXE = DEFAULT_EXE

OUT = ROOT / "data" / "output" / "validation"
RUNS = str(OUT / "runs")
OUT.mkdir(parents=True, exist_ok=True)
BASE_DIR = RUNS
MAX_WORKERS = 9
TOLERANCE = 0.02            # 2 % relative error accepted against the reference
RUN_OLD_COURANT = False     # True also runs the old uniform Courant scheme (very slow with corrected n)

SCHEMES = {
    'ref':    {'n_hyd': 1000, 'tsmult': 1.02, 'dt_max_hours': 1.0, 'flood_margin': 2.0, 'courant_sfr': 0.5},
    'fine':   {'n_hyd': 400,  'tsmult': 1.05, 'dt_max_hours': 3.0},
    'default': {},
    'coarse': {'n_hyd': 100,  'tsmult': 1.2,  'dt_max_hours': 12.0},
}
if RUN_OLD_COURANT:
    SCHEMES['old_courant'] = {'mode': 'courant', 'courant_factor': 1.0}

METRICS = ['Qexc_in', 'Qexc_out', 'T_in', 'T_out', 'delta_tp', 'H_ratio', 'delta_H']


def pick_rows(dt: pd.DataFrame) -> pd.DataFrame:
    """3 rows per main basin size: 10th, 50th and 90th percentile of tc = Sy*B^2/(Ka*d_bkf)."""
    dt = dt.copy()
    dt["tc"] = dt["Sy"] * dt["B"]**2 / (dt["Ka"] * dt["d_bkf"])
    picked = []
    for da in [1000, 10000, 100000]:
        sub = dt[dt["da_km2"] == da].sort_values("tc")
        for q in [0.1, 0.5, 0.9]:
            picked.append(sub.iloc[int(q * (len(sub) - 1))])
    return pd.DataFrame(picked)


if __name__ == "__main__":
    rows = pick_rows(load_inputs("main"))
    print(f"Testing IDMERGE {rows['IDMERGE'].tolist()} with schemes {list(SCHEMES)}")

    jobs = []
    for scheme, cfg in SCHEMES.items():
        for _, row in rows.iterrows():
            sim_id = int(row["IDMERGE"])
            jobs.append((scheme, (sim_id, row, BASE_DIR, f"time_step_test/{scheme}", MODFLOW_EXE, cfg)))

    results = []
    start = time.time()
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(run_single_simulation, params): scheme for scheme, params in jobs}
        for future in as_completed(futures):
            res = future.result()
            res["scheme"] = futures[future]
            results.append(res)
            print(f"  {res['scheme']:12s} IDMERGE {res['simulation_id']}: "
                  f"{'OK' if res['success'] else 'FAILED ' + res.get('error', '')} "
                  f"({res['execution_time']/60:.1f} min)")
    print(f"Total wall time: {(time.time() - start)/60:.1f} min")

    res = pd.DataFrame(results)
    res = res.merge(rows[["IDMERGE", "da_km2", "tc"]], left_on="simulation_id", right_on="IDMERGE")
    ref = res[res["scheme"] == "ref"].set_index("simulation_id")
    for m in METRICS:
        res[f"err_{m}"] = [
            abs(r[m] - ref.loc[r["simulation_id"], m]) / abs(ref.loc[r["simulation_id"], m])
            if r["simulation_id"] in ref.index and r["success"] else np.nan
            for _, r in res.iterrows()
        ]
    res["max_err"] = res[[f"err_{m}" for m in METRICS]].max(axis=1)
    res.to_csv(OUT / "test_time_steps_results.csv", index=False)

    summary = res.groupby(["scheme", "da_km2"]).agg(
        runs_ok=("success", "sum"),
        nstp=("nstp_total", "median"),
        minutes=("execution_time", lambda s: s.median() / 60),
        disk_mb=("disk_mb", "median"),
        max_err=("max_err", "max"),
        return_closed=("return_closed", "all"),
        discrepancy_pct=("PERCENT_DISCREPANCY", "max"),
    )
    print("\n" + summary.round(3).to_string())

    ok = summary.reset_index()
    ok = ok[ok["scheme"] != "ref"].groupby("scheme").agg(max_err=("max_err", "max"), minutes=("minutes", "sum"))
    ok = ok[ok["max_err"] <= TOLERANCE].sort_values("minutes")
    if len(ok):
        print(f"\nCheapest scheme within {TOLERANCE:.0%}: '{ok.index[0]}' -> {SCHEMES[ok.index[0]] or 'DEFAULT_TIME_CFG'}")
    else:
        print(f"\nNo scheme within {TOLERANCE:.0%} of the reference: refine n_hyd / dt_max_hours.")
    print(f"Details per run in {OUT / 'test_time_steps_results.csv'}")
