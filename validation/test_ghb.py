"""
Sensitivity of the results to the GHB conductance (audit finding A5: C = Ka*w_bkf*delc/D).
Same rows as test_mass_balance.py, conductance x0.1, x1 and x10, 'fine' time scheme, corrected n.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "simulations"))

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from parallel_modflow import DEFAULT_EXE, load_inputs, run_single_simulation


OUT = ROOT / "data" / "output" / "validation"
RUNS = str(OUT / "runs")
OUT.mkdir(parents=True, exist_ok=True)
IDS = [41, 96, 612]
FINE = {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0, 'save_cbc': False}
FACTORS = [0.1, 1.0, 10.0]
METRICS = ['Qexc_in', 'Qexc_out', 'T_in', 'T_out', 't_r', 'delta_tp', 'delta_H',
           'STORAGE_IN', 'STORAGE_OUT', 'HEAD_DEP_BOUNDS_IN', 'HEAD_DEP_BOUNDS_OUT',
           'STREAM_LEAKAGE_IN', 'STREAM_LEAKAGE_OUT', 'PERCENT_DISCREPANCY']

if __name__ == "__main__":
    dt = load_inputs("main").set_index("IDMERGE", drop=False)
    jobs = [(f, (i, dt.loc[i], RUNS, f"ghb_test/x{f:g}", DEFAULT_EXE,
                 {**FINE, 'ghb_cond_factor': f})) for f in FACTORS for i in IDS]
    results = []
    with ProcessPoolExecutor(max_workers=len(jobs)) as ex:
        futures = {ex.submit(run_single_simulation, p): f for f, p in jobs}
        for fut in as_completed(futures):
            r = fut.result(); r["ghb_factor"] = futures[fut]; results.append(r)
    res = pd.DataFrame(results).sort_values(["simulation_id", "ghb_factor"])
    res.to_csv(OUT / "test_ghb_results.csv", index=False)

    rows = []
    for sim_id, g in res.groupby("simulation_id"):
        base = g[g["ghb_factor"] == 1.0].iloc[0]
        for _, r in g.iterrows():
            row = {"ID": sim_id, "ghb_factor": r["ghb_factor"], "ok": r["success"]}
            for m in METRICS:
                row[m] = r[m]
                if m != 'PERCENT_DISCREPANCY':
                    row[f"chg_{m}_%"] = 100 * (r[m] - base[m]) / abs(base[m])
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "test_ghb_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(t[["ID", "ghb_factor", "ok", "PERCENT_DISCREPANCY"] + [c for c in t if c.startswith("chg_")]]
          .round(2).to_string(index=False))
