"""
Does the current script reproduce the master Excel when run with the ORIGINAL Manning n?
Rows from campaign 10092025 (run with the September 2025 script: nlay=4, 90 d, uniform Courant step).
Run with the new adaptive 'fine' scheme and duration 90 d as in 16f, compared column by column with the Excel.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "simulations"))

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from parallel_modflow import DATA_INPUT, DEFAULT_EXE, run_single_simulation

# master Excel (not in the repository): python test_reproduce_excel.py path/to/Appendix_2_...xlsx
EXCEL_PATH = sys.argv[1] if len(sys.argv) > 1 else "Appendix_2_sim_Master_sample_points_merge_20251013.xlsx"

OUT = ROOT / "data" / "output" / "validation"
RUNS = str(OUT / "runs")
OUT.mkdir(parents=True, exist_ok=True)
IDS = [341, 377, 330, 458, 511, 418]
OLD_CTE = {1000: 0.55, 10000: 0.5}
CUR_CTE = {1000: 0.5, 10000: 0.35}
COMMON = {'save_cbc': False, 'duration_days': 90}
SCHEMES = {
    'fine': lambda da: {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0},
}
# Excel column -> exported column
COMPARE = {'Qexc_in': 'Qexc_in', 'Qexc_out': 'Qexc_out', 'T_retorno': 'T_out',
           'delta_tp': 'delta_tp_reach', 'H_ratio': 'H_ratio_reach', 'delta_H': 'delta_H_reach',
           'STORAGE_IN': 'STORAGE_IN', 'STORAGE_OUT': 'STORAGE_OUT',
           'HEAD_DEP_BOUNDS_IN': 'HEAD_DEP_BOUNDS_IN', 'HEAD_DEP_BOUNDS_OUT': 'HEAD_DEP_BOUNDS_OUT',
           'STREAM_LEAKAGE_IN': 'STREAM_LEAKAGE_IN', 'STREAM_LEAKAGE_OUT': 'STREAM_LEAKAGE_OUT',
           'PERCENT_DISCREPANCY': 'PERCENT_DISCREPANCY'}

if __name__ == "__main__":
    dt = pd.read_csv(DATA_INPUT / "inputs_n_original.csv").set_index("IDMERGE", drop=False)
    jobs = [(s, (i, dt.loc[i], RUNS, f"reproduce_excel/{s}", DEFAULT_EXE,
                 {**COMMON, **f(dt.loc[i, "da_km2"])})) for s, f in SCHEMES.items() for i in IDS]
    results = []
    with ProcessPoolExecutor(max_workers=len(jobs)) as ex:
        futures = {ex.submit(run_single_simulation, p): s for s, p in jobs}
        for fut in as_completed(futures):
            r = fut.result(); r["scheme"] = futures[fut]; results.append(r)
            print(f"{r['scheme']:7s} {r['simulation_id']}: {'OK' if r['success'] else r.get('error')} "
                  f"({r['execution_time']/60:.1f} min, {r.get('nstp_total')} steps)", flush=True)
    res = pd.DataFrame(results)
    res.to_csv(OUT / "test_reproduce_excel_results.csv", index=False)

    xl = pd.read_excel(EXCEL_PATH, sheet_name="sim_outputs").set_index("IDMERGE")
    rows = []
    for _, r in res.sort_values(["simulation_id", "scheme"]).iterrows():
        e = xl.loc[r["simulation_id"]]
        row = {"ID": r["simulation_id"], "da_km2": e["da_km2"], "scheme": r["scheme"], "nstp": r["nstp_total"]}
        for xc, sc in COMPARE.items():
            row[f"{xc}_excel"] = e[xc]
            row[f"{xc}_sim"] = r[sc]
            row[f"{xc}_dif%"] = 100 * (r[sc] - e[xc]) / abs(e[xc]) if xc != 'PERCENT_DISCREPANCY' else np.nan
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "test_reproduce_excel_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(t[["ID", "da_km2", "scheme", "nstp"] + [c for c in t if c.endswith("dif%") and not c.startswith("PERCENT")]]
          .round(2).to_string(index=False))
