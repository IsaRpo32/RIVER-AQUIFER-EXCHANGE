"""
Compare the old uniform Courant time step (as used for the reference Excel) against the new
adaptive schemes, on two cheap rows with the ORIGINAL Manning n.
Output folders keep the .sfr.out so the depth series can be plotted afterwards.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "simulations"))

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from parallel_modflow import DATA_INPUT, DEFAULT_EXE, run_single_simulation

INPUT_CSV = DATA_INPUT / "inputs_n_original.csv"

OUT = ROOT / "data" / "output" / "validation"
RUNS = str(OUT / "runs")
OUT.mkdir(parents=True, exist_ok=True)
IDS = [96, 612]
SCHEMES = {
    'old_courant': {'mode': 'courant', 'courant_factor': 1.0},
    'coarse': {'n_hyd': 100, 'tsmult': 1.2, 'dt_max_hours': 12.0},
    'default': {},
    'fine': {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0},
}

if __name__ == "__main__":
    dt = pd.read_csv(INPUT_CSV).set_index("IDMERGE", drop=False)
    jobs = [(scheme, (i, dt.loc[i], RUNS, f"scheme_examples/{scheme}",
                      DEFAULT_EXE, cfg))
            for scheme, cfg in SCHEMES.items() for i in IDS]
    results = []
    with ProcessPoolExecutor(max_workers=len(jobs)) as ex:
        futures = {ex.submit(run_single_simulation, p): s for s, p in jobs}
        for f in as_completed(futures):
            r = f.result(); r["scheme"] = futures[f]; results.append(r)
    res = pd.DataFrame(results).sort_values(["simulation_id", "scheme"])
    res.to_csv(OUT / "test_schemes_examples_results.csv", index=False)
    cols = ["simulation_id", "scheme", "success", "nstp_total", "dt_min_s", "execution_time",
            "Qexc_in", "Qexc_out", "T_in", "T_out", "delta_tp", "H_ratio", "delta_H"]
    print(res[[c for c in cols if c in res]].to_string(index=False))
