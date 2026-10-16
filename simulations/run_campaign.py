"""
Run a campaign (or one part of it) on any machine, resumable, and pack the results to send back.

    python simulations/run_campaign.py inputs_n_corregido_mitad_A.csv     # one half of the main ensemble
    python simulations/run_campaign.py inputs_n_corregido_mitad_B.csv 8   # the other half, 8 parallel simulations
    python simulations/run_campaign.py inputs_high_connectivity.csv 8     # the high-connectivity ensemble

The CSV is read from ../data/input/. Rows already finished (data/output/<name>/rows_done.jsonl) are
skipped, so the run can be stopped and started again. At the end it writes
data/output/<name>/results.csv with every finished row (the file the analysis uses), and in data/output/:
    enviar_resultados_<name>.zip   small: results CSV and rows_done.jsonl (all metrics and budget terms)
    sfr_out_<name>.zip             large: every .sfr.out and .list (only needed to recompute the series)
"""
import json
import multiprocessing as mp
import sys
import zipfile
from pathlib import Path

import pandas as pd

from parallel_modflow import DATA_INPUT, DATA_OUTPUT, DEFAULT_EXE, consolidate_results, run_parallel_simulations

INPUT_NAME = sys.argv[1] if len(sys.argv) > 1 else "inputs_n_corregido.csv"
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else max(1, mp.cpu_count() - 2)
NAME = Path(INPUT_NAME).stem
NAME = NAME[len("inputs_"):] if NAME.startswith("inputs_") else NAME
# Time-step scheme validated in validation/test_time_steps.py and test_mass_balance.py ('fine'); no .cbc / .hds
TIME_CFG = {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0, 'save_cbc': False}


def finished_ids(out_dir: Path) -> set:
    done = set()
    jl = out_dir / "rows_done.jsonl"
    if jl.exists():
        for line in jl.read_text().splitlines():
            r = json.loads(line)
            if str(r.get("success")) == "True":
                done.add(int(r["IDMERGE"]))
    return done


def write_results(out_dir: Path, rows: pd.DataFrame) -> Path:
    """All finished rows in one CSV (same columns as complete_sample_points_*.csv), from rows_done.jsonl.
    Needed when the run was stopped and resumed: each call writes a complete_sample_points_*.csv with
    only the rows it ran."""
    last = {}
    for line in (out_dir / "rows_done.jsonl").read_text().splitlines():
        r = json.loads(line)
        if str(r.get("success")) == "True":
            last[int(r["IDMERGE"])] = {k: v for k, v in r.items() if k != "IDMERGE"}
    out = out_dir / "results.csv"
    consolidate_results(list(last.values()), rows.reset_index(drop=True)).to_csv(out, index=False)
    return out


def pack(out_dir: Path):
    small = DATA_OUTPUT / f"enviar_resultados_{NAME}.zip"
    with zipfile.ZipFile(small, "w", zipfile.ZIP_DEFLATED) as z:
        for f in out_dir.glob("*.csv"):
            z.write(f, f.relative_to(DATA_OUTPUT))
        z.write(out_dir / "rows_done.jsonl", (out_dir / "rows_done.jsonl").relative_to(DATA_OUTPUT))
    big = DATA_OUTPUT / f"sfr_out_{NAME}.zip"
    with zipfile.ZipFile(big, "w", zipfile.ZIP_DEFLATED) as z:
        for f in list(out_dir.glob("run_*/*.sfr.out")) + list(out_dir.glob("run_*/*.list")):
            z.write(f, f.relative_to(DATA_OUTPUT))
    print(f"\nDone. Send: {small} ({small.stat().st_size/1e6:.1f} MB)"
          f"\nOptional (full series): {big} ({big.stat().st_size/1e9:.2f} GB)")


if __name__ == "__main__":
    out_dir = DATA_OUTPUT / NAME
    rows = pd.read_csv(DATA_INPUT / INPUT_NAME)
    done = finished_ids(out_dir)
    todo = rows[~rows["IDMERGE"].isin(done)]
    print(f"{NAME}: {len(rows)} rows, {len(done)} already finished, {len(todo)} to run, {WORKERS} in parallel")

    if len(todo):
        out_dir.mkdir(parents=True, exist_ok=True)
        pending_csv = out_dir / "_pending.csv"
        todo.to_csv(pending_csv, index=False)
        run_parallel_simulations(max_workers=WORKERS, input_csv_path=str(pending_csv),
                                 base_results_dir=str(DATA_OUTPUT), simulation_subdir=NAME,
                                 modflow_exe_path=DEFAULT_EXE, time_cfg=TIME_CFG)

    done = finished_ids(out_dir)
    missing = sorted(set(rows["IDMERGE"]) - done)
    if missing:
        print(f"\nWARNING: {len(missing)} rows failed: {missing}\nRun again to retry them.")
    print(f"All finished rows: {write_results(out_dir, rows)}")
    pack(out_dir)
