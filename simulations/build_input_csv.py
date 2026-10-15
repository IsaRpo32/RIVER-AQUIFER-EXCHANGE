"""
Build the input CSVs for parallel_modflow.py from the master Excel
(Appendix_2_sim_Master_sample_points_merge_20251013.xlsx, sheet 'sim_outputs').

- Keeps only the input columns the model reads (+ IDMERGE / date for traceability).
- Fills the missing 'nrow' (campaigns 10092025 and 30082025) with L / dx, dx = C*dr/(Kr*w_bkf).
- Sets repetir = 1 so every row is run.
- Writes two versions: original Manning n and n/10 (corrected range 0.025-0.045).
  The master Excel has n = 0.25-0.45, ten times a realistic channel value; the corrected
  file keeps the same draws divided by 10 and the original value in 'n_original'.

Usage:
    python build_input_csv.py path/to/Appendix_2_sim_Master_sample_points_merge_20251013.xlsx
"""
import sys
from pathlib import Path

import pandas as pd

DATA_INPUT = Path(__file__).resolve().parents[1] / "data" / "input"
INPUT_COLUMNS = ['IDMERGE', 'date', 'da_km2', 'Sy', 'Ka', 'L', 'B', 'D', 'dr', 'Kr', 'S', 'n',
                 't_d_days', 'H_0', 'eta_0', 'H_p', 'w_bkf', 'd_bkf', 'C', 'nrow']


def build_inputs(excel_path: str) -> pd.DataFrame:
    dt = pd.read_excel(excel_path, sheet_name="sim_outputs")[INPUT_COLUMNS].copy()
    nrow_from_C = dt["L"] / (dt["C"] * dt["dr"] / (dt["Kr"] * dt["w_bkf"]))
    dt["nrow"] = dt["nrow"].fillna(nrow_from_C).round().astype(int)
    dt["repetir"] = 1
    return dt


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    dt = build_inputs(sys.argv[1])
    DATA_INPUT.mkdir(parents=True, exist_ok=True)
    dt.to_csv(DATA_INPUT / "inputs_n_original.csv", index=False)

    dt_corr = dt.copy()
    dt_corr["n_original"] = dt_corr["n"]
    dt_corr["n"] = dt_corr["n"] / 10
    dt_corr.to_csv(DATA_INPUT / "inputs_n_corregido.csv", index=False)

    print(f"{len(dt)} rows written to {DATA_INPUT}/inputs_n_original.csv and inputs_n_corregido.csv")
    print(f"n corrected range: {dt_corr['n'].min():.4f} - {dt_corr['n'].max():.4f}")
