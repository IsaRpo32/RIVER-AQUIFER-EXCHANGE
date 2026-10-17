"""
Rebuild the inputs of the main ensemble (758 configurations) from the master Excel
(Appendix_2_sim_Master_sample_points_merge_20251013.xlsx, sheet 'sim_outputs'; not in the repository).

- Keeps only the input columns the model reads (+ IDMERGE / date for traceability).
- Fills the missing 'nrow' (campaigns 10092025 and 30082025) with L / dx, dx = C*dr/(Kr*w_bkf).
- Sets run = 1 so every row is run.
- Divides Manning n by 10 (corrected range 0.025-0.045) and keeps the original value in 'n_original':
  the master Excel has n = 0.25-0.45, ten times a realistic channel value.

Writes data/input/main_inputs_rebuilt.csv and compares it with the 'main' rows of
data/input/model_inputs.csv, the inputs that were simulated.

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
    dt["run"] = 1
    return dt


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    dt = build_inputs(sys.argv[1])
    dt["n_original"] = dt["n"]
    dt["n"] = dt["n"] / 10
    dt.insert(1, "ensemble", "main")
    out = DATA_INPUT / "main_inputs_rebuilt.csv"
    dt.to_csv(out, index=False)
    print(f"{len(dt)} rows written to {out}")
    print(f"n corrected range: {dt['n'].min():.4f} - {dt['n'].max():.4f}")

    # compare with the main-ensemble rows that were simulated
    sim = pd.read_csv(DATA_INPUT / "model_inputs.csv")
    sim = sim[sim.ensemble == "main"].set_index("IDMERGE")
    new = dt.set_index("IDMERGE").reindex(sim.index)
    num = [c for c in new.select_dtypes("number").columns if c in sim]
    rel = ((new[num] - sim[num]).abs() / sim[num].abs().replace(0, 1)).max()
    print(f"Against data/input/model_inputs.csv: same rows {set(dt.IDMERGE) == set(sim.index)}; "
          f"max relative difference {rel.max():.1e} ({rel.idxmax()})")
