"""
Build the input CSV of the high-connectivity ensemble (300 configurations) from the main ensemble.

The main ensemble (data/input/inputs_n_corregido.csv, 758 rows) has low streambed conductances
(C_riv = 1e-4 to 0.2 m2/s). The high-connectivity ensemble keeps the geometry, Sy, forcing and
Manning n of those rows and reassigns Kr and Ka for coarse-sand and gravel beds:

  1. r = Ka/Kr log-uniform in [1, 10]; Kr log-uniform in [1e-4, 1e-2/r] m/s; Ka = r*Kr (so Ka <= 1e-2 m/s);
     C = Kr*w_bkf*(L/nrow)/dr recomputed with the same cell length.
  2. Keep the rows with Kr <= 1e-3 m/s (sand and fine gravel; no clean gravel) and draw 100 per
     drainage-area class (1e3, 1e4, 1e5 km2). Result sorted by C, descending.

Both steps use the seed 20261012. The simulated file was built from the main ensemble as saved from
the master Excel: rows sorted by Kr (ties by IDMERGE), with Kr, Ka and S rounded to 3 significant
figures. This script does the same, so the simulated configurations carry S with 3 significant figures
(relative change <= 3e-4). The original Kr and Ka are kept in Kr_old and Ka_old.

The file that was simulated is data/input/inputs_high_connectivity.csv. This script writes
inputs_high_connectivity_rebuilt.csv next to it and compares both: same 300 rows in the same order,
identical Kr, Ka, S and Ka_Kr, and the other model inputs within 2e-8 (digits kept by Excel). Only the
informative column Ka_old differs more, for the same reason.

The IDMERGE values are those of the main ensemble rows the configurations come from: they are NOT unique
across the two ensembles. Use (ensemble, IDMERGE) as the key when combining results.

Usage:
    python simulations/build_high_connectivity_inputs.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

DATA_INPUT = Path(__file__).resolve().parents[1] / "data" / "input"
SEED = 20261012
KR_MIN, KA_MAX, R_MIN, R_MAX = 1e-4, 1e-2, 1.0, 10.0
KR_MAX_SAMPLE = 1e-3            # m/s (~86 m/d)
PER_CLASS = 100
AREAS = [1000.0, 10000.0, 100000.0]


def round_sig(x: pd.Series, sig: int = 3) -> pd.Series:
    return x.map(lambda v: float(f"{v:.{sig}g}"))


def reassign_conductivities(d: pd.DataFrame) -> pd.DataFrame:
    """Step 1: new Kr and Ka for every row (coarse-sand and gravel beds)."""
    rng = np.random.default_rng(SEED)
    n = len(d)
    r = 10**rng.uniform(np.log10(R_MIN), np.log10(R_MAX), n)
    kr = 10**rng.uniform(np.log10(KR_MIN), np.log10(KA_MAX / r), n)
    out = d.copy()
    out["Kr_old"], out["Ka_old"] = d["Kr"], d["Ka"]
    out["Ka"], out["Kr"], out["Ka_Kr"] = r * kr, kr, r
    out["C"] = out.Kr * out.w_bkf * (out.L / out.nrow) / out.dr
    return out


def sample_per_class(d: pd.DataFrame) -> pd.DataFrame:
    """Step 2: 100 rows per drainage-area class with Kr <= KR_MAX_SAMPLE, sorted by C (descending)."""
    pool = d[d.Kr <= KR_MAX_SAMPLE]
    m = pd.concat([pool[pool.da_km2 == a].sample(n=PER_CLASS, random_state=SEED) for a in AREAS])
    return m.sort_values("C", ascending=False).reset_index(drop=True)


def build(main_csv: Path = DATA_INPUT / "inputs_n_corregido.csv") -> pd.DataFrame:
    d = pd.read_csv(main_csv)
    d = d.sort_values(["Kr", "IDMERGE"], kind="mergesort").reset_index(drop=True)
    for col in ("Kr", "Ka", "S"):
        d[col] = round_sig(d[col])
    m = sample_per_class(reassign_conductivities(d))

    r = m.Ka / m.Kr
    assert len(m) == 3 * PER_CLASS and m.IDMERGE.is_unique
    assert (r >= R_MIN).all() and (r <= R_MAX).all() and (m.Ka <= KA_MAX).all()
    assert (m.Kr >= KR_MIN).all() and (m.Kr <= KR_MAX_SAMPLE).all()
    assert ((m.n / m.n_original - 0.1).abs() < 1e-6).all()
    return m


if __name__ == "__main__":
    m = build()
    out = DATA_INPUT / "inputs_high_connectivity_rebuilt.csv"
    m.to_csv(out, index=False)
    print(f"{len(m)} rows written to {out}")
    print(m.groupby("da_km2").size().to_string())
    for c in ["Kr", "Ka", "C"]:
        print(f"{c}: min {m[c].min():.3g}  median {m[c].median():.3g}  max {m[c].max():.3g}")

    simulated = pd.read_csv(DATA_INPUT / "inputs_high_connectivity.csv")
    same_rows = (m.IDMERGE.to_numpy() == simulated.IDMERGE.to_numpy()).all()
    num = simulated.select_dtypes("number").columns.drop(["Ka_old"])
    rel = ((m[num] - simulated[num]).abs() / simulated[num].abs().replace(0, 1)).max()
    print(f"\nAgainst the simulated file: same rows and order: {same_rows}; "
          f"max relative difference {rel.max():.1e} ({rel.idxmax()}); Kr, Ka, S: {rel[['Kr', 'Ka', 'S']].max():.1e}")
