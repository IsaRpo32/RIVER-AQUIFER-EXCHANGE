"""
Post-processing of the 1,058 simulations: characteristic scales, pi-groups and dimensionless metrics.

Reads the three result files of the campaign (data/output/campaign/) and writes data/output/metrics.csv,
which every other script in analysis/ and figures/ reads.

    python analysis/postprocess_metrics.py

Peak delay and return time are scaled with the kinematic travel time t_L = L / c_k and not with t_c,
because t_c also enters pi1, pi3 and pi4 and induces spurious correlation (see time_scale_test.py).
The t_c-scaled versions are kept as delta_tp*_tc and t_r*_tc for comparison.

The IDMERGE values of the high-connectivity ensemble repeat those of the main ensemble (same base
configurations, different Kr and Ka): use (ensemble, IDMERGE) as the key.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "data" / "output" / "campaign"
OUT = ROOT / "data" / "output" / "metrics.csv"
SOURCES = {  # file -> ensemble
    "main_half_A.csv": "main (half A)",
    "main_half_B.csv": "main (half B)",
    "high_connectivity.csv": "high connectivity",
}

parts = []
for f, ensemble in SOURCES.items():
    d = pd.read_csv(CAMPAIGN / f)
    d.insert(1, "ensemble", ensemble)
    d.insert(2, "source_file", f)
    parts.append(d)
m = pd.concat(parts, ignore_index=True)

# Characteristic scales
m["td_s"] = m.t_d_days * 86400
m["tc"] = m.Sy * m.B**2 / (m.Ka * m.d_bkf)                       # aquifer response time (enters pi1, pi3, pi4)
m["V_bkf"] = m.d_bkf**(2/3) * m.S**0.5 / m.n                       # Manning velocity at bankfull (R ~ d_bkf)
m["tL"] = m.L / (5/3 * m.V_bkf)                                    # kinematic travel time, t_L = L / c_k
V_ref = m.B * m.L * m.d_bkf                                        # bankfull control volume

# pi-groups
m["pi1"] = m.tc * m.Kr / m.dr                                      # vertical exchange efficiency
m["pi2"] = m.w_bkf / (m.Sy * m.B)                                  # storage partitioning
m["pi3"] = m.d_bkf**(2/3) * m.S**0.5 * m.tc / (m.n * m.B)          # kinematic contrast (corrected n)
m["pi4"] = m.tc / m.td_s                                           # event persistence (forcing)
m["pi3_n_original"] = m.pi3 * m.n / m.n_original                   # pi3 with the original n, for comparison only

# Dimensionless metrics
m["delta_tp*"] = m.delta_tp / m.tL                                 # peak delay
m["delta_H*"] = m.delta_H / m.d_bkf                                # stage attenuation
m["delta_leakage*"] = (m.STREAM_LEAKAGE_IN - m.STREAM_LEAKAGE_OUT) / V_ref   # net river-to-aquifer exchange
m["t_r*"] = m.t_r / m.tL                                           # return time
m["delta_storage*"] = (m.STORAGE_OUT - m.STORAGE_IN) / V_ref       # net storage gain

# Auxiliary columns
m["delta_tp*_tc"] = m.delta_tp / m.tc                              # earlier normalization by t_c, comparison only
m["t_r*_tc"] = m.t_r / m.tc
m["delta_H*_reach"] = m.delta_H_reach / m.d_bkf                    # earlier reach-based definition, comparison only
m["Ka_Kr"] = m.Ka / m.Kr
m["L_B"] = m.L / m.B
m["group"] = np.where(m.ensemble == "high connectivity", "high_connectivity", "main")
m["area_class"] = np.select([m.da_km2 < 10**3.5, m.da_km2 < 10**4.5], ["1e3", "1e4"], "1e5")  # nearest in log10

m.to_csv(OUT, index=False)

cols = ["tc", "tL", "pi1", "pi2", "pi3", "pi4", "delta_tp*", "delta_H*", "delta_leakage*", "t_r*", "delta_storage*"]
print(f"{len(m)} rows -> {OUT.relative_to(ROOT)}")
print("values <= 0 or NaN:", {c: int((~(m[c] > 0)).sum()) for c in cols})
print("\nlog10 by group (min / median / max):")
print(np.log10(m[cols]).groupby(m.group).agg(["min", "median", "max"]).T.round(2).to_string())
print("\nrows by group and drainage-area class:\n", pd.crosstab(m.area_class, m.group).to_string())
