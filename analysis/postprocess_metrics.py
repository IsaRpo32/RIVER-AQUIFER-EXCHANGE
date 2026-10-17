"""
Post-processing of the 1,058 simulations: characteristic scales, pi-groups and dimensionless metrics.

Reads data/output/simulation_results.csv (raw results of every run) and writes data/output/model_results.csv,
the results table of the paper (columns described in model_results_dictionary.csv), which every other script
in analysis/ and figures/ reads.

    python analysis/postprocess_metrics.py

Peak delay and return time are scaled with the kinematic travel time t_L = L / c_k and not with t_c,
because t_c also enters pi1, pi3 and pi4 and induces spurious correlation (see time_scale_test.py).

The IDMERGE (config_id) values of the high-connectivity ensemble repeat those of the main ensemble (same base
configurations, different Kr and Ka): use (ensemble, config_id) as the key.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from results_table import RESULTS, SHORT

ROOT = Path(__file__).resolve().parents[1]
SIMULATIONS = ROOT / "data" / "output" / "simulation_results.csv"

m = pd.read_csv(SIMULATIONS)

# Characteristic scales
td_s = m.t_d_days * 86400
m["tc"] = m.Sy * m.B**2 / (m.Ka * m.d_bkf)                       # aquifer response time (enters pi1, pi3, pi4)
m["V_bkf"] = m.d_bkf**(2/3) * m.S**0.5 / m.n                       # Manning velocity at bankfull (R ~ d_bkf)
m["tL"] = m.L / (5/3 * m.V_bkf)                                    # kinematic travel time, t_L = L / c_k
V_ref = m.B * m.L * m.d_bkf                                        # bankfull control volume

# pi-groups
m["pi1"] = m.tc * m.Kr / m.dr                                      # vertical exchange efficiency
m["pi2"] = m.w_bkf / (m.Sy * m.B)                                  # storage partitioning
m["pi3"] = m.d_bkf**(2/3) * m.S**0.5 * m.tc / (m.n * m.B)          # kinematic contrast
m["pi4"] = m.tc / td_s                                             # event persistence (forcing)

# Dimensionless metrics
m["delta_tp*"] = m.delta_tp / m.tL                                 # peak delay
m["delta_H*"] = m.delta_H / m.d_bkf                                # stage attenuation
m["delta_leakage*"] = (m.STREAM_LEAKAGE_IN - m.STREAM_LEAKAGE_OUT) / V_ref   # net river-to-aquifer exchange
m["t_r*"] = m.t_r / m.tL                                           # return time
m["delta_storage*"] = (m.STORAGE_OUT - m.STORAGE_IN) / V_ref       # net storage gain

m["Ka_Kr"] = m.Ka / m.Kr
m["area_class_km2"] = np.select([m.da_km2 < 10**3.5, m.da_km2 < 10**4.5], [1e3, 1e4], 1e5)  # nearest in log10

table = m[list(SHORT.values())].set_axis(list(SHORT.keys()), axis=1)
table.to_csv(RESULTS, index=False)

cols = ["t_c_s", "t_L_s", "pi1", "pi2", "pi3", "pi4", "delta_tp_star", "delta_H_star", "delta_leakage_star",
        "t_r_star", "delta_storage_star"]
print(f"{len(table)} rows -> {RESULTS.relative_to(ROOT)}")
print("values <= 0 or NaN:", {c: int((~(table[c] > 0)).sum()) for c in cols})
print("\nlog10 by ensemble (min / median / max):")
print(np.log10(table[cols]).groupby(table.ensemble).agg(["min", "median", "max"]).T.round(2).to_string())
print("\nrows by ensemble and drainage-area class:\n", pd.crosstab(table.area_class_km2, table.ensemble).to_string())
