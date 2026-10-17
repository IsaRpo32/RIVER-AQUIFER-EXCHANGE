"""
The Nechí corridor in pi-space: pi-groups, metrics (peak delay and return time scaled with t_L), percentiles
among the 1,058 simulations, predictions of the pi1-pi4 models and comparison of its parameters with the
sampled ranges.

Nechí parameters (manuscript v14, Section 4 and Table 10): B = 4 km (total width), d_bkf = 3.66 m,
S = 1.5e-4, n = 0.04, Ka = 1.2e-3 m/s, Kr = 1e-6 m/s, Sy = 0.25, t_d = 7 d, L = 50 km. W_bkf = 200 m and
d_r = 1.25 m follow from pi2 = 0.20 and pi1 = 728.6 (Table 10). The dimensional metrics come from Table 9 of
v14 (computed with t_c = 9.107e8 s): delta_tp = 11.64 h, t_r = 7.21 d.

Some Nechí values lie outside the sampled ranges (d_bkf, and its distance to the nearest simulation in
pi-space is larger than the 95th percentile of the ensemble): its predictions are extrapolations.

    python analysis/nechi_position.py      (after regressions.py)  ->  data/output/tables/nechi_position.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

from results_table import load_results

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "data" / "output" / "tables"; TABLES.mkdir(parents=True, exist_ok=True)
d = load_results()
tp = pd.read_csv(TABLES / "predictive_models_pi1_pi4.csv")
N = dict(B=4000.0, d_bkf=3.66, S=1.5e-4, n=0.04, Ka=1.2e-3, Kr=1e-6, Sy=0.25, td_s=7 * 86400.0, L=50000.0,
         w_bkf=200.0, dr=1.25)
N["tc"] = N["Sy"] * N["B"] ** 2 / (N["Ka"] * N["d_bkf"])
N["V_bkf"] = N["d_bkf"] ** (2 / 3) * N["S"] ** 0.5 / N["n"]
N["tL"] = N["L"] / (5 / 3 * N["V_bkf"])
N["pi1"] = N["tc"] * N["Kr"] / N["dr"]
N["pi2"] = N["w_bkf"] / (N["Sy"] * N["B"])
N["pi3"] = N["V_bkf"] * N["tc"] / N["B"]
N["pi4"] = N["tc"] / N["td_s"]
N["delta_tp*"] = 4.6e-5 * N["tc"] / N["tL"]
N["t_r*"] = 6.84e-4 * N["tc"] / N["tL"]
N.update({"delta_H*": 1.72e-2, "delta_leakage*": 2.14e-3, "delta_storage*": 1.79e-3})
MET = ["delta_tp*", "delta_H*", "t_r*", "delta_leakage*", "delta_storage*"]

rows = []
for v in ["pi1", "pi2", "pi3", "pi4"] + MET:
    ref = d          # all 1,058 simulations (the restriction to the 758 of the main ensemble with pi1 is for figures only)
    pct = 100 * (ref[v] < N[v]).mean()
    row = dict(variable=v, nechi=N[v], percentile=pct, min=ref[v].min(), median=ref[v].median(), max=ref[v].max())
    if v in MET:
        c = tp[tp.metric == v].set_index("term").coef
        pred = 10 ** (c["a0"] + sum(c[p] * np.log10(N[p]) for p in ["pi1", "pi2", "pi3", "pi4"]))
        rmse = tp[tp.metric == v].rmse_log10.iloc[0]
        row.update(prediction=pred, error_log10=np.log10(N[v] / pred), rmse_log10=rmse)
    rows.append(row)
t = pd.DataFrame(rows)
t.to_csv(TABLES / "nechi_position.csv", index=False)
pd.set_option("display.width", 200)
print(t.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print(f"\ntc {N['tc']:.4g} s | tL {N['tL']:.4g} s ({N['tL']/3600:.2f} h) | V_bkf {N['V_bkf']:.3f} m/s")
print("Parameters against the sampled ranges:")
for k, v in [("S", N["S"]), ("Ka/Kr", N["Ka"] / N["Kr"]), ("L/B", N["L"] / N["B"]), ("V_bkf", N["V_bkf"]),
             ("tL", N["tL"]), ("d_bkf", N["d_bkf"]), ("Sy", N["Sy"])]:
    col = {"Ka/Kr": d.Ka / d.Kr, "L/B": d.L / d.B}.get(k, d[k] if k in d else None)
    print(f"  {k:6s} Nechí {v:.3g} | ensemble {col.min():.3g}–{col.max():.3g}")
# nearest neighbour in log pi1-pi4 space, normalized by the standard deviation of each axis
X = np.log10(d[["pi1", "pi2", "pi3", "pi4"]].to_numpy()); s = X.std(0)
x0 = np.log10([N[p] for p in ["pi1", "pi2", "pi3", "pi4"]])
dist = np.sqrt((((X - x0) / s) ** 2).sum(1))
nn = np.sort(np.sqrt((((X[:, None, :] - X[None, :, :]) / s) ** 2).sum(2)) + np.eye(len(X)) * 1e9, axis=1)[:, 0]
print(f"distance to the nearest simulation: Nechí {dist.min():.2f} | ensemble median {np.median(nn):.2f}, "
      f"p95 {np.percentile(nn, 95):.2f}")
