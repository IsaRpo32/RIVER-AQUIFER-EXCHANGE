"""
Choice of the time scale of the peak delay and the return time: spurious correlation by a shared scale
(Kenney, 1982).

For each candidate scale T, Y* = Y / T (Y = delta_tp or t_r, dimensional) and compares:
  - rho(log Y*, log pi_k) observed;
  - the "spurious" rho: the same with Y shuffled at random among simulations, which measures how much
    correlation appears only because T is shared with pi_k;
  - R2 of log Y* ~ log pi1..pi3 (observed and with Y shuffled).
Scales: t_c, t_d (event duration), t_p (= t_d/2), t_L = L / c_k (kinematic travel time, c_k = 5/3 V_bkf).

    python analysis/time_scale_test.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

from results_table import load_results

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(1)
d = load_results()
V = d.d_bkf ** (2 / 3) * d.S ** 0.5 / d.n               # Manning velocity at bankfull
d["tL"] = d.L / (5 / 3 * V)
d["tp_s"] = d.td_s / 2
SCALES = {"t_c": d.tc, "t_d": d.td_s, "t_p": d.tp_s, "t_L": d.tL}
Y = {"delta_tp": d.delta_tp, "t_r": d.t_r}
PI = ["pi1", "pi2", "pi3", "pi4"]
L = {p: np.log10(d[p]).to_numpy() for p in PI}
rk = lambda a: np.argsort(np.argsort(a))
sp = lambda a, b: np.corrcoef(rk(a), rk(b))[0, 1]


def r2_123(y):
    X = np.column_stack([np.ones(len(y))] + [L[p] for p in PI[:3]])
    res = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return 1 - res.var() / y.var()


print("variance of log10 of each scale:", {k: round(float(np.log10(v).var()), 3) for k, v in SCALES.items()},
      "| var log10 delta_tp:", round(float(np.log10(d.delta_tp).var()), 3), "| var log10 t_r:",
      round(float(np.log10(d.t_r).var()), 3))
print("median delta_tp / t_L:", round(float((d.delta_tp / d.tL).median()), 3),
      "| median t_r / t_d:", round(float((d.t_r / d.td_s).median()), 3),
      "| median t_r / t_p:", round(float((d.t_r / d.tp_s).median()), 3))
rows = []
for yn, y in Y.items():
    ly = np.log10(y.to_numpy())
    for sn, s in SCALES.items():
        ls = np.log10(s.to_numpy())
        obs = ly - ls
        perm = [ly[rng.permutation(len(ly))] - ls for _ in range(1000)]
        row = {"Y": yn, "scale": sn, "R2_pi123": r2_123(obs), "R2_spurious": np.mean([r2_123(p) for p in perm[:200]])}
        for p in PI:
            row[f"rho_{p}"] = sp(obs, L[p])
            row[f"spurious_{p}"] = np.mean([sp(q, L[p]) for q in perm[:200]])
        rows.append(row)
t = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(t.round(2).to_string(index=False))
