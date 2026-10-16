"""
Same shuffling test as time_scale_test.py for delta_H* (scale d_bkf) and delta_leakage*, delta_storage*
(scale B L d_bkf). These metrics do not show spurious correlation from their scale.

    python analysis/other_scales_test.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(2)
d = pd.read_csv(ROOT / "data" / "output" / "metrics.csv")
PI = ["pi1", "pi2", "pi3", "pi4"]
L = {p: np.log10(d[p]).to_numpy() for p in PI}
rk = lambda a: np.argsort(np.argsort(a))
sp = lambda a, b: np.corrcoef(rk(a), rk(b))[0, 1]


def r2(y):
    X = np.column_stack([np.ones(len(y))] + [L[p] for p in PI[:3]])
    res = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return 1 - res.var() / y.var()


V = d.B * d.L * d.d_bkf
cases = {"delta_H*  (d_bkf)": (d.delta_H, d.d_bkf),
         "delta_leakage* (BLd)": (d.STREAM_LEAKAGE_IN - d.STREAM_LEAKAGE_OUT, V),
         "delta_storage* (BLd)": (d.STORAGE_OUT - d.STORAGE_IN, V)}
for n, (y, s) in cases.items():
    ly, ls = np.log10(y.to_numpy()), np.log10(s.to_numpy())
    obs = ly - ls
    perm = [ly[rng.permutation(len(ly))] - ls for _ in range(200)]
    print(f"{n:22s} R2(pi1-3) obs {r2(obs):.2f} spurious {np.mean([r2(p) for p in perm]):.2f} | rho pi3 obs "
          f"{sp(obs, L['pi3']):+.2f} spurious {np.mean([sp(p, L['pi3']) for p in perm]):+.2f} | var log scale "
          f"{ls.var():.2f} var log Y {ly.var():.2f}")
