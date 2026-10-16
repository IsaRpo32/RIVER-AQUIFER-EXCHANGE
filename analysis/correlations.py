"""
Pearson and Spearman correlations between the metrics and pi1-pi4 (log10 values), and other numbers quoted
in the Results section. Correlations with pi1 use the main ensemble only (n = 758), as in the figures.

    python analysis/correlations.py      ->  data/output/tables/correlations_pi_metrics.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "data" / "output" / "tables"; TABLES.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(ROOT / "data" / "output" / "metrics.csv")
PI = ["pi1", "pi2", "pi3", "pi4"]
MET = ["delta_tp*", "delta_H*", "t_r*", "delta_leakage*", "delta_storage*"]
rk = lambda s: s.rank()
rows = []
for m in MET:
    for p in PI:
        g = d[d.group == "main"] if p == "pi1" else d
        x, y = np.log10(g[p]), np.log10(g[m])
        rows.append(dict(metric=m, pi=p, n=len(g), pearson=x.corr(y), spearman=rk(x).corr(rk(y))))
c = pd.DataFrame(rows)
c.to_csv(TABLES / "correlations_pi_metrics.csv", index=False)
print(c.pivot(index="metric", columns="pi", values="spearman").reindex(MET).round(2).to_string())
print(c.pivot(index="metric", columns="pi", values="pearson").reindex(MET).round(2).to_string())


def r2a(cols, m):
    X = np.column_stack([np.ones(len(d))] + [np.log10(d[c]) for c in cols]); y = np.log10(d[m])
    res = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]; r2 = 1 - res.var() / y.var()
    return 1 - (1 - r2) * (len(y) - 1) / (len(y) - X.shape[1])


print("\nR2adj pi1-pi3 | pi1-pi4 | gain")
for m in MET:
    a, b = r2a(PI[:3], m), r2a(PI, m)
    print(f"  {m:15s} {a:.3f} {b:.3f} {b - a:+.3f}")
L = np.log10(d[PI]); print("\ncorrelation of log pi:\n", L.corr().round(2).to_string())
print("\nranges (1058):", {p: (f"{d[p].min():.2g}", f"{d[p].max():.2g}") for p in PI + ["tc"]})
print("median tc by class:", d.groupby("area_class").tc.median().map(lambda v: f"{v:.2g}").to_dict())
print("fraction dH* >= 0.10:", (d["delta_H*"] >= 0.1).mean().round(3), "| >= 0.05:", (d["delta_H*"] >= 0.05).mean().round(3),
      "| median dH*:", d["delta_H*"].median().round(4))
print("median dtp (h):", (d.delta_tp / 3600).median().round(2), "| median t_r (d):", (d.t_r / 86400).median().round(2))
