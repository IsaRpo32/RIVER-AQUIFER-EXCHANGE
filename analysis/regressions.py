"""
Log-log regressions on the 1,058 simulations (data/output/metrics.csv).

Predictive models: log10(Y*) = a0 + a1 log10(pi1) + ... + a4 log10(pi4)
    -> coefficients, 95 % CI (OLS), p-values, R2adj, 10-fold cross-validated R2.
Hierarchy of controls: model with pi1-pi3 only; fully standardized beta (x and y), analytical 95 % CI (OLS)
    and percentile bootstrap 95 % CI (3,000 replicates). pi4 is left out because it describes the forcing
    and is collinear with pi3 through t_c (r = 0.96).

Cross-validation uses the same partition as sklearn's KFold(10, shuffle=True, random_state=42).
No scipy: p-values use the normal approximation (more than 1,000 degrees of freedom).

    python analysis/regressions.py
        -> data/output/tables/predictive_models_pi1_pi4.csv, hierarchy_beta_pi1_pi3.csv
"""
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "data" / "output" / "tables"; TABLES.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(ROOT / "data" / "output" / "metrics.csv")
MET = ["delta_tp*", "delta_H*", "t_r*", "delta_leakage*", "delta_storage*"]
PI4, PI3 = ["pi1", "pi2", "pi3", "pi4"], ["pi1", "pi2", "pi3"]
Z = 1.96
N_BOOT, SEED = 3000, 12345


def ols(X, y):
    X1 = np.column_stack([np.ones(len(y)), X])
    coef, *_ = np.linalg.lstsq(X1, y, rcond=None)
    res = y - X1 @ coef
    n, k = X1.shape
    s2 = res @ res / (n - k)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X1.T @ X1)))
    r2 = 1 - res.var() / y.var()
    return coef, se, r2, 1 - (1 - r2) * (n - 1) / (n - k), float(np.sqrt(np.mean(res**2)))


def r2_cv(X, y, k=10, seed=42):
    """Same partition as sklearn KFold(10, shuffle=True, random_state=42)."""
    from numpy.random import RandomState
    idx = np.arange(len(y)); RandomState(seed).shuffle(idx)
    sizes = np.full(k, len(y) // k); sizes[: len(y) % k] += 1
    yt, yp, start = [], [], 0
    for s in sizes:
        te = idx[start:start + s]; tr = np.setdiff1d(idx, te); start += s
        c = ols(X[tr], y[tr])[0]
        yt.append(y[te]); yp.append(np.column_stack([np.ones(len(te)), X[te]]) @ c)
    yt, yp = np.concatenate(yt), np.concatenate(yp)
    return 1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2)


def pval(t):
    return math.erfc(abs(t) / math.sqrt(2))


rows_pred, rows_beta = [], []
rng = np.random.default_rng(SEED)
for m in MET:
    y = np.log10(d[m].to_numpy())
    # --- predictive model with pi1-pi4
    X = np.log10(d[PI4].to_numpy())
    c, se, r2, r2a, rmse = ols(X, y)
    cv = r2_cv(X, y)
    for i, term in enumerate(["a0"] + PI4):
        rows_pred.append(dict(metric=m, term=term, coef=c[i], ci_low=c[i] - Z * se[i], ci_high=c[i] + Z * se[i],
                              p_value=pval(c[i] / se[i]), n=len(y), R2=r2, R2adj=r2a, R2cv=cv, rmse_log10=rmse))
    # --- hierarchy with pi1-pi3 (fully standardized beta)
    X3 = np.log10(d[PI3].to_numpy())
    Zx = (X3 - X3.mean(0)) / X3.std(0, ddof=1); zy = (y - y.mean()) / y.std(ddof=1)
    b, seb, r2b, r2ab, _ = ols(Zx, zy)
    boot = np.empty((N_BOOT, 3))
    for r in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        Xi, yi = X3[i], y[i]
        boot[r] = ols((Xi - Xi.mean(0)) / Xi.std(0, ddof=1), (yi - yi.mean()) / yi.std(ddof=1))[0][1:]
    vif = []
    for j in range(3):
        rj = ols(np.delete(X3, j, axis=1), X3[:, j])[2]
        vif.append(1 / (1 - rj))
    order = np.argsort(-np.abs(b[1:]))
    for j, p in enumerate(PI3):
        rows_beta.append(dict(metric=m, pi=p, beta=b[1 + j], ci_ols_low=b[1 + j] - Z * seb[1 + j],
                              ci_ols_high=b[1 + j] + Z * seb[1 + j], ci_boot_low=np.percentile(boot[:, j], 2.5),
                              ci_boot_high=np.percentile(boot[:, j], 97.5), rank=int(np.where(order == j)[0][0]) + 1,
                              VIF=vif[j], R2adj_pi123=r2ab, n=len(y)))

tp, tb = pd.DataFrame(rows_pred), pd.DataFrame(rows_beta)
tp.to_csv(TABLES / "predictive_models_pi1_pi4.csv", index=False)
tb.to_csv(TABLES / "hierarchy_beta_pi1_pi3.csv", index=False)

print("PREDICTIVE MODELS (coefficients, pi1-pi4)")
w = tp.pivot(index="metric", columns="term", values="coef").reindex(MET)[["a0"] + PI4]
w["R2adj"] = tp.groupby("metric").R2adj.first(); w["R2cv"] = tp.groupby("metric").R2cv.first()
print(w.round(3).to_string())
print("\nHIERARCHY (beta, pi1-pi3)  beta [bootstrap CI] (rank)")
for m in MET:
    s = tb[tb.metric == m]
    print(f"  {m:15s} R2adj={s.R2adj_pi123.iloc[0]:.3f} | " + " | ".join(
        f"{r.pi} {r.beta:+.3f} [{r.ci_boot_low:+.3f},{r.ci_boot_high:+.3f}] (#{r.rank})" for r in s.itertuples()))
print("\nmean CI width: OLS %.3f  bootstrap %.3f" % ((tb.ci_ols_high - tb.ci_ols_low).mean(),
                                                    (tb.ci_boot_high - tb.ci_boot_low).mean()))
