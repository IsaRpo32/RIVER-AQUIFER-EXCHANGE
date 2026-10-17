"""
Robustness of the correlations and of the standardized coefficients to sample size.

Bootstrap with replacement (100 %) and subsampling without replacement (90, 80, 70 and 50 %), 3,000
replicates, for the pooled ensemble (ALL) and each drainage-area class (nearest class centre in log10 Ad),
plus convergence curves of Spearman's rho with sample size.
  - correlations with pi1: main ensemble only (758), as in the figures;
  - standardized beta of the pi1-pi3 model (the hierarchy leaves pi4 out);
  - vectorized with numpy (ranks by argsort; no scipy).

Takes several minutes.

    python analysis/robustness.py      ->  data/output/robustness/01_ ... 04_*.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

from results_table import load_results

ROOT = Path(__file__).resolve().parents[1]
D = load_results()
OUT = ROOT / "data" / "output" / "robustness"; OUT.mkdir(parents=True, exist_ok=True)
PI = ["pi1", "pi2", "pi3", "pi4"]; PI3 = ["pi1", "pi2", "pi3"]
MET = ["delta_tp*", "delta_H*", "t_r*", "delta_leakage*", "delta_storage*"]
N_REP, SEED = 3000, 12345
FRAC = [1.0, 0.9, 0.8, 0.7, 0.5]
SIZES_ALL = [150, 300, 450, 600, 750, 900, 1050]
SIZES_CLASS = [50, 100, 150, 200, 250, 300]
rng = np.random.default_rng(SEED)
D["area"] = np.select([D.da_km2 < 10**3.5, D.da_km2 < 10**4.5], ["1e3", "1e4"], "1e5")
GROUPS = {"ALL": D, **{f"area={c}": D[D.area == c] for c in ("1e3", "1e4", "1e5")}}


def ranks(a):
    r = np.empty_like(a, dtype=float); r[np.argsort(a, kind="mergesort")] = np.arange(len(a)); return r


def corr(x, y):
    x = x - x.mean(); y = y - y.mean()
    return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum()))


def corr_data(g, pi, met):
    g = g[g.group == "main"] if pi == "pi1" else g
    return np.log10(g[pi].to_numpy()), np.log10(g[met].to_numpy())


def sample(n, frac=None, k=None):
    if frac == 1.0:
        return rng.integers(0, n, n)
    k = k or int(np.ceil(frac * n))
    return rng.choice(n, min(k, n), replace=False)


rows_c, rows_conv = [], []
for gn, g in GROUPS.items():
    for pi in PI:
        for met in MET:
            x, y = corr_data(g, pi, met)
            n = len(x)
            full = {"pearson": corr(x, y), "spearman": corr(ranks(x), ranks(y))}
            for frac in FRAC:
                v = {"pearson": [], "spearman": []}
                for _ in range(N_REP):
                    i = sample(n, frac)
                    v["pearson"].append(corr(x[i], y[i])); v["spearman"].append(corr(ranks(x[i]), ranks(y[i])))
                for mth, vals in v.items():
                    vals = np.array(vals)
                    rows_c.append(dict(group=gn, pi=pi, metric=met, method=mth, fraction=frac, n=n, full=full[mth],
                                       median=np.median(vals), q025=np.percentile(vals, 2.5),
                                       q975=np.percentile(vals, 97.5),
                                       sign_pres=float(np.mean(np.sign(vals) == np.sign(full[mth])))))
            for k in (SIZES_ALL if gn == "ALL" else SIZES_CLASS):
                if k > n:
                    continue
                vals = np.array([corr(ranks(x[i]), ranks(y[i])) for i in (sample(n, k=k) for _ in range(N_REP))])
                rows_conv.append(dict(group=gn, pi=pi, metric=met, n_sample=k, n=n, full=full["spearman"],
                                      median=np.median(vals), abs_diff=abs(np.median(vals) - full["spearman"]),
                                      q025=np.percentile(vals, 2.5), q975=np.percentile(vals, 97.5)))
    print("correlations done:", gn, flush=True)
pd.DataFrame(rows_c).to_csv(OUT / "01_correlation_stability.csv", index=False)
pd.DataFrame(rows_conv).to_csv(OUT / "02_spearman_convergence.csv", index=False)


def betas(X, y):
    Zx = (X - X.mean(0)) / X.std(0, ddof=1); zy = (y - y.mean()) / y.std(ddof=1)
    return np.linalg.lstsq(np.column_stack([np.ones(len(zy)), Zx]), zy, rcond=None)[0][1:]


rows_b, rows_r = [], []
for gn, g in GROUPS.items():
    X = np.log10(g[PI3].to_numpy())
    for met in MET:
        y = np.log10(g[met].to_numpy()); n = len(y)
        bf = betas(X, y); dom = PI3[int(np.argmax(np.abs(bf)))]
        for frac in FRAC:
            B = np.array([betas(X[i], y[i]) for i in (sample(n, frac) for _ in range(N_REP))])
            top = pd.Series([PI3[j] for j in np.argmax(np.abs(B), 1)]).value_counts(normalize=True)
            full_order = tuple(np.argsort(-np.abs(bf)))
            same_order = float(np.mean([tuple(np.argsort(-np.abs(b))) == full_order for b in B]))
            for j, p in enumerate(PI3):
                rows_b.append(dict(group=gn, metric=met, pi=p, fraction=frac, n=n, beta_full=bf[j],
                                   median=np.median(B[:, j]), q025=np.percentile(B[:, j], 2.5),
                                   q975=np.percentile(B[:, j], 97.5),
                                   sign_pres=float(np.mean(np.sign(B[:, j]) == np.sign(bf[j])))))
            rows_r.append(dict(group=gn, metric=met, fraction=frac, n=n, dominant_full=dom,
                               p_dominant=float(top.get(dom, 0.0)), p_full_ranking=same_order))
    print("betas done:", gn, flush=True)
pd.DataFrame(rows_b).to_csv(OUT / "03_beta_stability.csv", index=False)
pd.DataFrame(rows_r).to_csv(OUT / "04_hierarchy_stability.csv", index=False)
print("done")
