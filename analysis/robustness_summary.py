"""
Summary of data/output/robustness/ (the numbers quoted in the text).

    python analysis/robustness_summary.py      (after robustness.py)
"""
from pathlib import Path

import pandas as pd

R = Path(__file__).resolve().parents[1] / "data" / "output" / "robustness"
c = pd.read_csv(R / "01_correlation_stability.csv")
v = pd.read_csv(R / "02_spearman_convergence.csv")
b = pd.read_csv(R / "03_beta_stability.csv")
j = pd.read_csv(R / "04_hierarchy_stability.csv")
c70 = c[(c.fraction == 0.7) & (c.method == "spearman")]
strong = c70[c70.full.abs() >= 0.3]
print("Spearman, 70 %: pairs with |rho|>=0.3 ->", len(strong), "| min sign preservation:", strong.sign_pres.min().round(4),
      "| max |median-full|:", (strong["median"] - strong.full).abs().max().round(4))
print("  all pairs with sign preservation < 0.95:")
print(c70[c70.sign_pres < 0.95][["group", "pi", "metric", "full", "sign_pres"]].round(3).to_string(index=False))
print("\nconvergence (Spearman): |median - full| by sample size")
print(v.groupby(["group", "n_sample"]).abs_diff.max().unstack(0).round(4).to_string())
b70 = b[b.fraction == 0.7]
print("\nbeta 70 %: min sign preservation by group:", b70.groupby("group").sign_pres.min().round(3).to_dict())
print(b70[b70.sign_pres < 0.95][["group", "metric", "pi", "beta_full", "sign_pres"]].round(3).to_string(index=False))
print("\nhierarchy: P(dominant = that of the full fit) and P(same full ranking), 70 % and 50 %")
print(j[j.fraction.isin([0.7, 0.5])].pivot_table(index=["group", "metric"], columns="fraction",
                                                   values=["p_dominant", "p_full_ranking"]).round(3).to_string())
