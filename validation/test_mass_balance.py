"""
Time-step test focused on the MODFLOW volumetric budget (.list), with corrected Manning n.

Reference = 11 000 uniform steps (the full old Courant step would need 42 000 - 91 000 steps for these
rows, 0.5-1.5 h each). Compared against the adaptive schemes fine / default / coarse.
The budget is printed at every step so the per-step discrepancy and the rate series can be compared.

    python validation/test_mass_balance.py            # run the simulations (~10-15 min)
    python validation/test_mass_balance.py --analyze  # tables + figures from the finished runs
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "simulations"))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
from concurrent.futures import ProcessPoolExecutor, as_completed

from parallel_modflow import DATA_INPUT, DEFAULT_EXE, run_single_simulation

INPUT_CSV = DATA_INPUT / "inputs_main.csv"

OUT_DIR = ROOT / "data" / "output" / "validation"
RUNS = str(OUT_DIR / "runs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
IDS = [41, 96, 612]            # 1000, 10000, 100000 km2: cheapest rows for the old Courant scheme
BASE = RUNS
SUBDIR = "mass_balance_test"
OUT = str(OUT_DIR / "figuras_balance_masas")
COMMON = {'save_cbc': False, 'print_budget': True}
SCHEMES = {
    'uniform_11k': {'mode': 'courant', 'nstp_uniform': 11000},
    'fine': {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0},
    'default': {},
    'coarse': {'n_hyd': 100, 'tsmult': 1.2, 'dt_max_hours': 12.0},
}
TERMS = ['STORAGE', 'HEAD_DEP_BOUNDS', 'STREAM_LEAKAGE']
COLORS = {'uniform_11k': '#2a78d6', 'fine': '#eb6834', 'default': '#1baf7a', 'coarse': '#eda100'}


def list_path(scheme, sim_id):
    return os.path.join(BASE, SUBDIR, scheme, f"run_{sim_id:03d}", f"water_sim_bagre_{sim_id}.list")


def read_budget(scheme, sim_id):
    """Rate and cumulative budgets at every step, indexed by elapsed time [days] (steady period excluded)."""
    rate, cum = flopy.utils.MfListBudget(list_path(scheme, sim_id)).get_dataframes(start_datetime=None)
    return rate.iloc[1:], cum.iloc[1:]


def run():
    dt = pd.read_csv(INPUT_CSV).set_index("IDMERGE", drop=False)
    jobs = [(s, (i, dt.loc[i], BASE, f"{SUBDIR}/{s}", DEFAULT_EXE, {**COMMON, **cfg}))
            for s, cfg in SCHEMES.items() for i in IDS]
    results = []
    with ProcessPoolExecutor(max_workers=len(jobs)) as ex:
        futures = {ex.submit(run_single_simulation, p): s for s, p in jobs}
        for f in as_completed(futures):
            r = f.result(); r["scheme"] = futures[f]; results.append(r)
            print(f"{r['scheme']:12s} {r['simulation_id']}: {'OK' if r['success'] else r.get('error')} "
                  f"({r['execution_time']/60:.1f} min)", flush=True)
    pd.DataFrame(results).to_csv(OUT_DIR / "test_mass_balance_results.csv", index=False)


def analyze():
    os.makedirs(OUT, exist_ok=True)
    res = pd.read_csv(OUT_DIR / "test_mass_balance_results.csv").set_index(["simulation_id", "scheme"])
    rows = []
    for sim_id in IDS:
        rate_ref, cum_ref = read_budget('uniform_11k', sim_id)
        for scheme in SCHEMES:
            rate, cum = read_budget(scheme, sim_id)
            row = {'ID': sim_id, 'scheme': scheme, 'nstp': len(rate),
                   'minutes': res.loc[(sim_id, scheme), 'execution_time'] / 60,
                   'max_step_discrepancy_%': rate['PERCENT_DISCREPANCY'].abs().max(),
                   'final_cum_discrepancy_%': cum['PERCENT_DISCREPANCY'].iloc[-1]}
            for term in TERMS:
                for side in ('IN', 'OUT'):
                    col = f"{term}_{side}"
                    # final cumulative volume vs reference
                    row[f"errV_{col}_%"] = 100 * abs(cum[col].iloc[-1] - cum_ref[col].iloc[-1]) / cum_ref[col].iloc[-1]
                    # rate series: reference interpolated at this scheme's times, relative L2 error
                    ref_at = np.interp(rate.index, rate_ref.index, rate_ref[col])
                    row[f"errQ_{col}_%"] = 100 * np.linalg.norm(rate[col] - ref_at) / np.linalg.norm(ref_at)
            for m in ['Qexc_in', 'Qexc_out', 'T_in', 'T_out', 'delta_tp', 'delta_H']:
                ref_v = res.loc[(sim_id, 'uniform_11k'), m]
                row[f"err_{m}_%"] = 100 * abs(res.loc[(sim_id, scheme), m] - ref_v) / abs(ref_v)
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT_DIR / "test_mass_balance_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(t[['ID', 'scheme', 'nstp', 'minutes', 'max_step_discrepancy_%', 'final_cum_discrepancy_%']].round(3).to_string(index=False))
    print("\nFinal cumulative volumes, % error vs uniform_11k")
    print(t[['ID', 'scheme'] + [c for c in t if c.startswith('errV_')]].round(3).to_string(index=False))
    print("\nRate time series, relative L2 error % vs uniform_11k")
    print(t[['ID', 'scheme'] + [c for c in t if c.startswith('errQ_')]].round(3).to_string(index=False))
    print("\nExported metrics, % error vs uniform_11k")
    print(t[['ID', 'scheme'] + [c for c in t if c.startswith('err_')]].round(3).to_string(index=False))

    plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e4e0",
                         "savefig.dpi": 300, "savefig.bbox": "tight"})
    # Fig A: per-step discrepancy over time
    fig, axes = plt.subplots(1, len(IDS), figsize=(10, 3), sharey=True)
    for ax, sim_id in zip(axes, IDS):
        for scheme in SCHEMES:
            rate, _ = read_budget(scheme, sim_id)
            ax.plot(rate.index, rate['PERCENT_DISCREPANCY'], lw=1 if scheme != 'uniform_11k' else .6,
                    color=COLORS[scheme], label=f"{scheme} ({len(rate)} steps)")
        ax.set_title(f"Row {sim_id}"); ax.set_xlabel("Time (days)")
        ax.legend(frameon=False, fontsize=6)
    axes[0].set_ylabel("Budget discrepancy in the step (%)")
    fig.suptitle("Mass-balance error at every time step (.list)", y=1.03)
    fig.tight_layout(); fig.savefig(f"{OUT}/figA_discrepancy_per_step.png"); plt.close(fig)

    # Fig B: rate series of each budget term (one row per term), first two weeks + full window inset
    for sim_id in IDS:
        fig, axes = plt.subplots(len(TERMS), 2, figsize=(9, 7))
        for k, term in enumerate(TERMS):
            for j, side in enumerate(('IN', 'OUT')):
                ax = axes[k, j]
                for scheme in SCHEMES:
                    rate, _ = read_budget(scheme, sim_id)
                    ax.plot(rate.index, rate[f"{term}_{side}"], color=COLORS[scheme],
                            lw=.7 if scheme == 'uniform_11k' else 1.1, label=scheme,
                            ls='-' if scheme == 'uniform_11k' else '--')
                ax.set_title(f"{term}_{side}", fontsize=8); ax.set_ylabel("m³/s")
                if k == len(TERMS) - 1:
                    ax.set_xlabel("Time (days)")
        axes[0, 0].legend(frameon=False, fontsize=6)
        fig.suptitle(f"Row {sim_id}: budget terms (rates) by scheme", y=1.01)
        fig.tight_layout(); fig.savefig(f"{OUT}/figB_budget_terms_{sim_id}.png"); plt.close(fig)
    print(f"\nFigures in {OUT}/")


if __name__ == "__main__":
    analyze() if "--analyze" in sys.argv else run()
