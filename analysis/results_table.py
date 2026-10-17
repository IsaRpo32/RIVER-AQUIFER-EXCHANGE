"""
The results table of the paper, data/output/model_results.csv (one row per simulation; columns described in
data/output/model_results_dictionary.csv).

postprocess_metrics.py writes it; every other script in analysis/ and figures/ reads it with load_results(),
which returns the short column names used in the code (pi1, delta_tp*, tc, ...).
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "output" / "model_results.csv"

# column of model_results.csv -> name used in the code (also the name in simulation_results.csv)
SHORT = {
    "config_id": "IDMERGE", "ensemble": "ensemble", "area_class_km2": "area_class_km2",
    "drainage_area_km2": "da_km2", "W_bkf_m": "w_bkf", "d_bkf_m": "d_bkf", "L_m": "L", "B_m": "B", "D_m": "D",
    "d_r_m": "dr", "S": "S", "n_manning": "n", "S_y": "Sy", "K_a_m_s": "Ka", "K_r_m_s": "Kr",
    "Ka_over_Kr": "Ka_Kr", "C_riv_m2_s": "C", "n_rows": "nrow", "t_d_days": "t_d_days", "H0_m": "H_0",
    "Hp_m": "H_p", "sim_duration_days": "duration_days", "n_time_steps": "nstp_total",
    "mass_balance_error_pct": "PERCENT_DISCREPANCY", "t_p_in_s": "t_p_in", "t_p_out_s": "t_p_out",
    "H_p_in_m": "H_p_in", "H_p_out_m": "H_p_out", "delta_tp_s": "delta_tp", "delta_H_m": "delta_H",
    "t_r_s": "t_r", "stream_leakage_in_m3": "STREAM_LEAKAGE_IN", "stream_leakage_out_m3": "STREAM_LEAKAGE_OUT",
    "storage_in_m3": "STORAGE_IN", "storage_out_m3": "STORAGE_OUT", "ghb_in_m3": "HEAD_DEP_BOUNDS_IN",
    "ghb_out_m3": "HEAD_DEP_BOUNDS_OUT", "total_in_m3": "TOTAL_IN", "total_out_m3": "TOTAL_OUT",
    "exchange_end_ratio": "Qexc_end_ratio", "t_c_s": "tc", "V_bkf_m_s": "V_bkf", "t_L_s": "tL",
    "pi1": "pi1", "pi2": "pi2", "pi3": "pi3", "pi4": "pi4",
    "delta_tp_star": "delta_tp*", "delta_H_star": "delta_H*", "delta_leakage_star": "delta_leakage*",
    "t_r_star": "t_r*", "delta_storage_star": "delta_storage*",
}
AREA_LABEL = {1000.0: "1e3", 10000.0: "1e4", 100000.0: "1e5"}


def load_results() -> pd.DataFrame:
    """model_results.csv with short column names, plus td_s (event duration in s), group (= ensemble:
    'main' or 'high_connectivity') and area_class ('1e3', '1e4', '1e5')."""
    d = pd.read_csv(RESULTS).rename(columns=SHORT)
    d["td_s"] = d.t_d_days * 86400
    d["group"] = d.ensemble
    d["area_class"] = d.area_class_km2.map(AREA_LABEL)
    return d
