# Validation of `simulations/parallel_modflow.py`

Scripts used to choose and check the time discretization, the transient duration and the exported metrics.
They run the model through `run_single_simulation` with the inputs in `data/input/` and write their runs,
tables and figures to `data/output/validation/`. Run them from the repository root, e.g.
`python validation/test_time_steps.py`.

| Script | What it checks | Result (October 2026) |
|---|---|---|
| `test_reproduce_excel.py <master.xlsx>` | Current script with the **original** Manning n against the master Excel (6 rows of campaign 10092025) | Qexc_in, t_r, H_ratio, ΔH and budget terms within 0–4 %. Qexc_out in the Excel is up to 18 % low because the old runs stopped before the return finished |
| `test_time_steps.py` | Adaptive time-step schemes against a very fine reference (3 rows per basin size) | `fine` (`n_hyd=400`, `tsmult=1.05`, `dt_max_hours=3`) is the cheapest scheme within tolerance |
| `test_mass_balance.py` (`--analyze`) | Volumetric budget of the `.list` against 11 000 uniform steps | `fine`, about 1 700 steps: volumes ≤ 0.6 %, rates ≤ 1.9 %, metrics ≤ 3.4 %; cumulative discrepancy ≤ 0.04 % |
| `test_ghb.py` | GHB conductance ×0.1 / ×1 / ×10 | Exchange changes ≤ 5 %; boundary flow up to 61 % |
| `test_duration.py` (`--analyze`) | 730-day runs: when the net exchange falls below 1 % of its peak | ≤ 19 / 77 / 322 d for 10³ / 10⁴ / 10⁵ km², hence durations of 90 / 150 / 400 d |
| `test_schemes_examples.py` | Old uniform Courant step against the adaptive schemes on two rows with the original n | Used for the figures comparing schemes |

`simulations/sfr_series_utils.py` rebuilds the time vector of a finished run and reads the depth and
flow series of the inlet and outlet reaches.
