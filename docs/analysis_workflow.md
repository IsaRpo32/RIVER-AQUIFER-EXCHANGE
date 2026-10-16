# Analysis workflow (manuscript v15)

From the input files to the tables and figures of the paper. Model details are in `methodology.md`.

## 1. Ensembles

**Main ensemble (758 configurations)**, `data/input/inputs_n_corregido.csv`. Basin-dependent geometry from
drainage area (Bieger et al., 2015; Flint, 1974; Williams, 1978); stochastic aquifer and streambed
properties (Ka, Kr with Ka/Kr log-uniform in [10, 1000], Sy, nB, nD); three flood durations per
drainage-area class (td = {0.7, 1.0, 1.3} × 0.6 Tc, Kirpich). Streambed conductance Criv = 10⁻⁴–0.2 m²/s;
100–200 rows along the reach. A kernel-density subset in π₁–π₂–π₃ space gives the 758 configurations.
Manning coefficients n = 0.025–0.045 (an earlier campaign with values ten times larger was discarded).
The ensemble was run in two halves of 379 rows (`inputs_n_corregido_mitad_A.csv`, `..._mitad_B.csv`).

**High-connectivity ensemble (300 configurations)**, `data/input/inputs_high_connectivity.csv`. Same base
configurations (geometry, Sy, forcing, n); Kr and Ka reassigned for coarse-sand and gravel beds: Kr
log-uniform in [10⁻⁴, 10⁻²/r], r = Ka/Kr log-uniform in [1, 10], Ka ≤ 10⁻² m/s, Criv recomputed with the
same Δx. Configurations with Kr > 10⁻³ m/s (clean gravel) are excluded and 100 configurations are drawn
per drainage-area class (seed 20261012). Criv = 0.2–30 m²/s. Built by
`simulations/build_high_connectivity_inputs.py`.

The `IDMERGE` of a high-connectivity row is that of the main row it comes from, so `IDMERGE` alone is not
a unique key: use (`ensemble`, `IDMERGE`).

## 2. Simulations

`simulations/run_campaign.py <input CSV> [workers]` runs `simulations/parallel_modflow.py` on every row
(MODFLOW-NWT + SFR2), resumable, with the validated time-step scheme. Main features of the model:

- four layers (lowest inactive), 65 columns across B (channel in the central column), 100–200 rows along L;
- one-day steady-state period, then a transient period of 90, 150 or 400 days for 10³, 10⁴ and 10⁵ km²
  (interpolated in log Ad), until the exchange flux decays to about 1 % of its peak;
- adaptive time steps: Δt = min(td/400, 50Δx/ck) during the flood and one further event duration, then
  growing by a factor of 1.05 up to 3 h;
- Cooper–Rorabaugh stage hydrograph (tp = 0.5 td; asymmetry parameter at its lower limit 0.1ω) converted
  to inflow with Manning's equation for a wide channel;
- longitudinal boundaries: GHB with C = Ka Wbkf Δx / D;
- exported per run: peak times and stages at inlet and outlet, return time, cumulative budget terms
  (STORAGE, STREAM LEAKAGE, HEAD DEP BOUNDS) and mass-balance discrepancy.

The results used in the paper are in `data/output/campaign/`:

| File | Ensemble | Original name |
|---|---|---|
| `main_half_A.csv` | main ensemble, half A (379) | `complete_sample_points_parallel_20261011_230012.csv` |
| `main_half_B.csv` | main ensemble, half B (379) | `complete_sample_points_parallel_20261011_185710.csv` |
| `high_connectivity.csv` | high-connectivity ensemble (300) | `complete_sample_points_parallel_20261012_134235.csv` |

All 1,058 runs finished (`success = True`). Half A was run on Linux and half B on Windows with the same
script; both executables give the same results within 0.00004 %. To use a new run, copy the `results.csv`
that `run_campaign.py` writes in `data/output/<input name>/` over the corresponding file above.

Validation of the numerical model (time step, mass balance, GHB, duration, reproduction of the master
Excel) is in `validation/`.

## 3. Post-processing (`analysis/postprocess_metrics.py`)

Reads the three result files and writes `data/output/metrics.csv` (1,058 rows):

- t_c = Sy B² / (Ka d_bkf); V_bkf = d_bkf^⅔ S^½ / n; t_L = L / (5/3 V_bkf);
- π₁ = t_c Kr / dr; π₂ = Wbkf / (Sy B); π₃ = V_bkf t_c / B; π₄ = t_c / td;
- Δtₚ* = Δtₚ / t_L; ΔH* = ΔH / d_bkf; ΔLeakage* = (LEAKAGE_IN − LEAKAGE_OUT) / (B L d_bkf); tᵣ* = tᵣ / t_L;
  ΔStorage* = (STORAGE_OUT − STORAGE_IN) / (B L d_bkf);
- `ensemble`, `group` (`main` / `high_connectivity`) and `area_class` (nearest of 10³, 10⁴, 10⁵ km² in log).

Columns `delta_tp*_tc` and `t_r*_tc` keep the earlier normalization by t_c for comparison only.

**Choice of t_L.** t_c enters π₁, π₃ and π₄. With Δtₚ permuted at random among the simulations, Δtₚ/t_c
still gives R² = 0.87 against π₁–π₃ (0.84 for tᵣ), against 0.07 with t_L (0.05 for tᵣ)
(`time_scale_test.py`; Kenney, 1982). ΔH*, ΔLeakage* and ΔStorage* do not show this problem
(`other_scales_test.py`).

## 4. Statistical analysis

| Script | What it does | Output (`data/output/`) |
|---|---|---|
| `correlations.py` | Pearson and Spearman on log₁₀ values; correlations with π₁ use the main ensemble only (n = 758), as in the figures | `tables/correlations_pi_metrics.csv` |
| `regressions.py` | predictive models log₁₀ Y* = a₀ + Σ aᵢ log₁₀ πᵢ with π₁–π₄ (OLS; 10-fold cross-validated R²); hierarchy of controls: standardized model with π₁–π₃ only (π₄ describes the forcing and is collinear with π₃ through t_c, r = 0.96), 95 % bootstrap intervals (3,000 replicates) | `tables/predictive_models_pi1_pi4.csv`, `tables/hierarchy_beta_pi1_pi3.csv` |
| `robustness.py` | bootstrap and subsampling (100, 90, 80, 70 and 50 %), 3,000 replicates, pooled and by drainage-area class (several minutes) | `robustness/01_ … 04_*.csv` |
| `robustness_summary.py` | numbers quoted in the text from the robustness tables | (printed) |
| `nechi_position.py` | Nechí corridor: π-groups and metrics, percentiles in the ensemble, nearest-neighbour distance in π-space, model predictions | `tables/nechi_position.csv` |

The Nechí corridor lies at the edge of the ensemble (d_bkf = 3.66 m against a sampled maximum of 3.48 m;
distance to the nearest simulation in π-space 0.70 against a 95th percentile of 0.47 within the
ensemble), so its predictions are extrapolations.

## 5. Figures (`figures/`)

`style.py` sets the common style (Arial, single grey for data). Figures with π₁ on an axis use the main
ensemble only (n = 758). `fig03_histograms.py` also writes `tables/pi_statistics.csv` (statistics of t_c,
t_L and π₁–π₄ by drainage-area class). `fig10_11_pi_space.py` also writes a 2D alternative of Fig. 10
(`fig10_deltaH_pi_space_2D`), not used in the paper. The files in `figures/output/` were produced on
Windows with Arial; on a machine without Arial the same figures are drawn with DejaVu Sans.

## 6. Paper tables and figures

Numbering of manuscript v15.

| Paper | Produced by | File |
|---|---|---|
| Fig. 3 | `figures/fig03_histograms.py` | `figures/output/fig03_pi_histograms.*` |
| Table 5 | `figures/fig03_histograms.py` | `data/output/tables/pi_statistics.csv` |
| Figs. 4–5 | `figures/fig04_scatter_pi.py`, `fig05_scatter_processes.py` | `figures/output/fig04_*`, `fig05_*` |
| Table 6 | `analysis/regressions.py` | `data/output/tables/predictive_models_pi1_pi4.csv` |
| Table 7 | `analysis/regressions.py` | `data/output/tables/hierarchy_beta_pi1_pi3.csv` |
| Figs. 6–7 | `figures/fig06_07_obs_pred.py` | `figures/output/fig06_*`, `fig07_*` |
| Fig. 8 | `figures/fig08_beta.py` | `figures/output/fig08_beta_pi1_pi3.*` |
| Tables 8–9, Fig. 9 | `analysis/nechi_position.py`, `figures/fig09_nechi.py` | `data/output/tables/nechi_position.csv`, `figures/output/fig09_*` |
| Figs. 10–11a | `figures/fig10_11_pi_space.py` | `figures/output/fig10_deltaH_pi_space_3D.*`, `fig11a_*` |
| Fig. 11b | `figures/fig11b_schematic.py` | `figures/output/fig11b_regulation_schematic.*` |
| Supplementary: correlations, robustness | `analysis/correlations.py`, `analysis/robustness.py` | `data/output/tables/correlations_pi_metrics.csv`, `data/output/robustness/` |
