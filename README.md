# River–Aquifer Exchange: flood regulation in π-space

Code and data for:

> Restrepo-Correa, I., Gómez-Vélez, J. D., Jaramillo, M., Álvarez-Villa, O., & Cuéllar, E. *Timing matters: flood
> regulation by river–aquifer exchange across scales.* Manuscript in preparation for *Journal of Hydrology*.

The repository couples MODFLOW-NWT and the Streamflow-Routing package (SFR2) through FloPy to simulate transient
river–aquifer exchange during flood events. A nondimensional formulation of the governing equations yields three
system groups (π₁–π₃) and one forcing group (π₄). A Monte Carlo ensemble of 1,058 simulations, spanning drainage
areas of 10³–10⁵ km², is used to rank the controls of flood regulation and to derive log–log scaling relationships
for five dimensionless response metrics.

Everything needed to reproduce the paper from the simulation results is included: the 1,058 results, the
post-processing, the statistics and the figure scripts. Re-running the simulations needs MODFLOW-NWT.

---

## Repository structure

```
RIVER-AQUIFER-EXCHANGE/
├── src/                          reusable modules
│   ├── geomorphic_scaling.py     channel geometry from drainage area
│   ├── dimensionless_numbers.py  characteristic time and π-groups
│   └── modflow_utils.py          topography, river, hydrograph, SFR2 output
│
├── simulations/
│   ├── montecarlo_sampling.ipynb        Monte Carlo sampling of the configurations
│   ├── build_input_csv.py               main-ensemble inputs from the master Excel
│   ├── build_high_connectivity_inputs.py  high-connectivity ensemble from the main one
│   ├── parallel_modflow.py              the model: one MODFLOW-NWT + SFR2 run per row, in parallel
│   ├── run_campaign.py                  resumable run of one ensemble
│   └── sfr_series_utils.py              time vector and reach series of a finished run
│
├── validation/                   tests of the numerical model (time step, mass balance, GHB, duration)
├── analysis/                     post-processing and statistics
├── figures/                      figure scripts and their output (PDF, PNG 600 dpi, SVG)
│
├── data/
│   ├── input/
│   │   └── model_inputs.csv      model inputs of the 1,058 configurations (both ensembles)
│   └── output/
│       ├── simulation_results.csv      raw results of the 1,058 simulations
│       ├── model_results.csv           results table: scales, π-groups and metrics (post-processing)
│       ├── model_results_dictionary.csv  symbol, description and units of each column of model_results.csv
│       ├── tables/               paper tables
│       └── robustness/           bootstrap and subsampling results
│
├── docs/
│   ├── methodology.md            model, sampling and metric definitions
│   └── analysis_workflow.md      step by step, and which script makes each paper table and figure
├── executables/                  MODFLOW-NWT executable (not tracked; see executables/README.md)
└── legacy/                       archived code (reference only)
```

---

## Installation

Python ≥ 3.10 (tested with 3.12).

```bash
git clone https://github.com/IsaRpo32/RIVER-AQUIFER-EXCHANGE.git
cd RIVER-AQUIFER-EXCHANGE
python -m venv venv
venv\Scripts\activate            # Linux/macOS: source venv/bin/activate
pip install -r requirements.txt
```

Only the simulations need MODFLOW-NWT (Niswonger et al., 2011): see `executables/README.md`.

---

## Reproducing the paper

Every script can be run from any directory; the commands below assume the repository root.

**1. Post-processing** (seconds)

```bash
python analysis/postprocess_metrics.py
```

Reads `data/output/simulation_results.csv` and writes `data/output/model_results.csv` (1,058 rows): inputs,
results, characteristic scales, π-groups and the five metrics.

**2. Statistics**

```bash
python analysis/time_scale_test.py      # choice of the time scale of Δtₚ* and tᵣ*
python analysis/other_scales_test.py
python analysis/correlations.py
python analysis/regressions.py
python analysis/nechi_position.py
python analysis/robustness.py           # several minutes
python analysis/robustness_summary.py
```

**3. Figures**

```bash
python figures/fig03_histograms.py
python figures/fig04_scatter_pi.py
python figures/fig05_scatter_processes.py
python figures/fig06_07_obs_pred.py
python figures/fig08_beta.py
python figures/fig09_nechi.py
python figures/fig10_11_pi_space.py
python figures/fig11b_schematic.py
```

Which script produces each table and figure of the manuscript: `docs/analysis_workflow.md`, Section 7.

**4. Simulations** (optional; hours to days depending on the machine)

```bash
python simulations/run_campaign.py main 8                 # main ensemble, 8 in parallel
python simulations/run_campaign.py high_connectivity 8    # high-connectivity ensemble
```

Each call takes the rows of that ensemble from `data/input/model_inputs.csv`, skips the rows already finished (so
it can be stopped and run again) and writes `data/output/<name>/results.csv` with every finished row, plus a zip to
send the results. Put those rows in `data/output/simulation_results.csv` (same columns, with the `ensemble`
column) and go back to step 1.

---

## Ensembles

Every data file holds both ensembles, told apart by the `ensemble` column:

| `ensemble` | Configurations |
|---|---|
| `main` | 758 |
| `high_connectivity` | 300 |

Manning coefficients are n = 0.025–0.045 (Chow, 1959); the column `n_original` keeps the earlier values, ten
times larger (`run_campaign.py main_n_original` runs the main ensemble with them). The high-connectivity ensemble reuses the main configurations with Kr and Ka of
coarse-sand and gravel beds (C_riv up to 30 m²/s). Its `IDMERGE` values repeat those of the main ensemble, so
rows are keyed by (`ensemble`, `IDMERGE`). The main ensemble was run in two halves of 379 rows; the results
file joins them. The raw SFR2 and listing files of the runs (several GB) are not
included.

---

## Supplementary data

The data files published as supplementary material of the paper are:

| File | Content |
|---|---|
| `data/input/model_inputs.csv` | model inputs of the 1,058 configurations, in the column names the model reads (`simulations/parallel_modflow.py`) |
| `data/output/model_results.csv` | inputs, results, characteristic scales, π-groups and dimensionless metrics of the 1,058 simulations, one row per simulation, with self-explanatory column names |
| `data/output/model_results_dictionary.csv` | symbol, description and units of every column of `model_results.csv` |

The three files have the same rows in the same order. In `model_results.csv` the key is (`ensemble`,
`config_id`); `config_id` is the `IDMERGE` of the other files.

---

## Key definitions

| Quantity | Definition |
|---|---|
| Aquifer response time | t_c = S_y B² / (K_a d_bkf) |
| Kinematic travel time | t_L = L / c_k, with c_k = (5/3)(1/n) S^½ d_bkf^⅔ (bankfull stage) |
| Vertical exchange efficiency | π₁ = t_c K_r / d_r |
| Storage partitioning | π₂ = W_bkf / (S_y B) |
| Kinematic contrast | π₃ = (d_bkf^⅔ S^½ / n) / (B / t_c) |
| Event persistence (forcing) | π₄ = t_c / t_d |
| Peak delay | Δtₚ* = (t_p,out − t_p,in) / t_L |
| Stage attenuation | ΔH* = (H_p,in − H_p,out) / d_bkf |
| Net exchange | ΔLeakage* = ∫(Q_exc,in − Q_exc,out) dt / (B L d_bkf) |
| Return time | tᵣ* = tᵣ / t_L |
| Storage change | ΔStorage* = (V_storage,out − V_storage,in) / (B L d_bkf) |

B is the total aquifer width, with the channel at the centre line. Δtₚ* and tᵣ* are scaled with t_L and not with
t_c, because t_c also enters π₁, π₃ and π₄ and would induce spurious correlation (`analysis/time_scale_test.py`).
Details in `docs/methodology.md` and `docs/analysis_workflow.md`.

---

## References

- Bakker, M., et al. (2016). Scripting MODFLOW model development using Python and FloPy. *Groundwater*, 54(5), 733–739.
- Bieger, K., Rathjens, H., Allen, P. M., & Arnold, J. G. (2015). Development and evaluation of bankfull hydraulic
  geometry relationships for the physiographic regions of the United States. *Journal of the American Water
  Resources Association*, 51(3), 842–858.
- Chow, V. T. (1959). *Open-channel hydraulics*. McGraw-Hill.
- Cooper, H. H., Jr., & Rorabaugh, M. I. (1963). Ground-water movements and bank storage due to flood stages in surface
  streams. U.S. Geological Survey Water-Supply Paper 1536-J.
- Flint, J. J. (1974). Stream gradient as a function of order, magnitude, and discharge. *Water Resources Research*,
  10(5), 969–973.
- Gomez-Velez, J. D., Wilson, J. L., Cardenas, M. B., & Harvey, J. W. (2017). Flow and residence times of dynamic
  riverbank storage and sinuosity-driven hyporheic exchange. *Water Resources Research*, 53(10), 8572–8595.
- Kenney, B. C. (1982). Beware of spurious self-correlations! *Water Resources Research*, 18(4), 1041–1048.
- Niswonger, R. G., & Prudic, D. E. (2005). Documentation of the Streamflow-Routing (SFR2) Package. U.S. Geological
  Survey Techniques and Methods 6-A13.
- Niswonger, R. G., Panday, S., & Ibaraki, M. (2011). MODFLOW-NWT, a Newton formulation for MODFLOW-2005. U.S.
  Geological Survey Techniques and Methods 6-A37.

---

## Citation

```bibtex
@software{river_aquifer_exchange,
  author = {Restrepo-Correa, Isabel},
  title  = {River-Aquifer Exchange: flood regulation in pi-space},
  year   = {2026},
  url    = {https://github.com/IsaRpo32/RIVER-AQUIFER-EXCHANGE}
}
```

## License

MIT License. See [LICENSE](LICENSE).
