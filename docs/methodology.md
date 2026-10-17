# Methodology

## Overview

This document describes the methodology for modeling transient river-aquifer exchange during flood events using a combination of geomorphic scaling, dimensionless analysis, and numerical groundwater modeling.

## 1. Geomorphic Scaling

### 1.1 Regional Curves

Channel geometry is estimated from drainage area using regional curves based on empirical relationships from the literature.

#### Bankfull Width (Bieger et al., 2015)

```
w_bkf = 2.7 × A_d^0.352
```

Where:
- `w_bkf`: Bankfull width (m)
- `A_d`: Drainage area (km²)

This is the relation of the manuscript (Table 2) and of `src/geomorphic_scaling.get_w_bkf`.
Note: 300 of the 758 rows in `data/input/` (campaigns `25082025` and `30082025`) were generated
with an earlier coefficient, `w_bkf = 2.9 × A_d^0.352`; the other 458 rows (`10092025`,
`13092025_big_basins`) use 2.7. Each row is run with the `w_bkf` stored in the CSV, so every
simulation is internally consistent, but the ensemble does not follow a single width curve.

#### Bankfull Depth (Bieger et al., 2015)

```
d_bkf = 0.3 × A_d^0.213
```

Where:
- `d_bkf`: Bankfull depth (m)

#### Channel Slope (Lane et al., 2021; Zhou et al., 2024)

```
S = 0.001 × A_d^-0.35
```

Where:
- `S`: Channel slope (m/m); 1.78×10⁻⁵ to 8.91×10⁻⁵ for 10³ to 10⁵ km²

### 1.2 Hack's Law

Hydraulic length is estimated using Hack's Law:

```
L_h = 1.4 × A_d^0.6 × 1000
```

Where:
- `L_h`: Hydraulic length (m)
- `A_d`: Drainage area (km²)

## 2. Hydrograph Generation

### 2.1 Cooper-Rorabaugh Model

The flood hydrograph follows the Cooper and Rorabaugh (1963) analytical solution for stream stage variation:

```
h(t) = H_0 + H_p × exp(-β(t - t_p)) × [(1 - cos(ωt)) / (1 - cos(ωt_p))]
```

Where:
- `H_0`: Initial (base) stage (m)
- `H_p`: Peak stage increment (m)
- `t_p`: Time to peak (s)
- `ω = 2π/t_d`: Angular frequency
- `β = ω·cot(ω·t_p/2)`: Decay coefficient (asymmetry), with a lower limit `β ≥ 0.1ω`

Outside `0 ≤ t ≤ t_d` the stage is `H_0`. In the simulations `H_0 = 0.3·d_bkf`,
`H_p = 0.7·d_bkf` (so the peak stage is about `d_bkf`) and `t_p = t_d/2`.

**Lower limit of β.** With `t_p = t_d/2`, `ω·t_p/2 = π/2` and `cot(π/2) = 0`, so the limit
`β = 0.1ω` applies in every run. The hydrograph is therefore slightly skewed: its maximum lies
at `0.468·t_d` instead of `t_p = t_d/2`, and reaches `1.007·d_bkf`. For this reason `H_p,in` and `t_p,in` (Section 6) are
taken from the hydrograph as injected, not from the nominal `H_p` and `t_p`. The limit was kept
on purpose.

### 2.2 Inflow to SFR2

SFR2 receives a discharge, not a stage. The stage hydrograph is converted with Manning's
equation for a wide rectangular channel, `Q = (1/n)·h^(5/3)·S^(1/2)·w_bkf`, sampled every
600 s and written to the tabfile `QHidrom.dat`. The flood starts one day after the steady
period (`t_delay = 1 d`).

**Base flow after the event.** After `t_d` the table continues at the base flow `Q(H_0)` until
one day past the end of the run. Without these rows, SFR2 averaged the table over each time
step and took zero flow after the last entry, which gave a one-step drop of the inflow at
`t_d`. The extra rows change the exchange by ≤ 0.01 % and do not change the hydrograph.

### 2.3 Event Duration Scaling

Event duration is scaled with drainage area using the NRCS method:

```
t_d = 0.6 × T_c
```

Where `T_c` is the time of concentration estimated using Kirpich's formula:

```
T_c = 0.0195 × L^0.77 × S^-0.385 × 60
```

## 3. Dimensionless Analysis

Nondimensionalizing the governing equations with the length scales `B` and `d_bkf` gives three system
groups and one forcing group (computed in `analysis/postprocess_metrics.py`):

| Group | Definition | Meaning |
|-------|------------|---------|
| Aquifer response time | `t_c = S_y·B² / (K_a·d_bkf)` | characteristic time scale |
| π₁ | `t_c·K_r / d_r` | vertical exchange efficiency |
| π₂ | `w_bkf / (S_y·B)` | storage partitioning |
| π₃ | `(d_bkf^(2/3)·S^(1/2) / n) / (B / t_c)` | kinematic contrast (bankfull Manning velocity over aquifer response velocity) |
| π₄ | `t_c / t_d` | event persistence (forcing) |

`B` is the total aquifer width, with the channel in the central column. `t_c` enters π₁, π₃ and π₄, so
π₃ and π₄ are strongly collinear (r = 0.96 in log space); the hierarchy of controls is therefore ranked
with π₁–π₃ only, and π₄ is used in the predictive models (`analysis/regressions.py`).

The kinematic travel time along the reach, `t_L = L / c_k` with `c_k = (5/3)·d_bkf^(2/3)·S^(1/2)/n`,
is the scale of the peak delay and the return time (Section 6).

## 4. Numerical Modeling

### 4.1 MODFLOW-NWT Configuration

The groundwater model uses MODFLOW-NWT with the following packages:

| Package | Purpose |
|---------|---------|
| DIS | Discretization: 4 layers, 65 columns across `B`, 100–200 rows along `L` |
| BAS | Active cells and initial heads |
| UPW | Aquifer properties (`K_a`, `K_v = K_a`, `S_y`, `S_s = 10⁻⁷ m⁻¹`, convertible layers) |
| NWT | Newton solver (`headtol = 1e-5`, `fluxtol = 1000`) |
| SFR2 | Streamflow routing with kinematic-wave routing (`IRTFLG = 3`, `NUMTIM = 50`, `WEIGHT = 0.75`) |
| GHB | Upstream and downstream aquifer boundaries |
| OC | Output control |

### 4.2 Model Domain

```
        B (aquifer width)
    ←───────────────────────→
    ┌───────────────────────┐ ↑
    │                       │ │
    │     AQUIFER           │ │ L (length)
    │        ║ RIVER        │ │
    │        ║              │ │
    │        ║              │ │
    └───────────────────────┘ ↓
```

### 4.3 Boundary Conditions

- **River (SFR2)**: one segment along the central column; inflow hydrograph at the upstream end
  (Section 2.2), streambed conductance `K_r·w_bkf·Δx / d_r`
- **Upstream and downstream boundaries** (first and last rows): GHB in every layer, with head at the
  base-flow river level (minus `S·L` downstream) and conductance `K_a·w_bkf·Δx / D`. A sensitivity test
  with the conductance ×0.1 and ×10 changes the exchange by ≤ 5 % (`validation/test_ghb.py`)
- **Lateral boundaries** (first and last columns): inactive cells, i.e. no-flow
- **Bottom**: the lowest of the four layers is inactive, so the base is no-flow and the active
  thickness is somewhat less than `D`
- **Top**: water table (unconfined)

### 4.4 Time Discretization

The model has a steady-state stress period followed by transient periods built for each row by
`build_adaptive_periods` (`simulations/parallel_modflow.py`):

| Period | Length | Time step |
|--------|--------|-----------|
| SP1 | delay + flood + one more `t_d` (`t_delay + 2·t_d`) | uniform fine step `Δt_fine` |
| SP2 | until the step reaches `Δt_max` | grows geometrically (`tsmult`) |
| SP3 | rest of the run | uniform `Δt_max` |

`Δt_fine = min(t_d / n_hyd, numtim·Δy / c)`, where `c = (5/3)·v_max` is the kinematic-wave
celerity and `numtim = 50` the SFR2 routing sub-steps. The campaign uses the scheme validated as
`fine`: `n_hyd = 400`, `tsmult = 1.05`, `Δt_max = 3 h`, which gives about 1 700–3 700 steps per run.

This replaces the earlier uniform step set by the Courant criterion (4 000–250 000 steps per run).
Against a reference with 11 000 uniform steps, the `fine` scheme reproduces the metrics within
3.4 %, and the cumulative mass-balance discrepancy is ≤ 0.04 % (`validation/test_time_steps.py`,
`validation/test_mass_balance.py`).

### 4.5 Transient Length

The transient simulation lasts 90, 150 and 400 days for 10³, 10⁴ and 10⁵ km², interpolated in
`log10(A_d)` for intermediate areas (`duration_by_basin`). In 730-day test runs the net exchange
fell below 1 % of its peak after at most 19, 77 and 322 days for the three basin sizes
(`validation/test_duration.py`). Earlier runs used 215 days (or 90/150), which in some rows cut
the return flow before it finished.

## 5. Monte Carlo Simulation

### 5.1 Parameter Sampling

Geometry scales with drainage area (Section 1); aquifer and streambed properties are sampled
(`simulations/montecarlo_sampling.ipynb`). The 758 configurations of the main ensemble come from the
master Excel through `simulations/build_input_csv.py`, and span:

| Parameter | Main ensemble (758) | Distribution |
|-----------|---------------------|--------------|
| K_a | 5.8×10⁻⁷ – 2.6×10⁻² m/s | K_a = r·K_r, r = K_a/K_r log-uniform in [10, 1000] |
| K_r | 10⁻⁸ – 9.5×10⁻⁵ m/s | Log-uniform |
| S_y | 0.10 – 0.30 | Uniform |
| n (Manning) | 0.025 – 0.045 | Uniform (Section 5.4) |
| B/w_bkf | 13 – 200 | Uniform |
| D/d_bkf | 1 – 20 | Uniform |
| d_r/d_bkf, H_0/d_bkf | 0.3 | Fixed |
| t_d | {0.7, 1.0, 1.3} × 0.6·T_c (1.3 – 38 d) | Three durations per area class |

### 5.2 Constraint-Based Filtering

Samples are filtered to ensure numerical stability:

```
100 < L/dx < 200  (grid resolution constraint)
C = K_r × w_bkf × dx / d_r  (streambed conductance)
```

### 5.3 Density-Based Sampling

Final samples are selected using kernel density estimation (KDE) in dimensionless parameter space (π₁, π₂, π₃) to ensure coverage of the parameter space.

### 5.4 Manning Coefficient

The first sampling drew `n` from 0.25–0.45. Those values are about ten times too high for a
channel: they match the sheet-flow values of NRCS TR-55 (overland flow), not natural
streams, for which `n ≈ 0.025–0.045` (Chow, 1959; Arcement & Schneider, 1989). The corrected
inputs (rows `main` of `data/input/model_inputs.csv`) keep the same 758 draws divided by 10 and store the
original value in `n_original`, so the earlier campaign can be run for comparison
(`run_campaign.py main_n_original`). They are rebuilt from the master Excel by `simulations/build_input_csv.py`.

### 5.5 High-Connectivity Ensemble

The main ensemble has low streambed conductances (`C_riv` = 10⁻⁴ – 0.2 m²/s). A second ensemble of 300
configurations (rows `high_connectivity` of `data/input/model_inputs.csv`) keeps the geometry, `S_y`, forcing and `n`
of the main configurations and reassigns `K_r` and `K_a` for coarse-sand and gravel beds: `K_a/K_r`
log-uniform in [1, 10], `K_r` log-uniform in [10⁻⁴, 10⁻²·K_r/K_a] m/s; rows with `K_r > 10⁻³ m/s`
(clean gravel) are dropped and 100 configurations are drawn per drainage-area class. This gives
`C_riv` = 0.2 – 30 m²/s. The file is rebuilt by `simulations/build_high_connectivity_inputs.py` (see its
docstring for how it reproduces the simulated file). Its `IDMERGE` values repeat those of the main rows
the configurations come from, so results are keyed by (ensemble, `IDMERGE`).

## 6. Response Metrics

The five response metrics. `parallel_modflow.py` exports the dimensional quantities for each row; the
normalization (the metrics marked with `*`) is done by `analysis/postprocess_metrics.py`. All times are
measured from the start of the transient period.

| Metric | Definition | Exported columns |
|--------|------------|------------------|
| Peak delay | `Δt_p* = (t_p,out − t_p,in) / t_L` | `delta_tp`, `t_p_in`, `t_p_out` |
| Stage attenuation | `ΔH* = (H_p,in − H_p,out) / d_bkf` | `delta_H`, `H_p_in`, `H_p_out`, `H_ratio` |
| Flood-exchange efficiency | `ΔLeakage* = ∫(Q_exc,in − Q_exc,out) dt / (B·L·d_bkf)` | `STREAM_LEAKAGE_IN`, `STREAM_LEAKAGE_OUT` |
| Return time | `t_r* = t_r / t_L` | `t_r` |
| Storage change | `ΔStorage* = (V_storage,out − V_storage,in) / (B·L·d_bkf)` | `STORAGE_OUT`, `STORAGE_IN` |

**Time scale of `Δt_p*` and `t_r*`.** Earlier versions of the manuscript (up to v14) divided both by
`t_c`. Because `t_c` also enters π₁, π₃ and π₄, that normalization creates correlation by itself: with
`Δt_p` shuffled at random among the simulations, `Δt_p/t_c` still gives R² = 0.87 against π₁–π₃ (0.84
for `t_r`), against 0.07 with `t_L` (0.05 for `t_r`) (`analysis/time_scale_test.py`; Kenney, 1982).
Since manuscript v15 both are scaled with the kinematic travel time `t_L`; the `t_c` versions are computed
only in `time_scale_test.py`. `ΔH*`, `ΔLeakage*` and `ΔStorage*` do not show this
problem (`analysis/other_scales_test.py`). Note that `Δt_p` is close to `t_L` itself (`Δt_p*` ≈ 1–3), so
the peak delay mostly reflects the travel time of the wave along the reach.

How each quantity is obtained (`flood_wave_metrics`, `return_time`):

- **`t_p,in`, `H_p,in`**: maximum of the prescribed upstream stage hydrograph, as injected (with
  the β limit), refined with a parabola through the maximum and its two neighbours.
- **`t_p,out`**: peak of the inflow to the outlet reach. SFR2 computes the depth from that flow
  with Manning, so both peak at the same time, but the flow is printed with five significant
  digits and the depth only to 1 mm.
- **`H_p,out`**: maximum simulated depth of the outlet reach (resolution 1 mm).
- **`t_r`**: first time after the positive exchange peak at which the net river–aquifer exchange
  changes sign, linearly interpolated (`NaN` if it never changes).
- **Leakage and storage**: cumulative volumes of the MODFLOW budget in the `.list` file at the
  end of the run. In MODFLOW's sign convention, `STREAM_LEAKAGE_IN` is water entering the aquifer
  from the river and `STREAM_LEAKAGE_OUT` the reverse; `STORAGE_OUT` is water taken into storage
  and `STORAGE_IN` water released from it. So `STREAM_LEAKAGE_IN − STREAM_LEAKAGE_OUT` is the
  net river-to-aquifer volume of `ΔLeakage*`, and `STORAGE_OUT − STORAGE_IN` the net storage gain of
  `ΔStorage*`. The other budget terms (`CONSTANT_HEAD_*`, `HEAD_DEP_BOUNDS_*` for the GHB, `TOTAL_*`,
  `PERCENT_DISCREPANCY`) are exported too.

Earlier versions took `Δt_p` and `ΔH` from the time of maximum simulated depth at the second and
second-to-last reaches. Because the depth is printed to 1 mm, its maximum is often a plateau; this
gave errors of up to 46 % and made the values depend on the time step. Those values are still exported as
`delta_tp_reach`, `delta_H_reach` and `H_ratio_reach`, only to compare with the master Excel.

Also exported: `Qexc_in` and `Qexc_out` (integrals of the positive and negative net exchange),
`T_in` and `T_out` (times of the exchange maximum and minimum), `V_hid` (hydrograph volume),
`Qexc_end_ratio` (net exchange at the end relative to its peak, to check that the run was long
enough), `nstp_total`, `dt_min_s` and `duration_days`.

## References

1. Bieger, K., Rathjens, H., Allen, P. M., & Arnold, J. G. (2015). Development and evaluation of bankfull hydraulic geometry relationships for the physiographic regions of the United States. *Journal of the American Water Resources Association*, 51(3), 842-858.

2. Cooper, H.H., & Rorabaugh, M.I. (1963). Ground-water movements and bank storage due to flood stages in surface streams. *USGS Water Supply Paper*, 1536-J.

3. Zhou, Y., et al. (2024). Regional bankfull geometry relationships. *Geomorphology*, 108897.

4. Niswonger, R.G., & Prudic, D.E. (2005). Documentation of the Streamflow-Routing (SFR2) Package. *USGS Techniques and Methods*, 6-A13.

5. Chow, V.T. (1959). *Open-Channel Hydraulics*. McGraw-Hill, New York.

6. Arcement, G.J., & Schneider, V.R. (1989). Guide for selecting Manning's roughness coefficients for natural channels and flood plains. *USGS Water-Supply Paper* 2339.

7. Kenney, B.C. (1982). Beware of spurious self-correlations! *Water Resources Research*, 18(4), 1041-1048.

