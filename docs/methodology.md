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

To analyze the system, we purpose something pi numbers that are dimensionless number that allow understand the system. These are explained in paper.
## 4. Numerical Modeling

### 4.1 MODFLOW-NWT Configuration

The groundwater model uses MODFLOW-NWT with the following packages:

| Package | Purpose |
|---------|---------|
| DIS | Model discretization |
| BAS | Basic package (boundary conditions) |
| LPF/UPW | Layer property flow |
| NWT | Newton solver |
| SFR2 | Streamflow routing |
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

- **River (SFR2)**: Time-varying stage from hydrograph
- **Lateral boundaries**: No-flow (symmetry assumption)
- **Bottom**: No-flow (impermeable base)
- **Top**: Water table (unconfined)

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

Parameters are sampled within physically realistic ranges:

| Parameter | Range | Distribution |
|-----------|-------|--------------|
| K_a | 10⁻⁷ - 10⁻⁴ m/s | Log-uniform |
| K_r/K_a | 10⁻³ - 10⁻¹ | Log-uniform |
| S_y | 0.10 - 0.30 | Uniform |
| n (Manning) | 0.025 - 0.045 | Uniform |
| B/w_bkf | 10 - 200 | Uniform |
| D/d_bkf | 1 - 20 | Uniform |

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
input `data/input/inputs_n_corregido.csv` keeps the same 758 draws divided by 10 and stores the
original value in `n_original`. `inputs_n_original.csv` has the same rows with the original `n`,
for comparison. Both are built from the master Excel by `simulations/build_input_csv.py`.

## 6. Response Metrics

The five response metrics of the manuscript (Eq. 16–20). `parallel_modflow.py` exports the
dimensional quantities for each row; the normalization (the metrics marked with `*`) is done in
post-processing. All times are measured from the start of the transient period.

| Metric | Definition | Exported columns |
|--------|------------|------------------|
| Peak delay (Eq. 16) | `Δt_p* = (t_p,out − t_p,in) / t_c` | `delta_tp`, `t_p_in`, `t_p_out` |
| Stage attenuation (Eq. 17) | `ΔH* = (H_p,in − H_p,out) / d_bkf` | `delta_H`, `H_p_in`, `H_p_out`, `H_ratio` |
| Flood-exchange efficiency (Eq. 18) | `ΔLeakage* = ∫(Q_exc,in − Q_exc,out) dt / (B·L·d_bkf)` | `STREAM_LEAKAGE_IN`, `STREAM_LEAKAGE_OUT` |
| Return time (Eq. 19) | `t_r* = t_r / t_c` | `t_r` |
| Storage change (Eq. 20) | `ΔStorage* = (V_storage,out − V_storage,in) / (B·L·d_bkf)` | `STORAGE_OUT`, `STORAGE_IN` |

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
  net river-to-aquifer volume of Eq. 18, and `STORAGE_OUT − STORAGE_IN` the net storage gain of
  Eq. 20. The other budget terms (`CONSTANT_HEAD_*`, `HEAD_DEP_BOUNDS_*` for the GHB, `TOTAL_*`,
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

