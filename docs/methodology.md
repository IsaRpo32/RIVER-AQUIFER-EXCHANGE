# Methodology

## Overview

This document describes the methodology for modeling transient river-aquifer exchange during flood events using a combination of geomorphic scaling, dimensionless analysis, and numerical groundwater modeling.

## 1. Geomorphic Scaling

### 1.1 Regional Curves

Channel geometry is estimated from drainage area using regional curves based on empirical relationships from the literature.

#### Bankfull Width (Bieger et al., 2015)

```
w_bkf = a × A_d^b
```

Where:
- `w_bkf`: Bankfull width (m)
- `A_d`: Drainage area (km²)
- `a, b`: Regional coefficients

#### Bankfull Depth (Bieger et al., 2015)

```
d_bkf = c × A_d^d
```

Where:
- `d_bkf`: Bankfull depth (m)
- `c, d`: Regional coefficients

#### Channel Slope

```
S = e × A_d^f
```

Where:
- `S`: Channel slope (m/m)
- `e, f`: Regional coefficients (typically f ≈ -0.35)

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
- `β`: Decay coefficient

### 2.2 Event Duration Scaling

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

## 6. Analysis Outputs

### 6.1 Exchange Flow Metrics

- **Q_in**: Stream-to-aquifer exchange (infiltration)
- **Q_out**: Aquifer-to-stream exchange (exfiltration)
- **Net exchange**: Q_in - Q_out

### 6.2 Integral Metrics

- **Bank storage volume**: ∫Q_in dt during rising limb
- **Return flow volume**: ∫Q_out dt during falling limb
- **Hysteresis area**: Difference between bank storage and return flow

## References

1. Bieger, K., et al. (2015). Development and evaluation of bankfull hydraulic geometry relationships. *Journal of Hydrology*, 530, 761-771.

2. Cooper, H.H., & Rorabaugh, M.I. (1963). Ground-water movements and bank storage due to flood stages in surface streams. *USGS Water Supply Paper*, 1536-J.

3. Zhou, Y., et al. (2024). Regional bankfull geometry relationships. *Geomorphology*, 108897.

4. Niswonger, R.G., & Prudic, D.E. (2005). Documentation of the Streamflow-Routing (SFR2) Package. *USGS Techniques and Methods*, 6-A13.

