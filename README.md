# River-Aquifer Exchange Modeling

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A Python framework for modeling groundwater-surface water interactions using MODFLOW-NWT and dimensionless analysis. This repository contains tools for geomorphic scaling, dimensionless parameter computation, and Monte Carlo simulation of river-aquifer exchange dynamics.

## Overview

This project implements methodologies for analyzing transient river-aquifer exchange during flood events. It combines:

- **Geomorphic scaling relationships** based on regional curves (Bieger et al., 2015; Zhou et al., 2024)
- **Dimensionless Pi numbers** for hydraulic parameter analysis
- **MODFLOW-NWT** groundwater modeling with Streamflow Routing (SFR2) package
- **Monte Carlo simulation** for parameter space exploration

## Project Structure

```
river-aquifer-exchange/
├── README.md                           # This file
├── LICENSE                             # MIT License
├── requirements.txt                    # Python dependencies
├── .gitignore                          # Git ignore rules
│
├── src/                                # Active source code
│   ├── __init__.py
│   ├── dimensionless_numbers.py        # Dimensionless Pi number calculations
│   ├── geomorphic_scaling.py           # Geomorphic scaling relationships
│   └── modflow_utils.py                # MODFLOW utilities and helpers
│   └── signal_processing.py                # MODFLOW utilities and helpers
│
├── simulations/                        # Simulation scripts and notebooks
│   ├── parallel_modflow.py             # Parallel MODFLOW execution
│   └── montecarlo_sampling.ipynb       # Monte Carlo parameter sampling
│
├── executables/                        # MODFLOW executables (not tracked)
│   └── README.md                       # Instructions for obtaining executables
│
├── data/                               # Data directories
│   ├── input/                          # Input data files
│   └── output/                         # Simulation outputs
│
├── docs/                               # Documentation
│   └── methodology.md                  # Detailed methodology
│
└── legacy/                             # Archived legacy code (reference only)
    ├── README.md
    └── ...
```

## Installation

### Prerequisites

- Python 3.8 or higher
- MODFLOW-NWT executable (see `executables/README.md`)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/river-aquifer-exchange.git
   cd river-aquifer-exchange
   ```

2. Create a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Download MODFLOW-NWT executable and place in `executables/` directory.

## Usage

### Dimensionless Numbers

Calculate dimensionless parameters for hydraulic analysis:

```python
from src import dimensionless_numbers as dn

# Calculate characteristic time
tc = dn.get_tc(Sy=0.15, Ka=1e-5, B=500, d_bkf=2.5)

# Calculate dimensionless gamma (event duration ratio)
gamma_d = dn.get_gamma_d(tc, t_d=86400*5)

# Calculate lambda_v (vertical flux parameter)
lambda_v = dn.get_lambda_v(tc, Kr=1e-6, dr=0.2)
```

### Geomorphic Scaling

Estimate channel geometry from drainage area:

```python
from src import geomorphic_scaling as gs

# Get bankfull dimensions for a 1000 km² watershed
da_km2 = 1000
w_bkf = gs.get_w_bkf(da_km2)  # Bankfull width (m)
d_bkf = gs.get_d_bkf(da_km2)  # Bankfull depth (m)
slope = gs.get_S(da_km2)       # Channel slope (m/m)
```

### Running Simulations

For parallel MODFLOW simulations:

```python
from simulations.parallel_modflow import run_simulation

# Configure and run simulation
results = run_simulation(
    drainage_area=1000,
    event_duration=5,  # days
    n_workers=4
)
```

## Key Modules

### `src/dimensionless_numbers.py`

Calculates dimensionless Pi numbers for river-aquifer exchange analysis:

| Parameter | Symbol | Pi Number | Description |
|-----------|--------|-----------|-------------|
| `gamma_d` | γ_d | Π₁ | Dimensionless event duration |
| `lambda_v` | λ_v | Π₂ | Vertical exchange parameter |
| `theta_r` | θ_r | Π₃ | River width ratio |
| `chi` | χ | Π₄ | Flow resistance parameter |
| `rho_0` | ρ_0 | Π₅ | Initial head ratio |
| `zeta` | ζ | Π₆ | Longitudinal slope parameter |

### `src/geomorphic_scaling.py`

Regional geomorphic relationships based on:

- **Bieger et al. (2015)**: Bankfull width and depth scaling
- **Zhou et al. (2024)**: Channel geometry relationships
- **Hack's Law**: Hydraulic length scaling

### `src/modflow_utils.py`

Utilities for MODFLOW model construction:

- Cooper-Rorabaugh (1963) hydrograph generation
- Topography and river network creation
- SFR2 package configuration
- Results extraction and analysis

## Scientific Background

### Cooper-Rorabaugh Hydrograph

The flood hydrograph follows the Cooper and Rorabaugh (1963) analytical solution:

```
h(t) = H_0 + H_p × exp(-β(t - t_p)) × (1 - cos(ωt)) / (1 - cos(ωt_p))
```

Where:
- `H_0`: Base stage
- `H_p`: Peak stage increment
- `t_p`: Time to peak
- `ω`: Angular frequency (2π/t_d)

### Dimensionless Analysis

The river-aquifer exchange is characterized by dimensionless groups:

1. **Π₁ (γ_d) = t_c / t_d**: Ratio of aquifer response time to event duration
2. **Π₂ (λ_v) = K_r × d_r / (t_c × K_a)**: Vertical exchange intensity
3. **Π₃ (θ_r) = S_y × w_bkf / B**: Effective river width ratio

## References

- Bieger, K., et al. (2015). Development and evaluation of bankfull hydraulic geometry relationships. *Journal of Hydrology*, 530, 761-771.
- Cooper, H.H., & Rorabaugh, M.I. (1963). Ground-water movements and bank storage due to flood stages in surface streams. *USGS Water Supply Paper*, 1536-J.
- Zhou, Y., et al. (2024). Regional bankfull geometry relationships. *Geomorphology*, 108897.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{river_aquifer_exchange,
  author = {Isabel Restrepo},
  title = {River-Aquifer Exchange Modeling Framework},
  year = {2025},
  url = {https://github.com/IsaRpo32/RIVER-AQUIFER-EXCHANGE.git}
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

For questions or collaboration, please open an issue or contact the maintainer.

