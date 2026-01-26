# Legacy Code

This folder contains legacy code preserved for reference purposes only.

## Important Notes

- This code is **NOT actively maintained**
- It may not function with current dependencies
- It is preserved for historical reference and potential future use

## Active Code Location

The active functionality has been migrated to the `src/` directory:

| Legacy File | Active Equivalent |
|-------------|-------------------|
| `utils/inputs.py` | `src/modflow_utils.py` |
| `phi_numbers2.py` | `src/dimensionless_numbers.py` |
| `scalar_functions_v5.py` | `src/geomorphic_scaling.py` |

## Contents

### Core Modules

- **modflow_io.py**: MODFLOW model setup and execution class
- **settings.py**: Base configuration parameters for simulations
- **geology.py**: Topography generation and river placement
- **hydrology.py**: Hydrograph generation (Cooper-Rorabaugh method)
- **results.py**: Results analysis and SFR data extraction
- **visualization.py**: Comprehensive plotting functions

### Utility Modules

- **signal_processing.py**: Signal smoothing and filtering utilities

### Dashboard (subdirectory)

Interactive visualization dashboard components:
- **core.py**: Main dashboard application
- **components.py**: Reusable UI components
- **config.py**: Dashboard configuration
- **views.py**: Dashboard view definitions

## Original Class Hierarchy

```
ModflowConfig (settings.py)
    ├── Geology (geology.py)
    │   └── Hydrology (hydrology.py)
    │       └── Results (results.py)
    ├── Visualization (visualization.py)
    └── ModflowIO (modflow_io.py)
```

## Why Legacy?

This code was part of an earlier object-oriented design that proved overly complex for the research workflow. The current approach in `src/` uses simpler, standalone functions that are:

1. Easier to test and debug
2. More flexible for parallel simulations
3. Better suited for Monte Carlo sampling
4. Compatible with multiprocessing

The legacy code is preserved in case specific OOP patterns are needed in the future.
