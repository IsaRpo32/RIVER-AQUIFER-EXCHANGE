"""
Legacy Code Module

NOTE: This folder contains legacy code preserved for reference purposes only.
This code is NOT actively maintained and may not function with current dependencies.

The active functionality has been moved to the src/ directory:
- src/dimensionless_numbers.py - Dimensionless Pi numbers (from phi_numbers2.py)
- src/geomorphic_scaling.py - Geomorphic scaling functions (from scalar_functions_v5.py)
- src/modflow_utils.py - MODFLOW utilities (from utils/inputs.py)

Legacy modules preserved here:
- modflow_io.py - MODFLOW model setup class
- settings.py - Configuration settings
- geology.py - Topography and river generation
- hydrology.py - Hydrograph generation
- results.py - Results analysis
- visualization.py - Plotting functions
- signal_processing.py - Signal filtering utilities
- dashboard/ - Interactive dashboard components
"""

__all__ = [
    'ModflowIO',
    'ModflowConfig',
    'Geology',
    'Hydrology',
    'Results',
    'Visualization',
    'SignalProcessor',
]

# Note: Imports are commented out as this legacy code may have dependencies
# that are no longer available or compatible.
#
# from .modflow_io import ModflowIO
# from .settings import ModflowConfig
# from .geology import Geology
# from .hydrology import Hydrology
# from .results import Results
# from .visualization import Visualization
# from .signal_processing import SignalProcessor
