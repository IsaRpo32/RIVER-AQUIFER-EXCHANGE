"""
MODFLOW Professional Dashboard - Utils Submodule (Legacy)
=========================================================

NOTE: This is legacy code preserved for reference.
This dashboard module was part of the original project but is not actively maintained.

Elite-level interactive dashboard for MODFLOW groundwater-surface water interaction analysis.
Seamlessly integrated as a utils submodule for notebook use.

Usage Examples:
--------------

1. **From Jupyter Notebook (Automatic Integration):**
   ```python
   from utils.dashboard import launch_dashboard_from_notebook

   # Launch dashboard with current notebook data
   launch_dashboard_from_notebook(globals())
   ```

2. **From Jupyter Notebook (Manual Integration):**
   ```python
   from utils.dashboard import MODFLOWDashboard

   # Create dashboard from existing objects
   dashboard = MODFLOWDashboard.from_modflow_objects(
       geology=geology,
       hydrology=hydrology,
       model_io=model_io,
       results=results
   )
   dashboard.launch()
   ```

Features:
---------
- Professional UI with Nature journal styling
- Advanced time series analysis with spectral methods
- Comprehensive river flow and exchange analysis
- Dynamic animations and 3D visualizations
- Statistical analysis with distribution testing
- Publication-ready export tools
- Interactive 3D topography visualization
- Real-time performance optimization

Technology Stack:
----------------
- Streamlit + Plotly + PyDeck for visualization
- SciPy + NumPy + Pandas for analysis
- Professional CSS styling and animations
- FloPy integration for MODFLOW data
"""

# Note: This legacy dashboard requires additional dependencies that may not be installed.
# Original imports are commented out to prevent import errors.

__version__ = "1.0.0"
__author__ = "Claude Code"
__description__ = "Professional MODFLOW Dashboard - Utils Submodule (Legacy)"

# Imports commented out - this is legacy code
# import os
# import sys
# import time
# import warnings
# from pathlib import Path
# import traceback

DASHBOARD_AVAILABLE = False
MISSING_DEPS = ['streamlit', 'plotly', 'pydeck', 'folium', 'altair']
AVAILABLE_DEPS = []

# Note: The original code performed extensive dependency checking and
# dynamic imports. This legacy version is preserved for reference only.

__all__ = [
    # Main classes (not available - legacy)
    'MODFLOWDashboard',
    'Dashboard',

    # Launch functions (not available - legacy)
    'launch_dashboard',
    'launch_dashboard_from_notebook',
    'launch_dashboard_with_data',
    'dashboard_from_notebook',
    'run_dashboard_in_notebook',

    # Utility functions (not available - legacy)
    'create_dashboard_data',
    'check_dashboard_requirements',
    'install_dashboard_requirements',
    'get_dashboard_info',
    'run_dashboard_diagnostics',
    'test_dashboard_imports',

    # Components (not available - legacy)
    'TimeSeriesAnalyzer',
    'RiverAnalyzer',
    'AnimationCreator',
    'StatisticalAnalyzer',

    # Configuration (not available - legacy)
    'THEMES',
    'ANALYSIS_DEFAULTS',
    'get_theme',
    'get_analysis_defaults',

    # Status
    'DASHBOARD_AVAILABLE'
]

print("Legacy dashboard module - not actively maintained")
print("See active code in src/ directory")
