"""
MODFLOW Configuration Settings (Legacy)

NOTE: This is legacy code preserved for reference.
All measurements are in SI system (m, s, m3/s, m2, etc)
"""

import os
import numpy as np


class ModflowConfig:
    """
    Base configuration class for MODFLOW simulations.

    This class contains all the fundamental parameters needed for groundwater-surface water
    interaction modeling using MODFLOW with SFR2 package.

    All units are in SI system unless otherwise specified.
    """

    def __init__(self, verbose=True):
        """
        Initialize MODFLOW configuration with default parameters.

        Parameters
        ----------
        verbose : bool, optional
            Control output verbosity (default: True)
        """
        self.verbose = verbose

        # MODFLOW EXECUTION SETTINGS
        self.path_results = "../results/simple_model/"
        self.name_model = 'groundwater_surface_interaction'
        self.version = "mfnwt"
        self.exe_name = "MODFLOW-NWT.exe"

        # UNITS
        self.itmuni = 1  # Time unit: 1=seconds
        self.lenuni = 2  # Length unit: 2=meters

        # SPATIAL DISCRETIZATION
        self.L = 5000       # Domain length (m)
        self.B = 5000       # Domain width (m)
        self.D = 35         # Aquifer thickness (m)

        self.nrow = 20
        self.ncol = 11
        self.nlay = 10

        self.delr = self.B / self.ncol
        self.delc = self.L / self.nrow

        self.slope_x = 1e-12
        self.slope_y = 0.0001

        # AQUIFER PROPERTIES
        self.ka = 100/86400         # Horizontal hydraulic conductivity (m/s)
        self.kr = 0.1/86400         # Streambed hydraulic conductivity (m/s)
        self.sy = 0.15              # Specific yield
        self.ss = 1e-7              # Specific storage (1/m)
        self.kv_kh = 1.0            # Vertical to horizontal K ratio

        # TEMPORAL DISCRETIZATION
        self.cte_courant = .5
        self.perlen = [1*24*3600, 20*24*3600]
        self.nstp = [1, 3000]
        self.steady = [True, False]

        # STREAM PROPERTIES
        self.d_bkf = 3              # Bankfull depth (m)
        self.w_bkf = 50             # Bankfull width (m)
        self.dr = 1                 # Streambed thickness (m)
        self.n_manning = 0.03       # Manning's roughness

        # HYDROGRAPH PARAMETERS
        self.fname = "discharge_hydrograph.dat"
        self.H_0 = 1
        self.eta_0 = (self.d_bkf - self.H_0 - 1) / self.H_0
        self.t_d = 2.5 * 24 * 3600
        self.eta_d = 0.5
        self.t_delay = 1 * 24 * 3600
        self.H_p = self.eta_0 * self.H_0

        # SOLVER PARAMETERS
        self.maxiterout = 20000
        self.fluxtol = 1000
        self.headtol = 1e-6
        self.numtim = 50
        self.weight = 0.8
        self.flwtol = 1e-5

        # BOUNDARY CONDITIONS
        self.min_topo = self.D + self.d_bkf + 1

    def _get_transmivity(self):
        return self.ka * self.D
