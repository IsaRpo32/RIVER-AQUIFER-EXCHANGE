"""
MODFLOW Input/Output and Model Setup Module (Legacy)

NOTE: This is legacy code preserved for reference.
The active functionality has been moved to src/modflow_utils.py

This module handles MODFLOW model creation, configuration, and execution.
It inherits from ModflowConfig for base parameters and Results for analysis.
"""

import os
import time
import flopy
import numpy as np
from .settings import ModflowConfig
from .results import Results
from .geology import Geology


class ModflowIO(ModflowConfig):
    """
    MODFLOW model setup and execution class.

    This class handles the complete MODFLOW model workflow including:
    - Model initialization and configuration
    - Package setup (DIS, UPW, SFR, BAS, OC, NWT)
    - Model execution and results handling

    Inherits from:
    - ModflowConfig: Base configuration parameters

    Composition with:
    - Results: Results analysis and visualization (accessed via self.results)
    """

    def __init__(self, verbose=True):
        """
        Initialize MODFLOW model setup.

        Parameters
        ----------
        verbose : bool, optional
            Control output verbosity (default: True)
        """
        super().__init__()
        self.results = Results()
        self.verbose = verbose
        self.model = None
        self.success = False
        self.elapsed_time = None
        self.i_riv = None
        self.j_riv = None
        self.Q = None

        if self.verbose:
            print("ModflowIO initialized successfully")

    def _setup_modflow_model(self):
        """
        Initialize the MODFLOW model object.

        Creates the basic MODFLOW model with specified parameters.
        """
        os.makedirs(self.path_results, exist_ok=True)

        self.model = flopy.modflow.Modflow(
            self.name_model,
            exe_name=self.exe_name,
            version=self.version,
            model_ws=self.path_results,
            verbose=self.verbose
        )

        if self.verbose:
            print(f"MODFLOW model '{self.name_model}' created")
            print(f"   Workspace: {self.path_results}")
            print(f"   Executable: {self.exe_name}")

    def setup_discretization(self, topography, botm):
        """
        Setup spatial and temporal discretization package (DIS).

        Parameters
        ----------
        topography : numpy.ndarray
            Top elevation array (nrow, ncol)
        botm : numpy.ndarray or list
            Bottom elevation array(s) for each layer
        """
        if topography.shape != (self.nrow, self.ncol):
            raise ValueError(f"Topography shape {topography.shape} doesn't match grid ({self.nrow}, {self.ncol})")

        self.topography = topography

        self.dis = flopy.modflow.ModflowDis(
            self.model,
            nlay=self.nlay,
            nrow=self.nrow,
            ncol=self.ncol,
            delr=self.delr,
            delc=self.delc,
            top=topography,
            botm=botm,
            itmuni=self.itmuni,
            lenuni=self.lenuni,
            nper=len(self.nstp),
            perlen=self.perlen,
            nstp=self.nstp,
            steady=self.steady
        )

        if self.verbose:
            print(f"Discretization package setup complete")
            print(f"   Grid: {self.nrow} x {self.ncol} x {self.nlay}")
            print(f"   Periods: {len(self.perlen)} ({sum(self.perlen)/86400:.1f} days total)")

    def setup_external_file(self, discharge_data=None):
        """
        Setup external file for SFR discharge data.

        Parameters
        ----------
        discharge_data : numpy.ndarray, optional
            Discharge time series data. If not provided, uses self.Q
        """
        if discharge_data is not None:
            self.Q = discharge_data

        if not hasattr(self, 'Q') or self.Q is None:
            raise ValueError("Discharge data (Q) must be provided before setting up external file")

        self.model.add_external(
            fname=self.fname,
            unit=70
        )

        self.tabfiles_dict = {
            1: {
                "numval": len(self.Q),
                "inuit": 70
            }
        }

        if self.verbose:
            print(f"External file '{self.fname}' configured")
            print(f"   Discharge points: {len(self.Q)}")

    def setup_flow_package(self):
        """
        Setup Upstream Weighting (UPW) package for flow properties.

        Configures aquifer hydraulic properties including conductivity,
        storage parameters, and layer types.
        """
        self.upw = flopy.modflow.ModflowUpw(
            self.model,
            hk=self.nlay * [self.ka],
            vka=self.nlay * [self.kv_kh * self.ka],
            sy=self.nlay * [self.sy],
            ss=self.nlay * [self.ss],
            laytyp=self.nlay * [1],
            hdry=-1e30,
            ipakcb=53,
            iphdry=0
        )

        if self.verbose:
            print(f"Flow package (UPW) setup complete")
            print(f"   Ka: {self.ka:.2e} m/s, Kv/Kh: {self.kv_kh}")
            print(f"   Sy: {self.sy}, Ss: {self.ss:.1e} 1/m}")
            print(f"   All layers: Unconfined (laytyp=1)")

    def setup_newton_solver(self):
        """
        Setup Newton-Raphson solver (NWT) package.
        """
        self.nwt = flopy.modflow.ModflowNwt(
            self.model,
            iprnwt=1,
            maxiterout=self.maxiterout,
            fluxtol=self.fluxtol,
            headtol=self.headtol
        )

        if self.verbose:
            print(f"Newton solver (NWT) setup complete")
            print(f"   Max iterations: {self.maxiterout}")
            print(f"   Tolerances - Flux: {self.fluxtol}, Head: {self.headtol}")

    def run_model(self, silent=False):
        """
        Execute MODFLOW simulation.

        Parameters
        ----------
        silent : bool, optional
            Run model silently (default: False)

        Returns
        -------
        bool
            True if model ran successfully
        """
        if self.model is None:
            raise ValueError("Model must be setup before running")

        try:
            if self.verbose:
                print("Running MODFLOW simulation...")

            start_time = time.time()

            self.model.write_input()
            self.success, _ = self.model.run_model(silent=silent, report=not silent)

            end_time = time.time()
            self.elapsed_time = end_time - start_time

            if self.verbose:
                if self.success:
                    print(f"Model completed successfully in {self.elapsed_time:.1f} seconds")
                else:
                    print(f"Model failed after {self.elapsed_time:.1f} seconds")

            return self.success

        except Exception as e:
            print(f"Error running model: {e}")
            self.success = False
            return False
