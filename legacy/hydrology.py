"""
Hydrology Module (Legacy)

NOTE: This is legacy code preserved for reference.
The active hydrograph functionality has been moved to src/modflow_utils.py

Handles hydrograph generation and flow calculations for MODFLOW simulations.
"""

import os
import numpy as np
from .settings import ModflowConfig
from .signal_processing import SignalProcessor


class Hydrology(ModflowConfig, SignalProcessor):
    """
    Hydrology class for generating hydrographs and flow time series.

    This class implements various hydrograph models for use in groundwater-surface
    water interaction modeling with signal processing capabilities.
    """

    def __init__(self, verbose=True):
        """
        Initialize Hydrology with inherited configuration.

        Parameters
        ----------
        verbose : bool, optional
            Control output verbosity (default: True)
        """
        super().__init__(verbose=verbose)

        self.time_array = None
        self.stage_hydrograph = None
        self.discharge_hydrograph = None
        self.velocity_series = None
        self.hydrograph_file_path = None
        self.filtered_discharge = None

    def _generate_time_array(self):
        """Generate time array for hydrograph based on configuration."""
        time_resolution = 600
        self.time_array = np.linspace(0, self.t_d, time_resolution)

        if self.verbose:
            print(f"Generated time array: {len(self.time_array)} points over {self.t_d/86400:.2f} days")

    def generate_cooper_rorabaugh_hydrograph(self):
        """
        Generate deterministic hydrograph using Cooper and Rorabaugh (1963) method.

        Returns
        -------
        tuple
            (time_array, stage_hydrograph) in seconds and meters
        """
        self._generate_time_array()

        H_p = self.eta_0 * self.H_0
        t_p = self.eta_d * self.t_d
        omega = 2 * np.pi / self.t_d

        angle = omega * t_p / 2
        if np.abs(np.sin(angle)) < 1e-10:
            beta = omega / (angle + 1e-10)
        else:
            beta = omega / np.tan(angle)

        t = np.asarray(self.time_array)
        cos_term = (1 - np.cos(omega * t)) / (1 - np.cos(omega * t_p))
        exp_term = np.exp(-beta * (t - t_p))

        self.stage_hydrograph = np.where(
            (t >= 0) & (t <= self.t_d),
            self.H_0 + H_p * exp_term * cos_term,
            self.H_0
        )

        if self.verbose:
            print(f"Generated Cooper-Rorabaugh hydrograph:")
            print(f"  Peak stage: {np.max(self.stage_hydrograph):.3f} m")
            print(f"  Time to peak: {t_p/3600:.2f} hours")
            print(f"  Event duration: {self.t_d/86400:.2f} days")

        return self.time_array, self.stage_hydrograph

    def stage_to_discharge(self, use_manning=True):
        """
        Convert stage hydrograph to discharge using Manning's equation.

        Parameters
        ----------
        use_manning : bool, optional
            Use Manning's equation (True) or simple rating curve (False)

        Returns
        -------
        numpy.ndarray
            Discharge time series (m3/s)
        """
        if self.stage_hydrograph is None:
            raise ValueError("Generate stage hydrograph first")

        if use_manning:
            self.discharge_hydrograph = (
                (1 / self.n_manning) *
                self.w_bkf *
                (self.stage_hydrograph ** (5/3)) *
                (self.slope_y ** 0.5)
            )
        else:
            a = 1.0
            b = 2.5
            h0 = 0.1

            depth_above_minimum = np.maximum(self.stage_hydrograph - h0, 0)
            self.discharge_hydrograph = a * (depth_above_minimum ** b)

        if self.verbose:
            method = "Manning's equation" if use_manning else "rating curve"
            print(f"Converted stage to discharge using {method}")
            print(f"  Peak discharge: {np.max(self.discharge_hydrograph):.3f} m3/s")
            print(f"  Base discharge: {np.min(self.discharge_hydrograph):.3f} m3/s")

        return self.discharge_hydrograph

    def calculate_flow_velocities(self):
        """
        Calculate flow velocities from discharge and channel geometry.

        Returns
        -------
        numpy.ndarray
            Velocity time series (m/s)
        """
        if self.discharge_hydrograph is None:
            raise ValueError("Generate discharge hydrograph first")

        cross_sectional_area = self.w_bkf * self.stage_hydrograph
        cross_sectional_area = np.maximum(cross_sectional_area, 1e-6)
        self.velocity_series = self.discharge_hydrograph / cross_sectional_area

        if self.verbose:
            print(f"Calculated flow velocities:")
            print(f"  Maximum velocity: {np.max(self.velocity_series):.3f} m/s")
            print(f"  Mean velocity: {np.mean(self.velocity_series):.3f} m/s")

        return self.velocity_series

    def save_hydrograph_file(self, include_delay=True):
        """
        Save hydrograph to file for MODFLOW input.

        Parameters
        ----------
        include_delay : bool, optional
            Include delay time before event starts (default: True)

        Returns
        -------
        str
            Path to saved hydrograph file
        """
        if self.discharge_hydrograph is None:
            raise ValueError("Generate discharge hydrograph first")

        if include_delay:
            time_with_delay = self.time_array + self.perlen[0] + self.t_delay
        else:
            time_with_delay = self.time_array + self.perlen[0]

        hydrograph_data = np.column_stack([time_with_delay, self.discharge_hydrograph])

        os.makedirs(self.path_results, exist_ok=True)
        self.hydrograph_file_path = os.path.join(self.path_results, self.fname)

        np.savetxt(self.hydrograph_file_path, hydrograph_data,
                  fmt="%.4f", header="Time(s) Discharge(m3/s)", comments='# ')

        if self.verbose:
            print(f"Saved hydrograph to: {self.hydrograph_file_path}")
            print(f"  Data points: {len(hydrograph_data)}")
            if include_delay:
                print(f"  Includes {self.t_delay/86400:.2f} day delay")

        return self.hydrograph_file_path
