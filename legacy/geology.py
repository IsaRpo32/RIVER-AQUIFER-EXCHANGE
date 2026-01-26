"""
Geology Module (Legacy)

NOTE: This is legacy code preserved for reference.
The active functionality has been moved to src/modflow_utils.py

Handles topography generation and river placement for MODFLOW simulations.
"""

import numpy as np
from .settings import ModflowConfig


class Geology(ModflowConfig):
    """
    Geology class for generating topography and river networks.

    This class handles the creation of synthetic topography and river systems
    for groundwater-surface water interaction modeling.
    """

    def __init__(self, verbose=True):
        """
        Initialize Geology with inherited configuration.

        Parameters
        ----------
        verbose : bool, optional
            Control output verbosity (default: True)
        """
        super().__init__(verbose=verbose)

        self.topography = None
        self.river_network = None
        self.river_cells = None
        self.i_riv = None
        self.j_riv = None
        self.nc_riv = None

    def generate_valley_topography(self, plot=False):
        """
        Generate synthetic valley topography.

        Creates a valley-shaped topography with longitudinal and transverse slopes.

        Parameters
        ----------
        plot : bool, optional
            Whether to plot the generated topography (default: False)

        Returns
        -------
        numpy.ndarray
            2D array representing elevation at each grid cell
        """
        if not isinstance(self.nrow, int) or not isinstance(self.ncol, int):
            raise ValueError("Grid dimensions must be integers")
        if self.nrow <= 0 or self.ncol <= 0:
            raise ValueError("Grid dimensions must be positive")

        delta_h = self.slope_y * self.L

        y_coords = np.linspace(0, self.L, self.nrow)
        x_coords = np.linspace(0, self.B, self.ncol)
        X, Y = np.meshgrid(x_coords, y_coords)

        longitudinal_profile = (self.slope_y * Y)

        valley_center = self.B / 2
        transverse_profile = self.slope_x * ((X - valley_center) ** 2)

        elevation = longitudinal_profile + transverse_profile
        elevation = elevation - np.min(elevation)
        elevation = elevation + self.min_topo

        self.topography = elevation

        if self.verbose:
            print(f"Generated valley topography: {self.nrow}x{self.ncol} grid")
            print(f"Elevation range: {elevation.min():.1f} - {elevation.max():.1f} m")

        return self.topography

    def generate_straight_river(self, orientation='v', width=1, start_position=None):
        """
        Generate a straight river along the valley.

        Parameters
        ----------
        orientation : str, optional
            River orientation: 'v' for vertical, 'h' for horizontal (default: 'v')
        width : int, optional
            River width in grid cells (default: 1)
        start_position : int, optional
            Starting position perpendicular to flow direction (default: center)

        Returns
        -------
        numpy.ndarray
            2D binary array where 1 indicates river cells
        """
        if orientation not in ['v', 'h']:
            raise ValueError("Orientation must be 'v' (vertical) or 'h' (horizontal)")

        river_grid = np.zeros((self.nrow, self.ncol), dtype=int)

        if start_position is None:
            start_position = (self.ncol // 2) if orientation == 'v' else (self.nrow // 2)

        if orientation == 'v':
            for i in range(self.nrow):
                width_start = max(0, start_position - width // 2)
                width_end = min(self.ncol, start_position + width // 2 + 1)

                for j in range(width_start, width_end):
                    river_grid[i, j] = 1
        else:
            for j in range(self.ncol):
                width_start = max(0, start_position - width // 2)
                width_end = min(self.nrow, start_position + width // 2 + 1)

                for i in range(width_start, width_end):
                    river_grid[i, j] = 1

        self.river_network = river_grid
        self.i_riv, self.j_riv = np.where(river_grid == 1)
        self.i_riv, self.j_riv = self._adjust_river_coordinates(self.i_riv, self.j_riv, orientation)
        self.nc_riv = np.sum(river_grid)

        if self.verbose:
            print(f"Generated straight river: {orientation} orientation, width={width}")
            print(f"River cells: {self.nc_riv}")

        return self.river_network

    def _adjust_river_coordinates(self, i_riv, j_riv, orientation):
        """
        Adjust river coordinates to ensure proper flow direction.
        """
        if orientation == 'h':
            sort_idx = np.lexsort((i_riv, j_riv))
            i_riv, j_riv = j_riv[sort_idx], i_riv[sort_idx]

        unique_values, unique_indices, unique_counts = np.unique(i_riv, return_counts=True, return_index=True)

        duplicate_starts = np.ones(len(i_riv), dtype=bool)
        if len(unique_indices) > 1:
            duplicate_starts[unique_indices[1:] - 1] = False

        for i in range(len(duplicate_starts) - 1):
            if duplicate_starts[i] and duplicate_starts[i + 1]:
                duplicate_starts[i + 1] = False

        start_indices = np.where(duplicate_starts)[0]
        counts_with_duplicates = unique_counts[unique_counts != 1]

        for idx, start_pos in enumerate(start_indices):
            if idx < len(counts_with_duplicates):
                end_pos = start_pos + counts_with_duplicates[idx]

                if start_pos > 0:
                    before_val = j_riv[start_pos - 1]
                    current_val = j_riv[start_pos]

                    if before_val > current_val:
                        j_riv[start_pos:end_pos] = j_riv[start_pos:end_pos][::-1]

        if orientation == 'h':
            i_riv, j_riv = j_riv, i_riv

        return i_riv, j_riv

    def get_river_elevations(self):
        """
        Get elevations along the river network.
        """
        if self.topography is None or self.river_network is None:
            raise ValueError("Generate topography and river network first")

        river_elevations = self.topography[self.i_riv, self.j_riv]
        return river_elevations
