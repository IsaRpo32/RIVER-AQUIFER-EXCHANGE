"""
MODFLOW Utilities for River-Aquifer Exchange Modeling

This module provides utility functions for setting up and analyzing MODFLOW
simulations of river-aquifer interaction, including topography generation,
river geometry, hydrograph generation, and results analysis.

References:
    - Cooper, H. H., & Rorabaugh, M. I. (1963). Ground-water movements and bank
      storage due to flood stages in surface streams.
    - Kinoshita, R. (1961). Investigation of channel deformation in Ishikari River.
    - Abad, J. D., & Garcia, M. H. (2009). Experiments in a high-amplitude
      Kinoshita meandering channel.
"""

import os
import warnings

import flopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import ipywidgets as widgets
    from IPython.display import display
    JUPYTER_AVAILABLE = True
except ImportError:
    JUPYTER_AVAILABLE = False

from flopy.utils.sfroutputfile import SfrFile


# =============================================================================
# Topography Generation
# =============================================================================

def generate_topo(nrow, ncol, L, delta_h, slope_t, min_topo=20, plot=False):
    """
    Generate a simulated valley topography for MODFLOW model.

    Creates a topographic surface with longitudinal slope and a V-shaped
    transversal valley profile.

    Parameters
    ----------
    nrow : int
        Number of rows in the mesh (longitudinal direction).
    ncol : int
        Number of columns in the mesh (transversal direction).
    L : float
        Total physical length of the model in the longitudinal direction (m).
    delta_h : float
        Total elevation difference along L (m).
    slope_t : float
        Slope in the transversal direction. Controls valley shape.
    min_topo : float, optional
        Minimum topography elevation (default: 20 m).
    plot : bool, optional
        If True, displays a plot with topography and profiles (default: False).

    Returns
    -------
    numpy.ndarray
        2D array representing the generated topography.

    Raises
    ------
    ValueError
        If nrow or ncol are not positive integers.
    """
    if not (isinstance(nrow, int) and isinstance(ncol, int)):
        raise ValueError("nrow and ncol must be integers")
    if nrow <= 0 or ncol <= 0:
        raise ValueError("nrow and ncol must be greater than zero")

    # Calculate physical distance between points in longitudinal direction
    dy = L / (nrow - 1)

    # Calculate longitudinal slope from total elevation change
    slope_l = delta_h / L

    # Create physical coordinate meshes
    y_phys = np.linspace(0, L, nrow)
    x_phys = np.linspace(0, ncol, ncol)
    x, y = np.meshgrid(x_phys, y_phys)

    # Longitudinal profile using physical y coordinate
    longitudinal_profile = slope_l * y

    # Transversal profile: inverted V-shaped valley centered in x direction
    transversal_profile = slope_t * (np.abs(x - ncol / 2) ** 2)

    # Combine profiles
    z = longitudinal_profile - transversal_profile

    # Normalization and topography adjustment
    topo = -z - np.min(-z)
    topo += min_topo - np.min(topo)

    if plot:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Topography map
        im = axes[0].imshow(topo, cmap="terrain", aspect="auto",
                           extent=[x_phys.min(), x_phys.max(),
                                  y_phys.max(), y_phys.min()])
        axes[0].axvline(ncol / 2, color='r', lw=2)
        axes[0].axhline(L / 2, color='k', lw=2)
        axes[0].set_title("Valley Topography")
        axes[0].set_xlabel("X (Transversal Direction)")
        axes[0].set_ylabel("Y (Longitudinal Direction)")
        fig.colorbar(im, ax=axes[0])

        # Longitudinal profile along central row
        axes[1].plot(topo[int(nrow / 2), :], color='k')
        axes[1].set_title("Longitudinal Profile (Central Row)")
        axes[1].set_xlabel("X (Transversal Direction)")
        axes[1].set_ylabel("Elevation (m)")
        axes[1].grid(True)

        # Transversal profile along central column
        axes[2].plot(topo[:, int(ncol / 2)], color='r')
        axes[2].set_title("Transversal Profile (Central Column)")
        axes[2].set_xlabel("Y (Longitudinal Direction)")
        axes[2].set_ylabel("Elevation (m)")
        axes[2].grid(True)

        plt.tight_layout()
        plt.show()

    return topo


# =============================================================================
# River Generation
# =============================================================================

def generate_river(nrow, ncol, A, B, C, start_cell, orientation="h",
                   width_riv=None, plot=False,
                   plot_settings=dict(figsize=(10, 6), fontsize=8)):
    """
    Generate a sinuous river in a landscape with a cosine function trajectory.

    Parameters
    ----------
    nrow : int
        Number of rows in the mesh (Y direction dimension).
    ncol : int
        Number of columns in the mesh (X direction dimension).
    A : float
        Amplitude of the cosine function (defines river variability).
    B : float
        Frequency of the cosine function (controls oscillation density).
    C : float
        Initial phase of the cosine function (horizontal curve offset).
    start_cell : int
        Initial cell where the river begins.
        For orientation "h", it's the initial cell in Y.
        For orientation "v", it's the initial cell in X.
    orientation : str, optional
        River orientation: "h" (horizontal) or "v" (vertical). Default: "h".
    width_riv : int, optional
        River width in adjacent cells. If None, river has 1-cell width.
    plot : bool, optional
        If True, displays a plot of the generated river.
    plot_settings : dict
        Plot settings with figsize and fontsize.

    Returns
    -------
    numpy.ndarray
        2D array with river cells marked as 1 and rest as 0.

    Raises
    ------
    ValueError
        If orientation is not 'h' or 'v'.
    """
    # Validate orientation
    if orientation not in ["h", "v"]:
        raise ValueError("Orientation must be 'horizontal' (h) or 'vertical' (v).")

    # Create river trajectory
    if orientation == "v":
        x = np.arange(0, nrow)
        y = A * np.cos(B * x + C) + start_cell
        y = np.clip(np.round(y).astype(int), 0, ncol - 1)  # Ensure within bounds
        river_path = [(xi, yi) for xi, yi in zip(x, y)]
    elif orientation == "h":
        y = np.arange(0, ncol)
        x = A * np.cos(B * y + C) + start_cell
        x = np.clip(np.round(x).astype(int), 0, nrow - 1)  # Ensure within bounds
        river_path = [(xi, yi) for xi, yi in zip(x, y)]

    # Create empty landscape
    landscape = np.zeros((nrow, ncol), dtype=int)  # 0 indicates unoccupied area

    # Index for river flow
    river_indices = []

    # Add river to landscape
    for x, y in river_path:
        if orientation == "v":
            y_min = max(0, y - width_riv // 2 if width_riv else 0)
            y_max = min(ncol, y + width_riv // 2 + 1 if width_riv else 1)
            for yi in range(y_min, y_max):
                landscape[x, yi] = 1  # Assign value 1 to river cells

        elif orientation == "h":
            x_min = max(0, x - width_riv // 2 if width_riv else 0)
            x_max = min(nrow, x + width_riv // 2 + 1 if width_riv else 1)
            for xi in range(x_min, x_max):
                landscape[xi, y] = 1  # Assign value 1 to river cells

    if plot:
        plt.figure(figsize=plot_settings.get('figsize', (10, 6)))
        plt.imshow(landscape, cmap="Blues", origin="upper", interpolation="nearest")
        plt.title(f"River with A={A}, B={B}, C={C}, basic")
        plt.xlabel("X (Longitudinal Direction)")
        plt.ylabel("Y (Transversal Direction)")
        plt.colorbar(label="River (1 = Yes, 0 = No)")

        fontsize = plot_settings.get('fontsize', 8)
        for x, y in river_path:
            if orientation == "h":
                y_min = max(0, y - width_riv // 2 if width_riv else 0)
                y_max = min(ncol, y + width_riv // 2 + 1 if width_riv else 1)
                for yi in range(y_min, y_max):
                    plt.text(yi, x, f'({x},{yi})', fontsize=fontsize,
                            color='red', ha='center', va='center')
            elif orientation == "v":
                x_min = max(0, x - width_riv // 2 if width_riv else 0)
                x_max = min(nrow, x + width_riv // 2 + 1 if width_riv else 1)
                for xi in range(x_min, x_max):
                    plt.text(y, xi, f'({xi},{y})', fontsize=fontsize,
                            color='red', ha='center', va='center')

    return landscape


def kinoshita_curve(theta_0, lambda_value, j_s, j_f, n, m_points=1000, ds=None):
    """
    Generate a Kinoshita Curve for meandering river simulation.

    The Kinoshita curve is based on Kinoshita (1961). Equations are from
    Abad and Garcia (2009).

    Parameters
    ----------
    theta_0 : float
        Maximum angular amplitude in radians.
    lambda_value : float
        Arc wavelength.
    j_s : float
        Skewness parameter.
    j_f : float
        Flatness (or "fatness") parameter.
    n : int
        Number of meander loops.
    m_points : int, optional
        Number of points describing the meander (default: 1000).
        Overwritten if ds is provided.
    ds : float, optional
        Delta of distance. If None, uses m_points for calculation.

    Returns
    -------
    tuple
        (X coordinates, Y coordinates, data dictionary)
        Data dictionary contains: 'curve', 'theta', 's', 'lmax', 'ymax', 'sinuosity'.
    """
    # Total path length
    smax = n * lambda_value
    if ds is None:
        s = np.linspace(0, smax, m_points)
    else:
        s = np.arange(0, smax + ds, ds)
    deltas = s[1]

    k = 2 * np.pi / lambda_value

    theta_rad = theta_0 * np.sin(k * s) + theta_0**3 * (
        j_s * np.cos(3 * k * s) - j_f * np.sin(3 * k * s)
    )
    theta_rad = theta_rad[:-1]
    theta = theta_rad * 180 / np.pi

    curve = (
        k * theta_0 * np.cos(k * s)
        - 3 * k * theta_0**3 * (j_s * np.sin(3 * k * s))
        + j_f * np.cos(3 * k * s)
    )

    # Generate coordinates by integrating displacements
    deltax = deltas * np.cos(theta_rad)
    deltay = deltas * np.sin(theta_rad)

    x = np.array([0] + list(np.cumsum(deltax)))
    y = np.array([0] + list(np.cumsum(deltay)))

    lmax = x[-1]
    ymax = np.max(np.abs(y))
    sinuosity = smax / lmax

    data = {
        "curve": curve,
        "theta": theta,
        "s": s,
        "lmax": lmax,
        "ymax": ymax,
        "sinuosity": sinuosity,
    }

    return x, y, data


def generate_river_kinoshita(nrow, ncol, theta_0, lambda_value, j_s, j_f, n,
                             m_points=1000, ds=None, width_riv=1, cell_start=(0, 0),
                             base_matrix=None, orientation="h", plot=False):
    """
    Generate a river based on the Kinoshita curve and intersect with a base matrix.

    Parameters
    ----------
    nrow : int
        Number of rows in the matrix.
    ncol : int
        Number of columns in the matrix.
    theta_0, lambda_value, j_s, j_f, n, m_points, ds : various
        Parameters for Kinoshita curve generation.
    width_riv : int, optional
        River width in cells (default: 1).
    cell_start : tuple, optional
        Starting cell coordinates (default: (0, 0)).
    base_matrix : numpy.ndarray, optional
        Base matrix for intersection. If None, uses matrix of ones.
    orientation : str, optional
        River orientation "h" or "v" (default: "h").
    plot : bool, optional
        If True, generates a plot of the result.

    Returns
    -------
    numpy.ndarray
        Matrix retaining base matrix values only where river exists; 0 elsewhere.
    """
    # Generate Kinoshita curve trajectory
    if orientation == "h":
        x_curve, y_curve, _ = kinoshita_curve(theta_0, lambda_value, j_s, j_f, n, m_points, ds)
        start_row, start_col = cell_start
    else:
        y_curve, x_curve, _ = kinoshita_curve(theta_0, lambda_value, j_s, j_f, n, m_points, ds)
        start_row, start_col = cell_start

    x_curve = x_curve + start_col
    y_curve = y_curve + start_row

    # Create river mask in an empty grid
    # Assumes horizontal direction (columns) corresponds to x and vertical (rows) to y
    river_mask = np.zeros((nrow, ncol), dtype=int)

    # Traverse each curve point and mark in the mask
    for xi, yi in zip(x_curve, y_curve):
        # Convert coordinates to integer indices for the matrix
        col = int(round(xi))
        row = int(round(yi))

        # Calculate the range of cells spanning the river width
        row_min = max(0, row - width_riv // 2)
        row_max = min(nrow, row + width_riv // 2 + 1)
        col_min = max(0, col - width_riv // 2)
        col_max = min(ncol, col + width_riv // 2 + 1)

        # Mark cells corresponding to the river
        river_mask[row_min:row_max, col_min:col_max] = 1

    # If no base matrix provided, use a matrix of ones
    if base_matrix is None:
        base_matrix = np.ones((nrow, ncol), dtype=int)

    # Perform intersection: retain base matrix values where river is present
    intersection = np.where(river_mask == 1, base_matrix, 0)

    if plot:
        plt.figure(figsize=(10, 6))
        plt.imshow(intersection, cmap="Blues", origin="upper")
        plt.title("River (Kinoshita)")
        plt.xlabel("Col")
        plt.ylabel("Row")
        plt.colorbar(label="River (1 = Yes, 0 = No)")
        i_row, i_col = np.where(intersection != 0)
        for x, y in zip(i_col, i_row):
            plt.text(x, y, f'({x},{y})', fontsize=8, color='red', ha='center', va='center')

    return intersection


# =============================================================================
# Coordinate Utilities
# =============================================================================

def adjust_coordinates(i_riv, j_riv, orientation="v"):
    """
    Adjust river coordinates for proper flow direction and connectivity.

    Parameters
    ----------
    i_riv : numpy.ndarray
        River row coordinates.
    j_riv : numpy.ndarray
        River column coordinates.
    orientation : str
        River orientation ('v' or 'h').

    Returns
    -------
    tuple
        (adjusted_i_riv, adjusted_j_riv)

    Warnings
    --------
    This function is deprecated. Consider using RiverGeometry class instead.
    """
    if orientation == 'h':
        idx_sort = np.lexsort((i_riv, j_riv))
        i_riv, j_riv = j_riv[idx_sort], i_riv[idx_sort]

    idx_place, idx_unique_counts = np.unique(i_riv, return_counts=True, return_index=True)[1:]
    idx_start_repeats = np.ones(i_riv.shape)
    idx_start_repeats[idx_place - 1] = 0

    for n in range(len(idx_start_repeats) - 1):
        if idx_start_repeats[n] == 1 and idx_start_repeats[n + 1] == 1:
            idx_start_repeats[n + 1] = 0

    idx_start_repeats = np.where(idx_start_repeats)[0]
    idx_counts_riv = idx_unique_counts[idx_unique_counts != 1]

    for n, idx_P in enumerate(idx_start_repeats):
        before_duplicate = j_riv[idx_P - 1]
        start_duplicate = j_riv[idx_P]
        cond = before_duplicate > start_duplicate

        if cond:
            j_riv[idx_P:idx_P + idx_counts_riv[n]] = j_riv[idx_P:idx_P + idx_counts_riv[n]][::-1]

    if orientation == 'h':
        i_riv, j_riv = j_riv, i_riv

    return i_riv, j_riv


# =============================================================================
# MODFLOW Stress Period Data
# =============================================================================

def generate_river_stress_period_data(nrow, ncol, topo, del_rbot, nper,
                                      stage_per_time, river_params, type_gen='basic'):
    """
    Generate stress period data for a river in a MODFLOW model.

    Parameters
    ----------
    nrow : int
        Number of model rows.
    ncol : int
        Number of model columns.
    topo : numpy.ndarray
        Topography matrix.
    del_rbot : float
        Depth below topography for river bottom.
    nper : int
        Number of stress periods.
    stage_per_time : list
        List of stage values per period.
    river_params : dict
        River parameters for generate_river function.
    type_gen : str
        Type of river generation: 'basic' or 'kinoshita'.

    Returns
    -------
    dict
        Dictionary with river stress period data.

    Warnings
    --------
    This function is deprecated. Consider using RiverGeometry class instead.
    """
    warnings.warn("generate_river_stress_period_data is deprecated.",
                  DeprecationWarning, stacklevel=2)

    if type_gen == "basic":
        river = generate_river(nrow, ncol, **river_params)
    elif type_gen == "kinoshita":
        river = generate_river_kinoshita(nrow, ncol, **river_params)

    # Extract river topography
    extract_topo_riv_raster = np.where(river == 1, topo, np.nan)
    riv_topo = extract_topo_riv_raster[~np.isnan(extract_topo_riv_raster)]
    nc_riv = riv_topo.size
    steady = np.repeat(np.arange(nper), int(nc_riv / nper))

    # River coordinates
    i_riv, j_riv = np.where(river == 1)

    # Calculate stage, rbot and cond
    rbot = riv_topo - del_rbot
    lay = np.repeat(0, i_riv.size)

    i_riv, j_riv = adjust_coordinates(i_riv, j_riv, river_params["orientation"])

    river_SFR = {
        "i": i_riv,
        "j": j_riv,
        "top": extract_topo_riv_raster,
        "rbot": rbot,
        "lay": lay,
        "steady": steady,
        "river": river,
        "nc_river": nc_riv
    }
    return river_SFR


# =============================================================================
# Hydrograph Generation
# =============================================================================

def hydrograph_cooper_rorabaugh(t, H_0, eta_0, eta_d, t_d):
    """
    Calculate the deterministic Cooper and Rorabaugh (1963) hydrograph.

    Parameters
    ----------
    t : array_like
        Time(s) to evaluate the hydrograph.
    H_0 : float
        River base level.
    eta_0 : float
        Event intensity (H_p = eta_0 * H_0).
    eta_d : float
        Hydrograph asymmetry (t_p = eta_d * t_d).
    t_d : float
        Total duration of the flood event.

    Returns
    -------
    numpy.ndarray
        River level at each time t.

    Warnings
    --------
    This function is deprecated. Consider using HydrogramGenerator class instead.
    """
    warnings.warn("hydrograph_cooper_rorabaugh is deprecated.",
                  DeprecationWarning, stacklevel=2)

    t = np.asarray(t)
    H_p = eta_0 * H_0
    t_p = eta_d * t_d
    omega = 2 * np.pi / t_d
    angle = omega * t_p / 2

    # Enhanced numerical stability
    if np.abs(np.tan(angle)) > 1e-10:
        beta = omega / np.tan(angle)
    else:
        beta = omega  # Fallback for small angles

    # Ensure beta is reasonable
    beta = max(beta, omega * 0.1)

    # Calculate H_b(t)
    H_b = np.where(
        (t >= 0) & (t <= t_d),
        H_0 + H_p * np.exp(-beta * (t - t_p)) * (1 - np.cos(omega * t)) / (1 - np.cos(omega * t_p)),
        H_0  # Outside the range [0, t_d]
    )

    return H_b


# =============================================================================
# Results Extraction and Analysis
# =============================================================================

def extract_sfr_data(model, reach_numbers, labels, orientation='v'):
    """
    Extract SFR data for specific reaches or all reaches if reach_numbers=[-1].

    Parameters
    ----------
    model : flopy.modflow.Modflow
        MODFLOW model object.
    reach_numbers : int or list
        Reach number(s) to extract (0-indexed). Use [-1] for all reaches.
    labels : list
        SFR variables to extract ['Qin', 'Qaquifer', 'Qout', 'stage', 'depth', 'gradient'].
    orientation : str
        River orientation ('v' or 'h').

    Returns
    -------
    dict
        Dictionary with data structured by reach and variable.
    """
    # Read SFR file
    sfr_df = SfrFile(f"{model.model_ws}\\{model.name}.sfr.out").df

    # Handle special case: extract all reaches if reach_numbers=[-1]
    if reach_numbers == [-1]:
        unique_reaches = sorted(sfr_df['reach'].unique())
        reach_numbers = [rn - 1 for rn in unique_reaches]  # Convert to 0-indexed

    # Normalize reach_numbers to list
    if isinstance(reach_numbers, int):
        reach_numbers = [reach_numbers]

    # Convert to file indices (1-indexed)
    file_reaches = [rn + 1 for rn in reach_numbers]

    # Filter data
    filtered_df = sfr_df[sfr_df['reach'].isin(file_reaches)]

    # Output data structure
    data_dict = {}

    for file_rn, orig_rn in zip(file_reaches, reach_numbers):
        reach_df = filtered_df[filtered_df['reach'] == file_rn]
        data_dict[orig_rn] = {}

        for label in labels:
            if label in reach_df.columns:
                values = reach_df[label].to_numpy()
                data_dict[orig_rn][label] = values
            else:
                print(f"Warning: Variable '{label}' not found in SFR data")

    return data_dict


def analyze_peak_integrals(t, y, eps=0.0):
    """
    Analyze peak integrals by automatically finding integration bounds.

    Searches where the signal stops being >eps (for positive part) or
    < -eps (for negative part).

    Parameters
    ----------
    t : numpy.ndarray
        Time array.
    y : numpy.ndarray
        Signal array.
    eps : float, optional
        Threshold for signal detection (default: 0.0).

    Returns
    -------
    dict
        Dictionary containing:
        - 't_max', 'y_max': Time and value of positive peak
        - 't_min', 'y_min': Time and value of negative peak
        - 'bounds_pos': (i0, i1) indices for positive peak bounds
        - 'bounds_neg': (i0, i1) indices for negative peak bounds
        - 'int_pos': Integral of positive region
        - 'int_neg': Integral of negative region
    """
    # 1) Locate peaks
    i_max = np.argmax(y)
    i_min = np.argmin(y)
    t_max, y_max = t[i_max], y[i_max]
    t_min, y_min = t[i_min], y[i_min]

    n = len(y)

    # 2) Search bounds by scanning
    def scan_bounds(idx_peak, sign):
        """
        sign = +1 for positive peak (searches y[idx] <= eps),
        sign = -1 for negative peak (searches y[idx] >= -eps).
        Returns (i0, i1) left and right boundary indices.
        """
        # Left
        i0 = idx_peak
        while i0 > 0 and sign * y[i0] > eps:
            i0 -= 1
        # Right
        i1 = idx_peak
        while i1 < n - 1 and sign * y[i1] > eps:
            i1 += 1
        return i0, i1

    # For positive peak use sign=+1, for negative use sign=-1
    i0_pos, i1_pos = scan_bounds(i_max, +1)
    i0_neg, i1_neg = scan_bounds(i_min, -1)

    # 3) Integrals
    int_pos = np.trapz(y[i0_pos:i1_pos + 1], t[i0_pos:i1_pos + 1])
    int_neg = np.trapz(y[i0_neg:i1_neg + 1], t[i0_neg:i1_neg + 1])

    return {
        't_max': t_max, 'y_max': y_max,
        't_min': t_min, 'y_min': y_min,
        'bounds_pos': (i0_pos, i1_pos),
        'bounds_neg': (i0_neg, i1_neg),
        'int_pos': int_pos,
        'int_neg': int_neg
    }


def compute_exchange_metrics(t, y_raw, t_delay, window=11):
    """
    Compute integration bounds, areas, and peak indices/values for exchange flow.

    Parameters
    ----------
    t : array-like
        Time vector.
    y_raw : array-like
        Raw signal vector (same length as t).
    t_delay : float
        Time at which to start integration.
    window : int, optional
        Odd-sized moving-average window for smoothing (default: 11).

    Returns
    -------
    dict
        Dictionary containing smoothed signal, indices, areas, and peak information.
    """
    # 1) Smooth
    y_smooth = pd.Series(y_raw) \
                   .rolling(window, center=True, min_periods=1) \
                   .mean() \
                   .to_numpy()

    # 2) Start index
    start_idx = np.searchsorted(t, t_delay)

    # 3) Negative crossing - find first point where signal becomes negative
    neg_idx = None
    for i in range(start_idx, len(y_smooth)):
        if y_smooth[i] < 0:
            neg_idx = i
            break

    if neg_idx is None:
        print("Warning: Signal never goes negative after t_delay. Using end of signal.")
        neg_idx = len(y_smooth) - 1

    # 4) Return-to-zero crossing
    zero_idx = None
    for i in range(neg_idx + 1, len(y_smooth)):
        if y_smooth[i] >= 0:
            zero_idx = i
            break

    if zero_idx is None:
        print("Warning: Signal never returns to positive. Using end of signal.")
        zero_idx = len(y_smooth) - 1

    # 5) Areas (trapezoid rule)
    if neg_idx > start_idx:
        area_pos = np.trapezoid(y_smooth[start_idx:neg_idx + 1], t[start_idx:neg_idx + 1])
    else:
        area_pos = 0.0

    if zero_idx > neg_idx:
        area_neg = np.trapezoid(y_smooth[neg_idx:zero_idx + 1], t[neg_idx:zero_idx + 1])
    else:
        area_neg = 0.0

    # 6) Peaks with safety checks
    if neg_idx > start_idx:
        sub_pos = y_smooth[start_idx:neg_idx + 1]
        rel_peak_pos = np.argmax(sub_pos)
        peak_pos_idx = start_idx + rel_peak_pos
        peak_pos_val = y_smooth[peak_pos_idx]
    else:
        peak_pos_idx = start_idx
        peak_pos_val = y_smooth[start_idx]

    if zero_idx > neg_idx:
        sub_neg = y_smooth[neg_idx:zero_idx + 1]
        rel_peak_neg = np.argmin(sub_neg)
        peak_neg_idx = neg_idx + rel_peak_neg
        peak_neg_val = y_smooth[peak_neg_idx]
    else:
        peak_neg_idx = neg_idx
        peak_neg_val = y_smooth[neg_idx]

    return {
        'y_smooth': y_smooth,
        'start_idx': start_idx,
        'neg_idx': neg_idx,
        'zero_idx': zero_idx,
        'Qin_exc': area_pos,
        'Qout_exc': area_neg,
        'peak_pos_idx': peak_pos_idx,
        'peak_pos_val': peak_pos_val,
        'peak_pos_time': t[peak_pos_idx],
        'peak_neg_idx': peak_neg_idx,
        'peak_neg_val': peak_neg_val,
        'peak_neg_time': t[peak_neg_idx],
    }


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_results_hds_cbc(path_results, model, int_p, layers=0, **kwargs):
    """
    Plot hydraulic heads and storage from MODFLOW output files.

    Parameters
    ----------
    path_results : str
        Path to results directory.
    model : flopy.modflow.Modflow
        MODFLOW model object.
    int_p : int
        Time step index to plot.
    layers : int, optional
        Layer to plot. Use -1 for sum of all layers (default: 0).
    **kwargs
        Additional arguments passed to imshow.
    """
    n_hds = model.name + ".hds"
    n_cbc = model.name + ".cbc"
    hds = flopy.utils.HeadFile(path_results + n_hds)
    heads = hds.get_data(kstpkper=(int_p, 1))
    cbb = flopy.utils.binaryfile.CellBudgetFile(path_results + n_cbc)
    frf = cbb.get_data(text="FLOW RIGHT FACE")[int_p + 1]
    fff = cbb.get_data(text="FLOW FRONT FACE")[int_p + 1]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge((frf, fff, None), model)

    fig, ax = plt.subplots(1, 2, figsize=(10, 7.5))

    ax[0].set_title("Volumetric discharge (" + r"$L^3/T$" + ")")
    mapview = flopy.plot.PlotMapView(model=model)
    mapview.plot_ibound(ax=ax[0])
    quadmesh = mapview.plot_array(heads, alpha=1, cmap='turbo', ax=ax[0])

    plt.colorbar(quadmesh, ax=ax[0], orientation="horizontal", pad=.05, label='Hydraulic head [m]')

    if layers == -1:
        storage_data = cbb.get_data(text="STORAGE")[int_p].sum(axis=0)
    else:
        storage_data = cbb.get_data(text="STORAGE")[int_p][layers]

    ax[1].set_title("Storage")
    s = ax[1].imshow(storage_data, cmap='turbo', **kwargs)
    plt.colorbar(s, orientation="horizontal", ax=ax[1], pad=.05, label="")
    plt.tight_layout()


def plot_exchange_flow(t, y_raw, metrics, t_delay):
    """
    Plot raw vs smoothed exchange flow, shaded integrals, and mark peaks.

    Parameters
    ----------
    t : array-like
        Time vector.
    y_raw : array-like
        Raw signal.
    metrics : dict
        Output from compute_exchange_metrics.
    t_delay : float
        Vertical reference line at t_delay.
    """
    y_smooth = metrics['y_smooth']
    si, ni, zi = metrics['start_idx'], metrics['neg_idx'], metrics['zero_idx']
    ap, an = metrics['Qin_exc'], metrics['Qout_exc']
    ppi, ppi_val = metrics['peak_pos_idx'], metrics['peak_pos_val']
    pni, pni_val = metrics['peak_neg_idx'], metrics['peak_neg_val']
    ppt, pnt = metrics['peak_pos_time'], metrics['peak_neg_time']

    # Style
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'font.size': 12,
        'axes.linewidth': 1.2,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 11,
        'legend.frameon': False,
        'lines.linewidth': 1.5,
        'figure.dpi': 300,
    })

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot signals
    ax.plot(t, y_raw, color='gray', alpha=0.4, label='Raw signal')
    ax.plot(t, y_smooth, color='black', label='Smoothed signal')

    # Shaded integrals
    ax.fill_between(t[si:ni], y_smooth[si:ni],
                    where=y_smooth[si:ni] >= 0,
                    color='orange', alpha=0.3,
                    label=f"Positive area\n({ap:.2e})")
    ax.fill_between(t[ni:zi], y_smooth[ni:zi],
                    where=y_smooth[ni:zi] <= 0,
                    color='purple', alpha=0.3,
                    label=f"Negative area\n({an:.2e})")

    # Reference lines
    ax.axvline(t_delay, color='green', linestyle='--', label='t_delay')

    # Peak markers
    ax.scatter([ppt], [ppi_val],
               color='orange', edgecolor='black', zorder=5,
               label=f"Max +: {ppi_val:.2e} @ {ppt:.0f}")
    ax.scatter([pnt], [pni_val],
               color='purple', edgecolor='black', zorder=5,
               label=f"Max -: {pni_val:.2e} @ {pnt:.0f}")

    # Labels & legend
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Exchange flow (Q)')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.show()


def plot_sfr_results(model, river, reach_number, labels, orientation='h',
                     color=["b", "r"], marker=["+", '.'], T_range=None,
                     scalex='linear', scaley='linear', integral=False,
                     range_integral=[5, 5], middle_in_out=20,
                     dict_values=False, figsize=(15, 7.5), xlim=None, **kwargs):
    """
    Filter DataFrame by specified reach numbers and plot selected columns.

    Parameters
    ----------
    model : flopy.modflow.Modflow
        Object containing model_ws and name.
    river : numpy.ndarray
        Array with river information.
    reach_number : int or list
        Reach number(s) to filter (0-indexed).
    labels : list
        Labels to plot, e.g. ['Qin', 'Qaquifer', 'Qout', 'stage', 'depth', 'gradient'].
    orientation : str
        "v" or "h" for river orientation.
    color : list
        Colors for plotting.
    marker : list
        Markers for plotting.
    T_range : float
        Total time range for conversion.
    scalex, scaley : str
        Axis scales ('linear' or 'log').
    integral : bool
        If True, calculate and display integrals.
    range_integral : list
        Integration range indices.
    middle_in_out : int
        Middle point separating inflow from outflow.
    dict_values : bool
        If True, return values dictionary.
    figsize : tuple
        Figure size.
    xlim : tuple
        X-axis limits.
    **kwargs
        Additional arguments passed to imshow.

    Returns
    -------
    dict or None
        Values dictionary if dict_values=True.
    """
    # Read DataFrame from SFR file
    sfr_df = SfrFile(model.model_ws + f"\\{model.name}.sfr.out").df

    # Convert reach_number to list if received as integer
    if isinstance(reach_number, int):
        reach_numbers = [reach_number]
    else:
        reach_numbers = reach_number

    # Filter DataFrame for specified reaches
    # Add 1 since file reaches are indexed from 1
    file_reaches = [rn + 1 for rn in reach_numbers]
    filtered_df = sfr_df[sfr_df['reach'].isin(file_reaches)]

    # Setup figure with two subplots: one for river image and one for data
    fig, ax = plt.subplots(1, 2, figsize=figsize)
    fig.suptitle('Reach Data Visualization')

    # Subplot 0: Show river image and mark points corresponding to each reach
    i_riv, j_riv = np.where(river != 0)
    i_riv, j_riv = adjust_coordinates(i_riv, j_riv, orientation=orientation)
    for rn, cp in zip(reach_numbers, color):
        ax[0].scatter(j_riv[rn], i_riv[rn], label=f'Reach {rn + 1}', marker='s', s=60, c=cp)
    ax[0].imshow(river, **kwargs)
    ax[0].legend()

    values_dict = {}

    # Subplot 1: Plot series and mark maximum of each variable
    for file_rn, orig_rn, c in zip(file_reaches, reach_numbers, color):
        reach_df = filtered_df[filtered_df['reach'] == file_rn]
        n_steps = len(reach_df)
        conversion_factor = T_range / n_steps

        t = np.arange(len(reach_df))
        # Create a sub-dictionary for this reach
        values_dict[orig_rn] = {}
        for col, m in zip(labels, marker):
            # Store ALL values for this variable in current reach
            values = reach_df[col].to_numpy()
            values_dict[orig_rn][col] = values
            idx_max = np.argmax(values)
            max_val = values.max()

            ax[1].plot(t, values, label=f'{col} (Reach {orig_rn + 1})', marker=m, c=c)
            ax[1].scatter(t[idx_max], max_val, c='k', marker='s', zorder=10)

    ax[1].set_xscale(scalex)
    ax[1].set_yscale(scaley)
    ax[1].set_xlabel('Step')

    if xlim is not None:
        ax[1].set_xlim(xlim)

    secax = ax[1].secondary_xaxis('top',
                                   functions=(lambda x: x * conversion_factor,
                                             lambda x: x / conversion_factor))
    secax.set_xlabel('Time')

    if "q" in labels[-1].lower():
        ax[1].set_ylabel('Q (m^3)/s')
    elif "s" in labels[-1].lower():
        ax[1].set_ylabel("Stage (m)")

    ax[1].legend()

    if integral:
        for rn in reach_numbers:
            Q_Aquifer = values_dict[rn]["Qaquifer"]
            t_aquifer = conversion_factor * t
            idx_max_in = np.argmax(Q_Aquifer)
            idx_min_out = np.argmin(Q_Aquifer)

            t_in = t_aquifer[idx_max_in]
            t_out = t_aquifer[idx_min_out]

            Q_Aquifer_in = Q_Aquifer[range_integral[0]:middle_in_out]
            Q_Aquifer_out = Q_Aquifer[middle_in_out:range_integral[1]]

            Q_exc_in = np.trapezoid(Q_Aquifer_in, t_aquifer[range_integral[0]:middle_in_out])
            Q_exc_out = np.trapezoid(Q_Aquifer_out, t_aquifer[middle_in_out:range_integral[1]])

            print(f"Diff {t_out - t_in}", end='\n')
            print(f"Qexchange in for reach {rn + 1}: {Q_exc_in:.4f}", end='\n')
            print(f"Qexchange out for reach {rn + 1}: {Q_exc_out:.4f}", end='\n')

        ax[1].axhline(0, color='purple', linestyle='-')
        ax[1].axvline(t[range_integral[0]], color='r', linestyle='--')
        ax[1].axvline(t[middle_in_out], color='k', linestyle='--')
        ax[1].axvline(t[range_integral[1]], color='r', linestyle='--')

    if len(reach_numbers) == 2:
        diff_val = []
        r1 = reach_numbers[0]
        r2 = reach_numbers[1]

        for col in labels:
            max1 = np.max(values_dict[r1][col])
            idx_max1 = np.argmax(values_dict[r1][col])

            max2 = np.max(values_dict[r2][col])
            idx_max2 = np.argmax(values_dict[r2][col])

            diff = max1 - max2
            diff_T = (t[idx_max1] - t[idx_max2]) * conversion_factor

            diff_val.append(f"{col}: {diff:.4f} and T_diff {diff_T * conversion_factor:.4f} days")

        print(diff_val)
        ax[1].set_title(f"reach {r1 + 1} - reach {r2 + 1}\n" + "\n".join(diff_val) + "\n")

    if dict_values:
        return values_dict


def plot_sfr_data(data_dict, river, orientation='v', operation='sum',
                  colors=["C1", "C2", "C3", "C4", "C5"],
                  markers=[".", "+", "o", "s", "d"],
                  total_simulation_time=None, time_unit='auto', flow_units='m3/s',
                  river_width=None, figsize=(15, 7.5),
                  xlim=None, scalex='linear', scaley='linear', aspect=1):
    """
    Plot extracted SFR data with support for aggregation when all reaches are selected.

    Parameters
    ----------
    data_dict : dict
        Data extracted by extract_sfr_data.
    river : numpy.ndarray
        River array for spatial visualization.
    orientation : str
        River orientation ('v' or 'h').
    operation : str
        Aggregation operation ('sum' or 'mean') when plotting all reaches.
    colors, markers : list
        Plotting styles.
    total_simulation_time : float, optional
        Total simulation time in seconds.
    time_unit : str, optional
        Time unit for x-axis: 'auto', 'seconds', 'hours', 'days', 'steps'.
    flow_units : str, optional
        Flow units for y-axis: 'm3/s', 'm3/day', 'm2/s', 'm2/day'.
    river_width : float, optional
        River width in meters (required for m2 units).
    figsize : tuple
        Figure size.
    xlim : tuple
        X-axis limits.
    scalex, scaley : str
        Axis scales.
    aspect : float
        Aspect ratio for river visualization.

    Returns
    -------
    tuple
        (fig, ax) matplotlib figure and axes objects.
    """
    reach_numbers = list(data_dict.keys())
    labels = list(data_dict[reach_numbers[0]].keys())

    # Determine number of time steps
    num_time_steps = len(data_dict[reach_numbers[0]][labels[0]])

    # Time conversion
    if total_simulation_time is not None:
        time_factor = total_simulation_time / num_time_steps
        if time_unit == 'auto':
            if total_simulation_time < 3600:
                time_symbol = 'seconds'
            elif total_simulation_time < 86400:
                time_factor /= 3600
                time_symbol = 'hours'
            else:
                time_factor /= 86400
                time_symbol = 'days'
        elif time_unit == 'hours':
            time_factor /= 3600
            time_symbol = 'hours'
        elif time_unit == 'days':
            time_factor /= 86400
            time_symbol = 'days'
        else:
            time_symbol = 'seconds'
    else:
        time_factor = 1.0
        time_symbol = 'steps'

    # Create figure
    fig, ax = plt.subplots(1, 2, figsize=figsize)

    if len(reach_numbers) > 5:
        # Aggregation mode
        fig.suptitle(f'SFR Data Visualization - {operation.capitalize()} of All Reaches')

        i_riv, j_riv = np.where(river != 0)
        i_riv, j_riv = adjust_coordinates(i_riv, j_riv, orientation=orientation)

        ax[0].imshow(river, aspect=aspect, cmap='Blues')
        ax[0].scatter(j_riv, i_riv, c='red', s=20, alpha=0.7,
                     label=f'All Reaches ({len(reach_numbers)})')
        ax[0].legend()
        ax[0].set_title('River Layout - All Reaches')

        for label in labels:
            all_data = np.array([data_dict[rn][label] for rn in reach_numbers])

            if operation == 'sum':
                aggregated_data = np.sum(all_data, axis=0)
            else:
                aggregated_data = np.mean(all_data, axis=0)

            t_steps = np.arange(len(aggregated_data)) + 1
            t_real = (t_steps - 1) * time_factor

            ax[1].plot(t_real, aggregated_data,
                      label=f'{label} ({operation})',
                      marker=markers[0], linewidth=2)

        ax[1].set_xlabel(f'Time ({time_symbol})')
        ax[1].set_ylabel('Q (m³/s)' if 'q' in labels[0].lower() else 'Value')
    else:
        # Individual reaches mode
        fig.suptitle('SFR Data Visualization - Individual Reaches')

        i_riv, j_riv = np.where(river != 0)
        i_riv, j_riv = adjust_coordinates(i_riv, j_riv, orientation=orientation)

        for idx, (rn, color) in enumerate(zip(reach_numbers, colors)):
            ax[0].scatter(j_riv[rn], i_riv[rn],
                         label=f'Reach {rn + 1}',
                         marker='s', s=60, c=color)

        ax[0].imshow(river, aspect=aspect)
        ax[0].legend()
        ax[0].set_title('River Layout with Selected Reaches')

        for idx, (rn, color) in enumerate(zip(reach_numbers, colors)):
            for nmarker, label in enumerate(labels):
                values = data_dict[rn][label]
                t_steps = np.arange(len(values)) + 1
                t_real = (t_steps - 1) * time_factor

                ax[1].plot(t_real, values,
                           label=f'{label} (Reach {rn + 1})',
                           marker=markers[nmarker % len(markers)], c=color)

        ax[1].set_xlabel(f'Time ({time_symbol})')
        ax[1].set_ylabel('Q (m³/s)' if 'q' in labels[0].lower() else 'Value')

    ax[1].set_xscale(scalex)
    ax[1].set_yscale(scaley)
    if xlim is not None:
        ax[1].set_xlim(xlim)
    ax[1].legend()
    ax[1].grid(True)

    plt.tight_layout()
    return fig, ax
