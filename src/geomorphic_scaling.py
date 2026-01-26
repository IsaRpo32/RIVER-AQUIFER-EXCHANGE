"""
Geomorphic Scaling Functions for River-Aquifer Systems

This module provides empirical scaling relationships to estimate river and aquifer
geometric properties based on drainage area. These relationships are commonly used
in regional hydrological modeling.

References:
    - Bieger, K., et al. (2015). Development and evaluation of bankfull hydraulic
      geometry relationships for the physiographic regions of the United States.
    - Zhou, Y., et al. (2024). Scaling relationships for river systems.
    - Lane, S. N., et al. (2021). Geomorphic scaling in river networks.
"""

import numpy as np


def get_w_bkf(da_km2):
    """
    Calculate the bankfull width (w_bkf) based on drainage area.

    Uses the empirical relationship from Bieger et al. (2015).

    Parameters
    ----------
    da_km2 : float
        Drainage area in square kilometers.

    Returns
    -------
    float
        Bankfull width in meters.
    """
    w_bkf_m = 2.9 * da_km2**0.352  # Bieger et al. 2015
    return w_bkf_m


def get_d_bkf(da_km2):
    """
    Calculate the bankfull depth (d_bkf) based on drainage area.

    Uses the empirical relationship from Bieger et al. (2015).

    Parameters
    ----------
    da_km2 : float
        Drainage area in square kilometers.

    Returns
    -------
    float
        Bankfull depth in meters.
    """
    d_bkf_m = 0.3 * da_km2**0.213  # Bieger et al. 2015
    # Alternative: d_bkf_m = 0.2156 * da_km2**0.3018  # Eash (1993)
    return d_bkf_m


def get_S(da_km2):
    """
    Calculate the channel slope (S) based on drainage area.

    Uses the empirical relationship from Zhou et al. (2024) and Lane et al. (2021).

    Parameters
    ----------
    da_km2 : float
        Drainage area in square kilometers.

    Returns
    -------
    float
        Slope (dimensionless).
    """
    # S = 0.14 * da_km2**-0.45  # Mean basin slope - results were not optimal
    S = 0.001 * da_km2**-0.35  # Zhou et al. 2024, Lane et al. 2021
    return S


def get_L(w_bkf_m):
    """
    Calculate the reach length (L) based on bankfull width.

    Parameters
    ----------
    w_bkf_m : float
        Bankfull width in meters.

    Returns
    -------
    float
        Reach length in meters.

    Notes
    -----
    NOTE: TWO VALUES recommended based on basin size:
    - For basins of 100 and 1000 km², use factor = 150
    - For basins of 100000 km² and larger, use factor = 100
    """
    L = 150 * w_bkf_m
    return L


def get_dx(L):
    """
    Calculate the cell size (dx) for numerical discretization.

    Parameters
    ----------
    L : float
        Reach length in meters.

    Returns
    -------
    float
        Cell size in meters.

    Notes
    -----
    dx is needed to estimate Criv in a segment and constrain random sampling.
    Current value is empirical (1% of L).
    """
    dx = 0.01 * L
    return dx


def get_dr(d_bkf):
    """
    Calculate the streambed thickness (dr) based on bankfull depth.

    Parameters
    ----------
    d_bkf : float
        Bankfull depth in meters.

    Returns
    -------
    float
        Streambed thickness in meters.

    Notes
    -----
    Values smaller than 0.3 were tested but resulted in very unstable simulations.
    """
    dr = 0.3 * d_bkf
    return dr


def compute_conductance(K, dr, w_bkf_m, dx):
    """
    Compute river-aquifer conductance.

    Parameters
    ----------
    K : float
        River hydraulic conductivity (m/s).
    dr : float
        Thickness of the streambed (m).
    w_bkf_m : float
        Width of the river (m).
    dx : float
        Length of the river segment (m), typically a fraction of L.

    Returns
    -------
    float
        Conductance C = (K * width * length) / dr (m²/s).
    """
    Criv = (K * w_bkf_m * dx) / dr
    return Criv


def get_A(w_bkf_m, d_bkf_m):
    """
    Calculate the river cross-sectional area.

    Parameters
    ----------
    w_bkf_m : float
        Bankfull width in meters.
    d_bkf_m : float
        Bankfull depth in meters.

    Returns
    -------
    float
        River section area in square meters.
    """
    A_m2 = w_bkf_m * d_bkf_m
    return A_m2


def get_D(d_bkf_m, nd=30):
    """
    Calculate the aquifer depth (D) based on bankfull depth.

    Parameters
    ----------
    d_bkf_m : float
        Bankfull depth in meters.
    nd : int, optional
        Multiplier factor (default: 30).

    Returns
    -------
    float
        Aquifer depth in meters.

    Notes
    -----
    Recommended to test values between 10 and 100.
    Factor of 10 is discarded as it generates unstable simulations.
    """
    D_m = nd * d_bkf_m
    return D_m


def get_B(w_bkf_m, nb=120):
    """
    Calculate the aquifer lateral extension (B) based on bankfull width.

    Parameters
    ----------
    w_bkf_m : float
        Bankfull width in meters.
    nb : int, optional
        Multiplier factor (default: 120).

    Returns
    -------
    float
        Aquifer lateral extension in meters.

    Notes
    -----
    Recommended to test values between 100 and 1000.
    Factor of 100 is discarded as it generates unstable simulations.
    """
    B_m = nb * w_bkf_m
    return B_m


def get_realistic_K_range(da_km2):
    """
    Return a realistic range of hydraulic conductivity K.

    Parameters
    ----------
    da_km2 : float
        Watershed area in km² (currently unused but kept for future extensions).

    Returns
    -------
    tuple
        (K_min, K_max) hydraulic conductivity range in m/s.
    """
    return 1e-6, 1e-1


def get_K(area_list, C_threshold=0.2, num_K=None, dx_factor=0.025, L=None):
    """
    Sample combinations of hydraulic conductivity (K) for given watershed areas.

    For each watershed area, samples K values such that the computed conductance
    does not exceed C_threshold.

    Parameters
    ----------
    area_list : list of float
        List of watershed areas (km²).
    C_threshold : float
        Maximum allowable conductance (m²/s).
    num_K : int, optional
        Number of K values to sample logarithmically between K_min and K_max.
    dx_factor : float, optional
        Fraction of segment length L to use as dx (i.e., dx = dx_factor * L).
    L : float, optional
        Reach length. If None, computed from width.

    Returns
    -------
    dict
        A dictionary mapping each area A to results:
        {"C": [...], "K": [...], "dr": [...], "L": [...], "dx": [...], "w": [...]}
    """
    results = {}

    for A in area_list:
        C_all, K_all, dr_all, L_all, dx_all, w_all = [], [], [], [], [], []

        # Estimate streambed thickness dr (m) and river width w (m)
        d_bkf_aux = get_d_bkf(A)
        dr = get_dr(d_bkf_aux)
        w = get_w_bkf(A)

        # Compute segment length L (m) and cell size dx (m)
        L_val = L
        dx = dx_factor * L_val

        # Determine realistic K range for this watershed area
        K_min, K_max = get_realistic_K_range(A)

        # Generate logarithmically spaced K values
        K_values = np.logspace(np.log10(K_min), np.log10(K_max), num=num_K)

        np.random.shuffle(K_values)
        for K in K_values:
            C = compute_conductance(K, dr, w, dx)
            if C <= C_threshold:
                C_all.append(C)
                K_all.append(K)
                dr_all.append(dr)
                L_all.append(L_val)
                dx_all.append(dx)
                w_all.append(w)

        results[A] = {
            "C": C_all,
            "K": K_all,
            "dr": dr_all,
            "L": L_all,
            "dx": dx_all,
            "w": w_all,
        }

    return results
