"""
Dimensionless Numbers for River-Aquifer Exchange Analysis

This module provides functions to calculate dimensionless parameters (Pi numbers)
used in the analysis of river-aquifer interaction during flood events.

References:
    - Cooper, H. H., & Rorabaugh, M. I. (1963). Ground-water movements and bank
      storage due to flood stages in surface streams.
"""

import numpy as np


def get_tc(Sy, Ka, B, d_bkf):
    """
    Calculate the characteristic hydraulic time (tc).

    The hydraulic time represents the timescale of pressure propagation
    in the aquifer system.

    Parameters
    ----------
    Sy : float
        Specific yield of the aquifer (dimensionless).
    Ka : float
        Aquifer hydraulic conductivity (m/s).
    B : float
        Aquifer lateral extension from the river (m).
    d_bkf : float
        Bankfull depth of the river (m).

    Returns
    -------
    float
        Hydraulic time in seconds.
    """
    tc = (Sy * B**2) / (Ka * d_bkf)
    return tc


def get_gamma_d(tc, t_d):
    """
    Calculate the sensitivity to stage perturbation (Gamma_d).

    This parameter represents the ratio between the hydraulic propagation
    time and the flood event duration.

    Parameters
    ----------
    tc : float
        Hydraulic propagation time in seconds.
    t_d : float
        Event duration in seconds.

    Returns
    -------
    float
        Gamma_d (dimensionless).
    """
    gamma_d = tc / t_d
    return gamma_d


def get_lambda_v(tc, Kr, dr):
    """
    Calculate the streambed conductance parameter (lambda_v or pi1).

    This parameter characterizes the vertical exchange capacity through
    the streambed relative to the aquifer response time.

    Parameters
    ----------
    tc : float
        Hydraulic time in seconds.
    Kr : float
        Riverbed hydraulic conductivity (m/s).
    dr : float
        Streambed thickness (m).

    Returns
    -------
    float
        Lambda_v (dimensionless).
    """
    lambda_v = (Kr / dr) * tc
    return lambda_v


def get_theta_r(Sy, w_bkf, B):
    """
    Calculate the river-aquifer geometry parameter (theta_r or pi2).

    This parameter represents the geometric relationship between the river
    width and the aquifer storage capacity.

    Parameters
    ----------
    Sy : float
        Specific yield of the aquifer (dimensionless).
    w_bkf : float
        Bankfull width (m).
    B : float
        Aquifer lateral extension (m).

    Returns
    -------
    float
        Theta_r (dimensionless).
    """
    theta_r = w_bkf / (B * Sy)
    return theta_r


def get_chi(S, d_bkf, n, B, tc):
    """
    Calculate the flow advection parameter (chi or pi3).

    This parameter relates the hydraulic characteristics of the river
    flow to the aquifer response.

    Parameters
    ----------
    S : float
        Channel slope (dimensionless).
    d_bkf : float
        Bankfull depth (m).
    n : float
        Manning's roughness coefficient.
    B : float
        Aquifer lateral extension (m).
    tc : float
        Hydraulic time in seconds.

    Returns
    -------
    float
        Chi (dimensionless).
    """
    chi = (np.sqrt(S) * tc * (d_bkf ** (2/3))) / (B * n)
    return chi


def get_rho_0(Hp, H_0):
    """
    Calculate the flow intensity parameter (rho_0 or pi4).

    This parameter represents the relative magnitude of the flood peak
    compared to the base flow level.

    Parameters
    ----------
    Hp : float
        Peak stage above base level (m).
    H_0 : float
        Base flow level (m).

    Returns
    -------
    float
        Rho_0 (dimensionless).
    """
    rho_0 = Hp / H_0
    return rho_0


def get_zeta(S, L, H_0):
    """
    Calculate the longitudinal gradient parameter (zeta or pi5).

    This parameter represents the ratio between the elevation drop
    along the reach and the base flow depth.

    Parameters
    ----------
    S : float
        Channel slope (dimensionless).
    L : float
        Stream segment length (m).
    H_0 : float
        Base flow level (m).

    Returns
    -------
    float
        Zeta (dimensionless).
    """
    zeta = (S * L) / H_0
    return zeta
