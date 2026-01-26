"""
River-Aquifer Exchange Modeling Package

This package provides tools for modeling and analyzing river-aquifer
interaction during flood events using MODFLOW.

Modules
-------
dimensionless_numbers
    Dimensionless parameters (Pi numbers) for river-aquifer exchange analysis.
geomorphic_scaling
    Empirical scaling relationships for river and aquifer geometry.
modflow_utils
    Utilities for MODFLOW model setup and results analysis.
"""

from .dimensionless_numbers import (
    get_tc,
    get_gamma_d,
    get_lambda_v,
    get_theta_r,
    get_chi,
    get_rho_0,
    get_zeta,
)

from .geomorphic_scaling import (
    get_w_bkf,
    get_d_bkf,
    get_S,
    get_L,
    get_dx,
    get_dr,
    get_D,
    get_B,
    get_A,
    compute_conductance,
    get_realistic_K_range,
    get_K,
)

from .modflow_utils import (
    generate_topo,
    generate_river,
    generate_river_kinoshita,
    generate_river_stress_period_data,
    hydrograph_cooper_rorabaugh,
    extract_sfr_data,
    analyze_peak_integrals,
    compute_exchange_metrics,
    plot_exchange_flow,
    plot_sfr_results,
    plot_sfr_data,
    plot_results_hds_cbc,
)

__version__ = "1.0.0"
__author__ = "Isabel (PhD Candidate)"
