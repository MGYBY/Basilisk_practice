#!/usr/bin/env python3
"""User-editable settings for the two-layer power-law stability workflow.

Edit this file only for normal use.  The numerical/model files do not need to be
changed for ordinary parameter sweeps or dispersion-curve calculations.
"""

from pathlib import Path

import numpy as np

from two_layer_powerlaw_stability_core import Parameters


# =============================================================================
# A. BASE CASE USED FOR FINITE-WAVENUMBER DISPERSION TABLES
# =============================================================================

# These are dimensionless parameters of the normalized two-layer model.
BASE_PARAMETERS = Parameters(
    froude=0.8,                    # lower-layer Froude number Fr_l
    depth_ratio=1.0,               # h_r = H_u/H_l
    density_ratio=0.9,             # r_rho = rho_u/rho_l
    n_lower=0.6,                   # lower-layer power-law index n_l
    n_upper=0.8,                   # upper-layer power-law index n_u
    scaled_consistency_ratio=1.0,  # kappa_K = Lambda_u/Lambda_l
)

# A short label used to name the output directory.
CASE_NAME = "default_case"


# =============================================================================
# B. CRITICAL-FROUDE PARAMETER SWEEP (k -> 0+)
# =============================================================================

# Allowed names and common aliases:
#   "n_lower" or "n_l"
#   "depth_ratio" or "h_r"
#   "density_ratio" or "r_rho"
#   "n_upper" or "n_u"
#   "scaled_consistency_ratio" or "kappa_K"
SWEEP_PARAMETER = "n_lower"
SWEEP_VALUES = np.linspace(0.2, 1.0, 81)

# Usually the roll-wave threshold of interest is the free-surface branch.
# The alternative is "interfacial".
CRITICAL_BRANCH = "free_surface"

# Search interval used for the direct solution of d_j(Fr_l)=0.
FROUDE_SEARCH_MIN = 0.05
FROUDE_SEARCH_MAX = 3.0
FROUDE_SCAN_POINTS = 320
FROUDE_ROOT_TOLERANCE = 1.0e-8


# =============================================================================
# C. WAVENUMBER GRID FOR THE FOUR DISPERSION BRANCHES
# =============================================================================

WAVENUMBER_MIN = 1.0e-4
WAVENUMBER_MAX = 1.0e2
WAVENUMBER_POINTS = 101

# "log" is recommended when both the long-wave and finite-k ranges matter.
# "linear" is useful for a narrow prescribed interval.
# WAVENUMBER_SPACING = "log"
WAVENUMBER_SPACING = "linear"


# =============================================================================
# D. OUTPUT LOCATION
# =============================================================================

OUTPUT_ROOT = Path(__file__).resolve().parent / "results"
