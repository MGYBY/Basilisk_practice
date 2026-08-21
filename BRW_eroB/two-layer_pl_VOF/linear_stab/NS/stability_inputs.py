#!/usr/bin/env python3
"""Central user inputs for the full-NS threshold calculations.

Edit this file first.  The other scripts import the values below.
"""
from __future__ import annotations

import numpy as np

from fullns_stability import SolverOptions
from fullns_thresholds import (
    DimensionlessFamily,
    LongWaveOptions,
    NeutralCurveOptions,
    WaveSpecification,
)


# =============================================================================
# DIMENSIONLESS NORMAL-FLOW FAMILY
# =============================================================================

FAMILY = DimensionlessFamily(
    depth_ratio=1.0,             # h_r = H_u/H_l
    density_ratio=0.8,           # r_rho = rho_u/rho_l
    n_lower=0.4,
    n_upper=0.8,
    consistency_ratio=1.7336759276597786,  # kappa_K = Lambda_u/Lambda_l
    epsilon_lower=0.02233444771394924,
    epsilon_upper=0.02233444771394924,
)

# The regularized problem holds epsilon_i fixed as dimensionless transition
# shear rates.  Setting epsilon_upper=0 with n_upper<1 is not admissible in the
# full-NS eigenproblem because the ideal viscosity diverges at the free surface.


# =============================================================================
# SPECTRAL EIGENVALUE SOLVER
# =============================================================================

SOLVER = SolverOptions(
    chebyshev_order=44,
    residual_tolerance=2.0e-5,
    phase_speed_limit=100.0,
    seed_wavenumber=1.0e-4,
    base_samples=12001,
)

# Neutral maps require many spectra.  Order 32 is used for the production
# curves, while selected points are rechecked at higher order in the
# verification script.
NEUTRAL_SOLVER = SolverOptions(
    chebyshev_order=32,
    residual_tolerance=7.0e-5,
    phase_speed_limit=100.0,
    seed_wavenumber=1.0e-4,
    base_samples=6001,
)


# =============================================================================
# INFINITELY LONG-WAVE CRITICAL FROUDE NUMBER
# =============================================================================

LONG_WAVE_BRANCHES = (
    "surface_zigzag",
    "interface_varicose",
)

LONG_WAVE = LongWaveOptions(
    froude_min=0.10,
    froude_max=1.60,
    scan_points=25,
    reference_reynolds=20.0,
    wavenumbers=(3.0e-4, 5.0e-4, 8.0e-4, 1.2e-3, 2.0e-3),
    polynomial_degree=1,
)


# =============================================================================
# FINITE-WAVELENGTH NEUTRAL CURVES
# =============================================================================

# Exact analogue of the wavelength convention used in Figure 2 of Yu & Chu
# (2024): value = S0*lambda/H_l.  The script can also use
#   WaveSpecification("fixed_kH", k_H)
# or
#   WaveSpecification("fixed_lambda_over_H", lambda/H_l).
WAVE_SPECIFICATIONS = (
    WaveSpecification("fixed_slope_scaled_wavelength", 4.72),
    WaveSpecification("fixed_slope_scaled_wavelength", 0.472),
)

# The production Figure-2-style maps show the robust surface/zigzag neutral
# boundary.  Add "interface_varicose" here for an exploratory interface map;
# for the present family that weak branch has no converged positive long-wave
# threshold and is much more sensitive to regularization and spectral order.
NEUTRAL_BRANCHES = (
    "surface_zigzag",
)

# A logarithmic Reynolds grid gives a readable Fr-Re stability map.
REYNOLDS_VALUES = np.geomspace(2.5, 300.0, 12)

NEUTRAL_CURVE = NeutralCurveOptions(
    froude_min=0.10,
    froude_max=1.60,
    scan_points=16,
    root_xtol=1.0e-6,
    continuation_steps=6,
)
