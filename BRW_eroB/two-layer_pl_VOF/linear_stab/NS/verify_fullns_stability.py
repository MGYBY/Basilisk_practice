#!/usr/bin/env python3
"""Regression and mathematical checks for the readable full-NS solver."""
from __future__ import annotations

from pathlib import Path
import math
import sys

import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import eig

from fullns_stability import (
    CaseParameters,
    FluidLayer,
    SolverOptions,
    apparent_and_tangent_viscosity,
    assemble_phase_speed_pencil,
    build_base_flow,
    build_collocation_layer,
    chebyshev_grid,
    convergence_at_wavenumber,
    equilibrate_rows,
    generalized_residual,
    select_long_wave_pair,
    solve_spectrum,
)


RESULTS: list[tuple[str, bool, float, float]] = []


def check(name: str, error: float, tolerance: float) -> None:
    passed = bool(error <= tolerance)
    RESULTS.append((name, passed, float(error), float(tolerance)))
    if not passed:
        raise AssertionError(f"{name}: error={error:.6e} > tolerance={tolerance:.6e}")


# 1. Chebyshev differentiation on a mapped interval.
grid = chebyshev_grid(32, -0.3, 1.7)
z = grid.z
polynomial = z**4 - 2.0 * z**3 + 0.7 * z**2 + 3.0 * z - 1.0
first_exact = 4.0 * z**3 - 6.0 * z**2 + 1.4 * z + 3.0
second_exact = 12.0 * z**2 - 12.0 * z + 1.4
check(
    "Chebyshev first derivative",
    np.max(np.abs(grid.D1 @ polynomial - first_exact)),
    2.0e-9,
)
check(
    "Chebyshev second derivative",
    np.max(np.abs(grid.D2 @ polynomial - second_exact)),
    2.0e-7,
)


# 2. Default base-flow normalization and gravity traction balance.
case = CaseParameters()
base = build_base_flow(case)
check(
    "Lower mean velocity normalization",
    abs(float(trapezoid(base.lower.velocity, base.lower.z)) - 1.0),
    2.0e-8,
)

mu_lower, _ = apparent_and_tangent_viscosity(
    base.lower.shear,
    base.consistency_number_lower,
    case.lower.flow_index,
    base.epsilon_lower,
)
traction_lower = mu_lower * base.lower.shear
traction_lower_exact = (
    case.slope_tan
    / base.froude**2
    * (1.0 - base.lower.z + base.density_ratio * base.depth_ratio)
)
check(
    "Lower base traction",
    np.max(np.abs(traction_lower - traction_lower_exact)),
    2.0e-7,
)

mu_upper, _ = apparent_and_tangent_viscosity(
    base.upper.shear,
    base.consistency_number_upper,
    case.upper.flow_index,
    base.epsilon_upper,
)
traction_upper = mu_upper * base.upper.shear
traction_upper_exact = (
    base.density_ratio
    * case.slope_tan
    / base.froude**2
    * (1.0 + base.depth_ratio - base.upper.z)
)
check(
    "Upper base traction",
    np.max(np.abs(traction_upper - traction_upper_exact)),
    2.0e-7,
)


# 3. Exact base curvature from dT/dz = mu_t U'' versus spectral derivative.
collocation_lower = build_collocation_layer(base, 48, "lower")
curvature_spectral = np.real(collocation_lower.grid.D1 @ collocation_lower.shear)
interior = slice(2, -2)
check(
    "Base curvature identity",
    np.max(
        np.abs(
            curvature_spectral[interior]
            - collocation_lower.curvature[interior]
        )
    ),
    2.0e-4,
)


# 4. Newtonian reduction of the generalized viscous operator.
newtonian_case = CaseParameters(
    lower=FluidLayer(0.025, 1000.0, 1.0, 1.0),
    upper=FluidLayer(0.025, 1000.0, 1.0, 1.0),
    regularization_rate=0.0,
)
newtonian_base = build_base_flow(newtonian_case)
layer = build_collocation_layer(newtonian_base, 36, "lower")
k = 0.37
L_plus = layer.grid.D2 + k * k * layer.grid.identity
mu_matrix = np.diag(layer.apparent_viscosity)
mu_t_matrix = np.diag(layer.tangent_viscosity)
viscous = (
    L_plus @ mu_t_matrix @ L_plus
    - 4.0 * k * k * layer.grid.D1 @ mu_matrix @ layer.grid.D1
)
expected = (
    layer.apparent_viscosity[0]
    * (layer.grid.D2 - k * k * layer.grid.identity)
    @ (layer.grid.D2 - k * k * layer.grid.identity)
)
check(
    "Newtonian Orr-Sommerfeld reduction",
    np.linalg.norm(viscous - expected) / np.linalg.norm(expected),
    5.0e-10,
)


# 5. Equal-fluid passive-interface mode and total-film surface speed.
options = SolverOptions(
    chebyshev_order=64,
    residual_tolerance=1.0e-5,
    seed_wavenumber=5.0e-4,
)
modes, _, _ = solve_spectrum(newtonian_base, 0.05, options)
passive = min(
    modes,
    key=lambda mode: abs(mode.phase_speed - newtonian_base.interface_velocity),
)
check(
    "Equal-fluid passive-interface phase speed",
    abs(passive.phase_speed - newtonian_base.interface_velocity),
    2.0e-7,
)

modes, _, _ = solve_spectrum(newtonian_base, 0.005, options)
pair = select_long_wave_pair(modes)
check(
    "Equal-fluid long-wave surface celerity",
    abs(pair["surface_zigzag"].celerity - 4.8),
    0.12,
)


# 6. Row equilibration leaves the exact pencil unchanged and stabilizes the
# numerical solution.  Compare a moderate-wavenumber hydrodynamic root.
pencil = assemble_phase_speed_pencil(base, 0.1, 40)
A_scaled, B_scaled, scales = equilibrate_rows(pencil.A, pencil.B)
check(
    "Row scaling is nonzero",
    float(np.max(scales <= 0.0)),
    0.0,
)

values_raw, vectors_raw = eig(pencil.A, pencil.B, right=True, check_finite=False)
values_scaled, vectors_scaled = eig(A_scaled, B_scaled, right=True, check_finite=False)
finite_raw = values_raw[np.isfinite(values_raw.real) & np.isfinite(values_raw.imag)]
finite_scaled = values_scaled[
    np.isfinite(values_scaled.real) & np.isfinite(values_scaled.imag)
]
# Select the surface root near c ~ 2.7 in both spectra.
raw_root = min(finite_raw, key=lambda value: abs(value - (2.7 + 0.1j)))
scaled_root = min(finite_scaled, key=lambda value: abs(value - raw_root))
check(
    "Row-equilibrated eigenvalue invariance at moderate k",
    abs(raw_root - scaled_root),
    2.0e-5,
)


# 7. Long-wave convergence of the weak interface branch and finite-k
# convergence of the surface branch.
longwave_rows = convergence_at_wavenumber(
    case,
    0.01,
    orders=(32, 40, 52, 64),
)
interface_growth = np.array(
    [row["interface_varicose_growth"] for row in longwave_rows]
)
check(
    "Long-wave interface growth convergence",
    np.max(np.abs(interface_growth[1:] - interface_growth[-1])),
    3.0e-6,
)

finite_rows = convergence_at_wavenumber(
    case,
    0.4,
    orders=(32, 40, 52, 64),
)
surface_growth = np.array(
    [row["surface_zigzag_growth"] for row in finite_rows]
)
check(
    "Finite-k surface growth convergence",
    np.max(np.abs(surface_growth[1:] - surface_growth[-1])),
    7.0e-5,
)


# 8. Direct residual check on the selected default roots.
modes, pencil, _ = solve_spectrum(base, 0.01, options)
pair = select_long_wave_pair(modes)
for name, mode in pair.items():
    check(
        f"Original-pencil residual: {name}",
        generalized_residual(
            pencil.A,
            pencil.B,
            mode.phase_speed,
            mode.vector,
        ),
        1.0e-5,
    )


output = Path(__file__).resolve().parents[1] / "data" / "verification.txt"
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w", encoding="utf-8") as handle:
    for name, passed, error, tolerance in RESULTS:
        handle.write(
            f"{'PASS' if passed else 'FAIL'}\t{name}\t"
            f"error={error:.12e}\ttolerance={tolerance:.12e}\n"
        )
    handle.write("\nAll verification checks passed.\n")

print(output.read_text(encoding="utf-8"))
