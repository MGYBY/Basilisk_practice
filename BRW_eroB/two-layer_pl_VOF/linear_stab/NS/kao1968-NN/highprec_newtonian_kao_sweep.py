#!/usr/bin/env python3
"""High-precision Kao-type long-wave stability sweeps for two Newtonian layers.

This is a generalized, user-facing revision of ``highprec_newtonian.py``.
It calculates the infinitely-long-wave, branch-specific lower-layer Froude
thresholds for two immiscible Newtonian liquid layers flowing down an incline.
The two hydrodynamic branches are

* ``surface``: the faster, free-surface-dominated (zigzag) mode;
* ``interfacial``: the slower, liquid-interface-dominated (varicose) mode.

The calculation uses the long-wave expansion

    sigma_j(k) = -i*c_j*k + d_j(Fr_l)*k**2 + O(k**3),

with

    d_j(Fr_l) = A_j + B_j/Fr_l**2.

A positive ``d_j`` denotes temporal growth.  A positive neutral Froude number
exists only when ``A_j`` and ``B_j`` have opposite signs.  When no positive
root exists, the output states whether the branch is stable or unstable for
all positive Froude numbers.

The code uses an exact degree-six polynomial representation of the first
long-wave correction for Newtonian base profiles and high-precision ``mpmath``
linear algebra.  It computes results from the Kao-type full Navier--Stokes
long-wave formulation; it does not digitize values from Kao's published plots.

Dependencies
------------
    python -m pip install mpmath

Quick use
---------
1. Edit the USER SETTINGS block below and run

       python highprec_newtonian_kao_sweep.py

2. Or override settings from the command line, for example

       python highprec_newtonian_kao_sweep.py \
           --sweep viscosity_ratio --mode both \
           --minimum 0.01 --maximum 100 --points 41 --spacing log

The output is a self-explanatory tab-separated text file with comment lines
beginning with ``#``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import mpmath as mp


# =============================================================================
# USER SETTINGS
# =============================================================================

# Parameter to sweep: "depth_ratio", "density_ratio", or "viscosity_ratio".
SWEEP_PARAMETER = "depth_ratio"

# Mode to report: "surface", "interfacial", or "both".
MODE = "surface"

# Fixed parameters.  The value corresponding to SWEEP_PARAMETER is replaced by
# the sweep values.  Strings preserve the entered decimal values exactly.
FIXED_DEPTH_RATIO = "1.0"          # h_r = H_u/H_l
FIXED_DENSITY_RATIO = "0.8"        # r_rho = rho_u/rho_l
FIXED_VISCOSITY_RATIO = "1.0"      # mu_ratio = mu_u/mu_l

# Default sweep definitions.  Logarithmic spacing is useful for depth and
# viscosity ratios; linear spacing is convenient for the density ratio.
SWEEP_DEFINITIONS = {
    "depth_ratio": {
        "minimum": "0.01",
        "maximum": "100.0",
        "points": 41,
        "spacing": "log",
    },
    "density_ratio": {
        "minimum": "0.10",
        "maximum": "0.99",
        "points": 45,
        "spacing": "linear",
    },
    "viscosity_ratio": {
        "minimum": "0.01",
        "maximum": "100.0",
        "points": 41,
        "spacing": "log",
    },
}

# Optional explicit sweep values.  Set to a sequence of strings, for example
# ["0.01", "0.1", "1", "10", "100"], or leave as None.
CUSTOM_SWEEP_VALUES: Optional[Sequence[str]] = None

# Internal and printed precision.
MP_DECIMAL_DIGITS = 80
OUTPUT_SIGNIFICANT_DIGITS = 18

# None gives an automatic filename based on the selected sweep and mode.
OUTPUT_FILE: Optional[str] = None

# The checks reproduce the three equal-viscosity values from the original
# highprec_newtonian.py before the requested sweep begins.
RUN_LEGACY_REGRESSION_CHECKS = True


# =============================================================================
# Data structures and input utilities
# =============================================================================


@dataclass(frozen=True)
class Parameters:
    """Dimensionless parameters of the two-layer Newtonian normal flow."""

    depth_ratio: mp.mpf
    density_ratio: mp.mpf
    viscosity_ratio: mp.mpf

    def validate(self) -> None:
        values = {
            "depth_ratio h_r": self.depth_ratio,
            "density_ratio r_rho": self.density_ratio,
            "viscosity_ratio mu_u/mu_l": self.viscosity_ratio,
        }
        for name, value in values.items():
            if not mp.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite; got {value}")


@dataclass(frozen=True)
class LeadingMode:
    """Leading-order k -> 0 hydrodynamic mode."""

    label: str
    celerity: mp.mpf
    surface_to_interface: mp.mpf
    c_minus_interface_velocity: mp.mpf
    leading_residual: mp.mpf


@dataclass(frozen=True)
class CorrectionResult:
    """First long-wave correction at one Froude number."""

    growth_coefficient: mp.mpf
    relative_residual: mp.mpf


@dataclass(frozen=True)
class ModeResult:
    """Complete branch-specific result at one parameter point."""

    sweep_value: mp.mpf
    parameters: Parameters
    mode: str
    celerity: mp.mpf
    surface_to_interface: mp.mpf
    coefficient_at_infinite_froude: mp.mpf
    coefficient_of_inverse_froude_squared: mp.mpf
    d_at_froude_1: mp.mpf
    d_at_froude_2: mp.mpf
    critical_froude: mp.mpf
    onset_froude: mp.mpf
    stability_class: str
    leading_residual: mp.mpf
    correction_residual: mp.mpf
    neutral_residual: mp.mpf


@dataclass(frozen=True)
class ErrorResult:
    """Recoverable error at one sweep point."""

    sweep_value: mp.mpf
    parameters: Parameters
    mode: str
    error_type: str
    message: str


def _as_mpf(value: object, name: str) -> mp.mpf:
    try:
        result = mp.mpf(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Could not parse {name}={value!r} as a real number") from exc
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be finite; got {value!r}")
    return result


def _normalize_sweep_name(name: str) -> str:
    key = name.strip().lower().replace("-", "_")
    aliases = {
        "hr": "depth_ratio",
        "h_r": "depth_ratio",
        "depth": "depth_ratio",
        "depth_ratio": "depth_ratio",
        "rrho": "density_ratio",
        "r_rho": "density_ratio",
        "rho": "density_ratio",
        "density": "density_ratio",
        "density_ratio": "density_ratio",
        "mu": "viscosity_ratio",
        "m": "viscosity_ratio",
        "mu_ratio": "viscosity_ratio",
        "viscosity": "viscosity_ratio",
        "viscosity_ratio": "viscosity_ratio",
    }
    try:
        return aliases[key]
    except KeyError as exc:
        allowed = "depth_ratio, density_ratio, viscosity_ratio"
        raise ValueError(f"Unknown sweep parameter {name!r}; choose {allowed}") from exc


def _normalize_mode_name(name: str) -> str:
    key = name.strip().lower().replace("-", "_")
    aliases = {
        "surface": "surface",
        "surface_mode": "surface",
        "free_surface": "surface",
        "zigzag": "surface",
        "interface": "interfacial",
        "interfacial": "interfacial",
        "interfacial_mode": "interfacial",
        "varicose": "interfacial",
        "both": "both",
        "all": "both",
    }
    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(
            f"Unknown mode {name!r}; choose surface, interfacial, or both"
        ) from exc


def _generate_sweep_values(
    minimum: mp.mpf,
    maximum: mp.mpf,
    points: int,
    spacing: str,
) -> List[mp.mpf]:
    if points < 1:
        raise ValueError("points must be at least 1")
    if minimum <= 0 or maximum <= 0:
        raise ValueError("sweep endpoints must be positive")
    if maximum < minimum:
        raise ValueError("maximum must be greater than or equal to minimum")
    if points == 1:
        return [minimum]

    spacing_key = spacing.strip().lower()
    if spacing_key == "linear":
        return [
            minimum + (maximum - minimum) * mp.mpf(i) / mp.mpf(points - 1)
            for i in range(points)
        ]
    if spacing_key == "log":
        ratio = maximum / minimum
        return [
            minimum * mp.power(ratio, mp.mpf(i) / mp.mpf(points - 1))
            for i in range(points)
        ]
    raise ValueError("spacing must be 'linear' or 'log'")


def _parse_explicit_values(text: str) -> List[mp.mpf]:
    values = [token.strip() for token in text.split(",") if token.strip()]
    if not values:
        raise ValueError("--values must contain at least one comma-separated value")
    parsed = [_as_mpf(value, "sweep value") for value in values]
    if any(value <= 0 for value in parsed):
        raise ValueError("all explicit sweep values must be positive")
    return parsed


def _parameters_for_sweep_value(
    sweep_parameter: str,
    sweep_value: mp.mpf,
    fixed_depth_ratio: mp.mpf,
    fixed_density_ratio: mp.mpf,
    fixed_viscosity_ratio: mp.mpf,
) -> Parameters:
    values = {
        "depth_ratio": fixed_depth_ratio,
        "density_ratio": fixed_density_ratio,
        "viscosity_ratio": fixed_viscosity_ratio,
    }
    values[sweep_parameter] = sweep_value
    parameters = Parameters(**values)
    parameters.validate()
    return parameters


# =============================================================================
# Exact polynomial long-wave formulation
# =============================================================================


POLYNOMIAL_DEGREE = 6


def _polynomial_evaluation_row(
    degree: int,
    coordinate: mp.mpf,
    derivative: int = 0,
) -> List[mp.mpf]:
    row: List[mp.mpf] = []
    for power in range(degree + 1):
        if power < derivative:
            row.append(mp.mpf("0"))
        else:
            factor = mp.factorial(power) / mp.factorial(power - derivative)
            row.append(factor * coordinate ** (power - derivative))
    return row


def _multiply_polynomial_by_derivative_matrix(
    polynomial_coefficients: Sequence[mp.mpf],
    derivative: int,
    degree: int,
) -> List[List[mp.mpf]]:
    """Coefficient matrix of U(x) * d^derivative(x^j)/dx^derivative."""

    matrix = [
        [mp.mpf("0") for _ in range(degree + 1)]
        for _ in range(degree + 1)
    ]
    for column in range(degree + 1):
        if column < derivative:
            continue
        factor = mp.factorial(column) / mp.factorial(column - derivative)
        base_power = column - derivative
        for offset, coefficient in enumerate(polynomial_coefficients):
            row = base_power + offset
            if row <= degree:
                matrix[row][column] += coefficient * factor
    return matrix


def _base_quantities(parameters: Parameters) -> Tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Return chi, interface velocity, and free-surface velocity."""

    h = parameters.depth_ratio
    r = parameters.density_ratio
    viscosity_ratio = parameters.viscosity_ratio

    # The lower mean velocity is normalized to one.
    chi = 1 / (mp.mpf(1) / 3 + r * h / 2)
    interface_velocity = chi * (mp.mpf("0.5") + r * h)
    surface_velocity = (
        interface_velocity + (r / viscosity_ratio) * chi * h * h / 2
    )
    return chi, interface_velocity, surface_velocity


def _assemble_longwave_matrices(
    parameters: Parameters,
    froude: mp.mpf,
    degree: int = POLYNOMIAL_DEGREE,
) -> Tuple[mp.matrix, mp.matrix, mp.matrix, mp.matrix]:
    """Assemble A0, B0 and their first-wavenumber correction matrices.

    The generalized pencil is expanded in the small wavenumber.  ``A0`` and
    ``B0`` determine the leading celerities; ``A1`` and ``B1`` determine the
    O(k^2) growth coefficient.  The lower Reynolds normalization is one.  This
    changes the magnitude of d but not its neutral Froude number.
    """

    parameters.validate()
    if froude <= 0 or not mp.isfinite(froude):
        raise ValueError("froude must be positive and finite")

    h = parameters.depth_ratio
    r = parameters.density_ratio
    viscosity_ratio = parameters.viscosity_ratio
    chi, interface_velocity, _ = _base_quantities(parameters)

    lower_velocity = [
        mp.mpf("0"),
        chi * (1 + r * h),
        -chi / 2,
    ]
    upper_velocity = [
        interface_velocity,
        (r / viscosity_ratio) * chi * h,
        -(r / viscosity_ratio) * chi / 2,
    ]
    lower_velocity_derivative = [lower_velocity[1], 2 * lower_velocity[2]]
    upper_velocity_derivative = [upper_velocity[1], 2 * upper_velocity[2]]
    lower_velocity_second = 2 * lower_velocity[2]
    upper_velocity_second = 2 * upper_velocity[2]

    unknowns = 2 * (degree + 1) + 2
    A0 = mp.matrix(unknowns, unknowns)
    B0 = mp.matrix(unknowns, unknowns)
    A1 = mp.matrix(unknowns, unknowns)
    B1 = mp.matrix(unknowns, unknowns)

    eta_index = 2 * (degree + 1)
    zeta_index = eta_index + 1
    row = 0

    identity_polynomial = [
        [mp.mpf(int(i == j)) for j in range(degree + 1)]
        for i in range(degree + 1)
    ]
    second_derivative_matrix = _multiply_polynomial_by_derivative_matrix(
        [mp.mpf(1)], 2, degree
    )

    layer_data = (
        (
            0,
            mp.mpf(1),
            mp.mpf(1),
            lower_velocity,
            lower_velocity_second,
        ),
        (
            degree + 1,
            viscosity_ratio,
            r,
            upper_velocity,
            upper_velocity_second,
        ),
    )

    # Interior equations.  At the first correction, the exact Newtonian
    # polynomial has degree six, so degree=6 closes the system exactly.
    for offset, viscosity, density, velocity, velocity_second in layer_data:
        velocity_times_second = _multiply_polynomial_by_derivative_matrix(
            velocity, 2, degree
        )
        for power in range(degree - 3):
            fourth_power = power + 4
            A0[row, offset + fourth_power] = (
                viscosity
                * mp.factorial(fourth_power)
                / mp.factorial(power)
            )
            for column in range(degree + 1):
                A1[row, offset + column] = density * (
                    -velocity_times_second[power][column]
                    + velocity_second * identity_polynomial[power][column]
                )
                B1[row, offset + column] = (
                    -density * second_derivative_matrix[power][column]
                )
            row += 1

    # Rigid wall: psi_l=0 and psi_l'=0.
    A0[row, 0] = 1
    row += 1
    A0[row, 1] = 1
    row += 1

    # Interface continuity of normal velocity/streamfunction.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, mp.mpf(1), 0)
    ):
        A0[row, column] += value
    A0[row, degree + 1] -= 1
    row += 1

    # Interface continuity of tangential velocity, including displacement.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, mp.mpf(1), 1)
    ):
        A0[row, column] += value
    A0[row, degree + 2] -= 1
    lower_shear_at_interface = (
        lower_velocity_derivative[0] + lower_velocity_derivative[1]
    )
    upper_shear_at_interface = upper_velocity_derivative[0]
    A0[row, eta_index] = lower_shear_at_interface - upper_shear_at_interface
    row += 1

    # Interface tangential traction.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, mp.mpf(1), 2)
    ):
        A0[row, column] += value
    A0[row, degree + 3] -= 2 * viscosity_ratio
    A0[row, eta_index] -= (1 - r) * chi
    row += 1

    # Interface normal traction.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, mp.mpf(1), 3)
    ):
        A0[row, column] -= value
    A0[row, degree + 4] += 6 * viscosity_ratio

    lower_velocity_at_interface = sum(lower_velocity)
    upper_velocity_at_interface = upper_velocity[0]
    lower_velocity_prime_at_interface = lower_shear_at_interface
    upper_velocity_prime_at_interface = upper_shear_at_interface

    lower_prime_row = _polynomial_evaluation_row(degree, mp.mpf(1), 1)
    lower_value_row = _polynomial_evaluation_row(degree, mp.mpf(1), 0)
    upper_prime_row = _polynomial_evaluation_row(degree, mp.mpf(0), 1)
    upper_value_row = _polynomial_evaluation_row(degree, mp.mpf(0), 0)

    for column in range(degree + 1):
        A1[row, column] += (
            lower_velocity_at_interface * lower_prime_row[column]
            - lower_velocity_prime_at_interface * lower_value_row[column]
        )
        A1[row, degree + 1 + column] -= r * (
            upper_velocity_at_interface * upper_prime_row[column]
            - upper_velocity_prime_at_interface * upper_value_row[column]
        )
    A1[row, eta_index] += (1 - r) / (froude * froude)
    for column, value in enumerate(lower_prime_row):
        B1[row, column] += value
    for column, value in enumerate(upper_prime_row):
        B1[row, degree + 1 + column] -= r * value
    row += 1

    # Stress-free free surface: tangential traction.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, h, 2)
    ):
        A0[row, degree + 1 + column] += viscosity_ratio * value
    A0[row, zeta_index] -= r * chi
    row += 1

    # Stress-free free surface: normal traction.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, h, 3)
    ):
        A0[row, degree + 1 + column] -= viscosity_ratio * value

    upper_velocity_at_surface = sum(
        upper_velocity[power] * h**power
        for power in range(len(upper_velocity))
    )
    upper_velocity_prime_at_surface = (
        upper_velocity_derivative[0] + upper_velocity_derivative[1] * h
    )
    upper_surface_prime_row = _polynomial_evaluation_row(degree, h, 1)
    upper_surface_value_row = _polynomial_evaluation_row(degree, h, 0)
    for column in range(degree + 1):
        A1[row, degree + 1 + column] += r * (
            upper_velocity_at_surface * upper_surface_prime_row[column]
            - upper_velocity_prime_at_surface * upper_surface_value_row[column]
        )
    A1[row, zeta_index] += r / (froude * froude)
    for column, value in enumerate(upper_surface_prime_row):
        B1[row, degree + 1 + column] += r * value
    row += 1

    # Interface kinematic condition.
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, mp.mpf(1), 0)
    ):
        A0[row, column] += value
    A0[row, eta_index] += interface_velocity
    B0[row, eta_index] = 1
    row += 1

    # Free-surface kinematic condition.
    _, _, surface_velocity = _base_quantities(parameters)
    for column, value in enumerate(
        _polynomial_evaluation_row(degree, h, 0)
    ):
        A0[row, degree + 1 + column] += value
    A0[row, zeta_index] += surface_velocity
    B0[row, zeta_index] = 1
    row += 1

    if row != unknowns:
        raise RuntimeError(
            f"Internal equation-count error: assembled {row}, expected {unknowns}"
        )

    return A0, B0, A1, B1


def _leading_quadratic_coefficients(
    parameters: Parameters,
) -> Tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """Return the quadratic coefficients for x=c-U_I and base quantities."""

    parameters.validate()
    h = parameters.depth_ratio
    r = parameters.density_ratio
    viscosity_ratio = parameters.viscosity_ratio
    chi, interface_velocity, _ = _base_quantities(parameters)

    # The free-surface kinematic condition gives
    # A*x^2 + B*x + C = 0, where x=c-U_I.
    A = -2 / (chi * r)
    B = 2 * h * h / viscosity_ratio + 2 * h + 1 / r
    C = chi * h * h * (r - 1 / viscosity_ratio)
    return A, B, C, chi, interface_velocity


def _stable_quadratic_roots(
    A: mp.mpf,
    B: mp.mpf,
    C: mp.mpf,
) -> Tuple[mp.mpf, mp.mpf]:
    discriminant = B * B - 4 * A * C
    scale = max(abs(B * B), abs(4 * A * C), mp.mpf(1))
    tolerance = scale * mp.power(10, -max(20, mp.mp.dps // 2))
    if discriminant < -tolerance:
        raise ValueError(
            "The two leading long-wave celerities are complex for this "
            f"parameter set (discriminant={mp.nstr(discriminant, 12)})."
        )
    if discriminant < 0:
        discriminant = mp.mpf(0)

    square_root = mp.sqrt(discriminant)
    sign_B = mp.mpf(1) if B >= 0 else mp.mpf(-1)
    q = -(B + sign_B * square_root) / 2
    if q == 0:
        # Only a degenerate quadratic reaches this branch.
        root_1 = -B / (2 * A)
        root_2 = root_1
    else:
        root_1 = q / A
        root_2 = C / q
    return root_1, root_2


def _leading_vector(
    parameters: Parameters,
    celerity: mp.mpf,
    degree: int = POLYNOMIAL_DEGREE,
) -> mp.matrix:
    """Construct the exact leading-order polynomial eigenvector with eta=1."""

    h = parameters.depth_ratio
    r = parameters.density_ratio
    viscosity_ratio = parameters.viscosity_ratio
    chi, interface_velocity, _ = _base_quantities(parameters)
    x = celerity - interface_velocity

    lower_a2 = x
    upper_b0 = x
    upper_b1 = 2 * x + r * chi * h * (1 - 1 / viscosity_ratio)
    upper_b2 = (x - (1 - r) * chi / 2) / viscosity_ratio
    surface_to_interface = (2 * x - (1 - r) * chi) / (r * chi)

    unknowns = 2 * (degree + 1) + 2
    vector = mp.matrix(unknowns, 1)
    vector[2] = lower_a2
    vector[degree + 1] = upper_b0
    vector[degree + 2] = upper_b1
    vector[degree + 3] = upper_b2
    vector[2 * (degree + 1)] = 1
    vector[2 * (degree + 1) + 1] = surface_to_interface
    return vector


def _vector_max_abs(vector: mp.matrix) -> mp.mpf:
    return max((abs(vector[row]) for row in range(vector.rows)), default=mp.mpf(0))


def _relative_generalized_residual(
    A: mp.matrix,
    B: mp.matrix,
    eigenvalue: mp.mpf,
    vector: mp.matrix,
) -> mp.mpf:
    residual = A * vector - eigenvalue * (B * vector)
    scale = max(
        _vector_max_abs(A * vector),
        abs(eigenvalue) * _vector_max_abs(B * vector),
        mp.mpf(1),
    )
    return _vector_max_abs(residual) / scale


def leading_modes(parameters: Parameters) -> Dict[str, LeadingMode]:
    """Return the interface- and surface-dominated k -> 0 modes."""

    A, B, C, chi, interface_velocity = _leading_quadratic_coefficients(parameters)
    x_roots = _stable_quadratic_roots(A, B, C)

    A0, B0, _, _ = _assemble_longwave_matrices(
        parameters, mp.mpf(1), POLYNOMIAL_DEGREE
    )
    records: List[LeadingMode] = []
    for x in x_roots:
        celerity = interface_velocity + x
        surface_to_interface = (
            2 * x - (1 - parameters.density_ratio) * chi
        ) / (parameters.density_ratio * chi)
        vector = _leading_vector(parameters, celerity, POLYNOMIAL_DEGREE)
        residual = _relative_generalized_residual(A0, B0, celerity, vector)
        records.append(
            LeadingMode(
                label="unlabelled",
                celerity=celerity,
                surface_to_interface=surface_to_interface,
                c_minus_interface_velocity=x,
                leading_residual=residual,
            )
        )

    # Kao's two long-wave branches are labelled by their leading celerities:
    # the slower second mode is interface/varicose, while the faster first mode
    # is free-surface/zigzag.  The displacement ratio is retained as a useful
    # diagnostic but is not a globally reliable classifier at extreme density
    # or viscosity ratios.
    records.sort(key=lambda mode: mode.celerity)
    interface_mode, surface_mode = records[0], records[1]

    separation = abs(surface_mode.celerity - interface_mode.celerity)
    separation_scale = max(
        abs(surface_mode.celerity), abs(interface_mode.celerity), mp.mpf(1)
    )
    if separation <= separation_scale * mp.power(10, -max(20, mp.mp.dps // 2)):
        raise ValueError(
            "The two leading celerities are numerically degenerate, so the "
            "surface/interfacial labels are not uniquely defined."
        )

    return {
        "interfacial": LeadingMode(
            label="interfacial",
            celerity=interface_mode.celerity,
            surface_to_interface=interface_mode.surface_to_interface,
            c_minus_interface_velocity=interface_mode.c_minus_interface_velocity,
            leading_residual=interface_mode.leading_residual,
        ),
        "surface": LeadingMode(
            label="surface",
            celerity=surface_mode.celerity,
            surface_to_interface=surface_mode.surface_to_interface,
            c_minus_interface_velocity=surface_mode.c_minus_interface_velocity,
            leading_residual=surface_mode.leading_residual,
        ),
    }


def _first_longwave_correction(
    parameters: Parameters,
    froude: mp.mpf,
    leading_mode: LeadingMode,
    degree: int = POLYNOMIAL_DEGREE,
) -> CorrectionResult:
    """Calculate d in sigma=-i*c*k+d*k^2+... at one Froude number."""

    A0, B0, A1, B1 = _assemble_longwave_matrices(parameters, froude, degree)
    q0 = _leading_vector(parameters, leading_mode.celerity, degree)
    singular_matrix = A0 - leading_mode.celerity * B0
    right_hand_side = -(A1 - leading_mode.celerity * B1) * q0

    size = singular_matrix.rows
    augmented = mp.matrix(size + 1, size + 1)
    augmented_rhs = mp.matrix(size + 1, 1)
    Bq = B0 * q0

    for row in range(size):
        for column in range(size):
            augmented[row, column] = singular_matrix[row, column]
        augmented[row, size] = -Bq[row]
        augmented_rhs[row] = right_hand_side[row]

    # Fix the arbitrary first-order eigenvector normalization by setting eta1=0.
    eta_index = 2 * (degree + 1)
    augmented[size, eta_index] = 1
    augmented_rhs[size] = 0

    solution = mp.lu_solve(augmented, augmented_rhs)
    growth_coefficient = solution[size]

    residual_vector = augmented * solution - augmented_rhs
    residual_scale = max(
        _vector_max_abs(augmented * solution),
        _vector_max_abs(augmented_rhs),
        mp.mpf(1),
    )
    relative_residual = _vector_max_abs(residual_vector) / residual_scale
    return CorrectionResult(growth_coefficient, relative_residual)


def _coefficient_zero_tolerance(values: Iterable[mp.mpf]) -> mp.mpf:
    values_list = [abs(value) for value in values]
    scale = max(values_list, default=mp.mpf(1))
    if scale == 0:
        scale = mp.mpf(1)
    relative = scale * mp.power(10, -max(20, mp.mp.dps // 2))
    absolute = mp.power(10, -(mp.mp.dps - 10))
    return max(relative, absolute)


def _sign_with_tolerance(value: mp.mpf, tolerance: mp.mpf) -> int:
    if value > tolerance:
        return 1
    if value < -tolerance:
        return -1
    return 0


def _nan() -> mp.mpf:
    return mp.mpf("nan")


def evaluate_mode(
    parameters: Parameters,
    mode_name: str,
    sweep_value: mp.mpf,
) -> ModeResult:
    """Evaluate one branch at one sweep point."""

    modes = leading_modes(parameters)
    mode = modes[mode_name]
    correction_1 = _first_longwave_correction(
        parameters, mp.mpf(1), mode, POLYNOMIAL_DEGREE
    )
    correction_2 = _first_longwave_correction(
        parameters, mp.mpf(2), mode, POLYNOMIAL_DEGREE
    )

    d1 = correction_1.growth_coefficient
    d2 = correction_2.growth_coefficient

    # d(1)=A+B and d(2)=A+B/4.
    coefficient_B = mp.mpf(4) * (d1 - d2) / 3
    coefficient_A = d1 - coefficient_B

    tolerance = _coefficient_zero_tolerance(
        (coefficient_A, coefficient_B, d1, d2)
    )
    sign_A = _sign_with_tolerance(coefficient_A, tolerance)
    sign_B = _sign_with_tolerance(coefficient_B, tolerance)

    critical_froude = _nan()
    onset_froude = _nan()
    neutral_residual = _nan()

    if sign_A * sign_B < 0:
        critical_froude = mp.sqrt(-coefficient_B / coefficient_A)
        neutral_residual = abs(
            coefficient_A
            + coefficient_B / (critical_froude * critical_froude)
        )
        if sign_B < 0 and sign_A > 0:
            stability_class = "stable_below__unstable_above"
            onset_froude = critical_froude
        else:
            stability_class = "unstable_below__stable_above"
            # The branch is already unstable arbitrarily close to Fr=0.
            onset_froude = mp.mpf(0)
    elif sign_A > 0 and sign_B >= 0:
        stability_class = "unstable_for_all_positive_Fr"
        onset_froude = mp.mpf(0)
    elif sign_A < 0 and sign_B <= 0:
        stability_class = "stable_for_all_positive_Fr"
    elif sign_A == 0 and sign_B > 0:
        stability_class = "unstable_for_all_finite_positive_Fr__neutral_as_Fr_inf"
        onset_froude = mp.mpf(0)
    elif sign_A == 0 and sign_B < 0:
        stability_class = "stable_for_all_finite_positive_Fr__neutral_as_Fr_inf"
    elif sign_A > 0 and sign_B == 0:
        stability_class = "unstable_for_all_positive_Fr"
        onset_froude = mp.mpf(0)
    elif sign_A < 0 and sign_B == 0:
        stability_class = "stable_for_all_positive_Fr"
    elif sign_A == 0 and sign_B == 0:
        stability_class = "neutral_to_resolved_precision_for_all_Fr"
    else:
        stability_class = "degenerate__inspect_coefficients"

    return ModeResult(
        sweep_value=sweep_value,
        parameters=parameters,
        mode=mode_name,
        celerity=mode.celerity,
        surface_to_interface=mode.surface_to_interface,
        coefficient_at_infinite_froude=coefficient_A,
        coefficient_of_inverse_froude_squared=coefficient_B,
        d_at_froude_1=d1,
        d_at_froude_2=d2,
        critical_froude=critical_froude,
        onset_froude=onset_froude,
        stability_class=stability_class,
        leading_residual=mode.leading_residual,
        correction_residual=max(
            correction_1.relative_residual,
            correction_2.relative_residual,
        ),
        neutral_residual=neutral_residual,
    )


# =============================================================================
# Regression checks, output, and command-line interface
# =============================================================================


LEGACY_SURFACE_RESULTS = {
    "0.01": "0.52704608621000796464385520072206687115335746347577",
    "1": "0.48681309360276675176270305486349461983882643116700",
    "100": "0.078534749367984244165851388048387315469262517796689",
}


def run_legacy_regression_checks() -> List[str]:
    """Verify exact reproduction of the original equal-viscosity outputs."""

    messages: List[str] = []
    tolerance = mp.power(10, -min(40, max(20, mp.mp.dps // 2)))
    for depth_text, expected_text in LEGACY_SURFACE_RESULTS.items():
        parameters = Parameters(
            depth_ratio=mp.mpf(depth_text),
            density_ratio=mp.mpf("0.8"),
            viscosity_ratio=mp.mpf("1"),
        )
        result = evaluate_mode(parameters, "surface", parameters.depth_ratio)
        expected = mp.mpf(expected_text)
        error = abs(result.critical_froude - expected)
        if error > tolerance:
            raise RuntimeError(
                "Legacy regression failed for h_r="
                f"{depth_text}: calculated {mp.nstr(result.critical_froude, 30)}, "
                f"expected {mp.nstr(expected, 30)}, error={mp.nstr(error, 8)}"
            )
        messages.append(
            f"h_r={depth_text}: Fr_l,c,S={mp.nstr(result.critical_froude, 30)}, "
            f"absolute_error={mp.nstr(error, 6)}"
        )
    return messages


def _format_number(value: mp.mpf, digits: int) -> str:
    if mp.isnan(value):
        return "nan"
    if mp.isinf(value):
        return "+inf" if value > 0 else "-inf"
    if value == 0:
        return "0"
    return mp.nstr(value, n=digits, strip_zeros=False)


def _sanitize_text(text: str) -> str:
    return " ".join(str(text).replace("\t", " ").replace("\n", " ").split())


def _automatic_output_name(sweep_parameter: str, mode: str) -> str:
    return f"kao1968_newtonian_{sweep_parameter}_{mode}_sweep.txt"


def write_results(
    output_path: Path,
    sweep_parameter: str,
    mode_selection: str,
    sweep_spacing: str,
    sweep_values: Sequence[mp.mpf],
    fixed_parameters: Parameters,
    results: Sequence[ModeResult],
    errors: Sequence[ErrorResult],
    regression_messages: Sequence[str],
    output_digits: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    columns = [
        "sweep_value",
        "h_r",
        "r_rho",
        "mu_u_over_mu_l",
        "mode",
        "celerity_c_over_U_l",
        "free_surface_over_interface_amplitude_zeta_over_eta",
        "A_in_d_equals_A_plus_B_over_Fr2",
        "B_in_d_equals_A_plus_B_over_Fr2",
        "d_at_Fr_1",
        "d_at_Fr_2",
        "critical_Fr_l",
        "first_unstable_Fr_l",
        "stability_class",
        "leading_relative_residual",
        "correction_relative_residual",
        "neutral_equation_residual",
    ]

    with output_path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Kao (1968)-type two-layer Newtonian long-wave stability sweep\n")
        stream.write(f"# generated_utc = {timestamp}\n")
        stream.write("#\n")
        stream.write("# Scope and conventions\n")
        stream.write("# ---------------------\n")
        stream.write("# Both liquid layers are Newtonian: n_l = n_u = 1.\n")
        stream.write("# h_r = H_u/H_l.\n")
        stream.write("# r_rho = rho_u/rho_l.\n")
        stream.write("# mu_u_over_mu_l = dynamic-viscosity ratio mu_u/mu_l.\n")
        stream.write("# Fr_l = U_l/sqrt(g H_l cos(theta)), with U_l the lower-layer mean velocity.\n")
        stream.write("# Passive air, zero interfacial/surface tension, and k H_l -> 0 are used.\n")
        stream.write("# Normal mode: exp(i k x + sigma t).\n")
        stream.write("# Hydrodynamic branches: sigma_j = -i c_j k + d_j(Fr_l) k^2 + O(k^3).\n")
        stream.write("# The sign convention is d_j > 0 unstable and d_j < 0 stable.\n")
        stream.write("# In this normalization, d_j(Fr_l) = A_j + B_j/Fr_l^2.\n")
        stream.write("# A_j and B_j use the internal Re_l=1 normalization; a different positive\n")
        stream.write("#   Reynolds normalization multiplies both but leaves the neutral Fr_l unchanged.\n")
        stream.write("# critical_Fr_l is the positive neutral root sqrt(-B_j/A_j), when one exists.\n")
        stream.write("# first_unstable_Fr_l is the onset as Fr_l increases from zero: it is the\n")
        stream.write("#   neutral root for stable-below/unstable-above cases, 0 when the branch is\n")
        stream.write("#   already unstable arbitrarily close to Fr_l=0, and nan when no onset occurs.\n")
        stream.write("# A missing critical root does not mean stability; read stability_class.\n")
        stream.write("# Mode labels follow Kao's long-wave ordering: the slower second mode is\n")
        stream.write("#   interfacial/varicose and the faster first mode is surface/zigzag.\n")
        stream.write("# zeta/eta is printed as a mode-shape diagnostic, not used as the label.\n")
        stream.write("#\n")
        stream.write("# Run settings\n")
        stream.write("# ------------\n")
        stream.write(f"# sweep_parameter = {sweep_parameter}\n")
        stream.write(f"# mode_selection = {mode_selection}\n")
        stream.write(f"# sweep_spacing = {sweep_spacing}\n")
        stream.write(f"# sweep_points = {len(sweep_values)}\n")
        parameter_settings = {
            "depth_ratio": fixed_parameters.depth_ratio,
            "density_ratio": fixed_parameters.density_ratio,
            "viscosity_ratio": fixed_parameters.viscosity_ratio,
        }
        for parameter_name, parameter_value in parameter_settings.items():
            if parameter_name == sweep_parameter:
                stream.write(
                    f"# {parameter_name}_setting = swept "
                    f"(baseline input {_format_number(parameter_value, output_digits)})\n"
                )
            else:
                stream.write(
                    f"# {parameter_name}_setting = "
                    f"{_format_number(parameter_value, output_digits)}\n"
                )
        stream.write(f"# mpmath_decimal_digits = {mp.mp.dps}\n")
        stream.write(f"# printed_significant_digits = {output_digits}\n")
        stream.write(f"# exact_polynomial_degree = {POLYNOMIAL_DEGREE}\n")
        stream.write("#\n")
        stream.write("# Column definitions\n")
        stream.write("# ------------------\n")
        for index, column in enumerate(columns, start=1):
            stream.write(f"# {index:02d}. {column}\n")
        stream.write("#\n")
        stream.write("\t".join(columns) + "\n")

        for result in results:
            p = result.parameters
            values = [
                _format_number(result.sweep_value, output_digits),
                _format_number(p.depth_ratio, output_digits),
                _format_number(p.density_ratio, output_digits),
                _format_number(p.viscosity_ratio, output_digits),
                result.mode,
                _format_number(result.celerity, output_digits),
                _format_number(result.surface_to_interface, output_digits),
                _format_number(result.coefficient_at_infinite_froude, output_digits),
                _format_number(
                    result.coefficient_of_inverse_froude_squared, output_digits
                ),
                _format_number(result.d_at_froude_1, output_digits),
                _format_number(result.d_at_froude_2, output_digits),
                _format_number(result.critical_froude, output_digits),
                _format_number(result.onset_froude, output_digits),
                result.stability_class,
                _format_number(result.leading_residual, output_digits),
                _format_number(result.correction_residual, output_digits),
                _format_number(result.neutral_residual, output_digits),
            ]
            stream.write("\t".join(values) + "\n")

        root_count = sum(not mp.isnan(result.critical_froude) for result in results)
        conventional_onset_count = sum(
            result.stability_class == "stable_below__unstable_above"
            for result in results
        )
        unstable_from_zero_count = sum(
            not mp.isnan(result.onset_froude) and result.onset_froude == 0
            for result in results
        )
        stable_all_count = sum(
            result.stability_class.startswith("stable_for_all")
            for result in results
        )

        stream.write("#\n")
        stream.write("# Summary\n")
        stream.write("# -------\n")
        stream.write(f"# successful_rows = {len(results)}\n")
        stream.write(f"# calculation_error_rows = {len(errors)}\n")
        stream.write(f"# rows_with_positive_neutral_root = {root_count}\n")
        stream.write(
            f"# stable_below_unstable_above_rows = {conventional_onset_count}\n"
        )
        stream.write(
            f"# rows_unstable_arbitrarily_close_to_Fr_zero = {unstable_from_zero_count}\n"
        )
        stream.write(f"# stable_for_all_positive_Fr_rows = {stable_all_count}\n")

        if regression_messages:
            stream.write("#\n")
            stream.write("# Legacy regression checks\n")
            stream.write("# ------------------------\n")
            for message in regression_messages:
                stream.write(f"# PASS: {_sanitize_text(message)}\n")

        if errors:
            stream.write("#\n")
            stream.write("# Recoverable calculation errors\n")
            stream.write("# ------------------------------\n")
            for error in errors:
                stream.write(
                    "# ERROR: sweep_value="
                    f"{_format_number(error.sweep_value, output_digits)}, "
                    f"mode={error.mode}, type={error.error_type}, "
                    f"message={_sanitize_text(error.message)}\n"
                )


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "High-precision Kao-type Newtonian k->0 critical-Froude sweep. "
            "Command-line options override the editable USER SETTINGS block."
        )
    )
    parser.add_argument(
        "--sweep",
        "--sweep-parameter",
        dest="sweep_parameter",
        default=SWEEP_PARAMETER,
        help="depth_ratio (default), density_ratio, or viscosity_ratio",
    )
    parser.add_argument(
        "--mode",
        default=MODE,
        help="surface (default), interfacial, or both",
    )
    parser.add_argument("--minimum", default=None, help="sweep minimum")
    parser.add_argument("--maximum", default=None, help="sweep maximum")
    parser.add_argument("--points", type=int, default=None, help="number of points")
    parser.add_argument(
        "--spacing",
        choices=("linear", "log"),
        default=None,
        help="sweep spacing",
    )
    parser.add_argument(
        "--values",
        default=None,
        help="comma-separated explicit sweep values; overrides min/max/points",
    )
    parser.add_argument("--depth-ratio", default=FIXED_DEPTH_RATIO)
    parser.add_argument("--density-ratio", default=FIXED_DENSITY_RATIO)
    parser.add_argument("--viscosity-ratio", default=FIXED_VISCOSITY_RATIO)
    parser.add_argument("--dps", type=int, default=MP_DECIMAL_DIGITS)
    parser.add_argument("--digits", type=int, default=OUTPUT_SIGNIFICANT_DIGITS)
    parser.add_argument(
        "--output",
        default=OUTPUT_FILE,
        help="output text file; default is generated from sweep and mode",
    )
    parser.add_argument(
        "--skip-regression-checks",
        action="store_true",
        help="skip the three legacy equal-viscosity checks",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_argument_parser()
    arguments = parser.parse_args(argv)

    if arguments.dps < 40:
        parser.error("--dps must be at least 40 for this high-precision calculation")
    if arguments.digits < 8:
        parser.error("--digits must be at least 8")
    if arguments.digits > arguments.dps - 10:
        parser.error("--digits must be at least 10 smaller than --dps")

    mp.mp.dps = int(arguments.dps)

    try:
        sweep_parameter = _normalize_sweep_name(arguments.sweep_parameter)
        mode_selection = _normalize_mode_name(arguments.mode)
        fixed_parameters = Parameters(
            depth_ratio=_as_mpf(arguments.depth_ratio, "depth_ratio"),
            density_ratio=_as_mpf(arguments.density_ratio, "density_ratio"),
            viscosity_ratio=_as_mpf(arguments.viscosity_ratio, "viscosity_ratio"),
        )
        fixed_parameters.validate()

        if arguments.values is not None:
            sweep_values = _parse_explicit_values(arguments.values)
            spacing = "custom"
        elif (
            CUSTOM_SWEEP_VALUES is not None
            and arguments.minimum is None
            and arguments.maximum is None
            and arguments.points is None
            and arguments.spacing is None
        ):
            sweep_values = [
                _as_mpf(value, "CUSTOM_SWEEP_VALUES")
                for value in CUSTOM_SWEEP_VALUES
            ]
            if any(value <= 0 for value in sweep_values):
                raise ValueError("CUSTOM_SWEEP_VALUES must all be positive")
            spacing = "custom"
        else:
            definition = SWEEP_DEFINITIONS[sweep_parameter]
            minimum = _as_mpf(
                arguments.minimum
                if arguments.minimum is not None
                else definition["minimum"],
                "minimum",
            )
            maximum = _as_mpf(
                arguments.maximum
                if arguments.maximum is not None
                else definition["maximum"],
                "maximum",
            )
            points = (
                arguments.points
                if arguments.points is not None
                else int(definition["points"])
            )
            spacing = (
                arguments.spacing
                if arguments.spacing is not None
                else str(definition["spacing"])
            )
            sweep_values = _generate_sweep_values(
                minimum, maximum, points, spacing
            )
    except ValueError as exc:
        parser.error(str(exc))

    selected_modes = (
        ("interfacial", "surface")
        if mode_selection == "both"
        else (mode_selection,)
    )

    regression_messages: List[str] = []
    if RUN_LEGACY_REGRESSION_CHECKS and not arguments.skip_regression_checks:
        regression_messages = run_legacy_regression_checks()

    results: List[ModeResult] = []
    errors: List[ErrorResult] = []

    for sweep_value in sweep_values:
        parameters = _parameters_for_sweep_value(
            sweep_parameter=sweep_parameter,
            sweep_value=sweep_value,
            fixed_depth_ratio=fixed_parameters.depth_ratio,
            fixed_density_ratio=fixed_parameters.density_ratio,
            fixed_viscosity_ratio=fixed_parameters.viscosity_ratio,
        )
        for mode_name in selected_modes:
            try:
                results.append(evaluate_mode(parameters, mode_name, sweep_value))
            except Exception as exc:  # continue the sweep and document the point
                errors.append(
                    ErrorResult(
                        sweep_value=sweep_value,
                        parameters=parameters,
                        mode=mode_name,
                        error_type=type(exc).__name__,
                        message=str(exc),
                    )
                )

    output_name = arguments.output or _automatic_output_name(
        sweep_parameter, mode_selection
    )
    output_path = Path(output_name).expanduser().resolve()
    write_results(
        output_path=output_path,
        sweep_parameter=sweep_parameter,
        mode_selection=mode_selection,
        sweep_spacing=spacing,
        sweep_values=sweep_values,
        fixed_parameters=fixed_parameters,
        results=results,
        errors=errors,
        regression_messages=regression_messages,
        output_digits=int(arguments.digits),
    )

    print("Kao-type Newtonian long-wave sweep completed.")
    print(f"  sweep parameter : {sweep_parameter}")
    print(f"  mode selection  : {mode_selection}")
    print(f"  sweep points    : {len(sweep_values)}")
    print(f"  successful rows : {len(results)}")
    print(f"  error rows      : {len(errors)}")
    print(f"  output          : {output_path}")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
