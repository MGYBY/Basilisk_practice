#!/usr/bin/env python3
"""Readable linear-stability solver for the two-layer power-law model.

The purpose of this file is educational as well as computational.  The main
calculation is written with the same scalar coefficients that appear in the
report:

    h_l,t + q_l,x = 0
    h_u,t + q_u,x = 0

    q_l,t + a31 h_l,x + a32 h_u,x + a33 q_l,x + a34 q_u,x
        = c31 h_l + c32 h_u + c33 q_l + c34 q_u

    q_u,t + a41 h_l,x + a42 h_u,x + a43 q_l,x + a44 q_u,x
        = c41 h_l + c42 h_u + c43 q_l + c44 q_u.

For an infinitely long wave,

    sigma_j(k) = -i c_j k + d_j k**2 + O(k**3),   k -> 0+.

At exactly k=0 the two hydrodynamic roots are always zero because the two layer
masses are conserved.  The long-wave stability threshold is therefore d_j=0,
not sigma_j(0)=0.

For fixed depth ratio, density ratio, power-law indices and scaled consistency
ratio, the normalized normal-flow family satisfies

    d_j(Fr_l) = a_j Fr_l**2 + b_j.

A positive critical Froude number exists when a_j*b_j < 0:

    Fr_l,c,j = sqrt(-b_j/a_j).

The file also evaluates the complete finite-wavenumber dispersion relation by
solving the four eigenvalues of C0 - i*k*A0.  NumPy is the only requirement.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, replace
from typing import Callable, Sequence

import numpy as np


# -----------------------------------------------------------------------------
# 1. User parameters and result containers
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Parameters:
    """Dimensionless inputs of the normalized two-layer model."""

    froude: float = 0.8
    depth_ratio: float = 1.0
    density_ratio: float = 0.9
    n_lower: float = 0.6
    n_upper: float = 0.8
    scaled_consistency_ratio: float = 1.0

    def validated(self) -> "Parameters":
        values = {
            "froude": self.froude,
            "depth_ratio": self.depth_ratio,
            "density_ratio": self.density_ratio,
            "n_lower": self.n_lower,
            "n_upper": self.n_upper,
            "scaled_consistency_ratio": self.scaled_consistency_ratio,
        }
        for name, value in values.items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive; received {value!r}")
        return self


@dataclass(frozen=True)
class LinearCoefficients:
    """The sixteen scalar coefficients of the uniform-flow linearization."""

    # Base state and closure diagnostics.
    h_lower_0: float
    h_upper_0: float
    q_lower_0: float
    q_upper_0: float
    lambda_0: float
    lambda_derivatives: tuple[float, float, float, float]
    residual_lambda: float

    # Spatial-gradient coefficients in the two momentum equations.
    a31: float
    a32: float
    a33: float
    a34: float
    a41: float
    a42: float
    a43: float
    a44: float

    # Local source/relaxation coefficients.
    c31: float
    c32: float
    c33: float
    c34: float
    c41: float
    c42: float
    c43: float
    c44: float

    def A(self) -> np.ndarray:
        """Return A0 for finite-wavenumber calculations.

        The report develops the scalar equations first.  This compact array is
        constructed only when NumPy needs the full eigenproblem.
        """
        return np.array(
            [
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
                [self.a31, self.a32, self.a33, self.a34],
                [self.a41, self.a42, self.a43, self.a44],
            ],
            dtype=float,
        )

    def C(self) -> np.ndarray:
        """Return C0 for finite-wavenumber calculations."""
        return np.array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [self.c31, self.c32, self.c33, self.c34],
                [self.c41, self.c42, self.c43, self.c44],
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class LongWaveMode:
    """One hydrodynamic mode in the k -> 0+ expansion."""

    label: str
    celerity: float
    h_lower_amplitude: float
    h_upper_amplitude: float
    d_at_froude: float

    @property
    def free_surface_amplitude(self) -> float:
        return self.h_lower_amplitude + self.h_upper_amplitude


@dataclass(frozen=True)
class CriticalMode:
    """Affine long-wave coefficient and its positive critical Froude number."""

    label: str
    celerity: float
    h_lower_amplitude: float
    h_upper_amplitude: float
    slope_a: float
    intercept_b: float
    critical_froude: float
    stability_description: str

    def coefficient(self, froude: float) -> float:
        return self.slope_a * froude**2 + self.intercept_b


DEFAULT_PARAMETERS = Parameters()


# -----------------------------------------------------------------------------
# 2. Karman-Pohlhausen shape functions and exact normal flow
# -----------------------------------------------------------------------------


def lower_shape_functions(lam: complex, n_lower: float) -> tuple[complex, complex, complex]:
    """Return C_l, A_l and B_l for the lower-layer velocity profile.

    C_l relates interfacial velocity to basal traction, A_l relates mean
    velocity to interfacial velocity, and B_l gives the momentum moment.
    A short series is used close to lam=1 to avoid cancellation and remains
    compatible with complex-step differentiation.
    """
    m = (n_lower + 1.0) / n_lower
    epsilon = 1.0 - lam

    if abs(epsilon) < 1.0e-7:
        c_value = (
            1.0
            - 0.5 * (m - 1.0) * epsilon
            + (m - 1.0) * (m - 2.0) * epsilon**2 / 6.0
        )
        a_value = 0.5 + (m - 1.0) * epsilon / 12.0
        b_value = 1.0 / 3.0 + (m - 1.0) * epsilon / 12.0
        return c_value, a_value, b_value

    lambda_to_m = lam**m
    denominator = 1.0 - lambda_to_m
    j_m = (1.0 - lam ** (m + 1.0)) / ((m + 1.0) * (1.0 - lam))
    j_2m = (1.0 - lam ** (2.0 * m + 1.0)) / ((2.0 * m + 1.0) * (1.0 - lam))

    c_value = denominator / (m * (1.0 - lam))
    a_value = (1.0 - j_m) / denominator
    b_value = (1.0 - 2.0 * j_m + j_2m) / denominator**2
    return c_value, a_value, b_value


def upper_shape_constants(n_upper: float) -> tuple[float, float, float]:
    """Return m_u, A_u and B_u for the upper profile."""
    m_upper = (n_upper + 1.0) / n_upper
    a_upper = m_upper / (m_upper + 1.0)
    b_upper = 1.0 - 2.0 / (m_upper + 1.0) + 1.0 / (2.0 * m_upper + 1.0)
    return m_upper, a_upper, b_upper


def normal_flow_quantities(parameters: Parameters) -> dict[str, float]:
    """Return all quantities needed to define the normalized uniform state."""
    p = parameters.validated()

    lambda_0 = p.density_ratio * p.depth_ratio / (
        1.0 + p.density_ratio * p.depth_ratio
    )
    c_lower, a_lower, b_lower = lower_shape_functions(
        complex(lambda_0), p.n_lower
    )
    m_upper, a_upper, b_upper = upper_shape_constants(p.n_upper)

    # The normalized family is chosen so that h_l0=q_l0=1.  The rheological
    # coefficients therefore scale as Fr_l^{-2}.
    lambda_lower = (
        (1.0 + p.density_ratio * p.depth_ratio)
        * (a_lower.real * c_lower.real) ** p.n_lower
        / p.froude**2
    )
    lambda_upper = p.scaled_consistency_ratio * lambda_lower

    tau_b_0 = (1.0 + p.density_ratio * p.depth_ratio) / p.froude**2
    tau_i_0 = p.density_ratio * p.depth_ratio / p.froude**2
    interfacial_velocity_0 = 1.0 / a_lower.real
    upper_velocity_increment_0 = (
        p.depth_ratio
        / m_upper
        * (tau_i_0 / lambda_upper) ** (1.0 / p.n_upper)
    )
    q_upper_0 = p.depth_ratio * (
        interfacial_velocity_0 + a_upper * upper_velocity_increment_0
    )

    return {
        "lambda_0": lambda_0,
        "C_lower": c_lower.real,
        "A_lower": a_lower.real,
        "B_lower": b_lower.real,
        "m_upper": m_upper,
        "A_upper": a_upper,
        "B_upper": b_upper,
        "Lambda_lower": lambda_lower,
        "Lambda_upper": lambda_upper,
        "tau_b_0": tau_b_0,
        "tau_i_0": tau_i_0,
        "U_i_0": interfacial_velocity_0,
        "W_0": upper_velocity_increment_0,
        "q_upper_0": q_upper_0,
    }


def uniform_state(parameters: Parameters) -> np.ndarray:
    """Return [h_l0, h_u0, q_l0, q_u0]."""
    reference = normal_flow_quantities(parameters)
    return np.array(
        [1.0, parameters.depth_ratio, 1.0, reference["q_upper_0"]],
        dtype=float,
    )


# -----------------------------------------------------------------------------
# 3. Closure and total derivatives
# -----------------------------------------------------------------------------


def closure_with_lambda_held_fixed(
    state: np.ndarray,
    lam: complex,
    parameters: Parameters,
) -> dict[str, complex]:
    """Evaluate the local closure while treating lambda as independent.

    The implicit traction-continuity correction is applied later.  Separating
    these two steps makes the chain rule visible and avoids finite-differencing
    a nonlinear root solver.
    """
    h_lower, h_upper, q_lower, q_upper = state
    p = parameters

    c_lower, a_lower, b_lower = lower_shape_functions(lam, p.n_lower)
    m_upper, a_upper, b_upper = upper_shape_constants(p.n_upper)
    reference = normal_flow_quantities(p)

    mean_lower = q_lower / h_lower
    mean_upper = q_upper / h_upper
    interfacial_velocity = mean_lower / a_lower
    upper_increment = (mean_upper - interfacial_velocity) / a_upper

    tau_b = reference["Lambda_lower"] * (
        interfacial_velocity / (h_lower * c_lower)
    ) ** p.n_lower
    tau_i_from_upper = reference["Lambda_upper"] * (
        m_upper * upper_increment / h_upper
    ) ** p.n_upper

    residual = lam * tau_b - tau_i_from_upper
    tau_i = lam * tau_b
    momentum_lower = h_lower * interfacial_velocity**2 * b_lower
    momentum_upper = h_upper * (
        interfacial_velocity**2
        + 2.0 * a_upper * interfacial_velocity * upper_increment
        + b_upper * upper_increment**2
    )

    return {
        "residual": residual,
        "tau_b": tau_b,
        "tau_i": tau_i,
        "tau_i_from_upper": tau_i_from_upper,
        "momentum_lower": momentum_lower,
        "momentum_upper": momentum_upper,
        "interfacial_velocity": interfacial_velocity,
        "upper_increment": upper_increment,
    }


def flux_with_lambda_held_fixed(
    state: np.ndarray,
    lam: complex,
    parameters: Parameters,
) -> np.ndarray:
    h_lower, h_upper, q_lower, q_upper = state
    closure = closure_with_lambda_held_fixed(state, lam, parameters)
    return np.array(
        [
            q_lower,
            q_upper,
            closure["momentum_lower"] + h_lower**2 / (2.0 * parameters.froude**2),
            closure["momentum_upper"] + h_upper**2 / (2.0 * parameters.froude**2),
        ],
        dtype=complex,
    )


def source_with_lambda_held_fixed(
    state: np.ndarray,
    lam: complex,
    parameters: Parameters,
) -> np.ndarray:
    h_lower, h_upper = state[:2]
    closure = closure_with_lambda_held_fixed(state, lam, parameters)
    return np.array(
        [
            0.0,
            0.0,
            h_lower / parameters.froude**2
            + closure["tau_i"]
            - closure["tau_b"],
            h_upper / parameters.froude**2
            - closure["tau_i"] / parameters.density_ratio,
        ],
        dtype=complex,
    )


def complex_step_jacobian(
    function: Callable[[np.ndarray], np.ndarray],
    point: np.ndarray,
    step: float = 1.0e-30,
) -> np.ndarray:
    """Differentiate a real-valued vector function without subtraction error."""
    point = np.asarray(point, dtype=float)
    sample = np.asarray(function(point.astype(complex))).reshape(-1)
    jacobian = np.empty((sample.size, point.size), dtype=float)

    for column in range(point.size):
        perturbed = point.astype(complex)
        perturbed[column] += 1j * step
        jacobian[:, column] = np.imag(function(perturbed)).reshape(-1) / step

    return jacobian


def linear_coefficients(parameters: Parameters) -> LinearCoefficients:
    """Construct every scalar coefficient in the four linearized PDEs."""
    p = parameters.validated()
    state_0 = uniform_state(p)
    reference = normal_flow_quantities(p)
    lambda_0 = reference["lambda_0"]
    step = 1.0e-30

    # Partial derivatives of R(lambda,Q) with lambda held fixed.
    def residual_as_vector(state: np.ndarray) -> np.ndarray:
        value = closure_with_lambda_held_fixed(state, lambda_0, p)["residual"]
        return np.array([value], dtype=complex)

    residual_Q = complex_step_jacobian(residual_as_vector, state_0, step).reshape(4)
    residual_lambda = np.imag(
        closure_with_lambda_held_fixed(
            state_0.astype(complex), lambda_0 + 1j * step, p
        )["residual"]
    ) / step
    if abs(residual_lambda) < 1.0e-12:
        raise ValueError(
            "R_lambda is too small; the implicit traction-ratio closure is not regular."
        )
    lambda_Q = -residual_Q / residual_lambda

    # Fixed-lambda partial derivatives of flux and source.
    flux_Q_fixed = complex_step_jacobian(
        lambda state: flux_with_lambda_held_fixed(state, lambda_0, p),
        state_0,
        step,
    )
    source_Q_fixed = complex_step_jacobian(
        lambda state: source_with_lambda_held_fixed(state, lambda_0, p),
        state_0,
        step,
    )
    flux_lambda = np.imag(
        flux_with_lambda_held_fixed(state_0, lambda_0 + 1j * step, p)
    ) / step
    source_lambda = np.imag(
        source_with_lambda_held_fixed(state_0, lambda_0 + 1j * step, p)
    ) / step

    # Total derivatives: dG/dQ = G_Q + G_lambda * lambda_Q.
    flux_Q_total = flux_Q_fixed + np.outer(flux_lambda, lambda_Q)
    source_Q_total = source_Q_fixed + np.outer(source_lambda, lambda_Q)

    h_lower_0, h_upper_0, q_lower_0, q_upper_0 = state_0

    # Add the two cross-layer hydrostatic products to the momentum-gradient
    # coefficients.  This is the scalar equivalent of F_Q + B0.
    a31 = flux_Q_total[2, 0]
    a32 = flux_Q_total[2, 1] + p.density_ratio * h_lower_0 / p.froude**2
    a33 = flux_Q_total[2, 2]
    a34 = flux_Q_total[2, 3]

    a41 = flux_Q_total[3, 0] + h_upper_0 / p.froude**2
    a42 = flux_Q_total[3, 1]
    a43 = flux_Q_total[3, 2]
    a44 = flux_Q_total[3, 3]

    return LinearCoefficients(
        h_lower_0=float(h_lower_0),
        h_upper_0=float(h_upper_0),
        q_lower_0=float(q_lower_0),
        q_upper_0=float(q_upper_0),
        lambda_0=float(lambda_0),
        lambda_derivatives=tuple(float(value) for value in lambda_Q),
        residual_lambda=float(residual_lambda),
        a31=float(a31),
        a32=float(a32),
        a33=float(a33),
        a34=float(a34),
        a41=float(a41),
        a42=float(a42),
        a43=float(a43),
        a44=float(a44),
        c31=float(source_Q_total[2, 0]),
        c32=float(source_Q_total[2, 1]),
        c33=float(source_Q_total[2, 2]),
        c34=float(source_Q_total[2, 3]),
        c41=float(source_Q_total[3, 0]),
        c42=float(source_Q_total[3, 1]),
        c43=float(source_Q_total[3, 2]),
        c44=float(source_Q_total[3, 3]),
    )


# -----------------------------------------------------------------------------
# 4. Scalar long-wave reduction
# -----------------------------------------------------------------------------


def equilibrium_discharge_coefficients(
    coefficients: LinearCoefficients,
) -> tuple[float, float, float, float, float]:
    """Return Delta_c and e11,e12,e21,e22.

    At leading long-wave order the two momentum sources are in local
    equilibrium.  Solving those two scalar equations gives

        q_l_hat = e11 h_l_hat + e12 h_u_hat,
        q_u_hat = e21 h_l_hat + e22 h_u_hat.
    """
    c = coefficients
    delta_c = c.c33 * c.c44 - c.c34 * c.c43
    if abs(delta_c) < 1.0e-12:
        raise ValueError(
            "Delta_c = c33*c44-c34*c43 is too small; momentum relaxation is singular."
        )

    e11 = (c.c34 * c.c41 - c.c44 * c.c31) / delta_c
    e12 = (c.c34 * c.c42 - c.c44 * c.c32) / delta_c
    e21 = (c.c43 * c.c31 - c.c33 * c.c41) / delta_c
    e22 = (c.c43 * c.c32 - c.c33 * c.c42) / delta_c
    return delta_c, e11, e12, e21, e22


def longwave_celerities(
    coefficients: LinearCoefficients,
) -> tuple[float, float]:
    """Solve the explicit quadratic for the two equilibrium celerities."""
    _, e11, e12, e21, e22 = equilibrium_discharge_coefficients(coefficients)
    discriminant = (e11 - e22) ** 2 + 4.0 * e12 * e21
    tolerance = 1.0e-12 * max(1.0, abs(e11), abs(e12), abs(e21), abs(e22)) ** 2
    if discriminant < -tolerance:
        raise ValueError(
            "The two equilibrium celerities are complex; the long-wave equilibrium "
            "system is not hyperbolic for this parameter set."
        )
    discriminant = max(discriminant, 0.0)
    root = np.sqrt(discriminant)
    c_minus = 0.5 * (e11 + e22 - root)
    c_plus = 0.5 * (e11 + e22 + root)
    return float(c_minus), float(c_plus)


def _right_depth_amplitudes(
    celerity: float,
    e11: float,
    e12: float,
    e21: float,
    e22: float,
) -> tuple[float, float]:
    """Return a nonzero right eigenvector using only scalar operations."""
    candidate_1 = np.array([e12, celerity - e11], dtype=float)
    candidate_2 = np.array([celerity - e22, e21], dtype=float)
    candidate = candidate_1 if np.linalg.norm(candidate_1) >= np.linalg.norm(candidate_2) else candidate_2
    scale = np.max(np.abs(candidate))
    if scale < 1.0e-14:
        raise ValueError("Could not construct a nonzero long-wave depth eigenvector.")
    candidate /= scale
    # Use a repeatable sign: make the upper-depth component nonnegative.
    if candidate[1] < 0.0:
        candidate *= -1.0
    return float(candidate[0]), float(candidate[1])


def _left_depth_amplitudes(
    celerity: float,
    e11: float,
    e12: float,
    e21: float,
    e22: float,
) -> tuple[float, float]:
    """Return a left eigenvector, again without matrix notation."""
    candidate_1 = np.array([e21, celerity - e11], dtype=float)
    candidate_2 = np.array([celerity - e22, e12], dtype=float)
    candidate = candidate_1 if np.linalg.norm(candidate_1) >= np.linalg.norm(candidate_2) else candidate_2
    if np.linalg.norm(candidate) < 1.0e-14:
        raise ValueError("Could not construct a nonzero left long-wave eigenvector.")
    return float(candidate[0]), float(candidate[1])


def longwave_modes_at_froude(parameters: Parameters) -> list[LongWaveMode]:
    """Calculate c_j and d_j directly from the scalar coefficient equations."""
    coefficients = linear_coefficients(parameters)
    delta_c, e11, e12, e21, e22 = equilibrium_discharge_coefficients(coefficients)
    celerities = longwave_celerities(coefficients)

    raw_modes: list[LongWaveMode] = []
    for index, celerity in enumerate(celerities):
        r_lower, r_upper = _right_depth_amplitudes(
            celerity, e11, e12, e21, e22
        )
        l_lower, l_upper = _left_depth_amplitudes(
            celerity, e11, e12, e21, e22
        )

        # Leading equilibrium discharge amplitudes are c*r because the mass
        # equations require q_hat = c*h_hat at O(k).
        s_lower = celerity * r_lower
        s_upper = celerity * r_upper

        # The two scalar momentum equations supply the first departure from
        # local equilibrium.  These are the f_l and f_u in the report.
        f_lower = (
            coefficients.a31 * r_lower
            + coefficients.a32 * r_upper
            + coefficients.a33 * s_lower
            + coefficients.a34 * s_upper
            - celerity * s_lower
        )
        f_upper = (
            coefficients.a41 * r_lower
            + coefficients.a42 * r_upper
            + coefficients.a43 * s_lower
            + coefficients.a44 * s_upper
            - celerity * s_upper
        )

        # Solve the two relaxation equations for z_l and z_u explicitly.
        z_lower = (
            coefficients.c44 * f_lower - coefficients.c34 * f_upper
        ) / delta_c
        z_upper = (
            -coefficients.c43 * f_lower + coefficients.c33 * f_upper
        ) / delta_c

        denominator = l_lower * r_lower + l_upper * r_upper
        if abs(denominator) < 1.0e-13:
            raise ValueError("Left and right long-wave modes are orthogonal numerically.")
        d_value = (l_lower * z_lower + l_upper * z_upper) / denominator

        raw_modes.append(
            LongWaveMode(
                label=f"mode_{index + 1}",
                celerity=float(celerity),
                h_lower_amplitude=float(r_lower),
                h_upper_amplitude=float(r_upper),
                d_at_froude=float(d_value),
            )
        )

    # Identify the surface-dominated mode by the relative motion of the total
    # depth h_l+h_u.  This is more physical than labeling only by speed.
    surface_measure = [
        abs(mode.free_surface_amplitude)
        / (abs(mode.h_lower_amplitude) + abs(mode.h_upper_amplitude))
        for mode in raw_modes
    ]
    surface_index = int(np.argmax(surface_measure))

    labelled: list[LongWaveMode] = []
    for index, mode in enumerate(raw_modes):
        label = "free_surface" if index == surface_index else "interfacial"
        labelled.append(replace(mode, label=label))
    return labelled


def _stability_description(slope: float, intercept: float, critical: float) -> str:
    tolerance = 1.0e-12 * max(1.0, abs(slope), abs(intercept))
    if np.isfinite(critical):
        if slope > 0.0:
            return "stable below the root; unstable above it"
        return "unstable below the root; stable above it"
    if slope > tolerance and intercept >= -tolerance:
        return "unstable for every positive Froude number"
    if slope < -tolerance and intercept <= tolerance:
        return "stable for every positive Froude number"
    if abs(slope) <= tolerance:
        if intercept > tolerance:
            return "Froude-independent and unstable"
        if intercept < -tolerance:
            return "Froude-independent and stable"
        return "degenerate neutral coefficient"
    return "no positive root; inspect the sign of d(Fr)"


def critical_longwave_modes(parameters: Parameters) -> list[CriticalMode]:
    """Return d_j=a_j*Fr^2+b_j and the corresponding positive roots.

    Two Froude values are sufficient because the Fr^2 dependence is exact for
    this normalized normal-flow family.  Fr=1 and Fr=2 are used because they
    are well separated and keep the arithmetic transparent:

        d(1) = a+b,
        d(2) = 4a+b,
        a = [d(2)-d(1)]/3,
        b = d(1)-a.
    """
    base = parameters.validated()
    modes_at_1 = sorted(
        longwave_modes_at_froude(replace(base, froude=1.0)),
        key=lambda mode: mode.celerity,
    )
    modes_at_2 = sorted(
        longwave_modes_at_froude(replace(base, froude=2.0)),
        key=lambda mode: mode.celerity,
    )

    results: list[CriticalMode] = []
    for mode_1, mode_2 in zip(modes_at_1, modes_at_2):
        if abs(mode_1.celerity - mode_2.celerity) > 1.0e-9:
            raise AssertionError(
                "Long-wave celerity changed when only Fr was varied; the exact "
                "normal-family scaling is not satisfied."
            )
        slope = (mode_2.d_at_froude - mode_1.d_at_froude) / 3.0
        intercept = mode_1.d_at_froude - slope
        critical = np.nan
        if slope * intercept < 0.0:
            critical = float(np.sqrt(-intercept / slope))

        # The physical label is identical at Fr=1 and Fr=2 for the regular
        # branch; use the Fr=1 mode for amplitudes and celerity.
        results.append(
            CriticalMode(
                label=mode_1.label,
                celerity=mode_1.celerity,
                h_lower_amplitude=mode_1.h_lower_amplitude,
                h_upper_amplitude=mode_1.h_upper_amplitude,
                slope_a=float(slope),
                intercept_b=float(intercept),
                critical_froude=critical,
                stability_description=_stability_description(slope, intercept, critical),
            )
        )

    # Return interfacial first and free-surface second for predictable output.
    return sorted(results, key=lambda mode: 0 if mode.label == "interfacial" else 1)


# -----------------------------------------------------------------------------
# 5. Complete finite-wavenumber dispersion relation
# -----------------------------------------------------------------------------


def dispersion_roots(wavenumber: float, parameters: Parameters) -> np.ndarray:
    """Return the four temporal roots sigma(k) of the complete linear system."""
    coefficients = linear_coefficients(parameters)
    return np.linalg.eigvals(coefficients.C() - 1j * wavenumber * coefficients.A())


def _best_permutation(reference: np.ndarray, candidates: np.ndarray) -> np.ndarray:
    """Match four roots by the smallest total complex-plane distance."""
    best_cost = np.inf
    best: tuple[int, ...] | None = None
    for permutation in itertools.permutations(range(len(candidates))):
        ordered = candidates[list(permutation)]
        cost = float(np.sum(np.abs(ordered - reference)))
        if cost < best_cost:
            best_cost = cost
            best = permutation
    assert best is not None
    return candidates[list(best)]


def tracked_dispersion(
    wavenumbers: Sequence[float],
    parameters: Parameters,
) -> dict[str, np.ndarray]:
    """Track interface, free-surface and two relaxation roots over k.

    The matching is intentionally simple: roots are paired by continuity in
    the complex plane.  This is easier to read than a general assignment
    package and is reliable for the smooth grids used in the supplied examples.
    """
    k_values = np.asarray(wavenumbers, dtype=float)
    if k_values.ndim != 1 or k_values.size == 0 or np.any(k_values <= 0.0):
        raise ValueError("wavenumbers must be a nonempty one-dimensional positive array")
    if np.any(np.diff(k_values) <= 0.0):
        raise ValueError("wavenumbers must be strictly increasing")

    long_modes = critical_longwave_modes(parameters)
    at_reference = {
        mode.label: mode for mode in longwave_modes_at_froude(parameters)
    }
    coefficients = linear_coefficients(parameters)
    relaxation = np.linalg.eigvals(
        np.array(
            [[coefficients.c33, coefficients.c34], [coefficients.c43, coefficients.c44]],
            dtype=float,
        )
    )
    relaxation = relaxation[np.argsort(relaxation.real)]

    first_k = k_values[0]
    expected = np.array(
        [
            -1j * at_reference["interfacial"].celerity * first_k
            + at_reference["interfacial"].d_at_froude * first_k**2,
            -1j * at_reference["free_surface"].celerity * first_k
            + at_reference["free_surface"].d_at_froude * first_k**2,
            relaxation[0],
            relaxation[1],
        ],
        dtype=complex,
    )

    branches = np.empty((k_values.size, 4), dtype=complex)
    branches[0] = _best_permutation(expected, dispersion_roots(first_k, parameters))
    for index, k_value in enumerate(k_values[1:], start=1):
        roots = dispersion_roots(float(k_value), parameters)
        branches[index] = _best_permutation(branches[index - 1], roots)

    return {
        "wavenumber": k_values,
        "interfacial": branches[:, 0],
        "free_surface": branches[:, 1],
        "relaxation_fast": branches[:, 2],
        "relaxation_slow": branches[:, 3],
    }


