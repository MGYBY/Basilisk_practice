#!/usr/bin/env python3
"""Readable full-Navier--Stokes stability solver for two power-law layers.

The mathematical problem is a two-domain generalized Orr--Sommerfeld
problem.  The lower liquid occupies 0 <= z <= 1 and the upper liquid
occupies 1 <= z <= 1 + h_r after nondimensionalization with the steady
lower-layer depth H_l.  The velocity scale is the steady lower-layer
mean velocity U_l.

This module contains only the reusable mathematics:

1. construct the dimensional steady parallel base flow;
2. nondimensionalize it;
3. build two Chebyshev collocation grids;
4. assemble the generalized matrix pencil A q = c B q, where c is the
   complex wave speed;
5. solve and track the two hydrodynamic branches.

The phase-speed formulation is deliberately used instead of solving
for sigma directly.  With the normal mode exp(i k x + sigma t),

    sigma = -i k c,

so Re(c) is the phase celerity and k Im(c) is the temporal growth rate.
This keeps the hydrodynamic eigenvalues finite as k -> 0 and makes the
kinematic boundary conditions particularly transparent.

Surface and interfacial tension are omitted.  A smooth regularized
power-law relation is used because an ideal shear-thinning upper layer
has infinite apparent viscosity at its stress-free surface.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, Sequence
import math
import warnings

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.linalg import eig
from scipy.optimize import brentq


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


# ---------------------------------------------------------------------------
# 1. User-facing physical and numerical data
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FluidLayer:
    """Dimensional properties of one power-law liquid layer."""

    depth: float
    density: float
    consistency: float
    flow_index: float

    def validate(self, name: str) -> None:
        values = {
            "depth": self.depth,
            "density": self.density,
            "consistency": self.consistency,
            "flow_index": self.flow_index,
        }
        for field_name, value in values.items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(
                    f"{name}.{field_name} must be positive and finite; got {value!r}"
                )
        if self.flow_index > 1.0:
            warnings.warn(
                f"{name}.flow_index={self.flow_index:g} describes shear thickening, "
                "not shear thinning.  The equations still remain valid.",
                RuntimeWarning,
            )


@dataclass(frozen=True)
class CaseParameters:
    """Dimensional definition of the two-layer normal-flow case."""

    gravity: float = 9.81
    slope_tan: float = 0.05
    lower: FluidLayer = FluidLayer(
        depth=0.025,
        density=1500.0,
        consistency=6.0,
        flow_index=0.40,
    )
    upper: FluidLayer = FluidLayer(
        depth=0.025,
        density=1200.0,
        consistency=3.0,
        flow_index=0.80,
    )
    # Dimensional smoothing shear rate [s^-1].
    regularization_rate: float = 0.5

    def validate(self) -> None:
        if not np.isfinite(self.gravity) or self.gravity <= 0.0:
            raise ValueError("gravity must be positive and finite")
        if not np.isfinite(self.slope_tan) or self.slope_tan <= 0.0:
            raise ValueError("slope_tan must be positive and finite")
        if not np.isfinite(self.regularization_rate) or self.regularization_rate < 0.0:
            raise ValueError("regularization_rate must be finite and nonnegative")
        self.lower.validate("lower")
        self.upper.validate("upper")
        if self.regularization_rate == 0.0 and self.upper.flow_index < 1.0:
            warnings.warn(
                "An ideal shear-thinning upper power-law layer is singular at the "
                "stress-free free surface.  Use a positive regularization_rate for "
                "the numerical eigenproblem.",
                RuntimeWarning,
            )


@dataclass(frozen=True)
class SolverOptions:
    """Numerical controls for the spectral generalized eigenproblem."""

    chebyshev_order: int = 64
    residual_tolerance: float = 2.0e-6
    phase_speed_limit: float = 100.0
    seed_wavenumber: float = 5.0e-4
    base_samples: int = 6001

    def validate(self) -> None:
        if self.chebyshev_order < 12:
            raise ValueError("chebyshev_order should be at least 12")
        if self.residual_tolerance <= 0.0:
            raise ValueError("residual_tolerance must be positive")
        if self.phase_speed_limit <= 0.0:
            raise ValueError("phase_speed_limit must be positive")
        if self.seed_wavenumber <= 0.0:
            raise ValueError("seed_wavenumber must be positive")
        if self.base_samples < 200:
            raise ValueError("base_samples should be at least 200")


# ---------------------------------------------------------------------------
# 2. Base-state data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LayerProfile:
    """Dense nondimensional base profile in one liquid layer."""

    z: FloatArray
    velocity: FloatArray
    shear: FloatArray
    pressure: FloatArray


@dataclass(frozen=True)
class BaseFlow:
    """Complete nondimensional steady parallel base flow."""

    parameters: CaseParameters
    velocity_scale: float
    normal_gravity: float
    downslope_gravity: float
    froude: float
    depth_ratio: float
    density_ratio: float
    consistency_number_lower: float
    consistency_number_upper: float
    epsilon_lower: float
    epsilon_upper: float
    interface_velocity: float
    surface_velocity: float
    lower: LayerProfile
    upper: LayerProfile

    def summary(self) -> dict[str, float]:
        return {
            "velocity_scale_m_per_s": self.velocity_scale,
            "normal_gravity_m_per_s2": self.normal_gravity,
            "downslope_gravity_m_per_s2": self.downslope_gravity,
            "froude_lower": self.froude,
            "slope_tan": self.parameters.slope_tan,
            "depth_ratio": self.depth_ratio,
            "density_ratio": self.density_ratio,
            "Lambda_lower": self.consistency_number_lower,
            "Lambda_upper": self.consistency_number_upper,
            "epsilon_lower": self.epsilon_lower,
            "epsilon_upper": self.epsilon_upper,
            "interface_velocity_dimensionless": self.interface_velocity,
            "surface_velocity_dimensionless": self.surface_velocity,
        }


# ---------------------------------------------------------------------------
# 3. Smooth regularized power-law rheology and steady base flow
# ---------------------------------------------------------------------------


def regularized_stress(
    shear_rate: FloatArray | float,
    consistency: float,
    flow_index: float,
    regularization_rate: float,
) -> FloatArray | float:
    """Return tau(gamma) for the smooth regularized dimensional law.

    tau = K (gamma^2 + gamma_r^2)^((n-1)/2) gamma.
    """

    gamma = np.asarray(shear_rate)
    return (
        consistency
        * np.power(gamma * gamma + regularization_rate * regularization_rate,
                   0.5 * (flow_index - 1.0))
        * gamma
    )


def invert_regularized_stress(
    stress: FloatArray | float,
    layer: FluidLayer,
    regularization_rate: float,
    *,
    tolerance: float = 2.0e-13,
    maximum_iterations: int = 80,
) -> FloatArray:
    """Invert the monotone regularized power-law stress relation.

    The positive-shear constitutive equation is

        tau/K = gamma (gamma**2 + gamma_r**2)**((n - 1)/2).

    Earlier versions solved one scalar Brent problem at every profile point.
    This implementation applies a vectorized safeguarded Newton iteration.
    A lower and upper bracket is retained at every point, so an unacceptable
    Newton step is replaced by bisection.  The result is both faster and easier
    to use in parameter sweeps while preserving the monotone constitutive root.
    """

    tau = np.asarray(stress, dtype=float)
    if np.any(~np.isfinite(tau)):
        raise ValueError("stress must be finite")
    if np.any(tau < -1.0e-13):
        raise ValueError("The selected base-flow branch requires nonnegative stress")
    if regularization_rate < 0.0 or not np.isfinite(regularization_rate):
        raise ValueError("regularization_rate must be finite and nonnegative")

    reduced = np.maximum(tau, 0.0) / layer.consistency
    n = layer.flow_index
    epsilon = regularization_rate

    if epsilon == 0.0:
        return np.power(reduced, 1.0 / n)

    # Small-shear and ideal-power-law asymptotes provide useful initial scales.
    small_shear = reduced * epsilon ** (1.0 - n)
    ideal_shear = np.power(reduced, 1.0 / n)

    lower = np.zeros_like(reduced)
    upper = np.maximum.reduce(
        [
            2.0 * small_shear + epsilon,
            2.0 * ideal_shear + epsilon,
            np.full_like(reduced, max(epsilon, 1.0e-14)),
        ]
    )

    def reduced_stress(gamma: FloatArray) -> FloatArray:
        return gamma * np.power(
            gamma * gamma + epsilon * epsilon,
            0.5 * (n - 1.0),
        )

    # Expand any bracket that is still too small.
    for _ in range(40):
        mask = reduced_stress(upper) < reduced
        if not np.any(mask):
            break
        upper[mask] *= 2.0
    else:
        raise RuntimeError("Could not bracket the regularized shear-rate root")

    gamma = np.clip(0.5 * (small_shear + ideal_shear), lower, upper)
    gamma = np.where(reduced == 0.0, 0.0, gamma)

    for _ in range(maximum_iterations):
        value = reduced_stress(gamma)
        above = value > reduced
        upper = np.where(above, gamma, upper)
        lower = np.where(~above, gamma, lower)

        derivative = np.power(
            gamma * gamma + epsilon * epsilon,
            0.5 * (n - 3.0),
        ) * (epsilon * epsilon + n * gamma * gamma)
        newton = gamma - (value - reduced) / derivative
        midpoint = 0.5 * (lower + upper)
        acceptable = (newton > lower) & (newton < upper) & np.isfinite(newton)
        gamma = np.where(acceptable, newton, midpoint)

        width = np.max(upper - lower)
        scale = 1.0 + np.max(np.abs(gamma))
        if width <= tolerance * scale:
            break

    gamma = np.where(reduced == 0.0, 0.0, gamma)
    return gamma


def apparent_and_tangent_viscosity(
    shear: FloatArray,
    consistency_coefficient: float,
    flow_index: float,
    epsilon: float,
) -> tuple[FloatArray, FloatArray]:
    """Return dimensionless apparent and tangent viscosities.

    mu     = Lambda (gamma^2 + epsilon^2)^((n-1)/2),
    mu_t   = d(mu gamma)/d gamma.
    """

    gamma = np.asarray(shear, dtype=float)
    squared = gamma * gamma + epsilon * epsilon
    mu = consistency_coefficient * np.power(
        squared,
        0.5 * (flow_index - 1.0),
    )

    if epsilon == 0.0:
        if flow_index < 1.0 and np.any(gamma == 0.0):
            mu = np.where(gamma == 0.0, np.inf, mu)
        mu_t = flow_index * mu
    else:
        mu_t = consistency_coefficient * np.power(
            squared,
            0.5 * (flow_index - 3.0),
        ) * (epsilon * epsilon + flow_index * gamma * gamma)

    return mu, mu_t


def build_base_flow(
    parameters: CaseParameters,
    samples: int = 6001,
) -> BaseFlow:
    """Construct the mechanically consistent regularized normal flow."""

    parameters.validate()
    if samples < 200:
        raise ValueError("samples should be at least 200")

    lower = parameters.lower
    upper = parameters.upper
    angle = math.atan(parameters.slope_tan)
    g_s = parameters.gravity * math.sin(angle)
    g_n = parameters.gravity * math.cos(angle)

    z_lower_dimensional = np.linspace(0.0, lower.depth, samples)
    y_upper_dimensional = np.linspace(0.0, upper.depth, samples)

    interface_traction = upper.density * g_s * upper.depth
    bed_traction = g_s * (
        lower.density * lower.depth + upper.density * upper.depth
    )

    traction_lower = bed_traction - lower.density * g_s * z_lower_dimensional
    traction_upper = interface_traction - upper.density * g_s * y_upper_dimensional
    traction_lower[-1] = interface_traction
    traction_upper[-1] = 0.0

    shear_lower_dimensional = invert_regularized_stress(
        traction_lower,
        lower,
        parameters.regularization_rate,
    )
    shear_upper_dimensional = invert_regularized_stress(
        traction_upper,
        upper,
        parameters.regularization_rate,
    )

    velocity_lower_dimensional = cumulative_trapezoid(
        shear_lower_dimensional,
        z_lower_dimensional,
        initial=0.0,
    )
    interface_velocity_dimensional = float(velocity_lower_dimensional[-1])
    velocity_upper_dimensional = interface_velocity_dimensional + cumulative_trapezoid(
        shear_upper_dimensional,
        y_upper_dimensional,
        initial=0.0,
    )

    velocity_scale = float(
        trapezoid(velocity_lower_dimensional, z_lower_dimensional) / lower.depth
    )
    if not np.isfinite(velocity_scale) or velocity_scale <= 0.0:
        raise RuntimeError("The computed lower-layer mean velocity is not positive")

    H = lower.depth
    U = velocity_scale
    depth_ratio = upper.depth / lower.depth
    density_ratio = upper.density / lower.density
    froude = U / math.sqrt(g_n * H)

    consistency_number_lower = (
        lower.consistency
        * U ** (lower.flow_index - 2.0)
        / (lower.density * H ** lower.flow_index)
    )
    consistency_number_upper = (
        upper.consistency
        * U ** (upper.flow_index - 2.0)
        / (lower.density * H ** upper.flow_index)
    )
    epsilon = parameters.regularization_rate * H / U

    z_lower = z_lower_dimensional / H
    z_upper = 1.0 + y_upper_dimensional / H
    velocity_lower = velocity_lower_dimensional / U
    velocity_upper = velocity_upper_dimensional / U
    shear_lower = shear_lower_dimensional * H / U
    shear_upper = shear_upper_dimensional * H / U

    pressure_upper = (
        upper.density
        * g_n
        * (upper.depth - y_upper_dimensional)
        / (lower.density * U * U)
    )
    pressure_lower = (
        upper.density * g_n * upper.depth
        + lower.density * g_n * (lower.depth - z_lower_dimensional)
    ) / (lower.density * U * U)

    return BaseFlow(
        parameters=parameters,
        velocity_scale=velocity_scale,
        normal_gravity=g_n,
        downslope_gravity=g_s,
        froude=froude,
        depth_ratio=depth_ratio,
        density_ratio=density_ratio,
        consistency_number_lower=consistency_number_lower,
        consistency_number_upper=consistency_number_upper,
        epsilon_lower=epsilon,
        epsilon_upper=epsilon,
        interface_velocity=interface_velocity_dimensional / U,
        surface_velocity=float(velocity_upper_dimensional[-1] / U),
        lower=LayerProfile(
            z=z_lower,
            velocity=velocity_lower,
            shear=shear_lower,
            pressure=pressure_lower,
        ),
        upper=LayerProfile(
            z=z_upper,
            velocity=velocity_upper,
            shear=shear_upper,
            pressure=pressure_upper,
        ),
    )


# ---------------------------------------------------------------------------
# 4. Chebyshev collocation in each layer
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpectralGrid:
    z: FloatArray
    D1: ComplexArray
    D2: ComplexArray
    D3: ComplexArray
    D4: ComplexArray
    identity: ComplexArray


def chebyshev_grid(order: int, lower_bound: float, upper_bound: float) -> SpectralGrid:
    """Return ascending Chebyshev--Gauss--Lobatto nodes and derivatives."""

    N = int(order)
    if N < 1:
        raise ValueError("order must be positive")
    if upper_bound <= lower_bound:
        raise ValueError("upper_bound must exceed lower_bound")

    j = np.arange(N + 1)
    x_descending = np.cos(np.pi * j / N)
    weights = np.ones(N + 1)
    weights[0] = 2.0
    weights[-1] = 2.0
    weights *= (-1.0) ** j

    X = np.tile(x_descending, (N + 1, 1))
    dX = X.T - X
    D_descending = np.outer(weights, 1.0 / weights) / (
        dX + np.eye(N + 1)
    )
    D_descending -= np.diag(np.sum(D_descending, axis=1))

    reversal = np.eye(N + 1)[::-1]
    x = x_descending[::-1]
    D_reference = reversal @ D_descending @ reversal

    scale = 2.0 / (upper_bound - lower_bound)
    D1 = scale * D_reference
    D2 = D1 @ D1
    D3 = D2 @ D1
    D4 = D2 @ D2
    z = (
        0.5 * (lower_bound + upper_bound)
        + 0.5 * (upper_bound - lower_bound) * x
    )

    return SpectralGrid(
        z=z,
        D1=D1.astype(complex),
        D2=D2.astype(complex),
        D3=D3.astype(complex),
        D4=D4.astype(complex),
        identity=np.eye(N + 1, dtype=complex),
    )


@dataclass(frozen=True)
class CollocationLayer:
    grid: SpectralGrid
    velocity: FloatArray
    shear: FloatArray
    curvature: FloatArray
    apparent_viscosity: FloatArray
    tangent_viscosity: FloatArray
    density_ratio: float


def _interpolate(source_z: FloatArray, source_values: FloatArray, target_z: FloatArray) -> FloatArray:
    return np.interp(target_z, source_z, source_values)


def build_collocation_layer(
    base: BaseFlow,
    order: int,
    which: str,
) -> CollocationLayer:
    """Sample one base profile on its Chebyshev grid.

    The base curvature U'' is evaluated from the exact traction gradient,

        mu_t U'' = dT/dz,

    rather than by differentiating an interpolated array.  This is both more
    accurate and easier to audit.
    """

    if which == "lower":
        profile = base.lower
        fluid = base.parameters.lower
        grid = chebyshev_grid(order, 0.0, 1.0)
        Lambda = base.consistency_number_lower
        epsilon = base.epsilon_lower
        density_ratio = 1.0
        traction_gradient = -base.parameters.slope_tan / base.froude**2
    elif which == "upper":
        profile = base.upper
        fluid = base.parameters.upper
        grid = chebyshev_grid(order, 1.0, 1.0 + base.depth_ratio)
        Lambda = base.consistency_number_upper
        epsilon = base.epsilon_upper
        density_ratio = base.density_ratio
        traction_gradient = (
            -base.density_ratio
            * base.parameters.slope_tan
            / base.froude**2
        )
    else:
        raise ValueError("which must be 'lower' or 'upper'")

    velocity = _interpolate(profile.z, profile.velocity, grid.z)
    shear = _interpolate(profile.z, profile.shear, grid.z)
    mu, mu_t = apparent_and_tangent_viscosity(
        shear,
        Lambda,
        fluid.flow_index,
        epsilon,
    )
    if not np.all(np.isfinite(mu)) or not np.all(np.isfinite(mu_t)):
        raise FloatingPointError(
            f"Non-finite viscosity in the {which} layer.  Use a positive "
            "regularization_rate for a shear-thinning upper layer."
        )

    curvature = traction_gradient / mu_t

    return CollocationLayer(
        grid=grid,
        velocity=velocity,
        shear=shear,
        curvature=curvature,
        apparent_viscosity=mu,
        tangent_viscosity=mu_t,
        density_ratio=density_ratio,
    )


# ---------------------------------------------------------------------------
# 5. Matrix pencil A q = c B q
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PencilLayout:
    points_per_layer: int
    lower_slice: slice
    upper_slice: slice
    eta_index: int
    zeta_index: int
    wall_psi_row: int
    wall_u_row: int
    interface_w_row: int
    interface_u_row: int
    interface_t_row: int
    interface_n_row: int
    free_t_row: int
    free_n_row: int
    interface_kinematic_row: int
    free_kinematic_row: int


def pencil_layout(order: int) -> PencilLayout:
    n = order + 1
    return PencilLayout(
        points_per_layer=n,
        lower_slice=slice(0, n),
        upper_slice=slice(n, 2 * n),
        eta_index=2 * n,
        zeta_index=2 * n + 1,
        wall_psi_row=0,
        wall_u_row=1,
        interface_w_row=n - 2,
        interface_u_row=n - 1,
        interface_t_row=n,
        interface_n_row=n + 1,
        free_t_row=2 * n - 2,
        free_n_row=2 * n - 1,
        interface_kinematic_row=2 * n,
        free_kinematic_row=2 * n + 1,
    )


@dataclass(frozen=True)
class MatrixPencil:
    A: ComplexArray
    B: ComplexArray
    layout: PencilLayout
    lower: CollocationLayer
    upper: CollocationLayer
    wavenumber: float


@dataclass(frozen=True)
class PreparedProblem:
    """Base-flow data sampled once for repeated wavenumber calculations."""

    base: BaseFlow
    order: int
    layout: PencilLayout
    lower: CollocationLayer
    upper: CollocationLayer


def prepare_problem(base: BaseFlow, order: int) -> PreparedProblem:
    """Precompute all k-independent spectral data."""

    return PreparedProblem(
        base=base,
        order=int(order),
        layout=pencil_layout(int(order)),
        lower=build_collocation_layer(base, int(order), "lower"),
        upper=build_collocation_layer(base, int(order), "upper"),
    )


def _interior_operators(
    layer: CollocationLayer,
    wavenumber: float,
) -> tuple[ComplexArray, ComplexArray]:
    """Return the interior A and B operators for phase speed c."""

    k = wavenumber
    D = layer.grid.D1
    D2 = layer.grid.D2
    I = layer.grid.identity
    L_plus = D2 + k * k * I
    L_minus = D2 - k * k * I

    mu = np.diag(layer.apparent_viscosity)
    mu_t = np.diag(layer.tangent_viscosity)
    U = np.diag(layer.velocity)
    Udd = np.diag(layer.curvature)
    r = layer.density_ratio

    viscous = L_plus @ mu_t @ L_plus - 4.0 * k * k * D @ mu @ D

    # Generalized Orr--Sommerfeld equation with sigma = -i k c:
    #
    # [V - i k r U L_- + i k r U''] psi = c [-i k r L_-] psi.
    A = viscous - 1j * k * r * U @ L_minus + 1j * k * r * Udd
    B = -1j * k * r * L_minus
    return A, B


def _ik_normal_traction_without_c(
    layer: CollocationLayer,
    wavenumber: float,
    endpoint: int,
) -> ComplexArray:
    """Return the row for i k N with the c-dependent term removed.

    For the streamfunction perturbation u = D psi, w = -i k psi,

      i k N = -i k c r D psi
              + i k r U D psi - i k r U' psi
              + 4 k^2 mu D psi - D[mu_t (D^2+k^2) psi].
    """

    k = wavenumber
    D = layer.grid.D1
    I = layer.grid.identity
    L_plus = layer.grid.D2 + k * k * I
    mu = np.diag(layer.apparent_viscosity)
    mu_t = np.diag(layer.tangent_viscosity)
    r = layer.density_ratio

    return (
        1j * k * r * layer.velocity[endpoint] * D[endpoint, :]
        - 1j * k * r * layer.shear[endpoint] * I[endpoint, :]
        + 4.0 * k * k * mu[endpoint, endpoint] * D[endpoint, :]
        - (D @ mu_t @ L_plus)[endpoint, :]
    )


def assemble_prepared_pencil(
    prepared: PreparedProblem,
    wavenumber: float,
) -> MatrixPencil:
    """Assemble A q = c B q using precomputed spectral data."""

    k = float(wavenumber)
    if not np.isfinite(k) or k <= 0.0:
        raise ValueError("wavenumber must be positive and finite")

    base = prepared.base
    lower = prepared.lower
    upper = prepared.upper
    layout = prepared.layout
    n = layout.points_per_layer
    size = 2 * n + 2

    A = np.zeros((size, size), dtype=complex)
    B = np.zeros((size, size), dtype=complex)

    # Interior rows.  The first two and last two rows of each domain are
    # reserved for boundary and interface conditions.
    for layer, column_slice, offset in (
        (lower, layout.lower_slice, 0),
        (upper, layout.upper_slice, n),
    ):
        A_interior, B_interior = _interior_operators(layer, k)
        for local_row in range(2, n - 2):
            global_row = offset + local_row
            A[global_row, column_slice] = A_interior[local_row, :]
            B[global_row, column_slice] = B_interior[local_row, :]

    lower_bottom = 0
    lower_top = n - 1
    upper_bottom = 0
    upper_top = n - 1

    # Rigid wall: psi_l(0)=0 and D psi_l(0)=0.
    A[layout.wall_psi_row, lower_bottom] = 1.0
    A[layout.wall_u_row, layout.lower_slice] = lower.grid.D1[lower_bottom, :]

    # Interface normal-velocity continuity: psi_l(1)=psi_u(1).
    A[layout.interface_w_row, lower_top] = 1.0
    A[layout.interface_w_row, n + upper_bottom] = -1.0

    # Interface tangential-velocity continuity at the displaced interface.
    A[layout.interface_u_row, layout.lower_slice] = lower.grid.D1[lower_top, :]
    A[layout.interface_u_row, layout.upper_slice] = -upper.grid.D1[upper_bottom, :]
    A[layout.interface_u_row, layout.eta_index] = (
        lower.shear[lower_top] - upper.shear[upper_bottom]
    )

    L_plus_lower = lower.grid.D2 + k * k * lower.grid.identity
    L_plus_upper = upper.grid.D2 + k * k * upper.grid.identity

    # Interface tangential traction:
    # T_l - T_u - (1-rho_r) S0/Fr^2 eta = 0.
    A[layout.interface_t_row, layout.lower_slice] = (
        lower.tangent_viscosity[lower_top] * L_plus_lower[lower_top, :]
    )
    A[layout.interface_t_row, layout.upper_slice] = (
        -upper.tangent_viscosity[upper_bottom] * L_plus_upper[upper_bottom, :]
    )
    A[layout.interface_t_row, layout.eta_index] = (
        -(1.0 - base.density_ratio)
        * base.parameters.slope_tan
        / base.froude**2
    )

    # Interface normal traction, multiplied by i k.
    A[layout.interface_n_row, layout.lower_slice] = (
        _ik_normal_traction_without_c(lower, k, lower_top)
    )
    A[layout.interface_n_row, layout.upper_slice] = (
        -_ik_normal_traction_without_c(upper, k, upper_bottom)
    )
    A[layout.interface_n_row, layout.eta_index] = (
        1j * k * (1.0 - base.density_ratio) / base.froude**2
    )
    B[layout.interface_n_row, layout.lower_slice] = (
        1j * k * lower.density_ratio * lower.grid.D1[lower_top, :]
    )
    B[layout.interface_n_row, layout.upper_slice] = (
        -1j * k * upper.density_ratio * upper.grid.D1[upper_bottom, :]
    )

    # Stress-free free surface.
    A[layout.free_t_row, layout.upper_slice] = (
        upper.tangent_viscosity[upper_top] * L_plus_upper[upper_top, :]
    )
    A[layout.free_t_row, layout.zeta_index] = (
        -base.density_ratio * base.parameters.slope_tan / base.froude**2
    )

    A[layout.free_n_row, layout.upper_slice] = (
        _ik_normal_traction_without_c(upper, k, upper_top)
    )
    A[layout.free_n_row, layout.zeta_index] = (
        1j * k * base.density_ratio / base.froude**2
    )
    B[layout.free_n_row, layout.upper_slice] = (
        1j * k * upper.density_ratio * upper.grid.D1[upper_top, :]
    )

    # Kinematic conditions after division by -i k:
    # psi_I + U_I eta = c eta,
    # psi_s + U_s zeta = c zeta.
    A[layout.interface_kinematic_row, lower_top] = 1.0
    A[layout.interface_kinematic_row, layout.eta_index] = base.interface_velocity
    B[layout.interface_kinematic_row, layout.eta_index] = 1.0

    A[layout.free_kinematic_row, n + upper_top] = 1.0
    A[layout.free_kinematic_row, layout.zeta_index] = base.surface_velocity
    B[layout.free_kinematic_row, layout.zeta_index] = 1.0

    return MatrixPencil(
        A=A,
        B=B,
        layout=layout,
        lower=lower,
        upper=upper,
        wavenumber=k,
    )


def assemble_phase_speed_pencil(
    base: BaseFlow,
    wavenumber: float,
    order: int,
) -> MatrixPencil:
    """Convenience wrapper for one isolated wavenumber calculation."""

    return assemble_prepared_pencil(prepare_problem(base, order), wavenumber)


# ---------------------------------------------------------------------------
# 6. Conditioning, eigenvalue filtering, and mode tracking
# ---------------------------------------------------------------------------


def equilibrate_rows(
    A: ComplexArray,
    B: ComplexArray,
) -> tuple[ComplexArray, ComplexArray, FloatArray]:
    """Scale every equation to unit combined row norm.

    Multiplying a row of both A and B by the same nonzero number leaves the
    exact generalized eigenvalues unchanged.  It substantially improves the
    numerical resolution of weak long-wave modes in this singular pencil.
    """

    scales = np.sqrt(
        np.sum(np.abs(A) ** 2, axis=1)
        + np.sum(np.abs(B) ** 2, axis=1)
    )
    scales = np.where(scales > 0.0, scales, 1.0)
    return A / scales[:, None], B / scales[:, None], scales


def generalized_residual(
    A: ComplexArray,
    B: ComplexArray,
    eigenvalue: complex,
    eigenvector: ComplexArray,
) -> float:
    """Return a scale-invariant residual for A v = lambda B v."""

    numerator = np.linalg.norm(A @ eigenvector - eigenvalue * (B @ eigenvector))
    denominator = (
        np.linalg.norm(A) + abs(eigenvalue) * np.linalg.norm(B)
    ) * np.linalg.norm(eigenvector)
    return float(numerator / max(denominator, np.finfo(float).tiny))


def _phase_normalize(vector: ComplexArray, anchor_index: int) -> ComplexArray:
    """Normalize a complex eigenvector and make one anchor real and positive."""

    v = np.asarray(vector, dtype=complex).copy()
    norm = np.linalg.norm(v)
    v /= max(norm, np.finfo(float).tiny)

    anchor = v[anchor_index]
    if abs(anchor) > 1.0e-14:
        v *= np.exp(-1j * np.angle(anchor))
        if v[anchor_index].real < 0.0:
            v *= -1.0
    return v


@dataclass(frozen=True)
class Eigenmode:
    phase_speed: complex
    sigma: complex
    growth_rate: float
    celerity: float
    residual: float
    eta: complex
    zeta: complex
    surface_to_interface: complex
    vector: ComplexArray

    @property
    def amplitude_ratio(self) -> float:
        return float(abs(self.surface_to_interface))

    @property
    def phase_difference_degrees(self) -> float:
        return float(np.degrees(np.angle(self.surface_to_interface)))


def solve_spectrum(
    base: BaseFlow,
    wavenumber: float,
    options: SolverOptions,
    prepared: PreparedProblem | None = None,
) -> tuple[list[Eigenmode], MatrixPencil, FloatArray]:
    """Solve the complete collocation spectrum at one wavenumber."""

    options.validate()
    if prepared is None:
        prepared = prepare_problem(base, options.chebyshev_order)
    elif prepared.order != options.chebyshev_order:
        raise ValueError("prepared spectral order does not match SolverOptions")
    pencil = assemble_prepared_pencil(prepared, wavenumber)
    A_scaled, B_scaled, row_scales = equilibrate_rows(pencil.A, pencil.B)

    # Homogeneous generalized eigenvalues alpha/beta expose beta=0 infinite
    # eigenvalues directly, rather than relying on very large floating values.
    homogeneous, vectors = eig(
        A_scaled,
        B_scaled,
        right=True,
        homogeneous_eigvals=True,
        check_finite=False,
    )
    alpha, beta = homogeneous

    modes: list[Eigenmode] = []
    machine = np.finfo(float).eps
    for column in range(alpha.size):
        scale = max(abs(alpha[column]), abs(beta[column]), 1.0)
        if abs(beta[column]) <= 100.0 * machine * scale:
            continue

        phase_speed = complex(alpha[column] / beta[column])
        if not np.isfinite(phase_speed.real) or not np.isfinite(phase_speed.imag):
            continue
        if abs(phase_speed) > options.phase_speed_limit:
            continue

        raw_vector = vectors[:, column]
        residual = generalized_residual(
            pencil.A,
            pencil.B,
            phase_speed,
            raw_vector,
        )
        if residual > options.residual_tolerance:
            continue

        layout = pencil.layout
        # Prefer the interface displacement as phase anchor.  If it vanishes,
        # fall back to the free-surface displacement and then to the largest
        # component of the eigenvector.
        if abs(raw_vector[layout.eta_index]) > 1.0e-12:
            anchor = layout.eta_index
        elif abs(raw_vector[layout.zeta_index]) > 1.0e-12:
            anchor = layout.zeta_index
        else:
            anchor = int(np.argmax(np.abs(raw_vector)))
        vector = _phase_normalize(raw_vector, anchor)

        eta = complex(vector[layout.eta_index])
        zeta = complex(vector[layout.zeta_index])
        if abs(eta) > 1.0e-13:
            ratio = zeta / eta
        else:
            ratio = complex(np.nan, np.nan)

        sigma = -1j * wavenumber * phase_speed
        modes.append(
            Eigenmode(
                phase_speed=phase_speed,
                sigma=sigma,
                growth_rate=float(sigma.real),
                celerity=float(phase_speed.real),
                residual=residual,
                eta=eta,
                zeta=zeta,
                surface_to_interface=complex(ratio),
                vector=vector,
            )
        )

    modes.sort(key=lambda mode: abs(mode.sigma))
    return modes, pencil, row_scales


def _surface_amplitude(mode: Eigenmode) -> float:
    return float(math.sqrt(abs(mode.eta) ** 2 + abs(mode.zeta) ** 2))


def select_long_wave_pair(modes: Sequence[Eigenmode]) -> dict[str, Eigenmode]:
    """Select and label the two hydrodynamic roots at a small wavenumber.

    Hydrodynamic roots satisfy sigma -> 0 as k -> 0, whereas viscous/shear
    modes retain finite damping.  Therefore the two downstream roots with the
    smallest |sigma| and non-negligible moving-boundary amplitudes are the
    appropriate pair.  The slower root is labeled interface/varicose and the
    faster root surface/zigzag.
    """

    candidates = [
        mode
        for mode in modes
        if mode.celerity > 0.0
        and _surface_amplitude(mode) > 1.0e-8
        and np.isfinite(mode.amplitude_ratio)
    ]
    candidates.sort(key=lambda mode: abs(mode.sigma))
    if len(candidates) < 2:
        raise RuntimeError("Could not identify two long-wave hydrodynamic roots")

    pair = sorted(candidates[:2], key=lambda mode: mode.celerity)
    return {
        "interface_varicose": pair[0],
        "surface_zigzag": pair[1],
    }


def _eigenvector_overlap(previous: Eigenmode, current: Eigenmode) -> float:
    overlap = abs(np.vdot(previous.vector, current.vector))
    return float(min(max(overlap, 0.0), 1.0))


def _continuation_cost(previous: Eigenmode, current: Eigenmode) -> float:
    """Distance combining phase speed, eigenfunction, and boundary geometry."""

    speed_distance = abs(current.phase_speed - previous.phase_speed) / (
        1.0 + abs(previous.phase_speed)
    )
    vector_distance = 1.0 - _eigenvector_overlap(previous, current)

    old_ratio = previous.amplitude_ratio
    new_ratio = current.amplitude_ratio
    if old_ratio > 0.0 and new_ratio > 0.0 and np.isfinite(new_ratio):
        ratio_distance = abs(math.log(new_ratio / old_ratio))
    else:
        ratio_distance = 0.0

    return float(speed_distance + 0.20 * vector_distance + 0.02 * ratio_distance)


def continue_pair(
    previous: dict[str, Eigenmode],
    modes: Sequence[Eigenmode],
) -> dict[str, Eigenmode]:
    """Continue both branches while enforcing distinct eigenvalues."""

    candidates = [
        mode
        for mode in modes
        if -5.0 < mode.celerity < 20.0
        and _surface_amplitude(mode) > 1.0e-9
        and np.isfinite(mode.amplitude_ratio)
    ]
    if len(candidates) < 2:
        raise RuntimeError("Too few candidate roots for branch continuation")

    # Only a small neighbourhood of each previous root is relevant.  Limiting
    # the candidate sets keeps continuation fast and makes the logic easier to
    # inspect than an all-against-all comparison of the complete spectrum.
    interface_candidates = sorted(
        candidates,
        key=lambda mode: _continuation_cost(
            previous["interface_varicose"], mode
        ),
    )[:12]
    surface_candidates = sorted(
        candidates,
        key=lambda mode: _continuation_cost(
            previous["surface_zigzag"], mode
        ),
    )[:12]

    best_cost = math.inf
    best_pair: tuple[Eigenmode, Eigenmode] | None = None
    for interface_mode in interface_candidates:
        for surface_mode in surface_candidates:
            if interface_mode is surface_mode:
                continue
            total = (
                _continuation_cost(previous["interface_varicose"], interface_mode)
                + _continuation_cost(previous["surface_zigzag"], surface_mode)
            )
            if total < best_cost:
                best_cost = total
                best_pair = (interface_mode, surface_mode)

    if best_pair is None:
        raise RuntimeError("Could not continue the two hydrodynamic branches")
    return {
        "interface_varicose": best_pair[0],
        "surface_zigzag": best_pair[1],
    }


def _tracking_grid(
    requested_wavenumbers: FloatArray,
    seed_wavenumber: float,
) -> tuple[FloatArray, NDArray[np.int64]]:
    """Add a short hidden continuation path before the first requested k."""

    requested = np.asarray(requested_wavenumbers, dtype=float)
    if requested.ndim != 1 or requested.size == 0:
        raise ValueError("requested_wavenumbers must be a nonempty one-dimensional array")
    if np.any(~np.isfinite(requested)) or np.any(requested <= 0.0):
        raise ValueError("All requested wavenumbers must be positive and finite")
    if np.any(np.diff(requested) <= 0.0):
        raise ValueError("requested_wavenumbers must be strictly increasing")

    if requested[0] <= seed_wavenumber * (1.0 + 1.0e-12):
        return requested, np.arange(requested.size, dtype=int)

    prefix = np.geomspace(seed_wavenumber, requested[0], 10, endpoint=False)
    combined = np.concatenate([prefix, requested])
    requested_indices = np.arange(prefix.size, combined.size, dtype=int)
    return combined, requested_indices


def track_hydrodynamic_modes(
    base: BaseFlow,
    requested_wavenumbers: FloatArray,
    options: SolverOptions,
) -> tuple[dict[str, list[Eigenmode]], list[dict[str, float]]]:
    """Track the interface and surface branches over increasing k."""

    all_wavenumbers, requested_indices = _tracking_grid(
        requested_wavenumbers,
        options.seed_wavenumber,
    )
    requested_index_set = set(int(i) for i in requested_indices)

    tracks: dict[str, list[Eigenmode]] = {
        "interface_varicose": [],
        "surface_zigzag": [],
    }
    diagnostics: list[dict[str, float]] = []
    previous: dict[str, Eigenmode] | None = None
    prepared = prepare_problem(base, options.chebyshev_order)

    for grid_index, wavenumber in enumerate(all_wavenumbers):
        modes, pencil, row_scales = solve_spectrum(
            base,
            float(wavenumber),
            options,
            prepared=prepared,
        )
        if previous is None:
            selected = select_long_wave_pair(modes)
        else:
            selected = continue_pair(previous, modes)
        previous = selected

        diagnostics.append(
            {
                "wavenumber": float(wavenumber),
                "row_norm_spread": float(np.max(row_scales) / np.min(row_scales)),
                "finite_modes_retained": float(len(modes)),
                "matrix_size": float(pencil.A.shape[0]),
            }
        )

        if grid_index in requested_index_set:
            for name in tracks:
                tracks[name].append(selected[name])

    expected = len(requested_wavenumbers)
    if any(len(track) != expected for track in tracks.values()):
        raise RuntimeError("Internal branch-tracking length mismatch")
    return tracks, diagnostics


# ---------------------------------------------------------------------------
# 7. Convenience calculations used by the runner and verification suite
# ---------------------------------------------------------------------------


def convergence_at_wavenumber(
    parameters: CaseParameters,
    wavenumber: float,
    orders: Iterable[int] = (32, 40, 52, 64, 76),
    seed_wavenumber: float = 5.0e-4,
) -> list[dict[str, float]]:
    """Evaluate both tracked branches at one k for several spectral orders."""

    base = build_base_flow(parameters)
    target = float(wavenumber)
    if target <= 0.0:
        raise ValueError("wavenumber must be positive")

    rows: list[dict[str, float]] = []
    for order in orders:
        options = SolverOptions(
            chebyshev_order=int(order),
            residual_tolerance=1.0e-5,
            seed_wavenumber=min(seed_wavenumber, 0.5 * target),
        )
        # A short geometric path is sufficient for reliable branch identity.
        ks = np.geomspace(options.seed_wavenumber, target, 18)
        tracks, _ = track_hydrodynamic_modes(base, ks, options)

        row: dict[str, float] = {"order": float(order)}
        for name, track in tracks.items():
            mode = track[-1]
            row[f"{name}_growth"] = mode.growth_rate
            row[f"{name}_celerity"] = mode.celerity
            row[f"{name}_residual"] = mode.residual
        rows.append(row)
    return rows


def neutral_crossings(
    wavenumbers: FloatArray,
    growth_rates: FloatArray,
) -> list[float]:
    """Return linearly interpolated zeros of a sampled growth-rate curve."""

    k = np.asarray(wavenumbers, dtype=float)
    growth = np.asarray(growth_rates, dtype=float)
    if k.shape != growth.shape:
        raise ValueError("wavenumbers and growth_rates must have the same shape")

    roots: list[float] = []
    for index in range(k.size - 1):
        g0 = growth[index]
        g1 = growth[index + 1]
        if g0 == 0.0:
            roots.append(float(k[index]))
        elif g0 * g1 < 0.0:
            roots.append(
                float(k[index] - g0 * (k[index + 1] - k[index]) / (g1 - g0))
            )
    if growth.size and growth[-1] == 0.0:
        roots.append(float(k[-1]))
    return roots
