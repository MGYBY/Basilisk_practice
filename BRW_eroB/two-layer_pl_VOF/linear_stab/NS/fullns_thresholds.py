#!/usr/bin/env python3
"""Neutral thresholds for the two-layer power-law full-NS eigenproblem.

This module builds on :mod:`fullns_stability`.  It adds the two pieces that
are needed for stability maps:

1. a dimensionless normal-flow family in which the lower-layer consistency
   coefficient is removed by the steady-uniform compatibility condition;
2. branch-specific neutral searches in the infinitely long-wave limit and at
   a finite prescribed wavelength.

The nondimensionalization uses the lower normal-flow depth ``H_l`` and lower
normal-flow mean velocity ``U_l``.  The lower generalized Reynolds number is

    Re_l = 1 / Lambda_l,

and the compatibility ratio

    chi_l = (S0/Fr_l**2) / Lambda_l = S0*Re_l/Fr_l**2

is fixed by the lower mean-velocity normalization.  Hence

    S0 = chi_l*Fr_l**2/Re_l.

For a regularized power law, epsilon_l and epsilon_u are held fixed as
*dimensionless* transition shear rates.  Holding a dimensional transition
rate fixed while varying U_l defines a different parameter family and is not
silently mixed with the calculations here.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Iterable, Literal, Sequence
import math

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq

from fullns_stability import (
    BaseFlow,
    CaseParameters,
    Eigenmode,
    FluidLayer,
    LayerProfile,
    SolverOptions,
    continue_pair,
    invert_regularized_stress,
    select_long_wave_pair,
    solve_spectrum,
    track_hydrodynamic_modes,
)


FloatArray = NDArray[np.float64]
BranchName = Literal["interface_varicose", "surface_zigzag"]
WaveKind = Literal[
    "fixed_kH",
    "fixed_lambda_over_H",
    "fixed_slope_scaled_wavelength",
]


@dataclass(frozen=True)
class DimensionlessFamily:
    """Independent parameters of one dimensionless normal-flow family.

    ``consistency_ratio`` is ``kappa_K = Lambda_u/Lambda_l``.  It is the
    dimensionless consistency contrast; when the two flow indices differ it
    is not simply the raw dimensional ratio ``K_u/K_l``.
    """

    depth_ratio: float = 1.0
    density_ratio: float = 0.8
    n_lower: float = 0.4
    n_upper: float = 0.8
    consistency_ratio: float = 1.7336759276597786
    epsilon_lower: float = 0.02233444771394924
    epsilon_upper: float = 0.02233444771394924

    def validate(self) -> None:
        positive = {
            "depth_ratio": self.depth_ratio,
            "density_ratio": self.density_ratio,
            "n_lower": self.n_lower,
            "n_upper": self.n_upper,
            "consistency_ratio": self.consistency_ratio,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite; got {value!r}")
        for name, value in {
            "epsilon_lower": self.epsilon_lower,
            "epsilon_upper": self.epsilon_upper,
        }.items():
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be nonnegative and finite")
        if self.n_upper < 1.0 and self.epsilon_upper == 0.0:
            raise ValueError(
                "An ideal shear-thinning upper layer is singular at its "
                "stress-free surface.  Use epsilon_upper > 0 for numerical "
                "full-NS spectra."
            )


@dataclass(frozen=True)
class NormalFlowShape:
    """Normal-flow profiles that are independent of Fr, Re, and S0."""

    family: DimensionlessFamily
    compatibility_ratio: float
    lower_z: FloatArray
    upper_z: FloatArray
    lower_shear: FloatArray
    upper_shear: FloatArray
    lower_velocity: FloatArray
    upper_velocity: FloatArray
    interface_velocity: float
    surface_velocity: float
    lower_mean_velocity: float

    @property
    def generalized_reynolds_relation(self) -> str:
        return "Re_l = chi_l*Fr_l^2/S0"


@dataclass(frozen=True)
class WaveSpecification:
    """Definition of the finite disturbance wavelength.

    ``fixed_kH``
        ``value`` is k_H = k_dimensional H_l.
    ``fixed_lambda_over_H``
        ``value`` is lambda_dimensional/H_l.
    ``fixed_slope_scaled_wavelength``
        ``value`` is S0*lambda_dimensional/H_l.  This is the convention used
        in Figure 2 of Yu & Chu (2024).
    """

    kind: WaveKind = "fixed_slope_scaled_wavelength"
    value: float = 4.72

    def validate(self) -> None:
        if self.kind not in {
            "fixed_kH",
            "fixed_lambda_over_H",
            "fixed_slope_scaled_wavelength",
        }:
            raise ValueError(f"Unsupported wave specification: {self.kind!r}")
        if not np.isfinite(self.value) or self.value <= 0.0:
            raise ValueError("WaveSpecification.value must be positive and finite")


@dataclass(frozen=True)
class LongWaveOptions:
    """Controls for extrapolating a finite-k neutral root to k=0."""

    froude_min: float = 0.10
    froude_max: float = 1.60
    scan_points: int = 25
    reference_reynolds: float = 20.0
    wavenumbers: tuple[float, ...] = (
        3.0e-4,
        5.0e-4,
        8.0e-4,
        1.2e-3,
        2.0e-3,
    )
    polynomial_degree: int = 1

    def validate(self) -> None:
        if not (0.0 < self.froude_min < self.froude_max):
            raise ValueError("Require 0 < froude_min < froude_max")
        if self.scan_points < 9:
            raise ValueError("scan_points should be at least 9")
        if self.reference_reynolds <= 0.0:
            raise ValueError("reference_reynolds must be positive")
        if len(self.wavenumbers) < 3:
            raise ValueError("At least three small wavenumbers are required")
        if any((not np.isfinite(k) or k <= 0.0) for k in self.wavenumbers):
            raise ValueError("All long-wave wavenumbers must be positive")
        if any(b <= a for a, b in zip(self.wavenumbers, self.wavenumbers[1:])):
            raise ValueError("Long-wave wavenumbers must be strictly increasing")
        if self.polynomial_degree < 1:
            raise ValueError("polynomial_degree must be positive")
        if self.polynomial_degree >= len(self.wavenumbers):
            raise ValueError("polynomial_degree must be less than number of samples")


@dataclass(frozen=True)
class LongWaveThreshold:
    branch: BranchName
    critical_froude: float
    finite_k_roots: tuple[float, ...]
    wavenumbers: tuple[float, ...]
    fit_coefficients_in_k2: tuple[float, ...]
    maximum_fit_residual: float
    reference_reynolds: float
    compatibility_ratio: float
    reduced_growth_at_root: float
    long_wave_celerity: float
    converged: bool
    quality_message: str


@dataclass(frozen=True)
class NeutralPoint:
    branch: BranchName
    reynolds: float
    froude: float
    slope: float
    wavenumber: float
    wavelength_over_H: float
    slope_scaled_wavelength: float
    growth_rate: float
    celerity: float
    residual: float
    root_number: int


@dataclass(frozen=True)
class NeutralCurveOptions:
    froude_min: float = 0.10
    froude_max: float = 1.60
    scan_points: int = 36
    root_xtol: float = 3.0e-7
    continuation_steps: int = 12

    def validate(self) -> None:
        if not (0.0 < self.froude_min < self.froude_max):
            raise ValueError("Require 0 < froude_min < froude_max")
        if self.scan_points < 10:
            raise ValueError("scan_points should be at least 10")
        if self.root_xtol <= 0.0:
            raise ValueError("root_xtol must be positive")
        if self.continuation_steps < 4:
            raise ValueError("continuation_steps should be at least 4")


# ---------------------------------------------------------------------------
# Compatibility and dimensionless base flow
# ---------------------------------------------------------------------------


def _ideal_compatibility_ratio(family: DimensionlessFamily) -> float:
    """Closed ideal-power-law value chi_l = I_l**(-n_l)."""

    z = np.linspace(0.0, 1.0, 20001)
    A = 1.0 + family.density_ratio * family.depth_ratio
    integral = float(
        trapezoid(
            (1.0 - z) * np.power(A - z, 1.0 / family.n_lower),
            z,
        )
    )
    return integral ** (-family.n_lower)


@lru_cache(maxsize=64)
def _cached_shape(family: DimensionlessFamily, samples: int) -> NormalFlowShape:
    family.validate()
    if samples < 401:
        raise ValueError("samples should be at least 401")

    z_lower = np.linspace(0.0, 1.0, int(samples))
    y_upper = np.linspace(0.0, family.depth_ratio, int(samples))
    A = 1.0 + family.density_ratio * family.depth_ratio

    lower_unit = FluidLayer(1.0, 1.0, 1.0, family.n_lower)
    upper_unit = FluidLayer(family.depth_ratio, family.density_ratio, 1.0, family.n_upper)

    if family.epsilon_lower == 0.0:
        chi = _ideal_compatibility_ratio(family)
    else:
        def mean_minus_one(chi_trial: float) -> float:
            reduced_traction = chi_trial * (A - z_lower)
            shear = invert_regularized_stress(
                reduced_traction,
                lower_unit,
                family.epsilon_lower,
            )
            mean_velocity = float(trapezoid((1.0 - z_lower) * shear, z_lower))
            return mean_velocity - 1.0

        upper_chi = max(_ideal_compatibility_ratio(family), 1.0e-10)
        while mean_minus_one(upper_chi) < 0.0:
            upper_chi *= 2.0
            if upper_chi > 1.0e12:
                raise RuntimeError("Could not bracket the compatibility ratio chi_l")
        chi = float(brentq(mean_minus_one, 0.0, upper_chi, xtol=2.0e-13, rtol=2.0e-13))

    reduced_lower_traction = chi * (A - z_lower)
    reduced_upper_traction = (
        family.density_ratio
        * chi
        / family.consistency_ratio
        * (family.depth_ratio - y_upper)
    )

    shear_lower = invert_regularized_stress(
        reduced_lower_traction,
        lower_unit,
        family.epsilon_lower,
    )
    shear_upper = invert_regularized_stress(
        reduced_upper_traction,
        upper_unit,
        family.epsilon_upper,
    )

    velocity_lower = cumulative_trapezoid(shear_lower, z_lower, initial=0.0)
    interface_velocity = float(velocity_lower[-1])
    velocity_upper_local = interface_velocity + cumulative_trapezoid(
        shear_upper,
        y_upper,
        initial=0.0,
    )
    lower_mean = float(trapezoid(velocity_lower, z_lower))

    return NormalFlowShape(
        family=family,
        compatibility_ratio=chi,
        lower_z=z_lower,
        upper_z=1.0 + y_upper,
        lower_shear=shear_lower,
        upper_shear=shear_upper,
        lower_velocity=velocity_lower,
        upper_velocity=velocity_upper_local,
        interface_velocity=interface_velocity,
        surface_velocity=float(velocity_upper_local[-1]),
        lower_mean_velocity=lower_mean,
    )


def build_normal_flow_shape(
    family: DimensionlessFamily,
    samples: int = 8001,
) -> NormalFlowShape:
    """Return the compatible normal-flow shape for a dimensionless family."""

    return _cached_shape(family, int(samples))


def compatibility_slope(
    shape: NormalFlowShape,
    froude: float,
    reynolds: float,
) -> float:
    """Return S0 required by the steady-uniform compatibility condition."""

    if not np.isfinite(froude) or froude <= 0.0:
        raise ValueError("froude must be positive and finite")
    if not np.isfinite(reynolds) or reynolds <= 0.0:
        raise ValueError("reynolds must be positive and finite")
    return shape.compatibility_ratio * froude * froude / reynolds


def compatibility_reynolds(
    shape: NormalFlowShape,
    froude: float,
    slope: float,
) -> float:
    """Return Re_l required by the steady-uniform compatibility condition."""

    if not np.isfinite(slope) or slope <= 0.0:
        raise ValueError("slope must be positive and finite")
    return shape.compatibility_ratio * froude * froude / slope


def build_dimensionless_base_flow(
    shape: NormalFlowShape,
    froude: float,
    reynolds: float,
) -> BaseFlow:
    """Build a BaseFlow for one compatible (Fr_l, Re_l, S0) state."""

    family = shape.family
    family.validate()
    if not np.isfinite(froude) or froude <= 0.0:
        raise ValueError("froude must be positive and finite")
    if not np.isfinite(reynolds) or reynolds <= 0.0:
        raise ValueError("reynolds must be positive and finite")

    slope = compatibility_slope(shape, froude, reynolds)
    angle = math.atan(slope)
    normal_gravity = 1.0 / (froude * froude)
    gravity = normal_gravity / math.cos(angle)
    lambda_lower = 1.0 / reynolds
    lambda_upper = family.consistency_ratio / reynolds

    # This pseudo-dimensional object has H_l=U_l=rho_l=1.  Its only purpose
    # is to carry the physical labels and coefficients required by the
    # reusable spectral assembly.  The BaseFlow arrays are already fully
    # dimensionless.
    parameters = CaseParameters(
        gravity=gravity,
        slope_tan=slope,
        lower=FluidLayer(1.0, 1.0, lambda_lower, family.n_lower),
        upper=FluidLayer(
            family.depth_ratio,
            family.density_ratio,
            lambda_upper,
            family.n_upper,
        ),
        regularization_rate=family.epsilon_lower,
    )

    pressure_lower = (
        1.0 - shape.lower_z + family.density_ratio * family.depth_ratio
    ) / (froude * froude)
    pressure_upper = (
        family.density_ratio
        * (1.0 + family.depth_ratio - shape.upper_z)
        / (froude * froude)
    )

    return BaseFlow(
        parameters=parameters,
        velocity_scale=1.0,
        normal_gravity=normal_gravity,
        downslope_gravity=slope / (froude * froude),
        froude=froude,
        depth_ratio=family.depth_ratio,
        density_ratio=family.density_ratio,
        consistency_number_lower=lambda_lower,
        consistency_number_upper=lambda_upper,
        epsilon_lower=family.epsilon_lower,
        epsilon_upper=family.epsilon_upper,
        interface_velocity=shape.interface_velocity,
        surface_velocity=shape.surface_velocity,
        lower=LayerProfile(
            z=shape.lower_z,
            velocity=shape.lower_velocity,
            shear=shape.lower_shear,
            pressure=pressure_lower,
        ),
        upper=LayerProfile(
            z=shape.upper_z,
            velocity=shape.upper_velocity,
            shear=shape.upper_shear,
            pressure=pressure_upper,
        ),
    )


def family_from_base_flow(base: BaseFlow) -> DimensionlessFamily:
    """Extract the fixed dimensionless family represented by a BaseFlow."""

    return DimensionlessFamily(
        depth_ratio=base.depth_ratio,
        density_ratio=base.density_ratio,
        n_lower=base.parameters.lower.flow_index,
        n_upper=base.parameters.upper.flow_index,
        consistency_ratio=(
            base.consistency_number_upper / base.consistency_number_lower
        ),
        epsilon_lower=base.epsilon_lower,
        epsilon_upper=base.epsilon_upper,
    )


# ---------------------------------------------------------------------------
# Branch evaluation
# ---------------------------------------------------------------------------


def wavenumber_from_specification(
    specification: WaveSpecification,
    shape: NormalFlowShape,
    froude: float,
    reynolds: float,
) -> float:
    """Return k_H = k_dimensional H_l for the requested wave convention."""

    specification.validate()
    if specification.kind == "fixed_kH":
        return specification.value
    if specification.kind == "fixed_lambda_over_H":
        return 2.0 * math.pi / specification.value

    slope = compatibility_slope(shape, froude, reynolds)
    return 2.0 * math.pi * slope / specification.value


def wavelength_measures(
    kH: float,
    slope: float,
) -> tuple[float, float]:
    """Return (lambda/H_l, S0*lambda/H_l)."""

    wavelength_over_H = 2.0 * math.pi / kH
    return wavelength_over_H, slope * wavelength_over_H


def _pair_at_target(
    base: BaseFlow,
    wavenumber: float,
    solver: SolverOptions,
    continuation_steps: int = 12,
) -> dict[str, Eigenmode]:
    """Identify both hydrodynamic modes at one target wavenumber."""

    target = float(wavenumber)
    if target <= 0.0 or not np.isfinite(target):
        raise ValueError("wavenumber must be positive and finite")

    seed = min(solver.seed_wavenumber, target / 5.0)
    local_solver = replace(solver, seed_wavenumber=max(seed, target * 1.0e-4))

    # In the asymptotic long-wave range the two physical roots are the two
    # downstream roots with the smallest |sigma|, so an isolated solve is
    # both reliable and much faster than constructing a hidden k path.
    if target <= max(0.02, 2.0 * local_solver.seed_wavenumber):
        modes, _, _ = solve_spectrum(base, target, local_solver)
        return select_long_wave_pair(modes)

    ks = np.geomspace(local_solver.seed_wavenumber, target, continuation_steps)
    tracks, _ = track_hydrodynamic_modes(base, ks, local_solver)
    return {name: modes[-1] for name, modes in tracks.items()}


def mode_at_state(
    shape: NormalFlowShape,
    froude: float,
    reynolds: float,
    wavenumber: float,
    branch: BranchName,
    solver: SolverOptions,
    continuation_steps: int = 12,
) -> Eigenmode:
    base = build_dimensionless_base_flow(shape, froude, reynolds)
    return _pair_at_target(base, wavenumber, solver, continuation_steps)[branch]


def reduced_long_wave_growth(
    shape: NormalFlowShape,
    froude: float,
    reynolds: float,
    wavenumber: float,
    branch: BranchName,
    solver: SolverOptions,
) -> tuple[float, Eigenmode]:
    """Return Im(c)/(k Re), whose k->0 limit controls long-wave onset."""

    mode = mode_at_state(
        shape,
        froude,
        reynolds,
        wavenumber,
        branch,
        solver,
        continuation_steps=6,
    )
    return mode.phase_speed.imag / (wavenumber * reynolds), mode


def extrapolate_long_wave_coefficient(
    shape: NormalFlowShape,
    froude: float,
    reynolds: float,
    branch: BranchName,
    solver: SolverOptions,
    wavenumbers: Sequence[float],
    degree: int = 2,
) -> dict[str, object]:
    """Extrapolate the reduced k^2 growth coefficient to k=0.

    For a hydrodynamic branch,

        sigma = -i c0 k + Re_l D(Fr_l) k^2 + O(k^3),

    so ``Im(c)/(k Re_l) -> D``.  Reality symmetry makes the leading
    extrapolation error even in k; the fit therefore uses powers of k^2.
    """

    ks = np.asarray(wavenumbers, dtype=float)
    if ks.ndim != 1 or ks.size < degree + 1:
        raise ValueError("Insufficient wavenumbers for the requested fit")

    values: list[float] = []
    celerities: list[float] = []
    residuals: list[float] = []
    for k in ks:
        value, mode = reduced_long_wave_growth(
            shape,
            froude,
            reynolds,
            float(k),
            branch,
            solver,
        )
        values.append(value)
        celerities.append(mode.celerity)
        residuals.append(mode.residual)

    x = ks * ks
    coeff = np.polynomial.polynomial.polyfit(x, np.asarray(values), degree)
    fitted = np.polynomial.polynomial.polyval(x, coeff)
    c_coeff = np.polynomial.polynomial.polyfit(
        x,
        np.asarray(celerities),
        min(degree, len(ks) - 1),
    )

    return {
        "reduced_growth": float(coeff[0]),
        "long_wave_celerity": float(c_coeff[0]),
        "coefficients_in_k2": tuple(float(v) for v in coeff),
        "sample_values": tuple(float(v) for v in values),
        "wavenumbers": tuple(float(v) for v in ks),
        "maximum_fit_residual": float(np.max(np.abs(fitted - values))),
        "maximum_eigen_residual": float(max(residuals)),
    }


# ---------------------------------------------------------------------------
# Infinite-wavelength threshold
# ---------------------------------------------------------------------------


def _all_sign_change_brackets(
    x: FloatArray,
    values: FloatArray,
) -> list[tuple[float, float]]:
    brackets: list[tuple[float, float]] = []
    for index in range(len(x) - 1):
        a, b = float(x[index]), float(x[index + 1])
        fa, fb = float(values[index]), float(values[index + 1])
        if not np.isfinite(fa) or not np.isfinite(fb):
            continue
        if fa == 0.0:
            brackets.append((a, a))
        elif fa * fb < 0.0:
            brackets.append((a, b))
    if len(x) and values[-1] == 0.0:
        brackets.append((float(x[-1]), float(x[-1])))
    return brackets


def neutral_froude_roots_at_k(
    shape: NormalFlowShape,
    branch: BranchName,
    reynolds: float,
    wavenumber: float,
    solver: SolverOptions,
    froude_min: float,
    froude_max: float,
    scan_points: int,
    root_xtol: float = 2.0e-7,
) -> list[float]:
    """Return all branch-specific roots Im(c)=0 in a Froude interval."""

    fr_grid = np.linspace(froude_min, froude_max, scan_points)
    indicators = np.empty_like(fr_grid)

    for i, fr in enumerate(fr_grid):
        mode = mode_at_state(
            shape,
            float(fr),
            reynolds,
            wavenumber,
            branch,
            solver,
            continuation_steps=6 if wavenumber < 0.02 else 10,
        )
        indicators[i] = mode.phase_speed.imag

    brackets = _all_sign_change_brackets(fr_grid, indicators)
    roots: list[float] = []

    def indicator(fr: float) -> float:
        return mode_at_state(
            shape,
            fr,
            reynolds,
            wavenumber,
            branch,
            solver,
            continuation_steps=6 if wavenumber < 0.02 else 10,
        ).phase_speed.imag

    for a, b in brackets:
        if a == b:
            root = a
        else:
            root = float(brentq(indicator, a, b, xtol=root_xtol, rtol=2.0e-10))
        if not roots or abs(root - roots[-1]) > 10.0 * root_xtol:
            roots.append(root)
    return roots


def critical_froude_infinite_wavelength(
    shape: NormalFlowShape,
    branch: BranchName,
    solver: SolverOptions,
    options: LongWaveOptions = LongWaveOptions(),
) -> list[LongWaveThreshold]:
    """Compute all positive infinite-wavelength critical Froude numbers.

    A neutral Froude root is first calculated at each small finite k.  Each
    root family is then extrapolated as a polynomial in k^2 to k=0.  This
    avoids treating the singular descriptor pencil at exactly k=0 as an
    ordinary generalized eigenproblem.
    """

    options.validate()
    roots_by_k: list[list[float]] = []
    for k in options.wavenumbers:
        roots = neutral_froude_roots_at_k(
            shape,
            branch,
            options.reference_reynolds,
            k,
            solver,
            options.froude_min,
            options.froude_max,
            options.scan_points,
        )
        roots_by_k.append(roots)

    if not roots_by_k or any(len(roots) == 0 for roots in roots_by_k):
        return []

    # In ordinary use, the roots remain ordered and the same number is found
    # at every small k.  Restrict to the common count and state that root order
    # is a continuation label.
    common_count = min(len(roots) for roots in roots_by_k)
    ks = np.asarray(options.wavenumbers, dtype=float)
    x = ks * ks
    thresholds: list[LongWaveThreshold] = []

    for root_number in range(common_count):
        finite_roots = np.asarray(
            [roots[root_number] for roots in roots_by_k],
            dtype=float,
        )
        degree = min(options.polynomial_degree, len(ks) - 1)
        coeff = np.polynomial.polynomial.polyfit(x, finite_roots, degree)
        fitted = np.polynomial.polynomial.polyval(x, coeff)
        critical = float(coeff[0])

        longwave = extrapolate_long_wave_coefficient(
            shape,
            critical,
            options.reference_reynolds,
            branch,
            solver,
            options.wavenumbers,
            degree=degree,
        )
        fit_residual = float(np.max(np.abs(fitted - finite_roots)))
        root_spread = float(np.max(np.abs(finite_roots - critical)))
        allowed_residual = max(5.0e-4, 0.02 * abs(critical))
        converged = bool(
            np.isfinite(critical)
            and critical > 0.0
            and fit_residual <= allowed_residual
            and root_spread <= max(5.0e-3, 0.10 * abs(critical))
        )
        if converged:
            quality_message = (
                "finite-k neutral roots form a consistent k^2 extrapolation"
            )
        else:
            quality_message = (
                "finite-k neutral roots do not converge consistently to a "
                "positive k=0 intercept; do not interpret this as a physical "
                "critical Froude number"
            )

        thresholds.append(
            LongWaveThreshold(
                branch=branch,
                critical_froude=critical,
                finite_k_roots=tuple(float(v) for v in finite_roots),
                wavenumbers=tuple(float(v) for v in ks),
                fit_coefficients_in_k2=tuple(float(v) for v in coeff),
                maximum_fit_residual=fit_residual,
                reference_reynolds=options.reference_reynolds,
                compatibility_ratio=shape.compatibility_ratio,
                reduced_growth_at_root=float(longwave["reduced_growth"]),
                long_wave_celerity=float(longwave["long_wave_celerity"]),
                converged=converged,
                quality_message=quality_message,
            )
        )

    return thresholds


# ---------------------------------------------------------------------------
# Finite-wavelength neutral curves
# ---------------------------------------------------------------------------


def _scan_pair_over_froude(
    shape: NormalFlowShape,
    reynolds: float,
    specification: WaveSpecification,
    solver: SolverOptions,
    curve_options: NeutralCurveOptions,
) -> tuple[FloatArray, dict[str, list[Eigenmode]], FloatArray]:
    """Track both hydrodynamic modes as Fr varies at fixed Re."""

    specification.validate()
    curve_options.validate()
    fr_grid = np.linspace(
        curve_options.froude_min,
        curve_options.froude_max,
        curve_options.scan_points,
    )
    tracks: dict[str, list[Eigenmode]] = {
        "interface_varicose": [],
        "surface_zigzag": [],
    }
    k_values = np.empty_like(fr_grid)
    previous: dict[str, Eigenmode] | None = None

    for index, fr in enumerate(fr_grid):
        base = build_dimensionless_base_flow(shape, float(fr), reynolds)
        k = wavenumber_from_specification(
            specification,
            shape,
            float(fr),
            reynolds,
        )
        k_values[index] = k

        if previous is None:
            selected = _pair_at_target(
                base,
                k,
                solver,
                curve_options.continuation_steps,
            )
        else:
            local_solver = replace(
                solver,
                seed_wavenumber=min(solver.seed_wavenumber, k / 5.0),
            )
            modes, _, _ = solve_spectrum(base, k, local_solver)
            try:
                selected = continue_pair(previous, modes)
            except RuntimeError:
                selected = _pair_at_target(
                    base,
                    k,
                    solver,
                    curve_options.continuation_steps,
                )
        previous = selected
        for branch in tracks:
            tracks[branch].append(selected[branch])

    return fr_grid, tracks, k_values


def _neutral_point_from_root(
    shape: NormalFlowShape,
    branch: BranchName,
    reynolds: float,
    specification: WaveSpecification,
    solver: SolverOptions,
    curve_options: NeutralCurveOptions,
    root: float,
    root_number: int,
    mode: Eigenmode | None = None,
) -> NeutralPoint:
    slope = compatibility_slope(shape, root, reynolds)
    k = wavenumber_from_specification(specification, shape, root, reynolds)
    if mode is None:
        mode = mode_at_state(
            shape,
            root,
            reynolds,
            k,
            branch,
            solver,
            curve_options.continuation_steps,
        )
    lambda_over_H, slope_scaled = wavelength_measures(k, slope)
    return NeutralPoint(
        branch=branch,
        reynolds=reynolds,
        froude=root,
        slope=slope,
        wavenumber=k,
        wavelength_over_H=lambda_over_H,
        slope_scaled_wavelength=slope_scaled,
        growth_rate=mode.growth_rate,
        celerity=mode.celerity,
        residual=mode.residual,
        root_number=root_number,
    )


def finite_wavelength_neutral_points_at_reynolds(
    shape: NormalFlowShape,
    reynolds: float,
    specification: WaveSpecification,
    branch: BranchName,
    solver: SolverOptions,
    options: NeutralCurveOptions = NeutralCurveOptions(),
) -> list[NeutralPoint]:
    """Find every neutral Froude root at one generalized Reynolds number.

    The coarse scan supplies continuation-labelled endpoint modes.  During
    root refinement, each trial spectrum is matched to the nearer endpoint
    pair.  This prevents a generic fresh classification from jumping to a
    different branch inside a finite-k neutral bracket.
    """

    fr_grid, tracks, _ = _scan_pair_over_froude(
        shape,
        reynolds,
        specification,
        solver,
        options,
    )
    indicators = np.asarray(
        [mode.phase_speed.imag for mode in tracks[branch]],
        dtype=float,
    )

    points: list[NeutralPoint] = []
    for index in range(len(fr_grid) - 1):
        a = float(fr_grid[index])
        b = float(fr_grid[index + 1])
        fa = float(indicators[index])
        fb = float(indicators[index + 1])
        if not np.isfinite(fa) or not np.isfinite(fb):
            continue
        if fa != 0.0 and fa * fb >= 0.0:
            continue

        left_pair = {name: tracks[name][index] for name in tracks}
        right_pair = {name: tracks[name][index + 1] for name in tracks}
        pair_cache: dict[float, dict[str, Eigenmode]] = {
            round(a, 13): left_pair,
            round(b, 13): right_pair,
        }

        def pair_at(fr: float) -> dict[str, Eigenmode]:
            key = round(float(fr), 13)
            if key in pair_cache:
                return pair_cache[key]
            base = build_dimensionless_base_flow(shape, float(fr), reynolds)
            k = wavenumber_from_specification(
                specification, shape, float(fr), reynolds
            )
            local_solver = replace(
                solver,
                seed_wavenumber=min(solver.seed_wavenumber, k / 5.0),
            )
            modes, _, _ = solve_spectrum(base, k, local_solver)
            reference = left_pair if abs(fr - a) <= abs(fr - b) else right_pair
            try:
                selected = continue_pair(reference, modes)
            except RuntimeError:
                selected = _pair_at_target(
                    base, k, solver, options.continuation_steps
                )
            pair_cache[key] = selected
            return selected

        def indicator(fr: float) -> float:
            return pair_at(fr)[branch].phase_speed.imag

        if fa == 0.0:
            root = a
        else:
            try:
                root = float(
                    brentq(
                        indicator,
                        a,
                        b,
                        xtol=options.root_xtol,
                        rtol=2.0e-10,
                    )
                )
            except ValueError:
                # A rare mode-crossing ambiguity can make the locally matched
                # indicator discontinuous enough to defeat Brent's endpoint
                # check.  The continuation-labelled scan still brackets the
                # neutral point, so use its secant interpolation and retain a
                # diagnostic growth value in the output.
                root = a - fa * (b - a) / (fb - fa)

        if points and abs(root - points[-1].froude) <= 10.0 * options.root_xtol:
            continue
        root_mode = pair_at(root)[branch]
        points.append(
            _neutral_point_from_root(
                shape,
                branch,
                reynolds,
                specification,
                solver,
                options,
                root,
                len(points),
                mode=root_mode,
            )
        )

    return points


def finite_wavelength_neutral_curve(
    shape: NormalFlowShape,
    reynolds_values: Iterable[float],
    specification: WaveSpecification,
    branch: BranchName,
    solver: SolverOptions,
    options: NeutralCurveOptions = NeutralCurveOptions(),
) -> list[NeutralPoint]:
    """Return branch-specific neutral points in Fr-Re and Fr-S0 coordinates."""

    points: list[NeutralPoint] = []
    for reynolds in reynolds_values:
        if not np.isfinite(reynolds) or reynolds <= 0.0:
            raise ValueError("All Reynolds numbers must be positive and finite")
        points.extend(
            finite_wavelength_neutral_points_at_reynolds(
                shape,
                float(reynolds),
                specification,
                branch,
                solver,
                options,
            )
        )
    return points


# ---------------------------------------------------------------------------
# Tabular output
# ---------------------------------------------------------------------------


def write_long_wave_thresholds(
    thresholds: Sequence[LongWaveThreshold],
    path: str,
) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            "# branch\troot_number\tcritical_froude\treference_Re\tchi_l\t"
            "reduced_growth_at_root\tlong_wave_celerity\tfit_residual\t"
            "converged\tquality_message\twavenumbers\tfinite_k_roots\n"
        )
        for index, threshold in enumerate(thresholds):
            handle.write(
                f"{threshold.branch}\t{index}\t{threshold.critical_froude:.14e}\t"
                f"{threshold.reference_reynolds:.14e}\t"
                f"{threshold.compatibility_ratio:.14e}\t"
                f"{threshold.reduced_growth_at_root:.14e}\t"
                f"{threshold.long_wave_celerity:.14e}\t"
                f"{threshold.maximum_fit_residual:.14e}\t"
                f"{int(threshold.converged)}\t{threshold.quality_message}\t"
                + ",".join(f"{v:.14e}" for v in threshold.wavenumbers)
                + "\t"
                + ",".join(f"{v:.14e}" for v in threshold.finite_k_roots)
                + "\n"
            )


def write_neutral_curve(points: Sequence[NeutralPoint], path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            "# branch\troot_number\tRe_lower\tFr_lower\tS0\tkH\t"
            "lambda_over_H\tS0_lambda_over_H\tgrowth\tcelerity\tresidual\n"
        )
        for point in points:
            handle.write(
                f"{point.branch}\t{point.root_number}\t"
                f"{point.reynolds:.14e}\t{point.froude:.14e}\t"
                f"{point.slope:.14e}\t{point.wavenumber:.14e}\t"
                f"{point.wavelength_over_H:.14e}\t"
                f"{point.slope_scaled_wavelength:.14e}\t"
                f"{point.growth_rate:.14e}\t{point.celerity:.14e}\t"
                f"{point.residual:.14e}\n"
            )
