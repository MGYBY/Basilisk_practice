#!/usr/bin/env python3
"""Generate a dimensionless steady base state for the Basilisk case.

The input is an INI file containing the independent dimensionless groups used
in the audited full-Navier--Stokes report.  The lower consistency coefficient
Lambda_l is *derived* from the steady-uniform compatibility condition

    integral_0^1 U_l(z) dz = 1,

and Lambda_u = kappa_K Lambda_l.  The script writes ordinary text profiles and
a compact generated C header consumed by ``roll_wave_amr_dimensionless.c``.

No dimensional quantity enters the numerical PDE.  The optional
``[dimensional_reference]`` section is metadata used only to print conversion
scales and dimensional equivalents of selected outputs.
"""

from __future__ import annotations

import argparse
import configparser
import hashlib
import math
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np


MODEL_CODES = {
    "smooth_power_law": 0,
    "bounded_power_law": 1,
    "papanastasiou_hb": 2,
    "bounded_hb": 3,
}


def trapezoid(values: np.ndarray, coordinates: np.ndarray) -> float:
    """NumPy 1.x/2.x compatible trapezoidal integration."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, coordinates))
    return float(np.trapz(values, coordinates))  # pragma: no cover - NumPy < 2


def cumulative_trapezoid(values: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    result = np.zeros_like(values, dtype=float)
    if len(values) > 1:
        result[1:] = np.cumsum(
            0.5 * (values[:-1] + values[1:]) * np.diff(coordinates)
        )
    return result


@dataclass(frozen=True)
class Rheology:
    model: str
    n: float
    epsilon: float
    tau0: float
    m: float
    eta_min: float
    eta_max: float
    gamma_floor: float

    @property
    def model_code(self) -> int:
        return MODEL_CODES[self.model]


@dataclass(frozen=True)
class DimensionalReference:
    enabled: bool
    gravity: float
    lower_depth: float
    lower_density: float


@dataclass(frozen=True)
class Case:
    config_path: Path
    froude: float
    slope: float
    depth_ratio: float
    density_ratio: float
    consistency_ratio: float
    air_density_ratio: float
    air_eta: float
    air_streamwise_gravity: bool
    sigma_internal: float
    sigma_free: float
    lower: Rheology
    upper: Rheology
    lx: float
    nominal_ceiling: float
    max_level: int
    min_level: int
    initial_level: int
    ceiling_fine_fraction: float
    ceiling_level_drop: int
    initial_adapt_passes: int
    interface_padding: int
    interface_epsilon: float
    end_time: float
    max_dt: float
    cfl: float
    solver_tolerance: float
    solver_nitermax: int
    perturbation_mode: int
    tracked_mode: int
    perturbation_phase: float
    velocity_amplitude: float
    lower_depth_amplitude: float
    upper_depth_amplitude: float
    pressure_mode: int
    initial_audit_warning_tolerance: float
    u_error_lower: float
    v_error_lower: float
    u_error_air: float
    v_error_air: float
    vorticity_error: float
    shear_error: float
    output_dt: float
    gfs_dt: float
    dump_dt: float
    enable_dump_output: bool
    wave_amplitude_every: int
    lower_samples: int
    upper_samples: int
    air_samples: int
    integration_points: int
    reference: DimensionalReference

    @property
    def normal_gravity(self) -> float:
        return 1.0 / self.froude**2

    @property
    def streamwise_gravity(self) -> float:
        return self.slope / self.froude**2

    @property
    def liquid_depth(self) -> float:
        return 1.0 + self.depth_ratio

    @property
    def delta_min(self) -> float:
        return self.lx / (1 << self.max_level)

    @property
    def ceiling(self) -> float:
        return (
            math.ceil(self.nominal_ceiling / self.delta_min - 1.0e-12)
            * self.delta_min
        )

    @property
    def air_depth(self) -> float:
        return self.ceiling - self.liquid_depth

    @property
    def dimensional_velocity(self) -> float:
        if not self.reference.enabled:
            return math.nan
        cos_theta = 1.0 / math.sqrt(1.0 + self.slope**2)
        return self.froude * math.sqrt(
            self.reference.gravity * self.reference.lower_depth * cos_theta
        )

    @property
    def dimensional_time_scale(self) -> float:
        if not self.reference.enabled:
            return math.nan
        return self.reference.lower_depth / self.dimensional_velocity

    @property
    def dimensional_pressure_scale(self) -> float:
        if not self.reference.enabled:
            return math.nan
        return self.reference.lower_density * self.dimensional_velocity**2

    @property
    def dimensional_viscosity_scale(self) -> float:
        if not self.reference.enabled:
            return math.nan
        return (
            self.reference.lower_density
            * self.dimensional_velocity
            * self.reference.lower_depth
        )

    @property
    def dimensional_surface_tension_scale(self) -> float:
        if not self.reference.enabled:
            return math.nan
        return self.dimensional_pressure_scale * self.reference.lower_depth


def _required_float(parser: configparser.ConfigParser, section: str, name: str) -> float:
    try:
        value = parser.getfloat(section, name)
    except (configparser.Error, ValueError) as exc:
        raise ValueError(f"Missing or invalid [{section}] {name}") from exc
    if not math.isfinite(value):
        raise ValueError(f"[{section}] {name} must be finite")
    return value


def _required_int(parser: configparser.ConfigParser, section: str, name: str) -> int:
    try:
        return parser.getint(section, name)
    except (configparser.Error, ValueError) as exc:
        raise ValueError(f"Missing or invalid [{section}] {name}") from exc


def _required_bool(parser: configparser.ConfigParser, section: str, name: str) -> bool:
    try:
        return parser.getboolean(section, name)
    except (configparser.Error, ValueError) as exc:
        raise ValueError(f"Missing or invalid [{section}] {name}") from exc


def read_rheology(parser: configparser.ConfigParser, section: str) -> Rheology:
    model = parser.get(section, "model").strip().lower()
    if model not in MODEL_CODES:
        allowed = ", ".join(MODEL_CODES)
        raise ValueError(f"[{section}] model must be one of: {allowed}")
    return Rheology(
        model=model,
        n=_required_float(parser, section, "power_index"),
        epsilon=_required_float(parser, section, "epsilon"),
        tau0=_required_float(parser, section, "yield_stress_star"),
        m=_required_float(parser, section, "papanastasiou_m_star"),
        eta_min=_required_float(parser, section, "eta_min"),
        eta_max=_required_float(parser, section, "eta_max"),
        gamma_floor=_required_float(parser, section, "gamma_floor"),
    )


def load_case(path: Path) -> Case:
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not parser.read(path):
        raise FileNotFoundError(path)

    reference = DimensionalReference(
        enabled=_required_bool(parser, "dimensional_reference", "enabled"),
        gravity=_required_float(parser, "dimensional_reference", "gravity_m_s2"),
        lower_depth=_required_float(parser, "dimensional_reference", "lower_depth_m"),
        lower_density=_required_float(
            parser, "dimensional_reference", "lower_density_kg_m3"
        ),
    )

    case = Case(
        config_path=path.resolve(),
        froude=_required_float(parser, "physical_problem", "froude_lower"),
        slope=_required_float(parser, "physical_problem", "slope_tan"),
        depth_ratio=_required_float(parser, "physical_problem", "depth_ratio"),
        density_ratio=_required_float(parser, "physical_problem", "density_ratio"),
        consistency_ratio=_required_float(
            parser, "physical_problem", "consistency_ratio"
        ),
        air_density_ratio=_required_float(
            parser, "physical_problem", "air_density_ratio"
        ),
        air_eta=_required_float(
            parser, "physical_problem", "air_dynamic_viscosity"
        ),
        air_streamwise_gravity=_required_bool(
            parser, "physical_problem", "air_streamwise_gravity"
        ),
        sigma_internal=_required_float(
            parser, "physical_problem", "sigma_internal_star"
        ),
        sigma_free=_required_float(parser, "physical_problem", "sigma_free_star"),
        lower=read_rheology(parser, "lower_rheology"),
        upper=read_rheology(parser, "upper_rheology"),
        lx=_required_float(parser, "domain_and_mesh", "lx"),
        nominal_ceiling=_required_float(
            parser, "domain_and_mesh", "nominal_ceiling"
        ),
        max_level=_required_int(parser, "domain_and_mesh", "max_level"),
        min_level=_required_int(parser, "domain_and_mesh", "min_level"),
        initial_level=_required_int(parser, "domain_and_mesh", "initial_level"),
        ceiling_fine_fraction=_required_float(
            parser, "domain_and_mesh", "ceiling_fine_fraction"
        ),
        ceiling_level_drop=_required_int(
            parser, "domain_and_mesh", "ceiling_level_drop"
        ),
        initial_adapt_passes=_required_int(
            parser, "domain_and_mesh", "initial_adapt_passes"
        ),
        interface_padding=_required_int(
            parser, "domain_and_mesh", "interface_padding"
        ),
        interface_epsilon=_required_float(
            parser, "domain_and_mesh", "interface_epsilon"
        ),
        end_time=_required_float(parser, "time_integration", "end_time"),
        max_dt=_required_float(parser, "time_integration", "max_dt"),
        cfl=_required_float(parser, "time_integration", "cfl"),
        solver_tolerance=_required_float(
            parser, "time_integration", "solver_tolerance"
        ),
        solver_nitermax=_required_int(
            parser, "time_integration", "solver_nitermax"
        ),
        perturbation_mode=_required_int(
            parser, "initial_condition", "perturbation_mode"
        ),
        tracked_mode=_required_int(parser, "initial_condition", "tracked_mode"),
        perturbation_phase=_required_float(
            parser, "initial_condition", "perturbation_phase"
        ),
        velocity_amplitude=_required_float(
            parser, "initial_condition", "velocity_amplitude"
        ),
        lower_depth_amplitude=_required_float(
            parser, "initial_condition", "lower_depth_amplitude"
        ),
        upper_depth_amplitude=_required_float(
            parser, "initial_condition", "upper_depth_amplitude"
        ),
        pressure_mode=_required_int(parser, "initial_condition", "pressure_mode"),
        initial_audit_warning_tolerance=_required_float(
            parser, "initial_condition", "initial_audit_warning_tolerance"
        ),
        u_error_lower=_required_float(parser, "adaptation", "u_error_lower"),
        v_error_lower=_required_float(parser, "adaptation", "v_error_lower"),
        u_error_air=_required_float(parser, "adaptation", "u_error_air"),
        v_error_air=_required_float(parser, "adaptation", "v_error_air"),
        vorticity_error=_required_float(parser, "adaptation", "vorticity_error"),
        shear_error=_required_float(parser, "adaptation", "shear_error"),
        output_dt=_required_float(parser, "output", "output_dt"),
        gfs_dt=_required_float(parser, "output", "gfs_dt"),
        dump_dt=_required_float(parser, "output", "dump_dt"),
        enable_dump_output=_required_bool(
            parser, "output", "enable_dump_output"
        ),
        wave_amplitude_every=_required_int(
            parser, "output", "wave_amplitude_every"
        ),
        lower_samples=_required_int(parser, "profile_sampling", "lower_samples"),
        upper_samples=_required_int(parser, "profile_sampling", "upper_samples"),
        air_samples=_required_int(parser, "profile_sampling", "air_samples"),
        integration_points=_required_int(
            parser, "profile_sampling", "integration_points"
        ),
        reference=reference,
    )
    validate_case(case)
    return case


def validate_rheology(name: str, r: Rheology) -> None:
    if r.n <= 0.0:
        raise ValueError(f"{name}: power_index must be positive")
    if r.epsilon < 0.0 or r.tau0 < 0.0 or r.m < 0.0:
        raise ValueError(f"{name}: epsilon, yield stress and m must be nonnegative")
    if r.eta_min < 0.0 or r.eta_max <= 0.0 or r.eta_max < r.eta_min:
        raise ValueError(f"{name}: invalid viscosity bounds")
    if r.gamma_floor <= 0.0:
        raise ValueError(f"{name}: gamma_floor must be positive")
    if r.model == "smooth_power_law" and r.n < 1.0 and r.epsilon <= 0.0:
        raise ValueError(f"{name}: shear-thinning smooth model requires epsilon > 0")
    if r.model == "papanastasiou_hb" and r.tau0 > 0.0 and r.m <= 0.0:
        raise ValueError(f"{name}: Papanastasiou HB with yield stress requires m > 0")


def validate_case(c: Case) -> None:
    positive = {
        "froude_lower": c.froude,
        "slope_tan": c.slope,
        "depth_ratio": c.depth_ratio,
        "density_ratio": c.density_ratio,
        "consistency_ratio": c.consistency_ratio,
        "air_density_ratio": c.air_density_ratio,
        "air_dynamic_viscosity": c.air_eta,
        "lx": c.lx,
        "nominal_ceiling": c.nominal_ceiling,
        "end_time": c.end_time,
        "max_dt": c.max_dt,
        "cfl": c.cfl,
        "solver_tolerance": c.solver_tolerance,
        "output_dt": c.output_dt,
        "gfs_dt": c.gfs_dt,
        "dump_dt": c.dump_dt,
        "interface_epsilon": c.interface_epsilon,
    }
    invalid = [k for k, v in positive.items() if not math.isfinite(v) or v <= 0.0]
    if invalid:
        raise ValueError("These values must be finite and positive: " + ", ".join(invalid))
    if c.sigma_internal < 0.0 or c.sigma_free < 0.0:
        raise ValueError("Surface-tension coefficients cannot be negative")
    if not (1 <= c.min_level <= c.initial_level <= c.max_level):
        raise ValueError("Require 1 <= min_level <= initial_level <= max_level")
    if c.nominal_ceiling <= c.liquid_depth:
        raise ValueError("The ceiling must lie above both liquid layers")
    if c.ceiling_fine_fraction <= 0.0 or c.ceiling_fine_fraction > 1.0:
        raise ValueError("ceiling_fine_fraction must be in (0,1]")
    if c.ceiling <= c.liquid_depth:
        raise ValueError("The aligned ceiling lies inside the liquid")
    if c.initial_adapt_passes < 0 or c.interface_padding < 0:
        raise ValueError("AMR pass and padding counts cannot be negative")
    if c.solver_nitermax <= 0 or c.wave_amplitude_every <= 0:
        raise ValueError("Iteration/output counts must be positive")
    if c.perturbation_mode <= 0 or c.tracked_mode <= 0:
        raise ValueError("Perturbation and tracked mode numbers must be positive")
    if abs(c.lower_depth_amplitude) >= 1.0 or abs(c.upper_depth_amplitude) >= 1.0:
        raise ValueError("Depth amplitudes must have absolute value < 1")
    if abs(c.velocity_amplitude) >= 1.0:
        raise ValueError("Velocity amplitude must have absolute value < 1")
    if c.pressure_mode not in (0, 1):
        raise ValueError("pressure_mode must be 0 or 1")
    if min(c.lower_samples, c.upper_samples, c.air_samples) < 2:
        raise ValueError("Each profile region needs at least two samples")
    if c.integration_points < 501:
        raise ValueError("integration_points must be at least 501")
    validate_rheology("lower_rheology", c.lower)
    validate_rheology("upper_rheology", c.upper)
    if c.reference.enabled:
        if (
            c.reference.gravity <= 0.0
            or c.reference.lower_depth <= 0.0
            or c.reference.lower_density <= 0.0
        ):
            raise ValueError("Enabled dimensional reference values must be positive")


def apparent_viscosity(gamma: np.ndarray, lam: float, r: Rheology) -> np.ndarray:
    gabs = np.abs(np.asarray(gamma, dtype=float))
    gsafe = np.maximum(gabs, r.gamma_floor)

    if r.model == "smooth_power_law":
        eta = lam * (gabs * gabs + r.epsilon * r.epsilon) ** (
            0.5 * (r.n - 1.0)
        )
    else:
        eta = lam * gsafe ** (r.n - 1.0)

    if r.model == "papanastasiou_hb" and r.tau0 > 0.0:
        x = r.m * gabs
        ratio = np.empty_like(gabs)
        small = np.abs(x) < 1.0e-5
        ratio[small] = r.m * (
            1.0 - x[small] / 2.0 + x[small] ** 2 / 6.0 - x[small] ** 3 / 24.0
        )
        ratio[~small] = -np.expm1(-x[~small]) / np.maximum(
            gabs[~small], r.gamma_floor
        )
        eta = eta + r.tau0 * ratio
    elif r.model == "bounded_hb" and r.tau0 > 0.0:
        eta = eta + r.tau0 / gsafe

    return np.clip(eta, r.eta_min, r.eta_max)


def stress_from_gamma(gamma: np.ndarray, lam: float, r: Rheology) -> np.ndarray:
    gamma = np.asarray(gamma, dtype=float)
    return apparent_viscosity(gamma, lam, r) * gamma


def invert_stress(tau: np.ndarray, lam: float, r: Rheology) -> np.ndarray:
    """Invert the monotone positive-shear constitutive relation by bisection."""
    tau = np.maximum(np.asarray(tau, dtype=float), 0.0)
    low = np.zeros_like(tau)
    # Large enough for either a power-law or plateau-dominated response.
    estimate_power = (tau / max(lam, 1.0e-300) + 1.0e-300) ** (1.0 / r.n)
    estimate_plateau = tau / max(r.eta_max, 1.0e-300)
    high = np.maximum(1.0, 2.0 * np.maximum(estimate_power, estimate_plateau))

    for _ in range(80):
        insufficient = stress_from_gamma(high, lam, r) < tau
        if not np.any(insufficient):
            break
        high[insufficient] *= 2.0
    else:
        raise RuntimeError("Could not bracket constitutive inversion")

    for _ in range(90):
        mid = 0.5 * (low + high)
        left = stress_from_gamma(mid, lam, r) < tau
        low[left] = mid[left]
        high[~left] = mid[~left]
    return 0.5 * (low + high)


def traction_fields(c: Case, z_lower: np.ndarray, y_upper: np.ndarray, y_air: np.ndarray):
    air_weight = (
        c.air_density_ratio * c.air_depth if c.air_streamwise_gravity else 0.0
    )
    g = c.streamwise_gravity
    tau_lower = g * (
        1.0 - z_lower + c.density_ratio * c.depth_ratio + air_weight
    )
    tau_upper = g * (
        c.density_ratio * (c.depth_ratio - y_upper) + air_weight
    )
    if c.air_streamwise_gravity:
        tau_air = g * c.air_density_ratio * (c.air_depth - y_air)
    else:
        tau_air = np.zeros_like(y_air)
    return tau_lower, tau_upper, tau_air


def lower_mean_velocity(c: Case, lam_lower: float) -> float:
    z = np.linspace(0.0, 1.0, c.integration_points)
    tau, _, _ = traction_fields(c, z, np.asarray([0.0]), np.asarray([0.0]))
    gamma = invert_stress(tau, lam_lower, c.lower)
    return trapezoid((1.0 - z) * gamma, z)


def derive_lambda_lower(c: Case) -> float:
    """Solve the steady-uniform compatibility condition for Lambda_l."""
    low, high = 1.0e-14, 1.0
    mean_low = lower_mean_velocity(c, low)
    while mean_low <= 1.0:
        low *= 0.1
        mean_low = lower_mean_velocity(c, low)
        if low < 1.0e-300:
            raise RuntimeError("Could not find low-Lambda compatibility bracket")

    mean_high = lower_mean_velocity(c, high)
    while mean_high >= 1.0:
        high *= 10.0
        mean_high = lower_mean_velocity(c, high)
        if high > 1.0e300:
            raise RuntimeError("Could not find high-Lambda compatibility bracket")

    for _ in range(90):
        mid = math.sqrt(low * high)
        if lower_mean_velocity(c, mid) > 1.0:
            low = mid
        else:
            high = mid
    return math.sqrt(low * high)


def solve_base_state(c: Case):
    lam_lower = derive_lambda_lower(c)
    lam_upper = c.consistency_ratio * lam_lower

    z1 = np.linspace(0.0, 1.0, c.integration_points)
    y2 = np.linspace(0.0, c.depth_ratio, c.integration_points)
    y3 = np.linspace(0.0, c.air_depth, c.integration_points)
    tau1, tau2, tau3 = traction_fields(c, z1, y2, y3)
    gamma1 = invert_stress(tau1, lam_lower, c.lower)
    gamma2 = invert_stress(tau2, lam_upper, c.upper)
    gamma3 = tau3 / c.air_eta

    u1 = cumulative_trapezoid(gamma1, z1)
    interface_velocity = float(u1[-1])
    u2 = interface_velocity + cumulative_trapezoid(gamma2, y2)
    surface_velocity = float(u2[-1])
    u3 = surface_velocity + cumulative_trapezoid(gamma3, y3)

    mean1 = trapezoid(u1, z1)
    mean2 = trapezoid(u2, y2) / c.depth_ratio
    mean3 = trapezoid(u3, y3) / c.air_depth

    return {
        "lambda_lower": lam_lower,
        "lambda_upper": lam_upper,
        "z1": z1,
        "y2": y2,
        "y3": y3,
        "tau1": tau1,
        "tau2": tau2,
        "tau3": tau3,
        "gamma1": gamma1,
        "gamma2": gamma2,
        "gamma3": gamma3,
        "u1": u1,
        "u2": u2,
        "u3": u3,
        "interface_velocity": interface_velocity,
        "surface_velocity": surface_velocity,
        "ceiling_velocity": float(u3[-1]),
        "mean1": mean1,
        "mean2": mean2,
        "mean3": mean3,
    }


def profile_coordinates(c: Case) -> np.ndarray:
    lower = np.linspace(0.0, 1.0, c.lower_samples)
    upper = np.linspace(1.0, c.liquid_depth, c.upper_samples)[1:]
    air = np.linspace(c.liquid_depth, c.ceiling, c.air_samples)[1:]
    z = np.concatenate((lower, upper, air))
    if np.any(np.diff(z) <= 0.0):
        raise RuntimeError("Generated profile coordinates are not increasing")
    return z


def interpolate_base(z: np.ndarray, c: Case, base: dict):
    velocity = np.empty_like(z)
    gamma = np.empty_like(z)
    eta = np.empty_like(z)
    tau = np.empty_like(z)

    lower = z <= 1.0
    upper = (z > 1.0) & (z <= c.liquid_depth)
    air = z > c.liquid_depth

    velocity[lower] = np.interp(z[lower], base["z1"], base["u1"])
    gamma[lower] = np.interp(z[lower], base["z1"], base["gamma1"])
    tau[lower] = np.interp(z[lower], base["z1"], base["tau1"])
    eta[lower] = apparent_viscosity(
        gamma[lower], base["lambda_lower"], c.lower
    )

    local_upper = z[upper] - 1.0
    velocity[upper] = np.interp(local_upper, base["y2"], base["u2"])
    gamma[upper] = np.interp(local_upper, base["y2"], base["gamma2"])
    tau[upper] = np.interp(local_upper, base["y2"], base["tau2"])
    eta[upper] = apparent_viscosity(
        gamma[upper], base["lambda_upper"], c.upper
    )

    local_air = z[air] - c.liquid_depth
    velocity[air] = np.interp(local_air, base["y3"], base["u3"])
    gamma[air] = np.interp(local_air, base["y3"], base["gamma3"])
    tau[air] = np.interp(local_air, base["y3"], base["tau3"])
    eta[air] = c.air_eta

    return velocity, gamma, eta, tau


def pressure_profile(z: np.ndarray, c: Case) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    p = np.empty_like(z)
    lower = z <= 1.0
    upper = (z > 1.0) & (z <= c.liquid_depth)
    air = z > c.liquid_depth
    gn = c.normal_gravity
    p_air_surface = c.air_density_ratio * gn * c.air_depth
    p_upper_interface = p_air_surface + c.density_ratio * gn * c.depth_ratio
    p[air] = c.air_density_ratio * gn * (c.ceiling - z[air])
    p[upper] = p_air_surface + c.density_ratio * gn * (c.liquid_depth - z[upper])
    p[lower] = p_upper_interface + gn * (1.0 - z[lower])
    return p


def phase_fractions(z: np.ndarray, c: Case):
    f1 = (z <= 1.0).astype(float)
    f2 = ((z > 1.0) & (z <= c.liquid_depth)).astype(float)
    f3 = (z > c.liquid_depth).astype(float)
    return f1, f2, f3


def write_two_column(path: Path, z: np.ndarray, values: np.ndarray, name: str) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write("# Generated by generate_dimensionless_base_state.py\n")
        stream.write(f"# columns: z_star    {name}_star\n")
        for zz, value in zip(z, values):
            stream.write(f"{zz:.16e}\t{value:.16e}\n")


def c_bool(value: bool) -> str:
    return "1" if value else "0"


def c_float(value: float) -> str:
    if math.isnan(value):
        return "NAN"
    return f"{value:.17e}"


def write_generated_header(path: Path, c: Case, base: dict) -> None:
    h_ref = c.reference.lower_depth if c.reference.enabled else math.nan
    u_ref = c.dimensional_velocity
    rho_ref = c.reference.lower_density if c.reference.enabled else math.nan
    t_ref = c.dimensional_time_scale
    p_ref = c.dimensional_pressure_scale
    eta_ref = c.dimensional_viscosity_scale
    sigma_ref = c.dimensional_surface_tension_scale

    entries: list[tuple[str, str]] = [
        ("CASE_CONFIG_FILE", f'"{c.config_path.name}"'),
        ("CASE_FROUDE", c_float(c.froude)),
        ("CASE_SLOPE_TAN", c_float(c.slope)),
        ("CASE_DEPTH_RATIO", c_float(c.depth_ratio)),
        ("CASE_DENSITY_RATIO", c_float(c.density_ratio)),
        ("CASE_CONSISTENCY_RATIO", c_float(c.consistency_ratio)),
        ("CASE_AIR_DENSITY_RATIO", c_float(c.air_density_ratio)),
        ("CASE_AIR_ETA", c_float(c.air_eta)),
        ("CASE_AIR_STREAMWISE_GRAVITY", c_bool(c.air_streamwise_gravity)),
        ("CASE_SIGMA_INTERNAL", c_float(c.sigma_internal)),
        ("CASE_SIGMA_FREE", c_float(c.sigma_free)),
        ("CASE_LOWER_MODEL", str(c.lower.model_code)),
        ("CASE_LOWER_N", c_float(c.lower.n)),
        ("CASE_LOWER_EPSILON", c_float(c.lower.epsilon)),
        ("CASE_LOWER_TAU0", c_float(c.lower.tau0)),
        ("CASE_LOWER_M", c_float(c.lower.m)),
        ("CASE_LOWER_ETA_MIN", c_float(c.lower.eta_min)),
        ("CASE_LOWER_ETA_MAX", c_float(c.lower.eta_max)),
        ("CASE_LOWER_GAMMA_FLOOR", c_float(c.lower.gamma_floor)),
        ("CASE_UPPER_MODEL", str(c.upper.model_code)),
        ("CASE_UPPER_N", c_float(c.upper.n)),
        ("CASE_UPPER_EPSILON", c_float(c.upper.epsilon)),
        ("CASE_UPPER_TAU0", c_float(c.upper.tau0)),
        ("CASE_UPPER_M", c_float(c.upper.m)),
        ("CASE_UPPER_ETA_MIN", c_float(c.upper.eta_min)),
        ("CASE_UPPER_ETA_MAX", c_float(c.upper.eta_max)),
        ("CASE_UPPER_GAMMA_FLOOR", c_float(c.upper.gamma_floor)),
        ("CASE_LAMBDA_LOWER", c_float(base["lambda_lower"])),
        ("CASE_LAMBDA_UPPER", c_float(base["lambda_upper"])),
        ("CASE_LX", c_float(c.lx)),
        ("CASE_NOMINAL_CEILING", c_float(c.nominal_ceiling)),
        ("CASE_ALIGNED_CEILING", c_float(c.ceiling)),
        ("CASE_MAXLEVEL", str(c.max_level)),
        ("CASE_MINLEVEL", str(c.min_level)),
        ("CASE_INITLEVEL", str(c.initial_level)),
        ("CASE_CEILING_FINE_FRACTION", c_float(c.ceiling_fine_fraction)),
        ("CASE_CEILING_LEVEL_DROP", str(c.ceiling_level_drop)),
        ("CASE_INIT_ADAPT_PASSES", str(c.initial_adapt_passes)),
        ("CASE_INTERFACE_PADDING", str(c.interface_padding)),
        ("CASE_INTERFACE_EPS", c_float(c.interface_epsilon)),
        ("CASE_END_TIME", c_float(c.end_time)),
        ("CASE_MAX_DT", c_float(c.max_dt)),
        ("CASE_CFL", c_float(c.cfl)),
        ("CASE_SOLVER_TOLERANCE", c_float(c.solver_tolerance)),
        ("CASE_SOLVER_NITERMAX", str(c.solver_nitermax)),
        ("CASE_PERTURBATION_MODE", str(c.perturbation_mode)),
        ("CASE_TRACKED_MODE", str(c.tracked_mode)),
        ("CASE_PERTURBATION_PHASE", c_float(c.perturbation_phase)),
        ("CASE_VELOCITY_AMPLITUDE", c_float(c.velocity_amplitude)),
        ("CASE_LOWER_DEPTH_AMPLITUDE", c_float(c.lower_depth_amplitude)),
        ("CASE_UPPER_DEPTH_AMPLITUDE", c_float(c.upper_depth_amplitude)),
        ("CASE_PRESSURE_MODE", str(c.pressure_mode)),
        ("CASE_INITIAL_AUDIT_WARNING_TOL", c_float(c.initial_audit_warning_tolerance)),
        ("CASE_U_ERROR_LOWER", c_float(c.u_error_lower)),
        ("CASE_V_ERROR_LOWER", c_float(c.v_error_lower)),
        ("CASE_U_ERROR_AIR", c_float(c.u_error_air)),
        ("CASE_V_ERROR_AIR", c_float(c.v_error_air)),
        ("CASE_VORTICITY_ERROR", c_float(c.vorticity_error)),
        ("CASE_SHEAR_ERROR", c_float(c.shear_error)),
        ("CASE_OUTPUT_DT", c_float(c.output_dt)),
        ("CASE_GFS_DT", c_float(c.gfs_dt)),
        ("CASE_DUMP_DT", c_float(c.dump_dt)),
        ("CASE_ENABLE_DUMP_OUTPUT", c_bool(c.enable_dump_output)),
        ("CASE_WAVE_AMPLITUDE_EVERY", str(c.wave_amplitude_every)),
        ("CASE_U_INTERFACE", c_float(base["interface_velocity"])),
        ("CASE_U_SURFACE", c_float(base["surface_velocity"])),
        ("CASE_U_CEILING", c_float(base["ceiling_velocity"])),
        ("CASE_P_BED", c_float(float(pressure_profile(np.asarray([0.0]), c)[0]))),
        ("CASE_REFERENCE_ENABLED", c_bool(c.reference.enabled)),
        ("CASE_HREF_M", c_float(h_ref)),
        ("CASE_UREF_MS", c_float(u_ref)),
        ("CASE_RHOREF_KGM3", c_float(rho_ref)),
        ("CASE_TREF_S", c_float(t_ref)),
        ("CASE_PREF_PA", c_float(p_ref)),
        ("CASE_ETAREF_PAS", c_float(eta_ref)),
        ("CASE_SIGMAREF_NPM", c_float(sigma_ref)),
    ]

    with path.open("w", encoding="utf-8") as stream:
        stream.write("/* Generated file. Edit case_parameters.ini, not this header. */\n")
        stream.write("#ifndef BASILISK_DIMENSIONLESS_GENERATED_CASE_H\n")
        stream.write("#define BASILISK_DIMENSIONLESS_GENERATED_CASE_H\n\n")
        for key, value in entries:
            stream.write(f"#define {key} {value}\n")
        stream.write("\n#endif\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_outputs(output_dir: Path, c: Case) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    base = solve_base_state(c)
    z = profile_coordinates(c)
    velocity, gamma, eta, tau = interpolate_base(z, c, base)
    pressure = pressure_profile(z, c)
    f1, f2, f3 = phase_fractions(z, c)

    velocity_file = output_dir / "velocity.dat"
    pressure_file = output_dir / "pressure.dat"
    combined_file = output_dir / "combined_profile.tsv"
    summary_file = output_dir / "base_state_summary.txt"
    header_file = output_dir / "generated_case.h"
    config_copy = output_dir / "used_case_parameters.ini"
    manifest_file = output_dir / "MANIFEST.sha256"

    write_two_column(velocity_file, z, velocity, "u_x")
    write_two_column(pressure_file, z, pressure, "p")
    np.savetxt(
        combined_file,
        np.column_stack(
            (z, f1, f2, f3, velocity, np.zeros_like(z), pressure, gamma, eta, tau)
        ),
        delimiter="\t",
        header=(
            "z_star\tf1\tf2\tf3\tu_x_star\tu_z_star\tp_star\t"
            "gamma_star\teta_star\ttau_xz_star"
        ),
        comments="# ",
        fmt="%.16e",
    )
    write_generated_header(header_file, c, base)
    shutil.copyfile(c.config_path, config_copy)

    # Constitutive and compatibility checks on the dense integration grids.
    residual1 = np.max(
        np.abs(stress_from_gamma(base["gamma1"], base["lambda_lower"], c.lower) - base["tau1"])
    )
    residual2 = np.max(
        np.abs(stress_from_gamma(base["gamma2"], base["lambda_upper"], c.upper) - base["tau2"])
    )
    tau_b = float(base["tau1"][0])
    tau_i = float(base["tau1"][-1])
    tau_surface = float(base["tau2"][-1])
    we_internal = math.inf if c.sigma_internal == 0.0 else 1.0 / c.sigma_internal
    we_free = math.inf if c.sigma_free == 0.0 else 1.0 / c.sigma_free

    dimensional = "disabled"
    if c.reference.enabled:
        dimensional = f"""enabled
H_l                  = {c.reference.lower_depth:.12g} m
U_l                  = {c.dimensional_velocity:.12g} m/s
rho_l                = {c.reference.lower_density:.12g} kg/m^3
time scale           = {c.dimensional_time_scale:.12g} s
pressure scale       = {c.dimensional_pressure_scale:.12g} Pa
viscosity scale      = {c.dimensional_viscosity_scale:.12g} Pa s
surface-tension scale= {c.dimensional_surface_tension_scale:.12g} N/m
end time             = {c.end_time*c.dimensional_time_scale:.12g} s
"""

    summary = f"""Dimensionless analytical base state for Basilisk
==================================================
Configuration         = {c.config_path.name}
Length scale          = H_l
Velocity scale        = U_l (lower steady mean)
Time scale            = H_l/U_l
Pressure scale        = rho_l U_l^2
Dynamic-viscosity scale = rho_l U_l H_l

Independent liquid groups
-------------------------
Fr_l                  = {c.froude:.12g}
S_0=tan(theta)        = {c.slope:.12g}
h_r                   = {c.depth_ratio:.12g}
r_rho                 = {c.density_ratio:.12g}
n_l, n_u              = {c.lower.n:.12g}, {c.upper.n:.12g}
kappa_K               = {c.consistency_ratio:.12g}
epsilon_l, epsilon_u  = {c.lower.epsilon:.12g}, {c.upper.epsilon:.12g}

Explicit-air and capillary groups
---------------------------------
r_a                   = {c.air_density_ratio:.12g}
eta_a                 = {c.air_eta:.12g}
air streamwise gravity= {int(c.air_streamwise_gravity)}
sigma*_internal       = {c.sigma_internal:.12g}   (We={we_internal:.12g})
sigma*_free           = {c.sigma_free:.12g}   (We={we_free:.12g})

Derived compatibility quantities
--------------------------------
Lambda_l              = {base['lambda_lower']:.12g}
Lambda_u              = {base['lambda_upper']:.12g}
mean U_l              = {base['mean1']:.12g}
mean U_u              = {base['mean2']:.12g}
mean U_a              = {base['mean3']:.12g}
U_interface           = {base['interface_velocity']:.12g}
U_surface             = {base['surface_velocity']:.12g}
U_ceiling             = {base['ceiling_velocity']:.12g}
tau_b                 = {tau_b:.12g}
tau_I                 = {tau_i:.12g}
tau_liquid_surface    = {tau_surface:.12g}
p_bed                 = {pressure_profile(np.asarray([0.0]), c)[0]:.12g}
p_ceiling             = {pressure_profile(np.asarray([c.ceiling]), c)[0]:.12g}
max constitutive residual lower = {residual1:.6e}
max constitutive residual upper = {residual2:.6e}
mean-normalization error         = {base['mean1'] - 1.0:.6e}

Domain and resolution
---------------------
Lx                    = {c.lx:.12g}
nominal ceiling       = {c.nominal_ceiling:.12g}
aligned ceiling       = {c.ceiling:.12g}
air depth             = {c.air_depth:.12g}
Delta_min             = {c.delta_min:.12g}
cells through H_l     = {1.0/c.delta_min:.12g}
cells through H_u     = {c.depth_ratio/c.delta_min:.12g}
profile rows          = {len(z)}

Optional dimensional conversion
-------------------------------
{dimensional}
"""
    summary_file.write_text(summary, encoding="utf-8")

    generated: Iterable[Path] = (
        velocity_file,
        pressure_file,
        combined_file,
        summary_file,
        header_file,
        config_copy,
    )
    manifest_file.write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in generated),
        encoding="utf-8",
    )

    print(summary, end="")
    print(f"Wrote dimensionless base state to: {output_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parent / "case_parameters.ini",
        help="dimensionless INI file (default: case_parameters.ini)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "base_state",
        help="output directory (default: base_state)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case = load_case(args.config)
    write_outputs(args.output, case)


if __name__ == "__main__":
    main()
