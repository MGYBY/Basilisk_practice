#!/usr/bin/env python3
"""Generate a dimensionless steady base state for the Basilisk case.

The input is an INI file containing the independent dimensionless groups used
in the audited full-Navier--Stokes report.  The lower consistency coefficient
Lambda_l is *derived* from the steady-uniform compatibility condition

    integral_0^1 U_l(z) dz = 1,

The upper coefficient Lambda_u is derived from the prescribed steady-
reference interfacial apparent-viscosity ratio

    R_eta,I = eta_lower,I/eta_upper,I.

The former kappa_K = Lambda_u/Lambda_l is a generated diagnostic.  The ambient Newtonian viscosity is then
derived from the configured upper-to-air dynamic-viscosity ratio.  The upper
reference viscosity is evaluated from the same constitutive model at the
characteristic dimensionless shear rate gamma=1.  The script writes ordinary
text profiles and a compact generated C header consumed by
``roll_wave_amr_dimensionless.c``.

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

# U_l/H_l is the reference shear-rate scale, so its dimensionless value is one.
UPPER_VISCOSITY_REFERENCE_SHEAR_RATE = 1.0


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
class FrontRunnerSettings:
    enabled: bool = False
    wavelength_slope_scaled: float = 2.0
    center_wavelengths: float = 0.75
    stop_before_wrap: bool = True
    boundary_guard_wavelengths: float = 0.25
    boundary_guard_height: float = 0.01


def read_front_runner(parser: configparser.ConfigParser) -> FrontRunnerSettings:
    section = "front_runner"
    if not parser.has_section(section):
        return FrontRunnerSettings()
    return FrontRunnerSettings(
        enabled=parser.getboolean(section, "enabled", fallback=False),
        wavelength_slope_scaled=parser.getfloat(section, "wavelength_slope_scaled", fallback=2.0),
        center_wavelengths=parser.getfloat(section, "center_wavelengths", fallback=0.75),
        stop_before_wrap=parser.getboolean(section, "stop_before_wrap", fallback=True),
        boundary_guard_wavelengths=parser.getfloat(section, "boundary_guard_wavelengths", fallback=0.25),
        boundary_guard_height=parser.getfloat(section, "boundary_guard_height", fallback=0.01),
    )


def read_domain_length(parser: configparser.ConfigParser) -> float:
    section = "domain_and_mesh"
    if parser.has_option(section, "lx_slope_scaled"):
        if parser.has_option(section, "lx"):
            raise ValueError("Specify only lx_slope_scaled (= S0 Lx/H_l) OR lx (= Lx/H_l), not both")
        slope = _required_float(parser, "physical_problem", "slope_tan")
        if slope <= 0.0:
            raise ValueError("slope_tan must be positive")
        return _required_float(parser, section, "lx_slope_scaled") / slope
    return _required_float(parser, section, "lx")


@dataclass(frozen=True)
class Case:
    config_path: Path
    front_runner: FrontRunnerSettings
    froude: float
    slope: float
    depth_ratio: float
    density_ratio: float
    interfacial_viscosity_ratio: float
    air_density_ratio: float
    upper_to_air_viscosity_ratio: float
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


def _required_air_viscosity_ratio(
    parser: configparser.ConfigParser,
) -> float:
    section = "physical_problem"
    legacy_name = "air_dynamic_viscosity"
    ratio_name = "upper_to_air_dynamic_viscosity_ratio"
    if parser.has_option(section, legacy_name):
        raise ValueError(
            f"[{section}] {legacy_name} is no longer an independent input; "
            f"remove it and set {ratio_name} = eta_upper(gamma=1)/eta_air"
        )
    return _required_float(parser, section, ratio_name)




def _required_interfacial_viscosity_ratio(
    parser: configparser.ConfigParser,
) -> float:
    """Read R_eta,I and reject the former independent kappa_K input."""
    section = "physical_problem"
    legacy_name = "consistency_ratio"
    ratio_name = "interfacial_apparent_viscosity_ratio"
    if parser.has_option(section, legacy_name):
        raise ValueError(
            f"[{section}] {legacy_name} is no longer an independent input; "
            f"remove it and set {ratio_name} = eta_lower,I/eta_upper,I"
        )
    return _required_float(parser, section, ratio_name)


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
        front_runner=read_front_runner(parser),
        froude=_required_float(parser, "physical_problem", "froude_lower"),
        slope=_required_float(parser, "physical_problem", "slope_tan"),
        depth_ratio=_required_float(parser, "physical_problem", "depth_ratio"),
        density_ratio=_required_float(parser, "physical_problem", "density_ratio"),
        interfacial_viscosity_ratio=(
            _required_interfacial_viscosity_ratio(parser)
        ),
        air_density_ratio=_required_float(
            parser, "physical_problem", "air_density_ratio"
        ),
        upper_to_air_viscosity_ratio=_required_air_viscosity_ratio(parser),
        air_streamwise_gravity=_required_bool(
            parser, "physical_problem", "air_streamwise_gravity"
        ),
        sigma_internal=_required_float(
            parser, "physical_problem", "sigma_internal_star"
        ),
        sigma_free=_required_float(parser, "physical_problem", "sigma_free_star"),
        lower=read_rheology(parser, "lower_rheology"),
        upper=read_rheology(parser, "upper_rheology"),
        lx=read_domain_length(parser),
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
        "interfacial_apparent_viscosity_ratio": (
            c.interfacial_viscosity_ratio
        ),
        "air_density_ratio": c.air_density_ratio,
        "upper_to_air_dynamic_viscosity_ratio": (
            c.upper_to_air_viscosity_ratio
        ),
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
    fr = c.front_runner
    if fr.enabled:
        vals = [fr.wavelength_slope_scaled, fr.center_wavelengths,
                fr.boundary_guard_wavelengths, fr.boundary_guard_height]
        if not all(math.isfinite(v) and v > 0.0 for v in vals):
            raise ValueError("All front_runner length/guard values must be finite and positive")
        wavelength = fr.wavelength_slope_scaled / c.slope
        center = fr.center_wavelengths * wavelength
        guard = fr.boundary_guard_wavelengths * wavelength
        if not (guard < center - wavelength/4.0 < center + wavelength/4.0 < c.lx - guard):
            raise ValueError("Localized support must be inside the domain and outside both boundary guards")
        if c.air_streamwise_gravity:
            raise ValueError("Front-runner air extension requires air_streamwise_gravity = false")
        if abs(c.lower_depth_amplitude-c.upper_depth_amplitude) > 1.e-14:
            raise ValueError("Constant-Fr homothetic mapping requires equal lower/upper relative depth amplitudes")
        if not (0.0 <= c.lower_depth_amplitude < 1.0):
            raise ValueError("The positive upper-half sinusoid requires 0 <= depth amplitude < 1")
        if c.velocity_amplitude != 0.0:
            raise ValueError("Set velocity_amplitude = 0: constant-Fr velocity is derived from depth, not separately forced")
        if c.perturbation_phase != 0.0:
            raise ValueError("Set perturbation_phase = 0; use center_wavelengths to position the localized hump")
        if c.pressure_mode != 1:
            raise ValueError("Front runner requires pressure_mode = 1 (hydrostatic pressure on the disturbed geometry)")
        if c.liquid_depth*(1+c.lower_depth_amplitude) >= c.ceiling*c.ceiling_fine_fraction:
            raise ValueError("Initial hump reaches the ceiling coarsening band")
        if wavelength/2.0/c.delta_min < 16:
            raise ValueError("Resolve the half-sinusoid support with at least 16 finest cells")
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




def unclipped_viscosity_components(
    gamma: float,
    r: Rheology,
) -> tuple[float, float]:
    """Return eta_raw = Lambda*basis + offset at a positive shear rate.

    The decomposition mirrors :func:`apparent_viscosity` before final clipping.
    It is used only to derive the fixed upper consistency coefficient from the
    prescribed steady-reference interfacial apparent-viscosity ratio.
    """
    gabs = abs(float(gamma))
    gsafe = max(gabs, r.gamma_floor)

    if r.model == "smooth_power_law":
        basis = (gabs * gabs + r.epsilon * r.epsilon) ** (
            0.5 * (r.n - 1.0)
        )
    else:
        basis = gsafe ** (r.n - 1.0)

    offset = 0.0
    if r.model == "papanastasiou_hb" and r.tau0 > 0.0:
        x = r.m * gabs
        if abs(x) < 1.0e-5:
            ratio = r.m * (1.0 - x / 2.0 + x * x / 6.0 - x**3 / 24.0)
        else:
            ratio = -math.expm1(-x) / gsafe
        offset = r.tau0 * ratio
    elif r.model == "bounded_hb" and r.tau0 > 0.0:
        offset = r.tau0 / gsafe

    return float(basis), float(offset)


def derive_lambda_upper_from_interfacial_ratio(
    c: Case,
    lam_lower: float,
    gamma_lower_i: float,
) -> dict[str, float]:
    """Derive Lambda_u from the requested steady interfacial viscosity ratio.

    At the compatible positive-shear interface,

        R_eta,I = eta_l,I/eta_u,I = gamma_u,I/gamma_l,I,

    because the tangential traction is continuous.  The material coefficient
    Lambda_u is fixed once at generation time; R_eta,I is not imposed as a
    time-dependent constraint during the CFD calculation.
    """
    ratio = c.interfacial_viscosity_ratio
    eta_lower_i = float(
        apparent_viscosity(
            np.asarray([gamma_lower_i]), lam_lower, c.lower
        )[0]
    )
    gamma_upper_i = ratio * gamma_lower_i
    eta_upper_target = eta_lower_i / ratio

    basis, offset = unclipped_viscosity_components(gamma_upper_i, c.upper)
    tol = 5.0e-12 * max(1.0, abs(eta_upper_target))

    if eta_upper_target <= c.upper.eta_min + tol:
        raise ValueError(
            "The requested R_eta,I puts the upper interfacial apparent "
            "viscosity on/below eta_min; the required Lambda_u is not uniquely "
            "defined. Change R_eta,I or the upper viscosity bound."
        )
    if eta_upper_target >= c.upper.eta_max - tol:
        raise ValueError(
            "The requested R_eta,I puts the upper interfacial apparent "
            "viscosity on/above eta_max; the required Lambda_u is not uniquely "
            "defined. Change R_eta,I or the upper viscosity bound."
        )
    if not math.isfinite(basis) or basis <= 0.0:
        raise ValueError("Upper rheology has a nonpositive viscosity basis")

    lam_upper = (eta_upper_target - offset) / basis
    if not math.isfinite(lam_upper) or lam_upper <= 0.0:
        raise ValueError(
            "The requested R_eta,I is incompatible with the upper yield/"
            "regularization contribution: it would require Lambda_u <= 0."
        )

    eta_upper_i = float(
        apparent_viscosity(
            np.asarray([gamma_upper_i]), lam_upper, c.upper
        )[0]
    )
    ratio_realized = eta_lower_i / eta_upper_i
    tau_lower_i = eta_lower_i * gamma_lower_i
    tau_upper_i = eta_upper_i * gamma_upper_i
    kappa_derived = lam_upper / lam_lower

    rel_scale = max(1.0, abs(ratio), abs(tau_lower_i), abs(tau_upper_i))
    if abs(ratio_realized - ratio) > 2.0e-10 * rel_scale:
        raise RuntimeError(
            "Derived Lambda_u does not realize the requested interfacial "
            "apparent-viscosity ratio"
        )
    if abs(tau_lower_i - tau_upper_i) > 2.0e-10 * rel_scale:
        raise RuntimeError("Interfacial tangential traction is not continuous")

    return {
        "lambda_upper": lam_upper,
        "derived_consistency_ratio": kappa_derived,
        "gamma_lower_interface": gamma_lower_i,
        "gamma_upper_interface_target": gamma_upper_i,
        "eta_lower_interface": eta_lower_i,
        "eta_upper_interface": eta_upper_i,
        "interfacial_ratio_realized": ratio_realized,
        "interfacial_traction_lower": tau_lower_i,
        "interfacial_traction_upper": tau_upper_i,
        "interfacial_traction_mismatch": tau_upper_i - tau_lower_i,
    }


def derive_air_viscosity(c: Case, lam_upper: float) -> tuple[float, float]:
    """Return upper reference and derived air dynamic viscosities."""
    upper_reference_eta = float(
        apparent_viscosity(
            np.asarray([UPPER_VISCOSITY_REFERENCE_SHEAR_RATE]),
            lam_upper,
            c.upper,
        )[0]
    )
    if not math.isfinite(upper_reference_eta) or upper_reference_eta <= 0.0:
        raise ValueError(
            "The upper reference dynamic viscosity must be finite and positive"
        )
    air_eta = upper_reference_eta / c.upper_to_air_viscosity_ratio
    if not math.isfinite(air_eta) or air_eta <= 0.0:
        raise ValueError(
            "The derived air dynamic viscosity must be finite and positive; "
            "check upper_to_air_dynamic_viscosity_ratio"
        )
    return upper_reference_eta, air_eta


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

    z1 = np.linspace(0.0, 1.0, c.integration_points)
    y2 = np.linspace(0.0, c.depth_ratio, c.integration_points)
    y3 = np.linspace(0.0, c.air_depth, c.integration_points)
    tau1, tau2, tau3 = traction_fields(c, z1, y2, y3)

    gamma1 = invert_stress(tau1, lam_lower, c.lower)
    mapping = derive_lambda_upper_from_interfacial_ratio(
        c,
        lam_lower,
        float(gamma1[-1]),
    )
    lam_upper = mapping["lambda_upper"]
    upper_reference_eta, air_eta = derive_air_viscosity(c, lam_upper)

    gamma2 = invert_stress(tau2, lam_upper, c.upper)
    gamma3 = tau3 / air_eta

    u1 = cumulative_trapezoid(gamma1, z1)
    interface_velocity = float(u1[-1])
    u2 = interface_velocity + cumulative_trapezoid(gamma2, y2)
    surface_velocity = float(u2[-1])
    u3 = surface_velocity + cumulative_trapezoid(gamma3, y3)

    mean1 = trapezoid(u1, z1)
    mean2 = trapezoid(u2, y2) / c.depth_ratio
    mean3 = trapezoid(u3, y3) / c.air_depth

    gamma_upper_actual = float(gamma2[0])
    eta_upper_actual = float(
        apparent_viscosity(
            np.asarray([gamma_upper_actual]), lam_upper, c.upper
        )[0]
    )
    ratio_actual = mapping["eta_lower_interface"] / eta_upper_actual
    tau_upper_actual = eta_upper_actual * gamma_upper_actual
    interface_scale = max(
        1.0,
        abs(mapping["interfacial_traction_lower"]),
        abs(tau_upper_actual),
    )
    if abs(ratio_actual - c.interfacial_viscosity_ratio) > 5.0e-9 * max(
        1.0, abs(c.interfacial_viscosity_ratio)
    ):
        raise RuntimeError(
            "Full upper constitutive inversion does not realize the requested "
            "R_eta,I"
        )
    if abs(tau_upper_actual - mapping["interfacial_traction_lower"]) > 5.0e-9 * interface_scale:
        raise RuntimeError(
            "Full upper constitutive inversion violates interfacial traction "
            "continuity"
        )

    return {
        "lambda_lower": lam_lower,
        "lambda_upper": lam_upper,
        "derived_consistency_ratio": mapping["derived_consistency_ratio"],
        "upper_reference_eta": upper_reference_eta,
        "air_eta": air_eta,
        "gamma_lower_interface": mapping["gamma_lower_interface"],
        "gamma_upper_interface": gamma_upper_actual,
        "eta_lower_interface": mapping["eta_lower_interface"],
        "eta_upper_interface": eta_upper_actual,
        "interfacial_ratio_realized": ratio_actual,
        "interfacial_traction_lower": mapping["interfacial_traction_lower"],
        "interfacial_traction_upper": tau_upper_actual,
        "interfacial_traction_mismatch": (
            tau_upper_actual - mapping["interfacial_traction_lower"]
        ),
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
    eta[air] = base["air_eta"]

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
        (
            "CASE_INTERFACIAL_APPARENT_VISCOSITY_RATIO",
            c_float(c.interfacial_viscosity_ratio),
        ),
        (
            "CASE_DERIVED_CONSISTENCY_RATIO",
            c_float(base["derived_consistency_ratio"]),
        ),
        ("CASE_CONSISTENCY_RATIO", "CASE_DERIVED_CONSISTENCY_RATIO"),
        ("CASE_AIR_DENSITY_RATIO", c_float(c.air_density_ratio)),
        (
            "CASE_UPPER_TO_AIR_VISCOSITY_RATIO",
            c_float(c.upper_to_air_viscosity_ratio),
        ),
        (
            "CASE_UPPER_VISCOSITY_REFERENCE_SHEAR",
            c_float(UPPER_VISCOSITY_REFERENCE_SHEAR_RATE),
        ),
        ("CASE_UPPER_REFERENCE_ETA", c_float(base["upper_reference_eta"])),
        ("CASE_AIR_ETA", c_float(base["air_eta"])),
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
        (
            "CASE_LOWER_INTERFACE_GAMMA",
            c_float(base["gamma_lower_interface"]),
        ),
        (
            "CASE_UPPER_INTERFACE_GAMMA",
            c_float(base["gamma_upper_interface"]),
        ),
        (
            "CASE_LOWER_INTERFACE_ETA",
            c_float(base["eta_lower_interface"]),
        ),
        (
            "CASE_UPPER_INTERFACE_ETA",
            c_float(base["eta_upper_interface"]),
        ),
        (
            "CASE_REALIZED_INTERFACE_VISCOSITY_RATIO",
            c_float(base["interfacial_ratio_realized"]),
        ),
        (
            "CASE_INTERFACE_TRACTION_MISMATCH",
            c_float(base["interfacial_traction_mismatch"]),
        ),
        ("CASE_FRONT_RUNNER_ENABLED", c_bool(c.front_runner.enabled)),
        ("CASE_FR_WAVELENGTH", c_float(c.front_runner.wavelength_slope_scaled/c.slope)),
        ("CASE_FR_CENTER", c_float(c.front_runner.center_wavelengths*c.front_runner.wavelength_slope_scaled/c.slope)),
        ("CASE_FR_STOP_BEFORE_WRAP", c_bool(c.front_runner.stop_before_wrap)),
        ("CASE_FR_GUARD_WIDTH", c_float(c.front_runner.boundary_guard_wavelengths*c.front_runner.wavelength_slope_scaled/c.slope)),
        ("CASE_FR_GUARD_HEIGHT", c_float(c.front_runner.boundary_guard_height)),
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
        stream.write(
            "/* Generated file. Edit case_parameters.ini, not this header.\n"
            "   R_eta,I is prescribed at the compatible steady interface.\n"
            "   Lambda_u and kappa_K are derived once and remain fixed. */\n"
        )
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
    if c.front_runner.enabled:
        from front_runner_initialization import write_initial_tables
        write_initial_tables(output_dir, c, z, velocity)
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
R_eta,I prescribed    = {c.interfacial_viscosity_ratio:.12g}
R_eta,I realized      = {base['interfacial_ratio_realized']:.12g}
epsilon_l, epsilon_u  = {c.lower.epsilon:.12g}, {c.upper.epsilon:.12g}

Explicit-air and capillary groups
---------------------------------
r_a                   = {c.air_density_ratio:.12g}
R_mu,u/a              = {c.upper_to_air_viscosity_ratio:.12g}
gamma_ref,u           = {UPPER_VISCOSITY_REFERENCE_SHEAR_RATE:.12g}
eta_u,ref             = {base['upper_reference_eta']:.12g}
eta_a (derived)       = {base['air_eta']:.12g}
air streamwise gravity= {int(c.air_streamwise_gravity)}
sigma*_internal       = {c.sigma_internal:.12g}   (We={we_internal:.12g})
sigma*_free           = {c.sigma_free:.12g}   (We={we_free:.12g})

Derived compatibility quantities
--------------------------------
Lambda_l              = {base['lambda_lower']:.12g}
Lambda_u              = {base['lambda_upper']:.12g}
kappa_K (derived)     = {base['derived_consistency_ratio']:.12g}
gamma_l,I             = {base['gamma_lower_interface']:.12g}
gamma_u,I             = {base['gamma_upper_interface']:.12g}
eta_l,I               = {base['eta_lower_interface']:.12g}
eta_u,I               = {base['eta_upper_interface']:.12g}
traction mismatch I   = {base['interfacial_traction_mismatch']:.6e}
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
