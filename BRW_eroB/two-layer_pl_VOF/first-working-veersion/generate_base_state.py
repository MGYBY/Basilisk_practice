#!/usr/bin/env python3
"""Generate the analytical one-dimensional base state used by roll_wave_amr.c.

Workflow
--------
1. Edit the clearly marked ``USER INPUT PARAMETERS`` block below.
2. Run ``python3 generate_base_state.py``.
3. Edit the matching physical/geometrical parameters in ``roll_wave_amr.c``.
4. Compile and run the Basilisk case with ``./Allrun``.

The script writes only ordinary text data into ``base_state/``.  It does not
create C headers and it does not modify the Basilisk source.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np


@dataclass(frozen=True)
class Parameters:
    """Inputs for the dimensional two-layer power-law normal-flow solution."""

    gravity: float
    slope_tan: float
    h1: float
    h2: float
    nominal_height: float
    lx_factor: float
    max_level: int
    rho1: float
    rho2: float
    rho3: float
    k1: float
    n1: float
    k2: float
    n2: float
    mu1_max: float
    mu2_max: float
    mu_air: float
    gamma_min: float
    air_profile: str
    lower_samples: int
    upper_samples: int
    air_samples: int

    @property
    def lx(self) -> float:
        return self.h1 * self.lx_factor

    @property
    def delta_min(self) -> float:
        return self.lx / (1 << self.max_level)

    @property
    def ceiling(self) -> float:
        # Same alignment rule as roll_wave_amr.c.
        return math.ceil(
            self.nominal_height / self.delta_min - 1.0e-12
        ) * self.delta_min

    @property
    def liquid_depth(self) -> float:
        return self.h1 + self.h2

    @property
    def sin_theta(self) -> float:
        return self.slope_tan / math.sqrt(1.0 + self.slope_tan**2)

    @property
    def cos_theta(self) -> float:
        return 1.0 / math.sqrt(1.0 + self.slope_tan**2)

    @property
    def gs(self) -> float:
        return self.gravity * self.sin_theta

    @property
    def gn(self) -> float:
        return self.gravity * self.cos_theta


# =============================================================================
# USER INPUT PARAMETERS
# =============================================================================
# Edit this block directly.  The matching entries in roll_wave_amr.c are
# marked [SYNC] and must be assigned the same values.
#
# Power-law convention:
#     tau = K |du/dy|^(n-1) du/dy
# so K has units Pa s^n.
PARAMETERS = Parameters(
    # Geometry and gravity
    gravity=9.81,                 # m/s^2
    slope_tan=0.06,               # tan(theta)
    h1=0.025,                      # lower-layer depth [m]
    h2=0.025,                      # upper-layer depth [m]
    nominal_height=0.11,         # requested embedded-ceiling height [m]
    lx_factor=100.0,              # Lx/H1
    max_level=11,                 # must match MAXLEVEL in roll_wave_amr.c

    # Phase densities
    rho1=1500.0,                  # lower liquid [kg/m^3]
    rho2=1200.0,                  # upper liquid [kg/m^3]
    rho3=10.00,                     # passive ambient density [kg/m^3]

    # Lower and upper power-law rheology
    k1=6.0,                       # lower consistency [Pa s^n]
    n1=0.4,                       # lower power-law index
    k2=3.0,                       # upper consistency [Pa s^n]
    n2=0.8,                       # upper power-law index
    mu1_max=40.0,                 # lower apparent-viscosity cap [Pa s]
    mu2_max=40.0,                 # upper apparent-viscosity cap [Pa s]
    gamma_min=1.0e-8,            # shear-rate floor [1/s]

    # Passive ambient: constant U=U_surface, with zero tangential stress.
    mu_air=0.01,                # dynamic viscosity [Pa s]
    air_profile="constant",       # only compatible option in this release

    # Text-profile resolution
    lower_samples=101,
    upper_samples=101,
    air_samples=101,
)

# Files are written beside this script, under cases/roll_wave_amr/base_state.
OUTPUT_DIRECTORY = Path(__file__).resolve().parent / "base_state"
# =============================================================================
# END USER INPUT PARAMETERS
# =============================================================================


def validate(p: Parameters) -> None:
    positive = {
        "gravity": p.gravity,
        "h1": p.h1,
        "h2": p.h2,
        "nominal_height": p.nominal_height,
        "lx_factor": p.lx_factor,
        "rho1": p.rho1,
        "rho2": p.rho2,
        "rho3": p.rho3,
        "k1": p.k1,
        "k2": p.k2,
        "n1": p.n1,
        "n2": p.n2,
        "mu1_max": p.mu1_max,
        "mu2_max": p.mu2_max,
        "mu_air": p.mu_air,
        "gamma_min": p.gamma_min,
    }
    invalid = [
        name for name, value in positive.items()
        if not math.isfinite(value) or value <= 0.0
    ]
    if invalid:
        raise ValueError(
            "These parameters must be finite and positive: "
            + ", ".join(invalid)
        )
    if p.max_level < 1:
        raise ValueError("max_level must be at least one")
    if p.ceiling <= p.liquid_depth:
        raise ValueError("The aligned ceiling must lie above h1 + h2")
    if p.air_profile != "constant":
        raise ValueError(
            "air_profile must be 'constant' for the passive-ambient, "
            "free-slip validation problem"
        )
    for name, count in {
        "lower_samples": p.lower_samples,
        "upper_samples": p.upper_samples,
        "air_samples": p.air_samples,
    }.items():
        if count < 2:
            raise ValueError(f"{name} must be at least two")


def traction_scales(p: Parameters) -> Tuple[float, float]:
    tau_interface = p.rho2 * p.gs * p.h2
    tau_bed = p.gs * (p.rho1 * p.h1 + p.rho2 * p.h2)
    return tau_bed, tau_interface


def lower_velocity(y: np.ndarray, p: Parameters) -> np.ndarray:
    tau_bed, _ = traction_scales(p)
    exponent = (p.n1 + 1.0) / p.n1
    tau = np.maximum(tau_bed - p.rho1 * p.gs * y, 0.0)
    coefficient = p.n1 / (p.n1 + 1.0)
    coefficient /= p.rho1 * p.gs * p.k1 ** (1.0 / p.n1)
    return coefficient * (tau_bed**exponent - tau**exponent)


def interface_velocity(p: Parameters) -> float:
    return float(lower_velocity(np.asarray([p.h1]), p)[0])


def upper_velocity(local_y: np.ndarray, p: Parameters) -> np.ndarray:
    ui = interface_velocity(p)
    exponent = (p.n2 + 1.0) / p.n2
    coefficient = p.n2 / (p.n2 + 1.0)
    coefficient *= (p.rho2 * p.gs / p.k2) ** (1.0 / p.n2)
    return ui + coefficient * (
        p.h2**exponent
        - np.maximum(p.h2 - local_y, 0.0) ** exponent
    )


def surface_velocity(p: Parameters) -> float:
    return float(upper_velocity(np.asarray([p.h2]), p)[0])


def air_velocity(y: np.ndarray, p: Parameters) -> np.ndarray:
    """Passive ambient translating uniformly with the liquid surface.

    With streamwise gravity disabled in the ambient and a free-slip ceiling,
    this constant profile has zero shear everywhere and is a steady solution
    compatible with the stress-free analytical liquid base state.
    """
    return np.full_like(np.asarray(y, dtype=float), surface_velocity(p))


def air_shear_rate(y: np.ndarray, p: Parameters) -> np.ndarray:
    del p
    return np.zeros_like(np.asarray(y, dtype=float))


def velocity(y: np.ndarray, p: Parameters) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    result = np.empty_like(y)
    lower = y <= p.h1
    upper = (y > p.h1) & (y <= p.liquid_depth)
    air = y > p.liquid_depth
    result[lower] = lower_velocity(np.maximum(y[lower], 0.0), p)
    result[upper] = upper_velocity(y[upper] - p.h1, p)
    result[air] = air_velocity(y[air], p)
    return result


def pressure(y: np.ndarray, p: Parameters) -> np.ndarray:
    """Hydrostatic mechanical pressure with p=0 at the embedded ceiling."""
    y = np.asarray(y, dtype=float)
    result = np.empty_like(y)
    lower = y <= p.h1
    upper = (y > p.h1) & (y <= p.liquid_depth)
    air = y > p.liquid_depth

    p_air_at_surface = p.rho3 * p.gn * (p.ceiling - p.liquid_depth)
    p_upper_at_interface = p_air_at_surface + p.rho2 * p.gn * p.h2

    result[air] = p.rho3 * p.gn * (p.ceiling - y[air])
    result[upper] = (
        p_air_at_surface
        + p.rho2 * p.gn * (p.liquid_depth - y[upper])
    )
    result[lower] = (
        p_upper_at_interface
        + p.rho1 * p.gn * (p.h1 - y[lower])
    )
    return result


def shear_rate(y: np.ndarray, p: Parameters) -> np.ndarray:
    """Return sqrt(2 D:D), equal to |du/dy| for the base flow."""
    y = np.asarray(y, dtype=float)
    result = np.zeros_like(y)
    tau_bed, _ = traction_scales(p)
    lower = y <= p.h1
    upper = (y > p.h1) & (y <= p.liquid_depth)
    air = y > p.liquid_depth
    tau1 = np.maximum(tau_bed - p.rho1 * p.gs * y[lower], 0.0)
    result[lower] = (tau1 / p.k1) ** (1.0 / p.n1)
    local_y = y[upper] - p.h1
    tau2 = np.maximum(p.rho2 * p.gs * (p.h2 - local_y), 0.0)
    result[upper] = (tau2 / p.k2) ** (1.0 / p.n2)
    result[air] = air_shear_rate(y[air], p)
    return result


def apparent_viscosity(y: np.ndarray, p: Parameters) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    gd = np.maximum(shear_rate(y, p), p.gamma_min)
    result = np.full_like(y, p.mu_air)
    lower = y <= p.h1
    upper = (y > p.h1) & (y <= p.liquid_depth)
    result[lower] = np.minimum(
        p.mu1_max, p.k1 * gd[lower] ** (p.n1 - 1.0)
    )
    result[upper] = np.minimum(
        p.mu2_max, p.k2 * gd[upper] ** (p.n2 - 1.0)
    )
    return result


def phase_fractions(
    y: np.ndarray,
    p: Parameters,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = np.asarray(y, dtype=float)
    f1 = (y <= p.h1).astype(float)
    f2 = ((y > p.h1) & (y <= p.liquid_depth)).astype(float)
    f3 = (y > p.liquid_depth).astype(float)
    return f1, f2, f3


def profile_coordinates(p: Parameters) -> np.ndarray:
    lower = np.linspace(0.0, p.h1, p.lower_samples)
    upper = np.linspace(p.h1, p.liquid_depth, p.upper_samples)[1:]
    air = np.linspace(p.liquid_depth, p.ceiling, p.air_samples)[1:]
    y = np.concatenate((lower, upper, air))
    if np.any(np.diff(y) <= 0.0):
        raise RuntimeError("Generated profile coordinates are not strictly increasing")
    return y


def integrate_piecewise(
    y: np.ndarray,
    values: np.ndarray,
    low: float,
    high: float,
) -> float:
    mask = (y >= low) & (y <= high)
    integrate = getattr(np, "trapezoid", np.trapz)
    return float(integrate(values[mask], y[mask]))


def write_two_column(
    path: Path,
    y: np.ndarray,
    values: np.ndarray,
    value_name: str,
    units: str,
) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write("# Generated by generate_base_state.py\n")
        stream.write(f"# columns: y [m]    {value_name} [{units}]\n")
        for yy, value in zip(y, values):
            stream.write(f"{yy:.16e}\t{value:.16e}\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_outputs(output_dir: Path, p: Parameters) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    y = profile_coordinates(p)
    u = velocity(y, p)
    hydrostatic_p = pressure(y, p)
    gd = shear_rate(y, p)
    mu = apparent_viscosity(y, p)
    f1, f2, f3 = phase_fractions(y, p)

    velocity_file = output_dir / "velocity.dat"
    pressure_file = output_dir / "pressure.dat"
    combined_file = output_dir / "combined_profile.tsv"
    summary_file = output_dir / "base_state_summary.txt"
    manifest_file = output_dir / "MANIFEST.sha256"

    write_two_column(velocity_file, y, u, "u_x", "m/s")
    write_two_column(pressure_file, y, hydrostatic_p, "p", "Pa")

    combined = np.column_stack(
        (y, f1, f2, f3, u, np.zeros_like(y), hydrostatic_p, gd, mu)
    )
    np.savetxt(
        combined_file,
        combined,
        delimiter="\t",
        header=(
            "y[m]\tf1\tf2\tf3\tu_x[m/s]\tu_y[m/s]\tp[Pa]"
            "\tgammaDot[1/s]\tmu[Pa.s]"
        ),
        comments="# ",
        fmt="%.16e",
    )

    tau_bed, tau_interface = traction_scales(p)

    # Use a dense integration grid so reported means do not depend materially
    # on the text-profile sampling selected above.
    integrate = getattr(np, "trapezoid", np.trapz)
    y1q = np.linspace(0.0, p.h1, 20001)
    y2q = np.linspace(p.h1, p.liquid_depth, 20001)
    y3q = np.linspace(p.liquid_depth, p.ceiling, 20001)
    q1 = float(integrate(velocity(y1q, p), y1q))
    q2 = float(integrate(velocity(y2q, p), y2q))
    q3 = float(integrate(velocity(y3q, p), y3q))
    ubar1 = q1 / p.h1
    ubar2 = q2 / p.h2
    fr1 = ubar1 / math.sqrt(p.gn * p.h1)
    fr2 = ubar2 / math.sqrt(p.gn * p.h2)
    ambient_values = air_velocity(y3q, p)

    # Both sides of the liquid--ambient interface carry zero tangential stress.
    tau_liquid_surface = 0.0
    tau_ambient_surface = p.mu_air * float(air_shear_rate(np.asarray([p.liquid_depth]), p)[0])
    tau_jump = tau_ambient_surface - tau_liquid_surface

    summary = f"""Dimensional analytical base state for Basilisk
================================================
Input source         = manually edited generate_base_state.py
Profile format       = plain two-column text
Ambient profile      = {p.air_profile}
Ambient mechanics    = passive translation; zero streamwise gravity; free-slip ceiling
MAXLEVEL             = {p.max_level}
Lx                   = {p.lx:.12g} m
Delta_min            = {p.delta_min:.12g} m
Nominal height       = {p.nominal_height:.12g} m
Aligned ceiling      = {p.ceiling:.12g} m
Cells through H1     = {p.h1/p.delta_min:.12g}
Cells through H2     = {p.h2/p.delta_min:.12g}

g                    = {p.gravity:.12g} m/s^2
tan(theta)           = {p.slope_tan:.12g}
sin(theta)           = {p.sin_theta:.12g}
cos(theta)           = {p.cos_theta:.12g}
g_s                  = {p.gs:.12g} m/s^2
g_n                  = {p.gn:.12g} m/s^2
H1                   = {p.h1:.12g} m
H2                   = {p.h2:.12g} m
rho1, rho2, rho3     = {p.rho1:.12g}, {p.rho2:.12g}, {p.rho3:.12g} kg/m^3
K1, n1               = {p.k1:.12g} Pa s^n, {p.n1:.12g}
K2, n2               = {p.k2:.12g} Pa s^n, {p.n2:.12g}
mu1_max, mu2_max     = {p.mu1_max:.12g}, {p.mu2_max:.12g} Pa s
mu_air               = {p.mu_air:.12g} Pa s
gamma_min            = {p.gamma_min:.12g} 1/s

tau_b                = {tau_bed:.12g} Pa
tau_I                = {tau_interface:.12g} Pa
tau_liquid(surface)  = {tau_liquid_surface:.12g} Pa
tau_ambient(surface) = {tau_ambient_surface:.12g} Pa
tangential jump      = {tau_jump:.12g} Pa
U_I                  = {interface_velocity(p):.12g} m/s
U_surface            = {surface_velocity(p):.12g} m/s
U_air(min,max)       = {ambient_values.min():.12g}, {ambient_values.max():.12g} m/s
Q1                   = {q1:.12g} m^2/s
Q2                   = {q2:.12g} m^2/s
Q_air                = {q3:.12g} m^2/s
Ubar1                = {ubar1:.12g} m/s
Ubar2                = {ubar2:.12g} m/s
Fr1                  = {fr1:.12g}
Fr2                  = {fr2:.12g}
p(ceiling)           = {pressure(np.asarray([p.ceiling]), p)[0]:.12g} Pa
p(bed)               = {pressure(np.asarray([0.0]), p)[0]:.12g} Pa
Profile rows         = {len(y)}
"""
    summary_file.write_text(summary, encoding="utf-8")

    generated: Iterable[Path] = (
        velocity_file,
        pressure_file,
        combined_file,
        summary_file,
    )
    manifest_file.write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in generated),
        encoding="utf-8",
    )
    print(summary, end="")


def main() -> None:
    validate(PARAMETERS)
    write_outputs(OUTPUT_DIRECTORY, PARAMETERS)
    print(f"Wrote analytical profiles to: {OUTPUT_DIRECTORY}")
    print("Next: verify the matching [SYNC] parameters in roll_wave_amr.c, then run ./Allrun")


if __name__ == "__main__":
    main()
