#!/usr/bin/env python3
"""Run one full-Navier--Stokes two-layer stability calculation.

Normal use:

    python run_dispersion_case.py

The parameters that most users need are collected in the clearly marked
USER INPUT section below.  Command-line options can override the output
folder and Chebyshev order without editing the file.
"""
from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
import argparse
import json
import math
import sys

import numpy as np

from fullns_stability import (
    CaseParameters,
    FluidLayer,
    SolverOptions,
    apparent_and_tangent_viscosity,
    build_base_flow,
    convergence_at_wavenumber,
    neutral_crossings,
    track_hydrodynamic_modes,
)


# =============================================================================
# USER INPUT
# =============================================================================

CASE = CaseParameters(
    gravity=9.81,
    slope_tan=0.05,
    lower=FluidLayer(
        depth=0.025,          # m
        density=1500.0,       # kg/m^3
        consistency=6.0,      # Pa s^n
        flow_index=0.40,
    ),
    upper=FluidLayer(
        depth=0.025,          # m
        density=1200.0,       # kg/m^3
        consistency=3.0,      # Pa s^n
        flow_index=0.80,
    ),
    regularization_rate=0.5,  # s^-1
)

SOLVER = SolverOptions(
    chebyshev_order=64,
    residual_tolerance=2.0e-6,
    phase_speed_limit=100.0,
    seed_wavenumber=5.0e-4,
    base_samples=6001,
)

# A hybrid grid resolves the weak long-wave interface branch and the finite
# cutoff of the surface branch without wasting all points at large k.
WAVENUMBERS = np.unique(
    np.concatenate(
        [
            np.geomspace(5.0e-4, 5.0e-2, 60),
            np.linspace(5.0e-2, 5.0e-1, 90),
            np.linspace(5.0e-1, 3.0, 56),
        ]
    )
)

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "data"

# =============================================================================
# END USER INPUT
# =============================================================================


def write_base_profiles(base, output_directory: Path) -> Path:
    """Write one row per dense base-state point."""

    lower_mu, lower_mu_t = apparent_and_tangent_viscosity(
        base.lower.shear,
        base.consistency_number_lower,
        base.parameters.lower.flow_index,
        base.epsilon_lower,
    )
    upper_mu, upper_mu_t = apparent_and_tangent_viscosity(
        base.upper.shear,
        base.consistency_number_upper,
        base.parameters.upper.flow_index,
        base.epsilon_upper,
    )

    path = output_directory / "base_state_profiles.tsv"
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# layer\tz\tvelocity\tpressure\tshear\tapparent_viscosity\t"
            "tangent_viscosity\n"
        )
        for name, profile, mu, mu_t in (
            ("lower", base.lower, lower_mu, lower_mu_t),
            ("upper", base.upper, upper_mu, upper_mu_t),
        ):
            for values in zip(
                profile.z,
                profile.velocity,
                profile.pressure,
                profile.shear,
                mu,
                mu_t,
            ):
                handle.write(
                    name
                    + "\t"
                    + "\t".join(f"{float(value):.14e}" for value in values)
                    + "\n"
                )
    return path


def write_dispersion_table(wavenumbers, tracks, output_directory: Path) -> Path:
    path = output_directory / "dispersion_branches.tsv"
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# k\tbranch\tc_real\tc_imag\tsigma_real\tsigma_imag\t"
            "growth_rate\tcelerity\teta_real\teta_imag\tzeta_real\t"
            "zeta_imag\tamplitude_ratio\tphase_difference_deg\tresidual\n"
        )
        for branch in ("interface_varicose", "surface_zigzag"):
            for k, mode in zip(wavenumbers, tracks[branch]):
                handle.write(
                    f"{float(k):.14e}\t{branch}\t"
                    f"{mode.phase_speed.real:.14e}\t{mode.phase_speed.imag:.14e}\t"
                    f"{mode.sigma.real:.14e}\t{mode.sigma.imag:.14e}\t"
                    f"{mode.growth_rate:.14e}\t{mode.celerity:.14e}\t"
                    f"{mode.eta.real:.14e}\t{mode.eta.imag:.14e}\t"
                    f"{mode.zeta.real:.14e}\t{mode.zeta.imag:.14e}\t"
                    f"{mode.amplitude_ratio:.14e}\t"
                    f"{mode.phase_difference_degrees:.14e}\t"
                    f"{mode.residual:.14e}\n"
                )
    return path


def write_conditioning_table(diagnostics, output_directory: Path) -> Path:
    path = output_directory / "conditioning_diagnostics.tsv"
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# k\trow_norm_spread_before_equilibration\t"
            "finite_modes_retained\tmatrix_size\n"
        )
        for row in diagnostics:
            handle.write(
                f"{row['wavenumber']:.14e}\t{row['row_norm_spread']:.14e}\t"
                f"{row['finite_modes_retained']:.0f}\t{row['matrix_size']:.0f}\n"
            )
    return path


def write_convergence_table(rows, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# order\tinterface_growth\tinterface_celerity\tinterface_residual\t"
            "surface_growth\tsurface_celerity\tsurface_residual\n"
        )
        for row in rows:
            handle.write(
                f"{row['order']:.0f}\t"
                f"{row['interface_varicose_growth']:.14e}\t"
                f"{row['interface_varicose_celerity']:.14e}\t"
                f"{row['interface_varicose_residual']:.14e}\t"
                f"{row['surface_zigzag_growth']:.14e}\t"
                f"{row['surface_zigzag_celerity']:.14e}\t"
                f"{row['surface_zigzag_residual']:.14e}\n"
            )


def branch_metrics(wavenumbers, branch_modes) -> dict[str, object]:
    growth = np.array([mode.growth_rate for mode in branch_modes], dtype=float)
    celerity = np.array([mode.celerity for mode in branch_modes], dtype=float)
    maximum_index = int(np.argmax(growth))
    roots = neutral_crossings(wavenumbers, growth)
    return {
        "maximum_growth": float(growth[maximum_index]),
        "wavenumber_at_maximum": float(wavenumbers[maximum_index]),
        "celerity_at_maximum": float(celerity[maximum_index]),
        "neutral_wavenumbers": roots,
    }


def dimensional_metrics(base, dimensionless: dict[str, object]) -> dict[str, object]:
    result = dict(dimensionless)
    k = float(dimensionless["wavenumber_at_maximum"])
    result["dimensional_wavelength_at_maximum_m"] = (
        2.0 * math.pi * base.parameters.lower.depth / k
    )
    result["dimensional_growth_at_maximum_per_s"] = (
        float(dimensionless["maximum_growth"])
        * base.velocity_scale
        / base.parameters.lower.depth
    )
    result["dimensional_celerity_at_maximum_m_per_s"] = (
        float(dimensionless["celerity_at_maximum"]) * base.velocity_scale
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--order", type=int, default=SOLVER.chebyshev_order)
    parser.add_argument(
        "--regularization-rate",
        type=float,
        default=CASE.regularization_rate,
    )
    parser.add_argument(
        "--with-convergence",
        action="store_true",
        help="also run the more expensive multi-order convergence tables",
    )
    args = parser.parse_args()

    output_directory = args.output.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    case = replace(CASE, regularization_rate=args.regularization_rate)
    solver = replace(SOLVER, chebyshev_order=args.order)
    case.validate()
    solver.validate()

    base = build_base_flow(case, samples=solver.base_samples)
    tracks, diagnostics = track_hydrodynamic_modes(base, WAVENUMBERS, solver)

    base_path = write_base_profiles(base, output_directory)
    dispersion_path = write_dispersion_table(WAVENUMBERS, tracks, output_directory)
    conditioning_path = write_conditioning_table(diagnostics, output_directory)

    if args.with_convergence:
        longwave_rows = convergence_at_wavenumber(
            case,
            wavenumber=0.01,
            orders=(32, 40, 52, 64, 76),
            seed_wavenumber=solver.seed_wavenumber,
        )
        finite_rows = convergence_at_wavenumber(
            case,
            wavenumber=0.4,
            orders=(32, 40, 52, 64, 76),
            seed_wavenumber=solver.seed_wavenumber,
        )
        write_convergence_table(
            longwave_rows,
            output_directory / "convergence_k0p01.tsv",
        )
        write_convergence_table(
            finite_rows,
            output_directory / "convergence_k0p4.tsv",
        )

    branch_summaries = {
        branch: dimensional_metrics(
            base,
            branch_metrics(WAVENUMBERS, modes),
        )
        for branch, modes in tracks.items()
    }

    summary = {
        "case_parameters": asdict(case),
        "solver_options": asdict(solver),
        "base_flow": base.summary(),
        "wavenumber_count": int(WAVENUMBERS.size),
        "wavenumber_min": float(WAVENUMBERS.min()),
        "wavenumber_max": float(WAVENUMBERS.max()),
        "branches": branch_summaries,
        "files": {
            "base_state_profiles": str(base_path.name),
            "dispersion_branches": str(dispersion_path.name),
            "conditioning_diagnostics": str(conditioning_path.name),
            "convergence_longwave": "convergence_k0p01.tsv",
            "convergence_finite": "convergence_k0p4.tsv",
        },
    }
    summary_path = output_directory / "case_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print(f"\nWrote results to: {output_directory}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # concise user-facing failure message
        print(f"ERROR: {error}", file=sys.stderr)
        raise
