#!/usr/bin/env python3
"""Write the four growth-rate and phase-celerity branches versus wavenumber."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from branch_tracking import BRANCHES, hydrodynamic_amplitude_ratios, track_roots
from conventional_dispersion_analysis_detailed import (
    direct_critical_froude,
    longwave_modes,
)
from stability_user_inputs import (
    BASE_PARAMETERS,
    CASE_NAME,
    OUTPUT_ROOT,
    WAVENUMBER_MAX,
    WAVENUMBER_MIN,
    WAVENUMBER_POINTS,
    WAVENUMBER_SPACING,
)
from two_layer_powerlaw_stability_core import linear_coefficients, uniform_state


OUTPUT_BRANCH_ORDER = (
    "free_surface",
    "interfacial",
    "momentum_fast",
    "momentum_slow",
)


def _wavenumbers() -> np.ndarray:
    if WAVENUMBER_POINTS < 2:
        raise ValueError("WAVENUMBER_POINTS must be at least two")
    if not (0.0 < WAVENUMBER_MIN < WAVENUMBER_MAX):
        raise ValueError("Require 0 < WAVENUMBER_MIN < WAVENUMBER_MAX")

    spacing = WAVENUMBER_SPACING.strip().lower()
    if spacing == "log":
        return np.geomspace(
            WAVENUMBER_MIN,
            WAVENUMBER_MAX,
            WAVENUMBER_POINTS,
        )
    if spacing == "linear":
        return np.linspace(
            WAVENUMBER_MIN,
            WAVENUMBER_MAX,
            WAVENUMBER_POINTS,
        )
    raise ValueError("WAVENUMBER_SPACING must be 'log' or 'linear'")


def _write_table(path: Path, data: np.ndarray, header_lines: list[str]) -> None:
    np.savetxt(
        path,
        data,
        delimiter="\t",
        fmt="%.16e",
        header="\n".join(header_lines),
        comments="# ",
    )


def main() -> None:
    parameters = BASE_PARAMETERS.validated()
    k_values = _wavenumbers()
    tracked = track_roots(k_values, parameters)
    ratios = hydrodynamic_amplitude_ratios(tracked, parameters)

    output_dir = OUTPUT_ROOT / CASE_NAME / "dispersion"
    output_dir.mkdir(parents=True, exist_ok=True)

    root_arrays = {
        branch: np.asarray(tracked[branch], dtype=complex)
        for branch in OUTPUT_BRANCH_ORDER
    }

    growth = np.column_stack(
        [k_values] + [root_arrays[branch].real for branch in OUTPUT_BRANCH_ORDER]
    )
    growth_path = output_dir / "growth_rates_vs_wavenumber.txt"
    _write_table(
        growth_path,
        growth,
        [
            "Temporal growth rates from the four roots of the quartic dispersion relation",
            "normal-mode convention: exp(sigma*t + i*k*x)",
            "growth_rate = Re(sigma); positive means temporal growth",
            "branch labels are continuations from the k -> 0+ modes",
            "free_surface: hydrodynamic mode with mainly in-phase layer-depth motion",
            "interfacial: hydrodynamic mode with mainly out-of-phase layer-depth motion",
            "momentum_fast/slow: roots connected to the more/less strongly damped source-relaxation eigenvalues at k=0",
            "columns: k growth_free_surface growth_interfacial growth_momentum_fast growth_momentum_slow",
        ],
    )

    celerity = np.column_stack(
        [k_values]
        + [
            -root_arrays[branch].imag / k_values
            for branch in OUTPUT_BRANCH_ORDER
        ]
    )
    celerity_path = output_dir / "phase_celerities_vs_wavenumber.txt"
    _write_table(
        celerity_path,
        celerity,
        [
            "Phase celerities from the four roots of the quartic dispersion relation",
            "phase_celerity = -Im(sigma)/k for exp(sigma*t + i*k*x)",
            "the two momentum-relaxation celerities belong to damped relaxation roots and are not additional free-surface modes",
            "columns: k celerity_free_surface celerity_interfacial celerity_momentum_fast celerity_momentum_slow",
        ],
    )

    full_columns = [k_values]
    full_names = ["k"]
    for branch in OUTPUT_BRANCH_ORDER:
        full_columns.extend([root_arrays[branch].real, root_arrays[branch].imag])
        full_names.extend([f"sigma_{branch}_real", f"sigma_{branch}_imag"])
    full_columns.extend(
        [
            ratios["free_surface"].real,
            ratios["free_surface"].imag,
            ratios["interfacial"].real,
            ratios["interfacial"].imag,
        ]
    )
    full_names.extend(
        [
            "h_upper_over_h_lower_free_surface_real",
            "h_upper_over_h_lower_free_surface_imag",
            "h_upper_over_h_lower_interfacial_real",
            "h_upper_over_h_lower_interfacial_imag",
        ]
    )
    full_path = output_dir / "sigma_branches_vs_wavenumber.txt"
    _write_table(
        full_path,
        np.column_stack(full_columns),
        [
            "Complete branch-resolved dispersion data",
            "the amplitude-ratio columns help confirm the physical identity of the two hydrodynamic modes",
            "columns: " + " ".join(full_names),
        ],
    )

    branch_path = output_dir / "branch_definitions.txt"
    branch_path.write_text(
        "Branch definitions\n"
        "==================\n\n"
        "free_surface\n"
        "  The hydrodynamic root continuously connected to the k -> 0+ mode\n"
        "  whose lower- and upper-layer depth perturbations are mainly in phase.\n\n"
        "interfacial\n"
        "  The hydrodynamic root continuously connected to the k -> 0+ mode\n"
        "  whose layer-depth perturbations are mainly out of phase and whose\n"
        "  total free-surface displacement is relatively small.\n\n"
        "momentum_fast\n"
        "  The root connected at k=0 to the more negative source-relaxation\n"
        "  eigenvalue; it is the faster-decaying momentum adjustment mode.\n\n"
        "momentum_slow\n"
        "  The root connected at k=0 to the less negative source-relaxation\n"
        "  eigenvalue; it is the slower-decaying momentum adjustment mode.\n\n"
        "The labels are continuation labels. They are assigned at very small k\n"
        "and then followed continuously as k increases.\n",
        encoding="utf-8",
    )

    long_modes = {str(item["label"]): item for item in longwave_modes(parameters)}
    free_critical = direct_critical_froude(parameters, "free_surface", 0.05, 2.0)
    interface_critical = direct_critical_froude(parameters, "interfacial", 0.05, 2.0)
    coefficients = linear_coefficients(parameters)
    state_0 = uniform_state(parameters)

    summary_path = output_dir / "run_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("Two-layer power-law linear-stability dispersion run\n")
        handle.write("================================================\n\n")
        handle.write(f"case_name = {CASE_NAME}\n")
        handle.write(f"parameters = {parameters}\n")
        handle.write(f"uniform_state_[h_l,h_u,q_l,q_u] = {state_0.tolist()}\n")
        handle.write(
            f"wavenumber_grid = {WAVENUMBER_SPACING}, "
            f"{WAVENUMBER_MIN:g} to {WAVENUMBER_MAX:g}, "
            f"{WAVENUMBER_POINTS} points\n"
        )
        handle.write(
            "dimensionless outputs: k uses the report's streamwise length scale; "
            "sigma uses the inverse time scale; celerity uses the lower normal-flow velocity scale\n\n"
        )
        for branch in ("free_surface", "interfacial"):
            mode = long_modes[branch]
            handle.write(
                f"{branch}_longwave_celerity = {float(mode['celerity']):.16e}\n"
            )
            handle.write(
                f"{branch}_longwave_growth_coefficient_d = "
                f"{float(mode['growth_coefficient']):.16e}\n"
            )
            handle.write(
                f"{branch}_depth_amplitudes_[h_l,h_u] = "
                f"[{float(mode['h_lower_amplitude']):.16e}, "
                f"{float(mode['h_upper_amplitude']):.16e}]\n"
            )
        handle.write("\n")
        handle.write(f"free_surface_critical_result = {free_critical}\n")
        handle.write(f"interfacial_critical_result = {interface_critical}\n")
        handle.write(
            "source_relaxation_eigenvalues_at_k0 = "
            f"{np.linalg.eigvals(coefficients.C()).tolist()}\n"
        )
        handle.write(
            "maximum_normalized_branch_matching_cost = "
            f"{float(tracked['maximum_normalized_matching_cost']):.16e}\n"
        )

    print(f"Wrote {growth_path}")
    print(f"Wrote {celerity_path}")
    print(f"Wrote {full_path}")
    print(f"Wrote {branch_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
