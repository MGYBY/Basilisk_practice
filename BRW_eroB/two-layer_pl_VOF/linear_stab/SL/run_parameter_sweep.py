#!/usr/bin/env python3
"""Critical Froude number versus one selected model parameter.

The primary output has exactly two numerical columns:

    swept parameter    branch-specific critical Froude number

A second diagnostic table explains missing roots and also reports the other
hydrodynamic branch.  The critical value is obtained directly from the neutral
condition d_j(Fr_l)=0 by scan plus bisection; the optional affine shortcut is not
needed by this script.
"""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import numpy as np

from conventional_dispersion_analysis_detailed import (
    bisect_critical_froude,
    branch_growth_coefficient,
    find_sign_change,
)
from stability_user_inputs import (
    BASE_PARAMETERS,
    CASE_NAME,
    CRITICAL_BRANCH,
    FROUDE_ROOT_TOLERANCE,
    FROUDE_SCAN_POINTS,
    FROUDE_SEARCH_MAX,
    FROUDE_SEARCH_MIN,
    OUTPUT_ROOT,
    SWEEP_PARAMETER,
    SWEEP_VALUES,
)


PARAMETER_ALIASES = {
    "n_lower": "n_lower",
    "n_l": "n_lower",
    "depth_ratio": "depth_ratio",
    "h_r": "depth_ratio",
    "density_ratio": "density_ratio",
    "r_rho": "density_ratio",
    "n_upper": "n_upper",
    "n_u": "n_upper",
    "scaled_consistency_ratio": "scaled_consistency_ratio",
    "kappa_k": "scaled_consistency_ratio",
    "kappa_K": "scaled_consistency_ratio",
}

SYMBOL_NAMES = {
    "n_lower": "n_lower",
    "depth_ratio": "h_r",
    "density_ratio": "r_rho",
    "n_upper": "n_upper",
    "scaled_consistency_ratio": "kappa_K",
}


def _critical_for_branch(parameters, branch: str) -> dict[str, object]:
    bracket = find_sign_change(
        parameters,
        branch,
        FROUDE_SEARCH_MIN,
        FROUDE_SEARCH_MAX,
        FROUDE_SCAN_POINTS,
    )
    d_min = branch_growth_coefficient(
        FROUDE_SEARCH_MIN, parameters, branch
    )
    d_max = branch_growth_coefficient(
        FROUDE_SEARCH_MAX, parameters, branch
    )

    if bracket is None:
        if d_min > 0.0 and d_max > 0.0:
            status = "no_root__unstable_at_both_search_limits"
        elif d_min < 0.0 and d_max < 0.0:
            status = "no_root__stable_at_both_search_limits"
        else:
            status = "no_root__no_sign_change_resolved"
        return {
            "critical": np.nan,
            "status": status,
            "d_min": d_min,
            "d_max": d_max,
            "bracket_left": np.nan,
            "bracket_right": np.nan,
        }

    critical = bisect_critical_froude(
        parameters,
        branch,
        bracket,
        tolerance=FROUDE_ROOT_TOLERANCE,
    )
    delta = max(1.0e-6, 1.0e-4 * critical)
    d_below = branch_growth_coefficient(
        max(FROUDE_SEARCH_MIN, critical - delta), parameters, branch
    )
    d_above = branch_growth_coefficient(
        min(FROUDE_SEARCH_MAX, critical + delta), parameters, branch
    )
    if d_below < 0.0 < d_above:
        status = "stable_below__unstable_above"
    elif d_below > 0.0 > d_above:
        status = "unstable_below__stable_above"
    else:
        status = "root_found__local_sign_requires_inspection"

    return {
        "critical": float(critical),
        "status": status,
        "d_min": d_min,
        "d_max": d_max,
        "bracket_left": float(bracket[0]),
        "bracket_right": float(bracket[1]),
    }


def main() -> None:
    canonical = PARAMETER_ALIASES.get(SWEEP_PARAMETER)
    if canonical is None:
        allowed = ", ".join(sorted(PARAMETER_ALIASES))
        raise ValueError(
            f"Unsupported SWEEP_PARAMETER={SWEEP_PARAMETER!r}. Allowed: {allowed}"
        )
    if CRITICAL_BRANCH not in {"free_surface", "interfacial"}:
        raise ValueError(
            "CRITICAL_BRANCH must be 'free_surface' or 'interfacial'"
        )

    values = np.asarray(SWEEP_VALUES, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("SWEEP_VALUES must be a nonempty one-dimensional array")

    output_dir = OUTPUT_ROOT / CASE_NAME / "parameter_sweep"
    output_dir.mkdir(parents=True, exist_ok=True)
    symbol = SYMBOL_NAMES[canonical]

    simple_rows: list[list[float]] = []
    detailed_rows: list[dict[str, object]] = []

    for value in values:
        row: dict[str, object] = {
            "parameter_value": float(value),
            "error": "",
        }
        try:
            parameters = replace(
                BASE_PARAMETERS,
                **{canonical: float(value)},
            ).validated()
            free = _critical_for_branch(parameters, "free_surface")
            interface = _critical_for_branch(parameters, "interfacial")

            selected = free if CRITICAL_BRANCH == "free_surface" else interface
            simple_rows.append([float(value), float(selected["critical"])])

            row.update(
                {
                    "critical_Fr_free_surface": free["critical"],
                    "status_free_surface": free["status"],
                    "critical_Fr_interfacial": interface["critical"],
                    "status_interfacial": interface["status"],
                    "selected_branch": CRITICAL_BRANCH,
                    "selected_critical_Fr": selected["critical"],
                    "d_selected_at_Fr_min": selected["d_min"],
                    "d_selected_at_Fr_max": selected["d_max"],
                    "selected_bracket_left": selected["bracket_left"],
                    "selected_bracket_right": selected["bracket_right"],
                }
            )
        except Exception as exc:  # continue the sweep and record the failure
            simple_rows.append([float(value), np.nan])
            row.update(
                {
                    "critical_Fr_free_surface": np.nan,
                    "status_free_surface": "calculation_failed",
                    "critical_Fr_interfacial": np.nan,
                    "status_interfacial": "calculation_failed",
                    "selected_branch": CRITICAL_BRANCH,
                    "selected_critical_Fr": np.nan,
                    "d_selected_at_Fr_min": np.nan,
                    "d_selected_at_Fr_max": np.nan,
                    "selected_bracket_left": np.nan,
                    "selected_bracket_right": np.nan,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
        detailed_rows.append(row)

    simple_path = output_dir / f"critical_froude_vs_{symbol}.txt"
    header_lines = [
        "Two-layer power-law model: infinitely-long-wave critical Froude sweep",
        f"swept_parameter = {canonical} ({symbol})",
        f"selected_branch = {CRITICAL_BRANCH}",
        "critical condition = d_branch(Fr_l)=0 with a verified sign change",
        f"Froude search interval = [{FROUDE_SEARCH_MIN}, {FROUDE_SEARCH_MAX}]",
        "NaN means no positive sign-changing root was found in the search interval",
        f"columns: {symbol}    critical_Fr_{CRITICAL_BRANCH}",
    ]
    np.savetxt(
        simple_path,
        np.asarray(simple_rows, dtype=float),
        delimiter="\t",
        fmt="%.16e",
        header="\n".join(header_lines),
        comments="# ",
    )

    detailed_path = output_dir / f"critical_froude_vs_{symbol}_diagnostics.txt"
    fieldnames = list(detailed_rows[0].keys())
    with detailed_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "# Detailed branch diagnostics. The two-column file is the main plotting table.\n"
        )
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(detailed_rows)

    print(f"Wrote {simple_path}")
    print(f"Wrote {detailed_path}")


if __name__ == "__main__":
    main()
