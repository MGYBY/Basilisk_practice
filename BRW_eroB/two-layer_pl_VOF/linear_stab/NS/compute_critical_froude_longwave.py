#!/usr/bin/env python3
"""Compute branch-specific infinite-wavelength critical Froude numbers."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import argparse
import json

from fullns_thresholds import (
    build_normal_flow_shape,
    critical_froude_infinite_wavelength,
    write_long_wave_thresholds,
)
from stability_inputs import FAMILY, LONG_WAVE, LONG_WAVE_BRANCHES, SOLVER


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    shape = build_normal_flow_shape(FAMILY, SOLVER.base_samples)
    all_thresholds = []
    summary: dict[str, object] = {
        "family": asdict(FAMILY),
        "solver": asdict(SOLVER),
        "long_wave_options": asdict(LONG_WAVE),
        "compatibility_ratio_chi_lower": shape.compatibility_ratio,
        "lower_mean_velocity_check": shape.lower_mean_velocity,
        "branches": {},
    }

    for branch in LONG_WAVE_BRANCHES:
        thresholds = critical_froude_infinite_wavelength(
            shape,
            branch,
            SOLVER,
            LONG_WAVE,
        )
        all_thresholds.extend(thresholds)
        summary["branches"][branch] = [asdict(item) for item in thresholds]

        reliable = [item for item in thresholds if item.converged]
        if reliable:
            print(f"{branch}: reliable long-wave critical Froude number(s)")
            for item in reliable:
                print(
                    f"  Fr_l,c = {item.critical_froude:.10f}  "
                    f"(fit residual {item.maximum_fit_residual:.3e})"
                )
        elif thresholds:
            print(f"{branch}: no converged positive long-wave threshold")
            for item in thresholds:
                print(f"  diagnostic intercept {item.critical_froude:.10f}: {item.quality_message}")
        else:
            print(f"{branch}: no sign-changing finite-k roots in the requested range")

    write_long_wave_thresholds(
        all_thresholds,
        str(output / "critical_froude_longwave.tsv"),
    )
    with (output / "critical_froude_longwave.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"\nWrote results to {output}")


if __name__ == "__main__":
    main()
