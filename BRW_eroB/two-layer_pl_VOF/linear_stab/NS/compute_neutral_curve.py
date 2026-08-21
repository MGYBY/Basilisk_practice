#!/usr/bin/env python3
"""Compute finite-wavelength neutral curves in Fr-Re and Fr-S0 coordinates."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import argparse
import json

from fullns_thresholds import (
    build_normal_flow_shape,
    finite_wavelength_neutral_curve,
    write_neutral_curve,
)
from stability_inputs import (
    FAMILY,
    NEUTRAL_BRANCHES,
    NEUTRAL_CURVE,
    REYNOLDS_VALUES,
    NEUTRAL_SOLVER,
    SOLVER,
    WAVE_SPECIFICATIONS,
)


def wave_tag(specification) -> str:
    text = f"{specification.value:g}".replace(".", "p")
    return f"{specification.kind}_{text}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="recompute tables even when the target file already exists",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    shape = build_normal_flow_shape(FAMILY, NEUTRAL_SOLVER.base_samples)
    manifest: dict[str, object] = {
        "family": asdict(FAMILY),
        "solver": asdict(NEUTRAL_SOLVER),
        "high_order_check_solver": asdict(SOLVER),
        "neutral_curve_options": asdict(NEUTRAL_CURVE),
        "reynolds_values": [float(value) for value in REYNOLDS_VALUES],
        "compatibility_ratio_chi_lower": shape.compatibility_ratio,
        "calculations": [],
    }

    for specification in WAVE_SPECIFICATIONS:
        for branch in NEUTRAL_BRANCHES:
            tag = f"{branch}__{wave_tag(specification)}"
            path = output / f"neutral_curve__{tag}.tsv"
            if path.exists() and not args.force:
                print(f"Using existing {path.name}", flush=True)
                points = []
                with path.open("r", encoding="utf-8") as handle:
                    point_count = sum(1 for line in handle if line.strip() and not line.startswith("#"))
            else:
                print(
                    f"Computing {branch}, {specification.kind}="
                    f"{specification.value:g} ...",
                    flush=True,
                )
                points = finite_wavelength_neutral_curve(
                    shape,
                    REYNOLDS_VALUES,
                    specification,
                    branch,
                    NEUTRAL_SOLVER,
                    NEUTRAL_CURVE,
                )
                write_neutral_curve(points, str(path))
                point_count = len(points)
            manifest["calculations"].append(
                {
                    "branch": branch,
                    "wave_specification": asdict(specification),
                    "point_count": point_count,
                    "file": path.name,
                }
            )
            print(f"  table contains {point_count} neutral points: {path.name}")

    with (output / "neutral_curve_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    print(f"\nWrote results to {output}")


if __name__ == "__main__":
    main()
