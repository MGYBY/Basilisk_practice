#!/usr/bin/env python3
"""Small regularization-sensitivity study for the representative case.

This is deliberately separate from the main runner because the smoothing shear
rate is a constitutive/numerical modelling choice, not a spectral resolution
parameter.  The table helps identify conclusions that are robust and those
that are near-neutral and regularization-sensitive.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import numpy as np

from fullns_stability import (
    CaseParameters,
    SolverOptions,
    build_base_flow,
    neutral_crossings,
    track_hydrodynamic_modes,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "regularization_sensitivity.tsv"
RATES = (0.25, 0.5, 1.0)  # dimensional smoothing rates [s^-1]
WAVENUMBERS = np.unique(
    np.concatenate(
        [
            np.geomspace(5.0e-4, 5.0e-2, 24),
            np.linspace(5.0e-2, 3.5e-1, 30),
        ]
    )
)
OPTIONS = SolverOptions(
    chebyshev_order=36,
    residual_tolerance=1.0e-5,
    seed_wavenumber=2.0e-4,
    base_samples=6001,
)


def metrics(k: np.ndarray, modes) -> tuple[float, float, float, float]:
    growth = np.asarray([mode.growth_rate for mode in modes])
    celerity = np.asarray([mode.celerity for mode in modes])
    index = int(np.argmax(growth))
    crossings = neutral_crossings(k, growth)
    first_crossing = crossings[0] if crossings else float("nan")
    return (
        float(growth[index]),
        float(k[index]),
        float(celerity[index]),
        float(first_crossing),
    )


def main() -> None:
    rows = []
    for rate in RATES:
        case = replace(CaseParameters(), regularization_rate=rate)
        base = build_base_flow(case, samples=OPTIONS.base_samples)
        tracks, _ = track_hydrodynamic_modes(base, WAVENUMBERS, OPTIONS)
        for branch in ("interface_varicose", "surface_zigzag"):
            maximum, k_max, c_max, k_neutral = metrics(
                WAVENUMBERS,
                tracks[branch],
            )
            rows.append(
                (
                    rate,
                    base.epsilon_lower,
                    branch,
                    maximum,
                    k_max,
                    c_max,
                    k_neutral,
                )
            )

    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(
            "# regularization_rate_per_s\tepsilon\tbranch\tmaximum_growth\t"
            "k_at_maximum\tcelerity_at_maximum\tfirst_neutral_k\n"
        )
        for row in rows:
            handle.write(
                f"{row[0]:.14e}\t{row[1]:.14e}\t{row[2]}\t"
                f"{row[3]:.14e}\t{row[4]:.14e}\t{row[5]:.14e}\t"
                f"{row[6]:.14e}\n"
            )
    print(OUTPUT)


if __name__ == "__main__":
    main()
