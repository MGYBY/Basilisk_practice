#!/usr/bin/env python3
"""Verification checks for compatibility and neutral-threshold calculations."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import argparse
import math

import numpy as np

from fullns_stability import CaseParameters, SolverOptions, build_base_flow
from fullns_thresholds import (
    DimensionlessFamily,
    LongWaveOptions,
    NeutralCurveOptions,
    WaveSpecification,
    build_dimensionless_base_flow,
    build_normal_flow_shape,
    compatibility_reynolds,
    compatibility_slope,
    critical_froude_infinite_wavelength,
    family_from_base_flow,
    finite_wavelength_neutral_points_at_reynolds,
    wavenumber_from_specification,
    wavelength_measures,
)


RESULTS: list[tuple[str, bool, float, float]] = []


def check(name: str, error: float, tolerance: float) -> None:
    passed = bool(np.isfinite(error) and error <= tolerance)
    RESULTS.append((name, passed, float(error), float(tolerance)))
    if not passed:
        raise AssertionError(
            f"{name}: error={error:.6e} > tolerance={tolerance:.6e}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "verification",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    # 1. Reconstruct the dimensionless family of the audited dimensional case.
    dimensional = build_base_flow(CaseParameters(), samples=12001)
    family = family_from_base_flow(dimensional)
    shape = build_normal_flow_shape(family, samples=12001)
    dimensional_chi = (
        dimensional.parameters.slope_tan
        / (dimensional.froude**2 * dimensional.consistency_number_lower)
    )
    check(
        "Compatibility ratio versus dimensional normal flow",
        abs(shape.compatibility_ratio - dimensional_chi),
        3.0e-8,
    )
    check(
        "Lower mean-velocity compatibility",
        abs(shape.lower_mean_velocity - 1.0),
        1.0e-7,
    )
    check(
        "Dimensionless lower velocity shape",
        float(np.max(np.abs(shape.lower_velocity - dimensional.lower.velocity))),
        2.0e-7,
    )
    check(
        "Dimensionless upper velocity shape",
        float(np.max(np.abs(shape.upper_velocity - dimensional.upper.velocity))),
        2.0e-7,
    )

    # 2. Newtonian one-layer limit: chi_l -> 3, so Fr^2 = Re*S0/3.
    almost_one_layer = DimensionlessFamily(
        depth_ratio=1.0e-8,
        density_ratio=1.0,
        n_lower=1.0,
        n_upper=1.0,
        consistency_ratio=1.0,
        epsilon_lower=0.0,
        epsilon_upper=0.0,
    )
    one_layer_shape = build_normal_flow_shape(almost_one_layer, samples=12001)
    check(
        "Newtonian one-layer compatibility chi_l=3",
        abs(one_layer_shape.compatibility_ratio - 3.0),
        2.0e-7,
    )

    # 3. Round-trip transformation among Fr, Re, and S0.
    fr_test, re_test = 0.73, 18.0
    slope = compatibility_slope(shape, fr_test, re_test)
    check(
        "Compatibility Fr-Re-S0 round trip",
        abs(compatibility_reynolds(shape, fr_test, slope) - re_test),
        2.0e-12,
    )
    rebuilt = build_dimensionless_base_flow(shape, fr_test, re_test)
    check(
        "Base-flow slope reconstructed from compatibility",
        abs(rebuilt.parameters.slope_tan - slope),
        2.0e-14,
    )

    # 4. Figure-2 wavelength convention is reproduced identically.
    specification = WaveSpecification("fixed_slope_scaled_wavelength", 4.72)
    kH = wavenumber_from_specification(specification, shape, fr_test, re_test)
    lambda_H, slope_lambda_H = wavelength_measures(kH, slope)
    check(
        "Slope-scaled wavelength mapping",
        abs(slope_lambda_H - 4.72),
        2.0e-12,
    )
    check(
        "Positive depth-scaled wavelength",
        0.0 if lambda_H > 0.0 else 1.0,
        0.0,
    )

    # 5. Long-wave surface threshold converges with spectral order.
    small_k = (3.0e-4, 5.0e-4, 8.0e-4, 1.2e-3)
    lw = LongWaveOptions(
        froude_min=0.298,
        froude_max=0.323,
        scan_points=9,
        reference_reynolds=20.0,
        wavenumbers=small_k,
        polynomial_degree=1,
    )
    roots: list[float] = []
    for order in (36, 44):
        solver = SolverOptions(
            chebyshev_order=order,
            residual_tolerance=4.0e-5,
            seed_wavenumber=1.0e-4,
            base_samples=12001,
        )
        thresholds = critical_froude_infinite_wavelength(
            shape,
            "surface_zigzag",
            solver,
            lw,
        )
        reliable = [threshold for threshold in thresholds if threshold.converged]
        if not reliable:
            raise AssertionError(f"No converged surface threshold at order {order}")
        roots.append(reliable[0].critical_froude)
    check(
        "Surface critical Froude spectral convergence",
        abs(roots[0] - roots[1]),
        8.0e-5,
    )
    check(
        "Surface critical Froude benchmark",
        abs(roots[-1] - 0.31026),
        1.5e-4,
    )

    # 6. k->0 threshold is independent of the reference Re used for the
    # extrapolation, within the numerical resolution of the small-k pencil.
    re_roots: list[float] = []
    solver36 = SolverOptions(
        chebyshev_order=36,
        residual_tolerance=4.0e-5,
        seed_wavenumber=1.0e-4,
        base_samples=12001,
    )
    for reference_re in (5.0, 40.0):
        thresholds = critical_froude_infinite_wavelength(
            shape,
            "surface_zigzag",
            solver36,
            replace(lw, reference_reynolds=reference_re),
        )
        reliable = [threshold for threshold in thresholds if threshold.converged]
        if not reliable:
            raise AssertionError("No reliable threshold in Reynolds-invariance check")
        re_roots.append(reliable[0].critical_froude)
    check(
        "Long-wave critical Froude Reynolds invariance",
        abs(re_roots[0] - re_roots[1]),
        1.5e-4,
    )

    # 7. The interface branch should not be reported as a converged positive
    # threshold when its finite-k roots do not approach a common intercept.
    interface_thresholds = critical_froude_infinite_wavelength(
        shape,
        "interface_varicose",
        solver36,
        LongWaveOptions(
            froude_min=0.10,
            froude_max=1.50,
            scan_points=11,
            reference_reynolds=20.0,
            wavenumbers=small_k,
            polynomial_degree=1,
        ),
    )
    reliable_interface = [item for item in interface_thresholds if item.converged]
    check(
        "No spurious converged interface threshold",
        float(len(reliable_interface)),
        0.0,
    )

    # 8. A finite-wavelength surface neutral point is reproduced at two
    # spectral orders and satisfies the generalized eigenproblem residual.
    curve_options = NeutralCurveOptions(
        froude_min=0.26,
        froude_max=0.38,
        scan_points=14,
        root_xtol=5.0e-7,
        continuation_steps=8,
    )
    finite_roots: list[float] = []
    residuals: list[float] = []
    for order in (32, 36):
        solver = SolverOptions(
            chebyshev_order=order,
            residual_tolerance=6.0e-5,
            seed_wavenumber=1.0e-4,
            base_samples=8001,
        )
        points = finite_wavelength_neutral_points_at_reynolds(
            shape,
            20.0,
            specification,
            "surface_zigzag",
            solver,
            curve_options,
        )
        if not points:
            raise AssertionError("No finite-wavelength surface neutral point found")
        finite_roots.append(points[0].froude)
        residuals.append(points[0].residual)
        check(
            f"Finite-wave neutral growth at order {order}",
            abs(points[0].growth_rate),
            2.0e-7,
        )
    check(
        "Finite-wave neutral Froude spectral convergence",
        abs(finite_roots[0] - finite_roots[1]),
        1.0e-3,
    )
    check(
        "Finite-wave original-pencil residual",
        max(residuals),
        1.0e-5,
    )

    path = output / "threshold_verification.txt"
    with path.open("w", encoding="utf-8") as handle:
        for name, passed, error, tolerance in RESULTS:
            handle.write(
                f"{'PASS' if passed else 'FAIL'}\t{name}\t"
                f"error={error:.12e}\ttolerance={tolerance:.12e}\n"
            )
        handle.write("\nAll threshold verification checks passed.\n")
        handle.write(f"Surface long-wave roots (N=36,44): {roots}\n")
        handle.write(f"Reference-Re roots (Re=5,40): {re_roots}\n")
        handle.write(f"Finite-wave roots (N=32,36): {finite_roots}\n")

    print(path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
