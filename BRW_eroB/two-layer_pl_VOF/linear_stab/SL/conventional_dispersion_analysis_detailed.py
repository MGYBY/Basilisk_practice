#!/usr/bin/env python3
"""Conventional dispersion-relation analysis of the two-layer power-law model.

The script follows the same sequence as the revised report:

1. evaluate the scalar coefficients a_ij and c_ij;
2. form the quartic dispersion polynomial for a prescribed real wavenumber k;
3. calculate its four temporal roots sigma_j(k) with ``numpy.roots``;
4. obtain growth rate Re(sigma_j) and phase celerity -Im(sigma_j)/k;
5. expand the same quartic at k -> 0+ to obtain the two long-wave
   celerities and the two k**2 growth coefficients;
6. calculate branch-specific critical Froude numbers directly from d_j(Fr)=0
   by scanning for a sign change and refining the root by bisection.

The main root calculation does not depend on two arbitrarily selected Froude
numbers.  The optional affine-law calculation is retained only as a transparent
verification of d_j=A_j*Fr**2+B_j for the present normalization.

There is deliberately no equilibrium Jacobian, projected relaxation system, or
left/right-eigenvector calculation in this file.  The only imported model task is
the evaluation of the scalar linearized coefficients.
"""

from __future__ import annotations

import itertools
from dataclasses import replace
from typing import Iterable

import numpy as np

from two_layer_powerlaw_stability_core import (
    DEFAULT_PARAMETERS,
    LinearCoefficients,
    Parameters,
    linear_coefficients,
)


# =============================================================================
# USER INPUTS
# =============================================================================

PARAMETERS = replace(
    DEFAULT_PARAMETERS,
    froude=0.8,
    depth_ratio=1.0,
    density_ratio=0.9,
    n_lower=0.6,
    n_upper=0.8,
    scaled_consistency_ratio=1.0,
)

# Wavenumbers at which detailed finite-k roots are printed.
PRINT_WAVENUMBERS = (1.0e-3, 1.0e-2, 1.0e-1, 1.0)


# =============================================================================
# 1. Quartic dispersion relation
# =============================================================================


def quartic_polynomial(k: float, a: LinearCoefficients) -> np.ndarray:
    """Return [1, chi3, chi2, chi1, chi0] for real wavenumber ``k``.

    The polynomial is

        sigma**4 + chi3*sigma**3 + chi2*sigma**2
                   + chi1*sigma + chi0 = 0,

    with z = i*k.  The formulas are the scalar expansion of
    P11*P22 - P12*P21 = 0 in the report.
    """
    if not np.isfinite(k) or k < 0.0:
        raise ValueError("k must be finite and non-negative")

    z = 1j * float(k)

    chi3 = z * (a.a33 + a.a44) - (a.c33 + a.c44)

    chi2 = (
        a.c33 * a.c44
        - a.c34 * a.c43
        + z
        * (
            a.c31
            + a.c42
            - a.a33 * a.c44
            - a.a44 * a.c33
            + a.a34 * a.c43
            + a.a43 * a.c34
        )
        + z**2 * (a.a33 * a.a44 - a.a34 * a.a43 - a.a31 - a.a42)
    )

    chi1 = (
        z
        * (
            -a.c31 * a.c44
            + a.c32 * a.c43
            - a.c33 * a.c42
            + a.c34 * a.c41
        )
        + z**2
        * (
            a.a31 * a.c44
            - a.a32 * a.c43
            + a.a33 * a.c42
            - a.a34 * a.c41
            - a.a41 * a.c34
            + a.a42 * a.c33
            - a.a43 * a.c32
            + a.a44 * a.c31
        )
        + z**3
        * (
            -a.a31 * a.a44
            + a.a32 * a.a43
            - a.a33 * a.a42
            + a.a34 * a.a41
        )
    )

    chi0 = (
        z**2 * (a.c31 * a.c42 - a.c32 * a.c41)
        + z**3
        * (
            -a.a31 * a.c42
            + a.a32 * a.c41
            + a.a41 * a.c32
            - a.a42 * a.c31
        )
        + z**4 * (a.a31 * a.a42 - a.a32 * a.a41)
    )

    return np.array([1.0 + 0.0j, chi3, chi2, chi1, chi0], dtype=complex)


def dispersion_roots(k: float, parameters: Parameters) -> np.ndarray:
    """Calculate the four roots sigma_j(k) directly from the quartic."""
    return np.roots(quartic_polynomial(k, linear_coefficients(parameters)))


def depth_amplitude_ratio(
    sigma: complex,
    k: float,
    a: LinearCoefficients,
) -> complex:
    """Return h_upper_hat / h_lower_hat for one dispersion root.

    The ratio is recovered from the two scalar depth-amplitude equations.  The
    better-conditioned of the two rows is used automatically.
    """
    if k <= 0.0:
        raise ValueError("The finite-k amplitude ratio requires k > 0")

    z = 1j * float(k)
    p11 = sigma**2 + sigma * (z * a.a33 - a.c33) + z * a.c31 - z**2 * a.a31
    p12 = sigma * (z * a.a34 - a.c34) + z * a.c32 - z**2 * a.a32
    p21 = sigma * (z * a.a43 - a.c43) + z * a.c41 - z**2 * a.a41
    p22 = sigma**2 + sigma * (z * a.a44 - a.c44) + z * a.c42 - z**2 * a.a42

    if abs(p12) >= abs(p22) and abs(p12) > 1.0e-14:
        return -p11 / p12
    if abs(p22) > 1.0e-14:
        return -p21 / p22
    raise FloatingPointError("Both scalar amplitude rows are singular or ill-conditioned")


# =============================================================================
# 2. Direct k -> 0+ expansion of the same quartic
# =============================================================================


def small_k_quartic_coefficients(a: LinearCoefficients) -> dict[str, float]:
    """Return only the quartic coefficients required through O(k**2 growth)."""
    return {
        "chi30": -(a.c33 + a.c44),
        "chi20": a.c33 * a.c44 - a.c34 * a.c43,
        "chi21": (
            a.c31
            + a.c42
            - a.a33 * a.c44
            - a.a44 * a.c33
            + a.a34 * a.c43
            + a.a43 * a.c34
        ),
        "chi11": (
            -a.c31 * a.c44
            + a.c32 * a.c43
            - a.c33 * a.c42
            + a.c34 * a.c41
        ),
        "chi12": (
            a.a31 * a.c44
            - a.a32 * a.c43
            + a.a33 * a.c42
            - a.a34 * a.c41
            - a.a41 * a.c34
            + a.a42 * a.c33
            - a.a43 * a.c32
            + a.a44 * a.c31
        ),
        "chi02": a.c31 * a.c42 - a.c32 * a.c41,
        "chi03": (
            -a.a31 * a.c42
            + a.a32 * a.c41
            + a.a41 * a.c32
            - a.a42 * a.c31
        ),
    }


def classify_mode(c: float, a: LinearCoefficients) -> tuple[str, float, float]:
    """Classify a long-wave root and return normalized depth amplitudes.

    The upper-layer amplitude is set to one.  A mode with a small total-depth
    displacement relative to the internal-interface displacement is labelled
    interfacial; the other is labelled free_surface.
    """
    denominator = c * a.c34 + a.c32
    if abs(denominator) < 1.0e-12:
        raise FloatingPointError("Long-wave amplitude ratio denominator is too small")

    ratio_upper_over_lower = -(c * a.c33 + a.c31) / denominator
    h_upper = 1.0
    h_lower = 1.0 / ratio_upper_over_lower

    free_surface_displacement = h_lower + h_upper
    label = (
        "interfacial"
        if abs(free_surface_displacement) < abs(h_lower)
        else "free_surface"
    )
    return label, float(h_lower), float(h_upper)


def longwave_modes(parameters: Parameters) -> list[dict[str, float | str]]:
    """Return c, d and depth amplitudes directly from the quartic expansion."""
    a = linear_coefficients(parameters)
    chi = small_k_quartic_coefficients(a)

    c_roots = np.roots([chi["chi20"], -chi["chi11"], chi["chi02"]])
    if np.max(np.abs(c_roots.imag)) > 1.0e-9:
        raise FloatingPointError(f"Long-wave celerities are not real: {c_roots}")

    modes: list[dict[str, float | str]] = []
    for c_complex in c_roots:
        c = float(c_complex.real)
        denominator = 2.0 * chi["chi20"] * c - chi["chi11"]
        if abs(denominator) < 1.0e-11:
            raise FloatingPointError("Degenerate long-wave root: d denominator is too small")

        d = (
            chi["chi30"] * c**3
            - chi["chi21"] * c**2
            + chi["chi12"] * c
            - chi["chi03"]
        ) / denominator

        label, h_lower, h_upper = classify_mode(c, a)
        modes.append(
            {
                "label": label,
                "celerity": c,
                "growth_coefficient": float(d),
                "h_lower_amplitude": h_lower,
                "h_upper_amplitude": h_upper,
            }
        )

    return sorted(modes, key=lambda item: 0 if item["label"] == "interfacial" else 1)


def affine_critical_froude_numbers(parameters: Parameters) -> list[dict[str, float | str]]:
    """Optional affine-law check for d=A*Fr**2+B.

    ``sample_froude_a`` and ``sample_froude_b`` below are simply two distinct
    evaluation points used to reconstruct a straight line in Fr**2.  They are
    not the Froude numbers of the two layers and they are not critical values.
    The direct critical-Froude calculation is implemented separately below.
    """
    sample_froude_a, sample_froude_b = 0.5, 1.0
    modes1 = {
        str(mode["label"]): mode
        for mode in longwave_modes(replace(parameters, froude=sample_froude_a))
    }
    modes2 = {
        str(mode["label"]): mode
        for mode in longwave_modes(replace(parameters, froude=sample_froude_b))
    }

    results: list[dict[str, float | str]] = []
    for label in ("interfacial", "free_surface"):
        d1 = float(modes1[label]["growth_coefficient"])
        d2 = float(modes2[label]["growth_coefficient"])
        slope = (d2 - d1) / (sample_froude_b**2 - sample_froude_a**2)
        intercept = d1 - slope * sample_froude_a**2

        critical = np.nan
        if slope * intercept < 0.0:
            critical = float(np.sqrt(-intercept / slope))

        if np.isfinite(critical):
            below = intercept
            above = slope * (1.25 * critical) ** 2 + intercept
            interpretation = (
                "stable below and unstable above"
                if below < 0.0 < above
                else "unstable below and stable above"
            )
        elif slope > 0.0 and intercept > 0.0:
            interpretation = "unstable for every positive Froude number"
        elif slope < 0.0 and intercept < 0.0:
            interpretation = "stable for every positive Froude number"
        else:
            interpretation = "no positive sign-changing root in the affine family"

        results.append(
            {
                "label": label,
                "A": float(slope),
                "B": float(intercept),
                "critical_froude": critical,
                "interpretation": interpretation,
            }
        )

    return results



# =============================================================================
# 3. Direct critical-Froude calculation from d_j(Fr)=0
# =============================================================================


def branch_growth_coefficient(
    froude: float,
    parameters: Parameters,
    branch: str,
) -> float:
    """Return the long-wave coefficient d_j for one named branch."""
    if branch not in {"interfacial", "free_surface"}:
        raise ValueError("branch must be 'interfacial' or 'free_surface'")
    if not np.isfinite(froude) or froude <= 0.0:
        raise ValueError("froude must be positive and finite")
    modes = {
        str(mode["label"]): mode
        for mode in longwave_modes(replace(parameters, froude=float(froude)))
    }
    return float(modes[branch]["growth_coefficient"])


def find_sign_change(
    parameters: Parameters,
    branch: str,
    froude_min: float,
    froude_max: float,
    scan_points: int = 300,
) -> tuple[float, float] | None:
    """Return the first interval in which d_j changes sign."""
    if not (0.0 < froude_min < froude_max):
        raise ValueError("Require 0 < froude_min < froude_max")
    if scan_points < 2:
        raise ValueError("scan_points must be at least two")

    values = np.linspace(froude_min, froude_max, scan_points)
    left = float(values[0])
    d_left = branch_growth_coefficient(left, parameters, branch)

    for value in values[1:]:
        right = float(value)
        d_right = branch_growth_coefficient(right, parameters, branch)
        if d_right == 0.0:
            return right, right
        if np.signbit(d_left) != np.signbit(d_right):
            return left, right
        left, d_left = right, d_right
    return None


def bisect_critical_froude(
    parameters: Parameters,
    branch: str,
    bracket: tuple[float, float],
    tolerance: float = 1.0e-12,
    maximum_iterations: int = 200,
) -> float:
    """Refine a sign-changing bracket for d_j(Fr)=0 by bisection."""
    left, right = map(float, bracket)
    if left == right:
        return left
    d_left = branch_growth_coefficient(left, parameters, branch)
    d_right = branch_growth_coefficient(right, parameters, branch)
    if np.signbit(d_left) == np.signbit(d_right):
        raise ValueError("The supplied bracket does not contain a sign change")

    for _ in range(maximum_iterations):
        middle = 0.5 * (left + right)
        d_middle = branch_growth_coefficient(middle, parameters, branch)
        if abs(d_middle) < tolerance or 0.5 * (right - left) < tolerance:
            return middle
        if np.signbit(d_left) != np.signbit(d_middle):
            right, d_right = middle, d_middle
        else:
            left, d_left = middle, d_middle
    raise RuntimeError("Bisection did not converge")


def direct_critical_froude(
    parameters: Parameters,
    branch: str,
    froude_min: float = 0.05,
    froude_max: float = 1.20,
    scan_points: int = 300,
) -> dict[str, object]:
    """Locate the first positive sign-changing root of d_j(Fr)=0."""
    bracket = find_sign_change(
        parameters,
        branch,
        froude_min,
        froude_max,
        scan_points,
    )
    d_min = branch_growth_coefficient(froude_min, parameters, branch)
    d_max = branch_growth_coefficient(froude_max, parameters, branch)

    if bracket is None:
        if d_min > 0.0 and d_max > 0.0:
            interpretation = (
                "positive at both interval endpoints; no sign-changing root "
                "was resolved by the scan"
            )
        elif d_min < 0.0 and d_max < 0.0:
            interpretation = (
                "negative at both interval endpoints; no sign-changing root "
                "was resolved by the scan"
            )
        else:
            interpretation = "no sign-changing root resolved in the interval"
        return {
            "branch": branch,
            "critical_froude": None,
            "bracket": None,
            "d_at_min": d_min,
            "d_at_max": d_max,
            "interpretation": interpretation,
        }

    critical = bisect_critical_froude(parameters, branch, bracket)
    delta = max(1.0e-6, 1.0e-4 * critical)
    d_below = branch_growth_coefficient(max(froude_min, critical - delta), parameters, branch)
    d_above = branch_growth_coefficient(min(froude_max, critical + delta), parameters, branch)
    if d_below < 0.0 < d_above:
        interpretation = "stable below and unstable above"
    elif d_below > 0.0 > d_above:
        interpretation = "unstable below and stable above"
    else:
        interpretation = "root found; inspect the local sign change"
    return {
        "branch": branch,
        "critical_froude": critical,
        "bracket": bracket,
        "d_below": d_below,
        "d_above": d_above,
        "interpretation": interpretation,
    }


# =============================================================================
# 4. Optional finite-k branch tracking for readable tables or plots
# =============================================================================


def _best_permutation(reference: np.ndarray, candidates: np.ndarray) -> np.ndarray:
    best_cost = np.inf
    best_order: tuple[int, ...] | None = None
    for order in itertools.permutations(range(4)):
        ordered = candidates[list(order)]
        cost = float(np.sum(np.abs(ordered - reference)))
        if cost < best_cost:
            best_cost = cost
            best_order = order
    assert best_order is not None
    return candidates[list(best_order)]


def tracked_roots(
    wavenumbers: Iterable[float],
    parameters: Parameters,
) -> dict[str, np.ndarray]:
    """Track the four quartic roots by continuity in the complex plane."""
    k_values = np.asarray(tuple(wavenumbers), dtype=float)
    if k_values.ndim != 1 or k_values.size == 0 or np.any(k_values <= 0.0):
        raise ValueError("wavenumbers must be a nonempty positive 1-D sequence")
    if np.any(np.diff(k_values) <= 0.0):
        raise ValueError("wavenumbers must be strictly increasing")

    long = {str(mode["label"]): mode for mode in longwave_modes(parameters)}
    a = linear_coefficients(parameters)
    relaxation = np.roots(
        [1.0, -(a.c33 + a.c44), a.c33 * a.c44 - a.c34 * a.c43]
    )
    relaxation = relaxation[np.argsort(relaxation.real)]

    k0 = float(k_values[0])
    expected = np.array(
        [
            -1j * float(long["interfacial"]["celerity"]) * k0
            + float(long["interfacial"]["growth_coefficient"]) * k0**2,
            -1j * float(long["free_surface"]["celerity"]) * k0
            + float(long["free_surface"]["growth_coefficient"]) * k0**2,
            relaxation[0],
            relaxation[1],
        ],
        dtype=complex,
    )

    roots = np.empty((k_values.size, 4), dtype=complex)
    roots[0] = _best_permutation(expected, dispersion_roots(k0, parameters))
    for index, k in enumerate(k_values[1:], start=1):
        roots[index] = _best_permutation(roots[index - 1], dispersion_roots(float(k), parameters))

    return {
        "wavenumber": k_values,
        "interfacial": roots[:, 0],
        "free_surface": roots[:, 1],
        "relaxation_fast": roots[:, 2],
        "relaxation_slow": roots[:, 3],
    }


# =============================================================================
# 4. Human-readable demonstration
# =============================================================================


def main() -> None:
    parameters = PARAMETERS.validated()
    print("Inputs")
    print(parameters)

    print("\nLong-wave modes from the quartic")
    for mode in longwave_modes(parameters):
        print(
            f"  {mode['label']:12s}: "
            f"c={float(mode['celerity']):.12g}, "
            f"d={float(mode['growth_coefficient']):.12g}, "
            f"(h_l,h_u)=({float(mode['h_lower_amplitude']):.9g}, "
            f"{float(mode['h_upper_amplitude']):.9g})"
        )

    print("\nDirect infinite-wavelength critical Froude numbers")
    for branch in ("interfacial", "free_surface"):
        result = direct_critical_froude(parameters, branch)
        critical = result["critical_froude"]
        critical_text = f"{critical:.12g}" if isinstance(critical, float) else "none"
        print(
            f"  {branch:12s}: Fr_c={critical_text}; "
            f"{result['interpretation']}"
        )
        if result["bracket"] is not None:
            print(f"      initial sign-change bracket = {result['bracket']}")

    print("\nOptional affine-law verification")
    for result in affine_critical_froude_numbers(parameters):
        critical = float(result["critical_froude"])
        critical_text = f"{critical:.12g}" if np.isfinite(critical) else "none"
        print(
            f"  {result['label']:12s}: "
            f"d = ({float(result['A']):.12g}) Fr_l^2 "
            f"+ ({float(result['B']):.12g}); "
            f"Fr_c={critical_text}; {result['interpretation']}"
        )

    coefficients = linear_coefficients(parameters)
    print("\nSelected finite-wavenumber roots")
    for k in PRINT_WAVENUMBERS:
        roots = dispersion_roots(float(k), parameters)
        roots = roots[np.argsort(-roots.imag)]
        print(f"\n  k={k:g}")
        for index, sigma in enumerate(roots, start=1):
            growth = sigma.real
            celerity = -sigma.imag / k
            print(
                f"    root {index}: sigma={sigma.real:+.8e}{sigma.imag:+.8e}j, "
                f"growth={growth:+.8e}, c={celerity:+.8e}"
            )

    # Independent consistency check: the direct quartic roots must agree with
    # the eigenvalues of the original four scalar amplitude system.
    max_error = 0.0
    for k in PRINT_WAVENUMBERS:
        quartic = dispersion_roots(float(k), parameters)
        matrix = np.linalg.eigvals(coefficients.C() - 1j * k * coefficients.A())
        matched = _best_permutation(quartic, matrix)
        max_error = max(max_error, float(np.max(np.abs(quartic - matched))))
    print(f"\nMaximum quartic-versus-four-equation root discrepancy: {max_error:.3e}")


if __name__ == "__main__":
    main()
