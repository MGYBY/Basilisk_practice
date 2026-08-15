#!/usr/bin/env python3
"""Robust continuation of the four roots of the quartic dispersion relation.

The four labels are defined at k -> 0+ and then continued in increasing k:

* free_surface: hydrodynamic mode with mainly in-phase layer-depth motion;
* interfacial: hydrodynamic mode with mainly out-of-phase layer-depth motion;
* momentum_fast: root connected to the more strongly damped source eigenvalue;
* momentum_slow: root connected to the less strongly damped source eigenvalue.

No new governing equation is introduced.  This module only keeps the roots in a
consistent order while the existing quartic is solved at successive wavenumbers.
"""

from __future__ import annotations

import itertools
from typing import Iterable

import numpy as np

from conventional_dispersion_analysis_detailed import (
    depth_amplitude_ratio,
    dispersion_roots,
    longwave_modes,
)
from two_layer_powerlaw_stability_core import Parameters, linear_coefficients


BRANCHES = (
    "free_surface",
    "interfacial",
    "momentum_fast",
    "momentum_slow",
)


def _source_relaxation_roots(parameters: Parameters) -> tuple[complex, complex]:
    coefficients = linear_coefficients(parameters)
    roots = np.roots(
        [
            1.0,
            -(coefficients.c33 + coefficients.c44),
            coefficients.c33 * coefficients.c44
            - coefficients.c34 * coefficients.c43,
        ]
    )
    roots = roots[np.argsort(roots.real)]
    return complex(roots[0]), complex(roots[1])


def _longwave_predictions(k: float, parameters: Parameters) -> np.ndarray:
    modes = {str(item["label"]): item for item in longwave_modes(parameters)}
    relaxation_fast, relaxation_slow = _source_relaxation_roots(parameters)

    def hydro(branch: str) -> complex:
        c = float(modes[branch]["celerity"])
        d = float(modes[branch]["growth_coefficient"])
        return -1j * c * k + d * k * k

    return np.array(
        [
            hydro("free_surface"),
            hydro("interfacial"),
            relaxation_fast,
            relaxation_slow,
        ],
        dtype=complex,
    )


def _match(predicted: np.ndarray, candidates: np.ndarray) -> tuple[np.ndarray, float]:
    """Match candidates to predictions by a normalized global assignment."""
    best_cost = np.inf
    best: np.ndarray | None = None
    scale = 1.0 + np.abs(predicted)

    for order in itertools.permutations(range(4)):
        ordered = candidates[list(order)]
        cost = float(np.sum(np.abs(ordered - predicted) / scale))
        if cost < best_cost:
            best_cost = cost
            best = ordered

    if best is None:  # pragma: no cover - defensive only
        raise RuntimeError("Could not match dispersion roots")
    return best, best_cost


def track_roots(
    wavenumbers: Iterable[float],
    parameters: Parameters,
) -> dict[str, np.ndarray | float]:
    """Return four consistently labelled root arrays.

    The first two points are anchored to the analytical k -> 0+ expansion.  At
    later points, a secant predictor is used before the four quartic roots are
    matched.  This avoids the label ambiguity that can occur if roots are sorted
    independently at every wavenumber or matched only to the immediately
    preceding value.
    """
    k_values = np.asarray(tuple(wavenumbers), dtype=float)
    if k_values.ndim != 1 or k_values.size < 2:
        raise ValueError("At least two wavenumbers are required")
    if np.any(~np.isfinite(k_values)) or np.any(k_values <= 0.0):
        raise ValueError("All wavenumbers must be finite and strictly positive")
    if np.any(np.diff(k_values) <= 0.0):
        raise ValueError("Wavenumbers must be strictly increasing")

    roots = np.empty((k_values.size, 4), dtype=complex)
    matching_costs = np.empty(k_values.size, dtype=float)

    for index, k in enumerate(k_values):
        candidates = dispersion_roots(float(k), parameters)

        if index < 2:
            predicted = _longwave_predictions(float(k), parameters)
        else:
            dk_previous = k_values[index - 1] - k_values[index - 2]
            dk_current = k_values[index] - k_values[index - 1]
            predicted = roots[index - 1] + (
                roots[index - 1] - roots[index - 2]
            ) * (dk_current / dk_previous)

        roots[index], matching_costs[index] = _match(predicted, candidates)

    return {
        "wavenumber": k_values,
        "free_surface": roots[:, 0],
        "interfacial": roots[:, 1],
        "momentum_fast": roots[:, 2],
        "momentum_slow": roots[:, 3],
        "maximum_normalized_matching_cost": float(np.max(matching_costs)),
    }


def hydrodynamic_amplitude_ratios(
    tracked: dict[str, np.ndarray | float],
    parameters: Parameters,
) -> dict[str, np.ndarray]:
    """Return h_u_hat/h_l_hat for the two hydrodynamic branches."""
    k_values = np.asarray(tracked["wavenumber"], dtype=float)
    coefficients = linear_coefficients(parameters)
    result: dict[str, np.ndarray] = {}

    for branch in ("free_surface", "interfacial"):
        sigma_values = np.asarray(tracked[branch], dtype=complex)
        ratios = np.empty_like(sigma_values)
        for index, (k, sigma) in enumerate(zip(k_values, sigma_values)):
            ratios[index] = depth_amplitude_ratio(
                complex(sigma), float(k), coefficients
            )
        result[branch] = ratios

    return result
