#!/usr/bin/env python3
"""Small regression test for the practical output workflow."""

from __future__ import annotations

import numpy as np

from branch_tracking import BRANCHES, track_roots
from conventional_dispersion_analysis_detailed import (
    direct_critical_froude,
    dispersion_roots,
)
from stability_user_inputs import BASE_PARAMETERS


def main() -> None:
    critical = direct_critical_froude(
        BASE_PARAMETERS,
        "free_surface",
        0.05,
        1.2,
    )["critical_froude"]
    if critical is None or abs(float(critical) - 0.3822992557) > 5.0e-9:
        raise AssertionError(f"Unexpected default critical Froude number: {critical}")

    k = np.geomspace(1.0e-4, 10.0, 301)
    tracked = track_roots(k, BASE_PARAMETERS)

    for index, value in enumerate(k):
        direct = dispersion_roots(float(value), BASE_PARAMETERS)
        labelled = np.array([tracked[name][index] for name in BRANCHES])
        # Compare the two root sets without assuming their order.
        remaining = list(direct)
        errors = []
        for root in labelled:
            position = int(np.argmin([abs(root - item) for item in remaining]))
            errors.append(abs(root - remaining.pop(position)))
        if max(errors) > 1.0e-9:
            raise AssertionError(
                f"Tracked roots do not reproduce the quartic at k={value}: {max(errors)}"
            )

    print("All practical-workflow checks passed.")
    print(f"Default free-surface critical Fr = {float(critical):.12f}")
    print(
        "Maximum normalized branch-matching cost = "
        f"{float(tracked['maximum_normalized_matching_cost']):.6e}"
    )


if __name__ == "__main__":
    main()
