#!/usr/bin/env python3
"""Generate the long-wave and finite-wavenumber spectral-convergence tables."""
from __future__ import annotations

from pathlib import Path
import argparse

from fullns_stability import convergence_at_wavenumber
from run_stability_case import CASE, SOLVER, DEFAULT_OUTPUT, write_convergence_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    orders = (32, 40, 52, 64, 76)
    longwave = convergence_at_wavenumber(
        CASE,
        0.01,
        orders=orders,
        seed_wavenumber=SOLVER.seed_wavenumber,
    )
    finite = convergence_at_wavenumber(
        CASE,
        0.4,
        orders=orders,
        seed_wavenumber=SOLVER.seed_wavenumber,
    )
    write_convergence_table(longwave, output / "convergence_k0p01.tsv")
    write_convergence_table(finite, output / "convergence_k0p4.tsv")
    print(f"Wrote convergence tables to {output}")


if __name__ == "__main__":
    main()
