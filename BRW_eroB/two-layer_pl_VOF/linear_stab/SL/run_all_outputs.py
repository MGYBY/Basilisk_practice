#!/usr/bin/env python3
"""Run the parameter sweep and the finite-wavenumber tables in sequence."""

from run_dispersion_tables import main as run_dispersion
from run_parameter_sweep import main as run_sweep


def main() -> None:
    run_sweep()
    run_dispersion()


if __name__ == "__main__":
    main()
