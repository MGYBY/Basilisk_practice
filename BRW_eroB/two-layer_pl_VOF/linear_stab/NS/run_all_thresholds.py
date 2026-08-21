#!/usr/bin/env python3
"""Run the long-wave threshold, finite-wave neutral curves, and plots."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
FIGURES = ROOT / "report" / "figures"


def run(script: str, *arguments: str) -> None:
    command = [sys.executable, str(HERE / script), *arguments]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    run("compute_critical_froude_longwave.py", "--output", str(DATA))
    run("compute_neutral_curve.py", "--output", str(DATA))
    run("make_threshold_figures.py", "--data", str(DATA), "--output", str(FIGURES))
    run("verify_thresholds.py", "--output", str(ROOT / "verification"))


if __name__ == "__main__":
    main()
