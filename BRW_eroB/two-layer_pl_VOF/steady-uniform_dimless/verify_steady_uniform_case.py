#!/usr/bin/env python3
"""Regression checks for the executed steady-uniform relaxation benchmark."""

from __future__ import annotations

import configparser
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "steady_uniform_case_parameters.ini"


def read_history(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, names=True, delimiter="\t", dtype=float)
    if data.shape == ():
        data = np.asarray([data], dtype=data.dtype)
    return data


def assert_close(value: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isfinite(value) or abs(value - expected) > tolerance:
        raise AssertionError(
            f"{label}: got {value:.16g}, expected {expected:.16g} "
            f"within {tolerance:.3g}"
        )


def main() -> None:
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not parser.read(CONFIG):
        raise FileNotFoundError(CONFIG)

    lx = parser.getfloat("domain_and_mesh", "lx")
    hr = parser.getfloat("physical_problem", "depth_ratio")
    levels = [
        parser.getint("domain_and_mesh", key)
        for key in ("min_level", "initial_level", "max_level")
    ]
    if len(set(levels)) != 1:
        raise AssertionError(f"Stationary benchmark requires equal levels, got {levels}")
    if not (0.5 <= lx / (1.0 + hr) <= 10.0):
        raise AssertionError("Streamwise domain is not of the same order as liquid depth")

    initial_u = parser.getfloat("steady_uniform", "initial_uniform_velocity")
    if not (0.0 < initial_u < 0.1):
        raise AssertionError("Expected a small positive initial uniform velocity")

    history = read_history(ROOT / "history" / "comparison_history.tsv")
    end_time = parser.getfloat("time_integration", "end_time")
    output_dt = parser.getfloat("steady_uniform", "comparison_output_dt")
    expected_times = np.arange(0.0, end_time + 0.5 * output_dt, output_dt)
    if history.size != expected_times.size:
        raise AssertionError(
            f"Expected {expected_times.size} output rows, found {history.size}"
        )
    if not np.allclose(history["t_star"], expected_times, atol=1.0e-11, rtol=0.0):
        raise AssertionError("Output times do not match the configured physical interval")

    first = history[0]
    for name in ("meanUlower", "meanUupper", "meanUair"):
        assert_close(float(first[name]), initial_u, 1.0e-12, f"initial {name}")

    final = history[-1]
    worst_liquid = max(float(final["relRmsUlower"]), float(final["relRmsUupper"]))
    if worst_liquid > 5.0e-3:
        raise AssertionError(f"Final liquid relative RMS error is too large: {worst_liquid}")
    if float(final["rmsPressure"]) > 1.0e-5:
        raise AssertionError("Gauge-corrected pressure RMS error is too large")
    if float(final["maxVerticalVelocity"]) > 1.0e-6:
        raise AssertionError("Vertical velocity is too large for the parallel benchmark")

    for volume_name in ("Vlower", "Vupper", "Vair"):
        drift = (float(final[volume_name]) - float(first[volume_name])) / float(first[volume_name])
        if abs(drift) > 1.0e-6:
            raise AssertionError(f"{volume_name} relative drift is too large: {drift}")

    expected_profiles = expected_times.size
    profiles = sorted((ROOT / "profiles").glob("profile_*_t*.tsv"))
    if len(profiles) != expected_profiles:
        raise AssertionError(
            f"Expected {expected_profiles} profiles, found {len(profiles)}"
        )

    combined = np.loadtxt(ROOT / "base_state" / "combined_profile.tsv", comments="#")
    z = combined[:, 0]
    u = combined[:, 1]
    lower = z <= 1.0 + 1.0e-12
    zl = z[lower]
    ul = u[lower]
    lower_mean = float(np.sum(0.5*(ul[1:] + ul[:-1])*np.diff(zl)))
    if abs(lower_mean - 1.0) > 5.0e-4:
        raise AssertionError(f"Generated lower mean is not normalized: {lower_mean}")

    for relative in (
        "analysis/comparison_report.txt",
        "analysis/final_profile_comparison.tsv",
        "analysis/final_velocity_profile.png",
        "analysis/mean_velocity_evolution.png",
        "analysis/velocity_error_decay.png",
    ):
        if not (ROOT / relative).is_file():
            raise FileNotFoundError(ROOT / relative)

    print("All steady-uniform benchmark checks passed.")
    print(f"Outputs: t*=0:{output_dt:g}:{end_time:g}")
    print(f"Final worst liquid relative RMS error: {worst_liquid:.8e}")
    print(f"Final gauge-corrected pressure RMS: {float(final['rmsPressure']):.8e}")


if __name__ == "__main__":
    main()
