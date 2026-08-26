#!/usr/bin/env python3
"""Summarize the relaxation run and compare it with the reference solution."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_named_tsv(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, names=True, delimiter="\t", dtype=float)
    if data.shape == ():
        data = np.asarray([data], dtype=data.dtype)
    return data


def latest_profile(profile_dir: Path) -> Path:
    files = sorted(profile_dir.glob("profile_*_t*.tsv"))
    if not files:
        raise FileNotFoundError(f"No profile files found in {profile_dir}")
    return files[-1]


def first_time_below(time: np.ndarray, values: np.ndarray, threshold: float) -> float:
    indices = np.nonzero(values <= threshold)[0]
    return float(time[indices[0]]) if indices.size else float("nan")


def parse_initial_velocity(path: Path) -> float:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"^initial uniform u_x\*\s*=\s*([0-9.eE+-]+)", text, re.M)
    if not match:
        raise ValueError(f"Could not parse initial velocity from {path}")
    return float(match.group(1))


def phase_rms(profile: np.ndarray, phase_name: str, numeric: str, exact: str) -> tuple[float, float]:
    phase = profile[phase_name]
    mask = phase > 0.999
    if not np.any(mask):
        return float("nan"), float("nan")
    error = profile[numeric][mask] - profile[exact][mask]
    rms = float(np.sqrt(np.mean(error * error)))
    scale = float(np.sqrt(np.mean(profile[exact][mask] ** 2)))
    return rms, rms / max(scale, 1.0e-300)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()

    root = args.root.resolve()
    analysis = root / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)

    history_path = root / "history" / "comparison_history.tsv"
    history = load_named_tsv(history_path)
    profile_path = latest_profile(root / "profiles")
    profile = load_named_tsv(profile_path)
    initial_velocity = parse_initial_velocity(
        root / "base_state" / "steady_uniform_case_summary.txt"
    )

    final = history[-1]
    worst_liquid_rel = np.maximum(history["relRmsUlower"], history["relRmsUupper"])

    gamma_l_rms, gamma_l_rel = phase_rms(
        profile, "alphaLower", "gammaNum_star", "gammaExact_star"
    )
    gamma_u_rms, gamma_u_rel = phase_rms(
        profile, "alphaUpper", "gammaNum_star", "gammaExact_star"
    )
    eta_l_rms, eta_l_rel = phase_rms(
        profile, "alphaLower", "etaNum_star", "etaExact_star"
    )
    eta_u_rms, eta_u_rel = phase_rms(
        profile, "alphaUpper", "etaNum_star", "etaExact_star"
    )

    volume_drift_lower = (final["Vlower"] - history[0]["Vlower"]) / history[0]["Vlower"]
    volume_drift_upper = (final["Vupper"] - history[0]["Vupper"]) / history[0]["Vupper"]
    volume_drift_air = (final["Vair"] - history[0]["Vair"]) / history[0]["Vair"]

    report = analysis / "comparison_report.txt"
    with report.open("w", encoding="utf-8") as stream:
        stream.write("STEADY-UNIFORM RELAXATION COMPARISON\n")
        stream.write("====================================\n\n")
        stream.write(f"Final profile file                  : {profile_path.name}\n")
        stream.write(f"Initial uniform velocity u_x*       : {initial_velocity:.12g}\n")
        stream.write(f"Final comparison time t*            : {final['t_star']:.12g}\n")
        stream.write("\nLayer-mean streamwise velocities\n")
        stream.write("--------------------------------\n")
        stream.write(
            f"Lower: numerical={final['meanUlower']:.12g}, "
            f"reference={final['exactMeanUlower']:.12g}, "
            f"relative difference={(final['meanUlower']/final['exactMeanUlower']-1):.6e}\n"
        )
        stream.write(
            f"Upper: numerical={final['meanUupper']:.12g}, "
            f"reference={final['exactMeanUupper']:.12g}, "
            f"relative difference={(final['meanUupper']/final['exactMeanUupper']-1):.6e}\n"
        )
        stream.write(
            f"Air  : numerical={final['meanUair']:.12g}, "
            f"reference={final['exactMeanUair']:.12g}, "
            f"relative difference={(final['meanUair']/final['exactMeanUair']-1):.6e}\n"
        )
        stream.write("\nVelocity-profile errors\n")
        stream.write("-----------------------\n")
        stream.write(
            f"Lower relative RMS error             : {final['relRmsUlower']:.6e}\n"
        )
        stream.write(
            f"Upper relative RMS error             : {final['relRmsUupper']:.6e}\n"
        )
        stream.write(
            f"Air relative RMS error               : {final['relRmsUair']:.6e}\n"
        )
        stream.write(
            f"Lower/upper worst relative RMS       : {max(final['relRmsUlower'], final['relRmsUupper']):.6e}\n"
        )
        stream.write(
            f"Interface velocity error             : {(final['Uinterface']-final['exactUinterface']):.6e}\n"
        )
        stream.write(
            f"Free-surface velocity error          : {(final['Usurface']-final['exactUsurface']):.6e}\n"
        )
        stream.write("\nPressure and transverse motion\n")
        stream.write("-------------------------------\n")
        stream.write(
            f"Gauge-corrected pressure RMS error   : {final['rmsPressure']:.6e}\n"
        )
        stream.write(
            f"Gauge-corrected pressure Linf error  : {final['linfPressure']:.6e}\n"
        )
        stream.write(
            f"Maximum vertical velocity            : {final['maxVerticalVelocity']:.6e}\n"
        )
        stream.write("\nPure-phase strain-rate/viscosity diagnostics\n")
        stream.write("--------------------------------------------\n")
        stream.write(
            f"Lower gamma RMS / relative RMS       : {gamma_l_rms:.6e} / {gamma_l_rel:.6e}\n"
        )
        stream.write(
            f"Upper gamma RMS / relative RMS       : {gamma_u_rms:.6e} / {gamma_u_rel:.6e}\n"
        )
        stream.write(
            f"Lower eta RMS / relative RMS         : {eta_l_rms:.6e} / {eta_l_rel:.6e}\n"
        )
        stream.write(
            f"Upper eta RMS / relative RMS         : {eta_u_rms:.6e} / {eta_u_rel:.6e}\n"
        )
        stream.write("\nPhase-volume conservation\n")
        stream.write("-------------------------\n")
        stream.write(f"Lower relative drift                : {volume_drift_lower:.6e}\n")
        stream.write(f"Upper relative drift                : {volume_drift_upper:.6e}\n")
        stream.write(f"Air relative drift                  : {volume_drift_air:.6e}\n")
        stream.write("\nRelaxation times from the sampled history\n")
        stream.write("-----------------------------------------\n")
        for threshold in (5.0e-2, 1.0e-2, 5.0e-3, 2.5e-3):
            stream.write(
                f"First t* with worst liquid relative RMS <= {threshold:.4g}: "
                f"{first_time_below(history['t_star'], worst_liquid_rel, threshold):.12g}\n"
            )
        stream.write("\nInterpretation\n")
        stream.write("--------------\n")
        stream.write(
            "The pressure comparison removes one spatially constant offset because the "
            "closed periodic pressure field is defined only up to a gauge constant.\n"
        )
        stream.write(
            "The remaining velocity error contains both incomplete temporal relaxation "
            "and finite-grid error. The supplied level-6 sample has 16 cells through "
            "each liquid layer; increase the common mesh level for a convergence study.\n"
        )

    # Preserve the exact final table as a convenient publication/postprocessing input.
    header = "\t".join(profile.dtype.names or ())
    np.savetxt(
        analysis / "final_profile_comparison.tsv",
        np.column_stack([profile[name] for name in profile.dtype.names or ()]),
        delimiter="\t",
        header=header,
        comments="",
        fmt="%.16e",
    )

    # Mean-velocity evolution.
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(history["t_star"], history["meanUlower"], label="lower numerical")
    ax.plot(history["t_star"], history["meanUupper"], label="upper numerical")
    ax.plot(history["t_star"], history["meanUair"], label="air numerical")
    ax.plot(history["t_star"], history["exactMeanUlower"], "--", label="lower reference")
    ax.plot(history["t_star"], history["exactMeanUupper"], "--", label="upper reference")
    ax.plot(history["t_star"], history["exactMeanUair"], "--", label="air reference")
    ax.set_xlabel(r"dimensionless time $t^*$")
    ax.set_ylabel(r"layer-mean velocity $\overline{u}^*$")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(analysis / "mean_velocity_evolution.png", dpi=220)
    plt.close(fig)

    # Final profile.
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.plot(profile["uNum_star"], profile["z_star"], label="Basilisk")
    ax.plot(profile["uExact_star"], profile["z_star"], "--", label="reference")
    ax.axhline(1.0, linewidth=0.8, linestyle=":")
    ax.axhline(2.0, linewidth=0.8, linestyle=":")
    ax.set_xlabel(r"streamwise velocity $u^*$")
    ax.set_ylabel(r"normal coordinate $z^*$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(analysis / "final_velocity_profile.png", dpi=220)
    plt.close(fig)

    # Error decay.
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.semilogy(history["t_star"], history["relRmsUlower"], label="lower")
    ax.semilogy(history["t_star"], history["relRmsUupper"], label="upper")
    ax.semilogy(history["t_star"], history["relRmsUair"], label="air")
    ax.set_xlabel(r"dimensionless time $t^*$")
    ax.set_ylabel("relative RMS velocity error")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(analysis / "velocity_error_decay.png", dpi=220)
    plt.close(fig)

    print(report.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
