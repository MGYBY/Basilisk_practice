#!/usr/bin/env python3
"""Create report figures from the readable full-NS stability outputs."""
from __future__ import annotations

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig

from fullns_stability import (
    CaseParameters,
    SolverOptions,
    assemble_phase_speed_pencil,
    build_base_flow,
    equilibrate_rows,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "report" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def save(figure, name: str) -> None:
    figure.tight_layout()
    figure.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    figure.savefig(FIGURES / f"{name}.png", dpi=240, bbox_inches="tight")
    plt.close(figure)


def read_base_profiles():
    rows = []
    with (DATA / "base_state_profiles.tsv").open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            rows.append((fields[0],) + tuple(float(value) for value in fields[1:]))
    dtype = [
        ("layer", "U8"),
        ("z", float),
        ("velocity", float),
        ("pressure", float),
        ("shear", float),
        ("mu", float),
        ("mu_t", float),
    ]
    return np.array(rows, dtype=dtype)


def read_dispersion():
    dtype = [
        ("k", float),
        ("branch", "U32"),
        ("c_real", float),
        ("c_imag", float),
        ("sigma_real", float),
        ("sigma_imag", float),
        ("growth", float),
        ("celerity", float),
        ("eta_real", float),
        ("eta_imag", float),
        ("zeta_real", float),
        ("zeta_imag", float),
        ("ratio", float),
        ("phase", float),
        ("residual", float),
    ]
    return np.genfromtxt(
        DATA / "dispersion_branches.tsv",
        comments="#",
        delimiter="\t",
        dtype=dtype,
    )


def read_convergence(name: str):
    return np.loadtxt(DATA / name, comments="#")




def read_regularization_sensitivity():
    dtype = [
        ("rate", float),
        ("epsilon", float),
        ("branch", "U32"),
        ("maximum_growth", float),
        ("k_max", float),
        ("celerity_max", float),
        ("neutral_k", float),
    ]
    return np.genfromtxt(
        DATA / "regularization_sensitivity.tsv",
        comments="#",
        delimiter="\t",
        dtype=dtype,
    )


def conditioning_comparison() -> np.ndarray:
    """Compare raw and row-equilibrated interface roots at k=0.01."""

    orders = (32, 40, 52, 64, 76)
    scaled = read_convergence("convergence_k0p01.tsv")
    base = build_base_flow(CaseParameters())
    rows = []

    for row, order in zip(scaled, orders):
        target_growth = row[1]
        target_celerity = row[2]
        target_c_imag = target_growth / 0.01
        target = target_celerity + 1j * target_c_imag

        pencil = assemble_phase_speed_pencil(base, 0.01, order)
        raw_values = eig(
            pencil.A,
            pencil.B,
            right=False,
            check_finite=False,
        )
        finite = raw_values[
            np.isfinite(raw_values.real)
            & np.isfinite(raw_values.imag)
            & (np.abs(raw_values) < 100.0)
        ]
        raw_root = min(finite, key=lambda value: abs(value - target))
        raw_growth = float((-1j * 0.01 * raw_root).real)

        _, _, row_scales = equilibrate_rows(pencil.A, pencil.B)
        rows.append(
            (
                float(order),
                raw_growth,
                float(target_growth),
                float(np.max(row_scales) / np.min(row_scales)),
            )
        )

    result = np.asarray(rows)
    np.savetxt(
        DATA / "conditioning_comparison_k0p01.tsv",
        result,
        delimiter="\t",
        header=(
            "order\tunscaled_interface_growth\t"
            "row_equilibrated_interface_growth\trow_norm_spread"
        ),
        comments="# ",
    )
    return result


def main() -> None:
    base_data = read_base_profiles()
    dispersion = read_dispersion()
    lower = base_data[base_data["layer"] == "lower"]
    upper = base_data[base_data["layer"] == "upper"]

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.plot(lower["velocity"], lower["z"], label="lower")
    ax.plot(upper["velocity"], upper["z"], label="upper")
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"base velocity $U_i(z)$")
    ax.set_ylabel(r"normal coordinate $z$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "base_velocity")

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.plot(lower["pressure"], lower["z"], label="lower")
    ax.plot(upper["pressure"], upper["z"], label="upper")
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"base pressure $P_i(z)$")
    ax.set_ylabel(r"$z$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "base_pressure")

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.plot(lower["shear"], lower["z"], label="lower")
    ax.plot(upper["shear"], upper["z"], label="upper")
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"base shear $U_i'(z)$")
    ax.set_ylabel(r"$z$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "base_shear")

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.semilogx(lower["mu"], lower["z"], label=r"$\mu_\ell$")
    ax.semilogx(lower["mu_t"], lower["z"], linestyle="--", label=r"$\mu_{t,\ell}$")
    ax.semilogx(upper["mu"], upper["z"], label=r"$\mu_u$")
    ax.semilogx(upper["mu_t"], upper["z"], linestyle="--", label=r"$\mu_{t,u}$")
    ax.axhline(1.0, linestyle=":", linewidth=1.0)
    ax.set_xlabel("dimensionless viscosity")
    ax.set_ylabel(r"$z$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(ncol=2, fontsize=9)
    save(fig, "base_viscosities")

    labels = {
        "interface_varicose": "interface / varicose",
        "surface_zigzag": "surface / zigzag",
    }

    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    for branch, label in labels.items():
        values = dispersion[dispersion["branch"] == branch]
        ax.plot(values["k"], values["growth"], label=label)
    ax.axhline(0.0, linestyle="--", linewidth=1.0)
    ax.set_xlim(0.0, 0.5)
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel(r"growth rate $\Re(\sigma)$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "dispersion_growth_main")

    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    for branch, label in labels.items():
        values = dispersion[dispersion["branch"] == branch]
        mask = values["k"] <= 0.08
        ax.plot(values["k"][mask], values["growth"][mask], label=label)
    ax.axhline(0.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel(r"growth rate $\Re(\sigma)$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "dispersion_growth_longwave")

    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    for branch, label in labels.items():
        values = dispersion[dispersion["branch"] == branch]
        ax.plot(values["k"], values["celerity"], label=label)
    ax.set_xlim(0.0, 0.5)
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel(r"phase celerity $\Re(c)$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "dispersion_celerity")

    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    for branch, label in labels.items():
        values = dispersion[dispersion["branch"] == branch]
        ax.plot(values["k"], values["ratio"], label=label)
    ax.set_xlim(0.0, 0.5)
    ax.set_yscale("log")
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel(r"$|\widehat\zeta/\widehat\eta|$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    save(fig, "mode_amplitude_ratio")

    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    for branch, label in labels.items():
        values = dispersion[dispersion["branch"] == branch]
        ax.plot(values["k"], values["phase"], label=label)
    ax.axhline(0.0, linestyle=":", linewidth=1.0)
    ax.axhline(180.0, linestyle=":", linewidth=1.0)
    ax.axhline(-180.0, linestyle=":", linewidth=1.0)
    ax.set_xlim(0.0, 0.5)
    ax.set_ylim(-190.0, 190.0)
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel("relative phase (degrees)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "mode_phase_relation")

    for filename, k_label, output_name in (
        ("convergence_k0p01.tsv", "0.01", "convergence_k0p01"),
        ("convergence_k0p4.tsv", "0.4", "convergence_k0p4"),
    ):
        values = read_convergence(filename)
        fig, ax = plt.subplots(figsize=(5.9, 4.2))
        ax.plot(values[:, 0], values[:, 1], marker="o", label="interface")
        ax.plot(values[:, 0], values[:, 4], marker="s", label="surface")
        ax.set_xlabel("Chebyshev order N")
        ax.set_ylabel(rf"growth rate at $k={k_label}$")
        ax.grid(True, alpha=0.25)
        ax.legend()
        save(fig, output_name)



    regularization = read_regularization_sensitivity()
    fig, ax = plt.subplots(figsize=(6.1, 4.2))
    for branch, label in labels.items():
        values = regularization[regularization["branch"] == branch]
        ax.semilogy(values["rate"], values["maximum_growth"], marker="o", label=label)
    ax.set_xlabel(r"dimensional smoothing rate $\widetilde{\dot\gamma}_r$ [s$^{-1}$]")
    ax.set_ylabel(r"maximum $\Re(\sigma)$ on the sampled band")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    save(fig, "regularization_maximum_growth")

    fig, ax = plt.subplots(figsize=(6.1, 4.2))
    for branch, label in labels.items():
        values = regularization[regularization["branch"] == branch]
        ax.plot(values["rate"], values["neutral_k"], marker="o", label=label)
    ax.set_xlabel(r"dimensional smoothing rate $\widetilde{\dot\gamma}_r$ [s$^{-1}$]")
    ax.set_ylabel(r"first neutral wavenumber $k_c$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "regularization_neutral_wavenumber")


    conditioning = conditioning_comparison()
    fig, ax = plt.subplots(figsize=(6.1, 4.2))
    ax.plot(
        conditioning[:, 0],
        conditioning[:, 1],
        marker="o",
        label="unscaled pencil",
    )
    ax.plot(
        conditioning[:, 0],
        conditioning[:, 2],
        marker="s",
        label="row-equilibrated pencil",
    )
    ax.axhline(0.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Chebyshev order N")
    ax.set_ylabel(r"interface growth at $k=0.01$")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "conditioning_effect")

    summary = json.loads((DATA / "case_summary.json").read_text(encoding="utf-8"))
    (DATA / "figure_metrics.json").write_text(
        json.dumps(
            {
                "case_summary": summary,
                "conditioning_rows": conditioning.tolist(),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote figures to {FIGURES}")


if __name__ == "__main__":
    main()
