#!/usr/bin/env python3
"""Create report figures from the saved threshold tables."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import argparse
import csv
import json
import math

import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update(
    {
        "font.size": 9.0,
        "axes.labelsize": 9.5,
        "axes.titlesize": 10.0,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
)


def read_neutral_table(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open("r", encoding="utf-8") as handle:
        header = None
        for line in handle:
            if line.startswith("#"):
                header = line[1:].strip().split("\t")
                continue
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if header is None:
                raise RuntimeError(f"Missing header in {path}")
            row: dict[str, float | str] = {}
            for name, value in zip(header, values):
                if name == "branch":
                    row[name] = value
                else:
                    row[name] = float(value)
            rows.append(row)
    return rows


def save_figure(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.pdf")
    fig.savefig(output / f"{stem}.png")
    plt.close(fig)


def long_wave_figure(data: Path, output: Path) -> None:
    source = data / "critical_froude_longwave.json"
    if not source.exists():
        return
    content = json.loads(source.read_text(encoding="utf-8"))

    fig, ax = plt.subplots(figsize=(5.7, 3.8))
    plotted = False
    for branch, entries in content["branches"].items():
        for index, entry in enumerate(entries):
            ks = np.asarray(entry["wavenumbers"], dtype=float)
            roots = np.asarray(entry["finite_k_roots"], dtype=float)
            coeff = np.asarray(entry["fit_coefficients_in_k2"], dtype=float)
            x = ks * ks
            label = branch.replace("_", " ")
            if not entry["converged"]:
                label += " (diagnostic only)"
            line, = ax.plot(x, roots, "o", label=label)
            xfit = np.linspace(0.0, 1.04 * float(np.max(x)), 200)
            yfit = np.polynomial.polynomial.polyval(xfit, coeff)
            ax.plot(xfit, yfit, linestyle="--", linewidth=1.1, color=line.get_color())
            ax.plot(
                [0.0],
                [entry["critical_froude"]],
                marker="s",
                color=line.get_color(),
            )
            plotted = True

    ax.set_xlabel(r"$k_H^2$")
    ax.set_ylabel(r"finite-$k$ neutral $Fr_{\ell,c}(k_H)$")
    ax.set_title("Extrapolation of branch-specific neutral roots to infinite wavelength")
    ax.grid(True, alpha=0.28)
    if plotted:
        ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output, "longwave_critical_extrapolation")


def _curve_label(path: Path) -> str:
    stem = path.stem
    if "_4p72" in stem:
        return r"$S_0\lambda/H_\ell=4.72$"
    if "_0p472" in stem:
        return r"$S_0\lambda/H_\ell=0.472$"
    return stem.replace("neutral_curve__", "").replace("_", " ")


def _group_by_root(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[int(row["root_number"])].append(row)
    for root in groups:
        groups[root].sort(key=lambda item: float(item["Re_lower"]))
    return groups


def neutral_maps(data: Path, output: Path) -> None:
    surface_files = sorted(data.glob("neutral_curve__surface_zigzag__*.tsv"))
    interface_files = sorted(data.glob("neutral_curve__interface_varicose__*.tsv"))
    if not surface_files and not interface_files:
        return

    longwave = None
    longwave_path = data / "critical_froude_longwave.json"
    if longwave_path.exists():
        content = json.loads(longwave_path.read_text(encoding="utf-8"))
        for entry in content.get("branches", {}).get("surface_zigzag", []):
            if entry.get("converged"):
                longwave = float(entry["critical_froude"])
                break
        chi = float(content["compatibility_ratio_chi_lower"])
    else:
        chi = float("nan")

    # Fr-Re map.
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for path in surface_files:
        groups = _group_by_root(read_neutral_table(path))
        for root_number, rows in groups.items():
            label = _curve_label(path)
            if root_number:
                label += f" (root {root_number + 1})"
            ax.plot(
                [float(row["Fr_lower"]) for row in rows],
                [float(row["Re_lower"]) for row in rows],
                marker="o",
                markersize=3.0,
                linewidth=1.25,
                label=label,
            )
    for path in interface_files:
        groups = _group_by_root(read_neutral_table(path))
        for root_number, rows in groups.items():
            if not rows:
                continue
            ax.plot(
                [float(row["Fr_lower"]) for row in rows],
                [float(row["Re_lower"]) for row in rows],
                marker="x",
                markersize=3.0,
                linewidth=0.9,
                linestyle="--",
                label="interface: " + _curve_label(path),
            )

    if longwave is not None:
        ax.axvline(longwave, linestyle="--", linewidth=1.2, label=r"surface, $k_H\to0$")

    # Compatibility lines, analogous to the grey lines in Yu & Chu (2024).
    fr_line = np.linspace(0.1, 1.6, 300)
    if np.isfinite(chi):
        for slope in (0.01, 0.02, 0.05, 0.10, 0.20, 0.50):
            re_line = chi * fr_line**2 / slope
            ax.plot(fr_line, re_line, linestyle=":", linewidth=0.65, alpha=0.55)
            valid = np.where((re_line > 2.5) & (re_line < 300.0))[0]
            if valid.size:
                idx = valid[-1]
                ax.text(
                    fr_line[idx],
                    re_line[idx],
                    rf" $S_0={slope:g}$",
                    fontsize=6.4,
                    rotation=24,
                    va="bottom",
                )

    ax.set_yscale("log")
    ax.set_xlim(0.1, 1.6)
    ax.set_ylim(2.5, 300.0)
    ax.set_xlabel(r"lower-layer Froude number $Fr_\ell$")
    ax.set_ylabel(r"lower-layer generalized Reynolds number $Re_\ell^{(n_\ell)}$")
    ax.set_title(r"Finite-wavelength neutral boundaries in the $Fr_\ell$-$Re_\ell$ plane")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    save_figure(fig, output, "neutral_map_fr_re")

    # Coordinate transformation of exactly the same neutral points.
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for path in surface_files:
        groups = _group_by_root(read_neutral_table(path))
        for root_number, rows in groups.items():
            label = _curve_label(path)
            if root_number:
                label += f" (root {root_number + 1})"
            ax.plot(
                [float(row["Fr_lower"]) for row in rows],
                [float(row["S0"]) for row in rows],
                marker="o",
                markersize=3.0,
                linewidth=1.25,
                label=label,
            )
    for path in interface_files:
        groups = _group_by_root(read_neutral_table(path))
        for rows in groups.values():
            if rows:
                ax.plot(
                    [float(row["Fr_lower"]) for row in rows],
                    [float(row["S0"]) for row in rows],
                    marker="x",
                    markersize=3.0,
                    linewidth=0.9,
                    linestyle="--",
                    label="interface: " + _curve_label(path),
                )

    if longwave is not None:
        ax.axvline(longwave, linestyle="--", linewidth=1.2, label=r"surface, $k_H\to0$")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 1.6)
    ax.set_xlabel(r"lower-layer Froude number $Fr_\ell$")
    ax.set_ylabel(r"inclination parameter $S_0=\tan\theta$")
    ax.set_title(r"The same neutral boundaries in the $Fr_\ell$-$S_0$ plane")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    save_figure(fig, output, "neutral_map_fr_slope")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = args.data.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    long_wave_figure(data, output)
    neutral_maps(data, output)
    print(f"Wrote threshold figures to {output}")


if __name__ == "__main__":
    main()
