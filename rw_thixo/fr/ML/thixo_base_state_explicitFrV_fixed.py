#!/usr/bin/env python3
"""
Steady base state for the explicit-FrV multilayer formulation.

Important:
Although FrV appears explicitly in the transient multilayer momentum equation,
it cancels out of the steady uniform base-state equations. Therefore the steady
Balmforth base state is identical to the original Balmforth-scaled one and does
not depend on FrV.

This script writes three columns:
    z   U(z)   Lambda(z)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp


def A_of_lambda(lmbda: np.ndarray, a: float) -> np.ndarray:
    return (1.0 - lmbda) * (((1.0 - lmbda) * (1.0 - a)) + a)


def flowing_branch_lambda(tau: np.ndarray, Gamma: float, a: float) -> np.ndarray:
    tau = np.asarray(tau, dtype=float)
    tau_safe = np.maximum(tau, 1.0e-15)
    disc = np.maximum(1.0 - 4.0 * (1.0 - a) / (Gamma * tau_safe), 0.0)
    return (1.0 - np.sqrt(disc)) / (2.0 * (1.0 - a))


def transition_height(a: float, Gamma: float) -> float:
    if a > 0.5:
        return 1.0 - 1.0 / (Gamma * a)
    return 1.0 - 9.0 * (1.0 - a) / (Gamma * (2.0 - a) * (1.0 + a))


def diffusionless_profile(z: np.ndarray, a: float, Gamma: float) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    tau = 1.0 - z
    lmbda = np.ones_like(z)

    if a > 0.5:
        tau_A = 1.0 / (Gamma * a)
        mask = tau >= tau_A
        lmbda[mask] = flowing_branch_lambda(tau[mask], Gamma, a)
    else:
        Y = transition_height(a, Gamma)
        mask = z < Y
        lmbda[mask] = flowing_branch_lambda(tau[mask], Gamma, a)

    # In the steady uniform base state for the explicit-FrV multilayer model,
    # the FrV^{-2} factors cancel from the momentum balance.
    Uz = (1.0 - z) * A_of_lambda(lmbda, a)

    U = np.zeros_like(z)
    if len(z) > 1:
        U[1:] = np.cumsum(0.5 * (Uz[:-1] + Uz[1:]) * np.diff(z))
    return lmbda, U


def build_initial_mesh(a: float, Gamma: float, n_uniform: int = 500, n_cluster: int = 900) -> np.ndarray:
    z0 = transition_height(a, Gamma)
    s = np.linspace(-7.0, 7.0, n_cluster)
    cluster = z0 + 0.08 * np.tanh(s)
    cluster = cluster[(cluster > 0.0) & (cluster < 1.0)]
    return np.unique(np.concatenate([np.linspace(0.0, 1.0, n_uniform), cluster]))


def solve_base_state(a: float, Gamma: float, kappa: float, tol: float, max_nodes: int):
    if kappa <= 0.0:
        raise ValueError(f"kappa must be positive; got {kappa}")

    z = build_initial_mesh(a, Gamma)
    l0, U0 = diffusionless_profile(z, a, Gamma)
    q0 = np.gradient(l0, z, edge_order=2)
    y0 = np.vstack([l0, q0, U0])

    def ode(z_, y):
        lmbda = np.clip(y[0], 0.0, 1.2)
        q = y[1]
        Uz = (1.0 - z_) * A_of_lambda(lmbda, a)
        qz = (Gamma * lmbda * Uz - (1.0 - lmbda)) / kappa
        return np.vstack([q, qz, Uz])

    def bc(ya, yb):
        return np.array([ya[1], yb[1], ya[2]])

    return solve_bvp(ode, bc, z, y0, tol=tol, max_nodes=max_nodes, verbose=0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Solve the steady base state and write z U Lambda to a text file.')
    p.add_argument('--kappa', type=float, required=True, help='Positive diffusion coefficient kappa.')
    p.add_argument('--a', type=float, required=True, help='Parameter a.')
    p.add_argument('--Gamma', type=float, required=True, help='Parameter Gamma.')
    p.add_argument('--FrV', type=float, default=None,
                   help='Accepted for workflow compatibility, but it does not enter the steady base-state equations.')
    p.add_argument('--output', type=str, default='base_state_profile.txt', help='Output text filename.')
    p.add_argument('--n-rows', type=int, default=128, help='Number of rows written to the output text file.')
    p.add_argument('--tol', type=float, default=1.0e-6, help='solve_bvp tolerance.')
    p.add_argument('--max-nodes', type=int, default=400000, help='solve_bvp max_nodes.')
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_rows < 2:
        raise ValueError('--n-rows must be at least 2')

    sol = solve_base_state(args.a, args.Gamma, args.kappa, args.tol, args.max_nodes)
    if not sol.success:
        raise RuntimeError(f"Base-state BVP failed: {sol.message}")

    z_out = np.linspace(0.0, 1.0, args.n_rows)
    y_out = sol.sol(z_out)
    lmbda_out = y_out[0]
    U_out = y_out[2]

    out = np.column_stack([z_out, U_out, lmbda_out])
    outpath = Path(args.output)
    np.savetxt(outpath, out, fmt='%g')

    print('Solved base state successfully.')
    print(f'a = {args.a:g}, Gamma = {args.Gamma:g}, kappa = {args.kappa:.6e}')
    if args.FrV is not None:
        print(f'FrV = {args.FrV:g} (not used in the steady base-state equations)')
    print(f'Rows written = {len(z_out)}')
    print(f'Output file = {outpath}')
    print(f'U(1) = {U_out[-1]:.12g}')
    print(f'Lambda(1) = {lmbda_out[-1]:.12g}')


if __name__ == '__main__':
    main()
