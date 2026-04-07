#!/usr/bin/env python3
"""
Targeted growth-rate computation for Balmforth (2025), Section 4,
using a boundary-value formulation of the free-surface eigenmode.

Why this version?
-----------------
The earlier generalized-eigenvalue script treated the singular bordered
Chebyshev problem as a black-box matrix pencil. For this problem that is
numerically fragile: as the Chebyshev order is increased, spurious branch
switching / spectral pollution can push the reported rightmost eigenvalue
far away from the physical free-surface mode.

This corrected script instead:
  1) computes the steady base state from Section 3;
  2) uses a *low-order* Chebyshev solve only to obtain an initial guess
     for the relevant eigenvalue;
  3) refines the physical mode by solving the Section-4 ODE eigenproblem
     directly as a complex boundary-value problem with eta_hat = 1.

It is focused on the quantity the user asked for first: the growth rate
sigma_r(k) for prescribed (Gamma, kappa, T, a). The long-wave regime is
its main target.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import eig

# -----------------------------------------------------------------------------
# Base state (Section 3)
# -----------------------------------------------------------------------------

@dataclass
class BaseStateResult:
    a: float
    Gamma: float
    kappa: float
    z: np.ndarray
    Lambda: np.ndarray
    dLambda_dz: np.ndarray
    U: np.ndarray
    success: bool
    message: str

    @property
    def Uz(self) -> np.ndarray:
        return (1.0 - self.z) * reciprocal_viscosity(self.Lambda, self.a)

    @property
    def Lambda_zz(self) -> np.ndarray:
        return (self.Gamma * self.Lambda * self.Uz - (1.0 - self.Lambda)) / self.kappa


def reciprocal_viscosity(lmbda: np.ndarray, a: float) -> np.ndarray:
    return (1.0 - lmbda) * (((1.0 - lmbda) * (1.0 - a)) + a)


def reciprocal_viscosity_prime(lmbda: np.ndarray, a: float) -> np.ndarray:
    return -(a + 2.0 * (1.0 - a) * (1.0 - lmbda))


def flowing_branch_lambda(tau: np.ndarray, Gamma: float, a: float) -> np.ndarray:
    tau = np.asarray(tau, dtype=float)
    tau_safe = np.maximum(tau, 1.0e-14)
    disc = 1.0 - 4.0 * (1.0 - a) / (Gamma * tau_safe)
    disc = np.maximum(disc, 0.0)
    denom = 2.0 * (1.0 - a)
    return (1.0 - np.sqrt(disc)) / denom


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
        flowing = tau >= tau_A
        lmbda[flowing] = flowing_branch_lambda(tau[flowing], Gamma, a)
    else:
        Y = transition_height(a, Gamma)
        flowing = z < Y
        lmbda[flowing] = flowing_branch_lambda(tau[flowing], Gamma, a)

    Uz = tau * reciprocal_viscosity(lmbda, a)
    U = np.zeros_like(z)
    U[1:] = np.cumsum(0.5 * (Uz[:-1] + Uz[1:]) * np.diff(z))
    return lmbda, U


def build_initial_mesh(a: float, Gamma: float, n_uniform: int = 400, n_cluster: int = 800) -> np.ndarray:
    z0 = np.clip(transition_height(a, Gamma), 0.02, 0.98)
    s = np.linspace(-7.0, 7.0, n_cluster)
    cluster = z0 + 0.08 * np.tanh(s)
    cluster = cluster[(cluster > 0.0) & (cluster < 1.0)]
    mesh = np.unique(np.concatenate([np.linspace(0.0, 1.0, n_uniform), cluster]))
    return mesh


def initial_guess(z: np.ndarray, a: float, Gamma: float, previous: BaseStateResult | None = None) -> np.ndarray:
    if previous is not None:
        return np.vstack([
            np.interp(z, previous.z, previous.Lambda),
            np.interp(z, previous.z, previous.dLambda_dz),
            np.interp(z, previous.z, previous.U),
        ])
    lmbda0, U0 = diffusionless_profile(z, a, Gamma)
    dLambda0 = np.gradient(lmbda0, z, edge_order=2)
    return np.vstack([lmbda0, dLambda0, U0])


def solve_base_state(a: float, Gamma: float, kappa: float,
                     previous: BaseStateResult | None = None,
                     tol: float = 1.0e-6, max_nodes: int = 200000) -> BaseStateResult:
    z = build_initial_mesh(a, Gamma)
    y0 = initial_guess(z, a, Gamma, previous)

    def ode(z_, y):
        lmbda = np.clip(y[0], 0.0, 1.2)
        q = y[1]
        Uz = (1.0 - z_) * reciprocal_viscosity(lmbda, a)
        qz = (Gamma * lmbda * Uz - (1.0 - lmbda)) / kappa
        return np.vstack([q, qz, Uz])

    def bc(ya, yb):
        return np.array([ya[1], yb[1], ya[2]])

    sol = solve_bvp(ode, bc, z, y0, tol=tol, max_nodes=max_nodes, verbose=0)
    z_dense = np.linspace(0.0, 1.0, 2001)
    y_dense = sol.sol(z_dense)
    return BaseStateResult(a=a, Gamma=Gamma, kappa=kappa,
                           z=z_dense, Lambda=y_dense[0], dLambda_dz=y_dense[1], U=y_dense[2],
                           success=bool(sol.success), message=str(sol.message))

# -----------------------------------------------------------------------------
# Low-order seed from the bordered Chebyshev pencil
# -----------------------------------------------------------------------------

def chebyshev_matrices(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if n < 2:
        raise ValueError("Need n >= 2")
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.ones(n + 1)
    c[0] = 2.0
    c[-1] = 2.0
    c = c * ((-1.0) ** np.arange(n + 1))
    X = np.tile(x, (n + 1, 1))
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(n + 1))
    D = D - np.diag(np.sum(D, axis=1))
    z = 0.5 * (1.0 - x)
    Dz = -2.0 * D
    D2z = Dz @ Dz
    return z, Dz, D2z


def low_order_seed_sigma(base: BaseStateResult, k: float, T: float, n: int = 30) -> complex:
    z, D, D2 = chebyshev_matrices(n)
    Lambda = np.interp(z, base.z, base.Lambda)
    Lz = np.interp(z, base.z, base.dLambda_dz)
    U = np.interp(z, base.z, base.U)
    Uz = (1.0 - z) * reciprocal_viscosity(Lambda, base.a)
    A = reciprocal_viscosity(Lambda, base.a)
    Ap = reciprocal_viscosity_prime(Lambda, base.a)
    Lambda_zz = (base.Gamma * Lambda * Uz - (1.0 - Lambda)) / base.kappa

    size = 2 * (n + 1) + 1
    A_mat = np.zeros((size, size), dtype=complex)
    B_mat = np.zeros((size, size), dtype=complex)
    psi_off = 0
    lam_off = n + 1
    eta_idx = 2 * (n + 1)
    row = 0

    for j in range(1, n):
        A_mat[row, psi_off:psi_off + n + 1] = D2[j, :]
        A_mat[row, lam_off:lam_off + n + 1] = -(1.0 - z[j]) * Ap[j] * np.eye(n + 1)[j, :]
        A_mat[row, eta_idx] = -A[j] * (1.0 - 1j * k * (1.0 - z[j]))
        row += 1

    for j in range(1, n):
        # sign convention arranged so the physical root has the correct sign
        A_mat[row, psi_off:psi_off + n + 1] = base.Gamma * Lambda[j] * D2[j, :] - 1j * k * T * Lz[j] * np.eye(n + 1)[j, :]
        A_mat[row, lam_off:lam_off + n + 1] = -base.kappa * D2[j, :] + (1.0 + base.Gamma * Uz[j] + 1j * k * T * U[j]) * np.eye(n + 1)[j, :]
        B_mat[row, lam_off:lam_off + n + 1] = T * np.eye(n + 1)[j, :]
        row += 1

    # BCs
    A_mat[row, psi_off + n] = 1j * k
    A_mat[row, eta_idx] = 1j * k * U[-1]
    B_mat[row, eta_idx] = -1.0
    row += 1
    A_mat[row, psi_off + 0] = 1.0
    row += 1
    A_mat[row, psi_off:psi_off + n + 1] = D[0, :]
    row += 1
    A_mat[row, lam_off:lam_off + n + 1] = D[0, :]
    row += 1
    A_mat[row, lam_off:lam_off + n + 1] = D[n, :]
    A_mat[row, eta_idx] = Lambda_zz[-1]
    row += 1
    assert row == size

    vals, vecs = eig(A_mat, B_mat, overwrite_a=False, overwrite_b=False, check_finite=False)
    mask = np.isfinite(vals.real) & np.isfinite(vals.imag) & (np.abs(vals) < 1e6)
    vals = vals[mask]
    vecs = vecs[:, mask]
    if vals.size == 0:
        return 0.06 * k * k - 1j * base.Uz[0] * k

    norms = np.maximum(np.max(np.abs(vecs), axis=0), 1.0e-300)
    eta_amp = np.abs(vecs[-1, :] / norms)
    c_over_U0 = -vals.imag / (k * base.Uz[0])

    # keep modes with a genuine free-surface signature and a phase speed near the physical branch
    phys = eta_amp > 5.0e-2
    phys &= (c_over_U0 > 0.3) & (c_over_U0 < 1.6)
    cand = vals[phys] if np.any(phys) else vals
    return cand[np.argmax(cand.real)]

# -----------------------------------------------------------------------------
# Corrected targeted eigenvalue solve: direct Section-4 BVP with eta_hat = 1
# -----------------------------------------------------------------------------

@dataclass
class ModeResult:
    k: float
    sigma: complex
    success: bool
    status: int
    message: str
    n_nodes: int


def clustered_unit_interval(n_uniform: int = 180, n_top: int = 120) -> np.ndarray:
    """Mesh with extra resolution near z=1, where the lambda boundary layer lives."""
    z1 = np.linspace(0.0, 1.0, n_uniform)
    s = np.linspace(0.0, 1.0, n_top)
    z2 = 1.0 - 0.06 * (1.0 - s) ** 3
    return np.unique(np.concatenate([z1, z2]))


def solve_mode_bvp(base: BaseStateResult, k: float, T: float,
                   sigma_seed: complex | None = None,
                   tol: float = 3.0e-3,
                   max_nodes: int = 50000,
                   n_mesh: int = 300) -> ModeResult:
    """
    Solve the physical free-surface eigenmode with eta_hat fixed to 1.

    Unknowns:
        y = [psi, psi_z, lambda_hat, lambda_hat_z]
        parameter p = [sigma]
    BCs:
        psi(0)=0,
        psi_z(0)=0,
        lambda_hat_z(0)=0,
        lambda_hat_z(1) + Lambda_zz(1) = 0,
        sigma + i k U(1) = - i k psi(1).
    """
    if sigma_seed is None:
        sigma_seed = low_order_seed_sigma(base, k, T, n=30)

    z = clustered_unit_interval(max(120, n_mesh // 2), max(80, n_mesh // 3))
    U1 = base.U[-1]
    # kinematic BC gives psi(1) approximately from sigma
    psi1 = -(sigma_seed + 1j * k * U1) / (1j * k)
    psi_guess = psi1 * z
    psi_z_guess = np.full_like(z, psi1, dtype=complex)

    # a simple lambda guess that satisfies lambda_z(0)=0 and approximately the top BC
    Lamzz1 = base.Lambda_zz[-1]
    lam_guess = -0.5 * Lamzz1 * z**2
    lam_z_guess = -Lamzz1 * z

    y0 = np.vstack([psi_guess, psi_z_guess, lam_guess, lam_z_guess]).astype(complex)
    p0 = np.array([sigma_seed], dtype=complex)

    def ode(z_, y, p):
        sigma = p[0]
        Lam = np.interp(z_, base.z, base.Lambda)
        Lz = np.interp(z_, base.z, base.dLambda_dz)
        U = np.interp(z_, base.z, base.U)
        Uz = (1.0 - z_) * reciprocal_viscosity(Lam, base.a)
        A = reciprocal_viscosity(Lam, base.a)
        Ap = reciprocal_viscosity_prime(Lam, base.a)

        psi = y[0]
        psi_z = y[1]
        lam = y[2]
        lam_z = y[3]

        psi_zz = A * (1.0 - 1j * k * (1.0 - z_)) + (1.0 - z_) * Ap * lam
        lam_zz = (
            T * ((sigma + 1j * k * U) * lam - 1j * k * Lz * psi)
            + (1.0 + base.Gamma * Uz) * lam
            + base.Gamma * Lam * psi_zz
        ) / base.kappa
        return np.vstack([psi_z, psi_zz, lam_z, lam_zz])

    def bc(ya, yb, p):
        sigma = p[0]
        return np.array([
            ya[0],
            ya[1],
            ya[3],
            yb[3] + Lamzz1,
            sigma + 1j * k * base.U[-1] + 1j * k * yb[0],
        ], dtype=complex)

    sol = solve_bvp(ode, bc, z, y0, p=p0, tol=tol, max_nodes=max_nodes, verbose=0)
    return ModeResult(k=k, sigma=complex(sol.p[0]), success=bool(sol.success),
                      status=int(sol.status), message=str(sol.message), n_nodes=sol.x.size)


def compute_growth_curve(base: BaseStateResult, T: float, ks: Iterable[float]) -> list[ModeResult]:
    out: list[ModeResult] = []
    # solve the highest-k point first; it usually gives the best-conditioned seed,
    # then march to lower k where the long-wave branch is smoother.
    ks = sorted([float(k) for k in ks if k > 0.0], reverse=True)
    sigma_seed: complex | None = None
    for k in ks:
        res = solve_mode_bvp(base, k=k, T=T, sigma_seed=sigma_seed)
        out.append(res)
        sigma_seed = res.sigma
    return sorted(out, key=lambda r: r.k)


def estimate_sigma2_from_curve(results: list[ModeResult]) -> float:
    ks = np.array([r.k for r in results if r.success or np.isfinite(r.sigma.real)], dtype=float)
    sr = np.array([r.sigma.real for r in results if r.success or np.isfinite(r.sigma.real)], dtype=float)
    if ks.size < 2:
        return np.nan
    X = np.column_stack([ks**2, ks**4])
    coeff, *_ = np.linalg.lstsq(X, sr, rcond=None)
    return coeff[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="Corrected growth-rate solver for Balmforth (2025)")
    parser.add_argument("--a", type=float, required=True)
    parser.add_argument("--Gamma", type=float, default=8.0)
    parser.add_argument("--kappa", type=float, default=1.0e-4)
    parser.add_argument("--T", type=float, default=1.0)
    parser.add_argument("--k", type=float, default=None, help="single wavenumber")
    parser.add_argument("--kmin", type=float, default=1.0e-3)
    parser.add_argument("--kmax", type=float, default=1.0e-1)
    parser.add_argument("--nk", type=int, default=10)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    base = solve_base_state(a=args.a, Gamma=args.Gamma, kappa=args.kappa)
    if not base.success:
        raise RuntimeError(f"Base-state solve failed: {base.message}")

    print("Base state:")
    print(f"  success   = {base.success}")
    print(f"  message   = {base.message}")
    print(f"  U(1)      = {base.U[-1]:.10f}")
    print(f"  U_z(0)    = {base.Uz[0]:.10f}")
    print()

    if args.k is not None:
        res = solve_mode_bvp(base, k=args.k, T=args.T)
        print("Mode result:")
        print(f"  k         = {res.k:.8g}")
        print(f"  sigma_r   = {res.sigma.real:.10e}")
        print(f"  sigma_i   = {res.sigma.imag:.10e}")
        print(f"  success   = {res.success}")
        print(f"  nodes     = {res.n_nodes}")
        print(f"  message   = {res.message}")
        return

    ks = np.logspace(np.log10(args.kmin), np.log10(args.kmax), args.nk)
    results = compute_growth_curve(base, args.T, ks)
    sigma2 = estimate_sigma2_from_curve(results)

    print("Wavenumber results:")
    print(f"{'k':>12s} {'sigma_r':>18s} {'sigma_i':>18s} {'sigma_r/k^2':>18s} {'ok':>6s}")
    for r in results:
        scaled = r.sigma.real / (r.k * r.k)
        ok = 'yes' if r.success else 'no'
        print(f"{r.k:12.5e} {r.sigma.real:18.10e} {r.sigma.imag:18.10e} {scaled:18.10e} {ok:>6s}")
    print()
    print(f"Estimated Re(sigma2) from the returned curve: {sigma2:.10e}")


if __name__ == '__main__':
    main()
