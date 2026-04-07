
#!/usr/bin/env python3
from __future__ import annotations
import argparse, time, sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional, Sequence
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.linalg import eig

class Logger:
    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self.wall0 = time.perf_counter()
    def log(self, msg: str) -> None:
        if self.enabled:
            stamp = datetime.now().strftime("%H:%M:%S")
            wall = time.perf_counter() - self.wall0
            print(f"[{stamp} | wall +{wall:8.1f}s] {msg}", flush=True)

def reciprocal_viscosity(lmbda: np.ndarray, a: float) -> np.ndarray:
    return (1.0 - lmbda) * (((1.0 - lmbda) * (1.0 - a)) + a)

def transition_height(a: float, Gamma: float) -> float:
    if a > 0.5:
        return 1.0 - 1.0 / (Gamma * a)
    return 1.0 - 9.0 * (1.0 - a) / (Gamma * (2.0 - a) * (1.0 + a))

def flowing_branch_lambda(tau: np.ndarray, Gamma: float, a: float) -> np.ndarray:
    tau = np.maximum(np.asarray(tau, dtype=float), 1.0e-14)
    disc = np.maximum(1.0 - 4.0 * (1.0 - a) / (Gamma * tau), 0.0)
    return (1.0 - np.sqrt(disc)) / (2.0 * (1.0 - a))

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
    Uz = tau * reciprocal_viscosity(lmbda, a)
    U = np.zeros_like(z)
    U[1:] = np.cumsum(0.5 * (Uz[:-1] + Uz[1:]) * np.diff(z))
    return lmbda, U

@dataclass
class BaseState:
    z: np.ndarray
    Lambda: np.ndarray
    U: np.ndarray
    success: bool
    message: str

def build_initial_mesh(a: float, Gamma: float, n_uniform: int = 300, n_cluster: int = 700) -> np.ndarray:
    z0 = np.clip(transition_height(a, Gamma), 0.02, 0.98)
    s = np.linspace(-7.0, 7.0, n_cluster)
    cluster = z0 + 0.08 * np.tanh(s)
    cluster = cluster[(cluster > 0.0) & (cluster < 1.0)]
    return np.unique(np.concatenate([np.linspace(0.0, 1.0, n_uniform), cluster]))

def solve_base_state(a: float, Gamma: float, kappa: float, tol: float = 1e-6, logger: Optional[Logger] = None) -> BaseState:
    if logger: logger.log("Solving steady base state ...")
    z = build_initial_mesh(a, Gamma)
    lam0, U0 = diffusionless_profile(z, a, Gamma)
    dlam0 = np.gradient(lam0, z, edge_order=2)
    y0 = np.vstack([lam0, dlam0, U0])
    def ode(z_, y):
        lam = np.clip(y[0], 0.0, 1.2)
        q = y[1]
        Uz = (1.0 - z_) * reciprocal_viscosity(lam, a)
        qz = (Gamma * lam * Uz - (1.0 - lam)) / kappa
        return np.vstack([q, qz, Uz])
    def bc(ya, yb):
        return np.array([ya[1], yb[1], ya[2]])
    sol = solve_bvp(ode, bc, z, y0, tol=tol, max_nodes=200000)
    z_dense = np.linspace(0.0, 1.0, 2001)
    y_dense = sol.sol(z_dense)
    if logger: logger.log(f"Base state done: success={sol.success}, message='{sol.message}'")
    return BaseState(z=z_dense, Lambda=np.clip(y_dense[0],0.0,1.05), U=y_dense[2], success=sol.success, message=sol.message)

def cheb_lobatto_ascending(M: int) -> tuple[np.ndarray, np.ndarray]:
    n = M - 1
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.ones(n + 1); c[0]=2.0; c[-1]=2.0; c = c*((-1.0)**np.arange(n+1))
    X = np.tile(x, (n+1,1)); dX = X - X.T
    D = np.outer(c,1.0/c)/(dX + np.eye(n+1)); D = D - np.diag(np.sum(D,axis=1))
    P = np.eye(M)[::-1]
    D = -2.0 * (P @ D @ P)
    x = x[::-1]
    zeta = 0.5*(x+1.0)
    return zeta, D

def clencurt_weights_01(M: int) -> np.ndarray:
    N = M - 1
    theta = np.pi * np.arange(M) / N
    w = np.zeros(M); ii = np.arange(1,N); v = np.ones(N-1)
    if N % 2 == 0:
        w[0] = 1.0/(N*N-1.0); w[-1]=w[0]
        for k in range(1,N//2):
            v -= 2.0*np.cos(2.0*k*theta[ii])/(4.0*k*k-1.0)
        v -= np.cos(N*theta[ii])/(N*N-1.0)
    else:
        w[0]=1.0/(N*N); w[-1]=w[0]
        for k in range(1,(N-1)//2+1):
            v -= 2.0*np.cos(2.0*k*theta[ii])/(4.0*k*k-1.0)
    w[ii] = 2.0*v/N
    return 0.5*w[::-1]

def build_integration_operator(D: np.ndarray) -> np.ndarray:
    M = D.shape[0]
    A = D.copy()
    A[0,:]=0.0; A[0,0]=1.0
    P0=np.eye(M); P0[0,:]=0.0
    return np.linalg.solve(A, P0)

def fourier_wavenumbers(N: int, L: float) -> np.ndarray:
    return 2.0*np.pi*np.fft.fftfreq(N, d=L/N)

class FourierDiff1D:
    def __init__(self, N: int, L: float, dealias: bool = True):
        self.N=N; self.L=L; self.k = fourier_wavenumbers(N,L)
        self.dealias = dealias
        self.filter = np.ones(N)
        if dealias:
            cutoff = N//3
            kk = np.fft.fftfreq(N)*N
            self.filter[np.abs(kk) >= cutoff] = 0.0
        self.kmax = np.max(np.abs(self.k))
    def spectral_filter(self, f: np.ndarray, axis: int = 0) -> np.ndarray:
        fhat = np.fft.fft(f, axis=axis)
        shape=[1]*f.ndim; shape[axis]=self.N
        g = np.fft.ifft(self.filter.reshape(shape)*fhat, axis=axis)
        return g.real if np.isrealobj(f) else g
    def dx(self, f: np.ndarray, axis: int = 0) -> np.ndarray:
        fhat = np.fft.fft(f, axis=axis)
        if self.dealias:
            shape=[1]*f.ndim; shape[axis]=self.N
            fhat = self.filter.reshape(shape)*fhat
        shape=[1]*f.ndim; shape[axis]=self.N
        deriv=np.fft.ifft((1j*self.k).reshape(shape)*fhat,axis=axis)
        return deriv.real if np.isrealobj(f) else deriv

@dataclass
class NonlinearSolution:
    x: np.ndarray
    zeta: np.ndarray
    times: np.ndarray
    h: np.ndarray
    lambda_full: np.ndarray
    mean_square: np.ndarray
    base_state: BaseState
    snapshot_indices: np.ndarray

class BalmforthNonlinearSSP2:
    def __init__(self, *, a: float, Gamma: float, T: float, kappa: float, k_wave: float,
                 N: int=64, M: int=64, amplitude: float=5e-4, logger: Optional[Logger]=None):
        self.a=float(a); self.Gamma=float(Gamma); self.T=float(T); self.kappa=float(kappa); self.k_wave=float(k_wave)
        self.L = 2*np.pi/self.k_wave
        self.N=int(N); self.M=int(M); self.amplitude=float(amplitude)
        self.logger=logger or Logger(enabled=False)
        self.x=np.linspace(0.0,self.L,self.N,endpoint=False)
        self.fd=FourierDiff1D(self.N,self.L,dealias=True)
        self.zeta,self.Dz = cheb_lobatto_ascending(self.M)
        self.Dzz = self.Dz @ self.Dz
        self.wz = clencurt_weights_01(self.M)
        self.Iz = build_integration_operator(self.Dz)
        self.one_minus_zeta = 1.0 - self.zeta
        self.min_dz = float(np.min(np.diff(self.zeta)))
        D=self.Dz
        A = np.array([[D[0,0],D[0,-1]],[D[-1,0],D[-1,-1]]],dtype=float)
        B = -np.vstack([D[0,1:-1], D[-1,1:-1]])
        bc_reconstruct = np.linalg.solve(A,B)
        self.Bmap=np.zeros((self.M,self.M-2))
        self.Bmap[0,:]=bc_reconstruct[0]; self.Bmap[1:-1,:]=np.eye(self.M-2); self.Bmap[-1,:]=bc_reconstruct[1]
        self.Lint = self.Dzz[1:-1,:] @ self.Bmap
        self.L_evals, self.L_evecs = eig(self.Lint)
        self.L_evecs_inv = np.linalg.inv(self.L_evecs)
        self.base = solve_base_state(self.a,self.Gamma,self.kappa,logger=self.logger)
        if not self.base.success:
            raise RuntimeError(self.base.message)
        self.base_on_grid = np.interp(self.zeta, self.base.z, self.base.Lambda)
    def assemble_lambda_full(self, lam_int):
        return lam_int @ self.Bmap.T
    def filter_x(self, arr):
        return self.fd.spectral_filter(arr, axis=0)
    def compute_u(self, h, lam_full, hx):
        A = reciprocal_viscosity(lam_full, self.a)
        tau_b = h * (1.0 - hx)
        uz = tau_b[:,None] * self.one_minus_zeta[None,:] * A
        uzeta = h[:,None] * uz
        u = uzeta @ self.Iz.T
        q = h * (u @ self.wz)
        return self.filter_x(u), self.filter_x(uz), self.filter_x(q)
    def compute_hnu(self, h, u, qx):
        hu_x = self.fd.dx(self.filter_x(h[:,None]*u), axis=0)
        rhs = qx[:,None] - hu_x
        hnu = rhs @ self.Iz.T
        return self.filter_x(hnu)
    def explicit_rhs(self, h, lam_int):
        h = self.filter_x(h)
        lam_int = self.filter_x(lam_int)
        lam_full = self.assemble_lambda_full(lam_int)
        lam_full = np.clip(self.filter_x(lam_full), -0.02, 1.02)
        hx = self.fd.dx(h)
        u, uz, q = self.compute_u(h, lam_full, hx)
        qx = self.fd.dx(q)
        h_t = -qx
        hnu = self.compute_hnu(h, u, qx)
        lam_x = self.fd.dx(lam_full, axis=0)
        lam_zeta = lam_full @ self.Dz.T
        rhs = ((1.0 - lam_full) - self.Gamma*uz*lam_full)/self.T
        rhs -= u * lam_x
        rhs -= (hnu / np.maximum(h[:,None],1e-10)) * lam_zeta
        rhs = self.filter_x(rhs)
        return h_t, rhs[:,1:-1], float(np.max(np.abs(u))), float(np.max(np.abs(hnu/np.maximum(h[:,None],1e-10))))
    def exact_diffuse_step(self, lam_int, h_frozen, dt):
        """Exact solve of the frozen-h vertical diffusion subproblem.

        For each x-node, with h frozen, the interior lambda variables satisfy
            d lambda_int / dt = alpha(x) * Lint * lambda_int,
        where alpha(x)=kappa/(T h(x)^2). Because Lint has been diagonalized once,
        this substep is applied exactly by exponentiating the eigenvalues.
        """
        if dt == 0.0:
            return lam_int.copy()
        alpha = self.kappa / (self.T * np.maximum(h_frozen, 1e-8)**2)
        tmp = (self.L_evecs_inv @ lam_int.T).T
        tmp = tmp * np.exp(dt * alpha[:, None] * self.L_evals[None, :])
        out = (self.L_evecs @ tmp.T).T.real
        return out
    def initial_state(self):
        h0 = 1.0 + self.amplitude*np.sin(self.k_wave*self.x)
        lam_int0 = np.tile(self.base_on_grid[1:-1], (self.N,1))
        return h0, lam_int0
    def suggested_dt(self, h, lam_int, user_dt):
        # robust explicit spectral CFL estimate for Fourier advection + mapped-zeta advection
        _, _, umax, numax = self.explicit_rhs(h, lam_int)
        dt_x = 0.08 / max(self.fd.kmax * max(umax, 1e-8), 1e-8)
        dt_z = 0.30 * self.min_dz / max(numax, 1e-8)
        dt_use = min(user_dt, dt_x, dt_z, 1.0)
        return max(dt_use, 1e-5), dt_x, dt_z, umax, numax
    def solve(self, *, t_end: float, dt_max: float, n_eval: int=31, checkpoint_dir: Optional[Path]=None, progress_every: int=100):
        h, lam_int = self.initial_state()
        times_out = np.linspace(0.0, t_end, n_eval)
        h_hist = np.full((n_eval, self.N), np.nan)
        lam_hist = np.full((n_eval, self.N, self.M), np.nan)
        ms = np.full(n_eval, np.nan)
        h_hist[0] = h
        lam_hist[0] = self.assemble_lambda_full(lam_int)
        ms[0] = np.mean((h - 1.0) ** 2)
        snapshot_indices = np.array(sorted(set(np.linspace(0, n_eval - 1, min(7, n_eval), dtype=int))))
        if self.logger:
            self.logger.log(f"Stored initial state at physical time t=0; <(h-1)^2>={ms[0]:.6e}")

        t = 0.0
        out_index = 1
        step = 0
        while t < t_end - 1e-14:
            dt, dt_x, dt_z, umax, numax = self.suggested_dt(h, lam_int, dt_max)
            if t + dt > t_end:
                dt = t_end - t

            # --- Second-order SSP-IMEX / Strang-split step ---
            # 1) stiff vertical diffusion half-step (exact with h frozen at time level n)
            lam_a = self.exact_diffuse_step(lam_int, h, 0.5 * dt)
            lam_a = self.filter_x(np.clip(lam_a, -0.02, 1.02))
            h_a = h.copy()

            # 2) explicit full-step using SSPRK2 (Heun / TVD RK2) on the non-stiff operator
            h_rhs0, lam_rhs0, _, _ = self.explicit_rhs(h_a, lam_a)
            h1 = self.filter_x(np.maximum(h_a + dt * h_rhs0, 1e-6))
            lam1 = self.filter_x(np.clip(lam_a + dt * lam_rhs0, -0.02, 1.02))

            h_rhs1, lam_rhs1, _, _ = self.explicit_rhs(h1, lam1)
            h_b = self.filter_x(np.maximum(0.5 * h_a + 0.5 * (h1 + dt * h_rhs1), 1e-6))
            lam_b = self.filter_x(np.clip(0.5 * lam_a + 0.5 * (lam1 + dt * lam_rhs1), -0.02, 1.02))

            # 3) stiff vertical diffusion half-step (exact with h frozen at n+1 explicit state)
            lam_new = self.exact_diffuse_step(lam_b, h_b, 0.5 * dt)
            lam_new = self.filter_x(np.clip(lam_new, -0.02, 1.02))
            h_new = h_b

            if (not np.all(np.isfinite(h_new))) or (not np.all(np.isfinite(lam_new))):
                raise RuntimeError(f"Non-finite values detected at physical time t={t+dt:.6g}. Reduce --dt-max.")
            if np.max(np.abs(h_new - 1.0)) > 2.0:
                raise RuntimeError(f"Solution blew up at physical time t={t+dt:.6g}. Use smaller --dt-max (current adaptive dt={dt:.3g}).")

            h, lam_int = h_new, lam_new
            t += dt
            step += 1
            if (step % progress_every == 0) or step == 1 or abs(t - t_end) < 1e-12:
                ms_now = np.mean((h - 1.0) ** 2)
                self.logger.log(
                    f"Step {step}: physical time t={t:.6g}/{t_end:.6g}, <(h-1)^2>={ms_now:.6e}, "
                    f"dt={dt:.3g}, dt_x={dt_x:.3g}, dt_z={dt_z:.3g}, umax={umax:.3g}"
                )

            while out_index < n_eval and t >= times_out[out_index] - 1e-14:
                lam_full = self.assemble_lambda_full(lam_int)
                h_hist[out_index] = h
                lam_hist[out_index] = lam_full
                ms[out_index] = np.mean((h - 1.0) ** 2)
                self.logger.log(f"Stored output at physical time t={times_out[out_index]:.6g}; <(h-1)^2>={ms[out_index]:.6e}")
                if checkpoint_dir is not None:
                    checkpoint_dir.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(
                        checkpoint_dir / 'solution_checkpoint.npz',
                        x=self.x,
                        zeta=self.zeta,
                        times=times_out[:out_index + 1],
                        h=h_hist[:out_index + 1],
                        lambda_full=lam_hist[:out_index + 1],
                        mean_square=ms[:out_index + 1],
                        snapshot_indices=snapshot_indices,
                    )
                out_index += 1

        self.logger.log(f"Time marching complete at physical time t={t:.6g}.")
        return NonlinearSolution(self.x, self.zeta, times_out, h_hist, lam_hist, ms, self.base, snapshot_indices)

def estimate_wave_speed(times,h_hist,x,L,frac=0.35):
    finite = np.isfinite(h_hist).all(axis=1)
    times = times[finite]; h_hist=h_hist[finite]
    n0 = max(1, int((1.0-frac)*len(times)))
    xs = x[np.nanargmax(h_hist[n0:], axis=1)]
    xs = np.unwrap(2.0*np.pi*xs/L)*L/(2.0*np.pi)
    coeff = np.polyfit(times[n0:], xs, 1)
    return float(coeff[0])

def reconstruct_u_nu(solver, h, lam_full):
    hx = solver.fd.dx(h)
    u, uz, q = solver.compute_u(h, lam_full, hx)
    qx = solver.fd.dx(q)
    hnu = solver.compute_hnu(h,u,qx)
    return u, hnu/np.maximum(h[:,None],1e-10), hx, uz

def streamfunction_wave_frame(solver,h,lam_full,c_wave):
    u,_,_,_ = reconstruct_u_nu(solver,h,lam_full)
    psi_zeta = h[:,None]*(u-c_wave)
    psi = psi_zeta @ solver.Iz.T
    return psi

def center_periodic_profile(x, y, L):
    i_max = int(np.nanargmax(y))
    x_shift = x[i_max]
    x_center = ((x - x_shift + 0.5*L)%L)-0.5*L
    order=np.argsort(x_center)
    return x_center[order], y[order], order

def make_fig9_like(solution, solver, outdir: Path, tag: str):
    outdir.mkdir(parents=True, exist_ok=True)
    finite = np.isfinite(solution.mean_square)
    times = solution.times[finite]
    mean_square = solution.mean_square[finite]
    h_hist = solution.h[finite]
    c_ref = estimate_wave_speed(times, h_hist, solution.x, solver.L)
    plt.figure(figsize=(6,4))
    plt.semilogy(times, mean_square, lw=2)
    plt.xlabel('t'); plt.ylabel(r'$\langle (h-1)^2\rangle$')
    # mark snapshot times
    snap_ids = [i for i in solution.snapshot_indices if i < len(solution.times) and finite[i]]
    snap_times = solution.times[snap_ids]
    plt.plot(snap_times, solution.mean_square[snap_ids], 'r*', ms=8)
    plt.tight_layout(); plt.savefig(outdir/f'{tag}_fig9b_timeseries.png', dpi=200); plt.close()

    plt.figure(figsize=(7,4))
    for j, idx in enumerate(snap_ids):
        xcen, y, _ = center_periodic_profile(solution.x - c_ref*solution.times[idx], solution.h[idx], solver.L)
        alpha = 0.35 + 0.65*j/max(1,len(snap_ids)-1)
        plt.plot(solver.k_wave*xcen, y, lw=1.4, alpha=alpha, label=f't={solution.times[idx]:.0f}')
    plt.xlabel(r'$k(x-c_n t)$'); plt.ylabel('h')
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout(); plt.savefig(outdir/f'{tag}_fig9c_h_profiles.png', dpi=200); plt.close()

def make_fig10_like(solution, solver, outdir: Path, tag: str):
    finite = np.isfinite(solution.mean_square)
    last = np.where(finite)[0][-1]
    h = solution.h[last]; lam = solution.lambda_full[last]
    c_n = estimate_wave_speed(solution.times[finite], solution.h[finite], solution.x, solver.L)
    psi = streamfunction_wave_frame(solver,h,lam,c_n)
    x_center, h_ord, order = center_periodic_profile(solution.x, h, solver.L)
    lam = lam[order]; psi = psi[order]
    X=np.tile(x_center[:,None], (1,solver.M)); Z=h_ord[:,None]*solver.zeta[None,:]
    plt.figure(figsize=(7,4.5))
    pc=plt.pcolormesh(solver.k_wave*X, Z, lam, shading='gouraud', cmap='viridis')
    plt.colorbar(pc,label=r'$\lambda$')
    levels=np.linspace(np.nanmin(psi), np.nanmax(psi), 13)
    plt.contour(solver.k_wave*X, Z, psi, levels=levels, colors='r', linewidths=0.6)
    plt.xlabel(r'$k(x-c_n t)$'); plt.ylabel('z'); plt.ylim(0.0,1.05*np.nanmax(h_ord))
    plt.title(f'final stored time t={solution.times[last]:.0f}')
    plt.tight_layout(); plt.savefig(outdir/f'{tag}_fig10_final_wave.png', dpi=200); plt.close()

def parse_args(argv=None):
    p=argparse.ArgumentParser(description='Second-order SSP-IMEX nonlinear solver for Balmforth (2025), Section 5 / Appendix F.')
    p.add_argument('--a', type=float, default=0.2)
    p.add_argument('--Gamma', type=float, default=8.0)
    p.add_argument('--T', type=float, default=100.0)
    p.add_argument('--kappa', type=float, default=1e-4)
    p.add_argument('--k', type=float, default=0.1, dest='k_wave')
    p.add_argument('--N', type=int, default=64)
    p.add_argument('--M', type=int, default=64)
    p.add_argument('--amplitude', type=float, default=5e-4)
    p.add_argument('--t-end', type=float, default=3000.0)
    p.add_argument('--dt-max', type=float, default=0.1, help='Maximum allowed time step; code may reduce it adaptively.')
    p.add_argument('--n-eval', type=int, default=31)
    p.add_argument('--progress-every', type=int, default=200)
    p.add_argument('--output-dir', type=Path, default=Path('balmforth_nonlinear_out'))
    p.add_argument('--make-plots', action='store_true')
    p.add_argument('--quiet', action='store_true')
    p.add_argument('--paper-fig9a', action='store_true')
    p.add_argument('--paper-fig9b', action='store_true')
    return p.parse_args(argv)

def apply_presets(args):
    if args.paper_fig9a:
        args.a=1/5; args.Gamma=8.0; args.T=100.0; args.kappa=1e-4; args.k_wave=0.1
    if args.paper_fig9b:
        args.a=3/5; args.Gamma=8.0; args.T=100.0; args.kappa=1e-4; args.k_wave=0.1

def main(argv=None):
    args=parse_args(argv); apply_presets(args)
    logger=Logger(enabled=not args.quiet)
    logger.log("Launching corrected Balmforth nonlinear-wave solver.")
    logger.log("This version removes the unstable AB2 march, uses adaptive IMEX-Euler with de-aliasing, aborts on NaN/blow-up, and labels snapshot times in the plots.")
    solver=BalmforthNonlinearSSP2(a=args.a,Gamma=args.Gamma,T=args.T,kappa=args.kappa,k_wave=args.k_wave,N=args.N,M=args.M,amplitude=args.amplitude,logger=logger)
    logger.log(f"Parameters: a={args.a}, Gamma={args.Gamma}, T={args.T}, kappa={args.kappa}, k={args.k_wave}")
    logger.log(f"Domain length L={solver.L:.10g}, N={args.N}, M={args.M}, dt_max={args.dt_max}")
    try:
        sol=solver.solve(t_end=args.t_end, dt_max=args.dt_max, n_eval=args.n_eval, checkpoint_dir=args.output_dir, progress_every=args.progress_every)
    except RuntimeError as e:
        logger.log(f"ERROR: {e}")
        logger.log("No misleading plots were generated. Reduce --dt-max and rerun.")
        return 2
    c_n = estimate_wave_speed(sol.times, sol.h, sol.x, solver.L)
    logger.log(f"Estimated final wave speed c_n ≈ {c_n:.8f}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_dir/'solution.npz', x=sol.x, zeta=sol.zeta, times=sol.times, h=sol.h, lambda_full=sol.lambda_full, mean_square=sol.mean_square, snapshot_indices=sol.snapshot_indices)
    logger.log(f"Final solution written to {args.output_dir/'solution.npz'}")
    if args.make_plots:
        tag=f'a_{args.a:.6f}_k_{args.k_wave:.6f}'.replace('.','p')
        make_fig9_like(sol, solver, args.output_dir, tag)
        make_fig10_like(sol, solver, args.output_dir, tag)
        logger.log(f"Plots written to {args.output_dir}")
    logger.log("Run complete.")
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
