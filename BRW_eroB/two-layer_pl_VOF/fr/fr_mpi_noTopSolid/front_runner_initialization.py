#!/usr/bin/env python3
"""Analytical inspection counterpart of front_runner_ic.h (not a PDE solver).
Uses exactly the same piecewise-linear base velocity and its exact integral.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np


class InitialField:
    def __init__(self, z, u, *, amplitude, wavelength, center, depth, ceiling):
        self.z=np.asarray(z,dtype=float); self.u0=np.asarray(u,dtype=float)
        if self.z.ndim!=1 or self.z.shape!=self.u0.shape or np.any(np.diff(self.z)<=0):
            raise ValueError('The base profile must be strictly ordered 1D arrays')
        self.prefix=np.r_[0.,np.cumsum(.5*(self.u0[:-1]+self.u0[1:])*np.diff(self.z))]
        self.a=float(amplitude);self.lam=float(wavelength);self.xc=float(center)
        self.H=float(depth);self.yc=float(ceiling)
        self.Us=float(np.interp(self.H,self.z,self.u0)); self.Q0=float(self.integral(self.H))
        self.Qt=self.Q0+self.Us*(self.yc-self.H)

    def integral(self, y):
        y=np.asarray(y,dtype=float)
        j=np.clip(np.searchsorted(self.z,y,side='right')-1,0,len(self.z)-2)
        dz=y-self.z[j]
        v=self.prefix[j]+self.u0[j]*dz+.5*(self.u0[j+1]-self.u0[j])/(self.z[j+1]-self.z[j])*dz**2
        return np.where(y>=self.z[-1],self.prefix[-1]+self.u0[-1]*(y-self.z[-1]),
                        np.where(y<=self.z[0],self.u0[0]*(y-self.z[0]),v))

    def scale(self, x):
        x=np.asarray(x,dtype=float);r=(x-self.xc)/self.lam
        active=np.abs(r)<.25
        s=1.+np.where(active,self.a*np.cos(2*np.pi*r),0.)
        sx=np.where(active,-self.a*2*np.pi/self.lam*np.sin(2*np.pi*r),0.)
        return s,sx

    def evaluate(self, x, y):
        x,y=np.broadcast_arrays(np.asarray(x,dtype=float),np.asarray(y,dtype=float))
        s,sx=self.scale(x);rs=np.sqrt(s);h=s*self.H
        z0=y/s;U=np.interp(z0,self.z,self.u0);P=self.integral(z0)
        u=rs*U;w=-rs*sx*(1.5*P-z0*U);psi=s*rs*P
        air=(y>h)&(y<self.yc)
        d=self.yc-h;hx=self.H*sx;r=(y-h)/d
        # Values away from the air interval are masked, but keep polynomials bounded.
        r=np.clip(r,0.,1.)
        A=1-3*r*r+2*r**3;Ai=r-r**3+.5*r**4
        C=30*r*r*(1-r)**2;Ci=10*r**3-15*r**4+6*r**5
        Q=s*rs*self.Q0;Qx=1.5*rs*sx*self.Q0
        du=self.Us*(rs-1);dux=self.Us*sx/(2*rs)
        B=(self.Qt-Q)/d-self.Us-.5*du
        Bx=-Qx/d+(self.Qt-Q)*hx/d**2-.5*dux
        F=self.Us*r+du*Ai+B*Ci;ua=self.Us+du*A+B*C
        wa=-Qx+hx*F-d*(dux*Ai+Bx*Ci)-hx*(r-1)*ua
        u=np.where(air,ua,u);w=np.where(air,wa,w);psi=np.where(air,Q+d*F,psi)
        u=np.where(y<=0.,0.,np.where(y>=self.yc,self.Us,u))
        w=np.where((y<=0.)|(y>=self.yc),0.,w)
        psi=np.where(y<=0.,0.,np.where(y>=self.yc,self.Qt+self.Us*(y-self.yc),psi))
        return u,w,psi


def from_case(c, z, velocity):
    lam=c.front_runner.wavelength_slope_scaled/c.slope
    return InitialField(z,velocity,amplitude=c.lower_depth_amplitude,
        wavelength=lam,center=c.front_runner.center_wavelengths*lam,
        depth=c.liquid_depth,ceiling=c.ceiling)


def write_initial_tables(output_dir:Path,c,z,velocity):
    f=from_case(c,z,velocity);lam=f.lam
    xx=np.unique(np.r_[0.,np.linspace(max(0.,f.xc-.5*lam),min(c.lx,f.xc+.5*lam),801),c.lx])
    s,sx=f.scale(xx);qlo0=float(f.integral(1.));qup0=f.Q0-qlo0
    h1=s;h2=c.depth_ratio*s;ql=qlo0*s**1.5;qu=qup0*s**1.5
    us,ws,_=f.evaluate(xx,h1+h2)
    table=np.column_stack((xx,c.slope*xx,s,sx,h1,h2,h1+h2,ql,qu,ql/h1,qu/h2,
        c.froude*ql/h1/np.sqrt(h1),c.froude*qu/h2/np.sqrt(h2),
        us,ws,f.Qt-ql-qu,np.full_like(xx,f.Qt)))
    np.savetxt(output_dir/'initial_columns.tsv',table,delimiter='\t',fmt='%.16e',
        header='x_star\tS0_x_star\tdepth_factor\tdfactor_dx_star\th_lower\th_upper\th_total\tq_lower\tq_upper\tUbar_lower\tUbar_upper\tFr_lower\tFr_upper\tu_surface\tw_surface\tq_air\tq_total',comments='')
    selected=np.array([0.,f.xc-.2*lam,f.xc,f.xc+.2*lam,c.lx])
    rows=[]
    for x in selected:
        fac,_=f.scale(x)
        yy=np.unique(np.r_[np.linspace(0.,c.ceiling,601),float(fac),float(fac)*f.H])
        u,w,psi=f.evaluate(x,yy)
        rows.append(np.column_stack((np.full_like(yy,x),yy,u,w,psi)))
    np.savetxt(output_dir/'initial_velocity_slices.tsv',np.vstack(rows),delimiter='\t',fmt='%.16e',
               header='x_star\tz_star\tu_star\tw_star\tpsi_star',comments='')
    info=dict(convention='R_eta,I = eta_lower,I / eta_upper,I; inherited, not inverted',
      length_scale='H_lower',velocity_scale='Ubar_lower',lx_star=c.lx,lx_slope_scaled=c.slope*c.lx,
      wavelength_star=lam,wavelength_slope_scaled=c.slope*lam,center_star=f.xc,
      support_star=[f.xc-lam/4,f.xc+lam/4],amplitude=f.a,
      ceiling_star=c.ceiling,delta_min=c.delta_min,cells_per_lower_depth=1/c.delta_min,
      q_lower_base=qlo0,q_upper_base=qup0,q_liquid_base=f.Q0,u_surface_base=f.Us,q_total=f.Qt,
      excess_lower_volume=f.a*lam/np.pi,excess_upper_volume=c.depth_ratio*f.a*lam/np.pi,
      Fr_lower_base_interpolated=c.froude*qlo0,
      Fr_upper_base_interpolated=c.froude*qup0/c.depth_ratio**1.5,
      lower_mean_normalization_error=qlo0-1.,
      initial_mean_lower_depth=1.+f.a*lam/(np.pi*c.lx),
      initial_mean_total_depth=f.H*(1.+f.a*lam/(np.pi*c.lx)),
      min_air_streamwise_velocity=float(np.min(rows[2][rows[2][:,1]>=f.H*(1+f.a),2])),
      sigma_internal=c.sigma_internal,sigma_free=c.sigma_free,
      verification_scope='Analytical initial field; not a nonlinear DNS validation')
    (output_dir/'front_runner_initial_summary.json').write_text(json.dumps(info,indent=2)+'\n')
