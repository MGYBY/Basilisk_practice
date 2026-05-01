/**
# Two-phase VOF solver for thixotropic liquid--air flows

This header extends Basilisk's two-phase VOF framework by

1. adding a thixotropic structural field lambda,
2. transporting lambda consistently with the VOF phase fraction using
   the tracer machinery, and
3. computing a spatially-varying liquid viscosity from the rheology

    mu_liquid = 1 / [ (1 - lambda) * ( (1 - lambda)(1 - a) + a ) ]

in dimensionless form.

The thixotropic transport equation implemented here is

  d_t lambda + div(u lambda) = (So/T) [ (1 - lambda) - Gamma lambda gammadot ]
                             + (So/T) kappa nabla^2 lambda

with the advection handled through VOF-consistent tracer transport,
source integration by SSPRK2, and diffusion through diffusion.h.
*/

#include "vof.h"
#include "diffusion.h"
#include "fractions.h"

#ifndef THIXO_FERR
#define THIXO_FERR 1e-12
#endif

/* VOF-interface safeguards for lambda. */
#ifndef THIXO_LAMBDA_RECON_FMIN
#define THIXO_LAMBDA_RECON_FMIN 1e-8
#endif
#ifndef THIXO_LAMBDA_SOURCE_FMIN
#define THIXO_LAMBDA_SOURCE_FMIN 5e-2
#endif
#ifndef THIXO_GAMMADOT_FMIN
#define THIXO_GAMMADOT_FMIN 5e-2
#endif
#ifndef THIXO_GAMMADOT_NEIGHBOR_FMIN
#define THIXO_GAMMADOT_NEIGHBOR_FMIN 5e-2
#endif
#ifndef THIXO_LAMBDA_PLOT_FMIN
#define THIXO_LAMBDA_PLOT_FMIN 0.50
#endif
#ifndef THIXO_DIFF_FMIN
#define THIXO_DIFF_FMIN 1e-6
#endif
#ifndef THIXO_ENFORCE_DIFFUSIVE_LAMBDA_MASS
#define THIXO_ENFORCE_DIFFUSIVE_LAMBDA_MASS 1
#endif
#ifndef THIXO_LAMBDA_MIN
#define THIXO_LAMBDA_MIN 0.
#endif
#ifndef THIXO_LAMBDA_MAX
#define THIXO_LAMBDA_MAX 1.
#endif
#ifndef THIXO_MU_MAX
#define THIXO_MU_MAX (3.0*(1.0/0.05)) // roughly estimated from Balmforth as 1/\epsilon
#endif
#ifndef THIXO_MU_AIR_MIN
#define THIXO_MU_AIR_MIN 1e-12
#endif
#ifndef THIXO_MU_DENOM_MIN
#define THIXO_MU_DENOM_MIN (1./THIXO_MU_MAX)
#endif

scalar f[], * interfaces = {f};

/* Conservative tracer attached to f. Physically, flambda = f*lambda. */
scalar flambda[];

/* Cell-centred diagnostic fields. */
scalar lambda_thixo[];
scalar gammadot_thixo[];
scalar mu_liquid_thixo[];

/* Standard phase properties (dimensionless). */
double rho1 = 1., mu1 = 1.;
double rho2 = 1., mu2 = 0.;

/* Thixotropic parameters (dimensionless). */
double thixo_T = 1.;
double thixo_Gamma = 0.;
double thixo_kappa = 0.;
double thixo_a = 0.2;
double thixo_So = 0.05;
double thixo_FrV = 1.;
double thixo_lambda_air = 0.;

face vector alphav[];
scalar rhov[];

static inline double thixo_clamp_lambda (double l)
{
  return max (THIXO_LAMBDA_MIN, min (THIXO_LAMBDA_MAX, l));
}

static inline double thixo_source_prefactor()
{
  return thixo_So/(thixo_T + 1e-30);
}

static inline double thixo_visc_prefactor()
{
  return thixo_So/(sq(thixo_FrV) + 1e-30);
}

/*
  Regularised liquid viscosity.

  The theoretical rheology diverges as lambda -> 1. Numerically, that makes the
  implicit viscous solve extremely stiff and can stall the very first timestep.
  We therefore evaluate the viscosity with a floor on the denominator, which is
  equivalent to capping the liquid viscosity at THIXO_MU_MAX.
*/
static inline double thixo_mu_fluid_from_lambda (double l)
{
  double lc = thixo_clamp_lambda (l);
  double om = 1. - lc;
  double denom = om*(om*(1. - thixo_a) + thixo_a);
  denom = max (denom, THIXO_MU_DENOM_MIN);
  return min (1./denom, THIXO_MU_MAX);
}

static inline double thixo_mu_mix_harmonic (double mu_liq, double mu_air, double ff)
{
  double fcl = clamp (ff, 0., 1.);
  return 1./(fcl/max(mu_liq, THIXO_MU_AIR_MIN) +
             (1. - fcl)/max(mu_air, THIXO_MU_AIR_MIN));
}

static inline void thixo_sync_lambda_from_tracer_cells()
{
  foreach() {
    double fc = clamp (f[], 0., 1.);
    lambda_thixo[] = (fc > THIXO_LAMBDA_RECON_FMIN ? thixo_clamp_lambda(flambda[]/fc) : thixo_lambda_air);
    mu_liquid_thixo[] = thixo_mu_fluid_from_lambda (lambda_thixo[]);
  }
}

static inline void thixo_sync_lambda_from_tracer()
{
  thixo_sync_lambda_from_tracer_cells();
  boundary ({lambda_thixo, mu_liquid_thixo});
}

static inline void thixo_sync_tracer_from_lambda_cells()
{
  foreach() {
    double fc = clamp (f[], 0., 1.);
    lambda_thixo[] = (fc > THIXO_FERR ? thixo_clamp_lambda(lambda_thixo[]) : thixo_lambda_air);
    flambda[] = fc*lambda_thixo[];
    mu_liquid_thixo[] = thixo_mu_fluid_from_lambda (lambda_thixo[]);
  }
}

static inline void thixo_sync_tracer_from_lambda()
{
  thixo_sync_tracer_from_lambda_cells();
  /* Avoid boundary() on flambda here: on TREE+EMBED this can trigger
     VOF-concentration prolongation/boundary logic too early during init.
     flambda is synchronised cell-wise and only lambda/mu diagnostics need
     ghost updates here. */
  boundary ({lambda_thixo, mu_liquid_thixo});
}


/* Liquid-side lambda inventory. The conserved transported quantity is
   flambda = f*lambda. Advection should conserve this inventory up to VOF
   errors. Reaction changes it physically; diffusion with zero normal flux at
   the free surface should not change it. */
static inline double thixo_total_flambda()
{
  double m = 0.;
  foreach (reduction(+:m))
    m += flambda[]*dv();
  return m;
}

static inline double thixo_total_lambda_from_intrinsic()
{
  double m = 0.;
  foreach (reduction(+:m)) {
    double fc = clamp (f[], 0., 1.);
    m += fc*thixo_clamp_lambda(lambda_thixo[])*dv();
  }
  return m;
}

static inline void thixo_correct_lambda_mass (double target_mass)
{
#if THIXO_ENFORCE_DIFFUSIVE_LAMBDA_MASS
  thixo_sync_tracer_from_lambda_cells();
  double current_mass = thixo_total_flambda();
  double liquid_volume = 0.;
  foreach (reduction(+:liquid_volume)) {
    double fc = clamp (f[], 0., 1.);
    if (fc > THIXO_DIFF_FMIN)
      liquid_volume += fc*dv();
  }
  if (liquid_volume > 0.) {
    double dlambda = (target_mass - current_mass)/liquid_volume;
    foreach() {
      double fc = clamp (f[], 0., 1.);
      if (fc > THIXO_DIFF_FMIN)
        lambda_thixo[] = thixo_clamp_lambda (lambda_thixo[] + dlambda);
      else
        lambda_thixo[] = thixo_lambda_air;
      flambda[] = fc*lambda_thixo[];
      mu_liquid_thixo[] = thixo_mu_fluid_from_lambda (lambda_thixo[]);
    }
    boundary ({lambda_thixo, mu_liquid_thixo});
  }
#endif
}

/*
  Interface-safe velocity derivatives for the thixotropic source term.

  The structural parameter is a liquid material variable. Near the VOF
  interface, centered differences can sample air cells and convert a numerical
  air/liquid velocity jump into artificial microstructure destruction. We
  therefore update the source only in mostly-liquid cells and use one-sided
  liquid derivatives when one neighbour lies outside the liquid.
*/
static inline bool thixo_source_cell (Point point)
{
  return clamp(f[],0.,1.) > THIXO_GAMMADOT_FMIN;
}

static inline double thixo_limited_derivative (double qc, double qp, double qm,
                                              double fp, double fm, double Delta)
{
  bool has_p = clamp(fp,0.,1.) > THIXO_GAMMADOT_NEIGHBOR_FMIN;
  bool has_m = clamp(fm,0.,1.) > THIXO_GAMMADOT_NEIGHBOR_FMIN;

  if (has_p && has_m)
    return (qp - qm)/(2.*Delta);
  if (has_p)
    return (qp - qc)/Delta;
  if (has_m)
    return (qc - qm)/Delta;
  return 0.;
}

static inline double thixo_face_lambda_from_conservative (Point point)
{
  double f0 = clamp(f[],   0., 1.);
  double f1 = clamp(f[-1], 0., 1.);
  double fs = f0 + f1;
  if (fs > 2.*THIXO_LAMBDA_RECON_FMIN)
    return thixo_clamp_lambda ((flambda[] + flambda[-1])/fs);
  return thixo_lambda_air;
}

/*
  Cell-centred strain-rate measure following the spirit of the tested
  stencil used in two-phasePL.h. The returned D2 is sqrt(D_ij D_ij),
  so that gammadot = sqrt(2)*D2.
*/
static inline double thixo_D2_point (Point point)
{
#if dimension == 1
  double dxx = thixo_limited_derivative (u.x[], u.x[1], u.x[-1],
                                         f[1], f[-1], Delta);
  return fabs(dxx);
#elif dimension == 2
// TODO: check the consistency here
  double dxx = thixo_limited_derivative (u.x[], u.x[1,0], u.x[-1,0],
                                         f[1,0], f[-1,0], Delta);
  double dyy = thixo_limited_derivative (u.y[], u.y[0,1], u.y[0,-1],
                                         f[0,1], f[0,-1], Delta);
  double duxdy = thixo_limited_derivative (u.x[], u.x[0,1], u.x[0,-1],
                                           f[0,1], f[0,-1], Delta);
  double duydx = thixo_limited_derivative (u.y[], u.y[1,0], u.y[-1,0],
                                           f[1,0], f[-1,0], Delta);
  double dxy = 0.5*(duxdy + duydx);
  return sqrt (sq(dxx) + sq(dyy) + 2.*sq(dxy));
#else
  double dxx = thixo_limited_derivative (u.x[], u.x[1,0,0], u.x[-1,0,0],
                                         f[1,0,0], f[-1,0,0], Delta);
  double dyy = thixo_limited_derivative (u.y[], u.y[0,1,0], u.y[0,-1,0],
                                         f[0,1,0], f[0,-1,0], Delta);
  double dzz = thixo_limited_derivative (u.z[], u.z[0,0,1], u.z[0,0,-1],
                                         f[0,0,1], f[0,0,-1], Delta);
  double duxdy = thixo_limited_derivative (u.x[], u.x[0,1,0], u.x[0,-1,0],
                                           f[0,1,0], f[0,-1,0], Delta);
  double duydx = thixo_limited_derivative (u.y[], u.y[1,0,0], u.y[-1,0,0],
                                           f[1,0,0], f[-1,0,0], Delta);
  double duxdz = thixo_limited_derivative (u.x[], u.x[0,0,1], u.x[0,0,-1],
                                           f[0,0,1], f[0,0,-1], Delta);
  double duzdx = thixo_limited_derivative (u.z[], u.z[1,0,0], u.z[-1,0,0],
                                           f[1,0,0], f[-1,0,0], Delta);
  double duydz = thixo_limited_derivative (u.y[], u.y[0,0,1], u.y[0,0,-1],
                                           f[0,0,1], f[0,0,-1], Delta);
  double duzdy = thixo_limited_derivative (u.z[], u.z[0,1,0], u.z[0,-1,0],
                                           f[0,1,0], f[0,-1,0], Delta);
  double dxy = 0.5*(duxdy + duydx);
  double dxz = 0.5*(duxdz + duzdx);
  double dyz = 0.5*(duydz + duzdy);
  return sqrt (sq(dxx) + sq(dyy) + sq(dzz) +
               2.*sq(dxy) + 2.*sq(dxz) + 2.*sq(dyz));
#endif
}

static inline void thixo_update_gammadot()
{
  foreach() {
    if (thixo_source_cell (point)) {
      double D2 = thixo_D2_point (point);
      gammadot_thixo[] = sqrt(2.)*D2;
    }
    else
      gammadot_thixo[] = 0.;
  }
  boundary ({gammadot_thixo});
}

event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;

  if (mu1 || mu2)
    mu = new face vector;

  /* Register the conservative lambda tracer with the VOF field. */
  f.tracers = list_append (f.tracers, flambda);
  flambda.inverse = false;

#if TREE
  flambda.restriction = restriction_volume_average;
  flambda.refine = flambda.prolongation = vof_concentration_refine;
  // flambda.dirty = true;
  flambda.c = f;

  lambda_thixo.refine = lambda_thixo.prolongation = refine_embed_linear;
  lambda_thixo.restriction = restriction_volume_average;
  // lambda_thixo.dirty = true;

  gammadot_thixo.refine = gammadot_thixo.prolongation = refine_embed_linear;
  gammadot_thixo.restriction = restriction_volume_average;
  // gammadot_thixo.dirty = true;

  mu_liquid_thixo.refine = mu_liquid_thixo.prolongation = refine_embed_linear;
  mu_liquid_thixo.restriction = restriction_volume_average;
  // mu_liquid_thixo.dirty = true;
#endif
}

#ifndef rho
# define rho(ff) (clamp((ff),0.,1.)*(rho1 - rho2) + rho2)
#endif

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
#ifdef FILTERED
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
            2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
            f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else
  foreach()
    sf[] = (8.*f[] +
            4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
            2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
                f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
                f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
            f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
            f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif

#if TREE
  sf.prolongation = refine_bilinear;
  // sf.dirty = true;
#endif
}

event properties (i++)
{
  thixo_sync_lambda_from_tracer();

  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      double lf = thixo_face_lambda_from_conservative (point);
      double mu_liq = thixo_mu_fluid_from_lambda (lf);
      muv.x[] = fm.x[]*thixo_visc_prefactor()*
                thixo_mu_mix_harmonic (mu_liq, mu2, ff);
    }
  }

  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  // sf.dirty = true;
#endif
}

/*
  Operator-split update of the thixotropic field after VOF/tracer
  advection:
  1. recover lambda from flambda,
  2. explicit SSPRK2 reaction/source term,
  3. implicit diffusion through diffusion().
*/
event tracer_diffusion (i++, last)
{
  thixo_sync_lambda_from_tracer();
  thixo_update_gammadot();

  scalar lambda0[], lambda1[];
  double pref = thixo_source_prefactor();

  foreach() {
    lambda0[] = lambda_thixo[];
    if (clamp(f[],0.,1.) > THIXO_LAMBDA_SOURCE_FMIN) {
      double rhs0 = pref*((1. - lambda0[]) - thixo_Gamma*lambda0[]*gammadot_thixo[]);
      lambda1[] = thixo_clamp_lambda (lambda0[] + dt*rhs0);
    }
    else
      lambda1[] = lambda0[];
  }
  boundary ({lambda0, lambda1});

  foreach() {
    if (clamp(f[],0.,1.) > THIXO_LAMBDA_SOURCE_FMIN) {
      double rhs1 = pref*((1. - lambda1[]) - thixo_Gamma*lambda1[]*gammadot_thixo[]);
      lambda_thixo[] = thixo_clamp_lambda (0.5*(lambda0[] + lambda1[] + dt*rhs1));
    }
    else
      lambda_thixo[] = lambda0[];
  }
  boundary ({lambda_thixo});

  if (thixo_kappa > 0.) {
    /* Diffuse intrinsic lambda as a liquid-phase variable:

         f d(lambda)/dt = div( f kappa grad(lambda) ).

       The theta=f weighting makes the implicit solve conservative for the
       liquid inventory. The face coefficient is switched off at liquid-air
       faces to impose approximately zero flux of microstructure through the
       free surface. */
    double mass_before_diffusion = thixo_total_lambda_from_intrinsic();

    scalar theta_lambda[];
    face vector D[];
    foreach()
      theta_lambda[] = max (clamp(f[], 0., 1.), THIXO_DIFF_FMIN);
    boundary ({theta_lambda});

    foreach_face() {
      double f0 = clamp(f[],   0., 1.);
      double f1 = clamp(f[-1], 0., 1.);
      if (f0 > THIXO_DIFF_FMIN && f1 > THIXO_DIFF_FMIN) {
        double ff = 0.5*(f0 + f1);
        D.x[] = fm.x[]*(pref*thixo_kappa)*ff;
      }
      else
        D.x[] = 0.;
    }

    diffusion (lambda_thixo, dt, D, theta = theta_lambda);

    foreach()
      lambda_thixo[] = clamp(f[],0.,1.) > THIXO_DIFF_FMIN ?
                       thixo_clamp_lambda(lambda_thixo[]) : thixo_lambda_air;
    boundary ({lambda_thixo});

    thixo_correct_lambda_mass (mass_before_diffusion);
  }

  thixo_sync_tracer_from_lambda();
  thixo_update_gammadot();
}