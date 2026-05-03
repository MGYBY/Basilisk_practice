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
#define THIXO_MU_MAX (2.5*(1.0/0.05)) // roughly estimated from Balmforth as 1/\epsilon
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
/* Face-centred Frobenius norm D2 = sqrt(D:D), computed in the same
   staggered style as the power-law/viscoplastic face-viscosity method. */
face vector D2f_thixo[];
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

  The structural parameter is a liquid material variable. The shear rate used
  in the rebuilding/destruction source term should therefore approximate the
  liquid-side strain rate rather than a numerical air--liquid velocity jump.

  The helper below returns a physical derivative dq/dx. If both neighbouring
  cells are sufficiently liquid-filled, it uses the standard centred gradient.
  If only one neighbour is liquid-filled, it falls back to a one-sided liquid
  gradient. If the local stencil is not supported by liquid cells, it returns
  zero.  This keeps the free-surface treatment conservative but avoids forcing
  an artificial gas-phase lambda dynamics.
*/
static inline bool thixo_source_cell (Point point)
{
  return clamp(f[],0.,1.) > THIXO_GAMMADOT_FMIN;
}

static inline bool thixo_liquid_enough (double ff, double threshold)
{
  return clamp(ff,0.,1.) > threshold;
}

static inline double thixo_limited_derivative (double qc, double qp, double qm,
                                              double fc, double fp, double fm,
                                              double Delta)
{
  bool has_c = thixo_liquid_enough (fc, THIXO_GAMMADOT_FMIN);
  bool has_p = thixo_liquid_enough (fp, THIXO_GAMMADOT_NEIGHBOR_FMIN);
  bool has_m = thixo_liquid_enough (fm, THIXO_GAMMADOT_NEIGHBOR_FMIN);

  if (!has_c)
    return 0.;
  if (has_p && has_m)
    return (qp - qm)/(2.*Delta);
  if (has_p)
    return (qp - qc)/Delta;
  if (has_m)
    return (qc - qm)/Delta;
  return 0.;
}

static inline double thixo_pair_derivative (double qR, double qL,
                                            double fR, double fL,
                                            double Delta)
{
  if (thixo_liquid_enough (fR, THIXO_GAMMADOT_FMIN) &&
      thixo_liquid_enough (fL, THIXO_GAMMADOT_FMIN))
    return (qR - qL)/Delta;
  return 0.;
}

/*
  Power-law-style face-centred D2 calculation.

  This is deliberately different from the older purely cell-centred formula.
  The viscous solver stores mu on faces, and the tested power-law/viscoplastic
  Basilisk formulation evaluates the strain-rate invariant on faces first,
  then forms a cell-centred diagnostic by averaging adjacent face values.

  In 2-D Cartesian form, at an x-face,

      D11 = du_x/dx,
      D22 = averaged du_y/dy from the two neighbouring cells,
      D12 = 0.5*(du_y/dx + averaged du_x/dy),
      D2  = sqrt(D11^2 + D22^2 + 2 D12^2).

  A symmetric construction is used at y-faces.  The thixotropic source uses
  gammadot = sqrt(2)*D2_cell, where D2_cell is the average of the surrounding
  face-centred D2 values.  The limited derivatives prevent pure gas cells from
  entering the liquid structural breakdown rate near the VOF interface.
*/
static inline void thixo_update_gammadot()
{
  face vector D2f = D2f_thixo;

#if dimension == 1
  foreach_face(x) {
    double D11 = thixo_pair_derivative (u.x[], u.x[-1], f[], f[-1], Delta);
    D2f.x[] = sqrt (sq(D11));
  }
  foreach()
    gammadot_thixo[] = thixo_source_cell(point) ?
      sqrt(2.)*0.5*(D2f.x[] + D2f.x[1]) : 0.;

#elif dimension == 2
  foreach_face(x) {
    double D11 = thixo_pair_derivative (u.x[], u.x[-1,0],
                                        f[], f[-1,0], Delta);

    double uy_y_R = thixo_limited_derivative (u.y[], u.y[0,1], u.y[0,-1],
                                              f[], f[0,1], f[0,-1], Delta);
    double uy_y_L = thixo_limited_derivative (u.y[-1,0], u.y[-1,1], u.y[-1,-1],
                                              f[-1,0], f[-1,1], f[-1,-1], Delta);
    double D22 = 0.5*(uy_y_R + uy_y_L);

    double uy_x = thixo_pair_derivative (u.y[], u.y[-1,0],
                                         f[], f[-1,0], Delta);
    double ux_y_R = thixo_limited_derivative (u.x[], u.x[0,1], u.x[0,-1],
                                              f[], f[0,1], f[0,-1], Delta);
    double ux_y_L = thixo_limited_derivative (u.x[-1,0], u.x[-1,1], u.x[-1,-1],
                                              f[-1,0], f[-1,1], f[-1,-1], Delta);
    double D12 = 0.5*(uy_x + 0.5*(ux_y_R + ux_y_L));

    D2f.x[] = sqrt (max(0., sq(D11) + sq(D22) + 2.*sq(D12)));
  }

  foreach_face(y) {
    double ux_x_U = thixo_limited_derivative (u.x[], u.x[1,0], u.x[-1,0],
                                              f[], f[1,0], f[-1,0], Delta);
    double ux_x_D = thixo_limited_derivative (u.x[0,-1], u.x[1,-1], u.x[-1,-1],
                                              f[0,-1], f[1,-1], f[-1,-1], Delta);
    double D11 = 0.5*(ux_x_U + ux_x_D);

    double D22 = thixo_pair_derivative (u.y[], u.y[0,-1],
                                        f[], f[0,-1], Delta);

    double ux_y = thixo_pair_derivative (u.x[], u.x[0,-1],
                                         f[], f[0,-1], Delta);
    double uy_x_U = thixo_limited_derivative (u.y[], u.y[1,0], u.y[-1,0],
                                              f[], f[1,0], f[-1,0], Delta);
    double uy_x_D = thixo_limited_derivative (u.y[0,-1], u.y[1,-1], u.y[-1,-1],
                                              f[0,-1], f[1,-1], f[-1,-1], Delta);
    double D12 = 0.5*(ux_y + 0.5*(uy_x_U + uy_x_D));

    D2f.y[] = sqrt (max(0., sq(D11) + sq(D22) + 2.*sq(D12)));
  }

  boundary ((scalar *){D2f});

  foreach() {
    if (thixo_source_cell (point)) {
      double D2c = 0.25*(D2f.x[] + D2f.x[1,0] +
                         D2f.y[] + D2f.y[0,1]);
      gammadot_thixo[] = sqrt(2.)*max(0., D2c);
    }
    else
      gammadot_thixo[] = 0.;
  }

#else // dimension == 3
  /* Fallback: the previous cell-centred invariant is kept for 3-D until a
     fully face-centred 3-D version is validated. */
  foreach() {
    if (thixo_source_cell (point)) {
      double ux_x = thixo_limited_derivative (u.x[], u.x[1,0,0], u.x[-1,0,0],
                                              f[], f[1,0,0], f[-1,0,0], Delta);
      double uy_y = thixo_limited_derivative (u.y[], u.y[0,1,0], u.y[0,-1,0],
                                              f[], f[0,1,0], f[0,-1,0], Delta);
      double uz_z = thixo_limited_derivative (u.z[], u.z[0,0,1], u.z[0,0,-1],
                                              f[], f[0,0,1], f[0,0,-1], Delta);
      double ux_y = thixo_limited_derivative (u.x[], u.x[0,1,0], u.x[0,-1,0],
                                              f[], f[0,1,0], f[0,-1,0], Delta);
      double uy_x = thixo_limited_derivative (u.y[], u.y[1,0,0], u.y[-1,0,0],
                                              f[], f[1,0,0], f[-1,0,0], Delta);
      double ux_z = thixo_limited_derivative (u.x[], u.x[0,0,1], u.x[0,0,-1],
                                              f[], f[0,0,1], f[0,0,-1], Delta);
      double uz_x = thixo_limited_derivative (u.z[], u.z[1,0,0], u.z[-1,0,0],
                                              f[], f[1,0,0], f[-1,0,0], Delta);
      double uy_z = thixo_limited_derivative (u.y[], u.y[0,0,1], u.y[0,0,-1],
                                              f[], f[0,0,1], f[0,0,-1], Delta);
      double uz_y = thixo_limited_derivative (u.z[], u.z[0,1,0], u.z[0,-1,0],
                                              f[], f[0,1,0], f[0,-1,0], Delta);
      double dxy = 0.5*(ux_y + uy_x);
      double dxz = 0.5*(ux_z + uz_x);
      double dyz = 0.5*(uy_z + uz_y);
      gammadot_thixo[] = sqrt(max(0., 2.*(sq(ux_x) + sq(uy_y) + sq(uz_z) +
                                             2.*sq(dxy) + 2.*sq(dxz) + 2.*sq(dyz))));
    }
    else
      gammadot_thixo[] = 0.;
  }
#endif

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

  foreach_dimension()
    D2f_thixo.x.refine = refine_face;

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
      /*
        Liquid-side face reconstruction of lambda from the conservative
        VOF tracer flambda = f*lambda.  This is written directly inside
        foreach_face() rather than as a separate Point function so that qcc
        can see the stencil accesses without any ambiguity.
      */
      double fR_lam = clamp (f[],   0., 1.);
      double fL_lam = clamp (f[-1], 0., 1.);
      double denom_lam = fR_lam + fL_lam;
      double lf = thixo_lambda_air;

      if (denom_lam > THIXO_FERR) {
        double qR_lam = fR_lam > THIXO_FERR ? flambda[]   : 0.;
        double qL_lam = fL_lam > THIXO_FERR ? flambda[-1] : 0.;
        lf = thixo_clamp_lambda ((qR_lam + qL_lam)/denom_lam);
      }

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