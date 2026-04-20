/**
# Two-phase VOF solver for thixotropic liquid--air flows

This header extends Basilisk's standard two-phase VOF framework so that
fluid 1 is no longer a constant-viscosity liquid. Instead, fluid 1 is a
thixotropic liquid whose viscosity depends on a transported structure
parameter `lambda_thixo`.

The implementation follows the coding strategy described in the uploaded
notes:

1. keep the usual VOF volume fraction `f` for the liquid,
2. transport a conservative structural tracer `flambda = f*lambda`,
3. reconstruct the cell-centred structure `lambda_thixo = flambda/f` in
   cells containing liquid,
4. compute the liquid viscosity from the rheology,
5. form the liquid/air mixture viscosity with the harmonic average, and
6. advance the source and diffusion terms of the thixotropic equation in
   an operator-split fashion.

The dimensionless thixotropic model implemented here is

```text
∂t λ + ∇·(u λ) = (So/T) [ (1 - λ) - Γ λ γdot + κ ∇²λ ]
μ_liquid(λ)     = 1 / [ (1 - λ) ((1 - λ)(1 - a) + a) ]
```

The VOF advection of `flambda` is delegated to Basilisk's conservative
VOF machinery. The reaction/source part is advanced with a two-stage
SSPRK2 update, and the diffusive part is advanced with `diffusion.h`.
*/

#include "vof.h"
#include "diffusion.h"
#include "fractions.h"

/* --------------------------------------------------------------------------
   Small numerical tolerances used repeatedly in the structural update.
   -------------------------------------------------------------------------- */
#ifndef THIXO_FERR
#define THIXO_FERR 1e-12
#endif
#ifndef THIXO_LAMBDA_MIN
#define THIXO_LAMBDA_MIN 0.
#endif
#ifndef THIXO_LAMBDA_MAX
#define THIXO_LAMBDA_MAX 1.
#endif
#ifndef THIXO_EPS
#define THIXO_EPS 1e-30
#endif

/* --------------------------------------------------------------------------
   Standard two-phase VOF field.

   Convention used throughout this solver:
     f = 1 : thixotropic liquid phase,
     f = 0 : air phase.
   -------------------------------------------------------------------------- */
scalar f[], * interfaces = {f};

/* Conservative structural tracer transported with the VOF field.
   In liquid cells, flambda = f*lambda_thixo. */
scalar flambda[];

/* Cell-centred diagnostic/auxiliary fields.
   - lambda_thixo    : reconstructed structure parameter λ
   - gammadot_thixo  : scalar strain-rate magnitude entering the source term
   - mu_liquid_thixo : liquid-phase viscosity implied by λ (before mixing with air)
*/
scalar lambda_thixo[];
scalar gammadot_thixo[];
scalar mu_liquid_thixo[];

/* --------------------------------------------------------------------------
   Standard two-phase properties.

   `rho1` and `rho2` are still the liquid and air densities used by the
   dimensionless two-phase solver.

   `mu1` is kept for consistency with Basilisk's default interface, but the
   actual liquid viscosity used in the properties event is `mu_liquid_thixo`,
   i.e. a spatially varying viscosity reconstructed from λ.
   -------------------------------------------------------------------------- */
double rho1 = 1., mu1 = 1.;
double rho2 = 1., mu2 = 0.;

/* Dimensionless thixotropic parameters. */
double thixo_T = 1.;
double thixo_Gamma = 0.;
double thixo_kappa = 0.;
double thixo_a = 0.2;
double thixo_So = 0.05;

double thixo_lambda_air = 0.;

/* Standard Basilisk auxiliary fields for variable-density flows. */
face vector alphav[];
scalar rhov[];

/* --------------------------------------------------------------------------
   Helper functions.
   -------------------------------------------------------------------------- */

/* Keep λ inside the physically meaningful interval. */
static inline double thixo_clamp_lambda (double l)
{
  return max (THIXO_LAMBDA_MIN, min (THIXO_LAMBDA_MAX, l));
}

/* Equation (12) in the uploaded draft uses (So/T) as the common prefactor. */
static inline double thixo_source_prefactor()
{
  return thixo_So/(thixo_T + THIXO_EPS);
}

/* Liquid-phase viscosity from the thixotropic rheology. */
static inline double thixo_mu_fluid_from_lambda (double l)
{
  double om = 1. - thixo_clamp_lambda (l);
  double denom = om*(om*(1. - thixo_a) + thixo_a);
  return 1./max (denom, THIXO_EPS);
}

/* Harmonic average recommended for large viscosity ratios. */
static inline double thixo_mu_mix_harmonic (double mu_liq, double mu_air,
                                            double ff)
{
  double fcl = clamp (ff, 0., 1.);
  return 1./(fcl/max(mu_liq, THIXO_EPS) +
             (1. - fcl)/max(mu_air, THIXO_EPS));
}

/* --------------------------------------------------------------------------
   Synchronisation between the conservative tracer `flambda` and the explicit
   cell-centred structural field `lambda_thixo`.

   The conservative field is what must be advected by the VOF machinery.
   The explicit λ field is what we want for diagnostics, source terms and the
   viscosity law.
   -------------------------------------------------------------------------- */

/* Recover λ from the conservative tracer after VOF advection. */
static inline void thixo_sync_lambda_from_tracer()
{
  foreach() {
    double fc = clamp (f[], 0., 1.);
    lambda_thixo[] = (fc > THIXO_FERR ?
                      thixo_clamp_lambda (flambda[]/(fc + THIXO_EPS)) :
                      thixo_lambda_air);
    mu_liquid_thixo[] = thixo_mu_fluid_from_lambda (lambda_thixo[]);
  }
  boundary ({lambda_thixo, mu_liquid_thixo});
}

/* Rebuild the conservative tracer after source/diffusion updates on λ. */
static inline void thixo_sync_tracer_from_lambda()
{
  foreach() {
    double fc = clamp (f[], 0., 1.);
    lambda_thixo[] = (fc > THIXO_FERR ?
                      thixo_clamp_lambda (lambda_thixo[]) :
                      thixo_lambda_air);
    flambda[] = fc*lambda_thixo[];
    mu_liquid_thixo[] = thixo_mu_fluid_from_lambda (lambda_thixo[]);
  }
  boundary ({flambda, lambda_thixo, mu_liquid_thixo});
}

/* --------------------------------------------------------------------------
   Strain-rate magnitude.

   The uploaded notes asked to follow the strain-rate norm used in the user's
   tested `two-phasePL.h` implementation. That implementation evaluates the
   tensor norm

       sqrt(D_ij D_ij)

   and then uses γdot = sqrt(2) * sqrt(D_ij D_ij).

   Here we use the same invariant, but evaluated at cell centres because the
   structural source term is cell centred.
   -------------------------------------------------------------------------- */
static inline double thixo_D2_point (Point point)
{
#if dimension == 1
  double dxx = (u.x[1] - u.x[-1])/(2.*Delta);
  return fabs (dxx);
#elif dimension == 2
  double dxx = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
  double dyy = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
  double dxy = 0.5*((u.x[0,1] - u.x[0,-1])/(2.*Delta) +
                    (u.y[1,0] - u.y[-1,0])/(2.*Delta));
  return sqrt (sq(dxx) + sq(dyy) + 2.*sq(dxy));
#else
  double dxx = (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta);
  double dyy = (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta);
  double dzz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
  double dxy = 0.5*((u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta) +
                    (u.y[1,0,0] - u.y[-1,0,0])/(2.*Delta));
  double dxz = 0.5*((u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta) +
                    (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta));
  double dyz = 0.5*((u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta) +
                    (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta));
  return sqrt (sq(dxx) + sq(dyy) + sq(dzz) +
               2.*sq(dxy) + 2.*sq(dxz) + 2.*sq(dyz));
#endif
}

static inline void thixo_update_gammadot()
{
  foreach() {
    double D2 = thixo_D2_point (point);
    gammadot_thixo[] = sqrt (2.)*D2;
  }
  boundary ({gammadot_thixo});
}

/* --------------------------------------------------------------------------
   Basilisk events.
   -------------------------------------------------------------------------- */

/* Allocate the variable-density fields and register `flambda` as a VOF tracer. */
event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;

  if (mu1 || mu2)
    mu = new face vector;

  /* `flambda` must be transported consistently with `f`. */
  f.tracers = list_append (f.tracers, flambda);
  flambda.inverse = false;

#if TREE
  /*
     `flambda` behaves like a VOF concentration: on refinement it must remain
     consistent with the parent `f` field, which is why `vof_concentration_refine`
     is used here.
  */
  flambda.restriction = restriction_volume_average;
  flambda.refine = flambda.prolongation = vof_concentration_refine;
  flambda.dirty = true;
  flambda.c = f;

  /* The diagnostic fields are reconstructed linearly on trees. */
  lambda_thixo.refine = lambda_thixo.prolongation = refine_linear;
  lambda_thixo.restriction = restriction_volume_average;
  lambda_thixo.dirty = true;

  gammadot_thixo.refine = gammadot_thixo.prolongation = refine_linear;
  gammadot_thixo.restriction = restriction_volume_average;
  gammadot_thixo.dirty = true;

  mu_liquid_thixo.refine = mu_liquid_thixo.prolongation = refine_linear;
  mu_liquid_thixo.restriction = restriction_volume_average;
  mu_liquid_thixo.dirty = true;
#endif

  display ("draw_vof (c = 'f');");
}

#ifndef rho
# define rho(ff) (clamp((ff),0.,1.)*(rho1 - rho2) + rho2)
#endif

/* Optional smearing of material properties, following Basilisk's two-phase solver. */
#if FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
#ifndef sf
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
  sf.dirty = true;
#endif
}

/* Material properties needed by centered.h: density and face viscosity. */
event properties (i++)
{
  /* Reconstruct λ after the conservative VOF transport. */
  thixo_sync_lambda_from_tracer();

  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);

    if (mu1 || mu2) {
      face vector muv = mu;

      /*
         The liquid viscosity is evaluated from the average λ at the face and then
         mixed harmonically with the air viscosity.
      */
      double lf = thixo_clamp_lambda ((lambda_thixo[] + lambda_thixo[-1])/2.);
      double mu_liq = thixo_mu_fluid_from_lambda (lf);
      muv.x[] = fm.x[]*thixo_mu_mix_harmonic (mu_liq, mu2, ff);
    }
  }

  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true;
#endif
}

/* --------------------------------------------------------------------------
   Operator-split update for the thixotropic structure field.

   Step 1: recover λ from flambda after the conservative VOF step.
   Step 2: advance the source term with SSPRK2.
   Step 3: advance diffusion implicitly with diffusion.h.
   Step 4: rebuild the conservative tracer flambda = f*lambda.
   -------------------------------------------------------------------------- */
event tracer_diffusion (i++, last)
{
  thixo_sync_lambda_from_tracer();
  thixo_update_gammadot();

  scalar lambda0[], lambda1[];
  double pref = thixo_source_prefactor();

  /* SSPRK2 stage 1. */
  foreach() {
    lambda0[] = lambda_thixo[];
    if (clamp(f[], 0., 1.) > THIXO_FERR) {
      double rhs0 = pref*((1. - lambda0[]) -
                          thixo_Gamma*lambda0[]*gammadot_thixo[]);
      lambda1[] = thixo_clamp_lambda (lambda0[] + dt*rhs0);
    }
    else
      lambda1[] = thixo_lambda_air;
  }
  boundary ({lambda0, lambda1});

  /* SSPRK2 stage 2. */
  foreach() {
    if (clamp(f[], 0., 1.) > THIXO_FERR) {
      double rhs1 = pref*((1. - lambda1[]) -
                          thixo_Gamma*lambda1[]*gammadot_thixo[]);
      lambda_thixo[] = thixo_clamp_lambda (0.5*(lambda0[] + lambda1[] + dt*rhs1));
    }
    else
      lambda_thixo[] = thixo_lambda_air;
  }
  boundary ({lambda_thixo});

  /* Diffusion only acts where liquid is present. */
  if (thixo_kappa > 0.) {
    face vector D[];
    foreach_face() {
      double ff = clamp ((f[] + f[-1])/2., 0., 1.);
      D.x[] = fm.x[]*(pref*thixo_kappa)*ff;
    }
    diffusion (lambda_thixo, dt, D);

    /* Clamp again after the implicit solve, which may overshoot slightly. */
    foreach()
      if (clamp(f[], 0., 1.) > THIXO_FERR)
        lambda_thixo[] = thixo_clamp_lambda (lambda_thixo[]);
      else
        lambda_thixo[] = thixo_lambda_air;
    boundary ({lambda_thixo});
  }

  thixo_sync_tracer_from_lambda();
  thixo_update_gammadot();
}
