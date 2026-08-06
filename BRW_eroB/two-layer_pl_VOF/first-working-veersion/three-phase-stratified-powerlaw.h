/**
# Stratified three-phase generalized-Newtonian VOF properties

AMR-oriented branch for the topology

    lower liquid | upper liquid | air.

Two nested cumulative VOF functions are transported:

* `f_lower`: lower-liquid fraction;
* `f`: total-liquid fraction (lower + upper).

The physical fractions are

    f1 = f_lower,
    f2 = f - f_lower,
    f3 = 1 - f.

When the internal and free-surface interfaces remain separated, this
representation gives unit occupancy algebraically and avoids the
three-independent-fraction projection which is unsuitable for repeated
AMR.  A small hierarchy repair (`f_lower <= f`) is retained only as a
monitored safety operation.

The generalized-Newtonian invariant is

    gammaDot = sqrt(2 D:D),
    D = (grad(u) + grad(u)^T)/2,

which reduces to |du/dy| in parallel simple shear.
*/

#ifndef BASILISK_THREE_PHASE_STRATIFIED_POWERLAW_H
#define BASILISK_THREE_PHASE_STRATIFIED_POWERLAW_H

#include <float.h>
#include <math.h>
#include "vof.h"

#if dimension != 2
# error "The supplied stratified roll-wave branch is currently Cartesian 2D only"
#endif
#ifdef AXI
# error "Axisymmetric strain-rate terms are not implemented"
#endif

#ifndef STRATIFIED_EPS
# define STRATIFIED_EPS 1.e-14
#endif
#ifndef STRATIFIED_TOPOLOGY_ABORT_TOL
# define STRATIFIED_TOPOLOGY_ABORT_TOL 1.e-3
#endif
#ifndef STRATIFIED_ABORT_ON_TOPOLOGY
# define STRATIFIED_ABORT_ON_TOPOLOGY 1
#endif

/* For large liquid/ambient viscosity ratios, Basilisk examples recommend the
   harmonic rather than arithmetic interfacial viscosity.  Set this to zero
   to reproduce the original arithmetic face mixture. */
#ifndef STRATIFIED_HARMONIC_FACE_VISCOSITY
# define STRATIFIED_HARMONIC_FACE_VISCOSITY 1
#endif

/** `f` is deliberately the total-liquid VOF so that Basilisk/GfsView
    conventions recognise the free surface. */
scalar f[], f_lower[];
scalar * interfaces = {f_lower, f};

/** Diagnostic physical phase fractions, not independently advected. */
scalar f1[], f2[], f3[];

/** Phase densities. */
double rho1 = 1., rho2 = 1., rho3 = 1.;

typedef struct {
  double K;          // dynamic consistency [Pa s^n]
  double n;          // power-law index
  double tau0;       // yield stress [Pa]
  double m;          // Papanastasiou parameter [s]
  double mu_min;     // lower apparent-viscosity bound [Pa s]
  double mu_max;     // upper apparent-viscosity bound [Pa s]
  double gamma_min;  // shear-rate floor [1/s]
  int regularized;   // nonzero => Papanastasiou regularisation
} StratifiedRheology;

StratifiedRheology rheology1 = {
  .K = 0., .n = 1., .tau0 = 0., .m = 0.,
  .mu_min = 0., .mu_max = DBL_MAX, .gamma_min = 1.e-12,
  .regularized = 0
};
StratifiedRheology rheology2 = {
  .K = 0., .n = 1., .tau0 = 0., .m = 0.,
  .mu_min = 0., .mu_max = DBL_MAX, .gamma_min = 1.e-12,
  .regularized = 0
};
StratifiedRheology rheology3 = {
  .K = 0., .n = 1., .tau0 = 0., .m = 0.,
  .mu_min = 0., .mu_max = DBL_MAX, .gamma_min = 1.e-12,
  .regularized = 0
};

face vector alphav[];
scalar rhov[], gammaDot[], muMix[];

/** Most recent hierarchy/boundedness diagnostics. */
double stratified_hierarchy_max = 0.;
double stratified_bound_max = 0.;
double stratified_volume_correction_lower = 0.;
double stratified_volume_correction_liquid = 0.;
long stratified_corrected_cells = 0;

static inline double st_clip01 (double a)
{
  return clamp (a, 0., 1.);
}

/** Convert cumulative fractions to physical fractions without changing
    the stored VOF fields. */
static inline void stratified_phase_values (double lower, double liquid,
                                             double * a1, double * a2,
                                             double * a3)
{
  liquid = st_clip01 (liquid);
  lower = st_clip01 (lower);
  if (lower > liquid)
    lower = liquid;
  *a1 = lower;
  *a2 = liquid - lower;
  *a3 = 1. - liquid;
}

static inline double stratified_density_values (double lower,
                                                  double liquid)
{
  double a1, a2, a3;
  stratified_phase_values (lower, liquid, &a1, &a2, &a3);
  return rho1*a1 + rho2*a2 + rho3*a3;
}

static inline double stratified_density_cell (double lower,
                                                double liquid)
{
  return stratified_density_values (lower, liquid);
}

static inline double stratified_apparent_viscosity
  (const StratifiedRheology * r, double gd)
{
  const double g = max (fabs(gd), r->gamma_min);
  double muv = r->K > 0. ? r->K*pow(g, r->n - 1.) : 0.;

  if (r->tau0 > 0.) {
    if (r->regularized) {
      const double x = r->m*g;
      double ratio;
      if (fabs(x) < 1.e-5)
        ratio = r->m*(1. - x/2. + x*x/6. - x*x*x/24.);
      else
        ratio = -expm1(-x)/g;
      muv += r->tau0*ratio;
    }
    else
      muv += r->tau0/g;
  }
  return clamp (muv, r->mu_min, r->mu_max);
}

static void stratified_update_phase_fields (void)
{
  foreach() {
    double a1, a2, a3;
    stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
    f1[] = a1;
    f2[] = a2;
    f3[] = a3;
  }
  boundary ({f1, f2, f3});
}

/**
Enforce 0 <= f_lower <= f <= 1 while preserving total-liquid `f` when
there is a hierarchy violation.  For the intended separated-interface
topology, corrections should remain at roundoff.  A large correction is
reported as a topology failure rather than silently accepted.
*/
static void stratified_enforce_hierarchy (void)
{
  double max_h = 0., max_b = 0.;
  double before_lower = 0., after_lower = 0.;
  double before_liquid = 0., after_liquid = 0.;
  long corrected = 0;

  foreach (reduction(max:max_h) reduction(max:max_b)
           reduction(+:before_lower) reduction(+:after_lower)
           reduction(+:before_liquid) reduction(+:after_liquid)
           reduction(+:corrected)) {
    const double vol = dv();
    const double old_lower = f_lower[];
    const double old_liquid = f[];
    before_lower += old_lower*vol;
    before_liquid += old_liquid*vol;

    const double bound_error = max(max(-old_lower, old_lower - 1.),
                                   max(-old_liquid, old_liquid - 1.));
    const double hierarchy_error = old_lower - old_liquid;
    max_b = max(max_b, max(0., bound_error));
    max_h = max(max_h, max(0., hierarchy_error));

    double liquid = st_clip01(old_liquid);
    double lower = st_clip01(old_lower);
    if (lower > liquid)
      lower = liquid; // preserve the free-surface/total-liquid volume

    if (fabs(lower - old_lower) > 10.*DBL_EPSILON ||
        fabs(liquid - old_liquid) > 10.*DBL_EPSILON)
      corrected++;

    f_lower[] = lower;
    f[] = liquid;
    after_lower += lower*vol;
    after_liquid += liquid*vol;
  }

  boundary ({f_lower, f});
  stratified_hierarchy_max = max_h;
  stratified_bound_max = max_b;
  stratified_volume_correction_lower = after_lower - before_lower;
  stratified_volume_correction_liquid = after_liquid - before_liquid;
  stratified_corrected_cells = corrected;

#if STRATIFIED_ABORT_ON_TOPOLOGY
  if (max(max_h, max_b) > STRATIFIED_TOPOLOGY_ABORT_TOL) {
    fprintf (ferr,
             "FATAL: nested VOF topology violation: hierarchy=%g bounds=%g "
             "(tolerance=%g, t=%g, i=%d)\n",
             max_h, max_b, (double) STRATIFIED_TOPOLOGY_ABORT_TOL, t, iter);
    exit (1);
  }
#endif

  stratified_update_phase_fields();
}

static void stratified_validate_rheology (const char * name,
                                           const StratifiedRheology * r)
{
  if (r->K < 0. || r->n <= 0. || r->tau0 < 0. ||
      r->gamma_min <= 0. || r->mu_min < 0. ||
      r->mu_max < r->mu_min || (r->regularized && r->m <= 0.)) {
    fprintf (ferr,
             "Invalid rheology for %s: K=%g n=%g tau0=%g m=%g "
             "mu_min=%g mu_max=%g gamma_min=%g regularized=%d\n",
             name, r->K, r->n, r->tau0, r->m, r->mu_min,
             r->mu_max, r->gamma_min, r->regularized);
    exit (1);
  }
}

event defaults (i = 0)
{
  if (rho1 <= 0. || rho2 <= 0. || rho3 <= 0.) {
    fprintf (ferr, "All three phase densities must be positive.\n");
    exit (1);
  }
  stratified_validate_rheology ("lower liquid", &rheology1);
  stratified_validate_rheology ("upper liquid", &rheology2);
  stratified_validate_rheology ("air", &rheology3);

  alpha = alphav;
  rho = rhov;
  mu = new face vector;
  reset ((scalar *){mu}, 0.);

#if TREE
  f1.prolongation = refine_bilinear;
  f2.prolongation = refine_bilinear;
  f3.prolongation = refine_bilinear;
#endif

  display ("draw_vof (c = 'f'); draw_vof (c = 'f_lower');");
}

/** Lagged apparent-viscosity update. */
event properties (i++)
{
  boundary ((scalar *){u});

  foreach() {
    const double dudx = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    const double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    const double dvdx = (u.y[1,0] - u.y[-1,0])/(2.*Delta);
    const double dvdy = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    const double gd2 = 2.*(sq(dudx) + sq(dvdy)) + sq(dudy + dvdx);
    gammaDot[] = sqrt(max(gd2, 0.));

    double a1, a2, a3;
    stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
    const double m1 = stratified_apparent_viscosity (&rheology1,
                                                       gammaDot[]);
    const double m2 = stratified_apparent_viscosity (&rheology2,
                                                       gammaDot[]);
    const double m3 = stratified_apparent_viscosity (&rheology3,
                                                       gammaDot[]);
    muMix[] = a1*m1 + a2*m2 + a3*m3;
    rhov[] = cm[]*(a1*rho1 + a2*rho2 + a3*rho3);
    f1[] = a1;
    f2[] = a2;
    f3[] = a3;
  }
  boundary ({gammaDot, muMix, rhov, f1, f2, f3});

  foreach_face() {
    const double lower = (f_lower[] + f_lower[-1])/2.;
    const double liquid = (f[] + f[-1])/2.;
    double a1, a2, a3;
    stratified_phase_values (lower, liquid, &a1, &a2, &a3);

    const double gd = max(0., (gammaDot[] + gammaDot[-1])/2.);
    const double m1 = stratified_apparent_viscosity (&rheology1, gd);
    const double m2 = stratified_apparent_viscosity (&rheology2, gd);
    const double m3 = stratified_apparent_viscosity (&rheology3, gd);
    const double rhof = a1*rho1 + a2*rho2 + a3*rho3;

    alphav.x[] = fm.x[]/rhof;
    face vector muv = mu;
#if STRATIFIED_HARMONIC_FACE_VISCOSITY
    const double inv_mu = a1/max(m1, DBL_MIN) +
      a2/max(m2, DBL_MIN) + a3/max(m3, DBL_MIN);
    muv.x[] = fm.x[]/(inv_mu + DBL_MIN);
#else
    muv.x[] = fm.x[]*(a1*m1 + a2*m2 + a3*m3);
#endif
  }
}

#endif
