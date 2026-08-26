/**
# Dimensionless stratified three-phase generalized-Newtonian properties

This header implements the one-fluid material fields for the maintained
dimensionless problem documented in `README.md` and the formulation-report
addendum. The reference scales are

    x~ = H_l x,          u~ = U_l u,
    t~ = H_l/U_l t,      p~ - p_a = rho_l U_l^2 p,
    eta~ = rho_l U_l H_l eta.

Consequently, `rho1`, `rho2`, and `rho3` are density ratios, while the face
field `mu` stores dimensionless *dynamic* apparent viscosity.  Basilisk's
`centered.h` therefore advances

    r Du/Dt = -grad(p) + div(2 eta D) + r a.

The stratified topology is represented by two nested cumulative VOF fields:

    f_lower = alpha_lower,
    f       = alpha_lower + alpha_upper.

The physical fractions are reconstructed algebraically as

    alpha_1 = f_lower,
    alpha_2 = f - f_lower,
    alpha_3 = 1 - f.

Let g = |gamma|, g_safe = max(g, gamma_floor), and

    clip(q) = min(max(q, eta_min), eta_max).

Four rheological options are available through the integer `model` member.
The final clipping is applied after all consistency and yield-stress terms:

    0  smooth power law
       eta = clip[Lambda (g^2 + epsilon^2)^((n-1)/2)]

    1  bounded ideal power law
       eta = clip[Lambda g_safe^(n-1)]

    2  Papanastasiou Herschel--Bulkley
       eta = clip[Lambda g_safe^(n-1) + tau0 P(g)]

       P(g) = [1 - exp(-m g)]/g for g > 0 and P(0) = m. The code uses
       a small-argument series and the configured floor for stable evaluation.

    3  bounded Herschel--Bulkley
       eta = clip[Lambda g_safe^(n-1) + tau0/g_safe].

All variables in these expressions are dimensionless. The strain-rate
invariant is gamma = sqrt(2 D:D), which reduces to |du/dz| in parallel shear.
The case driver assigns the air phase a constant Newtonian viscosity by setting
n=1 and Lambda=eta_min=eta_max=the generated CASE_AIR_ETA.
*/

#ifndef BASILISK_THREE_PHASE_STRATIFIED_POWERLAW_DIMENSIONLESS_H
#define BASILISK_THREE_PHASE_STRATIFIED_POWERLAW_DIMENSIONLESS_H

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
#ifndef STRATIFIED_HARMONIC_FACE_VISCOSITY
# define STRATIFIED_HARMONIC_FACE_VISCOSITY 1
#endif

#define STRATIFIED_SMOOTH_POWER_LAW 0
#define STRATIFIED_BOUNDED_POWER_LAW 1
#define STRATIFIED_PAPANASTASIOU_HB 2
#define STRATIFIED_BOUNDED_HB 3

/** `f` is total liquid, which also makes the free surface visible to GfsView. */
scalar f[], f_lower[];
scalar * interfaces = {f_lower, f};

/** Reconstructed physical fractions; these are diagnostic, not advected. */
scalar f1[], f2[], f3[];

/** Dimensionless densities r_i = rho_i/rho_l. */
double rho1 = 1., rho2 = 1., rho3 = 1.;

typedef struct {
  int model;
  double Lambda;       // dimensionless consistency coefficient
  double n;            // power-law index
  double epsilon;      // dimensionless low-shear transition rate
  double tau0;         // dimensionless yield stress, tau0~/(rho_l U_l^2)
  double m;            // dimensionless Papanastasiou parameter
  double eta_min;      // dimensionless lower viscosity bound
  double eta_max;      // dimensionless upper viscosity bound
  double gamma_floor;  // positive dimensionless shear-rate floor
} StratifiedRheology;

StratifiedRheology rheology1 = {
  .model = STRATIFIED_BOUNDED_POWER_LAW,
  .Lambda = 0., .n = 1., .epsilon = 0., .tau0 = 0., .m = 0.,
  .eta_min = 0., .eta_max = DBL_MAX, .gamma_floor = 1.e-12
};
StratifiedRheology rheology2 = {
  .model = STRATIFIED_BOUNDED_POWER_LAW,
  .Lambda = 0., .n = 1., .epsilon = 0., .tau0 = 0., .m = 0.,
  .eta_min = 0., .eta_max = DBL_MAX, .gamma_floor = 1.e-12
};
StratifiedRheology rheology3 = {
  .model = STRATIFIED_BOUNDED_POWER_LAW,
  .Lambda = 0., .n = 1., .epsilon = 0., .tau0 = 0., .m = 0.,
  .eta_min = 0., .eta_max = DBL_MAX, .gamma_floor = 1.e-12
};

face vector alphav[];
scalar rhov[], gammaDot[], muMix[];

double stratified_hierarchy_max = 0.;
double stratified_bound_max = 0.;
double stratified_volume_correction_lower = 0.;
double stratified_volume_correction_liquid = 0.;
long stratified_corrected_cells = 0;

static inline double st_clip01 (double value)
{
  return clamp (value, 0., 1.);
}

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

static inline double stratified_density_values (double lower, double liquid)
{
  double a1, a2, a3;
  stratified_phase_values (lower, liquid, &a1, &a2, &a3);
  return rho1*a1 + rho2*a2 + rho3*a3;
}

static inline double stratified_density_cell (double lower, double liquid)
{
  return stratified_density_values (lower, liquid);
}

/** Numerically stable evaluation of [1-exp(-m g)]/g. */
static inline double stratified_papanastasiou_ratio (double g, double m,
                                                      double floor)
{
  if (m <= 0.)
    return 0.;
  const double x = m*g;
  if (fabs(x) < 1.e-5)
    return m*(1. - x/2. + x*x/6. - x*x*x/24.);
  return -expm1(-x)/max(g, floor);
}

static inline double stratified_apparent_viscosity
  (const StratifiedRheology * r, double gd)
{
  const double gabs = fabs(gd);
  const double gsafe = max(gabs, r->gamma_floor);
  double eta;

  if (r->model == STRATIFIED_SMOOTH_POWER_LAW)
    eta = r->Lambda*pow(sq(gabs) + sq(r->epsilon),
                        0.5*(r->n - 1.));
  else
    eta = r->Lambda*pow(gsafe, r->n - 1.);

  if (r->model == STRATIFIED_PAPANASTASIOU_HB && r->tau0 > 0.)
    eta += r->tau0*stratified_papanastasiou_ratio
      (gabs, r->m, r->gamma_floor);
  else if (r->model == STRATIFIED_BOUNDED_HB && r->tau0 > 0.)
    eta += r->tau0/gsafe;

  return clamp (eta, r->eta_min, r->eta_max);
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
Enforce 0 <= f_lower <= f <= 1 while preserving the total-liquid field `f`
when a hierarchy violation occurs.  Corrections should remain at roundoff for
the intended separated-interface topology.
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
      lower = liquid;

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
  const bool valid_model =
    r->model >= STRATIFIED_SMOOTH_POWER_LAW &&
    r->model <= STRATIFIED_BOUNDED_HB;
  if (!valid_model || r->Lambda < 0. || r->n <= 0. ||
      r->epsilon < 0. || r->tau0 < 0. || r->m < 0. ||
      r->gamma_floor <= 0. || r->eta_min < 0. ||
      r->eta_max <= 0. || r->eta_max < r->eta_min ||
      (r->model == STRATIFIED_SMOOTH_POWER_LAW &&
       r->n < 1. && r->epsilon <= 0.) ||
      (r->model == STRATIFIED_PAPANASTASIOU_HB &&
       r->tau0 > 0. && r->m <= 0.)) {
    fprintf (ferr,
             "Invalid dimensionless rheology for %s: model=%d Lambda=%g "
             "n=%g epsilon=%g tau0=%g m=%g eta_min=%g eta_max=%g "
             "gamma_floor=%g\n",
             name, r->model, r->Lambda, r->n, r->epsilon, r->tau0,
             r->m, r->eta_min, r->eta_max, r->gamma_floor);
    exit (1);
  }
}

event defaults (i = 0)
{
  if (rho1 <= 0. || rho2 <= 0. || rho3 <= 0.) {
    fprintf (ferr, "All three dimensionless phase densities must be positive.\n");
    exit (1);
  }
  stratified_validate_rheology ("lower liquid", &rheology1);
  stratified_validate_rheology ("upper liquid", &rheology2);
  stratified_validate_rheology ("ambient", &rheology3);

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

/** Lagged generalized-Newtonian update in dimensionless variables. */
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
    const double eta1 = stratified_apparent_viscosity (&rheology1,
                                                        gammaDot[]);
    const double eta2 = stratified_apparent_viscosity (&rheology2,
                                                        gammaDot[]);
    const double eta3 = stratified_apparent_viscosity (&rheology3,
                                                        gammaDot[]);
    muMix[] = a1*eta1 + a2*eta2 + a3*eta3;
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
    const double eta1 = stratified_apparent_viscosity (&rheology1, gd);
    const double eta2 = stratified_apparent_viscosity (&rheology2, gd);
    const double eta3 = stratified_apparent_viscosity (&rheology3, gd);
    const double rhof = a1*rho1 + a2*rho2 + a3*rho3;

    alphav.x[] = fm.x[]/rhof;
    face vector muv = mu;
#if STRATIFIED_HARMONIC_FACE_VISCOSITY
    const double inv_eta = a1/max(eta1, DBL_MIN) +
      a2/max(eta2, DBL_MIN) + a3/max(eta3, DBL_MIN);
    muv.x[] = fm.x[]/(inv_eta + DBL_MIN);
#else
    muv.x[] = fm.x[]*(a1*eta1 + a2*eta2 + a3*eta3);
#endif
  }
}

#endif
