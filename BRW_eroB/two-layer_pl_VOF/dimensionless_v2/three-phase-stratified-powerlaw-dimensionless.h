/**
# Dimensionless stratified three-phase generalized-Newtonian properties

This header implements the one-fluid material fields for the dimensionless
problem described in the audited full-Navier--Stokes report.  The reference
scales are

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

Material-property jumps can optionally be smeared with Basilisk's standard
vertex-average filter. Set `[numerical_options] filtered = true` in
`case_parameters.ini` and regenerate `base_state/generated_case.h` to filter
both cumulative VOF fields before reconstructing the three material fractions.
The default `false` setting uses the unsmeared VOF fields. Geometric VOF and
conservative momentum advection always use the original, unsmeared fractions.

Four rheological options are available through the integer `model` member:

    0  smooth power law
       eta = Lambda (gamma^2 + epsilon^2)^((n-1)/2)

    1  bounded ideal power law
       eta = Lambda max(|gamma|,gamma_floor)^(n-1)

    2  Papanastasiou Herschel--Bulkley
       eta = Lambda gamma^(n-1)
             + tau0 [1-exp(-m gamma)]/gamma

    3  bounded Herschel--Bulkley
       eta = Lambda gamma^(n-1) + tau0/gamma.

All variables in these expressions are dimensionless.  The strain-rate
invariant is gamma = sqrt(2 D:D), which reduces to |du/dz| in parallel shear.
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
#ifndef CASE_FILTERED
# define CASE_FILTERED 0
#endif
#ifndef FILTERED
# define FILTERED CASE_FILTERED
#endif
#if FILTERED != 0 && FILTERED != 1
# error "[numerical_options] filtered must generate either 0 or 1"
#endif

#define STRATIFIED_SMOOTH_POWER_LAW 0
#define STRATIFIED_BOUNDED_POWER_LAW 1
#define STRATIFIED_PAPANASTASIOU_HB 2
#define STRATIFIED_BOUNDED_HB 3

/** `f` is total liquid, which also makes the free surface visible to GfsView. */
scalar f[], f_lower[];
scalar * interfaces = {f_lower, f};

/**
Smeared cumulative material fractions.  Filtering cumulative fields rather
than three independent phase fractions preserves the stratified construction:
the filtered physical fractions are `sf_lower`, `sf - sf_lower` and `1 - sf`.
The positive filter stencil preserves `sf_lower <= sf` whenever the advected
VOF fields satisfy `f_lower <= f`.
*/
#if FILTERED
scalar sf[], sf_lower[];
#else
# define sf f
# define sf_lower f_lower
#endif

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
face vector gammaDotf[];
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

/**
Apply the vertex-average filter used by Basilisk's `two-phase-generic.h` to
both cumulative VOF fields.  The filter is used for material properties only.
*/
static void stratified_filter_material_fractions (void)
{
#if FILTERED
  foreach() {
    sf_lower[] = (4.*f_lower[] +
                  2.*(f_lower[0,1] + f_lower[0,-1] +
                      f_lower[1,0] + f_lower[-1,0]) +
                  f_lower[-1,-1] + f_lower[1,-1] +
                  f_lower[1,1] + f_lower[-1,1])/16.;
    sf[] = (4.*f[] +
            2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
            f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
  }
#if TREE
  set_prolongation (sf_lower, refine_bilinear);
  set_prolongation (sf, refine_bilinear);
#endif
  boundary ({sf_lower, sf});
#endif
}

/**
Compute gamma = sqrt(2 D:D) directly at faces using Vatsal Sanjay's
face-centered stencil.  Normal derivatives are centered on the face and the
two tangential derivatives are averaged from the adjacent cell centers.  The
cell scalar is only a diagnostic average of the four surrounding face values;
the constitutive law for the implicit viscous solve uses `gammaDotf` directly.
*/
static void stratified_update_strain_rate (void)
{
  boundary ((scalar *){u});

  foreach_face(x) {
    const double Dxx = (u.x[] - u.x[-1,0])/Delta;
    const double Dyy = (u.y[0,1] - u.y[0,-1] +
                        u.y[-1,1] - u.y[-1,-1])/(4.*Delta);
    const double dvdx = (u.y[] - u.y[-1,0])/Delta;
    const double dudy = (u.x[0,1] - u.x[0,-1] +
                         u.x[-1,1] - u.x[-1,-1])/(4.*Delta);
    const double Dxy = 0.5*(dvdx + dudy);
    gammaDotf.x[] = sqrt(max(0., 2.*(sq(Dxx) + sq(Dyy) +
                                     2.*sq(Dxy))));
  }

  foreach_face(y) {
    const double Dyy = (u.y[] - u.y[0,-1])/Delta;
    const double Dxx = (u.x[1,0] - u.x[-1,0] +
                        u.x[1,-1] - u.x[-1,-1])/(4.*Delta);
    const double dudy = (u.x[] - u.x[0,-1])/Delta;
    const double dvdx = (u.y[1,0] - u.y[-1,0] +
                         u.y[1,-1] - u.y[-1,-1])/(4.*Delta);
    const double Dxy = 0.5*(dudy + dvdx);
    gammaDotf.y[] = sqrt(max(0., 2.*(sq(Dxx) + sq(Dyy) +
                                     2.*sq(Dxy))));
  }

  boundary ((scalar *){gammaDotf});
  foreach()
    gammaDot[] = (gammaDotf.x[] + gammaDotf.x[1,0] +
                  gammaDotf.y[] + gammaDotf.y[0,1])/4.;
  boundary ({gammaDot});
}

/** Lagged generalized-Newtonian update in dimensionless variables. */
event properties (i++)
{
  stratified_filter_material_fractions();
  stratified_update_strain_rate();

  foreach_face() {
    const double lower = (sf_lower[] + sf_lower[-1])/2.;
    const double liquid = (sf[] + sf[-1])/2.;
    double a1, a2, a3;
    stratified_phase_values (lower, liquid, &a1, &a2, &a3);

    const double gd = gammaDotf.x[];
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

  foreach() {
    double m1, m2, m3;
    stratified_phase_values (sf_lower[], sf[], &m1, &m2, &m3);
    const double eta1 = stratified_apparent_viscosity (&rheology1,
                                                        gammaDot[]);
    const double eta2 = stratified_apparent_viscosity (&rheology2,
                                                        gammaDot[]);
    const double eta3 = stratified_apparent_viscosity (&rheology3,
                                                        gammaDot[]);
    muMix[] = m1*eta1 + m2*eta2 + m3*eta3;
    rhov[] = cm[]*(m1*rho1 + m2*rho2 + m3*rho3);

    double a1, a2, a3;
    stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
    f1[] = a1;
    f2[] = a2;
    f3[] = a3;
  }
  boundary ({muMix, rhov, f1, f2, f3});

#if FILTERED && TREE
  set_prolongation (sf_lower, fraction_refine);
  set_prolongation (sf, fraction_refine);
#endif
}

#endif
