/**
# Dimensionless two-layer power-law roll waves with three-phase VOF and AMR

This case reproduces the functionality of the uploaded dimensional solver but
advances the nondimensional one-fluid equations. The lower steady-uniform
depth and lower steady depth-averaged velocity are the reference length and
velocity. All material properties, pressure, gravity, surface tension, time,
mesh dimensions and output intervals entering the PDE are dimensionless.

The cumulative VOF functions are

    f_lower = alpha_lower,
    f       = alpha_lower + alpha_upper,

so alpha_1=f_lower, alpha_2=f-f_lower and alpha_3=1-f. The base-state
generator reads `case_parameters.ini`, derives Lambda_lower from the
steady-uniform compatibility condition and writes both the analytical profiles
and `base_state/generated_case.h`.

The solver retains nested VOF transport, total-momentum-consistent phase
advection, quadtree AMR, an embedded free-slip ceiling, optional capillarity,
optional phase-selective streamwise gravity, a divergence-free velocity
perturbation, diagnostics, Fourier amplitude tracking, vertical slices,
interface facets, GFS output and restart dumps.
*/

#include "grid/quadtree.h"
#include "adapt_wavelet_leave_interface_limited.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "base_state/generated_case.h"
#include "three-phase-stratified-powerlaw-dimensionless.h"
#include "tension.h"
#include "conserving3f-stratified.h"
#include "profile1d.h"
#include "output.h"
#include "gfs_output_safe.h"
#include "navier-stokes/perfs.h"

/* =============================================================================
   GENERATED DIMENSIONLESS CASE PARAMETERS
   =============================================================================
   Edit `case_parameters.ini` and run `generate_dimensionless_base_state.py`.
   Do not edit the macros below or `base_state/generated_case.h` manually.
*/

#define H1 1.0
#define H2 CASE_DEPTH_RATIO
#define LIQUID_DEPTH (H1 + H2)
#define LX CASE_LX
#define NOMINAL_HEIGHT CASE_ALIGNED_CEILING

/* Equivalent gravity magnitude. Its bed-aligned components are
   G_ACCEL*SIN_THETA=S_0/Fr_l^2 and G_ACCEL*COS_THETA=1/Fr_l^2. */
#define G_ACCEL (sqrt(1. + sq(CASE_SLOPE_TAN))/sq(CASE_FROUDE))
#define SLOPE_TAN CASE_SLOPE_TAN

#define RHOLOWERLAYER 1.0
#define RHOUPPERLAYER CASE_DENSITY_RATIO
#define AIRRHO CASE_AIR_DENSITY_RATIO

#define K_LOWER CASE_LAMBDA_LOWER
#define N_LOWER CASE_LOWER_N
#define K_UPPER CASE_LAMBDA_UPPER
#define N_UPPER CASE_UPPER_N
#define MUMAXLOWER CASE_LOWER_ETA_MAX
#define MUMAXUPPER CASE_UPPER_ETA_MAX
#define AIRMU CASE_AIR_ETA
#define GAMMA_MIN min(CASE_LOWER_GAMMA_FLOOR, CASE_UPPER_GAMMA_FLOOR)

#define MAXLEVEL CASE_MAXLEVEL
#define MINLEVEL CASE_MINLEVEL
#define INITLEVEL CASE_INITLEVEL
#define END_TIME CASE_END_TIME
#define MAX_DT CASE_MAX_DT
#define CFL_NUMBER CASE_CFL
#define SOLVER_TOLERANCE CASE_SOLVER_TOLERANCE
#define SOLVER_NITERMAX CASE_SOLVER_NITERMAX

#define OUTPUT_DT CASE_OUTPUT_DT
#define GFS_DT CASE_GFS_DT
#define DUMP_DT CASE_DUMP_DT
#define ENABLE_DUMP_OUTPUT CASE_ENABLE_DUMP_OUTPUT
#define WAVE_AMPLITUDE_EVERY CASE_WAVE_AMPLITUDE_EVERY
#define PERTURBATION_MODE_NUMBER CASE_PERTURBATION_MODE
#define WAVE_MODE_NUMBER CASE_TRACKED_MODE
#define INITIAL_AUDIT_WARNING_TOL CASE_INITIAL_AUDIT_WARNING_TOL

#define VELOCITY_PERTURBATION_AMPLITUDE CASE_VELOCITY_AMPLITUDE
#define LOWER_DEPTH_PERTURBATION_AMPLITUDE CASE_LOWER_DEPTH_AMPLITUDE
#define UPPER_DEPTH_PERTURBATION_AMPLITUDE CASE_UPPER_DEPTH_AMPLITUDE
#define PERTURBATION_PHASE CASE_PERTURBATION_PHASE
#define INITIAL_PRESSURE_MODE CASE_PRESSURE_MODE

#define U_ERR CASE_U_ERROR_LOWER
#define V_ERR CASE_V_ERROR_LOWER
#define U_ERRAIR CASE_U_ERROR_AIR
#define V_ERRAIR CASE_V_ERROR_AIR
#define OMEGA_ERR CASE_VORTICITY_ERROR
#define SHEAR_ERR CASE_SHEAR_ERROR
#define INTERFACE_PADDING CASE_INTERFACE_PADDING
#define INTERFACE_EPS CASE_INTERFACE_EPS
#define CEILING_FINE_FRACTION CASE_CEILING_FINE_FRACTION
#define CEILING_LEVEL_DROP CASE_CEILING_LEVEL_DROP
#define INIT_ADAPT_PASSES CASE_INIT_ADAPT_PASSES
#define BED_FINE_HEIGHT (H1*(1. + fabs(LOWER_DEPTH_PERTURBATION_AMPLITUDE)))

#define SIGMA_INTERNAL CASE_SIGMA_INTERNAL
#define SIGMA_FREE CASE_SIGMA_FREE
#define AIR_STREAMWISE_GRAVITY CASE_AIR_STREAMWISE_GRAVITY

#define BASE_STATE_DIR "base_state"
#define BASE_VELOCITY_FILE "velocity.dat"
#define BASE_PRESSURE_FILE "pressure.dat"

/* Optional output conversion only. These scales never enter the PDE. */
#define DIMENSIONAL_REFERENCE CASE_REFERENCE_ENABLED
#define HREF_M CASE_HREF_M
#define UREF_MS CASE_UREF_MS
#define TREF_S CASE_TREF_S
#define PREF_PA CASE_PREF_PA
#define ETAREF_PAS CASE_ETAREF_PAS

/* =============================================================================
   END GENERATED CASE PARAMETERS
   ============================================================================= */

#define SIN_THETA (SLOPE_TAN/sqrt(1. + sq(SLOPE_TAN)))
#define COS_THETA (1./sqrt(1. + sq(SLOPE_TAN)))

face vector gravity_accel[];

Profile1D profile_u0, profile_p;

scalar omega_amr[], shear_amr[];
scalar speed_lower[], speed_upper[], speed_air[];
scalar level_field[];

static double ceiling_y = NOMINAL_HEIGHT;
static double u_reference = 1.;
static int last_nf = 0, last_nc = 0;
static FILE * wave_amplitude_fp = NULL;

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
uf.n[bottom] = 0.;
uf.n[top] = 0.;

/* Free-slip embedded ceiling, following the demonstrated Basilisk
   component-wise homogeneous-Neumann pattern for embedded boundaries. */
u.n[embed] = neumann(0.);
u.t[embed] = neumann(0.);
// The pressure boundary is left to centered.h's acceleration-consistent rule.

static inline double finest_delta (void)
{
  return LX/(double)(1 << MAXLEVEL);
}

static void read_profile_from_dir (Profile1D * p,
                                   const char * directory,
                                   const char * basename)
{
  char path[1024];
  snprintf (path, sizeof(path), "%s/%s", directory, basename);
  profile1d_read (p, path);
  fprintf (ferr, "# read %s: %d values from y=%g to y=%g\n",
           path, p->n, p->x[0], p->x[p->n - 1]);
}

static double expected_interface_velocity (void)
{
  return CASE_U_INTERFACE;
}

static double expected_surface_velocity (void)
{
  return CASE_U_SURFACE;
}

static double expected_bed_pressure (void)
{
  return CASE_P_BED;
}

static void validate_case_parameters (void)
{
  if (MAXLEVEL < MINLEVEL || INITLEVEL < MINLEVEL || INITLEVEL > MAXLEVEL) {
    fprintf (ferr,
             "Invalid mesh levels: MINLEVEL=%d INITLEVEL=%d MAXLEVEL=%d.\n",
             MINLEVEL, INITLEVEL, MAXLEVEL);
    exit (1);
  }
  if (G_ACCEL <= 0. || SLOPE_TAN <= 0. || H1 <= 0. || H2 <= 0. ||
      LX <= 0. || NOMINAL_HEIGHT <= H1 + H2 ||
      RHOLOWERLAYER <= 0. || RHOUPPERLAYER <= 0. || AIRRHO <= 0. ||
      K_LOWER <= 0. || K_UPPER <= 0. || N_LOWER <= 0. || N_UPPER <= 0. ||
      MUMAXLOWER <= 0. || MUMAXUPPER <= 0. || AIRMU <= 0. ||
      GAMMA_MIN <= 0. || SOLVER_TOLERANCE <= 0. || SOLVER_NITERMAX <= 0 ||
      MAX_DT <= 0. || CFL_NUMBER <= 0.) {
    fprintf (ferr,
             "One or more geometry/material/solver parameters are invalid. "
             "Check case_parameters.ini and regenerate the base state.\n");
    exit (1);
  }
  if (WAVE_AMPLITUDE_EVERY <= 0 || WAVE_MODE_NUMBER <= 0 ||
      PERTURBATION_MODE_NUMBER <= 0) {
    fprintf (ferr,
             "Wave-output interval and mode numbers must be positive: "
             "every=%d perturbationMode=%d trackedMode=%d.\n",
             WAVE_AMPLITUDE_EVERY, PERTURBATION_MODE_NUMBER,
             WAVE_MODE_NUMBER);
    exit (1);
  }
  if (fabs(VELOCITY_PERTURBATION_AMPLITUDE) >= 1.) {
    fprintf (ferr,
             "The optional streamfunction velocity amplitude must satisfy "
             "|A_u| < 1.\n");
    exit (1);
  }
}

static void validate_simple_text_profiles (void)
{
  if (profile_u0.n != profile_p.n) {
    fprintf (ferr,
             "Simple profiles have different row counts: velocity=%d "
             "pressure=%d\n",
             profile_u0.n, profile_p.n);
    exit (1);
  }

  const double xtol = 1.e-10*max(1., ceiling_y);
  for (int j = 0; j < profile_u0.n; j++)
    if (fabs(profile_u0.x[j] - profile_p.x[j]) > xtol) {
      fprintf (ferr,
               "Simple profiles use different y coordinates at row %d.\n",
               j + 1);
      exit (1);
    }

  if (fabs(profile_u0.x[0]) > xtol ||
      fabs(profile_u0.x[profile_u0.n - 1] - ceiling_y) > xtol) {
    fprintf (ferr,
             "Profiles must cover y=0 to the compiled ceiling_y=%g; "
             "got [%g,%g].\n"
             "Regenerate the dimensionless base state, then "
             "regenerate the base state.\n",
             ceiling_y, profile_u0.x[0], profile_u0.x[profile_u0.n - 1]);
    exit (1);
  }

  const double ubed = profile1d_eval (&profile_u0, 0.);
  const double ui_file = profile1d_eval (&profile_u0, H1);
  const double us_file = profile1d_eval (&profile_u0, LIQUID_DEPTH);
  const double uceil = profile1d_eval (&profile_u0, ceiling_y);
  const double pbed_file = profile1d_eval (&profile_p, 0.);
  const double pceil = profile1d_eval (&profile_p, ceiling_y);
  const double ui_expected = expected_interface_velocity();
  const double us_expected = expected_surface_velocity();
  const double pbed_expected = expected_bed_pressure();

  double air_umin = HUGE, air_umax = -HUGE;
  for (int j = 0; j < profile_u0.n; j++)
    if (profile_u0.x[j] >= LIQUID_DEPTH - xtol) {
      air_umin = min(air_umin, profile_u0.value[j]);
      air_umax = max(air_umax, profile_u0.value[j]);
    }

  const double utol = 2.e-8*max(1., fabs(us_expected));
  const double ptol = 2.e-8*max(1., fabs(pbed_expected));
  if (fabs(ubed) > utol || fabs(uceil - us_expected) > utol ||
      fabs(pceil) > ptol || fabs(ui_file - ui_expected) > utol ||
      fabs(us_file - us_expected) > utol ||
      fabs(pbed_file - pbed_expected) > ptol ||
      air_umax - air_umin > utol) {
    fprintf (ferr,
             "Analytical profile is not the compatible passive-ambient "
             "base state.\n"
             "  interface velocity: file=%g expected=%g\n"
             "  surface velocity:   file=%g expected=%g\n"
             "  ceiling velocity:   file=%g expected=%g\n"
             "  ambient spread:     %g\n"
             "  bed pressure:       file=%g expected=%g\n"
             "  bed velocity / ceiling pressure: %g / %g\n",
             ui_file, ui_expected, us_file, us_expected,
             uceil, us_expected, air_umax - air_umin,
             pbed_file, pbed_expected, ubed, pceil);
    exit (1);
  }

  u_reference = us_file;
  fprintf (ferr,
           "# compatible base state validated: U_interface=%g U_surface=%g "
           "U_ambient=%g p_bed=%g ceiling=%g\n",
           ui_file, us_file, uceil, pbed_file, ceiling_y);
  fprintf (ferr,
           "# groups: Fr_l=%g S0=%g h_r=%g r_rho=%g kappa_K=%g "
           "Lambda_l=%g Lambda_u=%g epsilon_l=%g epsilon_u=%g "
           "r_air=%g eta_air=%g\n",
           (double)CASE_FROUDE, (double)CASE_SLOPE_TAN,
           (double)CASE_DEPTH_RATIO, (double)CASE_DENSITY_RATIO,
           (double)CASE_CONSISTENCY_RATIO, (double)CASE_LAMBDA_LOWER,
           (double)CASE_LAMBDA_UPPER, (double)CASE_LOWER_EPSILON,
           (double)CASE_UPPER_EPSILON, (double)CASE_AIR_DENSITY_RATIO,
           (double)CASE_AIR_ETA);
}

/** Periodic disturbance shared by the depth and velocity perturbations. */
static inline double disturbance_sine (double xx)
{
  return sin(2.*pi*PERTURBATION_MODE_NUMBER*xx/LX +
             PERTURBATION_PHASE);
}

static inline double initial_lower_depth (double xx)
{
  return H1*(1. + LOWER_DEPTH_PERTURBATION_AMPLITUDE*
                    disturbance_sine(xx));
}

static inline double initial_upper_depth (double xx)
{
  return H2*(1. + UPPER_DEPTH_PERTURBATION_AMPLITUDE*
                    disturbance_sine(xx));
}

static inline double initial_internal_interface (double xx)
{
  return initial_lower_depth(xx);
}

static inline double initial_free_surface (double xx)
{
  return initial_lower_depth(xx) + initial_upper_depth(xx);
}

/** Compatible base velocity plus an optional divergence-free perturbation.
    The perturbation comes from

      psi' = A Uref yc/(2 pi) [1 - cos(2 pi y/yc)] sin(kx + phi),

    with u'_x = d psi'/dy and u'_y = -d psi'/dx.  It vanishes at the bed
    and ceiling and causes no sinusoidal variation of total discharge. */
static inline void initial_velocity_components (double xx, double yy,
                                                double * ux, double * uy)
{
  const double mode = PERTURBATION_MODE_NUMBER;
  const double phase = 2.*pi*mode*xx/LX + PERTURBATION_PHASE;
  const double sy = clamp(yy/(ceiling_y + SEPS), 0., 1.);
  const double amp = VELOCITY_PERTURBATION_AMPLITUDE*u_reference;

  /* psi' = A Uref yc/(2 pi) [1 - cos(2 pi y/yc)] sin(kx + phi).
     This gives zero normal velocity at y=0 and y=yc, zero modal total
     discharge, and exactly divergence-free analytical perturbations. */
  *ux = profile1d_eval(&profile_u0, yy) +
    amp*sin(2.*pi*sy)*sin(phase);
  *uy = -amp*mode*(ceiling_y/LX)*
    (1. - cos(2.*pi*sy))*cos(phase);
}

/** Piecewise hydrostatic pressure for the locally perturbed interfaces,
    with dimensionless p=0 at the embedded ceiling. */
static inline double local_hydrostatic_pressure (double xx, double yy)
{
  const double h1x = initial_lower_depth(xx);
  const double h2x = initial_upper_depth(xx);
  const double eta2 = h1x + h2x;
  const double gn = G_ACCEL*COS_THETA;

  if (yy > eta2)
    return AIRRHO*gn*(ceiling_y - yy);
  if (yy > h1x)
    return AIRRHO*gn*(ceiling_y - eta2) +
      RHOUPPERLAYER*gn*(eta2 - yy);
  return AIRRHO*gn*(ceiling_y - eta2) +
    RHOUPPERLAYER*gn*h2x + RHOLOWERLAYER*gn*(h1x - yy);
}

static inline void build_embedded_ceiling (void)
{
  solid (cs, fs, ceiling_y - y);
  fractions_cleanup (cs, fs);
}

/** Maximum level allowed by position.  The wavelet criterion can coarsen
    the interior, while mixed f_lower, f and cs cells are protected by the
    custom adaptation function. */
static int local_maxlevel_cell (double xx, double yy, double zz,
                                double delta)
{
  (void) xx; (void) zz;
  /* Apply the requested y <= 0.975*ceiling rule to the complete cell
     footprint, not only its centre.  This prevents a coarse cut cell whose
     centre lies just below the threshold from being refined back to
     MAXLEVEL merely to represent the horizontal ceiling. */
  if (yy + 0.5*delta <= ceiling_y*CEILING_FINE_FRACTION)
    return MAXLEVEL;
  return max(MINLEVEL, MAXLEVEL - CEILING_LEVEL_DROP);
}

static int local_minlevel (double xx, double yy, double zz)
{
  (void) xx; (void) zz;
  const double dmin = finest_delta();
  if (yy <= BED_FINE_HEIGHT + 1.5*dmin)
    return MAXLEVEL;
  return MINLEVEL;
}

static inline void force_required_refinement (void)
{
  const double dmin = finest_delta();

  /* Resolve the lower liquid/basal high-shear region at the finest level. */
  refine (y <= 4.5*dmin  && level < MAXLEVEL);
}

static void validate_initialization_parameters (void)
{
  if (fabs(LOWER_DEPTH_PERTURBATION_AMPLITUDE) >= 1. ||
      fabs(UPPER_DEPTH_PERTURBATION_AMPLITUDE) >= 1.) {
    fprintf (ferr,
             "Depth-perturbation amplitudes must have absolute value < 1.\n");
    exit (1);
  }
  const double max_surface =
    H1*(1. + fabs(LOWER_DEPTH_PERTURBATION_AMPLITUDE)) +
    H2*(1. + fabs(UPPER_DEPTH_PERTURBATION_AMPLITUDE));
  if (max_surface >= ceiling_y*CEILING_FINE_FRACTION) {
    fprintf (ferr,
             "The maximum initial free surface (%g) reaches the ceiling "
             "coarsening band (%g). Reduce depth amplitude or raise the "
             "ceiling.\n",
             max_surface, ceiling_y*CEILING_FINE_FRACTION);
    exit (1);
  }
  if (max(fabs(LOWER_DEPTH_PERTURBATION_AMPLITUDE),
          fabs(UPPER_DEPTH_PERTURBATION_AMPLITUDE)) > 0.02)
    fprintf (ferr,
             "WARNING: depth perturbation exceeds 2%%; this is not a linear "
             "stability initialization.\n");
}

static void set_initial_fields (void)
{
  /* Both cumulative VOF fields receive the sinusoidal depth disturbance. */
  fraction (f_lower, initial_internal_interface(x) - y);
  fraction (f, initial_free_surface(x) - y);

  foreach() {
    initial_velocity_components (x, y, &u.x[], &u.y[]);
#if INITIAL_PRESSURE_MODE == 1
    p[] = local_hydrostatic_pressure(x, y);
#else
    p[] = profile1d_eval(&profile_p, y);
#endif
  }

  stratified_enforce_hierarchy();
  boundary ((scalar *){u, p});
}

static void update_amr_indicators (void)
{
  stratified_update_strain_rate();
  vorticity (u, omega_amr);
  foreach() {
    shear_amr[] = gammaDot[];
    const double speed = sqrt(sq(u.x[]) + sq(u.y[]));
    double a1, a2, a3;
    stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
    speed_lower[] = a1*speed;
    speed_upper[] = a2*speed;
    speed_air[] = a3*speed;
    level_field[] = level;
  }
  boundary ({omega_amr, shear_amr, speed_lower, speed_upper,
             speed_air, level_field});
}

static astats perform_adaptation (void)
{
  update_amr_indicators();
  Adapt_limited controls = {
    .slist = (scalar *){omega_amr, speed_lower, speed_upper, speed_air},
    .vol_frac = (scalar *){f_lower, f, cs},
    .max = (double []){OMEGA_ERR, U_ERR, U_ERR, U_ERRAIR},
    .MLFun = NULL,
    .MLFunDelta = local_maxlevel_cell,
    .MinFun = local_minlevel,
    .minlevel = MINLEVEL,
    .padding = INTERFACE_PADDING,
    .interface_eps = INTERFACE_EPS,
    .list = all
  };
  astats st = adapt_wavelet_limited (controls);

  build_embedded_ceiling();
  stratified_enforce_hierarchy();
  last_nf = st.nf;
  last_nc = st.nc;
  return st;
}

typedef struct {
  double mean_internal, mean_free, mean_upper;
  double cos_internal, sin_internal, amplitude_internal;
  double cos_free, sin_free, amplitude_free;
  double cos_upper, sin_upper, amplitude_upper;
} WaveMetrics;

/** Fourier amplitudes from conservative VOF volume integrals. */
static WaveMetrics compute_wave_metrics (void)
{
  const double k = 2.*pi*WAVE_MODE_NUMBER/LX;
  double vl = 0., vf = 0.;
  double cl = 0., sl = 0., cf = 0., sfree = 0.;
  foreach (reduction(+:vl) reduction(+:vf)
           reduction(+:cl) reduction(+:sl)
           reduction(+:cf) reduction(+:sfree)) {
    const double vol = dv();
    const double c = cos(k*x), ss = sin(k*x);
    vl += f_lower[]*vol;
    vf += f[]*vol;
    cl += f_lower[]*c*vol;
    sl += f_lower[]*ss*vol;
    cf += f[]*c*vol;
    sfree += f[]*ss*vol;
  }
  const double fac = 2./LX;
  WaveMetrics w = {0};
  w.mean_internal = vl/LX;
  w.mean_free = vf/LX;
  w.mean_upper = (vf - vl)/LX;
  w.cos_internal = fac*cl;
  w.sin_internal = fac*sl;
  w.cos_free = fac*cf;
  w.sin_free = fac*sfree;
  w.cos_upper = w.cos_free - w.cos_internal;
  w.sin_upper = w.sin_free - w.sin_internal;
  w.amplitude_internal = hypot(w.cos_internal, w.sin_internal);
  w.amplitude_free = hypot(w.cos_free, w.sin_free);
  w.amplitude_upper = hypot(w.cos_upper, w.sin_upper);
  return w;
}

/** Audit the initial divergence and streamwise discharge mode. */
static void write_initial_condition_audit (void)
{
  boundary ((scalar *){u});
  const double k = 2.*pi*WAVE_MODE_NUMBER/LX;
  double divmax = 0., div2 = 0., volume = 0.;
  double qcos = 0., qsin = 0.;
  foreach (reduction(max:divmax) reduction(+:div2) reduction(+:volume)
           reduction(+:qcos) reduction(+:qsin)) {
    const double vol = dv();
    qcos += u.x[]*cos(k*x)*vol;
    qsin += u.x[]*sin(k*x)*vol;
    if (cs[] > 1. - 1.e-10 && y > 1.5*Delta &&
        y < ceiling_y - 1.5*Delta) {
      const double divu =
        (u.x[1,0] - u.x[-1,0] + u.y[0,1] - u.y[0,-1])/(2.*Delta);
      divmax = max(divmax, fabs(divu));
      div2 += sq(divu)*vol;
      volume += vol;
    }
  }
  qcos *= 2./LX;
  qsin *= 2./LX;
  const double divrms = sqrt(div2/(volume + SEPS));
  const double divscale = max(u_reference/H1, SEPS);
  const WaveMetrics w = compute_wave_metrics();

  if (pid() == 0) {
    FILE * fp = fopen ("initial_condition_audit.tsv", "w");
    if (!fp) {
      perror("initial_condition_audit.tsv");
      exit (1);
    }
    fprintf (fp,
             "divMax_star\tdivRMS_star\tdivMaxScaled\tQcos_star"
             "\tQsin_star\tAinternal_star\tAfree_star\tAupper_star\n");
    fprintf (fp, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             divmax, divrms, divmax/divscale, qcos, qsin,
             w.amplitude_internal, w.amplitude_free, w.amplitude_upper);
    fclose (fp);
  }
  fprintf (ferr,
           "# initial audit: divMax=%g divRMS=%g scaledDivMax=%g "
           "Qmode=(%g,%g) Ainternal=%g Afree=%g Aupper=%g\n",
           divmax, divrms, divmax/divscale, qcos, qsin,
           w.amplitude_internal, w.amplitude_free, w.amplitude_upper);
  if (divmax/divscale > INITIAL_AUDIT_WARNING_TOL)
    fprintf (ferr,
             "WARNING: initialized velocity divergence exceeds scaled tolerance %g.\n",
             (double)INITIAL_AUDIT_WARNING_TOL);
}

int main (void)
{
  validate_case_parameters();
  size (LX);
  origin (0., 0.);
  periodic (right);

  /* Align the embedded ceiling with a MAXLEVEL cell face, avoiding a row
     of arbitrarily small cut cells. */
  const double dmin = finest_delta();
  ceiling_y = ceil(NOMINAL_HEIGHT/dmin - 1.e-12)*dmin;

  rho1 = RHOLOWERLAYER;
  rho2 = RHOUPPERLAYER;
  rho3 = AIRRHO;

  rheology1 = (StratifiedRheology) {
    .model = CASE_LOWER_MODEL,
    .Lambda = CASE_LAMBDA_LOWER,
    .n = CASE_LOWER_N,
    .epsilon = CASE_LOWER_EPSILON,
    .tau0 = CASE_LOWER_TAU0,
    .m = CASE_LOWER_M,
    .eta_min = CASE_LOWER_ETA_MIN,
    .eta_max = CASE_LOWER_ETA_MAX,
    .gamma_floor = CASE_LOWER_GAMMA_FLOOR
  };
  rheology2 = (StratifiedRheology) {
    .model = CASE_UPPER_MODEL,
    .Lambda = CASE_LAMBDA_UPPER,
    .n = CASE_UPPER_N,
    .epsilon = CASE_UPPER_EPSILON,
    .tau0 = CASE_UPPER_TAU0,
    .m = CASE_UPPER_M,
    .eta_min = CASE_UPPER_ETA_MIN,
    .eta_max = CASE_UPPER_ETA_MAX,
    .gamma_floor = CASE_UPPER_GAMMA_FLOOR
  };
  rheology3 = (StratifiedRheology) {
    .model = STRATIFIED_BOUNDED_POWER_LAW,
    .Lambda = CASE_AIR_ETA,
    .n = 1.,
    .epsilon = 0.,
    .tau0 = 0.,
    .m = 0.,
    .eta_min = CASE_AIR_ETA,
    .eta_max = CASE_AIR_ETA,
    .gamma_floor = min(CASE_LOWER_GAMMA_FLOOR,
                       CASE_UPPER_GAMMA_FLOOR)
  };

  f_lower.sigma = SIGMA_INTERNAL;
  f.sigma = SIGMA_FREE;

  a = gravity_accel;
  CFL = CFL_NUMBER;
  DT = MAX_DT;
  TOLERANCE = SOLVER_TOLERANCE;
  NITERMAX = SOLVER_NITERMAX;

  init_grid (1 << INITLEVEL);
  run();
}

event init (i = 0)
{
  if (pid() == 0)
    system ("mkdir -p gfs dumps interfaces slices");
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  read_profile_from_dir (&profile_u0, BASE_STATE_DIR, BASE_VELOCITY_FILE);
  read_profile_from_dir (&profile_p, BASE_STATE_DIR, BASE_PRESSURE_FILE);
  validate_simple_text_profiles();

  build_embedded_ceiling();
  validate_initialization_parameters();
  force_required_refinement();
  build_embedded_ceiling();

  for (int pass = 0; pass < INIT_ADAPT_PASSES; pass++) {
    set_initial_fields();
    astats st = perform_adaptation();
    fprintf (ferr, "# init AMR pass %d: nf=%d nc=%d cells=%ld\n",
             pass + 1, st.nf, st.nc, grid->tn);
    if (st.nf == 0 && st.nc == 0)
      break;
  }
  set_initial_fields();
  update_amr_indicators();
  write_initial_condition_audit();

  const double we_free = SIGMA_FREE > 0. ? 1./SIGMA_FREE : HUGE;
  const double we_internal = SIGMA_INTERNAL > 0. ? 1./SIGMA_INTERNAL : HUGE;

  fprintf (ferr,
           "# dimensionless AMR roll wave: Lx=%g nominalHeight=%g embeddedCeiling=%g "
           "dmin=%g MAXLEVEL=%d MINLEVEL=%d FILTERED=%d activeCells=%ld\n",
           LX, (double)NOMINAL_HEIGHT, ceiling_y, finest_delta(),
           MAXLEVEL, MINLEVEL, (int)FILTERED, grid->tn);
  fprintf (ferr,
           "# AMR policy: full MAXLEVEL below y=%g; ceiling band max=%d; "
           "bed fine height=%g\n",
           ceiling_y*CEILING_FINE_FRACTION,
           max(MINLEVEL, MAXLEVEL - CEILING_LEVEL_DROP),
           (double)BED_FINE_HEIGHT);
  fprintf (ferr,
           "# cells through H1=%g; sigma12=%g sigma23=%g "
           "We12=%g We23=%g velocityAmp=%g lowerDepthAmp=%g "
           "upperDepthAmp=%g perturbationMode=%d trackedMode=%d "
           "amplitudeEvery=%d\n",
           H1/finest_delta(), (double)SIGMA_INTERNAL, (double)SIGMA_FREE,
           we_internal, we_free,
           (double)VELOCITY_PERTURBATION_AMPLITUDE,
           (double)LOWER_DEPTH_PERTURBATION_AMPLITUDE,
           (double)UPPER_DEPTH_PERTURBATION_AMPLITUDE,
           PERTURBATION_MODE_NUMBER, WAVE_MODE_NUMBER,
           WAVE_AMPLITUDE_EVERY);
  fprintf (ferr,
           "# columns (physical fields dimensionless): i t dt leaves minLevel maxLevel lowerMin interfaceMin "
           "ceilingMax mgp mgu pRes uRes "
           "V1 V2 V3 Px KE Umax hierarchy bounds dVlower dVliquid corrected "
           "nf nc\n");
}

event acceleration (i++)
{
  face vector av = a;
  const double gs = G_ACCEL*SIN_THETA;
  const double gn = G_ACCEL*COS_THETA;

  foreach_face(x) {
    const double lower = (f_lower[] + f_lower[-1])/2.;
    const double liquid = (f[] + f[-1])/2.;
    double a1, a2, a3;
    stratified_phase_values (lower, liquid, &a1, &a2, &a3);
    const double rhof = a1*rho1 + a2*rho2 + a3*rho3;
#if AIR_STREAMWISE_GRAVITY
    av.x[] += gs;
#else
    av.x[] += gs*(a1*rho1 + a2*rho2)/rhof;
#endif
  }
  foreach_face(y)
    av.y[] -= gn;
}

event adapt (i++)
{
  perform_adaptation();
}

event wave_amplitude (i = 0; i += WAVE_AMPLITUDE_EVERY)
{
  const WaveMetrics w = compute_wave_metrics();
  if (pid() == 0) {
    static double previous_t = -1.;
    static double previous_internal = 0., previous_free = 0., previous_upper = 0.;
    if (!wave_amplitude_fp) {
      wave_amplitude_fp = fopen ("wave_amplitude.tsv", "w");
      if (!wave_amplitude_fp) {
        perror("wave_amplitude.tsv");
        exit (1);
      }
      fprintf (wave_amplitude_fp,
               "i\tt_star\tdt_star\tt_s\tdt_s\tmode"
               "\tmeanInternal_star\tmeanFree_star\tmeanUpper_star"
               "\tcosInternal_star\tsinInternal_star\tAinternal_star"
               "\tAinternal/H1\tcosFree_star\tsinFree_star\tAfree_star"
               "\tAfree/(H1+H2)\tcosUpper_star\tsinUpper_star"
               "\tAupper_star\tAupper/H2\tsigmaInternal_star"
               "\tsigmaFree_star\tsigmaUpper_star"
               "\tAinternal_m\tAfree_m\tAupper_m"
               "\tsigmaInternal_1_s\tsigmaFree_1_s\tsigmaUpper_1_s\n");
    }
    double sigma_internal = NAN, sigma_free = NAN, sigma_upper = NAN;
    if (previous_t >= 0. && t > previous_t) {
      if (previous_internal > 0. && w.amplitude_internal > 0.)
        sigma_internal = log(w.amplitude_internal/previous_internal)/(t - previous_t);
      if (previous_free > 0. && w.amplitude_free > 0.)
        sigma_free = log(w.amplitude_free/previous_free)/(t - previous_t);
      if (previous_upper > 0. && w.amplitude_upper > 0.)
        sigma_upper = log(w.amplitude_upper/previous_upper)/(t - previous_t);
    }
    const double tdim = DIMENSIONAL_REFERENCE ? t*TREF_S : NAN;
    const double dtdim = DIMENSIONAL_REFERENCE ? dt*TREF_S : NAN;
    const double length_scale = DIMENSIONAL_REFERENCE ? HREF_M : NAN;
    const double rate_scale = DIMENSIONAL_REFERENCE ? 1./TREF_S : NAN;
    fprintf (wave_amplitude_fp,
             "%d\t%.16e\t%.16e\t%.16e\t%.16e\t%d"
             "\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\n",
             i, t, dt, tdim, dtdim, WAVE_MODE_NUMBER,
             w.mean_internal, w.mean_free, w.mean_upper,
             w.cos_internal, w.sin_internal, w.amplitude_internal,
             w.amplitude_internal/H1,
             w.cos_free, w.sin_free, w.amplitude_free,
             w.amplitude_free/LIQUID_DEPTH,
             w.cos_upper, w.sin_upper, w.amplitude_upper,
             w.amplitude_upper/H2,
             sigma_internal, sigma_free, sigma_upper,
             w.amplitude_internal*length_scale,
             w.amplitude_free*length_scale,
             w.amplitude_upper*length_scale,
             sigma_internal*rate_scale,
             sigma_free*rate_scale,
             sigma_upper*rate_scale);
    fflush (wave_amplitude_fp);
    previous_t = t;
    previous_internal = w.amplitude_internal;
    previous_free = w.amplitude_free;
    previous_upper = w.amplitude_upper;
  }
}

event diagnostics (i += 10)
{
  stratified_update_phase_fields();
  long leaves = 0;
  int levmin = MAXLEVEL, levmax = 0;
  int lower_min_level = MAXLEVEL, interface_min_level = MAXLEVEL;
  int ceiling_max_level = 0;
  double v1 = 0., v2 = 0., v3 = 0., px = 0., ke = 0., umax = 0.;
  foreach (reduction(+:leaves) reduction(min:levmin) reduction(max:levmax)
           reduction(min:lower_min_level) reduction(min:interface_min_level)
           reduction(max:ceiling_max_level)
           reduction(+:v1) reduction(+:v2) reduction(+:v3)
           reduction(+:px) reduction(+:ke) reduction(max:umax)) {
    leaves++;
    levmin = min(levmin, level);
    levmax = max(levmax, level);
    if (y <= BED_FINE_HEIGHT)
      lower_min_level = min(lower_min_level, level);
    if ((f_lower[] > INTERFACE_EPS && f_lower[] < 1. - INTERFACE_EPS) ||
        (f[] > INTERFACE_EPS && f[] < 1. - INTERFACE_EPS))
      interface_min_level = min(interface_min_level, level);
    if (y + 0.5*Delta > ceiling_y*CEILING_FINE_FRACTION &&
        cs[] > INTERFACE_EPS)
      ceiling_max_level = max(ceiling_max_level, level);
    const double vol = dv();
    const double rhoc = stratified_density_cell(f_lower[], f[]);
    v1 += f1[]*vol;
    v2 += f2[]*vol;
    v3 += f3[]*vol;
    px += rhoc*u.x[]*vol;
    ke += 0.5*rhoc*(sq(u.x[]) + sq(u.y[]))*vol;
    umax = max(umax, sqrt(sq(u.x[]) + sq(u.y[])));
  }

  fprintf (ferr,
           "%d %.12g %.12g %ld %d %d %d %d %d %d %d %.3e %.3e "
           "%.10e %.10e %.10e %.10e %.10e %.8e %.3e %.3e "
           "%.3e %.3e %ld %d %d\n",
           i, t, dt, leaves, levmin, levmax,
           lower_min_level, interface_min_level, ceiling_max_level,
           mgp.i, mgu.i, mgp.resa, mgu.resa, v1, v2, v3, px, ke, umax,
           stratified_hierarchy_max, stratified_bound_max,
           stratified_volume_correction_lower,
           stratified_volume_correction_liquid,
           stratified_corrected_cells, last_nf, last_nc);
}

static void write_vertical_slice (const char * name, double xx)
{
  FILE * fp = fopen (name, "w");
  if (!fp) {
    perror(name);
    exit (1);
  }
  fprintf (fp,
           "z_star\tfLower\tfLiquid\tsfLower\tsfLiquid\tf1\tf2\tf3\tu_x_star\tu_z_star\tp_star\tgamma_star\teta_star\tr_star\tlevel\tcs\n");
  const double ds = finest_delta();
  for (double yy = ds/2.; yy < ceiling_y - ds/4.; yy += ds)
    fprintf (fp,
             "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t"
             "%.12e\t%.12e\t"
             "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
             yy, interpolate(f_lower, xx, yy), interpolate(f, xx, yy),
             interpolate(sf_lower, xx, yy), interpolate(sf, xx, yy),
             interpolate(f1, xx, yy), interpolate(f2, xx, yy),
             interpolate(f3, xx, yy), interpolate(u.x, xx, yy),
             interpolate(u.y, xx, yy), interpolate(p, xx, yy),
             interpolate(gammaDot, xx, yy), interpolate(muMix, xx, yy),
             interpolate(rhov, xx, yy), interpolate(level_field, xx, yy),
             interpolate(cs, xx, yy));
  fclose (fp);
}

static void write_interface_facets (scalar c, const char * stem,
                                    double time)
{
  char local[256];
  snprintf (local, sizeof(local), "interfaces/%s-p%05d.dat", stem, pid());
  FILE * fp = fopen (local, "w");
  output_facets (c, fp);
  fclose (fp);
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  if (pid() == 0) {
    char command[1024];
    snprintf (command, sizeof(command),
              "LC_ALL=C cat interfaces/%s-p*.dat > "
              "interfaces/%s-t%g.dat; rm -f interfaces/%s-p*.dat",
              stem, stem, time, stem);
    system (command);
  }
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
}

event field_outputs (t = 0.; t += OUTPUT_DT; t <= END_TIME + 1.e-12)
{
  stratified_update_phase_fields();
  update_amr_indicators();
  char name[256];
  snprintf (name, sizeof(name), "slices/x0-t%g.tsv", t);
  write_vertical_slice (name, 0.);
  snprintf (name, sizeof(name), "slices/xquarter-t%g.tsv", t);
  write_vertical_slice (name, LX/4.);
  write_interface_facets (f_lower, "internal", t);
  write_interface_facets (f, "free", t);
}

event output_gfs_files (t = 0.; t += GFS_DT; t <= END_TIME + 1.e-12)
{
  stratified_update_phase_fields();
  update_amr_indicators();
  char name[128];
  snprintf (name, sizeof(name), "out-%g.gfs", t);

  /* Keep `f` in the output list because Basilisk output_gfs() declares
     `VariableTracerVOF f` unconditionally.  f_lower exposes the internal
     liquid--liquid interface to GfsView. */
  scalar * gfs_fields = {f, f_lower, f1, f2, f3, cs,
#if FILTERED
                         sf, sf_lower,
#endif
                         u.x, u.y, p, gammaDot, muMix, rhov,
                         omega_amr, shear_amr, level_field};
  output_gfs_safe (name, gfs_fields, true);
}

#if ENABLE_DUMP_OUTPUT
event snapshots (t = 0.; t += DUMP_DT; t <= END_TIME + 1.e-12)
{
  char name[256];
  snprintf (name, sizeof(name), "dumps/dump-%010.5f", t);
  dump (file = name);
}
#endif

event stop (t = END_TIME)
{
  if (wave_amplitude_fp) {
    fclose (wave_amplitude_fp);
    wave_amplitude_fp = NULL;
  }
  profile1d_free (&profile_u0);
  profile1d_free (&profile_p);
  return 1;
}
