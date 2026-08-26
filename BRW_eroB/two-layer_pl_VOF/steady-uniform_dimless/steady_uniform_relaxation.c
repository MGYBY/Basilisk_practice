/**
# Short-domain steady-uniform relaxation benchmark

This case uses the uploaded dimensionless three-phase two-layer power-law
solver without a roll-wave disturbance.  The two interfaces are flat and the
complete fluid domain starts from a small uniform streamwise velocity.  Gravity,
no slip at the bed and the selected generalized-Newtonian rheologies then drive
the solution toward the one-dimensional steady-uniform reference profile
written by `generate_steady_uniform_case.py`.

The physical domain is periodic in x and deliberately short.  A single uniform
quadtree level is used so the comparison is not contaminated by refinement or
coarsening.  Text profiles and integral error metrics are written at exact
intervals in dimensionless time.
*/

#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "three-phase-stratified-powerlaw-dimensionless.h"
#include "conserving3f-stratified.h"
#include "profile1d.h"
#include "navier-stokes/perfs.h"
#include "base_state/generated_case.h"
#include "base_state/steady_uniform_generated.h"

#define H1 1.0
#define H2 CASE_DEPTH_RATIO
#define LIQUID_DEPTH (H1 + H2)
#define LX CASE_LX
#define CEILING_HEIGHT CASE_ALIGNED_CEILING
#define LEVEL CASE_MAXLEVEL
#define END_TIME CASE_END_TIME
#define OUTPUT_DT CASE_COMPARISON_OUTPUT_DT
#define INITIAL_U CASE_INITIAL_UNIFORM_VELOCITY

#define G_ACCEL (sqrt(1. + sq(CASE_SLOPE_TAN))/sq(CASE_FROUDE))
#define SIN_THETA (CASE_SLOPE_TAN/sqrt(1. + sq(CASE_SLOPE_TAN)))
#define COS_THETA (1./sqrt(1. + sq(CASE_SLOPE_TAN)))

face vector gravity_accel[];

Profile1D exact_u, exact_p, exact_gamma, exact_eta, exact_tau;
scalar level_field[];

static FILE * history_fp = NULL;
static int profile_index = 0;
static int convergence_count = 0;

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
uf.n[bottom] = 0.;
uf.n[top] = 0.;

/* Preserve the uploaded solver's embedded free-slip ceiling treatment. */
u.n[embed] = neumann(0.);
u.t[embed] = neumann(0.);

static inline double grid_delta (void)
{
  return LX/(double)(1 << LEVEL);
}

static void read_profile (Profile1D * p, const char * basename)
{
  char path[1024];
  snprintf (path, sizeof(path), "base_state/%s", basename);
  profile1d_read (p, path);
  fprintf (ferr, "# read %s: n=%d range=[%g,%g]\n",
           path, p->n, p->x[0], p->x[p->n - 1]);
}

static void validate_profiles (void)
{
  const Profile1D * profiles[] = {
    &exact_u, &exact_p, &exact_gamma, &exact_eta, &exact_tau
  };
  const char * names[] = {
    "velocity", "pressure", "gamma", "viscosity", "traction"
  };
  const double tol = 1.e-10*max(1., CEILING_HEIGHT);

  for (int k = 0; k < 5; k++) {
    const Profile1D * p = profiles[k];
    if (p->n < 2 || fabs(p->x[0]) > tol ||
        fabs(p->x[p->n - 1] - CEILING_HEIGHT) > tol) {
      fprintf (ferr,
               "Invalid %s reference profile: n=%d range=[%g,%g], "
               "expected [0,%g]. Regenerate base_state.\n",
               names[k], p->n, p->x[0], p->x[p->n - 1],
               (double) CEILING_HEIGHT);
      exit (1);
    }
  }

  if (fabs(profile1d_eval(&exact_u, 0.)) > 1.e-10 ||
      fabs(profile1d_eval(&exact_p, CEILING_HEIGHT)) > 1.e-10) {
    fprintf (ferr,
             "Reference profile does not satisfy zero bed velocity and "
             "zero ceiling pressure.\n");
    exit (1);
  }
}

static inline void build_ceiling (void)
{
  solid (cs, fs, CEILING_HEIGHT - y);
  fractions_cleanup (cs, fs);
}

typedef struct {
  double v1, v2, v3;
  double mean_u1, mean_u2, mean_u3;
  double exact_mean_u1, exact_mean_u2, exact_mean_u3;
  double rms_u1, rms_u2, rms_u3;
  double rel_rms_u1, rel_rms_u2, rel_rms_u3;
  double linf_u1, linf_u2, linf_u3;
  double p_offset, rms_p, linf_p;
  double rms_v, max_v;
  double u_interface, u_surface;
} ComparisonMetrics;

static ComparisonMetrics comparison_metrics (void)
{
  stratified_update_phase_fields();
  boundary ((scalar *){u, p, gammaDot, muMix, rhov});

  double pressure_difference = 0., fluid_volume = 0.;
  foreach (reduction(+:pressure_difference) reduction(+:fluid_volume))
    if (cs[] > 0.) {
      const double vol = dv();
      pressure_difference +=
        (p[] - profile1d_eval(&exact_p, y))*vol;
      fluid_volume += vol;
    }
  const double p_offset = pressure_difference/(fluid_volume + SEPS);

  double v1 = 0., v2 = 0., v3 = 0.;
  double u1 = 0., u2 = 0., u3 = 0.;
  double ue1 = 0., ue2 = 0., ue3 = 0.;
  double ue21 = 0., ue22 = 0., ue23 = 0.;
  double e21 = 0., e22 = 0., e23 = 0.;
  double einf1 = 0., einf2 = 0., einf3 = 0.;
  double p2 = 0., pinf = 0.;
  double vv2 = 0., vmax = 0.;

  foreach (reduction(+:v1) reduction(+:v2) reduction(+:v3)
           reduction(+:u1) reduction(+:u2) reduction(+:u3)
           reduction(+:ue1) reduction(+:ue2) reduction(+:ue3)
           reduction(+:ue21) reduction(+:ue22) reduction(+:ue23)
           reduction(+:e21) reduction(+:e22) reduction(+:e23)
           reduction(max:einf1) reduction(max:einf2) reduction(max:einf3)
           reduction(+:p2) reduction(max:pinf)
           reduction(+:vv2) reduction(max:vmax))
    if (cs[] > 0.) {
      double a1, a2, a3;
      stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
      const double vol = dv();
      const double un = u.x[];
      const double ve = u.y[];
      const double ux = profile1d_eval(&exact_u, y);
      const double du = un - ux;
      const double dp = p[] - p_offset - profile1d_eval(&exact_p, y);

      v1 += a1*vol; v2 += a2*vol; v3 += a3*vol;
      u1 += a1*un*vol; u2 += a2*un*vol; u3 += a3*un*vol;
      ue1 += a1*ux*vol; ue2 += a2*ux*vol; ue3 += a3*ux*vol;
      ue21 += a1*sq(ux)*vol; ue22 += a2*sq(ux)*vol;
      ue23 += a3*sq(ux)*vol;
      e21 += a1*sq(du)*vol; e22 += a2*sq(du)*vol;
      e23 += a3*sq(du)*vol;
      einf1 = max(einf1, a1*fabs(du));
      einf2 = max(einf2, a2*fabs(du));
      einf3 = max(einf3, a3*fabs(du));
      p2 += sq(dp)*vol;
      pinf = max(pinf, fabs(dp));
      vv2 += sq(ve)*vol;
      vmax = max(vmax, fabs(ve));
    }

  ComparisonMetrics m = {0};
  m.v1 = v1; m.v2 = v2; m.v3 = v3;
  m.mean_u1 = u1/(v1 + SEPS);
  m.mean_u2 = u2/(v2 + SEPS);
  m.mean_u3 = u3/(v3 + SEPS);
  m.exact_mean_u1 = ue1/(v1 + SEPS);
  m.exact_mean_u2 = ue2/(v2 + SEPS);
  m.exact_mean_u3 = ue3/(v3 + SEPS);
  m.rms_u1 = sqrt(e21/(v1 + SEPS));
  m.rms_u2 = sqrt(e22/(v2 + SEPS));
  m.rms_u3 = sqrt(e23/(v3 + SEPS));
  m.rel_rms_u1 = m.rms_u1/(sqrt(ue21/(v1 + SEPS)) + SEPS);
  m.rel_rms_u2 = m.rms_u2/(sqrt(ue22/(v2 + SEPS)) + SEPS);
  m.rel_rms_u3 = m.rms_u3/(sqrt(ue23/(v3 + SEPS)) + SEPS);
  m.linf_u1 = einf1; m.linf_u2 = einf2; m.linf_u3 = einf3;
  m.p_offset = p_offset;
  m.rms_p = sqrt(p2/(fluid_volume + SEPS));
  m.linf_p = pinf;
  m.rms_v = sqrt(vv2/(fluid_volume + SEPS));
  m.max_v = vmax;
  m.u_interface = interpolate(u.x, LX/2., H1);
  m.u_surface = interpolate(u.x, LX/2., LIQUID_DEPTH);
  return m;
}

static void write_profile (const ComparisonMetrics * m)
{
  char name[256];
  snprintf (name, sizeof(name),
            "profiles/profile_%04d_t%010.4f.tsv", profile_index, t);
  FILE * fp = fopen (name, "w");
  if (!fp) {
    perror(name);
    exit (1);
  }

  fprintf (fp,
           "z_star\tfLower\tfLiquid\talphaLower\talphaUpper\talphaAir"
           "\tuNum_star\tuExact_star\tuError_star\tvNum_star"
           "\tpNum_star\tpNumGaugeCorrected_star\tpExact_star\tpError_star"
           "\tgammaNum_star\tgammaExact_star\tetaNum_star\tetaExact_star"
           "\ttauNumMagnitude_star\ttauExact_star\trhoNum_star\tlevel\tcs\n");

  const double ds = grid_delta();
  const double xx = LX/2.;
  for (double yy = ds/2.; yy < CEILING_HEIGHT - ds/4.; yy += ds) {
    const double lower = interpolate(f_lower, xx, yy);
    const double liquid = interpolate(f, xx, yy);
    double a1, a2, a3;
    stratified_phase_values (lower, liquid, &a1, &a2, &a3);
    const double un = interpolate(u.x, xx, yy);
    const double vn = interpolate(u.y, xx, yy);
    const double pn = interpolate(p, xx, yy);
    const double gn = interpolate(gammaDot, xx, yy);
    const double etan = interpolate(muMix, xx, yy);
    const double uex = profile1d_eval(&exact_u, yy);
    const double pex = profile1d_eval(&exact_p, yy);
    const double gex = profile1d_eval(&exact_gamma, yy);
    const double etaex = profile1d_eval(&exact_eta, yy);
    const double tauex = profile1d_eval(&exact_tau, yy);
    const double rhon = interpolate(rhov, xx, yy)/
      max(interpolate(cs, xx, yy), SEPS);
    const double lev = interpolate(level_field, xx, yy);
    const double csv = interpolate(cs, xx, yy);

    fprintf (fp,
             "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e"
             "\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             yy, lower, liquid, a1, a2, a3,
             un, uex, un - uex, vn,
             pn, pn - m->p_offset, pex, pn - m->p_offset - pex,
             gn, gex, etan, etaex, etan*gn, tauex,
             rhon, lev, csv);
  }
  fclose (fp);
}

static void write_history_row (const ComparisonMetrics * m)
{
  if (!history_fp) {
    history_fp = fopen ("history/comparison_history.tsv", "w");
    if (!history_fp) {
      perror("history/comparison_history.tsv");
      exit (1);
    }
    fprintf (history_fp,
             "outputIndex\ti\tt_star\tdt_star\tVlower\tVupper\tVair"
             "\tmeanUlower\texactMeanUlower\tmeanUupper\texactMeanUupper"
             "\tmeanUair\texactMeanUair\tUinterface\texactUinterface"
             "\tUsurface\texactUsurface"
             "\trmsUlower\trelRmsUlower\tlinfUlower"
             "\trmsUupper\trelRmsUupper\tlinfUupper"
             "\trmsUair\trelRmsUair\tlinfUair"
             "\tpressureOffset\trmsPressure\tlinfPressure"
             "\trmsVerticalVelocity\tmaxVerticalVelocity\tmgpIterations"
             "\tmguIterations\tmgpResidual\tmguResidual\n");
  }

  fprintf (history_fp,
           "%d\t%d\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e"
           "\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e"
           "\t%.16e\t%.16e\t%.16e\t%.16e"
           "\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e"
           "\t%.16e\t%.16e\t%.16e"
           "\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e"
           "\t%d\t%d\t%.16e\t%.16e\n",
           profile_index, iter, t, dt, m->v1, m->v2, m->v3,
           m->mean_u1, m->exact_mean_u1,
           m->mean_u2, m->exact_mean_u2,
           m->mean_u3, m->exact_mean_u3,
           m->u_interface, profile1d_eval(&exact_u, H1),
           m->u_surface, profile1d_eval(&exact_u, LIQUID_DEPTH),
           m->rms_u1, m->rel_rms_u1, m->linf_u1,
           m->rms_u2, m->rel_rms_u2, m->linf_u2,
           m->rms_u3, m->rel_rms_u3, m->linf_u3,
           m->p_offset, m->rms_p, m->linf_p,
           m->rms_v, m->max_v,
           mgp.i, mgu.i, mgp.resa, mgu.resa);
  fflush (history_fp);
}

int main (void)
{
  if (CASE_MINLEVEL != LEVEL || CASE_INITLEVEL != LEVEL) {
    fprintf (ferr,
             "This benchmark requires min_level=initial_level=max_level; "
             "got %d/%d/%d.\n",
             CASE_MINLEVEL, CASE_INITLEVEL, LEVEL);
    return 1;
  }
  if (CEILING_HEIGHT >= LX || LIQUID_DEPTH >= CEILING_HEIGHT) {
    fprintf (ferr,
             "Require liquid depth < embedded ceiling < square-domain size; "
             "got liquid=%g ceiling=%g Lx=%g.\n",
             (double)LIQUID_DEPTH, (double)CEILING_HEIGHT, (double)LX);
    return 1;
  }

  size (LX);
  origin (0., 0.);
  periodic (right);

  rho1 = 1.;
  rho2 = CASE_DENSITY_RATIO;
  rho3 = CASE_AIR_DENSITY_RATIO;

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
    .n = 1., .epsilon = 0., .tau0 = 0., .m = 0.,
    .eta_min = CASE_AIR_ETA,
    .eta_max = CASE_AIR_ETA,
    .gamma_floor = min(CASE_LOWER_GAMMA_FLOOR,
                       CASE_UPPER_GAMMA_FLOOR)
  };

  a = gravity_accel;
  CFL = CASE_CFL;
  DT = CASE_MAX_DT;
  TOLERANCE = CASE_SOLVER_TOLERANCE;
  NITERMAX = CASE_SOLVER_NITERMAX;

  init_grid (1 << LEVEL);
  run();
}

event init (i = 0)
{
  if (pid() == 0)
    system ("mkdir -p profiles history analysis verification");
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  read_profile (&exact_u, "velocity.dat");
  read_profile (&exact_p, "pressure.dat");
  read_profile (&exact_gamma, "gamma.dat");
  read_profile (&exact_eta, "viscosity.dat");
  read_profile (&exact_tau, "traction.dat");
  validate_profiles();

  build_ceiling();
  fraction (f_lower, H1 - y);
  fraction (f, LIQUID_DEPTH - y);

  foreach() {
    u.x[] = cs[] > 0. ? INITIAL_U : 0.;
    u.y[] = 0.;
    p[] = cs[] > 0. ? profile1d_eval(&exact_p, y) : 0.;
    level_field[] = level;
  }
  stratified_enforce_hierarchy();
  boundary ((scalar *){u, p, level_field});

  /* Populate rheology and diagnostic fields before the t=0 text output. */
  event ("properties");

  fprintf (ferr,
           "# steady-uniform relaxation: Lx=%g liquidDepth=%g ceiling=%g "
           "level=%d Delta=%g cells/Hl=%g\n",
           (double)LX, (double)LIQUID_DEPTH, (double)CEILING_HEIGHT,
           LEVEL, grid_delta(), H1/grid_delta());
  fprintf (ferr,
           "# initial uniform velocity=%g; target lower mean=1; "
           "output interval=%g; end time=%g\n",
           (double)INITIAL_U, (double)OUTPUT_DT, (double)END_TIME);
  fprintf (ferr,
           "# groups: Fr_l=%g S0=%g h_r=%g r_rho=%g n_l=%g n_u=%g "
           "Lambda_l=%g Lambda_u=%g eta_air=%g\n",
           (double)CASE_FROUDE, (double)CASE_SLOPE_TAN,
           (double)CASE_DEPTH_RATIO, (double)CASE_DENSITY_RATIO,
           (double)CASE_LOWER_N, (double)CASE_UPPER_N,
           (double)CASE_LAMBDA_LOWER, (double)CASE_LAMBDA_UPPER,
           (double)CASE_AIR_ETA);
}

event acceleration (i++)
{
  face vector av = a;
  const double gs = G_ACCEL*SIN_THETA;
  const double gn = G_ACCEL*COS_THETA;

  foreach_face(x) {
    const double lower = face_value(f_lower, 0);
    const double liquid = face_value(f, 0);
    double a1, a2, a3;
    stratified_phase_values (lower, liquid, &a1, &a2, &a3);
    const double rhof = a1*rho1 + a2*rho2 + a3*rho3;
#if CASE_AIR_STREAMWISE_GRAVITY
    av.x[] += gs;
#else
    av.x[] += gs*(a1*rho1 + a2*rho2)/(rhof + SEPS);
#endif
  }
  foreach_face(y)
    av.y[] -= gn;
}

event text_outputs (t = 0.; t += OUTPUT_DT; t <= END_TIME + 1.e-12)
{
  foreach()
    level_field[] = level;
  boundary ({level_field});

  const ComparisonMetrics m = comparison_metrics();
  write_history_row (&m);
  write_profile (&m);

  const double worst_liquid_rel = max(m.rel_rms_u1, m.rel_rms_u2);
  fprintf (ferr,
           "# compare %d t=%g: meanU=(%g,%g,%g) relRMS=(%.3e,%.3e,%.3e) "
           "pRMS=%.3e vMax=%.3e\n",
           profile_index, t, m.mean_u1, m.mean_u2, m.mean_u3,
           m.rel_rms_u1, m.rel_rms_u2, m.rel_rms_u3,
           m.rms_p, m.max_v);

  if (worst_liquid_rel < CASE_CONVERGENCE_RELATIVE_RMS)
    convergence_count++;
  else
    convergence_count = 0;

  profile_index++;

#if CASE_STOP_WHEN_CONVERGED
  if (convergence_count >= CASE_CONVERGENCE_HOLD_OUTPUTS) {
    fprintf (ferr,
             "# convergence criterion satisfied for %d consecutive outputs; "
             "stopping at t=%g.\n",
             convergence_count, t);
    return 1;
  }
#endif
}

event iteration_diagnostics (i += 50)
{
  fprintf (ferr,
           "# iter %d t=%g dt=%g mgp=%d/%g mgu=%d/%g cells=%ld\n",
           i, t, dt, mgp.i, mgp.resa, mgu.i, mgu.resa, grid->tn);
}

event stop (t = END_TIME)
{
  if (history_fp) {
    fclose (history_fp);
    history_fp = NULL;
  }
  profile1d_free (&exact_u);
  profile1d_free (&exact_p);
  profile1d_free (&exact_gamma);
  profile1d_free (&exact_eta);
  profile1d_free (&exact_tau);
  return 1;
}
