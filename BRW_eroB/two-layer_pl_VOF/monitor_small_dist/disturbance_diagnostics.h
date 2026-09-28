/** Read-only, conservative AMR column diagnostics. 2026-09-25.
    Define DISTURBANCE_CORE_ONLY to test the ordinary C quadrature core.
    In the solver include after the case macros and VOF field definitions. */
#ifndef DISTURBANCE_DIAGNOSTICS_H
#define DISTURBANCE_DIAGNOSTICS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define DG_PI 3.14159265358979323846264338327950288

typedef struct {
  double mean, rms, linf, minimum, maximum, creal, cimag, amplitude, phase;
  double harmonic2_amplitude, harmonic3_amplitude, higher_rms;
} DisturbanceMeasure;

/* Conservatively add a constant-cell VOF height contribution to uniform bins.
   Domain and cell abscissae have origin zero in this supplied solver. */
static void disturbance_add_cell (double * h, int n, double length,
                                  double center, double delta,
                                  double fraction, double solid_fraction)
{
  const double dx = length/n;
  const double left = fmax(0., center - delta/2.);
  const double right = fmin(length, center + delta/2.);
  int first = (int) floor(left/dx);
  int last = (int) ceil(right/dx) - 1;
  if (first < 0) first = 0;
  if (last >= n) last = n - 1;
  for (int j = first; j <= last; j++) {
    const double overlap = fmax(0., fmin(right, (j+1)*dx) - fmax(left, j*dx));
    h[j] += fraction*solid_fraction*delta*overlap/dx;
  }
}

/* Discrete Fourier coefficient on conservative column-average values at bin
   centres, matching the common PyPDE diagnostic. This converges to the exact
   continuous integral as the horizontal projection is refined. */
static DisturbanceMeasure disturbance_measure (const double * h, int n, int mode)
{
  DisturbanceMeasure m = {0};
  m.minimum = DBL_MAX; m.maximum = -DBL_MAX;
  for (int j = 0; j < n; j++) {
    m.mean += h[j]/n;
    if (h[j] < m.minimum) m.minimum = h[j];
    if (h[j] > m.maximum) m.maximum = h[j];
  }
  double variance = 0., re[3] = {0}, im[3] = {0};
  for (int j = 0; j < n; j++) {
    const double d = h[j] - m.mean;
    variance += d*d/n;
    if (fabs(d) > m.linf) m.linf = fabs(d);
    for (int a = 0; a < 3; a++) {
      const double angle = 2.*DG_PI*(a+1)*mode*(j+0.5)/n;
      re[a] += d*cos(angle)/n;
      im[a] -= d*sin(angle)/n;
    }
  }
  m.creal = re[0]; m.cimag = im[0];
  m.amplitude = 2.*hypot(re[0], im[0]);
  m.phase = atan2(im[0], re[0]);
  m.rms = sqrt(fmax(0., variance));
  m.harmonic2_amplitude = 2.*hypot(re[1], im[1]);
  m.harmonic3_amplitude = 2.*hypot(re[2], im[2]);
  m.higher_rms = sqrt(fmax(0., variance - m.amplitude*m.amplitude/2.));
  return m;
}

#ifndef DISTURBANCE_CORE_ONLY
#ifndef CASE_DISTURBANCE_ENABLED
#define CASE_DISTURBANCE_ENABLED 0
#define CASE_DISTURBANCE_TIME_INTERVAL 0.
#define CASE_DISTURBANCE_EVERY_STEPS 50
#define CASE_DISTURBANCE_COLUMNS 1024
#endif

static FILE * disturbance_fp = NULL;
static double disturbance_next_time = 0., disturbance_last_time = -1.;
static double disturbance_initial_mean[2] = {0., 0.};

static void disturbance_metadata (int columns)
{
  if (pid() != 0) return;
  FILE * fp = fopen("disturbance_metadata.json", "w");
  if (!fp) { perror("disturbance_metadata.json"); exit(1); }
  fprintf(fp,
    "{\n  \"schema\": \"roll-wave-disturbance-metadata\",\n  \"schema_version\": 2,\n"
    "  \"case_groups\": {\"Fr_l\": %.17g, \"S0\": %.17g, \"n_l\": %.17g, \"n_u\": %.17g, \"rho_r\": %.17g, \"h_r\": %.17g, \"R_eta_I\": %.17g},\n"
    "  \"geometry\": {\"k_star\": %.17g, \"mode_number\": %d, \"domain_length_star\": %.17g},\n"
    "  \"boundary_conditions\": {\"outer_top\": \"u.n=dirichlet(0);u.t=neumann(0);uf.n=0 at y=LX\", \"impermeable_slip_equivalence_validated\": null, \"legacy_flag_status\": \"deprecated: metadata records the implemented boundary conditions, not a validation of model equivalence\"},\n"
    "  \"normalization\": {\"depth\": \"H_l\", \"time\": \"T=S0*t_star\", \"t_star\": \"t*Ubar_l/H_l\", \"amplitude\": \"eta/H_l\"},\n"
    "  \"observation\": {\"projection\": \"conservative column integral of cumulative VOF\", \"fourier\": \"discrete Fourier coefficient of conservative uniform columns\", \"ncolumns\": %d, \"time_interval_star\": %.17g, \"every_steps\": %d, \"sampling\": \"actual existing solver step; time cadence takes precedence\"},\n"
    "  \"simulation\": {\"model\": \"Basilisk three-phase VOF\", \"air_density_ratio\": %.17g, \"air_eta\": %.17g, \"sigma_internal_star\": %.17g, \"sigma_free_star\": %.17g, \"lower_epsilon\": %.17g, \"upper_epsilon\": %.17g, \"ceiling_star\": %.17g, \"lower_eta_max\": %.17g, \"upper_eta_max\": %.17g},\n",
    (double)CASE_FROUDE, (double)SLOPE_TAN, (double)N_LOWER, (double)N_UPPER,
    (double)RHOUPPERLAYER, (double)H2, (double)CASE_INTERFACIAL_APPARENT_VISCOSITY_RATIO,
    2.*DG_PI*WAVE_MODE_NUMBER/LX, WAVE_MODE_NUMBER, (double)LX,
    columns, (double)CASE_DISTURBANCE_TIME_INTERVAL, CASE_DISTURBANCE_EVERY_STEPS,
    (double)AIRRHO, (double)AIRMU, (double)SIGMA_INTERNAL, (double)SIGMA_FREE,
    (double)CASE_LOWER_EPSILON, (double)CASE_UPPER_EPSILON, (double)LX,
    (double)MUMAXLOWER, (double)MUMAXUPPER);
  fprintf(fp,
    "  \"vertical_geometry\": {\"actual_top_height_star\": %.17g, \"nominal_near_air_matching_height_star\": %.17g, \"aligned_near_air_matching_height_star\": %.17g, \"near_air_matching_plane_is_boundary\": false, \"ceiling_star_legacy_alias\": \"simulation.ceiling_star is the actual top height, not the near-air matching plane\"},\n"
    "  \"rheology_interpretation\": {\"R_eta_I_definition\": \"eta_lower_I / eta_upper_I in the compatible steady base state\", \"R_eta_I_role\": \"base-state input used once to derive fixed constitutive coefficients; not an instantaneous ratio constraint\", \"lambda_lower\": %.17g, \"lambda_upper\": %.17g, \"derived_consistency_ratio\": %.17g},\n"
    "  \"observable_definitions\": {\"interface_height\": \"eta_I(x)=integral f_lower dy\", \"surface_height\": \"eta_S(x)=integral f dy\", \"height_units\": \"H_l\", \"column_centres\": \"x_j=(j+1/2)*L/N\", \"coefficient\": \"C_m=(1/N)*sum_j[(eta_j-mean(eta))*exp(-i*2*pi*m*x_j/L)]\", \"amplitude\": \"A_m=2*abs(C_m); no finite-bin sinc correction\", \"phase\": \"arg(C_m) in radians\", \"rms\": \"sqrt(mean((eta-mean(eta))^2))\", \"ln_amplitude\": \"natural logarithm ln(A_m); NaN if A_m=0\", \"initial_amplitude_source\": \"first CSV row measures the projected VOF; analytical IC amplitudes below describe the installed continuum shape\"},\n"
    "  \"growth_rate_definitions\": {\"csv_time\": \"T=S0*t_star\", \"csv_time_star\": \"t_star=t_dimensional*Ubar_l/H_l\", \"natural_log_amplitude_slope\": \"sigma_T=d ln(A_m)/dT\", \"log10_amplitude_slope\": \"sigma_T/ln(10)\", \"quadratic_disturbance_slope\": \"d ln(A_m^2)/dT=2*sigma_T; this is not a definition of total kinetic energy\", \"solver_time_growth_rate\": \"sigma_star=S0*sigma_T\"},\n"
    "  \"initial_condition\": {\n"
    "    \"declared\": {\"perturbation_mode\": %d, \"tracked_mode\": %d, \"phase_radians\": %.17g, \"lower_relative_thickness_amplitude\": %.17g, \"upper_relative_thickness_amplitude\": %.17g, \"velocity_amplitude_relative_to_base_surface_velocity\": %.17g, \"base_surface_velocity_star\": %.17g},\n",
    (double)LX, (double)CASE_NOMINAL_CEILING, near_air_y,
    (double)CASE_LAMBDA_LOWER, (double)CASE_LAMBDA_UPPER,
    (double)CASE_DERIVED_CONSISTENCY_RATIO,
    PERTURBATION_MODE_NUMBER, WAVE_MODE_NUMBER, (double)PERTURBATION_PHASE,
    (double)LOWER_DEPTH_PERTURBATION_AMPLITUDE,
    (double)UPPER_DEPTH_PERTURBATION_AMPLITUDE,
    (double)VELOCITY_PERTURBATION_AMPLITUDE, u_reference);
#if CASE_LINEAR_EIGENMODE
  fprintf(fp,
    "    \"type\": \"full_ns_eigenmode\",\n"
    "    \"declared_controls_usage\": \"modes must match the installed eigenmode; declared phase, layer amplitudes and velocity amplitude are overridden\",\n"
    "    \"shape\": \"eta_a=base_height_a+scale*Re(eta_hat_a*exp(i*k_star*x_star))\",\n"
    "    \"velocity\": \"loaded full-NS Eulerian velocity eigenfunctions plus compatible base velocity, evaluated with the first-order phase-mapping correction\",\n"
    "    \"installed_eigenmode\": {\"k_star\": %.17g, \"scale\": %.17g, \"interface_eta_hat_real\": %.17g, \"interface_eta_hat_imag\": %.17g, \"surface_eta_hat_real\": %.17g, \"surface_eta_hat_imag\": %.17g, \"interface_amplitude_star\": %.17g, \"surface_amplitude_star\": %.17g, \"upper_thickness_amplitude_star\": %.17g}\n",
    linear_mode_k, linear_mode_amplitude,
    linear_mode_eta_real[0], linear_mode_eta_imag[0],
    linear_mode_eta_real[1], linear_mode_eta_imag[1],
    linear_mode_amplitude*hypot(linear_mode_eta_real[0],linear_mode_eta_imag[0]),
    linear_mode_amplitude*hypot(linear_mode_eta_real[1],linear_mode_eta_imag[1]),
    linear_mode_amplitude*hypot(linear_mode_eta_real[1]-linear_mode_eta_real[0],
                                linear_mode_eta_imag[1]-linear_mode_eta_imag[0]));
#elif CASE_FRONT_RUNNER_ENABLED
  fprintf(fp,
    "    \"type\": \"front_runner\",\n"
    "    \"declared_controls_usage\": \"relative layer amplitudes multiply the localized hump; tracked mode selects only the diagnostic Fourier component; declared periodic phase, perturbation mode and velocity amplitude do not control this initializer\",\n"
    "    \"shape\": \"h_a=H_a*(1+relative_amplitude_a*bump); bump=cos(2*pi*(x-center)/wavelength) for abs(x-center)<wavelength/4, otherwise zero\",\n"
    "    \"velocity\": \"localized constant-Froude streamfunction initialization with near-air matching and far-air plug extension\",\n"
    "    \"front_runner\": {\"wavelength_star\": %.17g, \"center_star\": %.17g, \"localized_shape_is_single_fourier_mode\": false}\n",
    (double)CASE_FR_WAVELENGTH, (double)CASE_FR_CENTER);
#else
  fprintf(fp,
    "    \"type\": \"geometric_periodic\",\n"
    "    \"declared_controls_usage\": \"layer depths share sin(2*pi*perturbation_mode*x/L+phase); tracked_mode independently selects the diagnostic component\",\n"
    "    \"shape\": \"h_a=H_a*(1+relative_amplitude_a*sin(2*pi*perturbation_mode*x/L+phase))\",\n"
    "    \"velocity\": \"compatible flat-layer base profile plus optional divergence-free sinusoidal streamfunction perturbation with declared velocity amplitude and phase; perturbation vanishes above the near-air matching plane\",\n"
    "    \"analytic_shape_amplitudes\": {\"mode\": %d, \"interface_amplitude_star\": %.17g, \"surface_amplitude_star\": %.17g, \"upper_thickness_amplitude_star\": %.17g}\n",
    PERTURBATION_MODE_NUMBER,
    fabs(H1*LOWER_DEPTH_PERTURBATION_AMPLITUDE),
    fabs(H1*LOWER_DEPTH_PERTURBATION_AMPLITUDE+H2*UPPER_DEPTH_PERTURBATION_AMPLITUDE),
    fabs(H2*UPPER_DEPTH_PERTURBATION_AMPLITUDE));
#endif
  fprintf(fp, "  }\n}\n");
  fclose(fp);
}

static void disturbance_write_measure (FILE * fp, DisturbanceMeasure m, double mean0)
{
  fprintf(fp, ",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",
    m.mean, m.mean-mean0, m.rms, m.maximum-m.minimum, m.creal, m.cimag,
    m.amplitude, m.phase, m.harmonic2_amplitude, m.higher_rms,
    m.minimum, m.maximum, m.harmonic3_amplitude, m.amplitude > 0. ? log(m.amplitude) : NAN, m.linf);
}

static void disturbance_maybe_record (double current_time, int step, int force)
{
  if (!CASE_DISTURBANCE_ENABLED) return;
  const double tol = 32.*DBL_EPSILON*fmax(1., fabs(current_time));
  if (disturbance_last_time >= 0. && fabs(current_time-disturbance_last_time) <= tol) return;
  if (!force && disturbance_last_time >= 0.) {
    if (CASE_DISTURBANCE_TIME_INTERVAL > 0.) {
      if (current_time + tol < disturbance_next_time) return;
    }
    else if (step % CASE_DISTURBANCE_EVERY_STEPS) return;
  }
  const int columns = CASE_DISTURBANCE_COLUMNS;
  double * column = (double *) calloc((size_t)2*columns, sizeof(double));
  if (!column) { perror("disturbance columns"); exit(1); }
  /* Serial within each rank avoids OpenMP array write races. foreach iterates
     owned leaf cells; MPI sums disjoint rank contributions after projection. */
  foreach(serial) {
    if (cm[] > 0.) {
      disturbance_add_cell(column, columns, LX, x, Delta, f[], cm[]);
      disturbance_add_cell(column + columns, columns, LX, x, Delta, f_lower[], cm[]);
    }
  }
#if _MPI
  MPI_Allreduce(MPI_IN_PLACE, column, 2*columns, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (pid() == 0) {
    DisturbanceMeasure surface = disturbance_measure(column, columns, WAVE_MODE_NUMBER);
    DisturbanceMeasure inter = disturbance_measure(column + columns, columns, WAVE_MODE_NUMBER);
    if (!disturbance_fp) {
      disturbance_fp = fopen("disturbance_history.csv", "w");
      if (!disturbance_fp) { perror("disturbance_history.csv"); exit(1); }
      disturbance_initial_mean[0] = surface.mean;
      disturbance_initial_mean[1] = inter.mean;
      fprintf(disturbance_fp, "time,time_star,mode,k_star,step,ncolumns");
      const char * names[2] = {"surface", "interface"};
      for (int a = 0; a < 2; a++) {
        const char * prefix = names[a];
        fprintf(disturbance_fp, ",%s_mean,%s_mean_drift,%s_rms,%s_range,%s_creal,%s_cimag,%s_amplitude,%s_phase,%s_harmonic2_amplitude,%s_higher_rms,%s_min,%s_max,%s_harmonic3_amplitude,%s_ln_amplitude,%s_linf",
          prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix);
      }
      fprintf(disturbance_fp, "\n");
      disturbance_metadata(columns);
    }
    fprintf(disturbance_fp, "%.17g,%.17g,%d,%.17g,%d,%d",
      SLOPE_TAN*current_time, current_time, WAVE_MODE_NUMBER,
      2.*DG_PI*WAVE_MODE_NUMBER/LX, step, columns);
    disturbance_write_measure(disturbance_fp, surface, disturbance_initial_mean[0]);
    disturbance_write_measure(disturbance_fp, inter, disturbance_initial_mean[1]);
    fprintf(disturbance_fp, "\n"); fflush(disturbance_fp);
  }
  free(column);
  disturbance_last_time = current_time;
  if (CASE_DISTURBANCE_TIME_INTERVAL > 0.)
    disturbance_next_time = (floor((current_time+tol)/CASE_DISTURBANCE_TIME_INTERVAL)+1.)*CASE_DISTURBANCE_TIME_INTERVAL;
}

static void disturbance_close (void)
{
  if (disturbance_fp) { fclose(disturbance_fp); disturbance_fp = NULL; }
}
#endif
#endif
