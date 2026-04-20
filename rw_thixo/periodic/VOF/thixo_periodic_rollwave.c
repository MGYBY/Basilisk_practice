/**
# Periodic roll-wave case for the thixotropic two-phase VOF solver

This case implements the simple periodic roll-wave configuration described
in the uploaded derivation notes.

Main ideas of the setup:

- The domain is periodic in the streamwise direction `x`.
- The lower phase is a thixotropic liquid, the upper phase is air.
- The liquid is initialised from a precomputed steady-uniform base-state
  profile stored in `base_state_profile.txt`.
- A sinusoidal disturbance is imposed on both the free-surface position and
  the liquid velocity profile.
- The mapped coordinate `zeta = y/h(x)` is used so that the base-state
  profile on `0 <= zeta <= 1` can be reused under a disturbed free surface.
- Following the modelling note in the uploaded draft, the streamwise gravity
  component is applied only to the liquid phase, whereas the vertical gravity
  component is applied to both phases.

The base-state file is expected to contain three whitespace-separated columns:

```text
zeta   U(zeta)   Lambda(zeta)
```
*/

#include "grid/quadtree.h"
#include "adapt_wavelet_leave_interface_limited.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-thixo.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "navier-stokes/perfs.h"

/* --------------------------------------------------------------------------
   Compile-time parameters (all dimensionless).

   These defaults are only a starting point for a first periodic test. They can
   be overridden with `-D...` flags in the compile command used by runCase.sh.
   -------------------------------------------------------------------------- */
#ifndef MAXLEVEL
#define MAXLEVEL 11
#endif
#ifndef MINLEVEL
#define MINLEVEL 4
#endif
#ifndef INITLEVEL
#define INITLEVEL 8
#endif

#ifndef THIXO_T
#define THIXO_T 1.0
#endif
#ifndef THIXO_GAMMA
#define THIXO_GAMMA 5.0
#endif
#ifndef THIXO_KAPPA
#define THIXO_KAPPA 2e-6
#endif
#ifndef THIXO_A
#define THIXO_A 0.2
#endif
#ifndef THIXO_FRV
#define THIXO_FRV 0.75
#endif
#ifndef THIXO_SO
#define THIXO_SO 0.06
#endif
#ifndef THIXO_RHOR
#define THIXO_RHOR 0.01
#endif
#ifndef THIXO_MUR
#define THIXO_MUR 0.02
#endif
#ifndef THIXO_SIGMA
#define THIXO_SIGMA 0.0
#endif

/* Disturbance amplitudes and geometry. */
#ifndef DIST_AMP_H
#define DIST_AMP_H 1e-3
#endif
#ifndef DIST_AMP_U
#define DIST_AMP_U 1e-3
#endif
#ifndef WAVE_LENGTH
#define WAVE_LENGTH (2.0/THIXO_SO)
#endif
#ifndef TOP_HEIGHT
#define TOP_HEIGHT 4.0
#endif
#ifndef AIR_DECAY
#define AIR_DECAY 0.75
#endif
#ifndef OUTPUT_DT
#define OUTPUT_DT 0.25
#endif
#ifndef MAXTIME
#define MAXTIME 80.0
#endif
#ifndef PROFILE_FILE
#define PROFILE_FILE "base_state_profile.txt"
#endif

/* Adaptation/output tolerances. */
#ifndef F_ERR_OUT
#define F_ERR_OUT 1e-7
#endif
#ifndef OMEGA_ERR
#define OMEGA_ERR 0.15
#endif
#ifndef LAMBDA_ERR
#define LAMBDA_ERR 2e-3
#endif
#ifndef UFLUID_ERR
#define UFLUID_ERR 2e-2
#endif
#ifndef UAIR_ERR
#define UAIR_ERR 3e-2
#endif
#ifndef INTERFACE_PICK_TOL
#define INTERFACE_PICK_TOL 1e-10
#endif

/* --------------------------------------------------------------------------
   Storage for the base-state profile read from the external text file.
   -------------------------------------------------------------------------- */
static double * zprof = NULL;
static double * uprof = NULL;
static double * lprof = NULL;
static int nprof = 0;

/* Diagnostic/output fields. */
scalar omega_field[];
scalar velFluidNorm[], velAirNorm[];
scalar xpos_interface[], ypos_interface[];

/* --------------------------------------------------------------------------
   Boundary conditions.

   The actual computational fluid region is bounded above by an embedded wall at
   `y = TOP_HEIGHT`. The domain top boundary is therefore far away and mostly
   irrelevant, but we still assign standard consistent values.
   -------------------------------------------------------------------------- */
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top]    = dirichlet(0.);
u.t[top]    = neumann(0.);

lambda_thixo[bottom] = neumann(0.);
lambda_thixo[top]    = neumann(0.);
flambda[bottom]      = neumann(0.);
flambda[top]         = neumann(0.);

/* --------------------------------------------------------------------------
   Base-state profile utilities.
   -------------------------------------------------------------------------- */

/* Read a plain-text base-state profile.

   Supported format:
   - optional blank lines,
   - optional comment lines starting with '#',
   - data lines containing: zeta U lambda
*/
static void read_base_profile (const char * fname)
{
  FILE * fp = fopen (fname, "r");
  if (!fp) {
    fprintf (ferr, "Cannot open base-state profile file '%s'.\n", fname);
    exit (1);
  }

  int cap = 256;
  zprof = (double *) malloc (cap*sizeof(double));
  uprof = (double *) malloc (cap*sizeof(double));
  lprof = (double *) malloc (cap*sizeof(double));
  if (!zprof || !uprof || !lprof) {
    fprintf (ferr, "Memory allocation failure while reading '%s'.\n", fname);
    exit (1);
  }

  char line[512];
  while (fgets (line, sizeof(line), fp)) {
    double zeta, uu, lam;
    char * c = line;

    /* Skip leading spaces/tabs. */
    while (*c == ' ' || *c == '\t')
      c++;

    /* Ignore blank lines and comments. */
    if (*c == '\0' || *c == '\n' || *c == '#')
      continue;

    if (sscanf (c, "%lf %lf %lf", &zeta, &uu, &lam) != 3) {
      fprintf (ferr, "Malformed line in base-state profile '%s': %s", fname, line);
      exit (1);
    }

    if (nprof >= cap) {
      cap *= 2;
      zprof = (double *) realloc (zprof, cap*sizeof(double));
      uprof = (double *) realloc (uprof, cap*sizeof(double));
      lprof = (double *) realloc (lprof, cap*sizeof(double));
      if (!zprof || !uprof || !lprof) {
        fprintf (ferr, "Memory reallocation failure while reading '%s'.\n", fname);
        exit (1);
      }
    }

    if (nprof > 0 && zeta <= zprof[nprof - 1]) {
      fprintf (ferr,
               "The zeta coordinate in '%s' must be strictly increasing.\n",
               fname);
      exit (1);
    }

    zprof[nprof] = zeta;
    uprof[nprof] = uu;
    lprof[nprof] = lam;
    nprof++;
  }
  fclose (fp);

  if (nprof < 2) {
    fprintf (ferr, "Base-state profile '%s' must contain at least two data rows.\n", fname);
    exit (1);
  }
}

/* Simple monotone piecewise-linear interpolation on the profile. */
static inline double interp_profile (double zeta, const double * y)
{
  if (zeta <= zprof[0])
    return y[0];
  if (zeta >= zprof[nprof - 1])
    return y[nprof - 1];

  int lo = 0, hi = nprof - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi)/2;
    if (zprof[mid] <= zeta)
      lo = mid;
    else
      hi = mid;
  }

  double dz = zprof[hi] - zprof[lo];
  double t = (zeta - zprof[lo])/(dz + 1e-30);
  return y[lo] + t*(y[hi] - y[lo]);
}

static inline double base_u (double zeta)
{
  return interp_profile (zeta, uprof);
}

static inline double base_lambda (double zeta)
{
  return interp_profile (zeta, lprof);
}

/* --------------------------------------------------------------------------
   Initial condition helpers.
   -------------------------------------------------------------------------- */

/* Disturbed free surface. */
static inline double disturbed_depth (double x)
{
  return 1. + DIST_AMP_H*sin (2.*pi*x/WAVE_LENGTH);
}

/* Liquid velocity initial condition obtained from the mapped coordinate zeta=y/h(x). */
static inline double fluid_velocity_ic (double x, double y)
{
  double h = disturbed_depth (x);
  if (y > h)
    return 0.;

  double zeta = clamp (y/(h + 1e-30), 0., 1.);
  return base_u (zeta)*(1. + DIST_AMP_U*sin (2.*pi*x/WAVE_LENGTH));
}

/* Structure field initial condition using the same mapped coordinate. */
static inline double fluid_lambda_ic (double x, double y)
{
  double h = disturbed_depth (x);
  if (y > h)
    return thixo_lambda_air;

  double zeta = clamp (y/(h + 1e-30), 0., 1.);
  return base_lambda (zeta);
}

/* Exponentially decaying air velocity with continuity at the interface. */
static inline double air_velocity_ic (double x, double y)
{
  double h = disturbed_depth (x);
  double uI = base_u (1.)*(1. + DIST_AMP_U*sin (2.*pi*x/WAVE_LENGTH));

  if (y <= h)
    return fluid_velocity_ic (x, y);

  return uI*exp (-(y - h)/(AIR_DECAY + 1e-30));
}

/* Spatially varying target max level for mesh adaptation.
   The focus is the wall-bounded liquid region and the air-liquid interface. */
static inline int refRegion (double x, double y, double z)
{
  (void) x;
  (void) z;

  if (y <= 1.75 || fabs (y - 1.0) <= 0.8)
    return MAXLEVEL;
  if (y <= TOP_HEIGHT + 0.35)
    return max (MINLEVEL + 2, MAXLEVEL - 2);
  return MINLEVEL + 1;
}

/* --------------------------------------------------------------------------
   Main program.
   -------------------------------------------------------------------------- */
int main()
{
  /* Square periodic box. The physical flow is truncated above by the embedded top wall. */
  size (WAVE_LENGTH);
  origin (0., 0.);

  /* Dimensionless phase properties. */
  rho1 = 1.;
  rho2 = THIXO_RHOR;
  mu1 = 1.;
  mu2 = THIXO_MUR;

  /* Dimensionless thixotropic parameters. */
  thixo_T = THIXO_T;
  thixo_Gamma = THIXO_GAMMA;
  thixo_kappa = THIXO_KAPPA;
  thixo_a = THIXO_A;
  thixo_So = THIXO_SO;
  thixo_lambda_air = 0.;

  /* Surface tension is optional; for the first roll-wave case the default is zero. */
  f.sigma = THIXO_SIGMA;

  read_base_profile (PROFILE_FILE);
  init_grid (1 << INITLEVEL);
  periodic (right);

  NITERMAX = 200;
  CFL = 0.45;
  TOLERANCE = 5e-4;

  fprintf (ferr,
           "Thixotropic periodic roll-wave case: T=%g Gamma=%g kappa=%g a=%g FrV=%g So=%g rho_r=%g mu_r=%g\n",
           thixo_T, thixo_Gamma, thixo_kappa, thixo_a,
           THIXO_FRV, THIXO_SO, rho2, mu2);

  run();
}

/* --------------------------------------------------------------------------
   Initialise the disturbed liquid layer, the air phase above it, and the embedded
   top wall. The initialisation is iterated with mesh refinement so that the initial
   interface and wall region are reasonably resolved before time stepping starts.
   -------------------------------------------------------------------------- */
event init (i = 0)
{
  if (!restore ("restart")) {
    int changed;
    do {
      changed = 0;

      /* Pre-refine the liquid layer and the neighbourhood of the interface. */
      refine (y <= TOP_HEIGHT + 0.5 && level < MAXLEVEL - 1);
      refine (y <= 1.5 && level < MAXLEVEL);

      /* Free surface and top embedded wall. */
      fraction (f, disturbed_depth (x) - y);
      solid (cs, fs, TOP_HEIGHT - y);

      foreach() {
        double h = disturbed_depth (x);
        bool in_liquid = (y <= h);

        /* Reconstruct the initial structure field in the liquid only. */
        lambda_thixo[] = in_liquid ? fluid_lambda_ic (x, y) : thixo_lambda_air;

        /* Streamwise velocity in the liquid and vertically decaying air motion above it. */
        u.x[] = in_liquid ? fluid_velocity_ic (x, y) : air_velocity_ic (x, y);
        u.y[] = 0.;

        /* Hydrostatic pressure with zero datum at the disturbed interface. */
        p[] = in_liquid ? (h - y)/sq (THIXO_FRV)
                        : -rho2*(y - h)/sq (THIXO_FRV);

        /* Diagnostics used by mesh adaptation. */
        velFluidNorm[] = sqrt (sq (u.x[]) + sq (u.y[]))*f[];
        velAirNorm[]   = sqrt (sq (u.x[]) + sq (u.y[]))*(1. - f[]);
      }

      /* Build the conservative tracer flambda = f*lambda_thixo. */
      thixo_sync_tracer_from_lambda();

      boundary ((scalar *) {f, u.x, u.y, p, lambda_thixo, flambda,
                            velFluidNorm, velAirNorm});

      vorticity (u, omega_field);
      changed = adapt_wavelet ((scalar *) {f, u.x, u.y, lambda_thixo, omega_field},
                               (double[]) {F_ERR_OUT/10., UFLUID_ERR/3., UFLUID_ERR/3.,
                                           LAMBDA_ERR/3., OMEGA_ERR/3.},
                               MAXLEVEL, MINLEVEL).nf;
    } while (changed);

    fractions_cleanup (cs, fs);
    thixo_sync_tracer_from_lambda();
  }
}

/* --------------------------------------------------------------------------
   Body force.

   - x-component: applied only to the liquid phase, following the modelling note.
   - y-component: applied to both phases.
   -------------------------------------------------------------------------- */
event acceleration (i++)
{
  face vector av = a;

  foreach_face(x)
    av.x[] += THIXO_SO/sq (THIXO_FRV)*clamp ((f[] + f[-1])/2., 0., 1.);

  foreach_face(y)
    av.y[] -= 1./sq (THIXO_FRV);

  foreach() {
    velFluidNorm[] = sqrt (sq (u.x[]) + sq (u.y[]))*f[];
    velAirNorm[]   = sqrt (sq (u.x[]) + sq (u.y[]))*(1. - f[]);
  }
}

/* Iteration log. */
event logfile (i += 10)
{
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
}

/* Stop the run explicitly at MAXTIME and write a restart dump. */
event stop (t = MAXTIME)
{
  dump (file = "restart",
        list = {f, lambda_thixo, flambda, gammadot_thixo,
                mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p});
  return 1;
}

/* Full field dump. */
event snapshot (t += OUTPUT_DT)
{
  char name[128];
  sprintf (name, "dump-%0.6f", t);
  dump (file = name,
        list = {f, lambda_thixo, flambda, gammadot_thixo,
                mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p});
}

/* Gerris-format output is convenient for quadtree post-processing. */
event output_gfs_files (t += OUTPUT_DT)
{
  char name[128];
  sprintf (name, "out-%0.6f.gfs", t);
  FILE * fp = fopen (name, "w");
  output_gfs (fp, translate = true,
              list = {f, lambda_thixo, flambda, gammadot_thixo,
                      mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p});
  fclose (fp);
}

/* Output the reconstructed interface as a text file.
   The time is embedded in the file name so that snapshots do not overwrite one another. */
event output_interface (t += OUTPUT_DT)
{
  char partname[160];
  sprintf (partname, "interface-part-%0.6f-%d.dat", t, pid());
  FILE * fp = fopen (partname, "w");
  output_facets (f, fp);
  fclose (fp);

  char command[256];
  sprintf (command,
           "LC_ALL=C cat interface-part-%0.6f-*.dat > interface-all-%0.6f.dat",
           t, t);
  system (command);
}

/* Vertical profile at the streamwise midpoint of one period. */
event output_centerline_profile (t += OUTPUT_DT)
{
  char name[128];
  sprintf (name, "centerline-%0.6f.txt", t);
  FILE * fp = fopen (name, "w");

  double dy = L0/pow (2., MAXLEVEL);
  for (double yy = 0.; yy <= TOP_HEIGHT + 1e-12; yy += dy)
    fprintf (fp, "%g %g %g %g %g\n", yy,
             interpolate (u.x, WAVE_LENGTH/2., yy),
             interpolate (u.y, WAVE_LENGTH/2., yy),
             interpolate (lambda_thixo, WAVE_LENGTH/2., yy),
             interpolate (f, WAVE_LENGTH/2., yy));

  fclose (fp);
}

/* Log the maximum interface elevation in the period.
   The crest position is obtained in a second pass to avoid unreliable direct
   floating-point equality checks. */
event amplitude_log (i += 20)
{
  FILE * fp = fopen ("amplitude.txt", "a");

  position (f, ypos_interface, (coord) {0, 1});
  position (f, xpos_interface, (coord) {1, 0});

  double amp = statsf (ypos_interface).max;
  double xpos_leftmost = HUGE;

  foreach(reduction(min:xpos_leftmost))
    if (ypos_interface[] != nodata &&
        fabs (ypos_interface[] - amp) <= INTERFACE_PICK_TOL)
      xpos_leftmost = min (xpos_leftmost, xpos_interface[]);

  double xpos = (xpos_leftmost < HUGE/2. ? xpos_leftmost : nodata);
  fprintf (fp, "%g %g %g\n", t, xpos, amp);
  fclose (fp);
}

/* --------------------------------------------------------------------------
   Adaptive mesh refinement.

   The refinement criteria follow the spirit of the existing viscoplastic case:
   - keep the interface region fine,
   - focus on the lower liquid layer and the vicinity of the disturbed interface,
   - use flow diagnostics (vorticity, velocity magnitude and lambda) as indicators.
   -------------------------------------------------------------------------- */
event adapt (i++)
{
  vorticity (u, omega_field);
  boundary ((scalar *) {omega_field});

  adapt_wavelet_limited ((scalar *) {omega_field, velFluidNorm, velAirNorm,
                                     lambda_thixo},
                         (scalar *) {f},
                         (double[]) {OMEGA_ERR, UFLUID_ERR, UAIR_ERR,
                                     LAMBDA_ERR},
                         refRegion,
                         minlevel = MINLEVEL);

  /* Mild extra refinement below the top wall to avoid an abrupt jump of levels. */
  refine (y <= 1.25*TOP_HEIGHT && level < MAXLEVEL - 1);
}
