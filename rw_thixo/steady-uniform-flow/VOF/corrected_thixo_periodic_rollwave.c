/**
# Periodic roll-wave case for the thixotropic VOF solver

This source file is a cleaned and corrected version of the recent
`mu_capped_noview` case. The main practical fix is the restart logic:
restart is disabled by default and the code reports explicitly whether a
restart was used.
*/

#include "grid/quadtree.h"
#include "adapt_wavelet_leave_interface_limited.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#define FILTERED 0
#include "corrected_two-phase-thixo.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "tag.h"
#include "navier-stokes/perfs.h"
// #include "./my-perfs.h"

#ifndef MAXLEVEL
#define MAXLEVEL 7
#endif
#ifndef MINLEVEL
#define MINLEVEL 7
#endif
#ifndef INITLEVEL
#define INITLEVEL 7
#endif

#ifndef THIXO_T
#define THIXO_T 100.0
#endif
#ifndef THIXO_GAMMA
#define THIXO_GAMMA 8.0
#endif
#ifndef THIXO_KAPPA
#define THIXO_KAPPA 1.0e-4
#endif
#ifndef THIXO_A
#define THIXO_A 0.2
#endif
#ifndef THIXO_FRV
#define THIXO_FRV 4.00
#endif
#ifndef THIXO_SO
#define THIXO_SO 0.05
#endif
#ifndef THIXO_RHOR
#define THIXO_RHOR 0.01
#endif
#ifndef THIXO_MUR
#define THIXO_MUR 0.02
#endif
#ifndef THIXO_MU_MAX
#define THIXO_MU_MAX (4.0*(1.0/THIXO_SO)) // roughly estimated from Balmforth as 1/\epsilon
#endif
#ifndef THIXO_SIGMA
#define THIXO_SIGMA 0.001
#endif

#ifndef DIST_AMP_H
// #define DIST_AMP_H 1e-1
// modified for S-U flow
#define DIST_AMP_H 0.0
#endif
#ifndef DIST_AMP_U
// #define DIST_AMP_U 1e-1
// modified for S-U flow
#define DIST_AMP_U 0.0
#endif
#ifndef WAVE_LENGTH
// #define WAVE_LENGTH (2.50/THIXO_SO)
// modified for S-U flow
#define WAVE_LENGTH (4.0)
#endif
#ifndef TOP_HEIGHT
#define TOP_HEIGHT 4.0
#endif
#ifndef AIR_DECAY
#define AIR_DECAY 0.75
#endif
#ifndef OUTPUT_DT
#define OUTPUT_DT 5.0
#endif
#ifndef MAXTIME
#define MAXTIME 405.10
#endif
#ifndef PROFILE_FILE
#define PROFILE_FILE "profile.txt"
#endif

#ifndef USE_RESTART
#define USE_RESTART 0
#endif
#ifndef RESTART_FILE
#define RESTART_FILE "thixo_restart"
#endif
#ifndef INIT_MAX_PASSES
#define INIT_MAX_PASSES 12
#endif

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

static double * zprof = NULL, * uprof = NULL, * lprof = NULL;
static int nprof = 0;

scalar omega_field[];
scalar velFluidNorm[], velAirNorm[];
scalar xpos_interface[], ypos_interface[];

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top]    = dirichlet(0.);
u.t[top]    = neumann(0.);

/* Embedded-ceiling boundary conditions. */
u.n[embed] = dirichlet(0.);
// modified here
u.t[embed] = neumann(0.);
p[embed]   = neumann(0.);
lambda_thixo[embed] = neumann(0.);
/* No explicit embed boundary on flambda during init; the conservative
   tracer is rebuilt cell-wise from lambda and f. */

lambda_thixo[bottom] = neumann(0.);
lambda_thixo[top]    = neumann(0.);

static inline void build_embedded_ceiling()
{
  solid (cs, fs, TOP_HEIGHT - y);
  fractions_cleanup (cs, fs);
}

static void read_base_profile (const char * fname)
{
  FILE * fp = fopen (fname, "r");
  if (!fp) {
    fprintf (ferr, "Cannot open base-state profile file '%s'\n", fname);
    exit (1);
  }

  int cap = 256;
  zprof = (double *) malloc (cap*sizeof(double));
  uprof = (double *) malloc (cap*sizeof(double));
  lprof = (double *) malloc (cap*sizeof(double));
  if (!zprof || !uprof || !lprof) {
    fprintf (ferr, "Memory allocation failure while reading base-state profile.\n");
    exit (1);
  }

  while (1) {
    double z, u, lam;
    int ret = fscanf (fp, "%lf %lf %lf", &z, &u, &lam);
    if (ret == EOF)
      break;
    if (ret != 3) {
      fprintf (ferr, "Malformed line in base-state profile '%s'.\n", fname);
      exit (1);
    }
    if (nprof >= cap) {
      cap *= 2;
      zprof = (double *) realloc (zprof, cap*sizeof(double));
      uprof = (double *) realloc (uprof, cap*sizeof(double));
      lprof = (double *) realloc (lprof, cap*sizeof(double));
      if (!zprof || !uprof || !lprof) {
        fprintf (ferr, "Memory reallocation failure while reading base-state profile.\n");
        exit (1);
      }
    }
    zprof[nprof] = z;
    uprof[nprof] = u;
    lprof[nprof] = lam;
    nprof++;
  }
  fclose (fp);

  if (nprof < 2) {
    fprintf (ferr, "Base-state profile '%s' must contain at least two rows.\n", fname);
    exit (1);
  }
}

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

  double t = (zeta - zprof[lo])/(zprof[hi] - zprof[lo] + 1e-30);
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

static inline double disturbed_depth (double x)
{
  return 1. + DIST_AMP_H*sin(2.*pi*x/WAVE_LENGTH);
}

static inline double fluid_velocity_ic (double x, double y)
{
  double h = disturbed_depth (x);
  if (y > h)
    return 0.;
  double zeta = clamp (y/(h + 1e-30), 0., 1.);
  return base_u (zeta)*(1. + DIST_AMP_U*sin(2.*pi*x/WAVE_LENGTH));
}

static inline double fluid_lambda_ic (double x, double y)
{
  double h = disturbed_depth (x);
  if (y > h)
    return thixo_lambda_air;
  double zeta = clamp (y/(h + 1e-30), 0., 1.);
  return base_lambda (zeta);
}

static inline double air_velocity_ic (double x, double y)
{
  double h = disturbed_depth (x);
  double uI = base_u (1.)*(1. + DIST_AMP_U*sin(2.*pi*x/WAVE_LENGTH));
  if (y <= h)
    return fluid_velocity_ic (x, y);
  return uI*exp (-(y - h)/(AIR_DECAY + 1e-30));
}

// static inline int refRegion (double x, double y, double z)
// {
//   (void) x; (void) z;
//   if (y <= 1.75 || fabs(y - 1.) <= 0.8)
//     return MAXLEVEL;
//   if (y <= TOP_HEIGHT + 0.35)
//     return max (MINLEVEL + 2, MAXLEVEL - 2);
//   return MINLEVEL + 1;
// }

int main()
{
  size (WAVE_LENGTH);
  // origin (0., 0.);

  rho1 = 1.;
  rho2 = THIXO_RHOR;
  mu1 = 1.;
  mu2 = THIXO_MUR;

  thixo_T = THIXO_T;
  thixo_Gamma = THIXO_GAMMA;
  thixo_kappa = THIXO_KAPPA;
  thixo_a = THIXO_A;
  thixo_So = THIXO_SO;
  thixo_FrV = THIXO_FRV;
  thixo_lambda_air = 0.;

  f.sigma = THIXO_SIGMA;

  // start from "at rest" for steady-uniform flows
  read_base_profile (PROFILE_FILE);
  init_grid (1 << INITLEVEL);
  periodic (right);

  /** Body-force gravity. This defines the acceleration vector $\pmb{a}$ in $\texttt{centered.h}$ file.*/
  // const face vector gravity[] = {(CHANNELSLOPE)*GRAV, (-CHANNELCOS)*GRAV, 0.0};
  //   const face vector gravity[] = {(CHANNELSLOPE)*GRAV*f, (-CHANNELCOS)*GRAV*f};
  // a = gravity;

  NITERMAX = 128;
  CFL = 0.475;
  TOLERANCE = 2.5e-4;

  fprintf (ferr,
           "Thixotropic periodic roll-wave case: T=%g Gamma=%g kappa=%g a=%g FrV=%g So=%g rho_r=%g mu_r=%g\n",
           thixo_T, thixo_Gamma, thixo_kappa, thixo_a, THIXO_FRV, THIXO_SO, rho2, mu2);
  fprintf (ferr, "USE_RESTART=%d RESTART_FILE=%s INIT_MAX_PASSES=%d\n",
           USE_RESTART, RESTART_FILE, INIT_MAX_PASSES);
  fflush (ferr);

  run();
  return 0;
}

event init (i = 0)
{
  bool restarted = false;
#if USE_RESTART
  restarted = restore (RESTART_FILE);
#endif

  fprintf (ferr, "init: entered, restarted=%d\n", restarted);
  fflush (ferr);

  if (!restarted) {
    for (int pass = 0; pass < INIT_MAX_PASSES; pass++) {
      fprintf (ferr, "init: pass=%d stage=refine-start\n", pass + 1);
      fflush (ferr);

      // modified for S-U flows
      // refine (y <= TOP_HEIGHT + 0.5 && level < MAXLEVEL - 1);
      refine (y <= 1.5 && level < MAXLEVEL);

      fprintf (ferr, "init: pass=%d stage=fractions\n", pass + 1);
      fflush (ferr);
      fraction (f, disturbed_depth(x) - y);
      // build_embedded_ceiling();

      fprintf (ferr, "init: pass=%d stage=fields\n", pass + 1);
      fflush (ferr);
      foreach() {
        double h = disturbed_depth (x);
        bool in_liquid = (y <= h);
        lambda_thixo[] = in_liquid ? fluid_lambda_ic (x, y) : f[]; // modified here for lambda init
        // lambda_thixo[] = in_liquid ? 0.10 : f[]; // modified here for lambda init
        // modified for S-U flows
        // u.x[] = in_liquid ? pow(fluid_velocity_ic (x, y),(3./2.)) : pow(air_velocity_ic (x, y),(3./2.));
        u.x[] = f[]>0?0.10:0.20;
        u.y[] = 0.;
        p[] = in_liquid ? (h - y)/sq(THIXO_FRV) : -rho2*(y - h)/sq(THIXO_FRV);
        velFluidNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*f[];
        velAirNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*(1. - f[]);
      }

      fprintf (ferr, "init: pass=%d stage=sync\n", pass + 1);
      fflush (ferr);
      fprintf (ferr, "init: pass=%d stage=sync-cell\n", pass + 1); fflush (ferr);
      thixo_sync_tracer_from_lambda_cells();
      fprintf (ferr, "init: pass=%d stage=sync-boundary\n", pass + 1); fflush (ferr);
      boundary ((scalar *) {f, lambda_thixo, mu_liquid_thixo, u.x, u.y,
                            velFluidNorm, velAirNorm});

      /* Avoid vorticity during initial AMR. This keeps the first-pass
         adaptation away from embedded/ghost-value sensitivities. */
      foreach()
        omega_field[] = 0.;
      boundary ((scalar *) {omega_field});

      // fprintf (ferr, "init: pass=%d stage=adapt\n", pass + 1);
      // fflush (ferr);
      // astats st = adapt_wavelet ((scalar *) {f, u.x, u.y, lambda_thixo},
      //                           (double []) {F_ERR_OUT/10., UFLUID_ERR/3., UFLUID_ERR/3.,
      //                                        LAMBDA_ERR/3.},
      //                           MAXLEVEL, MINLEVEL);

      // fprintf (ferr, "init: pass=%d nf=%d nc=%d\n", pass + 1, st.nf, st.nc);
      // fflush (ferr);

      // if (st.nf == 0 && st.nc == 0)
      //   break;
    }

    // modified for S-U flows
    // build_embedded_ceiling();
    fprintf (ferr, "init: final-sync-cell\n"); fflush (ferr);
    thixo_sync_tracer_from_lambda_cells();
    fprintf (ferr, "init: final-sync-boundary\n"); fflush (ferr);
    boundary ((scalar *) {f, lambda_thixo, mu_liquid_thixo, u.x, u.y,
                          velFluidNorm, velAirNorm});
    vorticity (u, omega_field);
    boundary ((scalar *) {omega_field});

    fprintf (ferr, "init: fresh initialisation complete\n");
    fflush (ferr);
  }
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] += THIXO_SO/sq(THIXO_FRV)*clamp((f[] + f[-1])/2., 0., 1.);
    // a more reasonable one but with more time restriction
    // av.x[] += THIXO_SO/sq(THIXO_FRV);
  foreach_face(y)
    av.y[] -= 1./sq(THIXO_FRV);

  foreach() {
    velFluidNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*f[];
    velAirNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*(1. - f[]);
  }
}

event logfile (i += 5)
{
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
}

// event snapshot ( t += OUTPUT_DT)
// {
//   char name[128];
//   sprintf (name, "dump-%g", t);
//   dump (file = name,
//         list = {f, lambda_thixo, flambda, gammadot_thixo,
//                 mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p});
// }

// event output_gfs_files (t += OUTPUT_DT)
// {
//   char name[128];
//   sprintf (name, "out-%g.gfs", t);
//   FILE * fp = fopen (name, "w");
//   output_gfs (fp, translate = true,
//               list = {f, lambda_thixo, flambda, gammadot_thixo,
//                       mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p,velFluidNorm,velAirNorm});
//   fclose (fp);
// }

event output_interface (t += OUTPUT_DT)
{
  char localname[80];
  sprintf (localname, "interface-%d.dat", pid());
  FILE * fp = fopen (localname, "w");
  output_facets (f, fp);
  fclose (fp);

  if (pid() == 0) {
    char command[256];
    sprintf (command,
             "LC_ALL=C cat interface-[0-9]*.dat > interface-all-%g.dat; rm -f interface-[0-9]*.dat",
             t);
    system (command);
  }
}

event output_centerline_profile (t += OUTPUT_DT)
{
  char name[128];
  sprintf (name, "centerline-%g.txt", t);
  FILE * fp = fopen (name, "w");
  for (double yy = 0.; yy <= TOP_HEIGHT; yy += WAVE_LENGTH/pow(2., MAXLEVEL))
    fprintf (fp, "%g %g %g %g %g\n", yy,
             interpolate (u.x, WAVE_LENGTH/2., yy),
             interpolate (u.y, WAVE_LENGTH/2., yy),
             interpolate (lambda_thixo, WAVE_LENGTH/2., yy),
             interpolate (f, WAVE_LENGTH/2., yy));
  fclose (fp);
}

event amplitude_log (i += 25)
{
  FILE * fp = fopen ("amplitude.txt", "a");
  FILE * fpV = fopen ("max_velocity.txt", "a");
  position (f, ypos_interface, (coord) {0,1});
  position (f, xpos_interface, (coord) {1,0});

  double amp = statsf(ypos_interface).max;
  double ampV = statsf(velFluidNorm).max;
  double xpos = 0.;
  foreach (reduction(max:xpos))
    if (ypos_interface[] == amp)
      xpos = xpos_interface[];

  fprintf (fp, "%g %g %g\n", t, xpos, amp);
  fprintf (fpV, "%g %g\n", t, ampV);
  fclose (fp);
  fclose (fpV);
}

event stop (t = MAXTIME)
{
  fprintf (ferr, "stop: reached MAXTIME=%g\n", MAXTIME);
  fflush (ferr);
}
