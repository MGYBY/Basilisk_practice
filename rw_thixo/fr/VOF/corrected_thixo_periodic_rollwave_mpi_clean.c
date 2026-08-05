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
// #define FILTERED 1
#include "corrected_two-phase-thixo.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "tag.h"
// #include "navier-stokes/perfs.h"
#include "./my-perfs.h"

#ifndef MAXLEVEL
#define MAXLEVEL 16
#endif
#ifndef MINLEVEL
#define MINLEVEL 4
#endif
#ifndef INITLEVEL
#define INITLEVEL 10
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
#ifndef normalVel
#define normalVel 0.20727
#endif
#ifndef THIXO_FRREAL
#define THIXO_FRREAL 0.750
#endif
#ifndef THIXO_FRV
#define THIXO_FRV (THIXO_FRREAL/normalVel)
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
#define THIXO_MU_MAX (3.0*(1.0/THIXO_SO)) // roughly estimated from Balmforth as 1/\epsilon
#endif
#ifndef THIXO_SIGMA
#define THIXO_SIGMA 0.003333
#endif

#ifndef DIST_AMP_H
#define DIST_AMP_H 2.e-1
#endif
#ifndef DIST_AMP_U
#define DIST_AMP_U 2.e-1
#endif
#ifndef DIST_WAVE_LENGTH
#define DIST_WAVE_LENGTH (2.0/THIXO_SO)
#endif
#ifndef DOMAIN_LENGTH
#define DOMAIN_LENGTH (123.50/THIXO_SO)
#endif
#ifndef TOP_HEIGHT
#define TOP_HEIGHT 5.650
#endif
#ifndef AIR_DECAY
#define AIR_DECAY 0.75
#endif
#ifndef OUTPUT_DT
#define OUTPUT_DT (1.00/THIXO_SO)
#endif
#ifndef DUMP_DT
#define DUMP_DT (2.00/THIXO_SO)
#endif
#ifndef MAXTIME
#define MAXTIME (240.0/THIXO_SO)
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
#define OMEGA_ERR 0.16
#endif
#ifndef LAMBDA_ERR
#define LAMBDA_ERR 2e-3
#endif
#ifndef UFLUID_ERR
#define UFLUID_ERR 6.e-2
#endif
#ifndef UAIR_ERR
#define UAIR_ERR 1.0e-1
#endif

static double * zprof = NULL, * uprof = NULL, * lprof = NULL;
static int nprof = 0;

scalar omega_field[];
scalar velFluidNorm[], velAirNorm[];
scalar xpos_interface[], ypos_interface[];
// scalar lambda_plot[];   // visualization-only: lambda masked outside mostly-liquid cells

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top]    = dirichlet(0.);
u.t[top]    = neumann(0.);

/* Embedded-ceiling boundary conditions. */
u.n[embed] = dirichlet(0.);
// u.t[embed] = dirichlet(0.);
u.t[embed] = neumann(0.);
p[embed]   = neumann(0.);
lambda_thixo[embed] = neumann(0.);
/* No explicit embed boundary on flambda during init; the conservative
   tracer is rebuilt cell-wise from lambda and f. */

lambda_thixo[bottom] = neumann(0.);
lambda_thixo[top]    = neumann(0.);

static inline double perturbation (double xx, double phase, double wavelength)
{
  if (xx<=wavelength && xx>=wavelength/2.0)
    return sin (2.*pi*xx/wavelength + phase - pi);
  return 0.0; 
}

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

// redundant, this is for periodic dist.
static inline double disturbed_depth (double xx, double phase, double wavelength)
{
  if (xx<=wavelength && xx>=wavelength/2.0)
    return 1.0*(1.0+DIST_AMP_H*(sin (2.*pi*xx/wavelength + phase - pi)));
  return 1.0; 
}

// this was for periodic dist.
// we change it for localized dist.
static inline double fluid_velocity_ic (double x, double y)
{
  double h = disturbed_depth (x,0.0,DIST_WAVE_LENGTH);
  if (y > h)
    return 0.;
  double zeta = clamp (y/(h + 1e-30), 0., 1.);
  double pu = perturbation (x, 0.0, DIST_WAVE_LENGTH);
  // return base_u (zeta)*(1. + DIST_AMP_U*sin(2.*pi*x/WAVE_LENGTH));
  return (base_u (zeta)*sqrt((1. + DIST_AMP_U*pu)));
}

static inline double fluid_lambda_ic (double x, double y)
{
  double h = disturbed_depth (x,0.0,DIST_WAVE_LENGTH);
  if (y > h)
    return thixo_lambda_air;
  double zeta = clamp (y/(h + 1e-30), 0., 1.);
  return base_lambda (zeta);
}

static inline double air_velocity_ic (double x, double y)
{
  double h = disturbed_depth (x,0.0,DIST_WAVE_LENGTH);
  double pu = perturbation (x, 0.0, DIST_WAVE_LENGTH);
  double uI = base_u (1.)*sqrt((1. + DIST_AMP_U*pu));
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
  size (DOMAIN_LENGTH);
  origin (0., 0.);

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

  read_base_profile (PROFILE_FILE);

  /* For MPI runs it is safer and clearer to set periodicity before
     grid allocation/initialisation. Basilisk will make existing and
     subsequently allocated fields periodic in the right/left direction. */
  periodic (right);
  init_grid (1 << INITLEVEL);

  /** Body-force gravity. This defines the acceleration vector $\pmb{a}$ in $\texttt{centered.h}$ file.*/
  // const face vector gravity[] = {(CHANNELSLOPE)*GRAV, (-CHANNELCOS)*GRAV, 0.0};
  //   const face vector gravity[] = {(CHANNELSLOPE)*GRAV*f, (-CHANNELCOS)*GRAV*f};
  // a = gravity;

  NITERMAX = 128;
  CFL = 0.475;
  TOLERANCE = 4.e-4;

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

      refine (y <= TOP_HEIGHT + 0.5 && level < MAXLEVEL - 1);
      refine (y <= 1.5 && level < MAXLEVEL);

      fprintf (ferr, "init: pass=%d stage=fractions\n", pass + 1);
      fflush (ferr);
      fraction (f, disturbed_depth(x,0.0,DIST_WAVE_LENGTH) - y);
      build_embedded_ceiling();

      fprintf (ferr, "init: pass=%d stage=fields\n", pass + 1);
      fflush (ferr);
      foreach() {
        double h = disturbed_depth (x,0.0,DIST_WAVE_LENGTH);
        // double h = H0*(1.0+DIST_AMP_H*perturbation (x, 0.0, DIST_WAVE_LENGTH));
        bool in_liquid = (y <= h);
        // lambda_thixo[] = in_liquid ? fluid_lambda_ic (x, y) : thixo_lambda_air;
        lambda_thixo[] = in_liquid ? fluid_lambda_ic (x, y) : f[];
        // u.x[] = in_liquid ? fluid_velocity_ic (x, y) : air_velocity_ic (x, y);
        u.x[] = in_liquid ? fluid_velocity_ic (x, y) : air_velocity_ic (x, y);
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

      fprintf (ferr, "init: pass=%d stage=adapt\n", pass + 1);
      fflush (ferr);
      astats st = adapt_wavelet ((scalar *) {f, u.x, u.y, lambda_thixo},
                                (double []) {F_ERR_OUT/1., UFLUID_ERR/3., UFLUID_ERR/3.,
                                             LAMBDA_ERR/3.},
                                MAXLEVEL, MINLEVEL);

      fprintf (ferr, "init: pass=%d nf=%d nc=%d\n", pass + 1, st.nf, st.nc);
      fflush (ferr);

      if (st.nf == 0 && st.nc == 0)
        break;
    }

    build_embedded_ceiling();
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
  foreach_face(y)
    av.y[] -= 1./sq(THIXO_FRV);

  foreach() {
    velFluidNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*f[];
    velAirNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*(1. - f[]);
  }
}

event logfile (i += 50)
{
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
}


event lambda_mass_log (i += 1000)
{
  double m_f = 0., m_flambda = 0.;
  double lmin =  HUGE, lmax = -HUGE;
  foreach (reduction(+:m_f) reduction(+:m_flambda) reduction(min:lmin) reduction(max:lmax)) {
    double fc = clamp (f[], 0., 1.);
    m_f += fc*dv();
    m_flambda += flambda[]*dv();
    if (fc > THIXO_LAMBDA_RECON_FMIN) {
      lmin = min (lmin, lambda_thixo[]);
      lmax = max (lmax, lambda_thixo[]);
    }
  }
  FILE * fp = fopen ("lambda_mass.txt", "a");
  fprintf (fp, "%g %d %.17g %.17g %.17g %.17g\n", t, i, m_f, m_flambda, lmin, lmax);
  fclose (fp);
}

/* Visualization-only field: avoids interpreting air/interface slivers as physical lambda. */
// event make_lambda_plot (i++)
// {
//   foreach()
//     lambda_plot[] = clamp(f[],0.,1.) > THIXO_LAMBDA_PLOT_FMIN ? lambda_thixo[] : thixo_lambda_air;
//   boundary ({lambda_plot});
// }

event snapshot ( t += DUMP_DT)
{
  char name[128];
  sprintf (name, "dump-%g", t);
  dump (file = name);
}

event output_gfs_files (t += OUTPUT_DT)
{
  char name[128];
  sprintf (name, "out-%g.gfs", t);
  FILE * fp = fopen (name, "w");
  output_gfs (fp, translate = true,
              // list = {f, lambda_thixo, flambda, lambda_plot, gammadot_thixo,
              //         mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p, velFluidNorm, velAirNorm});
              list = {f, lambda_thixo, flambda, gammadot_thixo,
                      mu_liquid_thixo, u.x, u.y, uf.x, uf.y, p, velFluidNorm, velAirNorm});
  fclose (fp);
}

event output_interface (t += OUTPUT_DT)
{
  /* MPI-safe interface output.

     Do not let all MPI ranks write to the same file.  Each rank writes
     a temporary file with a time/frame identifier, then rank 0 merges the
     files only after an MPI barrier.  The filtered field ff also removes
     tiny 0/1 VOF roundoff slivers before PLIC reconstruction. */
  static int iframe = 0;
  scalar ff[];
  foreach()
    ff[] = f[] <= 1e-10 ? 0. : f[] >= 1. - 1e-10 ? 1. : f[];
  boundary ({ff});

  char localname[128];
  sprintf (localname, "interface-i%09d-p%05d.dat", iframe, pid());
  FILE * fp = fopen (localname, "w");
  if (!fp) {
    fprintf (ferr, "Cannot open %s for interface output.\n", localname);
    exit (1);
  }
  output_facets (ff, fp);
  fclose (fp);

#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  if (pid() == 0) {
    char command[512];
    sprintf (command,
             "LC_ALL=C cat interface-i%09d-p*.dat > interface-all-i%09d-t%.12g.dat",
             iframe, iframe, t);
    int ret = system (command);
    if (ret != 0)
      fprintf (ferr, "Warning: cat command failed for interface frame %d.\n", iframe);

    sprintf (command, "rm -f interface-i%09d-p*.dat", iframe);
    ret = system (command);
    if (ret != 0)
      fprintf (ferr, "Warning: rm command failed for interface frame %d.\n", iframe);
  }

#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  iframe++;
}

// event output_centerline_profile (t += OUTPUT_DT)
// {
//   char name[128];
//   sprintf (name, "centerline-%g.txt", t);
//   FILE * fp = fopen (name, "w");
//   for (double yy = 0.; yy <= TOP_HEIGHT; yy += WAVE_LENGTH/pow(2., MAXLEVEL))
//     fprintf (fp, "%g %g %g %g %g\n", yy,
//              interpolate (u.x, WAVE_LENGTH/2., yy),
//              interpolate (u.y, WAVE_LENGTH/2., yy),
//              interpolate (lambda_thixo, WAVE_LENGTH/2., yy),
//              interpolate (f, WAVE_LENGTH/2., yy));
//   fclose (fp);
// }

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

// event adapt (i++)
// {
//   vorticity (u, omega_field);
//   boundary ((scalar *) {omega_field});
//
//   struct Adapt_limited adaptp;
//   adaptp.slist = (scalar *) {omega_field, velFluidNorm, velAirNorm, lambda_thixo};
//   adaptp.vol_frac = (scalar *) {f};
//   adaptp.max = (double []) {OMEGA_ERR, UFLUID_ERR, UAIR_ERR, LAMBDA_ERR};
//   adaptp.MLFun = refRegion;
//   adaptp.minlevel = MINLEVEL;
//   adaptp.padding = 1;
//   adaptp.list = all;
//   adapt_wavelet_limited (adaptp);
//
//   // refine (y <= 1.25*TOP_HEIGHT && level < MAXLEVEL - 1);
//   refine(y<=(4.550*WAVE_LENGTH/pow(2, MAXLEVEL)) && level<MAXLEVEL);
//   build_embedded_ceiling();
// }

int refRegion(double x,double y, double z){
  int lev;
  if( y < TOP_HEIGHT*0.975 )
    lev = MAXLEVEL;
  else
    lev = MINLEVEL+1;

  return lev;
}

// mesh adaptation
event adapt (i++) {
  /* MPI-sensitive cleanup before adaptation/output.  In MPI, ghost values
     at sub-domain boundaries must be up to date before computing gradients,
     vorticity and PLIC facets. */
  boundary ((scalar *) {f, lambda_thixo, flambda, mu_liquid_thixo,
                        u.x, u.y, velFluidNorm, velAirNorm});

  foreach() {
    velFluidNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*clamp(f[], 0., 1.);
    velAirNorm[]   = sqrt(sq(u.x[]) + sq(u.y[]))*(1. - clamp(f[], 0., 1.));
  }
  boundary ((scalar *) {velFluidNorm, velAirNorm});

  scalar omega[], gltn[];
  vector gradLambda_thixo[];

  vorticity (u, omega);
  boundary ((scalar *) {omega});

  gradients ({lambda_thixo}, {gradLambda_thixo});
  foreach()
    gltn[] = norm(gradLambda_thixo);
  boundary ((scalar *) {gltn});

  /* Keep f in the adaptation list as a primary VOF variable.  This is more
     robust in MPI than adapting only on diagnostic fields, especially for a
     long flat interface crossing many MPI sub-domain boundaries. */
  adapt_wavelet_limited ((scalar *) {f, u.x, u.y, lambda_thixo, omega}, {f},
                         (double []) {F_ERR_OUT, UFLUID_ERR, UFLUID_ERR,
                                      LAMBDA_ERR, OMEGA_ERR},
                         refRegion, minlevel = MINLEVEL);

  refine (y <= (4.550*DOMAIN_LENGTH/pow(2, MAXLEVEL)) && level < MAXLEVEL);

  /* After TREE adaptation/refinement, rebuild embedded fractions and
     reconstruct diagnostic lambda from the conservative VOF tracer. */
  build_embedded_ceiling();
  thixo_sync_lambda_from_tracer();

  foreach() {
    velFluidNorm[] = sqrt(sq(u.x[]) + sq(u.y[]))*clamp(f[], 0., 1.);
    velAirNorm[]   = sqrt(sq(u.x[]) + sq(u.y[]))*(1. - clamp(f[], 0., 1.));
  }
  boundary ((scalar *) {f, lambda_thixo, flambda, mu_liquid_thixo,
                        u.x, u.y, velFluidNorm, velAirNorm});
}
