/**
# Periodic roll-wave case for the Balmforth-consistent 1-D thixotropic multilayer solver

The base state is read from a text file with three columns

~~~
  z   U(z)   lambda(z)
~~~

written in the **explicit-$Fr_V$ dimensionless variables of the multilayer
solver** and by the Python script `thixo_base_state_explicitFrV.py`.

A small sinusoidal perturbation is applied to the free surface only.  The layer
velocities and structure variables are initialized from the steady base-state
vertical profiles.
*/

#include "grid/multigrid1D.h"
#include "thixo_multilayer_basilisk_explicitFrV_fixed.h"

#ifndef NLAYERS
#define NLAYERS 64
#endif
#ifndef INITLEVEL
#define INITLEVEL 7
#endif
#ifndef AMP
#define AMP 0.10
#endif
#ifndef DOMAINLENGTH
#define DOMAINLENGTH 2.5
#endif
#ifndef SIMTIME
#define SIMTIME 40.
#endif
#ifndef OUTPUTINTERVAL
#define OUTPUTINTERVAL 1.
#endif
#ifndef PROFILE_FILE
#define PROFILE_FILE "profile.txt"
#endif

#define FrReal 0.70
#define normalVel (0.20727)
#define normalDepth 1.00
#define MAXLEVEL 11
#define MINLEVEL 3
#define yieldSurfThre 1.10e-8
#define colorbarMax (normalVel*3.00)
#define PLOTRANGEMAX (3.250*normalDepth)

// Material parameters in the explicit-FrV formulation.
double Trelax = 100.;
double Gammath = 8.;
double kappath = 1e-4;
double ath = 0.20;
double FrV_input = 4.00;
int relax_scheme = THIXO_RELAX_SSPRK2;

double H0 = 1.;

double * zprof = NULL, * uprof = NULL, * lprof = NULL;
int nprof = 0;

static int read_profile (const char * fname)
{
  FILE * fp = fopen (fname, "r");
  if (!fp)
    return 0;

  int cap = 128, nrows = 0;
  double * z = malloc (cap*sizeof(double));
  double * u = malloc (cap*sizeof(double));
  double * l = malloc (cap*sizeof(double));
  if (!z || !u || !l) {
    fclose (fp);
    free (z), free (u), free (l);
    return 0;
  }

  while (1) {
    double zz, uu, ll;
    int ret = fscanf (fp, "%lf %lf %lf", &zz, &uu, &ll);
    if (ret != 3)
      break;

    if (nrows == cap) {
      cap *= 2;
      double * zn = realloc (z, cap*sizeof(double));
      double * un = realloc (u, cap*sizeof(double));
      double * ln = realloc (l, cap*sizeof(double));
      if (!zn || !un || !ln) {
        fclose (fp);
        free (zn ? zn : z);
        free (un ? un : u);
        free (ln ? ln : l);
        return 0;
      }
      z = zn, u = un, l = ln;
    }

    z[nrows] = zz;
    u[nrows] = uu;
    l[nrows] = ll;
    nrows++;
  }
  fclose (fp);

  if (nrows < 2) {
    free (z), free (u), free (l);
    return 0;
  }

  // Replace any previous profile.
  free (zprof), free (uprof), free (lprof);
  zprof = z, uprof = u, lprof = l, nprof = nrows;
  fprintf (ferr, "Read %d profile rows from %s\n", nprof, fname);
  return 1;
}

static int read_profile_any (void)
{
  const char * candidates[] = {
    PROFILE_FILE,
    "./profile.txt",
    "../profile.txt",
    "../../profile.txt",
    NULL
  };
  for (int i = 0; candidates[i]; i++)
    if (read_profile (candidates[i])) {
      fprintf (ferr, "Using base profile file: %s\n", candidates[i]);
      return 1;
    }
  perror (PROFILE_FILE);
  return 0;
}

static double interp1 (double * x, double * y, int n, double xi)
{
  if (xi <= x[0])
    return y[0];
  if (xi >= x[n - 1])
    return y[n - 1];
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi)/2;
    if (x[mid] <= xi)
      lo = mid;
    else
      hi = mid;
  }
  double t = (xi - x[lo])/(x[hi] - x[lo]);
  return y[lo] + t*(y[hi] - y[lo]);
}

// Average a tabulated profile over the layer interval [a,b].
static double interp_avg (double * x, double * y, int n, double a, double b)
{
  int m = 8;
  double s = 0.;
  for (int i = 0; i <= m; i++) {
    double z = a + (b - a)*i/(double) m;
    double w = (i == 0 || i == m ? 0.5 : 1.);
    s += w*interp1 (x, y, n, z);
  }
  return s/m;
}

event init (i = 0)
{
  if (restore (file = "restart"))
    return 1;

  if (!read_profile_any()) {
    fprintf (ferr, "Failed to read the base-state profile.\n");
    exit (1);
  }

  foreach() {
    zb[] = 0.;
    h[] = H0*(1. + AMP*sin (2.*pi*x/L0));
    eta[] = zb[] + h[];

    for (int l = 0; l < nl; l++) {
      double z0 = 0.;
      for (int k = 0; k < l; k++)
        z0 += layer[k];
      double z1 = z0 + layer[l];

      vector uk = ul[l];
      scalar lam = lambdal[l];
      uk.x[] = interp_avg (zprof, uprof, nprof, z0, z1);
      lam[] = interp_avg (zprof, lprof, nprof, z0, z1);
    }
  }
}

event logfile (i += 20)
{
  double hmin = HUGE, hmax = -HUGE, umax = -HUGE, lmin = HUGE, lmax = -HUGE;
  foreach (reduction (min:hmin) reduction (max:hmax)
           reduction (max:umax) reduction (min:lmin) reduction (max:lmax)) {
    hmin = min (hmin, h[]);
    hmax = max (hmax, h[]);
    for (int l = 0; l < nl; l++) {
      vector uk = ul[l];
      scalar lam = lambdal[l];
      umax = max (umax, fabs (uk.x[]));
      lmin = min (lmin, lam[]);
      lmax = max (lmax, lam[]);
    }
  }
  fprintf (ferr, "i=%d t=%g dt=%g hmin=%g hmax=%g umax=%g lmin=%g lmax=%g\n",
           i, t, dt, hmin, hmax, umax, lmin, lmax);
}

event profiles (t = 0; t <= SIMTIME; t += OUTPUTINTERVAL)
{
  char name[80];
  sprintf (name, "depth-%g.txt", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g\n", x, h[]);
  fclose (fp);

  sprintf (name, "layers-%g.txt", t);
  fp = fopen (name, "w");
  foreach() {
    double z = zb[];
    for (int l = 0; l < nl; l++) {
      z += 0.5*layer[l]*h[];
      vector uk = ul[l];
      scalar lam = lambdal[l];
      fprintf (fp, "%g %d %g %g %g %g\n", x, l, z, uk.x[], lam[], muifield[]);
      z += 0.5*layer[l]*h[];
    }
    fputc ('\n', fp);
  }
  fclose (fp);
}

// event end (t = SIMTIME)
// {
//   dump (file = "final");
// }

void setupVel (FILE * fp)
{
  fprintf (fp,
           "set pm3d map interpolate 2,2\n"
           "set palette rgb 33,13,10\n"
           "unset key\n"
           "set cbrange [0:%g]\n"
           "set xlabel 'x'\n"
           "set ylabel 'height'\n"
           "set xrange [%g:%g]\n"
           "set yrange [0.0:%g]\n"
           "set lmargin at screen 0.1\n"
           "set rmargin at screen 0.9\n",
           colorbarMax, 0.0, DOMAINLENGTH, PLOTRANGEMAX);
}

void plotVel (FILE * fp)
{
  fprintf (fp,
           "set title 't = %g'\n"
           "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    for (int l = 0; l < nl; l++) {
      z += layer[l]*h[]*0.50;
      vector u = ul[l];
      fprintf (fp, "%g %g %g\n", x, z, u.x[]);
      z += layer[l]*h[]*0.50;
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

void setupLam (FILE * fp)
{
  fprintf (fp,
           "set pm3d map interpolate 2,2\n"
           "set palette rgb 33,13,10\n"
           "unset key\n"
           "set cbrange [0:%g]\n"
           "set xlabel 'x'\n"
           "set ylabel 'height'\n"
           "set xrange [%g:%g]\n"
           "set yrange [0.0:%g]\n"
           "set lmargin at screen 0.1\n"
           "set rmargin at screen 0.9\n",
           1.0, 0.0, DOMAINLENGTH, PLOTRANGEMAX);
}

void plotLam (FILE * fp)
{
  fprintf (fp,
           "set title 't = %g'\n"
           "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    for (int l = 0; l < nl; l++) {
      scalar lam = lambdal[l];
      z += layer[l]*h[]*0.50;
      fprintf (fp, "%g %g %g\n", x, z, lam[]);
      z += layer[l]*h[]*0.50;
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event gnuplotVelocity (t += OUTPUTINTERVAL; t <= SIMTIME)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
           "set term pngcairo font \",10\" size 1024,320\n"
           "set output 'plotVel-%g.png'\n", t);
  if (i == 0)
    setupVel (fp);
  plotVel (fp);
}

event gnuplotLambda (t += OUTPUTINTERVAL; t <= SIMTIME)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
           "set term pngcairo font \",10\" size 1024,320\n"
           "set output 'plotLam-%g.png'\n", t);
  if (i == 0)
    setupLam (fp);
  plotLam (fp);
}

int main()
{
  nl = NLAYERS;
  init_grid (1 << INITLEVEL);
  L0 = DOMAINLENGTH;
  periodic (right);

  thixo_T = Trelax;
  thixo_Gamma = Gammath;
  thixo_kappa = kappath;
  thixo_a = ath;
  FrV = FrV_input;
  G = 1./sq (FrV);
  thixo_relaxation_scheme = relax_scheme;
  thixo_use_exact_relaxation = false;

  // TODO: add CFL number here ~0.45

  run();
}
