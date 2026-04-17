/**
# Periodic roll-wave test for the 1-D thixotropic multilayer solver

This case reads a steady-uniform base profile from a text file with three columns

~~~
 z   U(z)   lambda(z)
~~~

interpolates the profile to the multilayer grid, adds a small periodic free-surface
perturbation and advances the one-dimensional thixotropic multilayer system.
*/

// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "thixo_multilayer_basilisk.h"

#ifndef NLAYERS
#define NLAYERS 32
#endif
#ifndef LEVEL
#define LEVEL 8
#endif
#ifndef PERIODS
#define PERIODS 8.
#endif
#ifndef AMP
#define AMP 0.1
#endif

// handy variables added by me
#define DOMAINLENGTH 2.50
#define FrReal 0.70
#define normalVel (0.20727) // this value from the integration of steady state
#define FrV (FrReal/normalVel) 
#define normalDepth 1.00
#define MAXLEVEL 11
#define MINLEVEL 3
#define INITLEVEL 11
#define simTime 40.0
#define outputInterval 1.00
#define grav (1./(FrV*FrV))
#define yieldSurfThre 1.10e-8
#define colorbarMax (normalVel*3.00)
#define PLOTRANGEMAX (3.250*normalDepth)

// Dimensionless parameters of the periodic test.

// double FrV = 1.35;
double Tth = 10.;
double Gammath = 8.;
double kappath = 1e-4;
double ath = 0.20;
int relax_scheme = THIXO_RELAX_SSPRK2;

// Base-state profile read from the text file.

static double * zprof = NULL, * uprof = NULL, * lprof = NULL;
static int nprof = 0;

static int read_profile (const char * fname)
{
  FILE * fp = fopen (fname, "r");
  if (!fp) {
    perror (fname);
    return 0;
  }

  int cap = 128, nrows = 0;
  double * x = malloc (cap*sizeof(double));
  double * ux = malloc (cap*sizeof(double));
  double * lambda = malloc (cap*sizeof(double));

  if (!x || !ux || !lambda) {
    fprintf (ferr, "Allocation failure while reading %s\n", fname);
    fclose (fp);
    free (x), free (ux), free (lambda);
    return 0;
  }

  double a, b, c;
  while (fscanf (fp, "%lf %lf %lf", &a, &b, &c) == 3) {
    if (nrows == cap) {
      cap *= 2;

      double * xn = realloc (x, cap*sizeof(double));
      double * uxn = realloc (ux, cap*sizeof(double));
      double * ln = realloc (lambda, cap*sizeof(double));

      if (!xn || !uxn || !ln) {
        fprintf (ferr, "Reallocation failure while reading %s\n", fname);
        fclose (fp);
        free (xn ? xn : x);
        free (uxn ? uxn : ux);
        free (ln ? ln : lambda);
        return 0;
      }

      x = xn;
      ux = uxn;
      lambda = ln;
    }

    x[nrows] = a;
    ux[nrows] = b;
    lambda[nrows] = c;
    nrows++;
  }

  fclose (fp);

  if (nrows < 2) {
    fprintf (ferr, "Profile file %s has too few rows\n", fname);
    free (x), free (ux), free (lambda);
    return 0;
  }

  zprof = x;
  uprof = ux;
  lprof = lambda;
  nprof = nrows;

  fprintf (ferr, "Read %d profile rows from %s\n", nprof, fname);
  return 1;
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

static inline double layer_center_fraction (int l)
{
  double z = 0.;
  for (int k = 0; k < l; k++)
    z += layer[k];
  return z + 0.5*layer[l];
}

scalar depthGrad[], uAve[], yieldSurf[];

event init (i = 0)
{
  if (restore (file = "restart"))
    return 1;

  if (!read_profile ("profile.txt")) {
    fprintf (ferr, "Failed to read base profile from %s\n", "profile.txt");
    exit (1);
  }

  foreach() {
    zb[] = 0.;
    // double kx = 2.*pi*PERIODS/L0; // no need to keep this variable
    // h[] = normalDepth*(1.0+AMP*sin(2. * pi * x / DOMAINLENGTH));
    // we first perturb just velocity
    h[] = normalDepth;
    eta[] = h[];

    for (int l = 0; l < nl; l++) {
      vector uk = ul[l];
      scalar lam = lambdal[l];
      double zc = layer_center_fraction (l);
      // uk.x[] = interp1 (zprof, uprof, nprof, zc);
      uk.x[] = interp1 (zprof, uprof, nprof, zc)*(1.0+AMP*sin(2. * pi * x / DOMAINLENGTH));
      lam[] = interp1 (zprof, lprof, nprof, zc);
    }
  }
}

event acceleration (i++) {
  foreach(){
    depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    uAve[] = 0.0;
    double zCoord = zb[];
    double uVertGrad = 0.0;
    yieldSurf[] = 0.0;

      for (int l = 0; l < nl; l++) {
        zCoord += layer[l]*h[]*0.50;
        u = ul[l];
        u.x[] += dt*grav;

        uAve[] += pow((u.x[]*u.x[]+u.y[]*u.y[]),0.50)*layer[l];

        if (l>0 && l<(nl-1)) {
        // middle layers
        vector um = ul[l-1] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l+1]*h[]*0.50+layer[l]*h[]);
      }
      else if (l==(nl-1))
      {
        // top layer
        vector um = ul[l-1] ;
        vector up = ul[l] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l-1]*h[]*0.50+layer[l]*h[]*1.50);
      }
      else
      {
        // bottom layer
        vector um = ul[l] ;
        vector up = ul[l+1] ;
        uVertGrad = (up.x[]-um.x[])/(layer[l+1]*h[]*0.50+layer[l]*h[]*1.50);
      }

      if (uVertGrad<yieldSurfThre && yieldSurf[]<=0)
      {
        yieldSurf[] = zCoord;
      }

      zCoord += layer[l]*h[]*0.50;
      }
  }
    boundary ((scalar *){u});
}

/*
 * AMR here
 *
 */
event adapt1 (i++) {
  // adapt_wavelet({h, depthGrad, uAve}, (double[]){normalDepth/200.0, 0.01, normalVel/200.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // a more reasonable AMR criteria
  adapt_wavelet({h, depthGrad, uAve, yieldSurf}, (double[]){normalDepth/300.0, 0.007, normalVel/300.0, normalDepth/200.50}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  refine(x<=10.0*L0/pow(2, MAXLEVEL) && level<MAXLEVEL);
  refine(x>=L0-10.0*L0/pow(2, MAXLEVEL) && level<MAXLEVEL);
}

event logfile (i += 20)
{
  double hmin = HUGE, hmax = - HUGE;
  foreach (reduction (min:hmin) reduction (max:hmax)) {
    hmin = min (hmin, h[]);
    hmax = max (hmax, h[]);
  }
  fprintf (stderr, "i=%d t=%g dt=%g hmin=%g hmax=%g\n", i, t, dt, hmin, hmax);
}

event profiles (t = 0; t <= simTime; t+=outputInterval)
{
  // static int np = 0;
  char name[80];

  sprintf (name, "depth-%g.txt", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g\n", x, h[]);
  fclose (fp);

  sprintf (name, "layers-%g.txt", t);
  fp = fopen (name, "w");
  foreach() {
    for (int l = 0; l < nl; l++) {
      vector uk = ul[l];
      scalar lam = lambdal[l];
      fprintf (fp, "%g %d %g %g %g %g\n",
               x, l, layer_center_fraction (l)*h[], uk.x[], lam[], muifield[]);
    }
    fputc ('\n', fp);
  }
  fclose (fp);

  // np++;
}

/**
Post-processing visualization module
*/
// TODO: contours for lambda
void setup (FILE * fp)
{
  // FIXME: other customized color schemes
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
// 	   "# jet colormap\n"
// 	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
// 	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
// 	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
//        "# rainbow colormap\n"
       "set palette rgb 33,13,10\n"
      //  "load 'turbo.pal'\n"
//        "# green-red-violet colormap\n"
//        "set palette rgb 3,11,6\n"
	   "unset key\n"
	   "set cbrange [0:%g]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [%g:%g]\n"
	   "set yrange [0.0:%g]\n"
     "set lmargin at screen 0.1\n"
     "set rmargin at screen 0.9\n", colorbarMax, 0.0, DOMAINLENGTH, PLOTRANGEMAX
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %g'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    // fprintf (fp, "%g %g %g\n", x, max(0.5, z), u.x[]);
    for (int l = 0; l < nl; l++) {
      z += layer[l]*h[]*0.50;
      u = ul[l];
      // z += h[]/2.0;
      fprintf (fp, "%g %g %g \n", x, z, u.x[]);
      // z += h[]/2.0;
      z += layer[l]*h[]*0.50;
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);
}

event gnuplot (t += outputInterval; t <= simTime)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,320\n"
	   "set output 'plot-%g.png'\n", t);
  if (i == 0)
    setup (fp);
  plot (fp);
}

// event snapshots (t += 10.)
//   dump (file = "restart");

// event end (t = simTime)
//   dump (file = "final");

// event cleanup_case (i = end, last)
// {
//   free (zprof), zprof = NULL;
//   free (uprof), uprof = NULL;
//   free (lprof), lprof = NULL;
// }

int main()
{
  nl = NLAYERS;
  // N = 1 << LEVEL;
  init_grid(1 << (INITLEVEL));
  // L0 = 2.*pi*PERIODS;
  L0 = DOMAINLENGTH;
  // origin (0.);
  periodic (right);

  G = grav;
  thixo_T = Tth;
  thixo_Gamma = Gammath;
  thixo_kappa = kappath;
  thixo_a = ath;
  thixo_relaxation_scheme = relax_scheme;
  thixo_use_exact_relaxation = false;

  run();
}
