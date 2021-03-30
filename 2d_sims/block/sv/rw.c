/**
# Periodic wave propagation over an ellipsoidal shoal

We follow [Lannes and Marche, 2014](/src/references.bib#lannes2014)
and try to reproduce the experiment of [Berkhoff et al,
1982](/src/references.bib#berkhoff1982). The numerical wave tank is
25^2^ metres and periodic waves are generated on the left-hand side
and are damped on the right-hand side. */

// #include "grid/multigrid.h"
// #include "grid/cartesian.h"
// #include "green-naghdi.h"
#include "saint-venant.h"

#define MAXLEVEL 10
#define MINLEVEL 8
#define MAXMAXLEVEL 13

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.05;
double disPeriod = 0.933;
double simTime = 40.0;
double Lx = 42.0;
double Ly = 3.0;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
// block geometry
double leftCoord = 40.0;
double rightCoord = 40.40;
double bottomCoord = 1.30;
double topCoord = 1.70;
double gaugeArr[7];
double dyGauge;
double gfCoord;
double gbCoord;

// homogeneous bc
bid block;

int main()
{
  size(Lx);
  G = 9.81;
  init_grid(1 << MAXLEVEL);
  //N = 1000;
  CFL = 0.45; // CFL number should be sufficiently small
  // N = 1024;

  /**
  We turn off limiting to try to reduce wave dissipation. */
  // gradient = NULL;

  run();
}

event init(i = 0)
{
  // for 2-D case 40x3m
  mask(y > Ly ? top : none);
  // block bc
  mask(x >= leftCoord ? (x <= rightCoord ? (y >= bottomCoord ? (y <= topCoord ? block : none) : none) : none) : none);
  // L0 = 40.;

  h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
  u.n[left] = dirichlet(normalVelocity);

  u.n[right] = neumann(0.);
  h[right] = neumann(0.);

  u.n[top] = neumann(0.);
  h[top] = neumann(0.);

  u.n[bottom] = neumann(0.);
  h[bottom] = neumann(0.);

  /**
  The bathymetry is an inclined and skewed plane combined with an
  ellipsoidal shoal.

  ~~~gnuplot Bathymetry
  set term @PNG enhanced size 640,640 font ",8"
  set output 'bathy.png'
  set pm3d map
  set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
                        0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392,	  \
                        0.625 1 0.9333 0, 0.75 1 0.4392 0,		  \
                        0.875 0.9333 0 0, 1 0.498 0 0 )
  set size ratio -1
  set xlabel 'x (m)'
  set ylabel 'y (m)'
  splot [-10:12][-10:10]'end' u 1:2:4 w pm3d t ''
  # we remove the large border left by gnuplot using ImageMagick
  ! mogrify -trim +repage bathy.png
  ~~~
  */
  foreach ()
  {
    zb[] = -So * x;
    h[] = normalDepth;
    u.x[] = normalVelocity;
    u.y[] = 0.;
  }

  // related to the gauges
  dyGauge = (topCoord - bottomCoord - (Lx / pow(2, MAXLEVEL))) / 6.;
  for (int l = 1; l <= 7; l++)
  {
    gaugeArr[l - 1] = bottomCoord + (Lx / pow(2, MAXLEVEL) / 2.) + (l - 1) * dyGauge;
  }
  gfCoord = (leftCoord - (Lx / pow(2, MAXLEVEL) / 2.));
  gbCoord = (rightCoord + (Lx / pow(2, MAXLEVEL) / 2.));
}

// not using guages; use interploation instead
// /**
// We define the wave gauges for drag coeff calc. */
// Gauge gauges[] = {
//     {"GF1", gfCoord, gaugeArr[0]},
//     {"GF2", gfCoord, gaugeArr[1]},
//     {"GF3", gfCoord, gaugeArr[2]},
//     {"GF4", gfCoord, gaugeArr[3]},
//     {"GF5", gfCoord, gaugeArr[4]},
//     {"GF6", gfCoord, gaugeArr[5]},
//     {"GF7", gfCoord, gaugeArr[6]},
//     {"GB1", gbCoord, gaugeArr[0]},
//     {"GB2", gbCoord, gaugeArr[1]},
//     {"GB3", gbCoord, gaugeArr[2]},
//     {"GB4", gbCoord, gaugeArr[3]},
//     {"GB5", gbCoord, gaugeArr[4]},
//     {"GB6", gbCoord, gaugeArr[5]},
//     {"GB7", gbCoord, gaugeArr[6]},
//     {NULL}};

/**
To implement an absorbing boundary condition, we add an area for $x >
12$ for which quadratic friction increases linearly with $x$. */
// only use Euler backward time int for now

event friction(i++)
{
  cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
  foreach ()
  {
    double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
    foreach_dimension()
        u.x[] /= a;
  }
}
event outputGauges(i++)
{
  char name1[80];
  char name2[80];
  sprintf(name1, "gf");
  sprintf(name2, "gb");
  FILE *fp1 = fopen(name1, "a+");
  FILE *fp2 = fopen(name2, "a+");

  fprintf(fp1, "%g ", t);
  fprintf(fp2, "%g ", t);

  // output for the desired locations
  for (int l = 1; l <= 7; l++)
  {
    fprintf(fp1, "%g ", interpolate(h, gfCoord, gaugeArr[l - 1]));
    fprintf(fp2, "%g ", interpolate(h, gbCoord, gaugeArr[l - 1]));
  }

  fprintf(fp1, "\n");
  fprintf(fp2, "\n");

  fclose(fp1);
  fclose(fp2);
}

/**
At the end of the simulation, we output the maximum wave amplitudes
along the cross-sections corresponding with the experimental data. */

static bool isFluid(double x, double y, double Ly, double leftCoord, double rightCoord, double bottomCoord, double topCoord)
{
  bool fluid = true;
  if (y > Ly)
  {
    fluid = false;
  }
  if (x >= leftCoord && x <= rightCoord && y >= bottomCoord && y <= topCoord)
  {
    fluid = false;
  }
  return fluid;
}

event fieldOutput(t = 0; t <= simTime; t += 2)
{
  char name1[80];
  char name2[80];
  char name3[80];
  sprintf(name1, "out-h-%.0f", t);
  sprintf(name2, "out-u-%.0f", t);
  sprintf(name3, "profile-%.0f", t);
  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");
  FILE *fp3 = fopen(name3, "w");

  // compact form
  foreach ()
  {
    (isFluid(x, y, Ly, leftCoord, rightCoord, bottomCoord, topCoord)) ? (fprintf(fp1, "%g %g %g\n", x, y, h[]), fprintf(fp2, "%g %g\n", u.x[], u.y[])) : none;
  }

  for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
    fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}

// adaptivity
event adapt(i++)
{
  astats s = adapt_wavelet({h}, (double[]){normalDepth / 200.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
