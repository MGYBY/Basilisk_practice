/**
# */

// #include "grid/multigrid.h"
// #include "green-naghdi.h"
#include "saint-venant.h"

#define MAXLEVEL 10
#define MINLEVEL 8
#define MAXMAXLEVEL 10

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
double gaugeArr[15];
double dyGauge;
double gfCoord;
double gbCoord;
double heightThrehold = 5.0;
double islandHeight = 2.0;

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

static bool isInBlock(double x, double y)
{
  double Ly = 3.0, leftCoord = 40.0, rightCoord = 40.40, bottomCoord = 1.30, topCoord = 1.70;
  bool fluid = false;
  if (y < Ly && x >= leftCoord && x <= rightCoord && y >= bottomCoord && y <= topCoord)
  {
    fluid = true;
  }
  return fluid;
}

static double islandSlope(double x, double y)
{
  double Ly = 3.0, leftCoord = 40.0, rightCoord = 40.40, bottomCoord = 1.30, topCoord = 1.70, islandHeight = 2.0, So = 0.05011;
  double elevation = 0.0;
  if (y < Ly && x >= leftCoord && x <= rightCoord && y >= bottomCoord && y <= topCoord)
    elevation = islandHeight;
  else
    elevation = -So * x;

  return elevation;
}

void refine_zb(Point point, scalar zb)
{
  foreach_child()
      zb[] = islandSlope(x, y);
}

event init(i = 0)
{
  // for 2-D case 40x3m
  mask(y > Ly ? top : none);
  //friction coeff
  cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

  h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
  u.n[left] = dirichlet(normalVelocity);

  u.n[right] = neumann(0.);
  h[right] = neumann(0.);

  u.n[top] = neumann(0.);
  h[top] = neumann(0.);

  u.n[bottom] = neumann(0.);
  h[bottom] = neumann(0.);

  // related to the gauges
  dyGauge = (topCoord - bottomCoord) / 14.;
  for (int l = 1; l <= 15; l++)
  {
    gaugeArr[l - 1] = bottomCoord + (l - 1) * dyGauge;
  }

  gfCoord = (leftCoord - (Lx / pow(2, MAXMAXLEVEL)));
  gbCoord = (rightCoord + (Lx / pow(2, MAXMAXLEVEL)));

  zb.refine = refine_zb;

  foreach ()
  {
    zb[] = islandSlope(x, y);
    h[] = isInBlock(x, y) ? 0. : normalDepth;
    u.x[] = isInBlock(x, y) ? 0. : normalVelocity;
    u.y[] = 0.;
  }
}

/**
To implement an absorbing boundary condition, we add an area for $x >
12$ for which quadratic friction increases linearly with $x$. */

event friction(i++)
{

  foreach ()
  {
    double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
    foreach_dimension()
        u.x[] /= a;
  }
}

event outputGauges(i++)
{
  //first check if h is valid
  for (int l = 1; l <= 15; l++)
  {
    while (interpolate(h, gfCoord, gaugeArr[l - 1]) > heightThrehold)
    {
      gfCoord -= (Lx / pow(2, MAXMAXLEVEL));
    }

    while (interpolate(h, gbCoord, gaugeArr[l - 1]) > heightThrehold)
    {
      gbCoord += (Lx / pow(2, MAXMAXLEVEL));
    }
  }

  char name1[80];
  char name2[80];
  sprintf(name1, "gf");
  sprintf(name2, "gb");
  FILE *fp1 = fopen(name1, "a+");
  FILE *fp2 = fopen(name2, "a+");

  fprintf(fp1, "%g ", t);
  fprintf(fp2, "%g ", t);

  // output for the desired locations
  for (int l = 1; l <= 15; l++)
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
    // (isFluid(x, y, Ly, leftCoord, rightCoord, bottomCoord, topCoord)) ? (fprintf(fp1, "%g %g %g\n", x, y, h[]), fprintf(fp2, "%g %g\n", u.x[], u.y[])) : none;
    fprintf(fp1, "%g %g %g\n", x, y, h[]);
    fprintf(fp2, "%g %g\n", u.x[], u.y[]);
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
  astats s = adapt_wavelet({h}, (double[]){normalDepth / 100.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
