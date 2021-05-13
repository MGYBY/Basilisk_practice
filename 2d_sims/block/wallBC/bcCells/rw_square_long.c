// #include "grid/multigrid.h"
//#include "grid/cartesian.h"
// #include "green-naghdi.h"
#include "saint-venant.h"

#define MAXLEVEL 12
#define MINLEVEL 10
#define MAXMAXLEVEL 14

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.20;
double disPeriod = 0.933;
double simTime = 100.0;
double Lx = 162.0;
double Ly = 3.0;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
// block geometry
double leftCoord = 160.0;
double rightCoord = 160.40;
double bottomCoord = 1.30;
double topCoord = 1.70;
double coordTol = 1e-7;
// drag force calc
double rhoWater = 1e3;

// for debugging
double blockXLower = 1e2;
double blockXUpper = 0.0;

// homogeneous bc
bid block;
//free-slip bc
//not working???
// u.n[block] = dirichlet(0);


int main()
{
  size(Lx);
  G = 9.81;
  init_grid(1 << 13);
  //N = 1000;
  CFL = 0.495; // CFL number should be sufficiently small
//   theta = 1.5; // use Sweby limiter
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

  cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

  foreach ()
  {
    zb[] = 0.0;
    h[] = normalDepth;
    u.x[] = normalVelocity;
    u.y[] = 0.;
  }

  foreach_boundary (block){
    if (x<blockXLower)
    {
      blockXLower = x;
    }

    if (x>blockXUpper)
    {
      blockXUpper = x;
    }
  }
}

/**
To implement an absorbing boundary condition, we add an area for $x >
12$ for which quadratic friction increases linearly with $x$. */
// only use Euler backward time int for now

event friction(i++)
{
  // double uMed;
  
  foreach ()
  {
    // double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
      double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
      u.x[] = (u.x[] + G * So * dt) / a;
      u.y[] /= a;
      // uMed = u.x[] + dt *  (-(cf / 2.) * u.x[] * norm(u) / h[] + gravityCoeff*So);
      // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * u.x[] * sqrt(sq(uMed)+sq(u.y[])) / h[] + gravityCoeff*So);
      // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * u.x[] * sqrt(sq(uMed)+sq(u.y[])) / h[] + gravityCoeff*So);    

      // uMed = u.y[] + dt *  (-(cf / 2.) * u.y[] * norm(u) / h[]);
      // uMed = (3. / 4.) * u.y[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * u.y[] * sqrt(sq(uMed)+sq(u.x[])) / h[]);
      // u.y[] = (1. / 3.) * u.y[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * u.y[] * sqrt(sq(uMed)+sq(u.x[])) / h[]);  
  }
}

event outputDrag(i++)
{
  // char name1[80];
  // drag force in x dir
  double fd = 0.0;
  // name1 = "dragForce";
  FILE *fp1 = fopen("dragForce", "a+");

  foreach_boundary (block){
    if (x==blockXLower)
    {
      fd += 0.5*rhoWater*gravityCoeff*(Delta*h[])*h[];
    }
    else if (x==blockXUpper)
    {
      fd -= 0.5*rhoWater*gravityCoeff*(Delta*h[])*h[];
    }
  }

  /**
  // for debugging
  // foreach_boundary (block){
  //   fprintf(fp1, "%g", h[]);
  // }
  {
    fprintf(fp1, "%g %g", blockXLower, blockXUpper);
  }
  // for debugging
  **/

  fprintf(fp1, "%g %g %g\n", t, fd, fd/(0.5*rhoWater*(topCoord-bottomCoord)*sq(normalVelocity)*normalDepth));

  fclose(fp1);
}

/**
At the end of the simulation, we output the maximum wave amplitudes
along the cross-sections corresponding with the experimental data. */

// static bool isFluid(double x, double y, double Ly, double leftCoord, double rightCoord, double bottomCoord, double topCoord)
// {
//   bool fluid = true;
//   if (y > Ly)
//   {
//     fluid = false;
//   }
//   if (x >= leftCoord && x <= rightCoord && y >= bottomCoord && y <= topCoord)
//   {
//     fluid = false;
//   }
//   return fluid;
// }

event fieldOutput(t = 2; t <= simTime; t += 2)
{
  // ignore I.C.
  char name1[24];
  char name2[24];
  char name3[24];
  sprintf(name1, "out-h-%.0f", t);
  sprintf(name2, "out-u-%.0f", t);
  sprintf(name3, "out-h-zoomIn-%.0f", t);
  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");
  FILE *fp3 = fopen(name3, "w");

  // compact form
  foreach ()
  {
    fprintf(fp1, "%g %g %g\n", x, y, h[]);
    fprintf(fp2, "%g %g\n", u.x[], u.y[]);
    if (x>150)
    {
      fprintf(fp3, "%g %g %g\n", x, y, h[]);
    }
  }

  // for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  //   fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}

// adaptivity
event adapt(i++)
{
  astats s = adapt_wavelet({h}, (double[]){normalDepth / 150.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
  //fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
