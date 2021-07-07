// #include "grid/multigrid.h"
//#include "grid/cartesian.h"
// #include "green-naghdi.h"
#include "saint-venant.h"

#include "adapt_wavelet_limited.h"

#define MAXLEVEL 12
#define MINLEVEL 10
#define MAXMAXLEVEL 15

// problem-sepcific parameters
double So = 0.1192;
double normalDepth = 0.005333999999;
double normalVelocity = 1.28099972459;
double gravityCoeff = 9.81;
double disMag = 0.20;
double disPeriod = 0.694454984;
double simTime = 110.0;
double Lx = 162.0;
double Ly = 0.320;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
// block geometry
double leftCoord = 160.0;
double rightCoord = 160.05334;
double bottomCoord = 0.13333;
double topCoord = 0.18667;
double coordTol = 1e-7;
double FrNum = 0.0;
// drag force calc
double rhoWater = 1e3;

// for debugging
double blockXLower = 1e4;
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
  init_grid(1 << 12);
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
  refine (((x >= 0) && (x <= rightCoord+0.6)) && level < 14);
  refine(((x>leftCoord-0.06) && (x<rightCoord+0.06) && (y>bottomCoord-0.06) && (y<topCoord+0.06)) && level < 17);
  // block bc
  mask(x >= leftCoord ? (x <= rightCoord ? (y >= bottomCoord ? (y <= topCoord ? block : none) : none) : none) : none);

  // L0 = 40.;

  FrNum = normalVelocity/sqrt(gravityCoeff*normalDepth);
  // use only 2 periods to save computational cost?
  // h[left] = dirichlet(t>2*disPeriod? (normalDepth): normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
  h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
  u.n[left] = dirichlet(FrNum*sqrt(gravityCoeff*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod))));

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
  double rMax = 0.0;
  // name1 = "dragForce";
  FILE *fp1 = fopen("dragForce", "a+");
  FILE *fp2 = fopen("rMax", "a+");
//   FILE *fp3 = fopen("gauge10m", "a+");
//   FILE *fp4 = fopen("gauge20m", "a+");
//   FILE *fp5 = fopen("gauge40m", "a+");
//   FILE *fp6 = fopen("gauge80m", "a+");
//   FILE *fp7 = fopen("gauge120m", "a+");
//   FILE *fp8 = fopen("gauge140m", "a+");
//   FILE *fp9 = fopen("gauge160m", "a+");

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
  if (interpolate(h, leftCoord-Lx/pow(2, MAXMAXLEVEL), Ly/2.0)<10.0)
  {
  rMax = interpolate(h, leftCoord-Lx/pow(2, MAXMAXLEVEL), Ly/2.0);
  
  } else if(interpolate(h, leftCoord-Lx/pow(2, MAXMAXLEVEL)*2.0, Ly/2.0)<10.0){
  rMax = interpolate(h, leftCoord-Lx/pow(2, MAXMAXLEVEL)*2.0, Ly/2.0);
  }
  else{
  rMax = interpolate(h, leftCoord-Lx/pow(2, MAXLEVEL), Ly/2.0);
  }
  fprintf(fp2, "%g %g %g\n", t, rMax, rMax/normalDepth);
//   fprintf(fp3, "%g %g %g\n", t, interpolate(h, 10.0, Ly / 2.), interpolate(u.x, 10.0, Ly / 2.));
//   fprintf(fp4, "%g %g %g\n", t, interpolate(h, 20.0, Ly / 2.), interpolate(u.x, 20.0, Ly / 2.));
//   fprintf(fp5, "%g %g %g\n", t, interpolate(h, 40.0, Ly / 2.), interpolate(u.x, 40.0, Ly / 2.));
//   fprintf(fp6, "%g %g %g\n", t, interpolate(h, 80.0, Ly / 2.), interpolate(u.x, 80.0, Ly / 2.));
//   fprintf(fp7, "%g %g %g\n", t, interpolate(h, 120.0, Ly / 2.), interpolate(u.x, 120.0, Ly / 2.));
//   fprintf(fp8, "%g %g %g\n", t, interpolate(h, 140.0, Ly / 2.), interpolate(u.x, 140.0, Ly / 2.));
//   fprintf(fp9, "%g %g %g\n", t, interpolate(h, 159.2, Ly / 2.), interpolate(u.x, 159.2, Ly / 2.));

  fclose(fp1);
  fclose(fp2);
//   fclose(fp3);
//   fclose(fp4);
//   fclose(fp5);
//   fclose(fp6);
//   fclose(fp7);
//   fclose(fp8);
//   fclose(fp9);
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

event fieldOutput1(t = 0; t <= simTime; t += 6.0)
{
  // ignore I.C.
  char name1[28];
  char name2[28];
//   char name3[28];
  char name4[32];
  sprintf(name1, "out-h-local-%.1f", t);
  sprintf(name2, "out-u-local-%.1f", t);
//   sprintf(name3, "out-h-zoomIn-%.0f", t);
  sprintf(name4, "profile-%.1f", t);
  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");
//   FILE *fp3 = fopen(name3, "w");
  FILE *fp4 = fopen(name4, "w");

  // compact form
  foreach ()
  {
    if (x>152)
    {
      fprintf(fp1, "%g %g %g\n", x, y, h[]);
      fprintf(fp2, "%g %g %g %g\n", x, y, u.x[], u.y[]);
    }
  }
  
  for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  {
      if (interpolate(h, x, Ly / 2.)< 10.0)
      {
    fprintf(fp4, "%g %g\n", x, interpolate(h, x, Ly / 2.));
      }
  }

  // for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  //   fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

  fclose(fp1);
  fclose(fp2);
//   fclose(fp3);
  fclose(fp4);
}

event fieldOutput2(t = 60.55; t <= 68.60; t += 0.025)
{
  // ignore I.C.
  char name1[28];
  char name2[28];
//   char name3[28];
  char name4[32];
  sprintf(name1, "out-h-local-%.3f", t);
  sprintf(name2, "out-u-local-%.3f", t);
//   sprintf(name3, "out-h-zoomIn-%.0f", t);
  sprintf(name4, "profile-%.3f", t);
  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");
//   FILE *fp3 = fopen(name3, "w");
  FILE *fp4 = fopen(name4, "w");

  // compact form
  foreach ()
  {
    if (x>152)
    {
      fprintf(fp1, "%g %g %g\n", x, y, h[]);
      fprintf(fp2, "%g %g %g %g\n", x, y, u.x[], u.y[]);
    }
  }
  
  for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  {
      if (interpolate(h, x, Ly / 2.)< 10.0)
      {
    fprintf(fp4, "%g %g\n", x, interpolate(h, x, Ly / 2.));
      }
  }

  // for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  //   fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

  fclose(fp1);
  fclose(fp2);
//   fclose(fp3);
  fclose(fp4);
}

int refRegion(double x,double y){
    int lev;
    if( (x>leftCoord-0.4) && (x<rightCoord+0.1) )
      lev = 17;
    else
      lev = 15;

    return lev;
}

// adaptivity
event adapt(i++)
{
  astats s = adapt_wavelet_limited({h}, (double[]){normalDepth / 150.0},  refRegion, minlevel = 12);
  //fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
