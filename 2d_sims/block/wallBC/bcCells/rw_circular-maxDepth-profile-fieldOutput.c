#include "saint-venant.h"

#define MAXLEVEL 12
#define MINLEVEL 10
#define MAXMAXLEVEL 14
#define IN_CYLINDER (sq(x-20.20) + sq(y-1.50) < sq(0.2))

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.20;
double disPeriod = 0.933;
double simTime = 120.0;
double Lx = 162.0;
double Ly = 3.0;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
// block geometry
double cylinderCenterX = 160.20;
double cylinderCenterY = 1.50;
double cylinderRadius = 0.2;
double coordTol = 1e-7;
// drag force calc
double rhoWater = 1e3;

// homogeneous bc
bid cylinder;
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
  // N = 1024;

  /**
  We turn off limiting to try to reduce wave dissipation. */
  // gradient = NULL;

  run();
}

event init(i = 0)
{
    double FrNum = 0.0;
//   scalar region[];
//   foreach(){
//     region[] = IN_CYLINDER*noise()
//   }
  refine (sq(x-cylinderCenterX) + sq(y-cylinderCenterY) < sq(cylinderRadius*1.10) && level < MAXMAXLEVEL);
  // for 2-D case 40x3m
  mask(y > Ly ? top : none);
  // circular bc
  mask(sq(x - cylinderCenterX) + sq(y - cylinderCenterY) <= sq(cylinderRadius) ? cylinder : none);

  // L0 = 40.;

  FrNum = normalVelocity/sqrt(gravityCoeff*normalDepth);
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

// event refine_boundary(i++){
//     mask(sq(x - cylinderCenterX) + sq(y - cylinderCenterY) <= sq(cylinderRadius) ? cylinder : none);
// }


event outputDrag(i++)
{
  // char name1[80];
  // drag force in x dir
  double fd = 0.0;
  double maxDepth = 0.0;
  double maxDepthLocX = 0.0;
  double maxDepthVel = 0.0;
  // double k1 = 0.0;
  // name1 = "dragForce";
  FILE *fp1 = fopen("dragForce", "a+");
  FILE *fp2 = fopen("maxDepth", "a+");

  foreach_boundary(cylinder)
  {
    // coefficient to conver to x-component
    // k1 = (cylinderCenterX - x) / cylinderRadius;
    fd += ((cylinderCenterX - x) / cylinderRadius) * 0.5 * rhoWater * gravityCoeff * (Delta * h[]) * h[];
  }
  
  foreach ()
  {
      if (x<159.4){
          if (h[]>maxDepth){
             maxDepth = h[];
             maxDepthLocX = x;
             maxDepthVel = u.x[];
          }
      }
  }

  /**
  // for debugging
  // foreach_boundary (block){
  //   fprintf(fp1, "%g", h[]);
  // }
  // for debugging
  **/

  fprintf(fp1, "%g %g %g\n", t, fd, fd / (0.5 * rhoWater * (2 * cylinderRadius) * sq(normalVelocity) * normalDepth));
  fprintf(fp2, "%g %g %g\n", t, maxDepthLocX, maxDepth);

  fclose(fp1);
  fclose(fp2);
}

/**
At the end of the simulation, we output the maximum wave amplitudes
along the cross-sections corresponding with the experimental data. */

// static bool isFluid(double x, double y, double Ly, double cylinderCenterX, double cylinderCenterY, double cylinderRadius)
// {
//   bool fluid = true;
//   if (y > Ly)
//   {
//     fluid = false;
//   }
//   if (sq(x - cylinderCenterX) + sq(y - cylinderCenterY) <= sq(cylinderRadius))
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
  char name4[32];
  sprintf(name1, "out-h-%.0f", t);
  sprintf(name2, "out-u-%.0f", t);
  sprintf(name3, "out-h-zoomIn-%.0f", t);
  sprintf(name4, "profile-%.0f", t);
  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");
  FILE *fp3 = fopen(name3, "w");
  FILE *fp4 = fopen(name1, "w");

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
  
  for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
    fprintf(fp4, "%g %g\n", x, interpolate(h, x, Ly / 2.));

  // for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  //   fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
}

// adaptivity
event adapt(i++)
{
  astats s = adapt_wavelet({h}, (double[]){normalDepth / 120.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
  // fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
