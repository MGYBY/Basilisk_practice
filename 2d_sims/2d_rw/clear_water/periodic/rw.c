// #include "grid/multigrid.h"
//#include "grid/cartesian.h"
// #include "green-naghdi.h"
// #include "saint-venant.h"
#include "./saint-venant_hllc.h"
#include "utils.h"

#define MAXLEVEL 12
#define MINLEVEL 5
#define MAXMAXLEVEL 15

// problem-sepcific parameters
#define So 0.05011
#define frn 3.71
#define normalDepth 0.00798
#define gravityCoeff 9.81
#define gp (gravityCoeff*pow((1.0/(So*So+1.0)), 0.50))
#define normalVelocity (frn*pow(gp*normalDepth, 0.50))

#define disMag 0.20
#define disPeriod 0.933
#define disWavelength 1.272

#define simTime 128.0
#define OUTPUTINTERVAL 1.00

#define Lx 162.0
#define distLength (50.0*(Lx/pow(2, MAXMAXLEVEL)))
#define Ly (distLength*32.0)
// #define Ly (distLength*1.0)
#define cf (gp * So * 2. * normalDepth / (normalVelocity*normalVelocity)) // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

#define hemax (normalDepth/500.0)
#define uemax (normalVelocity/500.0)
#define hmaxemax (hemax/2.0)

// double FrNum=0.0;

// homogeneous bc
// bid block;
//free-slip bc
//not working???
// u.n[block] = dirichlet(0);

scalar hmax;

int main()
{
  size(Lx);
  G = gp;
  init_grid(1 << MAXLEVEL);
  //N = 1000;
  CFL = 0.399; // CFL number should be sufficiently small
//   theta = 1.5; // use Sweby limiter
  // N = 1024;

  periodic (right);

  theta = 1.30;

  /**
  We turn off limiting to try to reduce wave dissipation. */
  // gradient = NULL;

  run();
}

// h[left] = dirichlet((t<=disPeriod/2.0 && y<=distLength) ? normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod) : normalDepth);
// u.n[left] = dirichlet((t<=disPeriod/2.0 && y<=distLength) ? frn*sqrt(gp*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod))) : normalVelocity);
// u.t[left] = dirichlet(0.0);

// h[left] = ((t<=disPeriod/2.0 && y<=distLength) ? normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod) : normalDepth);
// u.n[left] = ((t<=disPeriod/2.0 && y<=distLength) ? frn*sqrt(gp*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod))) : normalVelocity);
// u.t[left] = (0.0);
//
// u.n[right] = neumann(0.);
// h[right] = neumann(0.);

// u.n[bottom] = neumann(0.);
// h[bottom] = neumann(0.);

event init(i = 0)
{
  refine(y>=(Ly-1.1*Lx/pow(2, MAXLEVEL)) && y<=(Ly+1.1*Lx/pow(2, MAXLEVEL)) && level<MAXMAXLEVEL);
  refine(y<=(distLength*2.0) && x<=(disWavelength*1.2) && x>=(disWavelength/2.0/1.2) && level<MAXMAXLEVEL);

  // for 2-D case 40x3m
  mask(y > Ly ? top : none);
  // block bc
//   mask(x >= leftCoord ? (x <= rightCoord ? (y >= bottomCoord ? (y <= topCoord ? block : none) : none) : none) : none);

//   h[left] = ((t<=disPeriod/2.0 && y<=distLength) ? normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod) : normalDepth);
//   u.n[left] = ((t<=disPeriod/2.0 && y<=distLength) ? frn*sqrt(gp*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod))) : normalVelocity);
//   u.t[left] = (0.0);
//
//   u.n[right] = neumann(0.);
//   h[right] = neumann(0.);

  u.n[top] = neumann(0.);
  h[top] = neumann(0.);

  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = neumann(0.);
  h[bottom] = neumann(0.);
  
//   refine(((x >= 0.0) && (x <= rightCoord+0.50)) && level < 14);

  // L0 = 40.;

  foreach ()
  {
    zb[] = 0.0;
    h[] = ((x>=disWavelength || x<=disWavelength/2.0) || y>distLength) ? normalDepth : normalDepth*(1.0+disMag*sin(2. * pi * (x-disWavelength/2.0) / disWavelength));
    u.x[] = ((x>=disWavelength || x<=disWavelength/2.0) || y>distLength) ? normalVelocity : frn*sqrt(gp*normalDepth*(1.0+disMag*sin(2. * pi * (x-disWavelength/2.0) / disWavelength)));
    u.y[] = 0.;

    hmax[] = h[];
  }

//   foreach_boundary (block){
//     if (x<blockXLower)
//     {
//       blockXLower = x;
//     }
//
//     if (x>blockXUpper)
//     {
//       blockXUpper = x;
//     }
//   }
}

/**
To implement an absorbing boundary condition, we add an area for $x >
12$ for which quadratic friction increases linearly with $x$. */
// only use Euler backward time int for now

event friction(i++)
{
  double uMed;
  
  foreach ()
  {
    // double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
//       double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
//       u.x[] = (u.x[] + G * So * dt) / a;
//       u.y[] /= a;
//       uMed = u.x[] + dt *  (-(cf / 2.) * u.x[] * norm(u) / h[] + gp*So);
//       uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * uMed * sqrt(sq(uMed)+sq(u.y[])) / h[] + gp*So);
//       u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * uMed * sqrt(sq(uMed)+sq(u.y[])) / h[] + gp*So);
//
//       uMed = u.y[] + dt *  (-(cf / 2.) * u.y[] * norm(u) / h[]);
//       uMed = (3. / 4.) * u.y[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * uMed * sqrt(sq(uMed)+sq(u.x[])) / h[]);
//       u.y[] = (1. / 3.) * u.y[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * uMed * sqrt(sq(uMed)+sq(u.x[])) / h[]);

      uMed = u.x[] + dt *  (-(cf / 2.) * u.x[] * norm(u) / h[] + gp*So);
      u.x[] = (1. / 2.) * u.x[] + (1. / 2.) * uMed + (1. / 2.) * dt * (-(cf / 2.) * uMed * sqrt(sq(uMed)+sq(u.y[])) / h[] + gp*So);

      uMed = u.y[] + dt *  (-(cf / 2.) * u.y[] * norm(u) / h[]);
      u.y[] = (1. / 2.) * u.y[] + (1. / 2.) * uMed + (1. / 2.) * dt * (-(cf / 2.) * uMed * sqrt(sq(uMed)+sq(u.x[])) / h[]);

      if (h[]>hmax[]) {
        hmax[]=h[];
      }
  }
}

// record max depth
event hmax(i+=25)
{
     double maxDepth = 0.0;
     double maxVel = 0.0;
//      double maxDepthLocX = 0.0;
//      double maxDepthVel = 0.0;
     FILE *fp3 = fopen("maxDepthVel", "a+");

//      maxDepth = h.max;
//      maxVel = ???;

     foreach ()
     {
               if (h[]>maxDepth){
                    maxDepth = h[];
               }

               if (norm(u)>maxVel){
                    maxVel = norm(u);
               }
     }

//      fprintf(fp3, "%g %g %g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel);
    fprintf(fp3, "%g %g %g \n", t, maxDepth, maxVel );
    fclose(fp3);
}

// event outputDrag(i++)
// {
//   // char name1[80];
//   // drag force in x dir
//   double fd = 0.0;
//   double rmax = 0.0;
// //   double maxDepth = 0.0;
// //   double maxDepthLocX = 0.0;
// //   double maxDepthVel = 0.0;
//   // name1 = "dragForce";
//   FILE *fp1 = fopen("dragForce", "a+");
//   FILE *fp2 = fopen("rMax", "a+");
//   FILE *fp3 = fopen("gauge10m", "a+");
//   FILE *fp4 = fopen("gauge20m", "a+");
//   FILE *fp5 = fopen("gauge40m", "a+");
//   FILE *fp6 = fopen("gauge80m", "a+");
//   FILE *fp7 = fopen("gauge120m", "a+");
//   FILE *fp8 = fopen("gauge140m", "a+");
//   FILE *fp9 = fopen("gauge160m", "a+");
//
//   foreach_boundary (block){
//     if (x==blockXLower)
//     {
//       fd += 0.5*rhoWater*gravityCoeff*(Delta*h[])*h[];
//     }
//     else if (x==blockXUpper)
//     {
//       fd -= 0.5*rhoWater*gravityCoeff*(Delta*h[])*h[];
//     }
//   }
//
// //     foreach ()
// //   {
// //       if (x<159.4){
// //           if (h[]>maxDepth){
// //              maxDepth = h[];
// //              maxDepthLocX = x;
// //              maxDepthVel = u.x[];
// //           }
// //       }
// //   }
//   /**
//   // for debugging
//   // foreach_boundary (block){
//   //   fprintf(fp1, "%g", h[]);
//   // }
//   {
//     fprintf(fp1, "%g %g", blockXLower, blockXUpper);
//   }
//   // for debugging
//   **/
//
//   fprintf(fp1, "%g %g %g\n", t, fd, fd/(0.5*rhoWater*(topCoord-bottomCoord)*sq(normalVelocity)*normalDepth));
//   rmax = interpolate(h, leftCoord-Lx/(pow(2, 14)), Ly/2.0);
//   fprintf(fp2, "%g %g %g\n", t, rmax, rmax/normalDepth);
//   fprintf(fp3, "%g %g %g\n", t, interpolate(h, 10.0, Ly / 2.), interpolate(u.x, 10.0, Ly / 2.));
//   fprintf(fp4, "%g %g %g\n", t, interpolate(h, 20.0, Ly / 2.), interpolate(u.x, 20.0, Ly / 2.));
//   fprintf(fp5, "%g %g %g\n", t, interpolate(h, 40.0, Ly / 2.), interpolate(u.x, 40.0, Ly / 2.));
//   fprintf(fp6, "%g %g %g\n", t, interpolate(h, 80.0, Ly / 2.), interpolate(u.x, 80.0, Ly / 2.));
//   fprintf(fp7, "%g %g %g\n", t, interpolate(h, 120.0, Ly / 2.), interpolate(u.x, 120.0, Ly / 2.));
//   fprintf(fp8, "%g %g %g\n", t, interpolate(h, 140.0, Ly / 2.), interpolate(u.x, 140.0, Ly / 2.));
//   fprintf(fp9, "%g %g %g\n", t, interpolate(h, 159.2, Ly / 2.), interpolate(u.x, 159.2, Ly / 2.));
//
//   fclose(fp1);
//   fclose(fp2);
//   fclose(fp3);
//   fclose(fp4);
//   fclose(fp5);
//   fclose(fp6);
//   fclose(fp7);
//   fclose(fp8);
//   fclose(fp9);
// }

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

event fieldOutput(t = 0; t <= simTime; t += OUTPUTINTERVAL)
{
  // ignore I.C.
//   char name1[24];
  char name2[40];
//   char name3[24];
//   char name4[32];
//   sprintf(name1, "out-h-%.0f", t);
  sprintf(name2, "snapshot-%g.gfs", t);
//   FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");

  // compact form
//   foreach ()
//   {
// //     fprintf(fp1, "%g %g %g\n", x, y, h[]);
// //     fprintf(fp2, "%g %g\n", u.x[], u.y[]);
//     if (x>150.0)
//     {
//       fprintf(fp3, "%g %g %g\n", x, y, h[]);
//       fprintf(fp2, "%g %g %g %g\n", x, y, u.x[], u.y[]);
//     }
//   }
//
//   for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
//   { if (interpolate(h, x, Ly / 2.)<10.0)
//       {
//       fprintf(fp4, "%g %g\n", x, interpolate(h, x, Ly / 2.));
//       }
//
//   }

  output_gfs(fp2, translate = true, list = {h, u.x, u.y});

//   sprintf (name2, "dumpSnapshot-%g", t);
//   dump(file=name2);

  // for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  //   fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

//   fclose(fp1);
  fclose(fp2);

}

// adaptivity
event adapt(i++)
{
  astats s = adapt_wavelet({h, hmax, u}, (double[]){hemax, hmaxemax, uemax, uemax, uemax}, maxlevel = MAXMAXLEVEL, minlevel = MINLEVEL);
  //fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

  refine(t<=disPeriod*2.0 && y<=(disWavelength*1.6) && x<=(disWavelength*2.5) && level<MAXMAXLEVEL);
}
