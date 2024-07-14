// normal flow simulation of the open channel
// use multilayer model only

// #include "grid/multigrid1D.h"
// #include "green-naghdi.h"
#include "grid/bitree.h"
#include "./my-saint-venant.h"

#define MAXLEVEL 8
#define MINLEVEL 2
#define MAXMAXLEVEL 15

#define gradTol 2.50e-3
#define disPeriod 0.933
#define initRefineStage (disPeriod*16.0)

// problem-sepcific parameters
#define So 0.05011
#define sinTheta (So*1.0)
#define cosTheta (pow((1.0-sinTheta*sinTheta),0.50))
#define normalDepth 0.00798
#define normalVelocity 1.0377
#define gravityCoeff 9.81
#define gPrime (gravityCoeff*cosTheta)
#define froude (normalVelocity/sqrt(gPrime*normalDepth))
#define disMag 0.20
#define simTime 120.0
#define cf (gravityCoeff * So * 2. * normalDepth / (normalVelocity*normalVelocity))  // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
#define Qin (normalDepth*normalVelocity)

// concentration
#define concInit 1.0

/**
The basin needs to be long enough so as to minimise the influence of
wave reflection at the outlet. Relatively high resolution is needed to
capture the dynamics properly. */

int main()
{
     N = 4096;
     L0 = (1200.0*normalDepth/sinTheta);
     G = gPrime;
     CFL = 0.40; // CFL number should be sufficiently small
//      theta = 2.0;
//      gradient = NULL;
     run();
}

/**
We use ["radiation"
conditions](/src/elevation.h#radiation-boundary-conditions) at the
inlet and outlet. At the inlet (on the left), we try to impose the
desired sinusoidal wave form.*/

// h[left] = dirichlet(normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.));


event init(i = 0)
{
     // for 2-D case
     // mask (y > 0.05 ? top : none);
//     h[left] = dirichlet(normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.));
     h[left] = dirichlet(t<=disPeriod/2.0 ? normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod) : normalDepth);
//     h[left] = dirichlet((normalDepth + 1.0 * normalDepth));
     // u.n[left] = dirichlet(Qin/((normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.))));
     u.n[left] = dirichlet(t<=disPeriod/2.0 ? froude*sqrt(gPrime*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod))) : normalVelocity);
     conc[left] = dirichlet(concInit);

     // u.n[right] = neumann(0.);
     // // h[right] = neumann(0.);

     // u.n[right] = radiation(normalVelocity);
     // h[right] = radiation(normalDepth);
     u.n[right] = neumann(0.);
     h[right] = neumann(0.);
     conc[right] = neumann(0.);

     foreach ()
     {
          // zb[] = -So * x;
          zb[] = 0.;
          h[] = normalDepth;
          u.x[] = normalVelocity;
          conc[] = concInit;
     }
}

// static double chezyBedFriction(double u, double h, double cf, double g, double So)
// {
//      double rhs;
//      rhs = -(cf / 2.) * u * fabs(u) / h + g*So;
//      return rhs;
// }

// Quadratic bottom friction
event friction(i++)
{
//     double uMed;
     
     foreach ()
     {
//           double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
          // double a = 1. + (cf / (2.)) * dt * u.x[] / h[];
          foreach_dimension(){
//                   u.x[] /= a;
              u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * norm(u) / (h[]));
          
//           uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf, G, So);
//           uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf, G, So);
//           u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf, G, So);
              
        }
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

// record max depth
// event hmaxUmax (i+=50)
// {
//      double maxDepth = 0.0;
//      double maxVel = 0.0;
//      double maxDepthLocX = 0.0;
//      double maxDepthVel = 0.0;
//      double maxVelLocX = 0.0;
//      double maxVelDepth = 0.0;
//
//      FILE *fp2 = fopen("maxDepth", "a+");
//      FILE *fp3 = fopen("maxVel", "a+");
//
//      foreach ()
//      {
//            if (h[]>maxDepth){
//                 maxDepth = h[];
//                 maxDepthLocX = x;
//                 maxDepthVel = u.x[];
//            }
//
//            if (u.x[]>maxVel){
//                 maxVel = u.x[];
//                 maxVelLocX = x;
//                 maxVelDepth = h[];
//            }
//      }
//
//      fprintf(fp2, "%.10g %.10g %.10g %.10g \n", t, maxDepth, maxDepthLocX, maxDepthVel);
//      fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, maxVel, maxVelLocX, maxVelDepth);
//
//      fclose(fp2);
//      fclose(fp3);
// }
//
// event hFront1 (t=0; i+=8; t<initRefineStage)
// {
//      double aveDepth = 0.0;
//      double aveVel = 0.0;
//      double xf = 0.0;
//      FILE *fp3 = fopen("frontPosMod", "a+");
//      foreach ()
//      {
//           xf = h[] > dry ?  max(xf,x) :  xf ;
//           aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
//           aveVel += (h[]>=dry ? Delta*u.x[] : 0.0);
//      }
//      fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf) );
//      fclose(fp3);
// }
//
// event hFront2 (t=initRefineStage; t<=simTime; i+=80)
// {
//      double aveDepth = 0.0;
//      double aveVel = 0.0;
//      double xf = 0.0;
//      FILE *fp3 = fopen("frontPosMod", "a+");
//      foreach ()
//      {
//           xf = h[] > dry ?  max(xf,x) :  xf ;
//           aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
//           aveVel += (h[]>=dry ? Delta*u.x[] : 0.0);
//      }
//
//      fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf) );
//      fclose(fp3);
// }

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=40$. 

![Snapshot of waves. The top of the bar is seen in white.](bar/snapshot.png)
*/

void plot_hProfile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't-%g_h.png'\n"
             "set title 't = %.2f'\n"
             "set xlabel 'x(m)'\n"
             "set ylabel 'eta(m)'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
          fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

void plot_concProfile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't-%g_conc.png'\n"
             "set title 't = %.2f'\n"
             "set xlabel 'x(m)'\n"
             "set ylabel 'conc'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
          fprintf(fp, "%g %g\n", x, conc[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplotHeight(t = 0; t <= simTime; t += 0.5)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_hProfile(t, fp);
}

event gnuplotConc(t = 0; t <= simTime; t += 0.5)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_concProfile(t, fp);
}

event output(t = 0; t <= simTime; t += 0.5)
{
     char name[80];
     sprintf(name, "out-%.0f", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g %g \n", x, h[], u.x[], conc[]);
     fprintf(fp, "\n");
     fclose(fp);
     fprintf(stderr, "%g\n", dt);
}

/**
The location of the gauges is difficult to find in the litterature, we
used a combination of [Yamazaki et al,
2009](/src/references.bib#yamazaki2009) and [Dingemans,
1994](/src/references.bib#dingemans1994). */

// Gauge gauges[] = {
//     {"WG4", 10.5},
//     {"WG5", 12.5},
//     {"WG6", 13.5},
//     {"WG7", 14.5},
//     {"WG8", 15.7},
//     {"WG9", 17.3},
//     {"WG10", 19},
//     {"WG11", 21},
//     {NULL}};

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

// AMR features
event adapt (i++) {
  astats s = adapt_wavelet({h}, (double[]){normalDepth / 425.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
