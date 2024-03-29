// normal flow simulation of the open channel
// use multilayer model only

#include "grid/multigrid1D.h"
#include "saint-venant.h"

// problem-sepcific parameters
#define So 0.05011
#define normalDepth 0.00798
#define normalVelocity  1.0311
#define gravityCoeff 9.81
#define disMag 0.05
#define disPeriod 0.933
#define simTime 40.0
double cf = 0.0; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

#define domainLength (10.0*normalDepth/So)

/**
The basin needs to be long enough so as to minimise the influence of
wave reflection at the outlet. Relatively high resolution is needed to
capture the dynamics properly. */

int main()
{
     N = 800;
//      L0 = 40.;
     L0 = domainLength;
     G = gravityCoeff;
     periodic (right);
//      nl = 1;
//      nu = 0.;
     CFL = 0.40; // CFL number should be sufficiently small
     // follow "$example/tsunami.c" method

     // h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
     // u.n[left] = dirichlet(normalVelocity);

     // u.n[right] = neumann(0);
     // h[right] = neumann(0);
     // u.n[right] = radiation(normalVelocity);
     // h[right] = radiation(normalDepth);
     // breaking = 0.1;
     run();
}

/**
We use ["radiation"
conditions](/src/elevation.h#radiation-boundary-conditions) at the
inlet and outlet. At the inlet (on the left), we try to impose the
desired sinusoidal wave form.*/

event init(i = 0)
{
     foreach ()
     {
//           zb[] = -So * x;
          zb[] = 0.;
          h[] = normalDepth*(1.0+0.02*sin(2.0*pi*x/domainLength));
          u.x[] = normalVelocity;
     }
}

static double chezyBedFriction(double u, double h, double cf)
{
     double rhs;
     rhs = gravityCoeff*So-(cf / 2.) * u * fabs(u) / h;
     return rhs;
}

// Quadratic bottom friction
event friction(i++)
{
     double uMed;
     cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
     foreach ()
     {
          // double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
          // double a = 1. + (cf / (2.)) * dt * u.x[] / h[];
          // foreach_dimension()
          // to be explicit!

          // in an inelegant way
          //     u.x[] /= a;
          uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=40$. 

![Snapshot of waves. The top of the bar is seen in white.](bar/snapshot.png)
*/

void plot_profile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't%g.png'\n"
             "set title 't = %.2f'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
               fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplot(t = 0; t <= simTime; t += 0.20)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
     // fprintf(fp,
     //         "set term pngcairo enhanced size 800,600 font \",10\"\n"
     //         "set output 't%.0f.png'\n"
     //         "set title 't = %.2f'\n"
     //         "set xrange [0:40]\n"
     //         "plot u 1:2 w l t\n",
     //         t, t);
     // fprintf(fp, "\n");
     // foreach ()
     //      fprintf(fp, "%g %g\n", x, h[]);
     // fprintf(fp, "e\n\n");
     // fflush(fp);
     // fprintf(stderr, "%.3f %.3f\n", t, statsf(h).max); // uncomment if needed
}

event output(t = 0; t <= simTime; t += 0.20)
{
     char name[80];
     sprintf(name, "out-%g", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g %g\n", x, eta[], u.x[], zb[]);
     fprintf(fp, "\n");
     fclose(fp);
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

//TODO: add AMR features
// event adapt (i++) {
// }
