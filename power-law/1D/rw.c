#include "grid/multigrid1D.h"
#include "saint-venant-power-law.h"



double n_coeff = 0.40;
double alpha_coeff = 1.0;
double lx = 2.878;
double disMag = 0.025;
// double betaCoeff = 0.0;
double simTime = 100.0;

int main()
{
     L0 = 2.878;
     // origin (-L0/2., -L0/2.);
     periodic (right);
     N = 512;
     alphaCoeff = alpha_coeff;
     betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
     // nCoeff = n_coeff;
     
     CFL = 0.40; // CFL number should be sufficiently small

     run();
}

event init(i = 0)
{
     foreach ()
     {
               zb[] = 0.0;
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               h[] = 1.0 + disMag * sin(2. * pi * x / lx);
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               u.x[] = 1.0;
     }
}

static double powerLawFriction(double u, double h, double n)
{
     double rhs;
     rhs = h-pow((u/h), n);
     return rhs;
}

// Quadratic bottom friction
event friction(i++)
{
     double uMed=0.0;
     foreach ()
     {
          // rk2tvd
          // if (h[] > dry) {
          //      uMed = u.x[] + dt * powerLawFriction(u.x[], h[], n_coeff);
          //      u.x[] = 0.5*u.x[]+0.5*uMed+0.5*dt*powerLawFriction(uMed, h[], n_coeff);
          // }
          // else {
          //      u.x[] = 0.0;
          // }

          //linearized backeard Euler
          u.x[] = h[] > dry ? (u.x[] + dt)/(1.0+dt*(1.0/h[])*((pow(u.x[], (n_coeff-1.0)))/(pow(h[], n_coeff)))) : 0.0;
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

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
             "set output 't%.0f.png'\n"
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

event gnuplot(t = 0; t <= simTime; t += 2)
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

event output(t = 0; t <= simTime; t += 2)
{
     char name[80];
     sprintf(name, "out-%.0f", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g \n", x, h[], u.x[]);
     fprintf(fp, "\n");
     fclose(fp);
}

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

//TODO: add AMR features
// event adapt (i++) {
// }
