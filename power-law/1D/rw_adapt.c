// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "saint-venant-power-law.h"

#define MAXLEVEL 12
#define MINLEVEL 10
#define MAXMAXLEVEL 17

double n_coeff = 0.40;
double alpha_coeff = 1.0;
// double lx = 2.878;
double disMag = 0.10;
// double betaCoeff = 0.0;
double simTime = 400.0;
double distPeriod = pi;

int main()
{
     L0 = 3600.0;
     N = 32768;
     alphaCoeff = alpha_coeff;
     betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
     // nCoeff = n_coeff;
     
     CFL = 0.25; // CFL number should be sufficiently small

     run();
}

event init(i = 0)
{
     h[left] = dirichlet( t<distPeriod ? 1.0 + disMag * 1.0 * sin(2. * pi * t / distPeriod) : 1.0);
     u.n[left] = dirichlet(1.0);

     u.n[right] = neumann(0.);
     h[right] = neumann(0.);

     foreach ()
     {
               zb[] = 0.0;
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               h[] = 1.0;
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

// bottom friction
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

          // rk3tvd
          if (h[] > dry) {
               uMed = u.x[] + dt * powerLawFriction(u.x[], h[], n_coeff);
               uMed = (3.0/4.0)*u.x[] + (1.0/4.0)*uMed + (1.0/4.0)*dt*powerLawFriction(uMed, h[], n_coeff);
               u.x[] = (1.0/3.0)*u.x[]+(2.0/3.0)*uMed+(2.0/3.0)*dt*powerLawFriction(uMed, h[], n_coeff);
          }
          else {
               u.x[] = 0.0;
          }

          //linearized backeard Euler
          // u.x[] = h[] > dry ? (u.x[] + dt)/(1.0+dt*(1.0/h[])*((pow(u.x[], (n_coeff-1.0)))/(pow(h[], n_coeff)))) : 0.0;
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

// record max depth
event hmax(i+=10)
{
     double maxDepth = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     FILE *fp2 = fopen("maxDepth", "a+");
     foreach ()
     {
               if (h[]>maxDepth){
                    maxDepth = h[];
                    maxDepthLocX = x;
                    maxDepthVel = u.x[];
               }
     }
     fprintf(fp2, "%g %g %g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel);
     fclose(fp2);
}

/**
Visualise the wave profile
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

event gnuplot(t = 0; t <= simTime; t += 4)
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

event output(t = 0; t <= simTime; t += 4)
{
     char name[80];
     sprintf(name, "out-%.0f.txt", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g \n", x, h[], u.x[]);
     fprintf(fp, "\n");
     fclose(fp);
}

// event output(i += 2; t <= 40)
//     output_gauges(gauges, {eta});

//TODO: add AMR features
event adapt (i++) {
     astats s = adapt_wavelet({h}, (double[]){1.0 / 150.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
}
