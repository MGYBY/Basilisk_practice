// normal flow simulation of the open channel
// use multilayer model only

// #include "grid/multigrid1D.h"
// #include "green-naghdi.h"
#include "grid/bitree.h"
#include "saint-venant.h"

#define MAXLEVEL 8
#define MINLEVEL 8
#define MAXMAXLEVEL 13

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.10;
double disPeriod = 0.933;
double simTime = 50.0;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
double Qin;

/**
The basin needs to be long enough so as to minimise the influence of
wave reflection at the outlet. Relatively high resolution is needed to
capture the dynamics properly. */

int main()
{
     N = 4096*2;
     L0 = 40.;
     G = gravityCoeff;
     CFL = 0.450; // CFL number should be sufficiently small
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
    Qin = normalDepth*normalVelocity;
     // for 2-D case
     // mask (y > 0.05 ? top : none);
//     h[left] = dirichlet(normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.));
    h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
//     h[left] = dirichlet((normalDepth + 1.0 * normalDepth));
     u.n[left] = dirichlet(Qin/((normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.))));

     // u.n[right] = neumann(0.);
     // // h[right] = neumann(0.);

     // u.n[right] = radiation(normalVelocity);
     // h[right] = radiation(normalDepth);
     u.n[right] = neumann(0.);
     h[right] = neumann(0.);
     
     cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

     foreach ()
     {
          // zb[] = -So * x;
          zb[] = 0.;
          h[] = normalDepth;
          u.x[] = normalVelocity;
     }
}

static double chezyBedFriction(double u, double h, double cf, double g, double So)
{
     double rhs;
     rhs = -(cf / 2.) * u * fabs(u) / h + g*So;
     return rhs;
}

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
  astats s = adapt_wavelet({h}, (double[]){normalDepth / 300.0}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
