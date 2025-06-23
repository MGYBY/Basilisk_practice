// normal flow simulation of the open channel
// use multilayer model only

#include "grid/multigrid1D.h"
// #include "green-naghdi.h"
// #include "grid/bitree.h"
#include "saint-venant.h"

#define MAXLEVEL 8
#define MINLEVEL 8
#define MAXMAXLEVEL 13

// problem-sepcific parameters
#define sinTheta (0.15/pow(0.15*0.15+0.3*0.3, 0.50)) // 30 degree
#define cosTheta (pow(1.0-sinTheta*sinTheta, 0.50))
#define tanTheta (sinTheta/cosTheta)
#define simTime 50.0
#define cf 1.0e-3
#define hin 0.30
#define uin 2.00

#define STAIRNUM 20
#define STAIRHEIGHT 0.150
#define FRONTLEN 2.80
#define WIDTH 0.25
#define BACKLEN 1.50
#define BOTTOMELEV ((-1.0)*STAIRNUM*STAIRHEIGHT)

#define PLOTTIME 0.25
#define TEXTTIME 1.0

#define GRAV 9.81

// u.n[right] = neumann(0.);
// h[right] = neumann(0.);

/**
The basin needs to be long enough so as to minimise the influence of
wave reflection at the outlet. Relatively high resolution is needed to
capture the dynamics properly. */

int main()
{
     N = 512;
     L0 = (STAIRNUM*STAIRHEIGHT)/tanTheta+FRONTLEN+BACKLEN;
     G = GRAV;
     CFL = 0.450; // CFL number should be sufficiently small
//      theta = 2.0;
//      gradient = NULL;
     run();
}

event init(i = 0)
{
    int i;
    // fixed-water-head BC
    h[left] = dirichlet(hin);
//     u.n[left] = dirichlet(uin);
    u.n[left] = neumann(0.);
     u.n[right] = neumann(0.);
     h[right] = neumann(0.);

     foreach ()
     {
          // zb[] = -So * x;
       zb[] = 0.0;
       for (i = 0; i < STAIRNUM; i++){
         if (x>=i*STAIRHEIGHT/tanTheta+FRONTLEN && x<(i+1)*STAIRHEIGHT/tanTheta+FRONTLEN){
          zb[] = (-1.0)*(i+1)*STAIRHEIGHT;
         }
         if (x>8.8) {
          zb[] = BOTTOMELEV;
         }
       }
       // "dam-break" IC
          h[] = x<=FRONTLEN ? hin : 0.0;
          u.x[] = 0.0;
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
         double a = h[] < dry ? HUGE : 1. + cf*dt*norm(u)/h[];
         foreach_dimension()
            u.x[] /= a;
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
             "set title 't = %g'\n"
             "set xlabel 'x(m)'\n"
             "set ylabel 'h (m)'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
          fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplot(t = 0; t <= simTime; t += PLOTTIME)
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

event output(t = 0; t <= simTime; t += TEXTTIME)
{
     char name[80];
     sprintf(name, "out-%.0f", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g %g %g \n", x, h[], u.x[], eta[], zb[]);
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
// event adapt (i++) {
//   astats s = adapt_wavelet({h, u.x}, (double[]){hin / 300.0, uin/300}, maxlevel = MAXMAXLEVEL, minlevel = MAXLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
// }
