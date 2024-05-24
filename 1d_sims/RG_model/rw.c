/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
// #include "grid/cartesian1D.h"
#include "grid/bitree.h"
// #include "grid/quadtree.h"
// #include "./saint-venant-cb.h"
#include "saint-venant_rg_new.h"

#define gravity 9.81
#define sinTheta 0.05011
#define cosTheta (pow((1.0-sinTheta*sinTheta),0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta)

#define froude 3.71
#define normalDepth 0.00800
#define normalVel (froude*pow(gravity*normalDepth,0.50)) // follow Brock's definition
#define cfVal (gravity*normalDepth*sinTheta/pow(normalVel,2.0))
#define crVal 0.00035

// dissipation in the roller
#define bigPhiVal 0.00
// dissipation at wall
#define smallPhiVal 22.76

#define piVal 3.14159265
#define disAmp 0.20

#define DOMAINLENGTH 1.356

#define INITLEVEL 10
#define MAXLEVEL 10
#define MINLEVEL 4

#define simTime 50.0
#define outputInterval 0.50

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

char s[80];
scalar hmax[], hGrad[], uxGrad[], teGrad[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  // FIXME: for 2D: G is gravity
  G = gPrime;
  init_grid(1 << (INITLEVEL));

  theta = 1.30;
  CFL = 0.40;

  periodic (right);

  run();
}

/**
## Initialization  */
event init (i = 0) {
  /**
  We initialize *h*. */
  foreach(){
    double totalDepth = normalDepth*(1.0+disAmp*sin(2. * pi * x / DOMAINLENGTH));
    zb[] = 0.0;
    h[] = totalDepth;
    // simple velocity field for now
    u.x[] = normalVel;
    te[] = 0.50*pow(u.x[],2.0)+0.50*(G*h[]+(bigPhiVal+smallPhiVal)*pow(h[],2.0));
    }
}

static double momXFric(double u, double v, double h, double te)
{
     double rhs;
    // FIXME: maybe this not general enough for 2D.
     rhs = h>dry ? (G*So - cfVal*u*pow((u*u+v*v),0.50)/h) : 0.0;
     return rhs;
}

static double momYFric(double u, double v, double h, double te)
{
     double rhs;
    // FIXME: maybe this not general enough for 2D.
     rhs = h>dry ? (-1.0)*cfVal*v*pow((u*u+v*v),0.50)/h : 0.0;
     return rhs;
}

static double teFric(double u, double v, double h, double te)
{
     double rhs, Ce, phiTotal, bigPhi;
     phiTotal = -1.0*((G*h-2.0*te+u*u)/pow(h,2.0));
     bigPhi = phiTotal-smallPhiVal;

     Ce = cfVal+crVal*bigPhi/phiTotal;
     rhs = h>dry ? -1.0*Ce*pow((u*u+v*v),(3.0/2.0))/h + G*So*u : 0.0;
     return rhs;
}

/**
  Source term.
*/
event frictionTerm (i++) {
  foreach(){
    double uMed, teMed;
//     double vPrev, vMed; // for 2D
//     uPrev = u.x[];
//     vPrev = u.y[];

    // RK3TVD
    // x-mom AND y-mom AND total-energy eqn source term
    uMed = u.x[] + dt *  momXFric(u.x[], 0.0, h[], te[]); // FIXME: maybe not general enough for 2D in v
//     vMed = u.y[] + dt *  momYFric(u.x[], u.y[], h[], te[]);
    teMed = te[] + dt *  teFric(u.x[], 0.0, h[], te[]); // FIXME: maybe not general enough for 2D in v

    uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * momXFric(uMed, 0.0, h[], teMed); // FIXME: maybe not general enough for 2D in v
//     vMed = (3. / 4.) * u.y[] + (1. / 4.) * vMed + (1. / 4.) * dt * momYFric(uMed, vMed, h[], te[]);
    teMed = (3. / 4.) * te[] + (1. / 4.) * teMed + (1. / 4.) * dt * teFric(uMed, 0.0, h[], teMed); // FIXME: maybe not general enough for 2D in v

    u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * momXFric(uMed, 0.0, h[], teMed); // FIXME: maybe not general enough for 2D in v
//     u.y[] = (1. / 3.) * u.y[] + (2. / 3.) * vMed + (2. / 3.) * dt * momYFric(uMed, vMed, h[], teMed);
    te[] = (1. / 3.) * te[] + (2. / 3.) * teMed + (2. / 3.) * dt * teFric(uMed, 0.0, h[], teMed);

    }

    boundary ((scalar *){u, te});
}

/*
 * AMR here
 *
 */
event adapt1 (i++) {
  vector gh[], gux[], gte[];
  gradients ({h}, {gh});
  gradients ({u.x}, {gux});
  gradients ({te}, {gte});
  foreach()
  {
    hGrad[] = norm(gh);
    uxGrad[] = norm(gux);
    teGrad[] = norm(gte);
  }

  adapt_wavelet({h, hGrad, uxGrad, teGrad}, (double[]){normalDepth/300.0, 0.004, normalVel/300.0, (normalVel*normalVel)/300.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval) {
  sprintf (s, "slice-%g.txt", t);
  FILE * fp2 = fopen (s, "w"); 
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g \n", x, h[], u.x[], te[]);
  }
  fclose(fp2);
}

/**
Visualise the wave profile
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

event gnuplotOutput1(t = 0; t <= simTime; t += outputInterval)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
}

// event moviemaker (t = end) {
//     system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
// 	    "ppm2mp4 movie.mp4");
// }
