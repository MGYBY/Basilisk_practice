/* Gray & Edwards equations in dimensionless form. */
// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "./my-saint-venant.h"
/**The problem is solved in one dimension, but can be extended to two.

Problem-specific parameters consistent with Table 1 and Table 2 of Gray & Edwards (2014).

Perturbation frequency follow Edwards & Gray(2015).

*Note that the parameters are in SI units.*
*/
#define MU1 0.381863
#define MU2 0.643468
#define MU 0.554309 // Tan(Theta)
#define sinTheta pow((MU*MU/(1.0+MU*MU)), 0.50)
#define cosTheta pow((1.0/(1.0+MU*MU)), 0.50)
#define BETA 0.136
#define NORMALDEPTH 0.00319914
#define NORMALVEL 0.16899
#define ell 0.000825
#define FR 1.02
#define GAMMA (BETA*NORMALDEPTH/(ell*FR)) // 0.517 in G&E
#define gravCoeff 9.81
#define NUE ((2.0/9.0)*(ell*pow(gravCoeff, 0.50)/BETA)*(sinTheta/pow(cosTheta, 0.50))*((MU2-MU)/(MU-MU1)))
#define RE ((NORMALVEL*pow(NORMALDEPTH, 0.50))/NUE)
// #define RE 8.45
#define xextent_ 30.710
#define OUTPUTTIME 20.00
#define simTime 1200.0

#define disMag 0.20

#define MAXLEVEL 10
#define MINLEVEL 2
#define INITLEVEL 9

#define uThreshold 1.0E-10

scalar depthGrad[];

/**
Main and parameters
*/
int main()
{
/**
  The domain is 10 long (follow Edwards & Gray(2015)'s granular wave sim),

  tantheta is $\tan \theta$,

  Initial shape and velocity */
  X0 = 0.;
  L0 = xextent_;
  init_grid (1 << INITLEVEL);

  G = 1.0/(FR*FR);

  //use CFL number for $dt$ control
  CFL = 0.35;
  theta = 1.50; // use Sweby to slightly increase resolution

/**
If viscosity $\nu_e$  is not zero, the maximal time step is defined due to this viscosity:

*$dt$ will be updated later??*
*/
  DT = (L0/pow(2.0,MAXLEVEL))*(L0/pow(2.0,MAXLEVEL))/2.0/(1.0/RE)/2.0;

  // note that we are using the hump periodic configuration
  periodic (right);

  fprintf (ferr,"u0=%lf  h0=%lf nu=%lf Re=%lf \n",NORMALVEL,NORMALDEPTH,NUE,RE);

  run();
}

/**
The initial conditions */
event init (i = 0){

/**
Initial condition: steady-uniform flow
*/
   foreach(){
    h[] = 1.0+disMag*sin(2.0*M_PI*x/(xextent_));
    u.x[]= 1.0;

    depthGrad[] = fabs(h[-1]-h[1])/Delta;
    }
    boundary ({u.x,h});
}

event maxdt (t <= simTime; t += 0.99);

/** bed-slope and bed-friction: use a backward linearized Euler scheme for now
*/
event coulomb_friction (i++) {
//   double Froude, ff;
/**
We use a simple implicit scheme to implement coulomb bottom friction i.e.
(note the simplification by $h$)
$$\frac{d\mathbf{u}}{dt} = -\mu g \frac{\mathbf{u}}{|\mathbf{u}|}$$
  with $\mu$ fonction of $I$.
  Note that the good implementation preserving equilibrium balance is in Bouchut's book
*/
  foreach() {
//     Froude = u.x[]/sqrt(G*cos(zeta)*h[]);
    // ff = norm(u) > 0 ? (1. +  dt *(tan(zeta1)+(tan(zeta2)-tan(zeta1))/(beta*h[]/(ell*Froude)+1.0))*G*cos(zeta)) : HUGE ;

//   foreach_dimension()
      // u.x[] /= ff;
//       u.x[] = (dt*MU + FR*FR*u.x[])*(u.x[]+pow(h[], 3.0/2.0)*GAMMA)*norm(u)/(dt*((MU2-MU1)*u.x[]+MU1*(u.x[]+pow(h[], 3.0/2.0)*GAMMA))+FR*FR*(u.x[]+pow(h[], 3.0/2.0)*GAMMA)*norm(u));
    // rk3tvd
    double uMed = u.x[]>=uThreshold ? u.x[]+(1/pow(FR,2.0))*dt*(MU-(MU1+(MU2-MU1)/(1.0+GAMMA*pow(h[],3.0/2.0)/u.x[]))) : u.x[]+(1/pow(FR,2.0))*dt*(MU-MU1);
    uMed = u.x[]>=uThreshold ? (3./4.)*u.x[] + (1./4.)*uMed + (1/pow(FR,2.0))*(1./4.)*dt*(MU-(MU1+(MU2-MU1)/(1.0+GAMMA*pow(h[],3.0/2.0)/uMed))) : (3./4.)*u.x[] + (1./4.)*uMed + (1/pow(FR,2.0))*(1./ 4.)*dt*(MU-MU1);
    u.x[] = u.x[]>=uThreshold ? (1./3.)*u.x[] + (2./3.)*uMed + (1/pow(FR,2.0))*(2./3.)*dt*(MU-(MU1+(MU2-MU1)/(1.0+GAMMA*pow(h[],3.0/2.0)/uMed))) : (1./3.)*u.x[] + (1/pow(FR,2.0))*(2./3.)*uMed + (2./3.)*dt*(MU-MU1);

    depthGrad[] = fabs(h[-1]-h[1])/Delta;
  }
  boundary ({u.x,h});
}
/**
 The new viscous term from Gray & Edwards
$$ \frac{\partial }{\partial x}(\nu_e h^{3/2}\frac{\partial }{\partial x}u)$$
*/
event friction_long (i++) {
 foreach_dimension() {
    face vector g[];
//     scalar a = u.x;
    foreach_face()
      g.x[] = (1.0/RE)*(u.x[] - u.x[-1,0])/Delta*pow((h[0,0] + h[-1,0])/2.,3./2.);
      // central difference
      // g.x[] = (1.0/RE)*(u.x[1,0] - u.x[-1,0])/(2.0*Delta)*pow(h[], 3.0/2.0);
    foreach ()
      u.x[] += dt/Delta*(g.x[1,0] - g.x[] + g.y[0,1] - g.y[]);
      // central difference
      // u.x[] += dt/(2.0*Delta)*(g.x[1,0] - g.x[-1,0]);
  }
  boundary ((scalar *){u});
}

/**
# Output Control
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

event gnuplot(t = 0; t <= simTime; t += OUTPUTTIME)
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

event output(t = 0; t <= simTime; t += OUTPUTTIME)
{
     char name[80];
     sprintf(name, "out-%.1f", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g \n", x, h[], u.x[]);
     fprintf(fp, "\n");
     fclose(fp);
     fprintf(stderr, "%g\n", dt);
}

event outputMax(t = 0; t <= simTime; i += 40)
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

  fprintf(fp2, "%g %g %g %g\n", t, maxDepthLocX, maxDepth, maxDepthVel);

  fclose(fp2);
}

// AMR features
event adapt1 (i++) {
  astats s = adapt_wavelet({h, depthGrad}, (double[]){1/300.10, 0.0025}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  refine(x<=10.0*L0/pow(2, MAXLEVEL) && level<MAXLEVEL);
  refine(x>=L0-10.0*L0/pow(2, MAXLEVEL) && level<MAXLEVEL);
}
