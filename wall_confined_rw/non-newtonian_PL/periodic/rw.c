// #include "grid/cartesian1D.h"
#include "grid/bitree.h"
#include "predictor-corrector.h"
#include "utils.h"

// #define piVal 3.141592653589793 // not needed anymore
#define Q0 1.0
#define R0 1.0
#define disAmp 0.10

#define DOMAINLENGTH 64.0

// AMR parameters
#define INITLEVEL 9
#define MAXLEVEL 10
#define MINLEVEL 4

#define simTime 1000.0
#define outputInterval 20.0

#include "conservation.h"

/**
 * ## Variables
 */

scalar a[], q[];
scalar * scalars = {a,q};
vector * vectors = NULL;
scalar aGrad[], qGrad[];
char s[80];

/**
 *The other parameters are specific to the example. */

double nVal, FrVal, betaVal;

/**
 * ## Functions
 *
 * We define the *flux* function required by the [generic
 * solver](/src/conservation.h). */

void flux (const double * s, double * f, double * e)
{
  double a = s[0], q = s[1], u = q/a, radius = pow((a/pi),0.50);
  f[0] = q;
  f[1] = ((3.0*nVal+1.0)*q*q)/((2.0*nVal+1.0)*a)+1.0/3.0*betaVal*(1.0/pow(pi, 0.5))*pow(a, 1.5);
  // min/max eigenvalues
  double c = pow(4.0*nVal*(1.0+3.0*nVal)*pow(pi,0.5)*q*q+2.0*pow(a,2.50)*betaVal*pow((1.0+2.0*nVal),2.0),0.5);
  e[0] = (1.0/(2.0*a*(1.0+2.0*nVal)*pow(pi,0.25)))*(2.0*pow(pi,0.25)*q*(1.0+3.0*nVal) - c); // min
  e[1] = (1.0/(2.0*a*(1.0+2.0*nVal)*pow(pi,0.25)))*(2.0*pow(pi,0.25)*q*(1.0+3.0*nVal) + c); // max
}

/**
 * ## Parameters
 *
 * For small amplitudes $Amp = 0.01$ at the input boundary condition the
 * system has analytical solutions for $e1 < e2$, in this case the spatial
 * envelope of the flux rate behaves like $Q=Amp\times e^{-e2/2x}$ [Wang at
 * al., 2013]. */

int main() {
  L0 = DOMAINLENGTH;
  init_grid(1 << (INITLEVEL));

  nVal = 0.3;
  FrVal = 0.7;
  betaVal = 0.50;

  periodic (right);

  run();
}

/**
 * ## Initial conditions
 *
 * The initial conditions are $A=1$ and $Q=0$. */

event init (i = 0) {
  theta = 1.2; // tune limiting from the default minmod

  foreach(){
    a[] = pi*R0*R0*(1.0+disAmp*sin(2. * pi * x / DOMAINLENGTH));
    // simple velocity field for now
    q[] = Q0;
  }
}

/**
 * Source term.
 */
static double sourceTerm(double a, double q, double r, double n, double Fr)
{
  return (((1.0/Fr/Fr)*(1.0/pi)*a/pi)-(pow(q,n)/(pi*Fr*Fr*pow(r,(3.0*n-1.0)))));
}

// #if 0
event frictionTerm (i++) {
  double QMed, rv;
  foreach(){
    rv = pow(a[]/pi, 0.50);
    // for RK2TVD
    // RK2TVD
    QMed = q[] + dt*sourceTerm(a[], q[], rv, nVal, FrVal);

    q[] = (1. / 2.) *(q[]) + (1. / 2.) * QMed + (1. / 2.) * dt * sourceTerm(a[], QMed, rv, nVal, FrVal);
}
}

/**
 * ## Outputs
 *
 */
event printdata (t += outputInterval; t <= simTime) {
  sprintf (s, "out-%g.txt", t);
  FILE * fp2 = fopen (s, "w");
  foreach()
    fprintf (fp2, "%g %g %g \n", x, a[], q[]);
  fclose(fp2);
}

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
    fprintf(fp, "%g %g\n", x, a[]);
  }
  fprintf(fp, "e\n\n");
  fflush(fp);
}

event gnuplotOutput1(t = 0; t <= simTime; t += outputInterval)
{
  static FILE *fp = popen("gnuplot 2> /dev/null", "w");
  plot_profile(t, fp);
}

event adapt1 (i++) {
  vector ga[], gq[];
  gradients ({a}, {ga});
  gradients ({q}, {gq});
  foreach()
  {
    aGrad[] = norm(ga);
    qGrad[] = norm(gq);
  }

  adapt_wavelet({a, ga, gq}, (double[]){pi*R0*R0/500.0, 0.01, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);

    refine(x<=(17.510*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
    refine(x>=(DOMAINLENGTH-17.510*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
}

