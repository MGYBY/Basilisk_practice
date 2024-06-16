// #include "grid/cartesian1D.h"
#include "grid/bitree.h"
#include "predictor-corrector.h"
#include "utils.h"

#include "conservation.h"
// #include "./my-conservation.h"

#define gravity 9.81
#define sinTheta 0.05011
#define cosTheta (pow((1.0-sinTheta*sinTheta),0.50))
#define So (sinTheta/cosTheta)
#define gPrime (gravity*cosTheta) // the actual gravity on the slope
#define gHat (gravity*sinTheta)

#define froude 3.71
#define normalDepth 0.00798

// problem-specific parameters
#define kappa 0.412
#define nuVis (9.36*pow(10.0,-7.0))
// from the approximate function 4.9 and 4.10
#define RInf 2.0434
#define R1Inf 3.7826
#define alpha (R1Inf-RInf+1.0)
#define alpha2 2.47
#define alpha1 (alpha-alpha2)
#define crCoeff 0.48
// note that cf is actually a function of h here
#define cfNHelpCoeff (pow(gHat*normalDepth*normalDepth*normalDepth, 0.50))
#define cfNVal (kappa*kappa*pow((RInf-2.0+2.0*log(2.0)+log(kappa)+log(cfNHelpCoeff/nuVis)),-2.0))
// use this trick to maintain normal flow for cf
#define normalVel (pow((gHat*normalDepth/cfNVal), 0.5))

// dissipation in the roller
#define bigPhiNVal 0.00
// dissipation at wall
#define smallPhiNVal (gHat/normalDepth/pow(kappa,2.0))

#define piVal 3.14159265
#define disAmp 0.05
#define distPeriod 0.933
#define initRefineStage (distPeriod*12.0)

#define DOMAINLENGTH (1200.0*normalDepth/sinTheta)

#define INITLEVEL 12
#define MAXLEVEL 15
#define MINLEVEL 1

#define simTime 128.0
#define outputInterval 0.5

/**
## Variables

We define the conserved scalar fields $a$ and $q$ which are passed to the
generic solver through the *scalars* list. We don't have any conserved
vector field. */

scalar h[], q[], hte[], hsphi[];
scalar * scalars = {h,q,hte,hsphi};
vector * vectors = NULL;

scalar hGrad[], qGrad[], hteGrad[], hsphiGrad[], bigPhiGrad[];
char s[80];

/**
## Functions

We define the *flux* function required by the [generic
solver](/src/conservation.h). */

void flux (const double * s, double * f, double e[2])
{  
  double h = s[0], q = s[1], u = q/h, hte = s[2], te = hte/h, hsphi = s[3], sphi = hsphi/h;
  double pt = -1.0*((gPrime*h-2.0*te+u*u)/pow(h,2.0));
  double capPi = 0.50*gPrime*h*h+pt*pow(h, 3.);
  f[0] = q;
  f[1] = q*u + capPi;
  f[2] = q*te + capPi*u;
  f[3] = q*sphi;
  // min/max eigenvalues
  double c = sqrt(gPrime*h+3.*h*h*pt);
  e[0] = u - c; // min
  e[1] = u + c; // max
}

/**
## Parameters

For small amplitudes $Amp = 0.01$ at the input boundary condition the
system has analytical solutions for $e1 < e2$, in this case the spatial
envelope of the flux rate behaves like $Q=Amp\times e^{-e2/2x}$ [Wang at
al., 2013]. */

int main() {
  L0 = DOMAINLENGTH;
  init_grid(1 << (INITLEVEL));

  theta = 1.30;

  // periodic (right);

  fprintf (ferr, "U = %g \n", normalVel);
  fprintf (ferr, "cf = %g \n", cfNVal);
  fprintf (ferr, "Fr = %g \n", (normalVel/pow(gravity*cosTheta*normalDepth, 0.50)));

  run();
}

double distDepthInlet(double t)
{
//   "type-c" inlet dist
  if (t<=distPeriod/2.0)
    return (normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod));
  else
    return normalDepth;

  // "type-a" inlet dist
//   return (normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod));

  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

double distQInlet(double t)
{
  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
//   "type-c", const Fr inlet dist
  if (t<=distPeriod/2.0)
    return (froude*sqrt(gPrime*(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod)))*(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod)));
  else
    return (normalVel*normalDepth);

//   "type-c", const q inlet dist
//   if (t<=distPeriod/2.0)
//     return (normalDepth*normalVel/(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod)));
//   else
//     return normalVel;

  // "type-a" inlet dist
//   return (froude*sqrt(G*(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod))));

  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

double distHteInlet(double t)
{
  double hVar = 0.0;
  double vVar = 0.0;
  double phiVar = 0.0;
  hVar = normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod);
  phiVar= gHat/(kappa*kappa)/hVar;
  vVar = froude*sqrt(gPrime*hVar);

  if (t<=distPeriod/2.0)
    return (hVar*(pow(vVar,2.0)/2.0 + pow(hVar,2.0)*phiVar/2.0 + 0.0 + gPrime*hVar/2.0));
  else
    return (normalDepth*(pow(normalVel,2.0)/2.0 + pow(normalDepth,2.0)*(gHat/(kappa*kappa)/normalDepth)/2.0 + 0.0 + gPrime*normalDepth/2.0));

}

/**
## Initial conditions 

The initial conditions are $A=1$ and $Q=0$. */

event init (i = 0) {
  // "wavemaker" inlet BC
  h[left] = dirichlet(distDepthInlet(t));
  q[left] = dirichlet(distQInlet(t));
  hsphi[left] = dirichlet(t<=distPeriod/2.0 ? gHat/(kappa*kappa)*(1.0+disAmp * sin(2. * pi * t / distPeriod)) : gHat/(kappa*kappa));
  hte[left] = dirichlet(distHteInlet(t));

  // steady-uniform inlet BC
//   h[left] = dirichlet(normalDepth);
//   u.n[left] = dirichlet(normalVel);
//   te[left] = dirichlet(0.50*pow(normalVel,2.0)+0.50*(G*normalDepth+(bigPhiVal+smallPhiVal)*pow(normalDepth,2.0)));

  q[right] = neumann(0.);
  h[right] = neumann(0.);
  hte[right] = neumann(0.);
  hsphi[right] = neumann(0.);
//   u.n[right] = dirichlet(normalDepth);
//   h[right] = dirichlet(normalVel);

  foreach(){
    h[] = normalDepth;
    // simple velocity field for now
//     u.x[] = normalVel;
    q[] = normalVel*normalDepth;
//     sphi[] = gHat/totalDepth/pow(kappa,2.0);
    hsphi[] = gHat/pow(kappa,2.0);
    hte[] = normalDepth*(0.50*pow((normalVel),2.0)+0.5*gPrime*normalDepth+0.5*pow(normalDepth,2.0)*(bigPhiNVal+hsphi[]/normalDepth));
    }
}

static double momXFric(double u, double v, double h, double te, double sphi, double cf)
{
  double firstTerm, secondTerm;
  firstTerm = (1-alpha1/kappa*pow(cf,0.5))*(gHat*h-cf*pow(u,2.0));
  secondTerm = (kappa*alpha2-alpha*alpha1*pow(cf,0.5))*h*pow(cf,0.5)*(h*sphi-gHat/pow(kappa,2.0));
  return (firstTerm+secondTerm);
}

static double teFric(double u, double v, double h, double te, double sphi, double cf)
{
  double firstTerm, secondTerm, thirdTerm, bigPhi;
  bigPhi = (2.0*te-pow(u,2.0)-gPrime*h)/(pow(h,2.0))-sphi;
  firstTerm = (1.0-alpha*pow(cf,0.5)/kappa)*(gHat*h-cf*u*u)*u;
  secondTerm = alpha*alpha*h*cf*u*(h*sphi-gHat/pow(kappa,2.0));
//   thirdTerm = crCoeff/2.0*pow(h,3.0)*pow(abs(bigPhi), (3.0/2.0));
  thirdTerm = bigPhi>0 ? crCoeff/2.0*pow(h,3.0)*pow(bigPhi, (3.0/2.0)) : 0.0;
  return (firstTerm-secondTerm-thirdTerm);
}

static double sphiFric(double u, double v, double h, double te, double sphi, double cf)
{
    double firstTerm, secondTerm;
    // the first term
    firstTerm = 2.0*alpha2/kappa*pow(cf,0.5)/pow(h,2.0)*u*(cf*u*u-gHat*h);
    // the second term
    secondTerm = 2.0*alpha2*(kappa+alpha*pow(cf,0.5))*pow(cf,0.5)/h*u*(h*sphi-gHat/pow(kappa,2.0));
    return (firstTerm-secondTerm);
}

/**
  Source term.
*/
// #if 0
event frictionTerm (i++) {
  foreach(){
    // for RK2TVD
    double cf, huMed, hteMed, hsphiMed;
    double uMed, teMed, sphiMed;
    // for RK3TVD
//     double cf, huMed1, hteMed1, hsphiMed1;
//     double huMed2, hteMed2, hsphiMed2;
//     double uMed, teMed, sphiMed;

    cf = kappa*kappa*pow((RInf-2.0+2.0*log(2.0)+log(kappa)+log(pow((gHat*h[]*h[]*h[]),0.50)/nuVis)),-2.0);

    // NOTE: we do not consider 2D extension for now, and we assume the flow is uidirectional
    // RK2TVD
    uMed = q[]/h[];
    teMed = hte[]/h[];
    sphiMed = hsphi[]/h[];
    huMed = q[] + dt * momXFric(uMed, 0.0, h[], teMed, sphiMed, cf);
    hteMed = hte[] + dt * teFric(uMed, 0.0, h[], teMed, sphiMed, cf);
    hsphiMed = hsphi[] + dt * sphiFric(uMed, 0.0, h[], teMed, sphiMed, cf);

    uMed = huMed/h[];
    teMed = hteMed/h[];
    sphiMed = hsphiMed/h[];
    q[] = (1. / 2.) *(q[]) + (1. / 2.) * huMed + (1. / 2.) * dt * momXFric(uMed, 0.0, h[], teMed, sphiMed, cf);
    hte[] = (1. / 2.) *(hte[]) + (1. / 2.) * hteMed + (1. / 2.) * dt * teFric(uMed, 0.0, h[], teMed, sphiMed, cf);
    hsphi[] = (1. / 2.) *(hsphi[]) + (1. / 2.) * hsphiMed + (1. / 2.) * dt * sphiFric(uMed, 0.0, h[], teMed, sphiMed, cf);

    // RK3VD
//     huMed1 = u.x[]*h[] + dt *  momXFric(u.x[], 0.0, h[], te[], sphi[], cf);
//     hteMed1 = te[]*h[] + dt *  teFric(u.x[], 0.0, h[], te[], sphi[], cf);
//     hsphiMed1 = te[]*h[] + dt *  sphiFric(u.x[], 0.0, h[], te[], sphi[], cf);
//
//     uMed = huMed1/h[];
//     teMed = hteMed1/h[];
//     sphiMed = hsphiMed1/h[];
//     huMed2 = (3. / 4.) * (u.x[]*h[]) + (1. / 4.) * (huMed1) + (1. / 4.) * dt * momXFric(uMed, 0.0, h[], teMed, sphiMed, cf);
//     hteMed2 = (3. / 4.) * (te[]*h[]) + (1. / 4.) * (hteMed1) + (1. / 4.) * dt * teFric(uMed, 0.0, h[], teMed, sphiMed, cf);
//     hsphiMed2 = (3. / 4.) * (sphi[]*h[]) + (1. / 4.) * (hsphiMed1) + (1. / 4.) * dt * sphiFric(uMed, 0.0, h[], teMed, sphiMed, cf);
//
//     uMed = huMed2/h[];
//     teMed = hteMed2/h[];
//     sphiMed = hsphiMed2/h[];
//     u.x[] = ((1. / 3.) * (u.x[]*h[]) + (2. / 3.) * (huMed2) + (2. / 3.) * dt *  momXFric(uMed, 0.0, h[], teMed, sphiMed, cf))/h[];
//     te[] = ((1. / 3.) * (te[]*h[]) + (2. / 3.) * (hteMed2) + (2. / 3.) * dt * teFric(uMed, 0.0, h[], teMed, sphiMed, cf))/h[];
//     sphi[] = ((1. / 3.) * (sphi[]*h[]) + (2. / 3.) * (hsphiMed2) + (2. / 3.) * dt * sphiFric(uMed, 0.0, h[], teMed, sphiMed, cf))/h[];
  }

//     boundary ((scalar *){u, te, sphi});
}
// #endif

/*
 * AMR here
 *
 */
event adapt1 (i++) {
  scalar bigPhiField[];
  vector gh[], gq[], ghte[], ghsphi[], gbigPhi[];

  foreach()
  {
    bigPhiField[] = (2.0*hte[]/h[]-pow((q[]/h[]),2.0)-gPrime*h[])/(pow(h[],2.0))-hsphi[]/h[];
  }

  gradients ({h}, {gh});
  gradients ({q}, {gq});
  gradients ({hte}, {ghte});
  gradients ({hsphi}, {ghsphi});
  gradients ({bigPhiField}, {gbigPhi});

  foreach()
  {
    hGrad[] = norm(gh);
    qGrad[] = norm(gq);
    hteGrad[] = norm(ghte);
    hsphiGrad[] = norm(ghsphi);
    bigPhiGrad[] = norm(gbigPhi);
  }

  adapt_wavelet({h, hGrad, qGrad, hsphiGrad, hteGrad, bigPhiGrad}, (double[]){normalDepth/(300.0), (2.e-3), (2.e-3), (2.e-3), (2.e-3), (2.e-3)}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // adapt_wavelet({h, hGrad}, (double[]){normalDepth/500.0, 0.0025}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

  refine(x<=(2.0*distPeriod*normalVel) && t<=(distPeriod) && level<MAXLEVEL);
  refine(x<=(17.0*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
  refine(x>=(DOMAINLENGTH-17.0*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);

  // refine(x<=(2.0*distPeriod*normalVel) && t<=(distPeriod) && level<MAXLEVEL);
  // refine(x<=(DOMAINLENGTH) && t>(distPeriod) && level<MAXLEVEL);
}

// record max depth
event hmax1(t=0; t<4.0; i+=20)
{
     double maxDepth = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     double minDepth = 0.0;
     double minDepthLocX = 0.0;
	 double max2LocX = 0.0;

     FILE *fp2 = fopen("maxMinDepth", "a+");
     // FR: maxDepth, minDepth, H, waveLength
     foreach ()
     {
               if ( h[]>h[1] && h[]>h[-1] && x>maxDepthLocX ){
                    maxDepth = h[];
                    maxDepthLocX = x;
                    maxDepthVel = q[]/h[];
               }
     }
	 foreach ()
	 {
		if (h[]<h[1] && h[]<h[-1] && h[-1]<h[-2] && h[1]<h[2] && x>minDepthLocX && x<maxDepthLocX ) {
                    minDepthLocX = x;
                    minDepth = h[];
               }

		else if (h[]>h[1] && h[]>h[-1] && h[1]>h[2] && h[-1]>h[-2] && h[]>1.1 && x<maxDepthLocX && x>max2LocX) {
                    max2LocX = x;
               }
	 }
     fprintf(fp2, "%g %13.9f %13.9f %13.9f %13.9f %13.9f \n", t, maxDepthLocX, maxDepth, minDepth, (maxDepth-minDepth), (maxDepthLocX-max2LocX));
     fclose(fp2);
}

event hmax2(t=4.0; t<simTime; i+=20)
{
     double maxDepth = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     double minDepth = 0.0;
     double minDepthLocX = 0.0;
	 double max2LocX = 0.0;

     FILE *fp2 = fopen("maxMinDepth", "a+");
     // FR: maxDepth, minDepth, H, waveLength
     foreach ()
     {
               if ( h[]>h[1] && h[]>h[-1] && x>maxDepthLocX ){
                    maxDepth = h[];
                    maxDepthLocX = x;
                    maxDepthVel = q[]/h[];
               }
     }
	 foreach ()
	 {
		if (h[]<h[1] && h[]<h[-1] && h[-1]<h[-2] && h[1]<h[2] && x>minDepthLocX && x<maxDepthLocX ) {
                    minDepthLocX = x;
                    minDepth = h[];
               }

		else if (h[]>h[1] && h[]>h[-1] && h[1]>h[2] && h[-1]>h[-2] && h[]>1.1 && x<maxDepthLocX && x>max2LocX) {
                    max2LocX = x;
               }
	 }
     fprintf(fp2, "%g %13.9f %13.9f %13.9f %13.9f %13.9f \n", t, maxDepthLocX, maxDepth, minDepth, (maxDepth-minDepth), (maxDepthLocX-max2LocX));
     fclose(fp2);
}

/**
save the hight the flux and the yield surface as a function of time
*/
event output  (t = 0; t <= simTime; t+=outputInterval) {
  sprintf (s, "slice-%g.txt", t);
  FILE * fp2 = fopen (s, "w");
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g %g \n", x, h[], q[], hte[], hsphi[]);
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
