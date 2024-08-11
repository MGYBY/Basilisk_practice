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
#include "saint-venant_r_new.h"
#include "tracer.h"

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
#define disAmp 0.20
#define distPeriod 0.933

#define DOMAINLENGTH (2500.0*normalDepth/sinTheta)

#define INITLEVEL 12
#define MAXLEVEL 16
#define MINLEVEL 1

#define simTime 250.0
#define outputInterval 0.50

#define gradTol 2.50e-3
#define initRefineStage (distPeriod*16.0)
#define bigPhiDryTol 1.0e-30

#define concInit 1.00

char s[80];
scalar hmax[], hGrad[], uxGrad[], teGrad[], sphiGrad[], bigPhiGrad[];
scalar bigPhiField[];
scalar conc[];
scalar * tracers = {conc};
face vector uf[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  // FIXME: for 2D: G is gravity
  G = gPrime;
  init_grid(1 << (INITLEVEL));

  // theta = 1.30;
  theta = 1.530;
  CFL = 0.425;

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
    
  //return normalDepth;

  // "type-a" inlet dist
//   return (normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod));

  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

double distuInlet(double t)
{
  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
//   "type-c", const Fr inlet dist
  if (t<=distPeriod/2.0)
    return (froude*sqrt(gPrime*(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod))));
    //return (normalVel*(1.0+disAmp * sin(2. * pi * t / distPeriod)));
  else
    return (normalVel);

//   "type-c", const q inlet dist
//   if (t<=distPeriod/2.0)
//     return (normalDepth*normalVel/(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod)));
//   else
//     return normalVel;

  // "type-a" inlet dist
//   return (froude*sqrt(G*(normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod))));

  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

double distteInlet(double t)
{
  double hVar = 0.0;
  double vVar = 0.0;
  double phiVar = 0.0;
  hVar = normalDepth + disAmp * normalDepth * sin(2. * pi * t / distPeriod);
  //hVar = normalDepth;
  phiVar= gHat/(kappa*kappa)/hVar;
  vVar = froude*sqrt(gPrime*hVar);
  //vVar = normalVel*(1.0+disAmp * sin(2. * pi * t / distPeriod));

  if (t<=distPeriod/2.0)
    return ((pow(vVar,2.0)/2.0 + pow(hVar,2.0)*phiVar/2.0 + 0.0 + gPrime*hVar/2.0));
  else
    return ((pow(normalVel,2.0)/2.0 + pow(normalDepth,2.0)*(gHat/(kappa*kappa)/normalDepth)/2.0 + 0.0 + gPrime*normalDepth/2.0));

}

/**
## Initialization  */
event init (i = 0) {
  // "wavemaker" inlet BC
  h[left] = dirichlet(distDepthInlet(t));
  u.x[left] = dirichlet(distuInlet(t));
  sphi[left] = dirichlet(t<=distPeriod/2.0 ? gHat/(kappa*kappa)/distDepthInlet(t) : gHat/(kappa*kappa)/normalDepth);
  te[left] = dirichlet(distteInlet(t));
  conc[left] = dirichlet(concInit);

  // steady-uniform inlet BC
//   h[left] = dirichlet(normalDepth);
//   u.n[left] = dirichlet(normalVel);
//   te[left] = dirichlet(0.50*pow(normalVel,2.0)+0.50*(G*normalDepth+(bigPhiVal+smallPhiVal)*pow(normalDepth,2.0)));

  u.x[right] = neumann(0.);
  h[right] = neumann(0.);
  te[right] = neumann(0.);
  sphi[right] = neumann(0.);
  conc[right] = neumann(0.);
//   u.n[right] = dirichlet(normalDepth);
//   h[right] = dirichlet(normalVel);

  /**
  We initialize *h*. */
//   fprintf (ferr, "Starting initizatiation \n");
  foreach(){
    // double totalDepth = normalDepth*(1.0+disAmp*sin(2. * pi * x / DOMAINLENGTH));
    zb[] = 0.0;
    h[] = normalDepth;
    // simple velocity field for now
//     u.x[] = normalVel;
    u.x[] = normalVel;
//     sphi[] = gHat/totalDepth/pow(kappa,2.0);
    sphi[] = gHat/normalDepth/pow(kappa,2.0);
    te[] = 0.50*pow(u.x[],2.0)+0.5*G*h[]+0.5*pow(h[],2.0)*(bigPhiNVal+sphi[]);
    conc[] = concInit;
    }

  foreach_face()
    uf.x[] = (u.x[] + u.x[-1])/2.;
//     fprintf (ferr, "Initizatiation finishes \n");
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
  bigPhi = (2.0*te-pow(u,2.0)-G*h)/(pow(h,2.0))-sphi;
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
    huMed = u.x[]*h[] + dt * momXFric(u.x[], 0.0, h[], te[], sphi[], cf);
    hteMed = te[]*h[] + dt * teFric(u.x[], 0.0, h[], te[], sphi[], cf);
    hsphiMed = sphi[]*h[] + dt * sphiFric(u.x[], 0.0, h[], te[], sphi[], cf);

    uMed = huMed/h[];
    teMed = hteMed/h[];
    sphiMed = hsphiMed/h[];
    u.x[] = ((1. / 2.) *(u.x[]*h[]) + (1. / 2.) * huMed + (1. / 2.) * dt * momXFric(uMed, 0.0, h[], teMed, sphiMed, cf))/h[];
    te[] = ((1. / 2.) *(te[]*h[]) + (1. / 2.) * hteMed + (1. / 2.) * dt * teFric(uMed, 0.0, h[], teMed, sphiMed, cf))/h[];
    sphi[] = ((1. / 2.) *(sphi[]*h[]) + (1. / 2.) * hsphiMed + (1. / 2.) * dt * sphiFric(uMed, 0.0, h[], teMed, sphiMed, cf))/h[];

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

// tracer advection
event tracerTerm (i++) {
  // compute uf
  foreach_face()
    uf.x[] = (u.x[] + u.x[-1])/2.;
}

/*
 * AMR here
 *
 */
 //event adapt1 (i++) {
event adapt1 (t=0; t<=initRefineStage; i++) {
  vector gh[], gux[], gte[], gsphi[], gbigPhi[];

  foreach()
  {
    bigPhiField[] = (2.0*te[]-pow((u.x[]),2.0)-gPrime*h[])/(pow(h[],2.0))-sphi[];
    
    if (bigPhiField[])
    {
        te[] = bigPhiField[]<bigPhiDryTol ? 0.5*(pow(u.x[],2.0)+pow(h[],2.0)*sphi[]+0.0+gPrime*h[]) : te[];
    }
  }

  gradients ({h}, {gh});
  gradients ({u.x}, {gux});
  gradients ({te}, {gte});
  gradients ({sphi}, {gsphi});
  gradients ({bigPhiField}, {gbigPhi});

  foreach()
  {
    hGrad[] = norm(gh);
    uxGrad[] = norm(gux);
    teGrad[] = norm(gte);
    sphiGrad[] = norm(gsphi);
    bigPhiGrad[] = norm(gbigPhi);
  }

  adapt_wavelet({h, hGrad, uxGrad, teGrad, sphiGrad, bigPhiGrad}, (double[]){normalDepth/350.0, gradTol, gradTol, gradTol, gradTol, gradTol}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//     adapt_wavelet({h}, (double[]){normalDepth/450.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

  refine(x<=(2.50*distPeriod*normalVel) && t<=(distPeriod*1.50) && level<MAXLEVEL);
  refine(x<=(17.510*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
  refine(x>=(DOMAINLENGTH-17.510*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
}

event adapt2 (t=initRefineStage; i++) {
  //vector gh[], gux[], gte[], gsphi[], gbigPhi[];
//
  foreach()
  {
    bigPhiField[] = (2.0*te[]-pow((u.x[]),2.0)-gPrime*h[])/(pow(h[],2.0))-sphi[];
    
    if (bigPhiField[])
    {
        te[] = bigPhiField[]<bigPhiDryTol ? 0.5*(pow(u.x[],2.0)+pow(h[],2.0)*sphi[]+0.0+gPrime*h[]) : te[];
    }
  }
//
  //gradients ({h}, {gh});
  //gradients ({u.x}, {gux});
  //gradients ({te}, {gte});
  //gradients ({sphi}, {gsphi});
  //gradients ({bigPhiField}, {gbigPhi});
//
  //foreach()
  //{
    //hGrad[] = norm(gh);
    //uxGrad[] = norm(gux);
    //teGrad[] = norm(gte);
    //sphiGrad[] = norm(gsphi);
    //bigPhiGrad[] = norm(gbigPhi);
  //}

  adapt_wavelet({h}, (double[]){normalDepth/(400.00)}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
   //adapt_wavelet({h, hGrad, uxGrad, teGrad, sphiGrad, bigPhiGrad}, (double[]){normalDepth/250.0, (gradTol*2500.0), (gradTol*2500.0), (gradTol*2500.0), (gradTol*2500.0), (gradTol*2500.0)}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);

  refine(x<=(17.510*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
  refine(x>=(DOMAINLENGTH-17.510*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval) {
  //foreach()
  //{
    //bigPhiField[] = (2.0*te[]-pow((u.x[]),2.0)-gPrime*h[])/(pow(h[],2.0))-sphi[];
  //}
  
  sprintf (s, "slice-%g.txt", t);
  FILE * fp2 = fopen (s, "w"); 
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g %g %g %g \n", x, h[], u.x[], te[], sphi[], bigPhiField[], conc[]);
  }
  fclose(fp2);
}

// record max depth
event hmaxUmax (i+=50)
{
     double maxDepth = 0.0;
     double maxVel = 0.0;
     double maxDepthLocX = 0.0;
     double maxDepthVel = 0.0;
     double maxVelLocX = 0.0;
     double maxVelDepth = 0.0;
  double maxConc=0.0;

     FILE *fp2 = fopen("maxDepth", "a+");
     FILE *fp3 = fopen("maxVel", "a+");

     foreach ()
     {
           if (h[]>maxDepth){
                maxDepth = h[];
                maxDepthLocX = x;
                maxDepthVel = u.x[];
        maxConc = conc[];
           }

           if (u.x[]>maxVel){
                maxVel = u.x[];
                maxVelLocX = x;
                maxVelDepth = h[];
           }
     }

     fprintf(fp2, "%.10g %.10g %.10g %.10g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel, maxConc);
     fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, maxVel, maxVelLocX, maxVelDepth);

     fclose(fp2);
     fclose(fp3);
}

event hFront1 (t=0; i+=8; t<initRefineStage)
{
     double aveDepth = 0.0;
     double aveVel = 0.0;
     double xf = 0.0;
     FILE *fp3 = fopen("frontPosMod", "a+");
     foreach ()
     {
          xf = h[] > dry ?  max(xf,x) :  xf ;
          aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
          aveVel += (h[]>=dry ? Delta*u.x[] : 0.0);
     }
     fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf) );
     fclose(fp3);
}

event hFront2 (t=initRefineStage; t<=simTime; i+=80)
{
     double aveDepth = 0.0;
     double aveVel = 0.0;
     double xf = 0.0;
     FILE *fp3 = fopen("frontPosMod", "a+");
     foreach ()
     {
          xf = h[] > dry ?  max(xf,x) :  xf ;
          aveDepth += (h[]>=dry ? Delta*h[] : 0.0);
          aveVel += (h[]>=dry ? Delta*u.x[] : 0.0);
     }

     fprintf(fp3, "%.10g %.10g %.10g %.10g \n", t, xf, (aveDepth/xf), (aveVel/xf) );
     fclose(fp3);
}

/**
Visualise the wave profile
*/

void plot_profile_h(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't-%g_h.png'\n"
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

void plot_profile_conc(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't-%g_conc.png'\n"
             "set title 't = %.2f'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2 linecolor rgb 'green' \n",
             t, t);
     foreach ()
     {
               fprintf(fp, "%g %g\n", x, conc[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplotOutput1(t = 0; t <= simTime; t += outputInterval)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile_h(t, fp);
}

event gnuplotOutput2(t = 0; t <= simTime; t += outputInterval)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile_conc(t, fp);
}

// event moviemaker (t = end) {
//     system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
// 	    "ppm2mp4 movie.mp4");
// }
