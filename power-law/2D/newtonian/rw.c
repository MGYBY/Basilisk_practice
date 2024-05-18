/**
# Steady-uniform flow (normal flow) simulation
Based on dimensional variables.

## Problem

# Code
*/
// #include "grid/cartesian1D.h"
// #include "grid/bitree.h"
#include "grid/quadtree.h"
// #include "./saint-venantNN.h"
#include "./saint-venant-power-law.h"
// #include "./output_vts.h"
#include "vtk.h"

// dimless parameters defining the problem
#define FR 0.90
#define So 0.060
#define RE (3.0*FR*FR/So)
#define n_coeff 1.0

#define normalDepth 1.00
#define normalVel 1.00

#define piVal 3.14159265
#define disAmp 0.20
// #define hBC(t) (t<=distPeriod/2.0 ? (normalDepth*(1.0+disAmp*sin(2. * pi * t / distPeriod))) : (normalDepth))

#define Lambda 20.0 // the Gaussian-function parameter
#define DOMAINLENGTH (120.0*Lambda)

#define refLen (Lambda*10.0)
#define Ly (Lambda*55.0)

// center of initial disturbance
#define x0 (DOMAINLENGTH/2.0)
#define y0 (Ly/2.0)

#define MAXLEVEL 11
#define MINLEVEL 2
#define INITLEVEL 8

#define simTime (1200.0/So)
#define outputInterval (round(5.0/So))

#define hError (normalDepth/500.0)
#define hmaxError (normalDepth/700.0)
#define uError (normalVel/500.0)
#define hGradError (0.002)
#define uxGradError (0.002)

double hSurf (double xCoord, double yCoord)
{
  return (normalDepth*(1.0+disAmp*exp(-0.5*(1.0/pow(Lambda, 2.0))*(pow((xCoord-x0), 2.0)+pow((yCoord-y0), 2.0)))));
}

char s[40], vtkName[40], gerrisName[40];
scalar uAve[], hmax[], hGrad[], uxGrad[];

int main() {
  L0 = DOMAINLENGTH;
  // G  = gravity;
  // change to g' for rotated axis frame
  G = (1.0/(FR*FR));
  betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
  init_grid(1 << (INITLEVEL));
  nu = 1.; // dummy

  CFL = 0.449;

  // slightly higher resolution
  theta = 1.30;

  periodic (right);
//   periodic(top);

  run();
}

/**
## Initialization  */
event init (i = 0) {
  refine(y>=(Ly-1.1*DOMAINLENGTH/pow(2, MAXLEVEL)) && y<=(Ly+1.1*DOMAINLENGTH/pow(2, MAXLEVEL)) && level<MAXLEVEL);
  refine(fabs(y0-y)<=(refLen/2.0) && fabs(x0-x)<=(refLen/2.0) && level<MAXLEVEL);

  // for 2-D case 40x3m
  mask(y > Ly ? top : none);

  /**
  We initialize *h*. */
  foreach(){
    double totalDepth = hSurf(x, y);
    // periodic BC cannot use realistic topo
    // zb[] = -x*sinTheta;
    zb[] = 0.0;
    h[] = totalDepth;

    hGrad[] = 0.0;
    uxGrad[] = 0.0;

    hmax[] = h[];
    u.x[] = normalVel*pow(h[]/normalDepth, 2.0);
    u.y[] = 0.0;

    uAve[] = norm(u);

//     depthGrad[] = fabs(h[-1]-h[1])/(2.*Delta);
    }
}

  /**
## Output
  We print the elevation  */
event acceleration (i++) {
  foreach ()
  {
      double a = h[] < dry ? HUGE : 1.0 + (3.0/RE) * dt / pow(h[], 2.0);
      u.x[] = (u.x[] + 3.0*dt/RE) / a;
      u.y[] /= a;

      // RK3TVD below
      // uMed = u.x[] + dt *  (-(cf / 2.) * u.x[] * norm(u) / h[] + gravityCoeff*So);
      // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * u.x[] * sqrt(sq(uMed)+sq(u.y[])) / h[] + gravityCoeff*So);
      // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * u.x[] * sqrt(sq(uMed)+sq(u.y[])) / h[] + gravityCoeff*So);

      // uMed = u.y[] + dt *  (-(cf / 2.) * u.y[] * norm(u) / h[]);
      // uMed = (3. / 4.) * u.y[] + (1. / 4.) * uMed + (1. / 4.) * dt * (-(cf / 2.) * u.y[] * sqrt(sq(uMed)+sq(u.x[])) / h[]);
      // u.y[] = (1. / 3.) * u.y[] + (2. / 3.) * uMed + (2. / 3.) * dt * (-(cf / 2.) * u.y[] * sqrt(sq(uMed)+sq(u.x[])) / h[]);

      uAve[] = norm(u);

//       if (h[]>hmax[]) {
//         hmax[]=h[];
//       }
  }

  boundary ((scalar *){u});
}

/*
 * AMR here
 *
 */
event adapt1 (i++) {
  // adapt_wavelet({h, depthGrad, uAve}, (double[]){normalDepth/200.0, 0.01, normalVel/200.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // a more reasonable AMR criteria
//   adapt_wavelet({h, depthGrad, uAve, yieldSurf}, (double[]){normalDepth/300.0, 0.007, normalVel/300.0, normalDepth/150.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  vector gh[], gux[];
  gradients ({h}, {gh});
  gradients ({u.x}, {gux});
  foreach()
  {
    hGrad[] = norm(gh);
    uxGrad[] = norm(gux);
  }

//   adapt_wavelet({h, hmax, uAve, hGrad, uxGrad}, (double[]){hError, hmaxError, uError, hGradError, uxGradError}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  adapt_wavelet({h, uAve, hGrad, uxGrad}, (double[]){hError, uError, hGradError, uxGradError}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);

//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

  refine(fabs(y0-y)<=(refLen/2.0) && fabs(x0-x)<=(refLen/2.0) && t<=(0.95*refLen/2.0/(1+disAmp)) && level<MAXLEVEL);
//   refine(x<=(9.9*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
//   refine(x>=(DOMAINLENGTH-9.9*DOMAINLENGTH/pow(2,MAXLEVEL)) && level<MAXLEVEL);
}

event hmaxFR(i+=20)
{
     double maxDepthLocX = 0.0;
     double pM2X = 0.0;
     double pM1X = 0.0;
     double ppX = 0.0;
     double pPL1X = 0.0;
     double FRH = 0.0;
     double FRV = 0.0;
     double FRX = 0.0;
     double hM2 = 0.0;
     double hM1 = 0.0;
     double hp = 0.0;
     double hPL1 = 0.0;

     FILE *fp2 = fopen("maxMinDepth_global", "a+");
     FILE *fp3 = fopen("maxMinDepth_FR", "a+");

     // global max here
     stats s1 = statsf (h);
     stats s2 = statsf (uAve);

     foreach ()
     {
       if (h[]==s1.max)
       {
         maxDepthLocX = x;
      }
    }
    //     fprintf(fp2, "%g %13.9f %13.9f %13.9f \n", t, maxDepthLocX, s1.max, s2.max );
    fprintf(fp2, "%g %g %g %g \n", t, maxDepthLocX, s1.max, s2.max );
    fclose(fp2);

    // FR max here
    FRH = normalDepth;
    FRV = normalVel;
    for (int l = 3; l <= (pow(2,MAXLEVEL)-1); l++)
    {
      ppX = DOMAINLENGTH/pow(2,MAXLEVEL)*l-0.5*DOMAINLENGTH/pow(2,MAXLEVEL);
      pM2X = ppX-2.0*DOMAINLENGTH/pow(2,MAXLEVEL);
      pM1X = ppX-1.0*DOMAINLENGTH/pow(2,MAXLEVEL);
      pPL1X = ppX+1.0*DOMAINLENGTH/pow(2,MAXLEVEL);
      hp = interpolate(h, ppX, y0, 0.0);
      hM2 = interpolate(h, pM2X, y0, 0.0);
      hM1 = interpolate(h, pM1X, y0, 0.0);
      hPL1 = interpolate(h, pPL1X, y0, 0.0);;
      if (hp>hPL1 && hp>=hM1 && hM1>=hM2 && hp>normalDepth*1.005 && ppX>FRX)
      {
        FRX = ppX;
        FRH = hp;
        FRV = interpolate(uAve, ppX, y0, 0.0);
      }
    }
    fprintf(fp3, "%g %g %g %g \n", t, FRX, FRH, FRV);
    fclose(fp3);
}

/**
save the hight the flux and the yield surface as a function of time
*/ 
event output  (t = 0; t <= simTime; t+=outputInterval){
  // sprintf (s, "slice-%g.txt", t);
  // sprintf(vtkName, "outVTK-%g.vtk", t);
  sprintf(gerrisName, "snapshot-%g.gfs", t);
  // FILE * fp1 = fopen (s, "w");
  // FILE * fp3 = fopen(vtkName, "w");
  FILE * fp4 = fopen(gerrisName, "w");

  sprintf (s, "depth-vel_%g.txt", t);
  FILE * fp2 = fopen (s, "w");
  foreach (serial) {
    fprintf (fp2, "%g %g %g %g %g \n", x, y, h[], u.x[], u.y[]);
  }
  fclose(fp2);

  // output_vtk((scalar *) {h, uAve}, 1 << MAXLEVEL, (FILE *) fp3, false);
  // fclose(fp3);
  output_gfs(fp4, list = {h, u.x, u.y}, translate = true);
  fclose(fp4);
}
