// #include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "saint-venant-flexible.h"
#include "utils.h"

// problem-specific parameters
#define FR 0.55
#define Fre 0.80
#define n1Val 10.0
#define n2Val -1.50

// double alpha_coeff = 40.0;
// double lx = 2.878;
#define disMag 0.20
// double betaCoeff = 0.0;
#define simTime 48.0
// double distPeriod = 1.0;
#define wavelength 100.0
#define OUTPUTTIME 2.00
#define DOMAINLENGTH (wavelength)

// #define INITDEPTHFUNC(amp, len, xCoord) (1.0*(1.0+amp*sin(2.0*M_PI*(xCoord-len/2.0)/len)))
// #define INLETDEPTHBC(amp, time, period) (1.0*(1.0+amp*sin(2.0*M_PI*(time)/period)))

#define MAXLEVEL 10
#define MINLEVEL 2
#define MAXMAXLEVEL 10

#define WARMUP (0.00)

// h[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? INLETDEPTHBC(disMag, (t-WARMUP), distPeriod) : 1.0 );
// u.n[left] = dirichlet( (t-WARMUP)<=(distPeriod/2.0) && (t-WARMUP)>=0 ? pow(INLETDEPTHBC(disMag, (t-WARMUP), distPeriod), 0.5) : 1.0 );
//
// u.n[right] = neumann(0.);
// h[right] = neumann(0.);

int main()
{
     L0 = DOMAINLENGTH;
//      N = 2048;
     init_grid (1 << MAXLEVEL);
     G = 1.0/(Fre*Fre*(n1-n2));
     n1 = n1Val;
     n2 = n2Val;
     // nCoeff = n_coeff;
     
     CFL = 0.399; // CFL number should be sufficiently small

     theta = 1.205;

     // note that we are using the hump periodic configuration
     periodic (right);

     run();
}

scalar te[], re[];
// scalar inletRef[];
scalar depthGrad[];

event init(i = 0)
{
     foreach ()
     {
               zb[] = 0.0;
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               // h[] = 1.0;
//                h[] = (x>=wavelength/2.0 && x<=wavelength) ? INITDEPTHFUNC(disMag, wavelength, x) : 1.0;
               // h[] = 1.0;
               h[] = 1.0+disMag*sin(2.0*M_PI*x/(DOMAINLENGTH));
//                u.x[] = (x>=wavelength/2.0 && x<=wavelength) ? pow(INITDEPTHFUNC(disMag, wavelength, x), 0.5) : 1.0;
               // u.x[] = 1.0;
               u.x[] = 1.0;

               depthGrad[] = fabs(h[-1]-h[1])/Delta;
     }

}

event calcDepthGrad(i++)
{
     foreach()
     {
          depthGrad[] = fabs(h[-1]-h[1])/Delta;
     }
}

static double vesselFriction(double u, double h, double froude)
{
     double rhs;
     // rhs = 1.0/(froude*froude)*(h-1.0/(pow(h, 0.50))*u);
     rhs = 1.0/(froude*froude)*(1.0-1.0/(pow(h, 1.50))*u);
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
               uMed = u.x[] + dt * vesselFriction(u.x[], h[], Fr);
               uMed = (3.0/4.0)*u.x[] + (1.0/4.0)*uMed + (1.0/4.0)*dt*vesselFriction(uMed, h[], Fr);
               u.x[] = (1.0/3.0)*u.x[]+(2.0/3.0)*uMed+(2.0/3.0)*dt*vesselFriction(uMed, h[], Fr);
          }
          else {
               u.x[] = 0.0;
          }

          //linearized backeard Euler
          // u.x[] = h[]>dry && u.x[]>1.0e-3  ? (u.x[] + dt)/(1.0+dt*(1.0/h[])*((pow(u.x[], (n_coeff-1.0)))/(pow(h[], n_coeff)))) : 0.0;
          // uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf);
          // uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf);
          // u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf);

          // u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * u.x[] / (h[]));

//           re[] = pow(u.x[], (2.0-n_coeff))*pow(h[], n_coeff);
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
     stats s1 = statsf (re);
     fprintf(fp2, "%g %g %g %g %g \n", t, maxDepth, maxDepthLocX, maxDepthVel, s1.max);
     fclose(fp2);
}

// event energyContent(i+=10)
// {
//      FILE *fp2 = fopen("totalEnergy", "a+");
//      foreach ()
//      {
//           te[] = 0.50*u.x[]*u.x[]+1.0/(pow(FR, 2.0))*h[];
//      }
//      stats s2 = statsf (te);
//      fprintf(fp2, "%g %g \n", t, s2.sum);
//      fclose(fp2);
// }

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
     sprintf(name, "out-%g.txt", t);
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
// AMR features
event adapt1 (i++) {
  astats s = adapt_wavelet({h, depthGrad}, (double[]){1/400.0, 0.0006}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//      astats s = adapt_wavelet({h, depthGrad}, (double[]){1.0/300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  refine(x<=10.0*L0/pow(2, MAXLEVEL) && level<MAXLEVEL);
  refine(x>=L0-10.0*L0/pow(2, MAXLEVEL) && level<MAXLEVEL);
}
