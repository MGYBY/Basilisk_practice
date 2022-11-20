#include "grid/multigrid1D.h"
#include "saint-venant-power-law.h"
#include "utils.h"

#define FR (1.756)
#define n_coeff 0.30
#define ALPHA (1/pow(FR,2.0))
// double alpha_coeff = 40.0;
// double lx = 2.878;
#define disMag 0.20
// double betaCoeff = 0.0;
#define simTime 120.0
// double distPeriod = 1.0;
#define wavelength (3.0/pow(FR,2.0))
#define OUTPUTTIME 0.50
#define DOMAINLENGTH (200.0/pow(FR,2.0))

#define INITDEPTHFUNC(amp, len, xCoord) (1.0*(1.0+amp*sin(2.0*M_PI*(xCoord-len/2.0)/len)))

#define MAXLEVEL 15
#define MINLEVEL 1
#define MAXMAXLEVEL 13

int main()
{
     L0 = DOMAINLENGTH;
     N = 2048;
     alphaCoeff = ALPHA;
     G = ALPHA;
     betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
     // nCoeff = n_coeff;
     
     CFL = 0.30; // CFL number should be sufficiently small

     theta = 1.225;

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
               h[] = (x>=wavelength/2.0 && x<=wavelength) ? INITDEPTHFUNC(disMag, wavelength, x) : 1.0;
               u.x[] = (x>=wavelength/2.0 && x<=wavelength) ? pow(INITDEPTHFUNC(disMag, wavelength, x), 0.5) : 1.0;

//                inletRef[] = x<=inletRefLen ? 1.0 : 0.0;
               depthGrad[] = (h[-1]-h[1])/Delta;
     }

}

event calcDepthGrad(i++)
{
     foreach()
     {
          depthGrad[] = (h[-1]-h[1])/Delta;
//           inletRef[] = x<=inletRefLen ? 1.0 : 0.0;
     }
}

static double powerLawFriction(double u, double h, double n)
{
     double rhs;
     rhs = 1.0-pow((u/h), n)/h;
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

          re[] = pow(u.x[], (2.0-n_coeff))*pow(h[], n_coeff);
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

event energyContent(i+=10)
{
     FILE *fp2 = fopen("totalEnergy", "a+");
     foreach ()
     {
          te[] = 0.50*u.x[]*u.x[]+1.0/(pow(FR, 2.0))*h[];
     }
     stats s2 = statsf (te);
     fprintf(fp2, "%g %g \n", t, s2.sum);
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
          fprintf(fp, "%g %g %g %g \n", x, h[], u.x[], re[]);
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
//   astats s = adapt_wavelet({h, depthGrad, inletRef}, (double[]){normalDepth / 300.0, 0.00016, 0.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
     astats s = adapt_wavelet({h, depthGrad}, (double[]){normalDepth / 300.0, 0.00016}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // astats s = adapt_wavelet({ depthGrad}, (double[]){ 0.00010}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
