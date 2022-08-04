#include "grid/multigrid1D.h"
#include "saint-venant-power-law.h"

#define FR 0.250
double n_coeff = 0.40;
#define ALPHA (1/pow(FR,2.0))
double alpha_coeff = 1.0;
// double lx = 2.878;
double disMag = 0.10;
// double betaCoeff = 0.0;
double simTime = 400.0;
double distPeriod = pi;

int main()
{
     L0 = 3200.0;
     N = 32768;
//      alphaCoeff = alpha_coeff;
	/* Note that Kurganov RS uses G as alhpa! */
	G = ALPHA;
     betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
     // nCoeff = n_coeff;
	
	theta = 1.35;
	
	// sometimes it is necessary to tune the dry threshold
	dry = 1e-7;
     
     CFL = 0.25; // CFL number should be sufficiently small

     run();
}

event init(i = 0)
{
     h[left] = dirichlet( t<distPeriod ? 1.0 + disMag * 1.0 * sin(2. * pi * t / distPeriod) : 1.0);
     u.n[left] = dirichlet(1.0);

     u.n[right] = neumann(0.);
     h[right] = neumann(0.);

     foreach ()
     {
               zb[] = 0.0;
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               h[] = 1.0;
               // eta[] = 1.0 + disMag * sin(2. * pi * x / lx);
               u.x[] = 1.0;
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
//      double uMed=0.0;
     double uMed1 = 0.0;
     double uMed2 = 0.0;
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
//           if (h[] > dry) {
//                uMed = u.x[] + dt * powerLawFriction(u.x[], h[], n_coeff);
//                uMed = (3.0/4.0)*u.x[] + (1.0/4.0)*uMed + (1.0/4.0)*dt*powerLawFriction(uMed, h[], n_coeff);
//                u.x[] = (1.0/3.0)*u.x[]+(2.0/3.0)*uMed+(2.0/3.0)*dt*powerLawFriction(uMed, h[], n_coeff);
//           }
//           else {
//                u.x[] = 0.0;
//           }
          
          // rk104
          if (h[] > dry) {
               uMed1 = u.x[];
               uMed2 = u.x[];
               for (int l = 0; l <= 5; l++) {
                    uMed1 += dt * powerLawFriction(uMed1, h[], n_coeff) /6.0;
               }
               uMed2 = (1.0/25.0)*uMed2+(9.0/25.0)*uMed1;
               uMed1 = 15.0*uMed2-5.0*uMed1;
               for (int l = 6; l <= 9; l++) {
                    uMed1 += dt * powerLawFriction(uMed1, h[], n_coeff) /6.0;
               }
               u.x[] = uMed2 + (3.0/5.0)*uMed1+(1.0/10.0)*dt*powerLawFriction(uMed1, h[], n_coeff) /6.0;
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
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

// record max depth
event hmax1(t=0; t<20.0; t+= 0.01)
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
                    maxDepthVel = u.x[];
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

event hmax2(t=20.0; t<simTime; t+= 0.06)
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
                    maxDepthVel = u.x[];
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

event gnuplot(t = 0; t <= simTime; t += 4)
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

event output(t = 0; t <= simTime; t += 4)
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
