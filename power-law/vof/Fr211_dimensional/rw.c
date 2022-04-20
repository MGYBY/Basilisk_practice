/**
   This is a standard file to reproduce the simulations presented in: [Mostert and Deike, Inertial energy dissipation in shallow-water breaking waves. Journal of Fluid Mechanics, 890:A12, 2020](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/inertial-energy-dissipation-in-shallowwater-breaking-waves/8B818C95BBB1C8BF2723DB611AFB9DCD).

   We use Navier-Stokes in 2D with surface tension, and with the momentum conserving VOF scheme. */

#include "grid/quadtree.h"
#include "adapt_wavelet_leave_interface.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"
#include "tag.h"
#include "./output_vtu_foreach.h"
#include "./waveprobes.h"

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"


#define MAXLEVEL 13
#define MINLEVEL 5

// do not consider surface tension for now
// #define BO 1000.0 //Bond number
#define RE 335.0 //Reynolds number
#define FR 2.11
#define RATIO 850.0/1.0 //density ratio, water to air
#define MURATIO 8.9e-4/17.4e-6 //dynamic viscosity ratio, water to air
#define POWERLAWN 1.0

#define CHANNELSLOPE 0.040
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))
#define DISTAMP 0.20

#define MAXTIME 20.0 // Maximum runtime.

#define NORMALDEPTH 0.001366334
#define NORMALVEL 0.2441864

double g_ = 9.81;

double xextent_ = 0.6118;

//----------------------MAIN-----------------------//
//Begin main function
int main()
{
  size (xextent_);
//   origin (0.0, -h_-yOffset_);
//   rho1 = 1.0;
//   rho2 = RATIO;
    rho1 = RATIO;
    rho2 = 1.0;
  /**
     For calculating viscosities, interpret Reynolds number as defined at depth of unity. */
//   mu1 = 17.4e-6;
//   mu2 = mu1*MURATIO;
    mu1 = 8.9e-4;
    mu2 = 17.4e-6;
  /**
     Use depth length scale for Bond number as well. */
  // do not consider surface tension for now
  f.sigma = 0.072;
  init_grid(1 << (9));
  /**
     Acceleration using reduced gravity. */
//   G.y = (-CHANNELCOS)*g_;
//   G.x = (CHANNELSLOPE)*g_;

  // body-force gravity
//   const face vector gravity[] = {(CHANNELSLOPE)*g_, (-CHANNELCOS)*g_, 0.0};
//   a = gravity;

  periodic (right);

  run();
}

//---------------------INITIALIZATION------------------------//
double distSurf(double xCoord)
{
  return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

/*Write initialization event*/
event init (i=0)
{
  if (!restore("restart")){
    double femax = 1e-4;
    double uemax = 2e-2;

    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      refine(y<(1.25*NORMALDEPTH) && level < MAXLEVEL);
      fraction (f, y - distSurf(x));

      foreach() {
        // variation of x-component velocity to keep discharge the same
    	u.x[] = y<=distSurf(x) ? (NORMALDEPTH*NORMALVEL/distSurf(x))*((1.0+2.0*POWERLAWN)/(1.0+POWERLAWN))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWN)/POWERLAWN)) : 0.0;
        u.y[] = 0.0;
      }
      boundary ((scalar *){u});
    }
    while (adapt_wavelet ({f,u}, (double[]){femax,uemax, uemax, uemax}, MAXLEVEL, MINLEVEL).nf);
  }
}

// max time
event maxdt (t <= MAXTIME; t += 0.50);

event timingLog(i += 10) {
  fprintf (stderr, "%d %g %g \n", i, t, dt);
}

//-------------------ADAPTIVITY---------------------//
/*Adapt once on error in volume fraction, velocity field, and beach fraction*/
event adapt(i++) {
  //double uemax = 1e-5;

//   double femax = 1e-3;
  double uemax = NORMALVEL/120.0;
//   adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
  adapt_wavelet_leave_interface ((scalar *){u},{f}, (double[]){uemax,uemax,uemax}, MAXLEVEL, MINLEVEL, 1);
}

/**
## Dump/restore

To be able to restart, we dump the entire simulation at regular
intervals. */

event snapshot (i += 100) {
  char nameOut[80];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut);
}

event outputGfsFiles (t = 0; t += 0.0025) {
    char name[80];
//     sprintf(name, "out-%.3f.vtk", t);
    sprintf(name, "out-%.3f.gfs", t);
    FILE *fp = fopen(name, "w");
//     output_vtk((scalar *) {h}, 1 << 10, (FILE *) fp, false);
    output_gfs(fp, translate = true);
    fclose (fp);
}

event outputInterface(t += 0.0025) {
  char resultname[40];
  sprintf( resultname, "surf_%g.txt", t );
  FILE * fp = fopen(resultname, "w");
  scalar xpos[];
  scalar ypos[];
  position (f, xpos, {1, 0});
  position (f, ypos, {0, 1});
  foreach()
    {
      if (xpos[] != nodata){
	fprintf (fp, "%g %g %g %g\n", x, y, xpos[], ypos[]);
      }
    }
  fclose (fp);
}

event vtk_file (t += 0.001) {
			FILE * fp;
			char name[40];
			sprintf(name, "outVTK-%g.vtu", t);
			fp = fopen(name,"w");
			output_vtu_bin_foreach((scalar *) {f}, (vector *) {u},N,fp,false);
			fclose(fp);
}

event waveprobes (i += 1) {
  heights (f,h);
  static FILE * fp0 = fopen("waveprobe.dat", "w");
  double xcoords[2]  = {0.20*NORMALDEPTH,3.0*NORMALDEPTH};
  double xMax0 = wprobe(xextent_/2.0,xcoords,72);
  fprintf(fp0, "%g %g \n", t, xMax0);
  fflush(fp0);
}
