/**
## Case set-up

  Steady-uniform flow down an incline.

  * VOF Navier-Stokes in 2D with surface tension.
  * Gravity is modelled as body force.
  * Periodic boundary condition for left-and-right boundaries.
  * 32x32 cells in 2D domain. No adaptation.
  * x-axis aligns with the inclined plane and y-axis is perpendicular to the inclined plane. Gravity $\pmb{g}$ is decomposited to two components: $(g_x, g_y)$.
  * Initial condition: 3/4 of the vertical domain (32x24) is set up as fluid phase and 1/4 of the vertical domain (32x8) is set up as air phase; velocity field is zero $u_x=u_y=0$ initially.
  * The steady-state numerical solution is to be compared with analytical solutions.
*/
/**
The analytical solutions for such a problem are:
$$$$
*/

// #include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phasePL.h"
#include "myTension.h"
// #include "vof.h"
// alternatively, use momentum-conserving scheme
// #include "navier-stokes/conserving.h"
// #include "./myTension.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"
#include "adapt_wavelet_limited.h"

#define FILTERED
// #define NFILTERED 2

/**
   Include profiling information. */

// #include "navier-stokes/perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 14
#define MINLEVEL 3
#define INITLEVEL 8

/** Problem-related parameters */

#define FR 1.00 // Froude number for the expected steady-state solution, but this dimensionless parameter is not used in the codes here.
#define MUDRHO 1120.0 //density ratio, water to air
// #define MURATIO 8.9e-4/17.4e-6 //dynamic viscosity ratio, water to air
#define POWERLAWINDEX 0.40
#define MUMUD 0.140

#define AIRRHO 1.12
#define AIRMU (0.001/50.0)

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 0.002099552
#define NORMALVEL 0.1433858

#define WEBER 200.0
#define COEFFST (MUDRHO*pow(NORMALVEL,2.0)*NORMALDEPTH/WEBER)

// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELSLOPE 0.06
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))

#define MAXTIME 40.0 // Maximum runtime.
#define TOUTPUT 0.20

#define DISTAMP 0.250
#define DISTPERIOD (2.0*NORMALDEPTH/CHANNELSLOPE/NORMALVEL)

#define GRAV 9.81

// square domain size
#define xextent_ (750.0*NORMALDEPTH)
#define topExtent (NORMALDEPTH*3.610)

#define KAPPAErr (1e-3)
#define OmegaErr (1.250)
#define fErr (1e-5)
#define VelErr (NORMALVEL/60.0)
#define KErr (0.80e-4)

/**
## Main body of the current codes
*/
double inletDistDepth(double timeVal)
{
  return NORMALDEPTH*(1.0+DISTAMP*sin(2*pi*timeVal/DISTPERIOD));
}

double inletDistVel(double yCoord, double timeVal)
{
  double distDepth = NORMALDEPTH*(1.0+DISTAMP*sin(2*pi*timeVal/DISTPERIOD));
  return (yCoord<=(distDepth) ? FR*sqrt(CHANNELCOS*GRAV*distDepth)*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-yCoord/distDepth), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : FR*sqrt(CHANNELCOS*GRAV*distDepth)*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*1.0);
}

/**
  slip at the top, dummy
*/
  u.t[top] = neumann(0.);
  u.n[top] = dirichlet(0.);
/**
  no slip at the bottom
*/
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
/**
  Real top BC control
 * **/
  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

  /*
  Inlet
  */
  f[left] = inletDistDepth(t);
  // f[left] = dirichlet(inletDistDepth(t));
  u.n[left] = dirichlet(inletDistVel(y, t));
  u.t[left] = dirichlet(0.0);
  p[left] = neumann(0);
  pf[left] = neumann(0);

  /*
  Outlet
  */
  f[right] = neumann(0.);
  u.n[right] = neumann(0.);
  u.t[right] = neumann(0.);
  p[right] = neumann(0);
  pf[right] = neumann(0);

// scalar boundaryLayer[];

/** ### Main */
int main()
{
  size (xextent_);

  rho1 = MUDRHO;
  rho2 = AIRRHO;

  mu1 = MUMUD;
  mu2 = AIRMU;

  powerLawIndex = POWERLAWINDEX;
  muRef = MUMUD;
  mumax = 175.0;

  f.sigma = COEFFST;

  // iteration specs 
  NITERMAX=250;
  // TOLERANCE=1e-4;

  // Surface tension seems not to change the solution too much, since there is very little interface curvature.
  // f.sigma = 0.072;
  init_grid(1 << (INITLEVEL));
  // Acceleration using reduced gravity. But reduced gravity approach does not work for this case.
//   G.y = (-CHANNELCOS)*g_;
//   G.x = (CHANNELSLOPE)*g_;

  /** Body-force gravity. This defines the acceleration vector $\pmb{a}$ in $\texttt{centered.h}$ file.*/
  const face vector gravity[] = {(CHANNELSLOPE)*GRAV, (-CHANNELCOS)*GRAV, 0.0};
//   const face vector gravity[] = {(CHANNELSLOPE)*GRAV*f, (-CHANNELCOS)*GRAV*f};
  a = gravity;

//   TOLERANCE = 1e-2;

  // periodic BC
  // periodic (right);

  run();
}

//---------------------INITIALIZATION------------------------//
double normalFlowVel(double yCoord)
{
  return (FR*sqrt(CHANNELCOS*GRAV*NORMALDEPTH)*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-yCoord/NORMALDEPTH), (1.0+POWERLAWINDEX)/POWERLAWINDEX)));
}

/** ### Init event */
event init (i=0)
{
	if (!restore("restart")){
    solid (cs, fs, topExtent - y );
    refine(y<=(22.50*NORMALDEPTH) && level < (MAXLEVEL-3));
    refine(y<=(5.00*NORMALDEPTH) && level < (MAXLEVEL-1));

    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      // refine(y<(4.00*NORMALDEPTH) && level < (MAXLEVEL-3));
      // refine(y<(2.00*NORMALDEPTH) && level < (MAXLEVEL-1));
      // refine(y<(1.50*NORMALDEPTH) && level < (MAXLEVEL));
      fraction (f, NORMALDEPTH-y);

      foreach() {
        // variation of x-component velocity to keep discharge the same
        u.x[] = (y<=topExtent*1.10) ? (y<=NORMALDEPTH ? normalFlowVel(y) : normalFlowVel(NORMALDEPTH)) : 0.0;
        // u.x[] = 0.0;
        u.y[] = 0.0;
        // hydrostatic pressure, zero pressure datum at free-surface
        p[] = (y<=NORMALDEPTH) ? MUDRHO*CHANNELCOS*GRAV*(NORMALDEPTH-y) : (-1.0)*CHANNELCOS*GRAV*AIRRHO*(y-NORMALDEPTH);

        // boundaryLayer[] = (y<=xextent_/pow(2,MAXLEVEL)*4.0 ? 1.0 : 0.0);
      }
      boundary ((scalar *){u});
    }
    while (adapt_wavelet ((scalar *){f, cs, u.x, u.y}, (double[]){(fErr/4.0), fErr, VelErr, VelErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);
    unrefine(y>topExtent*1.20 && level>MAXLEVEL-5);
//     jdx = 0;
//     do {
//       jdx += 1;
//       fprintf(ferr, "Refining solid, %d\n", jdx );
//     }
//     while (adapt_wavelet ({cs}, (double[]){0.005}, MAXLEVEL, MINLEVEL).nf);
//     mask (y > (NORMALDEPTH*4.01) ? top : none);
  }
}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= MAXTIME; t += 0.50);

event drop_remove (i += 1) {
  remove_droplets (f, 2, 1e-4, false);
  remove_droplets (f, 2, 1e-4, true);

  // remove_droplets (f, 1, 1e-4, false);
  // remove_droplets (f, 1, 1e-4, true);
}

// event timingLog(i += 10) {
//   fprintf (stderr, "%d %g %g \n", i, t, dt);
//   fflush (stderr);
// }

event iterLog(i += 10) {
  FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (fp1, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
  fclose (fp1);
}

/** ### No AMR is used for now */
//-------------------ADAPTIVITY---------------------//
/*Adapt once on error in volume fraction, velocity field, and beach fraction*/
// event adapt(i++) {
//   //double uemax = 1e-5;
//
// //   double femax = 1e-3;
//   double uemax = NORMALVEL/120.0;
// //   adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
//   adapt_wavelet_leave_interface ((scalar *){u},{f}, (double[]){uemax,uemax,uemax}, MAXLEVEL, MINLEVEL, 1);
// }

/**
### Dump output, Gfs file output, free-surface output*/

event snapshot (t += TOUTPUT) {
  char nameOut[50];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut);

  // text output
  // sprintf(nameOutText, "slice-%g.txt", t);
  // FILE *fp2 = fopen(nameOutText, "w");
  // for (double yCoord = 0.0; yCoord <= xextent_; yCoord += (xextent_/ (pow(2, INITLEVEL))))
  // {
  //   fprintf(fp2, "%g %g %g\n", yCoord, interpolate(u.x, (xextent_/2.0), yCoord), interpolate(f, (xextent_/2.0), yCoord));
  // }
  // fclose(fp2);
}

event snapshotCheck (t=0; t += 0.025; t<=0.125) {
  char nameOut[50];
  sprintf (nameOut, "checkSnapshot-%g", t);
  dump(file=nameOut);
}

event outputGfsFiles (t += TOUTPUT) {
    char name[80];
    sprintf(name, "out-%g.gfs", t);
    FILE *fp = fopen(name, "w");
    output_gfs(fp, translate = true);
    fclose (fp);
}

event outputInterface(t += TOUTPUT) {
  char names[36];
  sprintf( names, "interfaceMed-%d.dat", pid() );
  FILE * fp = fopen (names, "w");
  output_facets (f,fp);
  fclose(fp);
  char command[80];
  sprintf(command, "LC_ALL=C  cat interface* > ALLINTER-%g.dat",t);
  system(command);
}

int refRegion(double x, double y, double z){
  return (y<=3.10*NORMALDEPTH ? MAXLEVEL : (MAXLEVEL-2));
}

// mesh adaptation
event adapt (i++) {
  // scalar KAPPA[];
  scalar omega[], KAPPA[];
  curvature(f, KAPPA);
  vorticity (u, omega);
  // boundary ((scalar *){omega});
//   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, omega}, (double[]){fErr, VelErr, VelErr, KAPPAErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // adapt_wavelet ((scalar *){f, u.x, u.y, omega, KAPPA}, (double[]){fErr, VelErr, VelErr, OmegaErr, KErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  adapt_wavelet_limited ((scalar *){f, omega,KAPPA}, (double[]){fErr, OmegaErr,KErr}, refRegion, minlevel = MINLEVEL);
  refine(y<=xextent_/pow(2,MAXLEVEL)*5.0 && level<MAXLEVEL);
//   adapt_wavelet ((scalar *){f, cs}, (double[]){fErr, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

/**## Numerical Results */
/**
 <span style="color:red"> **The fluid stays stationary ($u_x(t)=u_y(t)=0$), which is apparently different from the analytical solutions and it is not reasonable at all. Why??** </span>*/
// event movies (i += 50) {
//   view (quat = {0.000, 0.000, 0.000, 1.000},
//       fov = 30, near = 0.01, far = 1000,
//       tx = -0.477, ty = -0.040, tz = -1.116,
//       width = 1620, height = 873);
//   box ();
//   squares (color = "u.x", max = 4.5, linear = true);
// //   draw_vof (c = "cs", fc = {0.7843137254901961,0,0}, lc = {0.7843137254901961,0,0}, lw = 3.0);
//   draw_vof (c = "cs", fc = {0.7843137254901961,0,0}, lc = {0.7843137254901961,0,0});
//   cells ();
// //   draw_vof (c = "f", lw = 2.0);
//   draw_vof (c = "f");
//   save ("movie.mp4");
// }

