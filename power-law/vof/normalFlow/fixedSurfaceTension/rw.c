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
// #include "./myCentered.h"
#include "navier-stokes/centered.h"
#include "./two-phasePL.h"
#include "tension.h"
// alternatively, use momentum-conserving scheme
// #include "navier-stokes/conserving.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"

#define FILTERED

/**
   Include profiling information. */

// #include "navier-stokes/perfs.h"
// #include "profiling.h"

/** Grid levels */

// #define MAXLEVEL 11
// #define MINLEVEL 3
#define INITLEVEL 5

/** Problem-related parameters */

#define FR 1.00 // Froude number for the expected steady-state solution, but this dimensionless parameter is not used in the codes here.
#define MUDRHO 1120.0 //density ratio, water to air
// #define MURATIO 8.9e-4/17.4e-6 //dynamic viscosity ratio, water to air
#define POWERLAWINDEX 0.40
#define MUMUD 0.140

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 0.002099552
#define NORMALVEL 0.1433858

#define WEBER 300.0
#define COEFFST (MUDRHO*pow(NORMALVEL,2.0)*NORMALDEPTH/WEBER)

// inclination angle of the channel. \sin\theta and \cos\theta
#define channelSlope 0.06
#define CHANNELCOS (pow((1.0-pow(channelSlope,2.0)),0.50))

#define MAXTIME 400.0 // Maximum runtime.
#define TOUTPUT 1.00

#define DISTAMP 0.20

#define grav 9.81

// square domain size
#define xextent_ (6.0*NORMALDEPTH/5.0)

#define KAPPAErr (1e-3)
#define OmegaErr (5.75e-2)
#define fErr (1e-7)
#define VelErr (8e-2)

/**
## Main body of the current codes
*/

/**
  slip at the top
*/
  u.t[top] = neumann(0.);
  u.n[top] = dirichlet(0.);
/**
 no slip at the bottom
*/
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

/** ### Main */
int main()
{
  size (xextent_);

  rho1 = MUDRHO;
  rho2 = 1.12;

  mu1 = MUMUD;
  mu2 = 0.001/50.0;

  powerLawIndex = POWERLAWINDEX;
  muRef = MUMUD;
  mumax = 175.0;

  f.sigma = COEFFST;

  // Surface tension seems not to change the solution too much, since there is very little interface curvature.
  // f.sigma = 0.072;
  init_grid(1 << (INITLEVEL));

  // periodic BC
  periodic (right);

  run();
}

//---------------------INITIALIZATION------------------------//
// double distSurf(double xCoord)
// {
//   return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
// }

/** ### Init event */
event init (i=0)
{
	if (!restore("restart")){

      fraction (f, NORMALDEPTH-y);

      foreach() {
        u.x[] = 0.0;
        u.y[] = 0.0;
      }
  }
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] += grav*channelSlope;
  foreach_face(y)
    av.y[] -= grav*CHANNELCOS;
}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= MAXTIME; t += 0.50);

event iterLog(i += 10) {
  FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (fp1, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
  fclose (fp1);
}

/**
### Dump output, Gfs file output, free-surface output*/

event snapshot (t += TOUTPUT) {
  char nameOut[50], nameOutText[50];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut);

  // text output
  sprintf(nameOutText, "slice-%g.txt", t);
  FILE *fp2 = fopen(nameOutText, "w");
  for (double yCoord = 0.0; yCoord <= xextent_; yCoord += (xextent_/ (pow(2, INITLEVEL))))
  {
    fprintf(fp2, "%g %g %g\n", yCoord, interpolate(u.x, (xextent_/2.0), yCoord), interpolate(f, (xextent_/2.0), yCoord));
  }
  fclose(fp2);
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
  // sprintf(names, "interface%d", pid());
  sprintf(names, "surf-%g.dat", t);
  FILE * fp = fopen (names, "w");
  output_facets (f,fp);
  fclose(fp);
}

// mesh adaptation
// event adapt (i++) {
// //   scalar KAPPA[], omega[];
//   scalar omega[];
// //   curvature(f, KAPPA);
//   vorticity (u, omega);
//   boundary ((scalar *){omega});
// //   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, omega}, (double[]){fErr, VelErr, VelErr, KAPPAErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   adapt_wavelet ((scalar *){f, u.x, u.y, omega}, (double[]){fErr, VelErr, VelErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
// //   adapt_wavelet ((scalar *){f, cs}, (double[]){fErr, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
// }

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

