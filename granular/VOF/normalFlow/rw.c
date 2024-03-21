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

#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
// #define FILTERED
#include "two-phasePL.h"
// #include "./myTension.h"
#include "tension.h"
// #include "vof.h"
// alternatively, use momentum-conserving scheme
#include "navier-stokes/conserving.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 7
#define MINLEVEL 3
#define INITLEVEL 7

/** Problem-related parameters */
#define MU1 0.381863
#define MU2 0.643468
#define MU 0.554309 // Tan(Theta)
#define PHIVAL 0.60
#define NORMALDEPTHDIM 0.00312955
#define NORMALVELDIM 0.16899
#define DGDIM 0.50e-3
#define DGVAL (DGDIM/NORMALDEPTHDIM)
#define I0VAL 0.26602 // tuned I0 value

#define MUDRHO 1.00 //density ratio, water to air
#define rhoR 1000.10
#define muR 50.10

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 1.00
#define NORMALVEL 1.00

#define FR 1.02
#define MUMUD (1.0)

#define MUMAXVAL 50.0

// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELTAN MU
#define CHANNELSLOPE (pow((CHANNELTAN*CHANNELTAN/(1.0+CHANNELTAN*CHANNELTAN)), 0.50)) // sinTheta
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))
#define GRAV 1.0
#define GRAVRED (GRAV*CHANNELCOS)

#define AIRRHO (MUDRHO/rhoR) // 1.12
// #define AIRMU (MUMUD/48.0)
#define AIRMU ((MU1+(MU2-MU1)/(I0VAL/(DGVAL/pow(fabs(PHIVAL),0.50))+1.0))/muR)

#define WEBER 250.0
#define COEFFST (MUDRHO*pow(NORMALVEL,2.0)*NORMALDEPTH/WEBER)

#define MAXTIME 600.0 // Maximum runtime.
#define TOUTPUT 5.0

#define DISTAMP 0.0

// square domain size
#define xextent_ (NORMALDEPTH*6.0/5.0)
#define topExtent (NORMALDEPTH*3.500)

#define KAPPAErr (1e-3)
#define OmegaErr (1.20)
#define fErr (1e-7)
#define VelErr (NORMALVEL/64.0)
#define KErr (1e-4)

// #define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#define mu(muTemp, mu2, f)  (1./(clamp(f,0,1)/muTemp  + (1.-clamp(f,0,1))/mu2))

/**
## Main body of the current codes
*/

/**
  slip at the top
*/
  u.t[top] = neumann(0.);
  u.n[top] = dirichlet(0.);
  p[top] = dirichlet(-1.0*AIRRHO*(xextent_-NORMALDEPTH)*(1.0/(FR*FR)));
  pf[top] = dirichlet(-1.0*AIRRHO*(xextent_-NORMALDEPTH)*(1.0/(FR*FR)));

/**
 no slip at the bottom
*/
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

scalar fluidVel[];
int j = 0;

/** ### Main */
int main()
{
  size (xextent_);

  rho1 = MUDRHO;
  rho2 = AIRRHO;

  mu1 = MUMUD;
  mu2 = AIRMU;

  froude = FR;
  mus = MU1;
  deltaMu = MU2-MU1;
  dg = DGVAL;
  phiVF = PHIVAL;
  I0 = I0VAL;
  mumax = MUMAXVAL;

  f.sigma = COEFFST;

  // Surface tension seems not to change the solution too much, since there is very little interface curvature.
  // f.sigma = 0.072;
  init_grid(1 << (INITLEVEL));
  // Acceleration using reduced gravity. But reduced gravity approach does not work for this case.
//   G.y = (-CHANNELCOS)*g_;
//   G.x = (CHANNELSLOPE)*g_;

  /** Body-force gravity. This defines the acceleration vector $\pmb{a}$ in $\texttt{centered.h}$ file.*/
//   const face vector gravity[] = {(CHANNELSLOPE)*GRAV, (-CHANNELCOS)*GRAV, 0.0};
//   const face vector gravity[] = {(CHANNELSLOPE)*GRAV*f, (-CHANNELCOS)*GRAV*f};
//   a = gravity;

//   TOLERANCE = 1e-2;

  // periodic BC
  periodic (right);

  NITERMAX = 200;
  TOLERANCE = 2.50e-4;
  CFL = 0.455;

  run();
}

//---------------------INITIALIZATION------------------------//
double distSurf(double xCoord)
{
  return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

/** ### Init event */
event init (i=0)
{
//     if (!restore("restart")){
    // double femax = 1e-4;
    // double uemax = 2e-2;
//     solid (cs, fs, topExtent - y );
    int jdx = 0;
//     do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      fraction (f, distSurf(x)-y);

      foreach() {
        // variation of x-component velocity to keep discharge the same
        // u.x[] = y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX));
//         u.x[] = (y<=topExtent*1.10) ? y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX)) : 0.0;
        u.x[] = 0.0;
        // u.x[] = 0.0;
        u.y[] = 0.0;
        // hydrostatic pressure, zero pressure datum at free-surface
        p[] = (y<=distSurf(x)) ? (1.0/(FR*FR))*(distSurf(x)-y) : (-1.0)*AIRRHO*(1.0/(FR*FR))*(y-distSurf(x));
        fluidVel[] = f[]*u.x[];
      }
      boundary ((scalar *){u});
//     }
//     while (adapt_wavelet ((scalar *){f, u.x, u.y}, (double[]){fErr, VelErr, VelErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);
    // avoid excessively large init file.
    // unrefine(y>topExtent*1.10 && level>MAXLEVEL-4);

}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= MAXTIME; t += 0.50);

event acceleration (i++) {
  face vector av = a;
//   foreach_face(x)
//     av.x[] += CHANNELTAN*(1.0/(FR*FR))*f[];
//   foreach_face(y)
//     av.y[] -= 1.00*(1.0/(FR*FR))*f[];

  foreach_face(x)
    av.x[] += CHANNELTAN*(1.0/(FR*FR))*f[];
  foreach_face(y)
    av.y[] -= 1.00*(1.0/(FR*FR));

  foreach()
  {
    fluidVel[] = f[]*u.x[];
  }
}

event iterLog(i += 5) {
  // FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
  // fclose (fp1);
}

/**
### Dump output, Gfs file output, free-surface output*/

event snapshot (t += TOUTPUT) {
  char nameOut[50], nameOutText[50];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut);
}

event outputField (t += TOUTPUT) {
    char nameField[40], nameSlice[40];
    sprintf(nameField, "field-%g.txt", t);
//     sprintf(nameSlice, "slice-%g.txt", t);
    FILE *fp1 = fopen(nameField, "w");
//     FILE *fp2 = fopen(nameSlice, "w");
    output_field({f,u.x}, fp1, linear=true);
//     for (j=1; j<30; j++)
//       fprintf(fp2, "%g %g %g %g\n", (xextent_/2.0), (j*1.1*NORMALDEPTH/30.0), interpolate(u.x, (xextent_/2.0), (j*1.1*NORMALDEPTH/30.0), 0.0), interpolate(f, (xextent_/2.0), (j*1.1*NORMALDEPTH/30.0), 0.0));
    fclose (fp1);
//     fclose (fp2);
}

event globalMaxFluidVel(i+=25)
{
  FILE *fp2 = fopen("maxVel", "a+");
//   foreach()
//   {
//     double maxVelVal = 0.0;
//     if (f[]*u.x[]>maxVelVal)
//       maxVelVal = f[]*u.x[];
//   }
  stats s1 = statsf (fluidVel);
  fprintf(fp2, "%g %g \n", t, s1.max);
  fclose(fp2);
}

event outputInterface(t += TOUTPUT) {
  char names[36];
  sprintf( names, "interfaceMed-%d.dat", pid() );
  FILE * fp = fopen (names, "w");
  output_facets (f,fp);
  fclose(fp);
  char command[80];
  sprintf(command, "LC_ALL=C  cat interface* > ALLINTER-%g.dat",t);
  system(command);// allow to use linux command in the c code to concatenate our files

}

event outputGfsFiles (t += TOUTPUT) {
// event outputGfsFiles (i+=20) {
    char name[40];
    sprintf(name, "out-%g.gfs", t);
    FILE *fp1 = fopen(name, "w");
    output_gfs(fp1, translate = true, list={f, u.x, u.y, uf.x, uf.y, p});
    fclose (fp1);
}

/**## Numerical Results */
/**
 <span style="color:red"> **The fluid stays stationary ($u_x(t)=u_y(t)=0$), which is apparently different from the analytical solutions and it is not reasonable at all. Why??** </span>*/
// event movies (i += 5) {
//   view (quat = {0.000, 0.000, 0.000, 1.000},
//       fov = 30, near = 0.01, far = 1000,
//       tx = -0.493, ty = -0.012, tz = -1.106,
//       width = 1130, height = 618);
//   box ();
//   // cells ();
//   squares (color = "u.x");
//   // vectors (u = "u", scale = 0.00005);
//   draw_vof (c = "f", lw = 2.4, fc = {1.0, 0.0, 1.0});
//   save ("movie.mp4");
// }

