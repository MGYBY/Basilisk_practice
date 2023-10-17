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
#include "./adapt_wavelet_leave_interface_limited.h"
#include "embed.h"
#include "navier-stokes/centered.h"
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

#define FILTERED

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 12
#define MINLEVEL 3
#define INITLEVEL 10

/** Problem-related parameters */

#define MUDRHO 1437.25 //density ratio, water to air
// #define MURATIO 8.9e-4/17.4e-6 //dynamic viscosity ratio, water to air
#define POWERLAWINDEX 0.140
#define MUMUD 74.60
#define YIELDSTRESS 0.00

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 0.12986706
#define NORMALVEL 0.22553941

#define AIRRHO (MUDRHO/850.0) // 1.12
#define AIRMU ((MUMUD*pow((NORMALVEL/NORMALDEPTH), (POWERLAWINDEX-1.0))+YIELDSTRESS/(NORMALVEL/NORMALDEPTH))/2.0e2)

// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELSLOPE 0.060 // sinTheta
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))
#define CHANNELTAN (CHANNELSLOPE/CHANNELCOS)
#define GRAV 9.81
#define GRAVRED (GRAV*CHANNELCOS)

#define FR (NORMALVEL/pow((GRAVRED*NORMALDEPTH), 0.50)) // Froude number for the expected steady-state solution, but this dimensionless parameter is not used in the codes here.

#define WEBER 250.0
#define COEFFST (MUDRHO*pow(NORMALVEL,2.0)*NORMALDEPTH/WEBER)

#define MAXTIME 200.0 // Maximum runtime.
#define TOUTPUT 0.250

#define DISTAMP 0.20

// square domain size
#define xextent_ (6.0*NORMALDEPTH/CHANNELTAN)
#define topExtent (NORMALDEPTH*3.600)
#define initTransitCells 6

// #define KAPPAErr (1e-3)
#define OmegaErr (2.00)
#define fErr (1e-7)
#define VelErr (NORMALVEL/50.0)
#define KErr (1e-3)

// #define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
// harmonic mean
#define mu(muTemp, mu2, f)  (1./(clamp(f,0,1)/muTemp  + (1.-clamp(f,0,1))/mu2))

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

scalar volDroplet[];

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
  mumax = 3200.0;

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

  NITERMAX = 256;
  TOLERANCE = 2.00e-4;
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
    solid (cs, fs, topExtent - y );
    int jdx = 0;
//     do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      refine(y<(7.0*NORMALDEPTH) && level < (MAXLEVEL-3));
      refine(y<(5.0*NORMALDEPTH) && level < (MAXLEVEL-2));
      refine(y<(2.50*NORMALDEPTH) && level < MAXLEVEL);
      refine(y>(topExtent-2.5*xextent_/pow(2, INITLEVEL)) && y<(topExtent+2.5*xextent_/pow(2, INITLEVEL)) && level < MAXLEVEL);
      fraction (f, distSurf(x)-y);

      foreach() {
        // variation of x-component velocity to keep discharge the same
        // u.x[] = y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX));
        // u.x[] = (y<=topExtent*1.10) ? y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX)) : 0.0;
            u.x[] = (y<=topExtent*1.0) ? y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (y<=(distSurf(x)+initTransitCells*xextent_/pow(2, MAXLEVEL))) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(((distSurf(x)+initTransitCells*xextent_/pow(2, MAXLEVEL))-y+xextent_/pow(2, MAXLEVEL))/(initTransitCells*xextent_/pow(2, MAXLEVEL))) : 0.0 : 0.0;
  // u.x[] = 0.0;
        u.y[] = 0.0;
        // hydrostatic pressure, zero pressure datum at free-surface
        p[] = (y<=distSurf(x)) ? MUDRHO*CHANNELCOS*GRAV*(distSurf(x)-y) : (-1.0)*CHANNELCOS*GRAV*AIRRHO*(y-distSurf(x));
      }
      boundary ((scalar *){u});
//     }
//     while (adapt_wavelet ((scalar *){f, u.x, u.y}, (double[]){fErr, VelErr, VelErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);
    // avoid excessively large init file.
    // unrefine(y>topExtent*1.10 && level>MAXLEVEL-4);

    

//     jdx = 0;
//     do {
//       jdx += 1;
//       fprintf(ferr, "Refining solid, %d\n", jdx );
//     }
//     while (adapt_wavelet ({cs}, (double[]){0.005}, MAXLEVEL, MINLEVEL).nf);
//     mask (y > (NORMALDEPTH*4.01) ? top : none);
//   }
}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= MAXTIME; t += 0.50);

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] += GRAV*CHANNELSLOPE;
  foreach_face(y)
    av.y[] -= GRAV*CHANNELCOS;
}

// event drop_remove (i += 1) {
//   remove_droplets (f, 2, 1e-4, false);
//   remove_droplets (f, 2, 1e-4, true);
// }

// event timingLog(i += 10) {
//   fprintf (stderr, "%d %g %g \n", i, t, dt);
//   fflush (stderr);
// }

event iterLog(i += 5) {
  // FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
  // fclose (fp1);
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
  char nameOut[50], nameOutText[50];
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
  system(command);// allow to use linux command in the c code to concatenate our files

//   char resultname[32];
//   sprintf( resultname, "interface-%g.dat", t );
//   FILE * fp = fopen(resultname, "w");
//   foreach(serial)
//     {
//       scalar xpos[];
//       scalar ypos[];
//       position (f, xpos, {1, 0});
//       position (f, ypos, {0, 1});
//       if (xpos[] != nodata){
// // 	fprintf (fp, "%g %g %g %g\n", x, y, xpos[], ypos[]);
//         fprintf (fp, "%g %g \n", xpos[], ypos[]);
//         // fclose (fp);
//       }
//     }
//   fflush(fp);
}

// event outputCentVel(t += TOUTPUT) {
//   char resultname[40];
//   sprintf( resultname, "centVel_%g.txt", t );
//   FILE * fp = fopen(resultname, "w");
//   for (double y = 0.; y < xextent_; y += xextent_/pow(2.,LEVEL))
//         fprintf (fp, "%g %g %g %g \n", y, interpolate (u.x, xextent_/2, y), interpolate (u.y, xextent_/2, y), interpolate (f, xextent_/2, y));
//   fclose (fp);
// }

event dropletsAndhMax (i += 20)
{
  FILE *fp1 = fopen("totalDroplets", "a+");
  FILE *fp2 = fopen("dropletDetail", "a+");

  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach (serial)* to avoid doing a parallel traversal when
  using OpenMP. This is because we don't have reduction operations for
  the *v* and *b* arrays (yet). */

  double v[n];
//   coord b[n];
  for (int j = 0; j < n; j++)
//     v[j] = b[j].x = b[j].y = b[j].z = 0.;
    v[j] = 0.;
  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
//       coord p = {x,y,z};
//       foreach_dimension()
// 	b[j].x += dv()*f[]*p.x;
    }

  foreach(serial)
  {
    volDroplet[] = 0.0;
    if (m[] > 0) {
      int j = m[] - 1;
      volDroplet[] = v[j];
//       coord p = {x,y,z};
//       foreach_dimension()
// 	b[j].x += dv()*f[]*p.x;
    }
  }

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */

  fprintf (fp1, "%g %d \n", t, n);
  for (int j = 0; j < n; j++)
    fprintf (fp2, "%g %d %g \n", t, j, v[j]);
  fflush (fp1);
  fflush (fp2);
}

event depthAmplitude (i += 20) {
  FILE *fp2 = fopen("amplitude", "a+");
  double ampY = 0.0;
  foreach (reduction(max:ampY))
  {
    // double yCoordInt = volDroplet[]>(16.50*sq(xextent_/pow(2, MAXLEVEL))) ? y : 0.0;
    double yCoordInt = volDroplet[]>(16.10*sq(xextent_/pow(2, MAXLEVEL))) ? y+0.50*(xextent_/pow(2, MAXLEVEL))-(1-f[])*(xextent_/pow(2, MAXLEVEL)) : 0.0;
    if (ampY<yCoordInt)
    {
      ampY = yCoordInt;
    }
  }
  fprintf (fp2, "%g %g \n", t, ampY);
  fprintf (ferr, "%g %g \n", t, ampY);
  fflush (fp2);
}

int refRegion(double x,double y, double z){
    int lev;
    if( y < topExtent )
      lev = MAXLEVEL;
    else
      lev = MINLEVEL+1;

    return lev;
}

// mesh adaptation
event adapt (i++) {
  scalar omega[], KAPPA[];
  curvature(f, KAPPA);
  vorticity (u, omega);
  boundary ((scalar *){omega});
//   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, omega}, (double[]){fErr, VelErr, VelErr, KAPPAErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  adapt_wavelet_limited ((scalar *){ u.x, u.y, omega, KAPPA}, {f}, (double[]){ VelErr, VelErr, OmegaErr, KErr}, refRegion, minlevel = MINLEVEL);
  refine(y<=(15.510*xextent_/pow(2, MAXLEVEL)) && level<MAXLEVEL);
//   adapt_wavelet ((scalar *){f, cs}, (double[]){fErr, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

/**## Numerical Results */
/**
 <span style="color:red"> **The fluid stays stationary ($u_x(t)=u_y(t)=0$), which is apparently different from the analytical solutions and it is not reasonable at all. Why??** </span>*/
event movies (i += 25) {
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.493, ty = -0.012, tz = -1.106,
      width = 1200, height = 600);
  box ();
  // cells ();
  squares (color = "u.x");
  // vectors (u = "u", scale = 0.00005);
  draw_vof (c = "f", lw = 2.4, fc = {1.0, 0.0, 1.0});
  save ("movie.mp4");
}

