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
// #include "./adapt_wavelet_leave_interface_limited.h"
#include "./adapt_wavelet_limited_v4.h"
// #include "embed.h"
#include "navier-stokes/centered.h"
#define FILTERED
// not sure what is the harmonic mean in the three-phase model
// #define mu(f)  (clamp(f3,0.,1.)*rho3 + clamp(f2,0.,1.)*rho2 + clamp(f1,0.,1.)*rho1)
// #define mu(f)  (1.0/(clamp(f3,0.,1.)*rho3 + clamp(f2,0.,1.)*rho2 + clamp(f1,0.,1.)*rho1))
#include "./three-phase.h"
// #include "./conserving3f_tracer.h"
#include "./conserving3f.h"
#include "tension.h"
// #include "vof.h"
// alternatively, use momentum-conserving scheme
// #include "navier-stokes/conserving.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"

/**
   Include profiling information. */

#include "./my-perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 11
#define MINLEVEL 3
#define INITLEVEL 9

/** Problem-related parameters */
// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELTAN 0.06
#define CHANNELSLOPE (pow((CHANNELTAN*CHANNELTAN/(1.0+CHANNELTAN*CHANNELTAN)), 0.50)) // sinTheta
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))
// #define CHANNELTAN (CHANNELSLOPE/CHANNELCOS)
#define GRAV 1.0
#define GRAVRED (GRAV*CHANNELCOS)

// two-layer parameters
#define muR 0.33333
#define rhoR 0.9868
#define hR 1.0
#define FR1 0.60
#define Re1 (FR1*FR1/CHANNELTAN/((2.0+3.0*hR*rhoR)/6.0))

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH1 1.00 // layer 1's normal depth
#define NORMALVEL1 1.00 // layer 1's normal vel

// densities and viscosities
// layer 1: f3; layer 2: f2; air: f1.
#define MUDRHO 1.0 // layer 1's density
#define MUDMU (1.0/Re1) // layer 1's viscosity
#define FLUIDRHO (1.0*rhoR) // layer 2's density
#define FLUIDMU (1.0/Re1*muR) // layer 2's viscosity

#define rRho 100.0 // ratio between layer 2 and air
#define rMu 50.0 // ratio between layer 2 and air

#define AIRRHO (FLUIDRHO/rhoR)
// #define AIRMU (MUMUD/50.0)
#define AIRMU (FLUIDMU/rMu)

#define WEBER 150.0 // defined based on layer 2
#define COEFFST ((MUDRHO*muR)*pow(NORMALVEL1,2.0)*NORMALDEPTH1/WEBER)

// #define MAXTIME 200.0 // Maximum runtime.
#define DIMLESSTOUTPUT 0.10
#define TOUTPUT (DIMLESSTOUTPUT*NORMALDEPTH1/CHANNELTAN/NORMALVEL1)
#define MAXTIME (1250.01*TOUTPUT) // Maximum runtime.

#define DISTAMP 0.10

// square domain size
#define xextent_ (3.0*NORMALDEPTH1/CHANNELTAN)
#define topExtent (NORMALDEPTH1*(1.0+hR)*3.01)
#define initTransitCells 30

// #define KAPPAErr (1e-3)
#define OmegaErr (NORMALVEL1/NORMALDEPTH1/12.50)
#define fErr (1e-7)
#define VelFluidErr (NORMALVEL1/45.50)
#define VelAirErr (NORMALVEL1/12.50)
#define KErr (1e-3)

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

//   u.n[embed] = dirichlet(0.);
//   // u.t[embed] = neumann(0.);
//   u.t[embed] = dirichlet(0.);

scalar volDroplet[];
scalar omega[];
scalar velFluidNorm[], velAirNorm[];
scalar ypos[];
// scalar T[];

/** ### Main */
int main()
{
  size (xextent_);

  rho3 = MUDRHO;
  rho2 = FLUIDRHO;
  rho1 = AIRRHO;

  mu3 = MUDMU;
  mu2 = FLUIDMU;
  mu1 = AIRMU;

  f1.sigma = COEFFST;
  f2.sigma = COEFFST;
  f3.sigma = COEFFST;

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

  NITERMAX = 255;
  TOLERANCE = 9.999e-5;
  CFL = 0.4750;

  run();
}

//---------------------INITIALIZATION------------------------//
double distSurf(double xCoord, double hn)
{
  return (hn*(1.0+DISTAMP*sin(2.0*pi*xCoord/(xextent_/1.0))));
}

double velSU(double xCoord, double yCoord, double hTotal)
{
  // only for the fluid & mud phases
  double term, interY;
  term = 2.0+3.0*hR*rhoR;
  interY = hR*(hTotal*(1.0+DISTAMP*sin(2.0*pi*xCoord/(xextent_/1.0)))); // essentially "distSurf@ND1"*hR
  if (yCoord<=interY)
    return ((3.0*yCoord*(2.0+2.0*hR*rhoR))/term-3.0*(yCoord*yCoord)/term);
  else
    return ((3.0*muR-3.0*rhoR-6.0*hR*rhoR+hR*muR*rhoR)/(term*muR)+yCoord*(6.0*rhoR+6.0*hR*rhoR)/(term*muR)-3.0*(yCoord*yCoord)*(rhoR/(term*muR)));
}

/** ### Init event */
event init (i=0)
{
    if (!restore("dumpSnapshot-6.25")){
    // double femax = 1e-4;
    // double uemax = 2e-2;

    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
//       refine(y<(7.0*NORMALDEPTH) && level < (MAXLEVEL-3));
      refine(y<(9.0*NORMALDEPTH1*(1.0+hR)) && level < (MAXLEVEL-1));
      refine(y<(2.50*NORMALDEPTH1*(1.0+hR)) && level < MAXLEVEL);
      refine(y>(topExtent-1.51*xextent_/pow(2, INITLEVEL)) && y<(topExtent+1.51*xextent_/pow(2, INITLEVEL)) && level < MAXLEVEL);

      // fraction (f, distSurf(x)-y);
      fraction (f1, -distSurf(x, (NORMALDEPTH1*(1.0+hR)))+y);
      fraction (f2, distSurf(x, (NORMALDEPTH1*(1.0+hR)))-y);
      fraction (f3, distSurf(x, NORMALDEPTH1)-y);

//       solid (cs, fs, topExtent - y );

      foreach() {
        // variation of x-component velocity to keep discharge the same
        // u.x[] = y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX));
        // u.x[] = (y<=topExtent*1.10) ? y<=distSurf(x) ? (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-y/distSurf(x)), (1.0+POWERLAWINDEX)/POWERLAWINDEX)) : (FR*sqrt(CHANNELCOS*GRAV*distSurf(x)))*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX)) : 0.0;
        u.x[] = (y<=topExtent*1.0) ? velSU(x, y, NORMALDEPTH1) : 1.14*velSU(x, NORMALDEPTH1*(1.0+hR), NORMALDEPTH1);
        // u.x[] = 0.0;
        u.y[] = 0.0;
        // hydrostatic pressure, zero pressure datum at free-surface
        p[] = (y<=distSurf(x,NORMALDEPTH1)) ? (1.0/(FR1*FR1))*FLUIDRHO*(distSurf(x,(NORMALDEPTH1*hR)))+(1.0/(FR1*FR1))*MUDRHO*(distSurf(x,NORMALDEPTH1)-y) : ((y<(distSurf(x,NORMALDEPTH1)+distSurf(x,(NORMALDEPTH1*hR)))) ? (1.0/(FR1*FR1))*FLUIDRHO*(distSurf(x,NORMALDEPTH1)+distSurf(x,(NORMALDEPTH1*hR))-y) : (-1.0)*AIRRHO*(1.0/(FR1*FR1))*(y-(distSurf(x, NORMALDEPTH1)+distSurf(x,(NORMALDEPTH1*hR)))));

        // T[] = f1[];

        velFluidNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*(f2[]+f3[]);
        velAirNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*(f1[]);
      }
      // boundary ((scalar *){u, p, T});
      boundary ((scalar *){u, p});
      vorticity (u, omega);
    }
    while (adapt_wavelet ((scalar *){f1, f2, f3, u.x, u.y, omega}, (double[]){(fErr/1.0), (fErr/1.0), (fErr/1.0), VelFluidErr/5.0, VelFluidErr/5.0, OmegaErr/5.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);
    // avoid excessively large init file.
    // unrefine(y>topExtent*1.10 && level>MAXLEVEL-4);

//     fractions_cleanup (cs, fs);

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
event maxdt (t <= MAXTIME; t += 50.0);

event acceleration (i++) {
  face vector av = a;
//   foreach_face(x)
//     av.x[] += CHANNELTAN*(1.0/(FR*FR))*f[];
//   foreach_face(y)
//     av.y[] -= 1.00*(1.0/(FR*FR))*f[];
  
    foreach_face(x)
    av.x[] += CHANNELTAN*(1.0/(FR1*FR1))*(f2[]+f3[]);
  foreach_face(y)
    av.y[] -= 1.00*(1.0/(FR1*FR1));

  foreach()
  {
    velFluidNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*(f2[]+f3[]);
    velAirNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*(f1[]);
  }
}

// static scalar * interfaces1 = NULL;
//
// event vof (i++) {
//   /**
//    *   We allocate three temporary vector fields to store the three components
//    *   of the momentum and set the boundary conditions and prolongation
//    *   functions. */
//
//   vector q1[], q2[], q3[];
//   for (scalar s in {q1,q2,q3})
//     foreach_dimension()
//       s.v.x.i = -1; // not a vector
//       for (int i = 0; i < nboundary; i++)
//         foreach_dimension() {
//           q1.x.boundary[i] = boundary_q1_x;
//           q2.x.boundary[i] = boundary_q2_x;
//           q3.x.boundary[i] = boundary_q3_x;
//         }
//         #if TREE
//         foreach_dimension() {
//           q1.x.prolongation = prolongation_q1_x;
//           q2.x.prolongation = prolongation_q2_x;
//           q3.x.prolongation = prolongation_q3_x;
//         }
//         #endif
//
//         /**
//          *   We split the total momentum $q$ into its three components $q2$,$q3$ and $q1$
//          *   associated with $f2$, $f3$ and $f1$ respectively. */
//
//         foreach()
//           foreach_dimension() {
//             double fc1 = clamp(f1[],0,1);
//             double fc2 = clamp(f2[],0,1);
//             double fc3 = clamp(f3[],0,1);
//             q2.x[] = fc2*rho2*u.x[];
//             q3.x[] = fc3*rho3*u.x[];
//             q1.x[] = fc1*rho1*u.x[];
//           }
//           boundary ((scalar *){q1,q2,q3});
//
//         /**
//          *   We use the same slope-limiting as for the
//          *   velocity field. */
//
//         foreach_dimension() {
//           q1.x.gradient = q2.x.gradient = q3.x.gradient = u.x.gradient;
//         }
//
//         /**
//          *   Momentum $C2$ is associated with $1 - f$, so we set the *inverse*
//          *   attribute to *true*. We use the same slope-limiting as for the
//          *   velocity field.
//          *   We associate the transport of $q1$ and $q2$ with $f$ and transport
//          *   all fields consistently using the VOF scheme. */
//
//         f1.tracers = (scalar *){q1,T};
//         f2.tracers = (scalar *){q2};
//         f3.tracers = (scalar *){q3};
//         vof_advection ({f1,f2,f3},i);
//
//         /**
//          *   We recover the advected velocity field using the total momentum and
//          *   the density */
//
//         foreach()
//           foreach_dimension()
//             u.x[] = (q1.x[] + q2.x[]+ q3.x[])/rho(f1[],f2[],f3[]);
//           boundary ((scalar *){u});
//
//         /**
//          *   We set the list of interfaces to NULL so that the default *vof()*
//          *   event does nothing (otherwise we would transport $f$ twice). */
//
//         interfaces1 = interfaces, interfaces = NULL;
// }
//
// /**
//  * We set the list of interfaces back to its default value. */
//
// event tracer_advection (i++) {
//   interfaces = interfaces1;
// }

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

event snapshot (t += TOUTPUT*5.0) {
  char nameOut[50];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut, list={f1, f2, f3, u.x, u.y, uf.x, uf.y, p, omega});

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
    FILE *fp1 = fopen(name, "w");
    output_gfs(fp1, translate = true, list={f1, f2, f3, u.x, u.y, uf.x, uf.y, p, omega});
    fclose (fp1);
}

event outputInterface(t += TOUTPUT) {
  char names[36];
  sprintf( names, "interfaceMed-%d.dat", pid() );
  FILE * fp2 = fopen (names, "w");
  output_facets (f2,fp2);
  fclose(fp2);
  char command[80];
  sprintf(command, "LC_ALL=C  cat interface* > ALLINTERPhaseFluid-%g.dat",t);
  system(command);// allow to use linux command in the c code to concatenate our files

  sprintf(command, "rm interfaceMed-*");
  system(command);

  sprintf(names, "interfaceMed-%d.dat", pid() );
  output_facets (f3,fp2);
  fclose(fp2);
  sprintf(command, "LC_ALL=C  cat interface* > ALLINTERPhaseMud-%g.dat",t);
  system(command);// allow to use linux command in the c code to concatenate our files

  sprintf(command, "rm interfaceMed-*");
  system(command);
}

// TODO: implement the amplitude calculation later.
// event depthAmplitude (i += 25) {
//   double ampY = 0.0; // essentially fr depth
//   double frLoc = 0.0;
// //   double frFrVal = 0.0;
// //   double frReVal = 0.0;
// //   double frAveVel = 0.0;
// //   int np = 90;
//   // coord c[np];
//   // double depthAveArray[np];
//
//   FILE *fp5 = fopen("amplitude", "a+");
//   FILE *fp3 = fopen("totalDroplets", "a");
//   FILE *fp4 = fopen("dropletDetail", "a");
//
// //   FILE *fp6 = fopen("frFr", "a");
// //   FILE *fp7 = fopen("frRe", "a");
//
//   // first calculate droplet statistics
//   scalar m[];
//   foreach()
//     m[] = f[] > 1e-3;
//   int n = tag (m);
//
//   /**
//   Once each cell is tagged with a unique droplet index, we can easily
//   compute the volume *v* and position *b* of each droplet. Note that
//   we use *foreach (serial)* to avoid doing a parallel traversal when
//   using OpenMP. This is because we don't have reduction operations for
//   the *v* and *b* arrays (yet). */
//
//   double v[n];
// //   coord b[n];
//   for (int j = 0; j < n; j++)
// //     v[j] = b[j].x = b[j].y = b[j].z = 0.;
//     v[j] = 0.;
//   foreach (serial)
//     if (m[] > 0) {
//       int j = m[] - 1;
//       v[j] += dv()*f[];
// //       coord p = {x,y,z};
// //       foreach_dimension()
// // 	b[j].x += dv()*f[]*p.x;
//     }
//
//   foreach(serial)
//   {
//     volDroplet[] = 0.0;
//     if (m[] > 0) {
//       int j = m[] - 1;
//       volDroplet[] = v[j];
// //       coord p = {x,y,z};
// //       foreach_dimension()
// // 	b[j].x += dv()*f[]*p.x;
//     }
//   }
//
//  /**
//  When using MPI we need to perform a global reduction to get the
//  volumes and positions of droplets which span multiple processes. */
//
// #if _MPI
//   MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// //   MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// #endif
//
//   /**
//   Finally we output the volume and position of each droplet to
//   standard output. */
//
//   fprintf (fp3, "%g %d \n", t, n);
//   for (int j = 0; j < n; j++)
//     fprintf (fp4, "%g %d %g \n", t, j, v[j]);
//   fclose (fp3);
//   fclose (fp4);
//
//   position (f, ypos, {0, 1});
//
//   foreach()
//   {
//     ypos[] = volDroplet[]>(259.1*sq(xextent_/pow(2, MAXLEVEL))) ? ypos[] : 0.0;
//   }
//
//   ampY = statsf(ypos).max;
//
//   // amp and amp location
// //   foreach (reduction(max:ampY))
// //   {
// //     if (ampY<ypos[])
// //     {
// //       ampY = ypos[];
// //     }
// //   }
//
//   foreach (reduction(max:frLoc))
//   {
//     if (f[]>1.0E-3 && ampY*0.99<y)
//     {
//       frLoc = x;
//     }
//   }
//
//   fprintf (fp5, "%g %g %g \n", t, frLoc, ampY);
// //   fprintf (ferr, "%g %g \n", t, ampY);
//   fclose (fp5);
// }

// event outputCentVel(t += TOUTPUT) {
//   char resultname[40];
//   sprintf( resultname, "centVel_%g.txt", t );
//   FILE * fp = fopen(resultname, "w");
//   for (double y = 0.; y < xextent_; y += xextent_/pow(2.,LEVEL))
//         fprintf (fp, "%g %g %g %g \n", y, interpolate (u.x, xextent_/2, y), interpolate (u.y, xextent_/2, y), interpolate (f, xextent_/2, y));
//   fclose (fp);
// }

int refRegion(double x,double y, double z){
    int lev;
    if( y < topExtent*0.975 )
      lev = MAXLEVEL;
    else
      lev = MINLEVEL+1;

    return lev;
}

// mesh adaptation
event adapt (i++) {
  double femax = 8.e-3;
  scalar omega[] ;
  vector gf1[], gf2[], gf3[];
  gradients ({f1}, {gf1});
  gradients ({f2}, {gf2});
  gradients ({f3}, {gf3});
//   curvature(f, KAPPA);
  vorticity (u, omega);
  boundary ((scalar *){omega});
//   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, omega}, (double[]){fErr, VelErr, VelErr, KAPPAErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  // adapt_wavelet_limited ((scalar *){ omega, velFluidNorm, velAirNorm}, {f1,f2,f3}, (double[]){ OmegaErr, VelFluidErr, VelAirErr}, refRegion, minlevel = MINLEVEL);
  adapt_wavelet_limited((scalar *){gf1.x, gf1.y, gf2.x, gf2.y, gf3.x, gf3.y , omega, velFluidNorm, velAirNorm}, (double[]){femax, femax, femax, femax, femax, femax, OmegaErr, VelFluidErr, VelAirErr}, refRegion, minlevel = MINLEVEL);
  refine(y<=(7.510*xextent_/pow(2, MAXLEVEL)) && level<MAXLEVEL);
//   adapt_wavelet ((scalar *){f, cs}, (double[]){fErr, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

/**## Numerical Results */
/**
 <span style="color:red"> **The fluid stays stationary ($u_x(t)=u_y(t)=0$), which is apparently different from the analytical solutions and it is not reasonable at all. Why??** </span>*/
// event movies (i += 10) {
//   view (quat = {0.000, 0.000, 0.000, 1.000},
//       fov = 30, near = 0.01, far = 1000,
//       tx = -0.493, ty = -0.012, tz = -1.106,
//       width = 1200, height = 600);
//   box ();
//   cells ();
//   squares (color = "u.x");
//   vectors (u = "u", scale = 0.00005);
//   draw_vof (c = "f", lw = 2.4, fc = {1.0, 0.0, 1.0});
//   save ("movie.mp4");
// }

