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
// #define FILTERED
// #include "two-phase.h"
// #include "log-conform.h"
#include "log-conform-viscoelastic.h"
// #define FILTERED // Smear density and viscosity jumps
#include "two-phaseVE.h"
// #include "vof.h"
// alternatively, use momentum-conserving scheme
#include "navier-stokes/conserving.h"
// #include "./myTension.h"
#include "tension.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"

/**
   Include profiling information. */

// #include "./my-perfs.h"
#include "navier-stokes/perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 10
#define MINLEVEL 3
#define INITLEVEL 10

/** Problem-related parameters */

// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELTAN 0.08
#define CHANNELSLOPE (pow((CHANNELTAN*CHANNELTAN/(1.0+CHANNELTAN*CHANNELTAN)), 0.50)) // sinTheta
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))
#define GRAV 1.0
#define GRAVRED (GRAV*CHANNELCOS)

#define FR 0.650
// #define RE 32.00
#define RE (3.0*FR*FR/CHANNELTAN)
#define BETAPARAM 0.50
#define WI 1.0
// more comfortable variable names
#define MUTOT (1.0/RE)
#define MUSOL (BETAPARAM*MUTOT)
#define MUPOLY ((1.0-BETAPARAM)*MUTOT)
#define L1 (WI)
#define L2 (WI*BETAPARAM)
#define GF (1.0/RE*(1.0-BETAPARAM)/WI)

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 1.00
#define NORMALVEL 1.00

#define MUDRHO 1.00 //density ratio, water to air
#define rhoR 100.0
#define muR 60.0
#define AIRRHO (MUDRHO/rhoR) // 1.12
// #define AIRMU (MUMUD/48.0)
#define AIRMU ((1.0/RE)/muR)

#define WEBER 250.0
#define COEFFST (MUDRHO*pow(NORMALVEL,2.0)*NORMALDEPTH/WEBER)

#define MAXTIME 50.0 // Maximum runtime.
#define TOUTPUT 1.0
#define TOUTPUTDUMP 5.0

#define DISTAMP 0.10

// square domain size
// #define xextent_ (NORMALDEPTH*50.0)
#define DIMLESSLENGTH 3.20
#define xextent_ (DIMLESSLENGTH/CHANNELTAN)
#define topExtent (NORMALDEPTH*4.0)

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
// #define mu(muTemp, mu2, f)  (1./(clamp(f,0,1)/muTemp  + (1.-clamp(f,0,1))/mu2))

// #define KAPPAErr (1e-3)
#define OmegaErr (NORMALVEL/NORMALDEPTH/4.0)
#define fErr (1.0e-7)
#define VelFluidErr (NORMALVEL/40.50)
#define VelErr VelFluidErr
#define VelAirErr (NORMALVEL/12.50)
#define KErr (1.0e-3) // important for large Weber num
#define normTauDevFieldError (2.50e-3)

/**
## Main body of the current codes
*/

/**
  slip at the top
*/
u.t[top] = neumann(0.);
u.n[top] = dirichlet(0.);
// /**
//  no slip at the bottom
// */
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
// tau_qq[left] = dirichlet(0);
//
u.n[embed] = dirichlet(0.);
u.t[embed] = neumann(0.);

int j = 0;

// scalar lambdav[], mupv[];
scalar velFluidNorm[], velAirNorm[];
scalar omega[];
scalar volDroplet[];
scalar ypos[], xpos[];
scalar norm_tau_dev_field[]; // Storing the norm of tau_p_dev

/** ### Main */
// int main()
int main(int argc, char const *argv[])
{
  size (xextent_);

  rho1 = MUDRHO;
  rho2 = AIRRHO;

  mu1 = MUSOL; // the solvent visco
  mu2 = AIRMU;

  // VE3D parameters here
  G1 = GF;
  G2 = 0.0;
  lambda1 = L1;
  lambda2 = 0.0;

  f.sigma = COEFFST;

  init_grid(1 << (INITLEVEL));

  // periodic BC
  periodic (right);

  NITERMAX = 25;
  TOLERANCE = 3.0e-4;
  CFL = 0.465;
  // DT = 0.005;

  /**
  We set a maximum timestep. This is necessary for proper temporal
  resolution of the viscoelastic stresses. */
//   DT = 2.0e-4;

  run();

  return 0;
}

//---------------------INITIALIZATION------------------------//
double distSurf(double xCoord)
{
  return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}

/** ### Init event */
event init (i=0)
{
//   scalar strComp = tau_p.x.x;
// //   strComp[bottom] = dirichlet (0.);
//   strComp[bottom] = neumann (0.);
  scalar strComp = tau_p.y.y;
//   strComp[bottom] = dirichlet (0.);
  strComp[bottom] = dirichlet (0.);
//   strComp = tau_p.x.y;
//   // strComp[bottom] = dirichlet(0.);
//   strComp[bottom] = neumann (0.);
//   strComp = tau_p.y.x;
//   strComp[bottom] = neumann (0.);

//     if (!restore("restart")){
    // double femax = 1e-4;
    // double uemax = 2e-2;
//     solid (cs, fs, topExtent - y );
    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      fraction (f, distSurf(x)-y);
      solid (cs, fs, topExtent - y );

      foreach() {
        u.x[] = y<=distSurf(x) ? 3.0/2.0*distSurf(x)*(2.0*(y/distSurf(x))-(y/distSurf(x))*(y/distSurf(x))) : 1.02*3.0/2.0*(2.0*distSurf(x)-distSurf(x)*distSurf(x));
        // u.x[] = 0.0;
        u.y[] = 0.0;
        // hydrostatic pressure, zero pressure datum at free-surface
        p[] = (y<=distSurf(x)) ? (1.0/(FR*FR))*(distSurf(x)-y) : (-1.0)*AIRRHO*(1.0/(FR*FR))*(y-distSurf(x));
        tau_p.x.x[] = (y<=distSurf(x)) ? 9.0/2.0*(WI-BETAPARAM*WI)*pow((2.0-2.0*y),2.0) : 0.0;
        tau_p.y.y[] = 0.0;
        tau_p.x.y[] = (y<=distSurf(x)) ? 3.0/2.0*(2.0-2.0*y) : 0.0;
        tau_p.y.x[] = (y<=distSurf(x)) ? 3.0/2.0*(2.0-2.0*y) : 0.0;

        velFluidNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*f[];
        velAirNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*(1.-f[]);
      }
      boundary ((scalar *){u});
    }
    while (adapt_wavelet ((scalar *){f, u.x, u.y}, (double[]){fErr, VelErr, VelErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);
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
    velFluidNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*f[];
    velAirNorm[] = pow(u.x[]*u.x[]+u.y[]*u.y[], 0.50)*(1.0-f[]);
  }
}

// from https://basilisk.fr/sandbox/hugofranca/saramito_droplet_spreading/droplet_spreading.c
event propertiesStress (i++) {
  double one_third = 1.0/3.0;
  foreach() {
    /// === Norm of the polymeric stress tensor
    double trace_tau = tau_p.x.x[] + tau_p.y.y[];
    double tau_dev_xx = tau_p.x.x[] - one_third*trace_tau;
    double tau_dev_yy = tau_p.y.y[] - one_third*trace_tau;
    double tau_dev_qq = -1.0*one_third*trace_tau;
    double tau_dev_xy = tau_p.x.y[];
    double norm_tau_dev = sqrt( sq(tau_dev_xx) + sq(tau_dev_yy) + sq(tau_dev_qq) + 2.0*sq(tau_dev_xy) );
    norm_tau_dev_field[] = norm_tau_dev;
  }

  boundary(all);
}

event iterLog(i += 5) {
  // FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
  // fclose (fp1);
}

/**
### Dump output, Gfs file output, free-surface output*/

event snapshot (t += TOUTPUTDUMP*5.0) {
  char nameOut[50];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut);
}

event outputGfsFiles (t += TOUTPUT) {
  char name[80];
  sprintf(name, "out-%g.gfs", t);
  FILE *fp1 = fopen(name, "w");
  output_gfs(fp1, translate = true, list={f, u.x, u.y, uf.x, uf.y, p, omega, tau_p.x.x, tau_p.y.y, tau_p.x.y});
//   output_gfs(fp1, translate = true, list={f, u.x, u.y, uf.x, uf.y, p});
  fclose (fp1);
}

// event outputField (t += TOUTPUT) {
//     char nameField[40], nameSlice[40];
//     sprintf(nameField, "field-%g.txt", t);
// //     sprintf(nameSlice, "slice-%g.txt", t);
//     FILE *fp1 = fopen(nameField, "w");
// //     FILE *fp2 = fopen(nameSlice, "w");
//     output_field({f, u.x, tau_p.x.x, tau_p.y.y, tau_p.x.y}, fp1, linear=true);
// //     for (j=1; j<30; j++)
// //       fprintf(fp2, "%g %g %g %g\n", (xextent_/2.0), (j*1.1*NORMALDEPTH/30.0), interpolate(u.x, (xextent_/2.0), (j*1.1*NORMALDEPTH/30.0), 0.0), interpolate(f, (xextent_/2.0), (j*1.1*NORMALDEPTH/30.0), 0.0));
//     fclose (fp1);
// //     fclose (fp2);
// }

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

event depthAmplitude (i += 25) {
  double ampY = 0.0; // essentially fr depth
  double frLoc = 0.0;
//   double frFrVal = 0.0;
//   double frReVal = 0.0;
//   double frAveVel = 0.0;
//   int np = 90;
  // coord c[np];
  // double depthAveArray[np];

  FILE *fp5 = fopen("amplitude", "a+");
  FILE *fp3 = fopen("totalDroplets", "a");
  FILE *fp4 = fopen("dropletDetail", "a");

//   FILE *fp6 = fopen("frFr", "a");
//   FILE *fp7 = fopen("frRe", "a");

  // first calculate droplet statistics
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

  fprintf (fp3, "%g %d \n", t, n);
  for (int j = 0; j < n; j++)
    fprintf (fp4, "%g %d %g \n", t, j, v[j]);
  fclose (fp3);
  fclose (fp4);

  position (f, ypos, {0, 1});
  position (f, xpos, {1, 0});

  foreach()
  {
    ypos[] = volDroplet[]>(81.1*sq(xextent_/pow(2, MAXLEVEL))) ? ypos[] : 0.0;
  }

  ampY = statsf(ypos).max;

  // amp and amp location
//   foreach (reduction(max:ampY))
//   {
//     if (ampY<ypos[])
//     {
//       ampY = ypos[];
//     }
//   }

  foreach (reduction(max:frLoc))
  {
//     if (0.99*(ampY-0.5*xextent_/pow(2, MAXLEVEL))<y)
//     {
//       frLoc = x;
//     }

    if (ypos[]==ampY)
    {
      frLoc = xpos[];
    }
  }

  fprintf (fp5, "%g %g %g \n", t, frLoc, ampY);
//   fprintf (ferr, "%g %g \n", t, ampY);
  fclose (fp5);

  // depth-average frontal properties
//   for (int l = 1; l <= np; l++)
//   {
//     frAveVel += interpolate(f, frLoc, (ampY*1.10)/np*l, 0.0)*interpolate(u.x, frLoc, (ampY*1.10)/np*l, 0.0)/np;
//   }
//
//   frFrVal = frAveVel/(pow(GRAVRED*ampY, 0.50)+1.0E-25);
//   frReVal = 1.0/((YIELDSTRESS/(2.0*MUDRHO*pow(frAveVel, 2.0)+1.0E-26))+(MUMUD/(MUDRHO*pow(ampY, POWERLAWINDEX)*pow(frAveVel, (2.0-POWERLAWINDEX))+1.0E-26)));
//
//   fprintf (fp6, "%g %g \n", t, frFrVal);
//   fprintf (fp7, "%g %g \n", t, frReVal);
//
//   fclose (fp6);
//   fclose (fp7);
}

int refRegion(double x,double y, double z){
  int lev;
  if( y < topExtent*0.98 )
    lev = MAXLEVEL;
  else
    lev = MINLEVEL+1;

  return lev;
}

// mesh adaptation
event adapt (i++) {
  scalar omega[] ;
  //   curvature(f, KAPPA);
  vorticity (u, omega);
  boundary ((scalar *){omega});
  //   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, omega}, (double[]){fErr, VelErr, VelErr, KAPPAErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  adapt_wavelet_limited ((scalar *){ omega, velFluidNorm, velAirNorm, norm_tau_dev_field}, {f}, (double[]){ OmegaErr, VelFluidErr, VelAirErr, normTauDevFieldError}, refRegion, minlevel = MINLEVEL);
  // adapt_wavelet ((scalar *){  velFluidNorm, velAirNorm, f}, (double[]){ VelFluidErr, VelAirErr, 0.0001}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  refine(y<=(4.50*xextent_/pow(2, MAXLEVEL)) && level<MAXLEVEL);
  //   adapt_wavelet ((scalar *){f, cs}, (double[]){fErr, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}
