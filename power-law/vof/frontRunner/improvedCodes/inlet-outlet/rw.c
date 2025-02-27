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

// #define FILTERED

/**
   Include profiling information. */

// #include "navier-stokes/perfs.h"
#include "./my-perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 13
#define MINLEVEL 2
#define INITLEVEL 11

/** Problem-related parameters */

#define MUDRHO 1003.2 //density ratio, water to air
// #define MURATIO 8.9e-4/17.4e-6 //dynamic viscosity ratio, water to air
#define POWERLAWINDEX 0.390
#define MUMUD 5.73
#define YIELDSTRESS 9.56

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 0.0135272
#define NORMALVEL 0.2754396

#define AIRRHO (MUDRHO/100.10) // 1.12
// #define AIRMU (MUMUD/50.0)
#define AIRMU ((MUMUD*pow((NORMALVEL/NORMALDEPTH), (POWERLAWINDEX-1.0))+YIELDSTRESS/(NORMALVEL/NORMALDEPTH))/60.510)

// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELSLOPE 0.34202 // sinTheta
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))
#define CHANNELTAN (CHANNELSLOPE/CHANNELCOS)
#define GRAV 9.81
#define GRAVRED (GRAV*CHANNELCOS)

#define BPARAM (YIELDSTRESS/(MUDRHO*CHANNELSLOPE*NORMALDEPTH*GRAV))
#define NORMALVELMAX (NORMALVEL*1.0/(1.0-POWERLAWINDEX/(2.0*POWERLAWINDEX+1.0)*(1.0-BPARAM)))

#define FR (NORMALVEL/pow((GRAVRED*NORMALDEPTH), 0.50)) // Froude number for the expected steady-state solution, but this dimensionless parameter is not used in the codes here.

#define WEBER (0.0708*2.99)
#define COEFFST (MUDRHO*pow(NORMALVEL,2.0)*NORMALDEPTH*WEBER)

#define MAXTIME 200.0 // Maximum runtime.
// #define TOUTPUT 0.250
#define DIMLESSTOUTPUT 0.50
#define TOUTPUT (DIMLESSTOUTPUT*NORMALDEPTH/NORMALVEL/CHANNELTAN)

#define DISTAMP 0.21
#define DISTWAVELENGTH (2.0*NORMALDEPTH/CHANNELSLOPE)
#define DISTPERIOD 0.66667

// square domain size
#define xextent_ (277.5*NORMALDEPTH)
#define topExtent (NORMALDEPTH*3.10)
#define initTransitCells 31

// #define KAPPAErr (1e-3)
#define OmegaErr (0.25)
#define fErr (1e-7)
#define VelErr (NORMALVEL/101.0)
#define KErr (1e-3)

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

/**
## Main body of the current codes
*/

double distSurf(double xCoord)
{
  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
//   if (xCoord>=DISTWAVELENGTH/2.0 && xCoord<=DISTWAVELENGTH)
    return (NORMALDEPTH);


  // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
}


double distVel(double xCoord, double yCoord)
{
  // double distDepth = NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_));
  double distDepth = NORMALDEPTH;
  double umax = sqrt(distDepth/NORMALDEPTH)*NORMALVELMAX;
//   if (yCoord<=bParam*distDepth)
//     return (umax*(1-pow((1.0-yCoord/(bParam*distDepth)), (POWERLAWINDEX+1.0)/POWERLAWINDEX)));
//   else if (yCoord<=distDepth)
//     return (umax);
//   else if (yCoord<=distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL))
//     return (umax*(distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL)-yCoord)/(initTransitCells*xextent_/pow(2,MAXLEVEL)));
//   else
//     return 0.0;

  if (yCoord<=(1.0-BPARAM)*NORMALDEPTH)
    return (umax*(1-pow((1.0-yCoord/((1.0-BPARAM)*NORMALDEPTH)), (POWERLAWINDEX+1.0)/POWERLAWINDEX)));
  else if (yCoord<=distDepth)
    return (umax);
  else if (yCoord<=distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL))
//     return (umax*(distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL)-yCoord)/(initTransitCells*xextent_/pow(2,MAXLEVEL)));
    return (umax*1.05);
  else
    return (umax*1.05);
}

double distVelInlet(double xCoord, double yCoord, double depth)
{
  // double distDepth = NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_));
  double distDepth = depth;
  double umax = sqrt(distDepth/NORMALDEPTH)*NORMALVELMAX;
//   if (yCoord<=bParam*distDepth)
//     return (umax*(1-pow((1.0-yCoord/(bParam*distDepth)), (POWERLAWINDEX+1.0)/POWERLAWINDEX)));
//   else if (yCoord<=distDepth)
//     return (umax);
//   else if (yCoord<=distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL))
//     return (umax*(distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL)-yCoord)/(initTransitCells*xextent_/pow(2,MAXLEVEL)));
//   else
//     return 0.0;

  if (yCoord<=(1.0-BPARAM)*NORMALDEPTH)
    return (umax*(1-pow((1.0-yCoord/((1.0-BPARAM)*NORMALDEPTH)), (POWERLAWINDEX+1.0)/POWERLAWINDEX)));
  else if (yCoord<=distDepth)
    return (umax);
  else if (yCoord<=distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL))
//     return (umax*(distDepth+initTransitCells*xextent_/pow(2,MAXLEVEL)-yCoord)/(initTransitCells*xextent_/pow(2,MAXLEVEL)));
    return (umax*1.051);
  else
    return (umax*1.051);
}

double inletDistDepth(double timeVal, double yCoord, double del)
{
  // type-a disturbance
  double eta = NORMALDEPTH*(1.0+DISTAMP*sin(2*pi*timeVal/DISTPERIOD));
	//cout << wwelev << endl;
	if (eta < yCoord - (del / 2.)) {
		return 0.0;
	}
	else if (eta > yCoord + (del / 2.)) {
		return 1.0;
	}
	else {
		// Calculate volume fraction for the given cell with size del and position y
		return clamp((eta - (yCoord - (del / 2.))) / del,0,1);
	}

//   if (eta < yCoord) {
//     return 0.0;
//   }
//   else {
//     return 1.0;
//   }
}

double inletDistVel(double yCoord, double timeVal)
{
  double distDepth = NORMALDEPTH*(1.0+DISTAMP*sin(2*pi*timeVal/DISTPERIOD));
  return (yCoord<=(distDepth) ? distVelInlet(timeVal, yCoord, distDepth) : (yCoord<=(distDepth+initTransitCells*xextent_/pow(2, MAXLEVEL))) ? 1.051*NORMALVELMAX : 0.0);
}

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

   /*
Inlet
*/
// Oystelan's suggestion
// f[left] = dirichlet(inletDistDepth(t, y, Delta));
f[left] = inletDistDepth(t, y, Delta);
// // u.n[left]  = f[]*dirichlet(inletDistVel(y, t)) + (1.-f[])*neumann(0);
u.n[left]  = dirichlet(inletDistVel(y, t));
// u.n[left]  = dirichlet(0.0);
u.t[left] = dirichlet(0.0);

/*
Outlet
*/
// f[right] = neumann(0.);
// uf.n[right] = neumann(0.);
// u.n[right] = neumann(0.);
// uf.t[right] = neumann(0.);
// u.t[right] = neumann(0.);
// p[right] = neumann(0.0);
// pf[right] = neumann(0.0);

f[right] = neumann(0.);
u.n[right] = neumann(0.);
// // u.t[right] = neumann(0.);
// // p[right] = neumann(0);
// // pf[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

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
  tauP = YIELDSTRESS;
  muRef = MUMUD;
  mumax = 77.60;

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
//   periodic (right);

  NITERMAX = 255;
  TOLERANCE = 2.50e-4;
  CFL = 0.470;

  run();
}

//---------------------INITIALIZATION------------------------//
// double distSurf(double xCoord)
// {
//   // return (NORMALDEPTH*(1.0+DISTAMP*sin(2.0*pi*xCoord/xextent_)));
//   return (NORMALDEPTH);
// }

// double normalFlowVel(double yCoord, double depth)
// {
//   double bParam = YIELDSTRESS/(MUDRHO*GRAVRED*NORMALDEPTH*CHANNELSLOPE);
//   double umax = NORMALVEL*(1.0/(1.0-POWERLAWINDEX/(2.0*POWERLAWINDEX+1.0)*(1.0-bParam)));
//   if (yCoord<=bParam*depth)
//     return (umax*(1-pow((1.0-yCoord/(bParam*depth)), (POWERLAWINDEX+1.0)/POWERLAWINDEX)));
//   else
//     return (umax);
// //   return (FR*sqrt(CHANNELCOS*GRAV*NORMALDEPTH)*((1.0+2.0*POWERLAWINDEX)/(1.0+POWERLAWINDEX))*(1.0-pow((1.0-yCoord/NORMALDEPTH), (1.0+POWERLAWINDEX)/POWERLAWINDEX)));
// }

/** ### Init event */
event init (i=0)
{
    if (!restore("restart")){
    scalar omega[];
    // refine(y<=(22.50*NORMALDEPTH) && level < (MAXLEVEL-3));
    // refine(y<=(5.00*NORMALDEPTH) && level < (MAXLEVEL-1));
    // refine(y<=(1.2*NORMALDEPTH) && level < (MAXLEVEL));
    // solid (cs, fs, topExtent - y );

    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      // refine(y<(4.00*NORMALDEPTH) && level < (MAXLEVEL-3));
      // refine(y<(2.00*NORMALDEPTH) && level < (MAXLEVEL-1));
      // refine(y<(1.50*NORMALDEPTH) && level < (MAXLEVEL));
      fraction (f, distSurf(x)-y);
      solid (cs, fs, topExtent - y );

      foreach() {
        // variation of x-component velocity to keep discharge the same
        u.x[] = (y<=topExtent*1.10) ? (y<=NORMALDEPTH ? distVel(x, y) : (y<=(NORMALDEPTH+initTransitCells*xextent_/pow(2,MAXLEVEL)) ? 1.05*NORMALVELMAX : 0.0)) : 0.0;
        // u.x[] = 0.0;
        u.y[] = 0.0;
        // hydrostatic pressure, zero pressure datum at free-surface
        p[] = (y<=distSurf(x)) ? MUDRHO*CHANNELCOS*GRAV*(distSurf(x)-y) : (-1.0)*CHANNELCOS*GRAV*AIRRHO*(y-distSurf(x));

        // boundaryLayer[] = (y<=xextent_/pow(2,MAXLEVEL)*4.0 ? 1.0 : 0.0);
      }
      boundary ((scalar *){u, p});
      vorticity (u, omega);
    }
    while (adapt_wavelet ((scalar *){f, cs, u.x, omega}, (double[]){(fErr/25.0), 1.e-4, VelErr, OmegaErr}, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);

    fractions_cleanup (cs, fs);

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

event snapshot (t += TOUTPUT*5.0) {
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
    FILE *fp1 = fopen(name, "w");
    output_gfs(fp1, translate = true);
    fclose (fp1);
}

event outputInterface(t += TOUTPUT) {
  char names[36];
  sprintf( names, "interfaceMed-%d.dat", pid() );
  FILE * fp2 = fopen (names, "w");
  output_facets (f,fp2);
  fclose(fp2);
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

event dropletsAndhMax (i += 20)
{
  FILE *fp3 = fopen("totalDroplets", "a");
  FILE *fp4 = fopen("dropletDetail", "a");

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
}

event depthAmplitude (i += 20) {
  FILE *fp5 = fopen("amplitude", "a+");
  double ampY = 0.0;
  double frLoc = 0.0;
  foreach (reduction(max:ampY) reduction(max:frLoc))
  {
    double yCoordInt = volDroplet[]>(16.10*sq(xextent_/pow(2, MAXLEVEL))) ? y+0.50*(xextent_/pow(2, MAXLEVEL))-(1-f[])*(xextent_/pow(2, MAXLEVEL)) : 0.0;
    if (ampY<yCoordInt)
    {
      ampY = yCoordInt;
      frLoc = x;
    }
  }
  fprintf (fp5, "%g %g %g \n", t, frLoc, ampY);
//   fprintf (ferr, "%g %g \n", t, ampY);
  fclose (fp5);
}

// event maintainNormalDepth1 (i=1; i<=3; i++) {
//   foreach() {
//     if(x>5.0*xextent_/pow(2,MAXLEVEL))
//     {
//       if (y+xextent_/pow(2, MAXLEVEL)/2.0<=NORMALDEPTH)
//         f[] = 1.0;
//       else if (y+xextent_/pow(2, MAXLEVEL)/2.0>NORMALDEPTH && y-xextent_/pow(2, MAXLEVEL)/2.0<=NORMALDEPTH)
//         f[] = (NORMALDEPTH-(y-xextent_/pow(2, MAXLEVEL)/2.0))/(xextent_/pow(2, MAXLEVEL));
//       else
//         f[] = 0.0;
//     }
//   }
// }
//
// event maintainNormalDepth2 (i=3; i+=80) {
//   double ampY = 0.0;
//   double frLoc = 0.0;
//   foreach (reduction(max:ampY) reduction(max:frLoc))
//   {
//     double yCoordInt = volDroplet[]>(16.10*sq(xextent_/pow(2, MAXLEVEL))) ? y+0.50*(xextent_/pow(2, MAXLEVEL))-(1-f[])*(xextent_/pow(2, MAXLEVEL)) : 0.0;
//     if (ampY<yCoordInt)
//     {
//       ampY = yCoordInt;
//       frLoc = x;
//     }
//   }
//
//   foreach() {
//     if (t>1.55*DISTPERIOD/2.0)
//     {
//     if(x>=1.10*frLoc)
//     {
//       if (y+xextent_/pow(2, MAXLEVEL)/2.0<=NORMALDEPTH)
//         f[] = 1.0;
//       else if (y+xextent_/pow(2, MAXLEVEL)/2.0>NORMALDEPTH && y-xextent_/pow(2, MAXLEVEL)/2.0<=NORMALDEPTH)
//         f[] = (NORMALDEPTH-(y-xextent_/pow(2, MAXLEVEL)/2.0))/(xextent_/pow(2, MAXLEVEL));
//       else
//         f[] = 0.0;
//     }
//     }
//     else
//     {
//       if(x>=xextent_/20.0 && x>=3.0*frLoc)
//       {
//         if (y+xextent_/pow(2, MAXLEVEL)/2.0<=NORMALDEPTH)
//           f[] = 1.0;
//         else if (y+xextent_/pow(2, MAXLEVEL)/2.0>NORMALDEPTH && y-xextent_/pow(2, MAXLEVEL)/2.0<=NORMALDEPTH)
//           f[] = (NORMALDEPTH-(y-xextent_/pow(2, MAXLEVEL)/2.0))/(xextent_/pow(2, MAXLEVEL));
//         else
//           f[] = 0.0;
//       }
//     }
//   }
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
    if( y < topExtent )
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
//   adapt_wavelet_limited ((scalar *){ u.x, u.y, omega }, {f}, (double[]){ VelErr, VelErr, OmegaErr }, refRegion, minlevel = MINLEVEL);
  adapt_wavelet_limited ((scalar *){ omega }, {f}, (double[]){ OmegaErr }, refRegion, minlevel = MINLEVEL);
  refine(y<=(7.5510*xextent_/pow(2, MAXLEVEL)) && level<MAXLEVEL);
//   adapt_wavelet ((scalar *){f, cs}, (double[]){fErr, 0.01}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

/**## Numerical Results */
/**
 <span style="color:red"> **The fluid stays stationary ($u_x(t)=u_y(t)=0$), which is apparently different from the analytical solutions and it is not reasonable at all. Why??** </span>*/
// event movies (i += 10) {
//   view (quat = {0.000, 0.000, 0.000, 1.000},
//       fov = 30, near = 0.01, far = 1000,
//       tx = -0.038, ty = -0.006, tz = -0.096,
//       width = 1620, height = 899);
//   squares (color = "u.x", max = 0.5);
//   draw_vof (c = "f", lw = 1.6);
//   save ("movie.mp4");
// }

