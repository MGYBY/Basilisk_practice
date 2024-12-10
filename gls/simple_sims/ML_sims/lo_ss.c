#include <sys/stat.h>
// #include "grid/cartesian.h"
#include "grid/multigrid.h"
#include "spherical.h"

// ML headers
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

#if corforce
const double Omega = 7.292205e-5;
#define F0() (2.*Omega*sin(y*pi/180.))

#include "layered/coriolis.h"
# endif

// SV headers
// #include "saint-venant.h"
// #include "my-saint-venant.h" // use HLLC solver

// SGN headers
// # include "green-naghdi.h"

#include "terrain.h"

#define GRAVITY 9.81
#define LENGTH 3000
#define ENDTIME 500.

// terrain-related parameters
#define lakeLevel 0.00
#define lonBL -80.05
#define latBL 43.0
#define aspectRaTio 3
#define LY (44.4-43.0)
#define LX (LY*aspectRaTio)

// wind
#define windAmp 15.0
#define windSin ((1.0)*1.0/pow(2.0,0.50))
#define windCos ((1.0)*pow((1.0-windSin*windSin),0.50))
#define rho_air 1.2
#define rho0 1000.0

// ML setups
#define NL 5

// simulation time
#define day (24.0*3600.0)
#define tFinal (day*2.50)

// for output
#define tOutput (0.125*3600.0)

scalar etamax[];
scalar etamin[];
scalar etamax_time[];
scalar zbs[];

/**
In the experiment, three radials of wave gauges were set up southwards from
ground zero towards the shore. These are defined here as R1, R2 and R3.
The highest number is closest to ground zero, and descends with distance. */

int main()
{
//   struct stat st = {0};
//   if (stat("./R1", &st) == -1) { mkdir("./R1", 0755); }
//   if (stat("./R2", &st) == -1) { mkdir("./R2", 0755); }
//   if (stat("./R3", &st) == -1) { mkdir("./R3", 0755); }

  nl = NL;
  dimensions (nx = 3);
  G = GRAVITY;

  Radius = 6371220.;
  breaking = 0.07;
  CFL_H = 0.48;
//   CFL = 0.45; // CFL number should be sufficiently small
  theta_H = 0.55;

  size (LX);
  origin (lonBL,latBL);
  N = 800;
  
  run();
}


// Gauge gauges[] = {
//   {"R1/17.dat",	66.97,		-249.95},
//   {"R1/16.dat",	118.80,		-443.37},
//   {"R1/15.dat",	168.50,		-628.84},
//   {NULL}
// };


/**
The [DEM terrain]
*/

// scalar hstart1[];
// scalar hstart2[];

void laplacian_smoothing()
{
  for (int i = 0; i < 2; i++) {
    foreach() {
      if (zb[] < 0.)
	zbs[] = (zb[1] + zb[-1] + zb[0,1] + zb[0,-1] +
		 zb[1,1] + zb[-1,-1] + zb[-1,1] + zb[1,-1])/8.;
      else
	zbs[] = zb[];
    }
    foreach()
      zb[] = zbs[];
  }
}

event init (i = 0)
{
  // for debugging
  char namesZb[36];
  sprintf( namesZb, "zbInit.asc" );
  FILE * fp1 = fopen (namesZb, "w");

  terrain (zb, "../terrain/dem_2", NULL);
//   laplacian_smoothing(); // the deepest point is 244m
  conserve_elevation();

  foreach() {
//     zb[] = zb[]-lakeLevel;

    // initialization for ML method 1.
    foreach_layer()
    {
//       h[] = (max(0., (-1.0)*zb[])) / nl;
//       h[] = zb[]>=0 ? 0.0 : ((-1.0)*zb[]);
      h[] = zb[]>=0 ? 0.0 : fabs((-1.0)*zb[]/nl);
//       h[] = x>(lonBL+LX) ? 0.0 : h[];
    }

    etamax_time[] = 0.;
  }

  output_grd(zb, fp1);
  fclose(fp1);

  u.n[left]  = - radiation(0);
  u.n[right] = + radiation(0);
  u.n[top]   = + radiation(0);
  u.n[bottom]   = - radiation(0);
}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= tFinal; t += 1.e4);

event iterLog(i += 5) {
  // FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (ferr, "%d %g %g \n", i, t, dt);
  // fclose (fp1);
}

/**
From [examples/tohoku.c](examples/tohoku.c) */

event friction (i++)
{
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
    {
      u.x[] /= a;
    }

    h[] = h[]<dry ? 0.0 : h[];
    foreach_dimension()
    {
      u.x[] = h[]<dry ? 0.0 : u.x[];
    }

    if (h[] > dry && h[] + zb[] > etamax[])
      etamax[] = h[] + zb[];
  }
  boundary ({etamax, u});
}

vector windVec[];
scalar windLocalAmp[];
event acceleration (i++)
{
  /**
  The wind stress is added directly as an acceleration, only in the
  topmost layer and only if the fluid layer thickness is larger than
  10 metres. We also interpolate linearly in time, between the times
  associated with `wind1` and `wind2`. */

  foreach()
  {
    windVec.x[] = windAmp*windSin;
    windVec.y[] = windAmp*windCos;
    windLocalAmp[] = pow(windVec.x[]*windVec.x[]+windVec.y[]*windVec.y[], 0.50);
  }

  foreach_face() {
    point.l = nl - 1;
    // wind only act on top layer greater than 0.1m
//     if (hf.x[] > 0.10) {
    if (hf.x[] > (0.10)) {
      double tauw = (windVec.x[] + windVec.x[-1])/2.0;

      double n = (0.80+0.065*windLocalAmp[])*windLocalAmp[];

      ha.x[] += n*tauw/rho0;
    }
  }
}

// face vector haW;
// event windStress (i++)
// {
//   haW = new face vector[nl];
//
//   foreach()
//   {
//     windVec.x[] = windAmp*windSin;
//     windVec.y[] = windAmp*windCos;
//     windLocalAmp[] = pow(windVec.x[]*windVec.x[]+windVec.y[]*windVec.y[], 0.50);
//   }
//
//   foreach_face() {
//     point.l = nl - 1;
//     // wind only act on top layer greater than 0.1m
// //     if (hf.x[] > 0.10) {
//     if (hf.x[] < (1.0*0.10)) {
//       double tauw = (windVec.x[] + windVec.x[-1])/2.0;
//
//       double n = (0.80+0.065*windLocalAmp[])*windLocalAmp[];
//
//       haW.x[] = n*tauw/rho0;
//     }
//   }
//
//   foreach_face()
//     foreach_layer()
//       hu.x[] += dt*haW.x[];
//
//   foreach()
//     foreach_layer() {
//       foreach_dimension()
// 	u.x[] += dt*(haW.x[] + haW.x[1])/(hf.x[] + hf.x[1] + dry);
// #if dimension == 2
//       // metric terms
//       double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
//       double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
//       double ux = u.x[], uy = u.y[];
//       double fG = uy*dmdl - ux*dmdt;
//       u.x[] += dt*fG*uy;
//       u.y[] -= dt*fG*ux;
// #endif // dimension == 2
//     }
//   delete ((scalar *){haW});
// }

// event windStress (i++)
// {
//   foreach() {
//     // formula by Wu (1982) JGRC
//     double Cw = h[] > dry ? (0.80+0.065*windAmp)*(1.0e-3) : 0.0;
//     u.x[] += Cw*windAmp*(windAmp*windSin);
//     u.y[] += Cw*windAmp*(windAmp*windCos);
//   }
//   boundary ({h,u});
// }

double nu_H = 10; // m^2/s

event viscous_term (i++)
{
  if (nu_H > 0.) {
    vector d2u[];
    foreach_layer() {
      double dry = 1.;
      foreach()
	foreach_dimension()
	d2u.x[] = 2.*(sq(fm.x[1])/(cm[1] + cm[])*u.x[1]*(h[1] > dry) +
		      sq(fm.x[])/(cm[-1] + cm[])*u.x[-1]*(h[-1] > dry) +
		      sq(fm.y[0,1])/(cm[0,1] + cm[])*u.x[0,1]*(h[0,1] > dry) +
		      sq(fm.y[0,-1])/(cm[0,-1] + cm[])*u.x[0,-1]*(h[0,-1] > dry))
	/(sq(Delta)*cm[]);
      foreach()
	foreach_dimension() {
	double n = 2.*(sq(fm.x[1])/(cm[1] + cm[])*(1. + (h[1] <= dry)) +
		       sq(fm.x[])/(cm[-1] + cm[])*(1. + (h[-1] <= dry)) +
		       sq(fm.y[0,1])/(cm[0,1] + cm[])*(1. + (h[0,1] <= dry)) +
		       sq(fm.y[0,-1])/(cm[0,-1] + cm[])*(1. + (h[0,-1] <= dry)))
	  /(sq(Delta)*cm[]);
	u.x[] = (u.x[] + dt*nu_H*d2u.x[])/(1. + dt*nu_H*n);
      }
    }
  }
}

/**
For a simple way to supplement the numerical wave gauge array, fields
to store the maximum (and minimum) wave height along with the time it
was reached are implemented.  */

event maxmin (i++)
{
  scalar etainst[];
  foreach() {

    double H = 0.;
    foreach_layer() {
      H += h[];
    }
    etainst[] = H < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);

    if (etainst[] > etamax[]) {
      etamax[] = etainst[];
      etamax_time[] = t;
    }

    if (etainst[] < etamin[])
      etamin[] = etainst[];

  }
  boundary ({etamax, etamin, etamax_time});
}

event outputGrdField (t+=tOutput)
{
  char names[45];
  scalar totalDepth[], totalDepth2[], uSurf[];
//   scalar uAmp[];
  // FIXME: more general coding for parallelism
//   foreach(){
//     uAmp[] = 0.0;
//     foreach_layer {
//     uAmp[] += norm(u);
//     }
//     uAmp[] /= nl;
//   }

  foreach() {

    totalDepth[] = 0.;
    totalDepth2[] = eta[];
    foreach_layer() {
      totalDepth[] += h[];
    }

    totalDepth2[] -= zb[];

    point.l = 0;
    uSurf[] = pow(sq(u.x[]) + sq(u.y[]), 0.5);
  }

  sprintf( names, "h-Grd-%g.asc", t );
  FILE * fp1 = fopen (names, "w");
//   output_grd(h, fp1);
  output_grd(eta, fp1);
  fclose(fp1);

  sprintf( names, "totalDepth-Grd-%g.asc", t );
  FILE * fp3 = fopen (names, "w");
//   output_grd(h, fp1);
  output_grd(totalDepth, fp3);
  fclose(fp3);

  sprintf( names, "totalDepth2-Grd-%g.asc", t );
  FILE * fp4 = fopen (names, "w");
//   output_grd(h, fp1);
  output_grd(totalDepth2, fp4);
  fclose(fp4);

  sprintf( names, "uAmp-Grd-%g.asc", t );
  FILE * fp2 = fopen (names, "w");
  output_grd(uSurf, fp1);
  fclose(fp2);
}


// event logGauges (i++) {
//   if (i == 0) {
//     fprintf (ferr,
// 	"t dt mgp.i mgp.nrelax grid->tn perf.t perf.speed npe\n");
//   }
//   fprintf (ferr, "%g %g %d %d %ld %g %g %d\n",
//       t, dt, mgp.i, mgp.nrelax, grid->tn*nl, perf.t, perf.speed*nl, npe());
//
//   output_gauges (gauges, {eta});
// }


/**
While unnecessary, animations can be output. Useful primarily for sanity checks. */

// event movies (t += 0.5) {
//
//   scalar m[], etam[];
//   foreach()
//   {
//     double H = 0.;
//     foreach_layer() {
//       H += h[];
//     }
//     m[] = -zb[];
//     etam[] = H < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);
//   };
//   boundary ({m, etam});
//
//   output_ppm (eta, mask = m, min = -0.5, max = 0.5, n = 1024, linear = true, box = {{-1250.,-1930.},{1000.,500.}}, file = "eta.mp4");
//   output_ppm (etamax, mask = m, min = 0.0, max = 0.5, n = 1024, linear = true, box = {{-1250.,-1930.},{1000.,500.}}, file = "etamax.mp4");
// }


/**
At the end of the simulation, the maximum and minimum fields are output.  */

// event end (t = ENDTIME) {
//   char name4[60];
//   sprintf(name4, "etamax.dat");
//   FILE * fp4 = fopen (name4, "w");
//   output_field ({etamax}, fp4, n = N );
//   fclose(fp4);
//
//   char name5[60];
//   sprintf(name5, "etamin.dat");
//   FILE * fp5 = fopen (name5, "w");
//   output_field ({etamin}, fp5, n = N );
//   fclose(fp5);
//
//   char name6[60];
//   sprintf(name6, "etamax_time.dat");
//   FILE * fp6 = fopen (name6, "w");
//   output_field ({etamax_time}, fp6, n = N );
//   fclose(fp6);
// }

