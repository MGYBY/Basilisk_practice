#include <sys/stat.h>
#include "spherical.h"

// ML headers
// #include "grid/multigrid.h"
// #include "layered/hydro.h"
// #include "layered/nh.h"
// #include "layered/remap.h"
// #include "layered/perfs.h"

// SV headers
#include "saint-venant.h"

// SGN headers
// # include "green-naghdi.h"

#include "terrain.h"

#define GRAVITY 9.81
#define LENGTH 3000
#define ENDTIME 500.

// terrain-related parameters
#define lakeLevel 74.77
#define lonBL -79.895
#define latBL 43.095
#define aspectRaTio 4
#define LY (44.45-43.095)
#define LX (LY*aspectRaTio)

scalar etamax[];
scalar etamin[];
scalar etamax_time[];

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

//   nl = 3;
  G = GRAVITY;
//   breaking = 0.07;
//   CFL_H = 0.5;
  CFL = 0.46; // CFL number should be sufficiently small
  // FIXME: define variables and domain size related to spherical coordinate
  size (LENGTH);
  origin (lonBL,latBL);
  N = 128;
  
  run();
}


// Gauge gauges[] = {
//   {"R1/17.dat",	66.97,		-249.95},
//   {"R1/16.dat",	118.80,		-443.37},
//   {"R1/15.dat",	168.50,		-628.84},
//   {"R1/14.dat",	214.25,		-799.59},
//   {"R1/13.dat",	248.41,		-927.07},
//   {"R1/12.dat",	282.25,		-1053.36},
//   {"R1/11.dat",	300.39,		-1121.07},
//   {"R1/10.dat",	321.37,		-1199.38},
//   {"R1/9.dat",	330.84,		-1234.71},
//   {"R1/8.dat",	337.54,		-1259.74},
//   {"R1/7.dat",	344.33,		-1285.06},
//   {"R1/6.dat",	353.32,		-1318.62},
//   {"R1/5.dat",	362.55,		-1353.06},
//   {"R1/4.dat",	367.91,		-1373.08},
//   {"R1/3.dat",	371.86,		-1387.80},
//   {"R1/2.dat",	375.65,		-1401.93},
//   {"R1/1.dat",	377.07,		-1407.23},
//   {"R2/12.dat",	0.,		-934.17},
//   {"R2/11.dat",	0.,		-1061.87},
//   {"R2/10.dat",	0.,		-1197.20},
//   {"R2/9.dat",	0.,		-1284.97},
//   {"R2/8.dat",	0.,		-1345.63},
//   {"R2/7.dat",	0.,		-1406.89},
//   {"R2/6.dat",	0.,		-1449.86},
//   {"R2/5.dat",	0.,		-1491.92},
//   {"R2/4.dat",	0.,		-1514.78},
//   {"R2/3.dat",	0.,		-1538.56},
//   {"R2/2.dat",	0.,		-1557.45},
//   {"R2/1.dat",	0.,		-1566.60},
//   {"R3/14.dat",	-67.98,		-349.75},
//   {"R3/13.dat",	-106.19,	-546.31},
//   {"R3/12.dat",	-141.78,	-729.41},
//   {"R3/11.dat",	-164.58,	-846.69},
//   {"R3/10.dat",	-186.10,	-957.39},
//   {"R3/9.dat",	-248.03,	-1276.03},
//   {"R3/8.dat",	-264.72,	-1361.89},
//   {"R3/7.dat",	-286.42,	-1473.49},
//   {"R3/6.dat",	-297.93,	-1532.73},
//   {"R3/5.dat",	-311.25,	-1601.24},
//   {"R3/4.dat",	-317.07,	-1631.16},
//   {"R3/3.dat",	-323.93,	-1666.46},
//   {"R3/2.dat",	-326.37,	-1679.03},
//   {"R3/1.dat",	-327.48,	-1684.71},
//   {NULL}
// };


/**
The [DEM terrain](#raumann2002) is preprocessed to be up-positive and
have "zero" elevation set at lake level. The initial disturbance is
a cavity and is of the function:

$$ \eta_{0} = - D    \left (  \frac{1}{3}   \left (   \frac{r}{W}   \right ) ^{4}     -     \frac{4}{3}   \left (   \frac{r}{W}   \right ) ^{2}   + 1  \right ) $$

...where $\eta_{0}$ is the initial water level, $r$ is the horizontal
distance from ground zero (origin), $D$ and $W$ are parameters set by
empirical relations involving charge yield and depth, and physically
define the depth and radius of the cavity respectively.

The function is valid for $r < W\sqrt{3}$ and is a quadratic deformation
where, in 3D, the volume 'displaced' from  $r < W\sqrt{2}$ is equal to the
volume added as a 'lip' from $W\sqrt{2} < r < W\sqrt{3}$. At higher $r$
the water level remains at zero.

For efficiency, the computational domain is limited to the appropriate
area of the lake, and radiation boundaries are added to mitigate
reflections.  */

// scalar hstart1[];
// scalar hstart2[];

event init (i = 0)
{
  terrain (zb, "../terrain/dem", NULL);
  conserve_elevation();

  foreach() {
    zb[] = zb[]-lakeLevel;
    
    hstart1[] = max(0., -zb[]);
    if (sqrt(x*x + y*y) < (WIDTH * sqrt(3.0)) )
      hstart2[] = (   -DEPTH * ( (1./3.)*pow((pow((pow(x,2) + pow(y,2)),0.5)/WIDTH),4.) - (4./3.)*pow((pow((pow(x,2) + pow(y,2)),0.5)/WIDTH),2.) + 1. )   );
    
    foreach_layer()
      h[] = max(0., -zb[]);

//     etamax_time[] = 0.;
    // TODO: define adapt-related variables

  }

  u.n[left]  = - radiation(0);
  u.n[right] = + radiation(0);
  u.n[top]   = + radiation(0);
  u.n[bottom]   = - radiation(0);
}


/**
From [examples/tohoku.c](examples/tohoku.c) */

event friction (i++)
{
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
    if (h[] > dry && h[] + zb[] > etamax[])
      etamax[] = h[] + zb[];
  }
  boundary ({etamax, u});
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


event logGauges (i++) {
  if (i == 0) {
    fprintf (ferr,
	"t dt mgp.i mgp.nrelax grid->tn perf.t perf.speed npe\n");
  }
  fprintf (ferr, "%g %g %d %d %ld %g %g %d\n", 
      t, dt, mgp.i, mgp.nrelax, grid->tn*nl, perf.t, perf.speed*nl, npe());
  
  output_gauges (gauges, {eta});
}


/**
While unnecessary, animations can be output. Useful primarily for sanity checks. */

event movies (t += 0.5) {

  scalar m[], etam[];
  foreach()
  {
    double H = 0.;
    foreach_layer() {
      H += h[];
    }
    m[] = -zb[];
    etam[] = H < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);
  };
  boundary ({m, etam});

  output_ppm (eta, mask = m, min = -0.5, max = 0.5, n = 1024, linear = true, box = {{-1250.,-1930.},{1000.,500.}}, file = "eta.mp4");
  output_ppm (etamax, mask = m, min = 0.0, max = 0.5, n = 1024, linear = true, box = {{-1250.,-1930.},{1000.,500.}}, file = "etamax.mp4");
}


/**
At the end of the simulation, the maximum and minimum fields are output.  */

event end (t = ENDTIME) {
  char name4[60];
  sprintf(name4, "etamax.dat");
  FILE * fp4 = fopen (name4, "w");
  output_field ({etamax}, fp4, n = N );
  fclose(fp4);

  char name5[60];
  sprintf(name5, "etamin.dat");
  FILE * fp5 = fopen (name5, "w");
  output_field ({etamin}, fp5, n = N );
  fclose(fp5);

  char name6[60];
  sprintf(name6, "etamax_time.dat");
  FILE * fp6 = fopen (name6, "w");
  output_field ({etamax_time}, fp6, n = N );
  fclose(fp6);
}

