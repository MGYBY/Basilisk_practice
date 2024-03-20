#include "grid/quadtree.h"
#include "./adapt_wavelet_leave_interface_limited.h"
// #include "embed.h"
#include "navier-stokes/centered.h"
#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
// #include "two-phasePL.h"
#include "two-phase.h"
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

#include "my-perfs.h"
// #include "profiling.h"

/** Grid levels */

#define MAXLEVEL 7
#define MINLEVEL 2
#define INITLEVEL 7

// square domain size
#define xextent_ 5.00
#define diam 1.00

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

/** ### Main */
int main()
{
  size (xextent_);

  rho1 = 1.0;
  rho2 = 0.01;

  mu1 = 1.0;
  mu2 = 0.01;

  f.sigma = 0.10;

  init_grid(1 << (INITLEVEL));


  // periodic BC
  periodic (right);

  NITERMAX = 255;
  TOLERANCE = 9.9e-5;
  CFL = 0.4750;

  run();
}

/** ### Init event */
event init (i=0)
{
    if (!restore("dumpSnapshot-6.25")){

    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );

      fraction (f, pow((x-xextent_/2.0), 2.0)+pow((y-xextent_/2.0), 2.0)-pow((diam/2.0), 2.0));

      // solid (cs, fs, topExtent - y );

      foreach() {
        u.x[] = 0.0;

        u.y[] = f[]*(-1.00);

      }
      boundary ((scalar *){u, p});

    }
    while (adapt_wavelet ((scalar *){f, u.x, u.y}, (double[]){(0.001), 0.1/5.0, 0.1/5.0 }, maxlevel = MAXLEVEL, minlevel = MINLEVEL).nf);
  }
}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= 1.0E-5; t += 50.0);

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] += 0.00;
  foreach_face(y)
    // av.y[] -= 1.00*(1.0/(FR*FR))*f[];
    av.y[] -= 1.00*f[];
}

event iterLog(i += 5) {
  // FILE *fp1 = fopen("iterStat.txt", "a+");
  fprintf (ferr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
  // fclose (fp1);
}

/**
### Dump output, Gfs file output, free-surface output*/

event snapshot (i += 1) {
  char nameOut[50];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut, list={f, u.x, u.y, uf.x, uf.y, p});
}

event outputGfsFiles (i += 1) {
    char name[80];
    sprintf(name, "out-%g.gfs", t);
    FILE *fp1 = fopen(name, "w");
    output_gfs(fp1, translate = true, list={f, u.x, u.y, uf.x, uf.y, p});
    fclose (fp1);
}

event outputInterface(i += 1) {
  char names[36];
  sprintf( names, "interfaceMed-%d.dat", pid() );
  FILE * fp2 = fopen (names, "w");
  output_facets (f,fp2);
  fclose(fp2);
  char command[80];
  sprintf(command, "LC_ALL=C  cat interface* > ALLINTER-%g.dat",t);
  system(command);// allow to use linux command in the c code to concatenate our files
}
