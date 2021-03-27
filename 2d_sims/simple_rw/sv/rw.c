/**
# Periodic wave propagation over an ellipsoidal shoal

We follow [Lannes and Marche, 2014](/src/references.bib#lannes2014)
and try to reproduce the experiment of [Berkhoff et al,
1982](/src/references.bib#berkhoff1982). The numerical wave tank is
25^2^ metres and periodic waves are generated on the left-hand side
and are damped on the right-hand side. */

// #include "grid/multigrid.h"
// #include "grid/cartesian.h"
// #include "green-naghdi.h"
#include "saint-venant.h"

#define MAXLEVEL 10
#define MINLEVEL 0

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.05;
double disPeriod = 0.933;
double simTime = 40.0;
double Lx = 40.0;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

int main()
{
  size(Lx);
  G = 9.81;
  init_grid(1 << MAXLEVEL);
  CFL = 0.40; // CFL number should be sufficiently small
  // N = 1024;

  /**
  We turn off limiting to try to reduce wave dissipation. */
  // gradient = NULL;

  run();
}

event init(i = 0)
{
  // for 2-D case 40x3m
  mask(y > 3. ? top : none);
  // L0 = 40.;

  h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
  u.n[left] = dirichlet(normalVelocity);

  u.n[right] = neumann(0.);
  h[right] = neumann(0.);

  u.n[top] = neumann(0.);
  h[top] = neumann(0.);

  u.n[bottom] = neumann(0.);
  h[bottom] = neumann(0.);

  /**
  The bathymetry is an inclined and skewed plane combined with an
  ellipsoidal shoal. 
  
  ~~~gnuplot Bathymetry
  set term @PNG enhanced size 640,640 font ",8"
  set output 'bathy.png'
  set pm3d map
  set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
                        0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392,	  \
                        0.625 1 0.9333 0, 0.75 1 0.4392 0,		  \
                        0.875 0.9333 0 0, 1 0.498 0 0 )
  set size ratio -1
  set xlabel 'x (m)'
  set ylabel 'y (m)'
  splot [-10:12][-10:10]'end' u 1:2:4 w pm3d t ''
  # we remove the large border left by gnuplot using ImageMagick
  ! mogrify -trim +repage bathy.png
  ~~~
  */
  foreach ()
  {
    zb[] = -So * x;
    h[] = normalDepth;
    u.x[] = normalVelocity;
    u.y[] = 0.;
  }
}

/**
To implement an absorbing boundary condition, we add an area for $x >
12$ for which quadratic friction increases linearly with $x$. */
// only use Euler backward time int for now

event friction(i++)
{
  cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
  foreach ()
  {
    double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
    foreach_dimension()
        u.x[] /= a;
  }
}

/**
Optionally, we can make a "stroboscopic" movie of the wave field. This
is useful to check the amount of waves reflected from the outflow. */
// ppm output does not look good
#if 0
event movie(t += 1)
{
  output_ppm(h, linear = true);
}
#endif

// gnuplot related
#if 0 
event contourData(t = 0; t <= simTime; t += 2)
{
  char name[80];
  sprintf(name, "OC");
  FILE *fp = fopen(name, "w");
  foreach ()
  {
    fprintf(fp, "%g %g %g\n", x, y, h[]);
  }
}

void plot_profile(double t, FILE *fp)
{
  fprintf(fp,
          "set term pngcairo enhanced size 600,600 font \",10\"\n"
          "set output 'g1_%.0f.png'\n"
          "set pm3d map\n"
          "set view map\n"
          "unset surface\n"
          // "set contour\n"
          "set title 't = %.2f'\n"
          "set size ratio -1\n"
          "set xlabel 'x(m)'\n"
          "set ylabel 'y(m)'\n"
          "splot 'OC' \n", //"plot [0:40][0:]'-' u 1:2 w l lw 2\n"
          t, t);

  // second strategy

  // foreach ()
  // {
  //   fprintf(fp, "%g %g %g\n", x, y, h[]);
  // }
  fprintf(fp, "e\n\n");
  fflush(fp);
}

event gnuplot(t = 0; t <= simTime; t += 2)
{
  static FILE *fp = popen("gnuplot 2> /dev/null", "w");
  plot_profile(t, fp);

  // fprintf(fp,
  //         "set term pngcairo enhanced size 800,600 font \",10\"\n"
  //         "set output 't%.0f.png'\n"
  //         "set title 't = %.2f'\n"
  //         "set xrange [0:40]\n"
  //         "plot u 1:2 w l t\n",
  //         t, t);
  // fprintf(fp, "\n");
  // foreach ()
  //      fprintf(fp, "%g %g\n", x, h[]);
  // fprintf(fp, "e\n\n");
  // fflush(fp);
  // fprintf(stderr, "%.3f %.3f\n", t, statsf(h).max); // uncomment if needed
}
#endif

/**
At the end of the simulation, we output the maximum wave amplitudes
along the cross-sections corresponding with the experimental data. */

event fieldOutput(t = 0; t <= simTime; t += 2)
{
  char name1[80];
  char name2[80];
  sprintf(name1, "out-h-%.0f", t);
  sprintf(name2, "out-u-%.0f", t);
  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");
  // note the usage of output_field!!! (x,y coordinate is output at the same time)
  // output_field({h, u, zb}, fp, linear = true);
  foreach ()
  {
    y <= 3. ? (fprintf(fp1, "%g %g %g\n", x, y, h[]), fprintf(fp2, "%g %g %g %g\n", x, y, u.x[], u.y[])) : none;
  }
  
  // compact form
  foreach ()
  {
    y <= 3. ? (fprintf(fp1, "%g %g %g\n", x, y, h[]), fprintf(fp2, %g %g\n", u.x[], u.y[])) : none;
  }

  // cross sectional data may be useful in the next stage
  // FILE *fp2 = fopen("section2", "w");
  // FILE *fp5 = fopen("section5", "w");
  // FILE *fp7 = fopen("section7", "w");
  // for (double y = -10.; y <= 10.; y += 0.02)
  // {
  //   fprintf(fp2, "%g %g\n", y, interpolate(maxa, 3., y));
  //   fprintf(stderr, "%g %g\n", y, interpolate(maxa, 5., y));
  //   fprintf(fp5, "%g %g\n", y, interpolate(maxa, 9., y));
  // }
  // for (double x = -10.; x <= 10.; x += 0.02)
  //   fprintf(fp7, "%g %g\n", x, interpolate(maxa, x, 0.));
}

/**
~~~gnuplot Instantaneous wave field at $t=50$
set output 'snapshot.png'
a0 = 0.026
splot [-10:12][-10:10]'end' u 1:2:($3/a0) w pm3d t ''
! mogrify -trim +repage snapshot.png
~~~

~~~gnuplot Maximum wave amplitude
set output 'maxa.png'
set label 1 "section 2" at 2.5,-6,1 rotate front
set label 2 "section 3" at 4.5,-6,1 rotate front
set label 3 "section 5" at 8.5,-6,1 rotate front
set label 4 "section 7" at -3,0.5,1 front
splot [-10:12][-10:10]'end' u 1:2:($5/a0) w pm3d t '', \
 '-' w l lt -1 t ''
3 -10 1
3 10 1


5 -10 1
5 10 1


9 -10 1
9 10 1


-10 0 1
10 0 1
e
unset label
! mogrify -trim +repage maxa.png
~~~

The results are comparable to [Lannes and Marche,
2014](/src/references.bib#lannes2014), Figure 19, but we have to use a
higher resolution (1024^2^ instead of 300^2^) because our numerical
scheme is only second order (versus 4th order).

~~~gnuplot Comparison of the maximum wave height with the experimental data (symbols) along various cross-sections
reset
set term svg enhanced size 640,480 font ",10"
set multiplot layout 2,2
set key top left
set xlabel 'y (m)'
set ylabel 'a_{max}/a_{0}'
set yrange [0:2.5]
plot [-5:5] '../section-2' pt 7 lc -1 t '', \
          'section2' u (-$1):($2/a0) w l lc 1 lw 2 t 'section 2'
plot [-5:5] '../section-3' pt 7 lc -1 t '', \
     'log' u (-$1):($2/a0) w l lc 1 lw 2 t 'section 3'
plot [-5:5] '../section-5' pt 7 lc -1 t '', \
     'section5' u (-$1):($2/a0) w l lc 1 lw 2 t 'section 5'
set xlabel 'x (m)'
plot [0:10] '../section-7' pt 7 lc -1 t '', \
     'section7' u 1:($2/a0) w l lc 1 lw 2 t 'section 7'
unset multiplot
~~~
*/
