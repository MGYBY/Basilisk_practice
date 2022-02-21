//#include "grid/multigrid.h" // multigrid cannot be used for mask
//#include "grid/cartesian.h"
// #include "green-naghdi.h"
#include "saint-venant-power-law.h"

#define MAXLEVEL 8
#define MINLEVEL 6
#define MAXMAXLEVEL 14

// problem-sepcific parameters
double n_coeff = 0.40;
double alpha_coeff = 1.0;
double lx = 0.0;
double ly = 0.0;
double disMag = 0.025;
// double betaCoeff = 0.0;
double simTime = 25.0;

// for debugging
double blockXLower = 0.0;
double blockXUpper = 0.0;
double blockYLower = 0.0;
double blockYUpper = 0.0;
double blockWidth = 0.0;
// avoid DBZ for U
double epsilonU = 1e-10;

double blockXLowerDrag = 0.0;
double blockXUpperDrag = 0.0;

// homogeneous bc
bid block;
//free-slip bc
//not working???
// u.n[block] = dirichlet(0);


int main()
{
  blockWidth = 3.0*1.0;
  lx = blockWidth*7.0;
  // square domain
  ly = lx;
  // lx = 7.0;
  blockXLower = lx/2.0+2.0;
  blockXUpper = lx/2.0+2.0+blockWidth;
  blockYLower = lx/2.0-blockWidth/2.0;
  blockYUpper = lx/2.0+blockWidth/2.0;
  size(lx);
  alphaCoeff = alpha_coeff;
  betaCoeff = (2.0*(1.0+2.0*n_coeff))/(2.0+3.0*n_coeff);
  // G = 9.81;
  init_grid(1 << MAXLEVEL);
  CFL = 0.2; // CFL number should be sufficiently small
//   theta = 1.5; // use Sweby limiter
  // gradient = NULL;

  run();
}

event init(i = 0)
{

  char name[80];

  mask(x >= blockXLower ? (x <= blockXUpper ? (y >= blockYLower ? (y <= blockYUpper ? block : none) : none) : none) : none);

  // use this trick to save storage if needed
  // refine(((x >= 0.0) && (x <= rightCoord+0.50)) && level < 14);

  h[left] = dirichlet(1.0);
  u.n[left] = dirichlet(1.0);

  u.n[right] = neumann(0.);
  h[right] = neumann(0.);

  u.n[top] = neumann(0.);
  h[top] = neumann(0.);

  u.n[bottom] = neumann(0.);
  h[bottom] = neumann(0.);

  foreach ()
  {
    zb[] = 0.0;
    h[] = 1.0;
    u.x[] = 1.0;
    u.y[] = 0.0;
  }

   foreach_boundary (block){
     if (x<=blockXLower)
     {
       blockXLowerDrag = x;
     }

     if (x>=blockXUpper)
     {
       blockXUpperDrag = x;
     }
   }

  //     sprintf(name, "out-%.3f.vtk", t);
      sprintf(name, "out-%.1f.gfs", t);
      FILE *fp = fopen(name, "w");
  //     output_vtk((scalar *) {h}, 1 << 10, (FILE *) fp, false);
      output_gfs(fp, translate = true);
      sprintf(name, "outText-%.1f.txt", t);
      FILE *fp1 = fopen(name, "w");
      foreach ()
  {
    fprintf(fp1, "%g %g %g\n", x, y, h[]);
  }
  fclose(fp);
  fclose(fp1);
}

static double powerLawFrictionX(double u, double h, double n)
{
     double rhs;
     // rhs = h-pow((u/h), n);
     rhs = h-u*pow(fabs(u), (n-1))/pow(h, n);
     return rhs;
}

static double powerLawFrictionY(double v, double h, double n)
{
     double rhs;
     // rhs = -pow((v/h), n);
     rhs = -v*pow(fabs(v), (n-1))/pow(h, n);
     return rhs;
}

event friction(i++)
{

  foreach ()
  {

    //linearized backeard Euler
    // u.x[] = (h[] > dry && u.x[]>epsilonU) ? (u.x[] + dt)/(1.0+dt*(1.0/h[])*((pow(u.x[], (n_coeff-1.0)))/(pow(h[], n_coeff)))) : 0.0;
    // u.y[] = (h[] > dry && u.y[]>epsilonU) ? (u.y[])/(1.0+dt*(1.0/h[])*((pow(u.y[], (n_coeff-1.0)))/(pow(h[], n_coeff)))) : 0.0;

    // rk2tvd
    if (h[] > dry) {
        double uMed = 0.0;
        double vMed = 0.0;
        if (fabs(u.x[])>epsilonU) {
          uMed = u.x[] + dt * powerLawFrictionX(u.x[], h[], n_coeff);
          u.x[] = 0.5*u.x[]+0.5*uMed+0.5*dt*powerLawFrictionX(uMed, h[], n_coeff);
        }
        else {
          u.x[] = 0.0;
        }

        if (fabs(u.y[])>epsilonU) {
          vMed = u.y[] + dt * powerLawFrictionY(u.y[], h[], n_coeff);
          u.y[] = 0.5*u.y[]+0.5*vMed+0.5*dt*powerLawFrictionY(vMed, h[], n_coeff);
        }
        else {
          u.y[] = 0.0;
        }
    }
    else {
         u.x[] = 0.0;
         u.y[] = 0.0;
    }
  }

  boundary((scalar *){u.x});
  boundary((scalar *){u.y});
}

// event outputDrag(i++)
event outputDrag(i+=10) // for steady-state calculation
{

  double fd = 0.0;
  double rmax = 0.0;

  FILE *fp1 = fopen("dragForce", "a+");
  FILE *fp2 = fopen("rMax", "a+");

  foreach_boundary (block){
    if (x==blockXLowerDrag)
    {
      fd += 0.5*alpha_coeff*(Delta*h[])*h[];
    }
    else if (x==blockXUpperDrag)
    {
      fd -= 0.5*alpha_coeff*(Delta*h[])*h[];
    }
  }
  fprintf(fp1, "%g %g %g\n", t, fd, fd/(0.5*alpha_coeff*(blockYUpper-blockYLower)*sq(1.0)*1.0));
  rmax = interpolate(h, blockXLower-lx/(pow(2, 8)), ly/2.0);
  fprintf(fp2, "%g %g \n", t, rmax);

  fclose(fp1);
  fclose(fp2);
}

event fieldOutput(t = 0; t <= simTime; t += 1)
{

  char name[80];
//     sprintf(name, "out-%.3f.vtk", t);
    sprintf(name, "out-%.1f.gfs", t);
    FILE *fp = fopen(name, "w");
//     output_vtk((scalar *) {h}, 1 << 10, (FILE *) fp, false);
    output_gfs(fp, translate = true);
    // output_gfs(fp, translate = false);
    fclose(fp);

    sprintf(name, "outText-%.1f.txt", t);
    FILE *fp1 = fopen(name, "w");
    foreach ()
    {
      fprintf(fp1, "%g %g %g\n", x, y, h[]);
      }
    fclose(fp1);

    // text output below
//   // ignore I.C.
// //   char name1[24];
//   char name2[24];
//   char name3[24];
//   char name4[32];
// //   sprintf(name1, "out-h-%.0f", t);
//   sprintf(name2, "out-u-zoomIn-%.0f", t);
//   sprintf(name3, "out-h-zoomIn-%.0f", t);
//   sprintf(name4, "profile-%.0f", t);
// //   FILE *fp1 = fopen(name1, "w");
//   FILE *fp2 = fopen(name2, "w");
//   FILE *fp3 = fopen(name3, "w");
//   FILE *fp4 = fopen(name4, "w");
//// event adapt(i++)
// {
//   astats s = adapt_wavelet({h}, (double[]){1.0 / 150.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
// }
//   // compact form
//   foreach ()
//   {
// //     fprintf(fp1, "%g %g %g\n", x, y, h[]);
// //     fprintf(fp2, "%g %g\n", u.x[], u.y[]);
//     if (x>150.0)
//     {
//       fprintf(fp3, "%g %g %g\n", x, y, h[]);
//       fprintf(fp2, "%g %g %g %g\n", x, y, u.x[], u.y[]);
//     }
//   }
//
//   for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
//   { if (interpolate(h, x, Ly / 2.)<10.0)
//       {
//       fprintf(fp4, "%g %g\n", x, interpolate(h, x, Ly / 2.));
//       }

  // }

  // for (double x = 0.0; x <= Lx; x += (Lx / (pow(2, MAXMAXLEVEL))))
  //   fprintf(fp3, "%g %g\n", x, interpolate(h, x, Ly / 2.));

//   fclose(fp1);
  // fclose(fp2);
  // fclose(fp3);
  // fclose(fp4);
}

// adaptivity
// event adapt(i++)
// {
//   astats s = adapt_wavelet({h}, (double[]){1.0 / 150.0}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
//   fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
// }
