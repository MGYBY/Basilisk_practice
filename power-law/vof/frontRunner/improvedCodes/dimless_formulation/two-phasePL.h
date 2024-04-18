/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof.h"

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

double powerLawIndex = 1.0, muRef = 0.001, mumax = 1.0, tauP = 0.0;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */

  if (mu1 || mu2)
    mu = new face vector;

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
// # define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
# define mu(muTemp, mu2, f)  (clamp(f,0.,1.)*(muTemp - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

scalar yieldSurface[], strainRate[];

event tracer_advection (i++)
{

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

event properties (i++)
{
  // not for parallelism
  // double muTemp = mu1;
  foreach_face() {
    double muTemp = mu1;
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      #if dimension == 2
        double D11 = (u.x[] - u.x[-1,0]);
        double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
        double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
        double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
      #else // dimension == 3
        double D2 = 0.;
        double dxx = (u.x[1,0,0] - u.x[-1,0,0])/2.0;
        double dyy = ((u.y[0,1,0] - u.y[0,-1,0])/2.0+(u.y[-1,1,0] - u.y[-1,-1,0])/2.0)/2.0;
        double dzz = ((u.z[0,0,1] - u.z[0,0,-1])/2.0+(u.z[-1,0,1] - u.z[-1,0,-1])/2.0)/2.0;
        double dxy = 0.5*(((u.x[0,1,0]-u.x[0,-1,0])/2.0+(u.x[-1,1,0]-u.x[-1,-1,0])/2.0)/2.0+(u.y[]-u.y[-1,0,0]));
        double dxz = 0.5*(((u.x[0,0,1]-u.x[0,0,-1])/2.0+(u.x[-1,0,1]-u.x[-1,0,-1])/2.0)/2.0+(u.z[]-u.z[-1,0,0]));
        double dyz = 0.5*(((u.y[0,0,1]-u.y[0,0,-1])/2.0+(u.y[-1,0,1]-u.y[-1,0,-1])/2.0)/2.0+((u.z[0,1,0]-u.z[0,-1,0])/2.0+(u.z[-1,1,0]-u.z[-1,-1,0])/2.0)/2.0);
        D2 = sqrt(sq(dxx) + sq(dyy) +sq(dzz) +2.*sq(dxy) +2.*sq(dxz) +2.*sq(dyz))/Delta;
        // fprintf (ferr, "3D D2 has not been implemented yet");
        // exit (1);
      #endif
      strainRate[] = D2;
      if (D2 > 0.) {
        D2 = max(D2, 1.0e-45);
//         double temp = muRef * exp((powerLawIndex - 1.) * log(D2 * pow(2,0.5)));
//         double temp = muRef * exp((powerLawIndex - 1.) * log(D2 * pow(2,0.5))) + tauP/((D2 * pow(2,0.5))+1.0e-15);
        double temp = muRef * exp((powerLawIndex - 1.) * log(D2 * pow(2,0.5))) + tauP/((D2 * pow(2,0.5)));
//         double temp = muRef * exp((powerLawIndex - 1.) * log(D2));
//         m = MUREF*exp ((N - 1.)*log (d2*pow(2,0.5)));
        yieldSurface[] = temp>=mumax ? 1.0*f[] : 0.0*f[];
        muTemp = min(temp, mumax);
      } else {
        if (powerLawIndex<1.){
          muTemp = mumax;
        } else {
          muTemp = mu1;
        }
      }
      /**
      Note that only the heavier fluid is Viscoplastic.
      */
      muv.x[] = fm.x[]*mu(muTemp, mu2, ff);
    }
  }

  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
