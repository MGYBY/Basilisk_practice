#include "poisson.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
};

#if AXI
# define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
					    mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif
#endif

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {
    foreach_dimension()
      w.x[] = (dt/rho[]*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
               #if dimension > 1
			   + mu.y[0,1]*(u.x[0,1] +
					(u.y[1,0] + u.y[1,1])/4. -
					(u.y[-1,0] + u.y[-1,1])/4.)
			   - mu.y[]*(- u.x[0,-1] +
				     (u.y[1,-1] + u.y[1,0])/4. -
				     (u.y[-1,-1] + u.y[-1,0])/4.)
               #endif
	       #if dimension > 2
			   + mu.z[0,0,1]*(u.x[0,0,1] +
					  (u.z[1,0,0] + u.z[1,0,1])/4. -
					  (u.z[-1,0,0] + u.z[-1,0,1])/4.)
			   - mu.z[]*(- u.x[0,0,-1] +
				     (u.z[1,0,-1] + u.z[1,0,0])/4. -
				     (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
               #endif
			   ) + r.x[]*sq(Delta))/
    (sq(Delta)*lambda.x + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
                                    #if dimension > 1
				      + mu.y[0,1] + mu.y[]
                                    #endif
			            #if dimension > 2
				      + mu.z[0,0,1] + mu.z[]
			            #endif
			     ));
  }

#if JACOBI
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (u.x[] + 2.*w.x[])/3.;
#endif
  
#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  /**
  We manually apply boundary conditions, so that all components are
  treated simultaneously. Otherwise (automatic) BCs would be applied
  component by component before each foreach_face() loop. */
  
  boundary ({u});
  
  foreach_dimension() {
    face vector taux[];
    foreach_face(x)
      taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] + 
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] + 
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
    #endif
    foreach (reduction(max:maxres)) {
      double d = 0.;
      foreach_dimension()
	d += taux.x[1] - taux.x[];
      res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*d/Delta;
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres))
    foreach_dimension() {
      res.x[] = r.x[] - lambda.x*u.x[] +
        dt/rho[]*(2.*mu.x[1,0]*(u.x[1] - u.x[])
		  - 2.*mu.x[]*(u.x[] - u.x[-1])
        #if dimension > 1
		  + mu.y[0,1]*(u.x[0,1] - u.x[] +
			       (u.y[1,0] + u.y[1,1])/4. -
			       (u.y[-1,0] + u.y[-1,1])/4.)
		  - mu.y[]*(u.x[] - u.x[0,-1] +
			    (u.y[1,-1] + u.y[1,0])/4. -
			    (u.y[-1,-1] + u.y[-1,0])/4.)
	#endif
        #if dimension > 2
		  + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
				 (u.z[1,0,0] + u.z[1,0,1])/4. -
				 (u.z[-1,0,0] + u.z[-1,0,1])/4.)
		  - mu.z[]*(u.x[] - u.x[0,0,-1] +
			    (u.z[1,0,-1] + u.z[1,0,0])/4. -
			    (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
	#endif
		  )/sq(Delta);
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
#endif
  return maxres;
}

#undef lambda

trace
mgstats viscosity (vector u, face vector mu, scalar rho, double dt,
		   int nrelax = 4, scalar * res = NULL)
{
  vector r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];

  restriction ({mu,rho});
  struct Viscosity p = { mu, rho, dt };
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, nrelax, res);
}

trace
mgstats viscosity_explicit (vector u, face vector mu, scalar rho, double dt)
{
  vector r[];
  mgstats mg = {0};
  struct Viscosity p = { mu, rho, dt };
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      u.x[] += r.x[];
  return mg;
}

trace
mgstats viscosity_cn (vector u, face vector mu, scalar rho, double dt,
		   int nrelax = 4, scalar * res = NULL)
{
  vector r[], uexnext[];
  // explicit step
  mgstats mg = {0};
  struct Viscosity p = { mu, rho, dt };
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      uexnext.x[] = u.x[] + r.x[];
  // implicit step
  foreach()
    foreach_dimension()
      r.x[] = u.x[];
  restriction ({mu,rho});
  mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, nrelax, res);
  // CN
  foreach()
    foreach_dimension()
      u.x[] = 0.5*(uexnext.x[] + u.x[]);
  return mg;
}
