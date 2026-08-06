/**
# Total-momentum-consistent VOF advection for nested three phases

This is the AMR-compatible stratified counterpart of
`navier-stokes/conserving.h`.  With cumulative fractions `f_lower` and
`f`, total momentum is decomposed as

  rho1 f_lower u
+ rho2 f u
- rho2 f_lower u
+ rho3 (1 - f) u.

Each term is transported as a VOF-confined tracer of the corresponding
cumulative interface.  Their signed sum is the physical total momentum
of lower liquid, upper liquid and air.
*/

#ifndef BASILISK_CONSERVING3F_STRATIFIED_H
#define BASILISK_CONSERVING3F_STRATIFIED_H

#if TREE
static void stratified_momentum_refine (Point point, scalar uc)
{
  refine_bilinear (point, uc);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*stratified_density_cell(f_lower[], f[])*uc[];
  const double parent_rho = stratified_density_cell(f_lower[], f[]);
  const double du = uc[] - rhou/
    ((1 << dimension)*(cm[] + SEPS)*parent_rho);
  foreach_child()
    uc[] += du;
}

static void stratified_momentum_restriction (Point point, scalar uc)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*stratified_density_cell(f_lower[], f[])*uc[];
  uc[] = rhou/((1 << dimension)*(cm[] + SEPS)*
               stratified_density_cell(f_lower[], f[]));
}

static void stratified_move_to_front (scalar c)
{
  int j = 0;
  while (all[j].i >= 0 && all[j].i != c.i)
    j++;
  if (all[j].i < 0)
    return;
  while (j > 0)
    all[j] = all[j - 1], j--;
  all[0] = c;
}
#endif

event defaults (i = 0)
{
  stokes = true; // switch off centered.h's ordinary velocity advection

#if TREE
  stratified_move_to_front (f);
  stratified_move_to_front (f_lower);

  foreach_dimension() {
    u.x.refine = u.x.prolongation = stratified_momentum_refine;
    u.x.restriction = stratified_momentum_restriction;
    u.x.depends = list_add (u.x.depends, f_lower);
    u.x.depends = list_add (u.x.depends, f);
  }
#endif
}

event stability (i++)
  dtmax = timestep (uf, dtmax);

foreach_dimension()
static double boundary_q1_x (Point neighbor, Point point, scalar q,
                              bool * data)
{
  double a1, a2, a3;
  stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
  return a1*rho1*u.x[];
}

foreach_dimension()
static double boundary_q2lower_x (Point neighbor, Point point, scalar q,
                                   bool * data)
{
  double a1, a2, a3;
  stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
  return a1*rho2*u.x[];
}

foreach_dimension()
static double boundary_q2total_x (Point neighbor, Point point, scalar q,
                                   bool * data)
{
  const double liquid = st_clip01(f[]);
  return liquid*rho2*u.x[];
}

foreach_dimension()
static double boundary_q3_x (Point neighbor, Point point, scalar q,
                              bool * data)
{
  const double liquid = st_clip01(f[]);
  return (1. - liquid)*rho3*u.x[];
}

#if TREE
foreach_dimension()
static void prolongation_q1_x (Point point, scalar q)
{
  foreach_child() {
    double a1, a2, a3;
    stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
    q[] = a1*rho1*u.x[];
  }
}

foreach_dimension()
static void prolongation_q2lower_x (Point point, scalar q)
{
  foreach_child() {
    double a1, a2, a3;
    stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
    q[] = a1*rho2*u.x[];
  }
}

foreach_dimension()
static void prolongation_q2total_x (Point point, scalar q)
{
  foreach_child()
    q[] = st_clip01(f[])*rho2*u.x[];
}

foreach_dimension()
static void prolongation_q3_x (Point point, scalar q)
{
  foreach_child()
    q[] = (1. - st_clip01(f[]))*rho3*u.x[];
}
#endif

static scalar * stratified_interfaces_saved = NULL;

event vof (i++)
{
  vector q1[], q2lower[], q2total[], q3[];

  for (scalar s in {q1, q2lower, q2total, q3}) {
    s.depends = list_add (s.depends, f_lower);
    s.depends = list_add (s.depends, f);
    foreach_dimension()
      s.v.x.i = -1;
  }

  for (int b = 0; b < nboundary; b++)
    foreach_dimension() {
      q1.x.boundary[b] = boundary_q1_x;
      q2lower.x.boundary[b] = boundary_q2lower_x;
      q2total.x.boundary[b] = boundary_q2total_x;
      q3.x.boundary[b] = boundary_q3_x;
    }

#if TREE
  foreach_dimension() {
    q1.x.prolongation = prolongation_q1_x;
    q2lower.x.prolongation = prolongation_q2lower_x;
    q2total.x.prolongation = prolongation_q2total_x;
    q3.x.prolongation = prolongation_q3_x;
  }
#endif

  foreach()
    foreach_dimension() {
      double a1, a2, a3;
      stratified_phase_values (f_lower[], f[], &a1, &a2, &a3);
      q1.x[] = a1*rho1*u.x[];
      q2lower.x[] = a1*rho2*u.x[];
      q2total.x[] = (a1 + a2)*rho2*u.x[];
      q3.x[] = a3*rho3*u.x[];
    }

  foreach_dimension() {
    q3.x.inverse = true;
    q1.x.gradient = q2lower.x.gradient = q2total.x.gradient =
      q3.x.gradient = u.x.gradient;
  }

  scalar * old_lower = f_lower.tracers;
  scalar * old_total = f.tracers;
  f_lower.tracers = list_concat (old_lower,
                                  (scalar *){q1, q2lower});
  f.tracers = list_concat (old_total,
                            (scalar *){q2total, q3});

  vof_advection ({f_lower, f}, i);

  free (f_lower.tracers);
  free (f.tracers);
  f_lower.tracers = old_lower;
  f.tracers = old_total;

  stratified_enforce_hierarchy();

  foreach()
    foreach_dimension()
      u.x[] = (q1.x[] + q2total.x[] - q2lower.x[] + q3.x[])/
        stratified_density_cell(f_lower[], f[]);
  boundary ((scalar *){u});

  stratified_interfaces_saved = interfaces;
  interfaces = NULL;
}

event tracer_advection (i++)
{
  interfaces = stratified_interfaces_saved;
}

#endif
