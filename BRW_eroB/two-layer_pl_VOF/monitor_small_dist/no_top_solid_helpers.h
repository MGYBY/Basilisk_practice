/** Ordinary-C helpers for the no-embedded-ceiling geometry (v2p7, retained from v2p4/v2p5).
 * These functions prescribe initialization and a mesh-level ceiling only;
 * they never reset evolved velocity or pressure in the far air.
 */
#ifndef FRONT_RUNNER_NO_TOP_SOLID_HELPERS_H
#define FRONT_RUNNER_NO_TOP_SOLID_HELPERS_H
#include <math.h>

#ifndef CASE_FAR_AIR_CELL_GROWTH
# define CASE_FAR_AIR_CELL_GROWTH 0.25
#endif

/* Start coarse, but allow at least one initial leaf per rank when the
   requested INITLEVEL permits it. This is a Cartesian 2-D case. */
static inline int nts_initial_level (int minlevel, int initlevel, int ranks)
{
  int seed = minlevel;
  while (seed < initlevel && ldexp(1., 2*seed) < ranks)
    seed++;
  return seed;
}

static inline double nts_air_pressure (double y, double reference_height,
                                       double reference_pressure,
                                       double air_density, double normal_gravity)
{
  return reference_pressure - air_density*normal_gravity*(y-reference_height);
}

static inline int nts_maxlevel (double y, double delta, double near_height,
                                double fine_fraction, double finest,
                                int maxlevel, int minlevel, int level_drop,
                                double cell_growth)
{
  if (y + 0.5*delta <= near_height*fine_fraction)
    return maxlevel;
  int cap = maxlevel - level_drop;
  if (cap < minlevel) cap = minlevel;
  if (cap > maxlevel) cap = maxlevel;
  double cell_size = ldexp(finest, maxlevel-cap);
  const double distance = fmax(0., y-0.5*delta-near_height);
  const double target = cell_size + fmax(0., cell_growth)*distance;
  /* Integer powers of two avoid a log/rounding ambiguity at level changes. */
  while (cap > minlevel && 2.*cell_size <= target) {
    cell_size *= 2.;
    cap--;
  }
  return cap;
}
#endif
