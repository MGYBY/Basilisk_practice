/**
# One-dimensional hydrostatic multilayer Saint-Venant solver for thixotropic fluids

This header implements the user's dimensionless hydrostatic multilayer model for
thixotropic fluids,

$$
\partial_t h + \sum_{\alpha=0}^{nl-1} \partial_x (h_\alpha u_\alpha) = 0,
$$

$$
\partial_t (h_\alpha u_\alpha) + \partial_x\left(h_\alpha u_\alpha^2 +
\frac{G}{2} l_\alpha h^2\right)
= G h_\alpha + u_{\alpha + 1/2} G_{\alpha + 1/2} -
  u_{\alpha - 1/2} G_{\alpha - 1/2} + \text{vertical viscosity},
$$

$$
\partial_t (h_\alpha \lambda_\alpha) +
\partial_x (h_\alpha u_\alpha \lambda_\alpha)
= \lambda_{\alpha + 1/2}^{up} G_{\alpha + 1/2} -
  \lambda_{\alpha - 1/2}^{up} G_{\alpha - 1/2} +
  h_\alpha \left[\frac{1 - \lambda_\alpha}{T} -
  \frac{\Gamma}{T}\lambda_\alpha \left|\frac{\partial u_\alpha}{\partial z}\right|\right]
+ \text{vertical diffusion}.
$$

The design follows the hydrostatic/Riemann-solver pathway of Basilisk's
`saint-venant.h` and `multilayer.h`, but is restricted to the one-dimensional
periodic roll-wave setting.  The code keeps the structure of the Basilisk solver:

* hyperbolic update for `h` and `h u_\alpha` through `kurganov()`;
* conservative transport of `h \lambda_\alpha` using the same mass flux;
* identical donor-cell logic for the interlayer `G_{\alpha + 1/2}` terms in the
  momentum and structure equations;
* implicit vertical viscosity and implicit vertical diffusion of `\lambda`;
* optional local relaxation integrators, with SSP RK2 as the default.
*/

scalar zb[], h[], eta[];
vector u[];          // first layer horizontal velocity
scalar lambda[];     // first layer structure variable
scalar muifield[];   // diagnostic: max interface viscosity in a column

/** Gravity and dry tolerance follow the standard Saint-Venant solver. */

double G = 1.;
double dry = 1e-10;

#if !LAYERS
int nl = 1;
#endif

/** Lists of layer fields.  `ul[l]` and `lambdal[l]` are the velocity and
    structure variables of layer `l`.  `divl[l]` stores the layer divergence used
    to assemble the interlayer exchange velocity `G_{l + 1/2}`. */

vector * ul = NULL;
scalar * lambdal = NULL;
scalar * divl = NULL;
double * layer = NULL;

attribute {
  int l;
}

/** Dimensionless thixotropic parameters. */

double thixo_T = 1.;
double thixo_Gamma = 8.;
double thixo_kappa = 1e-4;
double thixo_a = 0.6;
double thixo_lambda_init = 1.;
double thixo_Aeps = 1e-10;
double thixo_mu_max = 1e10;

enum {
  THIXO_RELAX_EXACT = 0,
  THIXO_RELAX_EULER = 1,
  THIXO_RELAX_SSPRK2 = 2
};

int thixo_relaxation_scheme = THIXO_RELAX_SSPRK2;
bool thixo_use_exact_relaxation = false;

/** Predictor-corrector integrates these fields. */

scalar * evolving = NULL;

static inline double clamp01 (double x)
{
  return x < 0. ? 0. : (x > 1. ? 1. : x);
}

static inline double A_thixo (double lambda)
{
  double om = 1. - clamp01 (lambda);
  return max (om*((1. - thixo_a)*om + thixo_a), thixo_Aeps);
}

static inline double mu_thixo (double lambda)
{
  return min (1./A_thixo (lambda), thixo_mu_max);
}

static inline double interface_lambda (Point point, scalar * ll, int i)
{
  if (i <= 0) {
    scalar s = ll[0];
    return clamp01 (s[]);
  }
  if (i >= nl) {
    scalar s = ll[nl - 1];
    return clamp01 (s[]);
  }
  scalar sb = ll[i - 1], st = ll[i];
  return clamp01 (0.5*(sb[] + st[]));
}

/** The dominant long-wave strain rate is approximated by the vertical shear
    `|du/dz|`.  The estimate is centered in the interior and one-sided near the
    bottom and free surface. */

static inline double shear_rate_x (Point point, double H, vector * uu, int l)
{
  if (H <= dry)
    return 0.;
  if (nl == 1) {
    vector uk = uu[0];
    return fabs (uk.x[])/(0.5*max (layer[0]*H, dry));
  }
  if (l == 0) {
    vector u0 = uu[0], u1 = uu[1];
    double dz = max (0.5*(layer[0] + layer[1])*H, dry);
    return fabs ((u1.x[] - u0.x[])/dz);
  }
  if (l == nl - 1) {
    vector ub = uu[nl - 2], uk = uu[nl - 1];
    double dz = max (0.5*(layer[nl - 2] + layer[nl - 1])*H, dry);
    return fabs ((uk.x[] - ub.x[])/dz);
  }
  vector ub = uu[l - 1], ut = uu[l + 1];
  double dz = max ((0.5*(layer[l - 1] + layer[l + 1]) + layer[l])*H, dry);
  return fabs ((ut.x[] - ub.x[])/dz);
}

/**
## Implicit vertical viscosity for momentum

This is the thixotropic counterpart of `vertical_viscosity()` in Basilisk's
`multilayer.h`.  The interface viscosity is evaluated from the layer-averaged
`\lambda` values and the resulting tridiagonal system is solved by the Thomas
algorithm.
*/

trace
static void vertical_viscosity_thixo (Point point, double H,
                                      vector * uu, scalar * ll, double dt)
{
  if (H <= dry)
    return;

  double muif[nl + 1];
  for (int i = 0; i <= nl; i++)
    muif[i] = mu_thixo (interface_lambda (point, ll, i));

  double mumax = muif[0];
  for (int i = 1; i <= nl; i++)
    if (muif[i] > mumax)
      mumax = muif[i];
  muifield[] = mumax;

  double a[nl], b[nl], c[nl], rhs[nl];

  for (int l = 0; l < nl; l++) {
    vector uk = uu[l];
    rhs[l] = H*layer[l]*uk.x[];
    a[l] = c[l] = 0.;
    b[l] = H*layer[l];
  }

  if (nl == 1) {
    b[0] += 2.*muif[0]*dt/max (H*layer[0], dry); // no-slip at the bottom
  }
  else {
    for (int l = 1; l < nl - 1; l++) {
      a[l] = - 2.*muif[l]*dt/(H*(layer[l - 1] + layer[l]));
      c[l] = - 2.*muif[l + 1]*dt/(H*(layer[l] + layer[l + 1]));
      b[l] -= a[l] + c[l];
    }
    c[0] = - 2.*muif[1]*dt/(H*(layer[0] + layer[1]));
    b[0] += 2.*muif[0]*dt/max (H*layer[0], dry) - c[0];
    a[nl - 1] = - 2.*muif[nl - 1]*dt/(H*(layer[nl - 2] + layer[nl - 1]));
    b[nl - 1] -= a[nl - 1];
  }

  for (int l = 1; l < nl; l++) {
    double m = a[l]/b[l - 1];
    b[l] -= m*c[l - 1];
    rhs[l] -= m*rhs[l - 1];
  }

  vector uk = uu[nl - 1];
  uk.x[] = rhs[nl - 1]/b[nl - 1];
  double back = uk.x[];
  for (int l = nl - 2; l >= 0; l--) {
    vector uloc = uu[l];
    back = (rhs[l] - c[l]*back)/b[l];
    uloc.x[] = back;
  }
}

/**
## Implicit vertical diffusion of `lambda`

The diffusion term is written as the difference of interface diffusive fluxes,
which leads to another tridiagonal solve.  Since `h` is frozen over this substep,
solving for the layer averages `\lambda_\alpha` is equivalent to solving for the
conservative variables `h_\alpha \lambda_\alpha`.
*/

trace
static void vertical_diffusion_lambda (Point point, double H, scalar * ll, double dt)
{
  if (thixo_kappa <= 0. || H <= dry || nl == 1)
    return;

  double D = thixo_kappa/thixo_T;
  if (D <= 0.)
    return;

  double a[nl], b[nl], c[nl], rhs[nl];
  for (int l = 0; l < nl; l++) {
    scalar lam = ll[l];
    rhs[l] = H*layer[l]*lam[];
    a[l] = c[l] = 0.;
    b[l] = H*layer[l];
  }

  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(H*(layer[l - 1] + layer[l]));
    c[l] = - 2.*D*dt/(H*(layer[l] + layer[l + 1]));
    b[l] -= a[l] + c[l];
  }
  c[0] = - 2.*D*dt/(H*(layer[0] + layer[1]));
  b[0] -= c[0];
  a[nl - 1] = - 2.*D*dt/(H*(layer[nl - 2] + layer[nl - 1]));
  b[nl - 1] -= a[nl - 1];

  for (int l = 1; l < nl; l++) {
    double m = a[l]/b[l - 1];
    b[l] -= m*c[l - 1];
    rhs[l] -= m*rhs[l - 1];
  }

  scalar lam = ll[nl - 1];
  lam[] = clamp01 (rhs[nl - 1]/b[nl - 1]);
  double back = lam[];
  for (int l = nl - 2; l >= 0; l--) {
    scalar sl = ll[l];
    back = (rhs[l] - c[l]*back)/b[l];
    sl[] = clamp01 (back);
  }
}

static inline double lambda_rhs_local (double lambda, double gdot)
{
  return (1. - lambda - thixo_Gamma*lambda*gdot)/thixo_T;
}

/** Local relaxation/source update for `\lambda`. */

trace
static void relax_structure (Point point, double H, vector * uu, scalar * ll, double dt)
{
  if (H <= dry)
    return;

  int scheme = thixo_use_exact_relaxation ? THIXO_RELAX_EXACT :
    thixo_relaxation_scheme;

  for (int l = 0; l < nl; l++) {
    scalar lam = ll[l];
    double lam0 = clamp01 (lam[]);
    double gdot = shear_rate_x (point, H, uu, l);

    if (scheme == THIXO_RELAX_EXACT) {
      double rate = (1. + thixo_Gamma*gdot)/thixo_T;
      double leq = 1./(1. + thixo_Gamma*gdot);
      lam[] = clamp01 (leq + (lam0 - leq)*exp (- dt*rate));
    }
    else if (scheme == THIXO_RELAX_SSPRK2) {
      double k1 = lambda_rhs_local (lam0, gdot);
      double lam1 = clamp01 (lam0 + dt*k1);
      double k2 = lambda_rhs_local (lam1, gdot);
      lam[] = clamp01 (0.5*lam0 + 0.5*(lam1 + dt*k2));
    }
    else
      lam[] = clamp01 (lam0 + dt*lambda_rhs_local (lam0, gdot));
  }
}

/** Interlayer momentum fluxes using the standard donor-cell choice. */

trace
static void vertical_fluxes_momentum (vector * uu, vector * duu,
                                      scalar * div, scalar dh)
{
  foreach() {
    double Gi = 0.;
    for (int l = 0; l < nl - 1; l++) {
      scalar d = div[l];
      Gi += d[] + layer[l]*dh[];
      vector ub = uu[l], ut = uu[l + 1];
      double ui = Gi < 0. ? ub.x[] : ut.x[];
      double flux = Gi*ui;
      vector dub = duu[l], dut = duu[l + 1];
      dub.x[] += flux/layer[l];
      dut.x[] -= flux/layer[l + 1];
    }
  }
}

/** Interlayer conservative fluxes for `h lambda`. */

trace
static void vertical_fluxes_lambda (scalar * ll, scalar * dll,
                                    scalar * div, scalar dh)
{
  foreach() {
    double Gi = 0.;
    for (int l = 0; l < nl - 1; l++) {
      scalar d = div[l];
      Gi += d[] + layer[l]*dh[];
      scalar lb = ll[l], lt = ll[l + 1];
      double li = Gi < 0. ? clamp01 (lb[]) : clamp01 (lt[]);
      double flux = Gi*li;
      scalar dlb = dll[l], dlt = dll[l + 1];
      dlb[] += flux/layer[l];
      dlt[] -= flux/layer[l + 1];
    }
  }
}

#include "predictor-corrector.h"

/** Convert updates of `h`, `h u` and `h lambda` into updates of `h`, `u` and
    `lambda`.  Vertical viscosity, local relaxation and vertical diffusion are
    applied after the conservative hyperbolic update. */

trace
static void advance_saint_venant_thixo (scalar * soutput, scalar * sinput,
                                        scalar * updates, double dt)
{
  scalar hi = sinput[0], ho = soutput[0], dh = updates[0];
  vector * uo = (vector *) &soutput[1];
  vector * ui = (vector *) &sinput[1];
  scalar * llo = &soutput[1 + nl];
  scalar * lli = &sinput[1 + nl];
  scalar * dll = &updates[1 + nl];

  foreach() {
    double hold = hi[];
    ho[] = hold + dt*dh[];
    eta[] = zb[] + ho[];

    if (ho[] > dry) {
      for (int l = 0; l < nl; l++) {
        vector uol = uo[l], uil = ui[l];
        vector dhu = ((vector *) &updates[1])[l];
        uol.x[] = (hold*uil.x[] + dt*dhu.x[])/ho[];

        scalar lamo = llo[l], lami = lli[l], dql = dll[l];
        lamo[] = clamp01 ((hold*lami[] + dt*dql[])/ho[]);
      }
      if (nl > 1)
        vertical_viscosity_thixo (point, ho[], uo, llo, dt);
      else
        muifield[] = mu_thixo (interface_lambda (point, llo, 0));
      relax_structure (point, ho[], uo, llo, dt);
      if (nl > 1)
        vertical_diffusion_lambda (point, ho[], llo, dt);
    }
    else {
      for (int l = 0; l < nl; l++) {
        vector uloc = uo[l];
        scalar lam = llo[l];
        uloc.x[] = 0.;
        lam[] = thixo_lambda_init;
      }
      muifield[] = 0.;
    }
  }
}

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h[];
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[] + h[];
}
#endif

/**
## Hyperbolic update

The mass and momentum fluxes are those of the one-dimensional hydrostatic
multilayer Saint-Venant system.  The hyperbolic flux for `h lambda` is built from
exactly the same numerical mass flux as a conservative passive scalar.
*/

#include "riemann.h"

trace
static double update_saint_venant_thixo (scalar * evolving, scalar * updates,
                                         double dtmax)
{
  scalar hh = evolving[0], dh = updates[0];
  face vector Fh[], Fl[], S[];
  face vector Fq[];

  vector gh[], geta[];
  for (scalar s in {gh, geta}) {
    s.gradient = zero;
#if TREE
    s.prolongation = refine_linear;
#endif
  }
  gradients ({hh, eta}, {gh, geta});

  for (int l = 0; l < nl; l++) {
    vector uk = ((vector *) &evolving[1])[l];
    scalar lam = evolving[1 + nl + l];

    vector gu[], glam[];
    for (scalar s in {gu, glam}) {
      s.gradient = zero;
#if TREE
      s.prolongation = refine_linear;
#endif
    }
    gradients ((scalar *) {uk}, (vector *) {gu});
    gradients ((scalar *) {lam}, (vector *) {glam});

    foreach_face (x, reduction (min:dtmax)) {
      double hi = hh[], hn = hh[-1];
      if (hi > dry || hn > dry) {
        double dx = Delta/2.;
        double zi = eta[] - hi;
        double zl = zi - dx*(geta.x[] - gh.x[]);
        double zn = eta[-1] - hn;
        double zr = zn + dx*(geta.x[-1] - gh.x[-1]);
        double zlr = max (zl, zr);

        double hl = hi - dx*gh.x[];
        double hp = max (0., hl + zl - zlr);
        double up = uk.x[] - dx*gu.x[];

        double hr = hn + dx*gh.x[-1];
        double hm = max (0., hr + zr - zlr);
        double um = uk.x[-1] + dx*gu.x[-1];

        double fh, fu;
        kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);

#if TREE
        if (is_prolongation(cell)) {
          hi = coarse (hh);
          zi = coarse (zb);
        }
        if (is_prolongation(neighbor(-1))) {
          hn = coarse (hh, -1);
          zn = coarse (zb, -1);
        }
#endif

        double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
        double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));

        Fh.x[] = fm.x[]*fh;
        Fq.x[] = fm.x[]*(fu - sl);
        S.x[]  = fm.x[]*(fu - sr);

        double lamp = clamp01 (lam[]  - dx*glam.x[]);
        double lamm = clamp01 (lam[-1] + dx*glam.x[-1]);
        Fl.x[] = fm.x[]*fh*(fh > 0. ? lamm : lamp);
      }
      else
        Fh.x[] = Fq.x[] = S.x[] = Fl.x[] = 0.;
    }

    vector dhu = ((vector *) &updates[1])[l];
    scalar dql = updates[1 + nl + l];
    double layerl = layer[l];

    foreach() {
      double dhl = layerl*(Fh.x[1] - Fh.x[])/(cm[]*Delta);
      dh[] = - dhl + (l > 0 ? dh[] : 0.);
      dhu.x[] = (Fq.x[] - S.x[1])/(cm[]*Delta);
      dql[] = - (Fl.x[1] - Fl.x[])/(cm[]*Delta);

      if (l < nl - 1) {
        scalar d = divl[l];
        d[] = dhl;
      }
    }
  }

  if (nl > 1) {
    vertical_fluxes_lambda (&evolving[1 + nl], &updates[1 + nl], divl, dh);
    vertical_fluxes_momentum ((vector *) &evolving[1], (vector *) &updates[1],
                              divl, dh);
  }

  return dtmax;
}

/**
## Initialisation and cleanup
*/

event defaults (i = 0)
{
  assert (ul == NULL && lambdal == NULL && divl == NULL);
  assert (nl > 0);

  ul = vectors_append (ul, u);
  lambdal = list_append (lambdal, lambda);
  for (int l = 1; l < nl; l++) {
    vector uk = new vector;
    scalar lam = new scalar;
    uk.x.l = l;
    lam.l = l;
    ul = vectors_append (ul, uk);
    lambdal = list_append (lambdal, lam);
  }
  for (int l = 0; l < nl - 1; l++) {
    scalar d = new scalar;
    d.l = l;
    divl = list_append (divl, d);
  }

  evolving = list_concat ({h}, (scalar *) ul);
  evolving = list_concat (evolving, lambdal);

  foreach() {
    h[] = zb[] = eta[] = 0.;
    muifield[] = 0.;
    for (scalar s in (scalar *) ul)
      s[] = 0.;
    for (scalar lam in lambdal)
      lam[] = thixo_lambda_init;
    for (scalar d in divl)
      d[] = 0.;
  }

  if (!layer) {
    layer = qmalloc (nl, double);
    for (int l = 0; l < nl; l++)
      layer[l] = 1./nl;
  }

  advance = advance_saint_venant_thixo;
  update = update_saint_venant_thixo;

#if TREE
  for (scalar s in {h, zb, eta, u, lambda}) {
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
    s.dirty = true;
  }
  eta.refine = refine_eta;
  eta.restriction = restriction_eta;
  eta.depends = list_copy ({zb, h});
  eta.dirty = true;
#endif
}

event init (i = 0)
{
  foreach()
    eta[] = zb[] + h[];
}

event cleanup (i = end, last)
{
  free (evolving), evolving = NULL;
  free (layer), layer = NULL;
  free (ul), ul = NULL;
  free (lambdal), lambdal = NULL;
  free (divl), divl = NULL;
}
