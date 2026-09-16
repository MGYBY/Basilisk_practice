/** Localized constant-Froude initialization only: no time-stepping changes.
 * Coordinates x,y are x/H_l,z/H_l. The prefix integral is the EXACT integral
 * of the piecewise-linear velocity.dat interpolant, not another interpolant.
 * Include profile1d.h before this header. Pure C99; also tested with gcc.
 */
#ifndef FRONT_RUNNER_IC_H
#define FRONT_RUNNER_IC_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
  const Profile1D *base;
  double *prefix;
  double amplitude, wavelength, center, depth, ceiling;
  double q0, us, qt;
} FrontRunnerIC;

static double fr_base_integral(const FrontRunnerIC *ic, double yy)
{
  const Profile1D *p = ic->base;
  if (yy <= p->x[0]) return (yy-p->x[0])*p->value[0];
  if (yy >= p->x[p->n-1])
    return ic->prefix[p->n-1] + (yy-p->x[p->n-1])*p->value[p->n-1];
  int lo=0, hi=p->n-1;
  while (hi-lo>1) {
    const int mid=(lo+hi)/2;
    if (p->x[mid]<=yy) lo=mid; else hi=mid;
  }
  const double dz=yy-p->x[lo];
  const double slope=(p->value[hi]-p->value[lo])/(p->x[hi]-p->x[lo]);
  return ic->prefix[lo]+p->value[lo]*dz+0.5*slope*dz*dz;
}

static void fr_init(FrontRunnerIC *ic, const Profile1D *p,
                    double amplitude, double wavelength, double center,
                    double depth, double ceiling)
{
  if (!p || p->n<2 || fabs(p->x[0])>1.e-13 || depth<=0. ||
      ceiling<=depth*(1.+amplitude) || wavelength<=0. ||
      amplitude<0. || amplitude>=1.) {
    fprintf(stderr,"Invalid localized front-runner profile/geometry.\n");
    exit(1);
  }
  ic->base=p; ic->amplitude=amplitude; ic->wavelength=wavelength;
  ic->center=center; ic->depth=depth; ic->ceiling=ceiling;
  ic->prefix=(double *)calloc((size_t)p->n,sizeof(double));
  if (!ic->prefix) { perror("front-runner prefix allocation"); exit(1); }
  for (int j=1;j<p->n;j++)
    ic->prefix[j]=ic->prefix[j-1]+
      0.5*(p->value[j-1]+p->value[j])*(p->x[j]-p->x[j-1]);
  ic->q0=fr_base_integral(ic,depth);
  ic->us=profile1d_eval(p,depth);
  ic->qt=ic->q0+ic->us*(ceiling-depth);
  for (int j=0;j<p->n;j++)
    if (p->x[j]>=depth && fabs(p->value[j]-ic->us)>1.e-10*fmax(1.,fabs(ic->us))) {
      fprintf(stderr,"Front-runner extension requires the uploaded uniform, unforced-air base profile.\n");
      exit(1);
    }
}

static inline double fr_bump(const FrontRunnerIC *ic, double xx)
{
  const double xi=(xx-ic->center)/ic->wavelength;
  if (xi<=-0.25 || xi>=0.25) return 0.;
  return cos(2.*acos(-1.)*xi);
}

static inline void fr_scale(const FrontRunnerIC *ic,double xx,
                            double *s,double *sx)
{
  const double xi=(xx-ic->center)/ic->wavelength;
  *s=1.; *sx=0.;
  if (xi<=-0.25 || xi>=0.25) return;
  const double k=2.*acos(-1.)/ic->wavelength;
  *s+=ic->amplitude*cos(2.*acos(-1.)*xi);
  *sx=-ic->amplitude*k*sin(2.*acos(-1.)*xi);
}

/** u=psi_y, w=-psi_x. The optional psi output permits an independent curl
 * check and a conservative face-flux initialization in future developments.
 * This release does NOT add a projection or alter native face-velocity setup.
 */
static void fr_evaluate(const FrontRunnerIC *ic,double xx,double yy,
                        double *ux,double *wy,double *psi)
{
  double s,sx;
  fr_scale(ic,xx,&s,&sx);
  if (yy<=0.) { *ux=0.; *wy=0.; if(psi) *psi=0.; return; }
  if (yy>=ic->ceiling) {
    *ux=ic->us; *wy=0.;
    if(psi) *psi=ic->qt+ic->us*(yy-ic->ceiling);
    return;
  }
  const double rs=sqrt(s), h=s*ic->depth;
  if (yy<=h) {
    const double z0=yy/s;
    const double U=profile1d_eval(ic->base,z0);
    const double P=fr_base_integral(ic,z0);
    *ux=rs*U;
    *wy=-rs*sx*(1.5*P-z0*U);
    if(psi) *psi=s*rs*P;
    return;
  }
  /* Confined air: surface velocity and zero normal velocity at the fixed lid;
     a smooth return-flow correction maintains x-independent TOTAL discharge. */
  const double d=ic->ceiling-h, hx=ic->depth*sx;
  const double r=(yy-h)/d, r2=r*r, r3=r2*r, r4=r3*r;
  const double A=1.-3.*r2+2.*r3, Ai=r-r3+0.5*r4;
  const double C=30.*r2*(1.-r)*(1.-r), Ci=10.*r3-15.*r4+6.*r4*r;
  const double Q=s*rs*ic->q0, Qx=1.5*rs*sx*ic->q0;
  const double du=ic->us*(rs-1.), dux=ic->us*sx/(2.*rs);
  const double B=(ic->qt-Q)/d-ic->us-0.5*du;
  const double Bx=-Qx/d+(ic->qt-Q)*hx/(d*d)-0.5*dux;
  const double F=ic->us*r+du*Ai+B*Ci;
  *ux=ic->us+du*A+B*C;
  *wy=-Qx+hx*F-d*(dux*Ai+Bx*Ci)-hx*(r-1.)*(*ux);
  if(psi) *psi=Q+d*F;
}

static inline void fr_velocity(const FrontRunnerIC *ic,double xx,double yy,
                               double *ux,double *wy)
{
  fr_evaluate(ic,xx,yy,ux,wy,NULL);
}
#endif
