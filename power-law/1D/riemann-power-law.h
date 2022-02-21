/* A collection of Riemann solvers for the Saint-Venant system 
* power-law model
 */

#define epsilon 1e-30

void hllc (double hm, double hp, double um, double up, double Delta,
	   double * fh, double * fq, double * dtmax)
{
  // eigenvalues based on power-law system:
  // lambda1 = u beta + sqrt(alpha h + u^2 (beta^2-beta))
  // lambda2 = u beta
  // lambda3 = u beta - sqrt(alpha h + u^2 (beta^2-beta))
  // double cm = sqrt (G*hm), cp = sqrt (G*hp);
  double cm = sqrt (alphaCoeff*hm+um*um*(betaCoeff*betaCoeff-betaCoeff));
  double cp = sqrt (alphaCoeff*hp+up*up*(betaCoeff*betaCoeff-betaCoeff));
  double ustar = betaCoeff*(um + up)/2. + cm - cp;
  double cstar = (cm + cp)/2. + betaCoeff*(um - up)/4.;
  double SL = hm == 0. ? betaCoeff*up - 2.*cp : min (betaCoeff*um - cm, ustar - cstar);
  double SR = hp == 0. ? betaCoeff*um + 2.*cm : max (betaCoeff*up + cp, ustar + cstar);

  if (0. <= SL) {
    // *fh = um*hm;
    // *fq = hm*(um*um + G*hm/2.);
    *fh = um*hm;
    *fq = hm*(betaCoeff*um*um + alphaCoeff*hm/2.);
  }
  else if (0. >= SR) {
    // *fh = up*hp;
    // *fq = hp*(up*up + G*hp/2.);
    *fh = up*hp;
    *fq = hp*(betaCoeff*up*up + alphaCoeff*hp/2.);
  }
  else {
    // double fhm = um*hm;
    // double fum = hm*(um*um + G*hm/2.);
    double fhm = um*hm;
    double fum = hm*(betaCoeff*um*um + alphaCoeff*hm/2.);
    // double fhp = up*hp;
    // double fup = hp*(up*up + G*hp/2.);
    double fhp = up*hp;
    double fup = hp*(betaCoeff*up*up + alphaCoeff*hp/2.);
    *fh = (SR*fhm - SL*fhp + SL*SR*(hp - hm))/(SR - SL);
    *fq = (SR*fum - SL*fup + SL*SR*(hp*up - hm*um))/(SR - SL);
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > epsilon) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}
