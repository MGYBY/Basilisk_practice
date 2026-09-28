/** Optional full-NS eigenmode initial condition, disabled unless explicitly
    compiled with -DCASE_LINEAR_EIGENMODE=1. Solver time advancement unchanged.
    Fields are Eulerian modes on flat phases; mapping corrections maintain their
    first-order meaning when evaluated within sinusoidally displaced layers. */
#ifndef LINEAR_EIGENMODE_IC_H
#define LINEAR_EIGENMODE_IC_H
#ifndef CASE_LINEAR_EIGENMODE
#define CASE_LINEAR_EIGENMODE 0
#endif
#if CASE_LINEAR_EIGENMODE
#if CASE_FRONT_RUNNER_ENABLED
#error "A periodic eigenmode requires front_runner.enabled=false"
#endif

typedef struct { int n; double (*row)[11]; } LinearModeTable;
static LinearModeTable linear_mode_table[3];
static double linear_mode_k, linear_mode_amplitude;
static double linear_mode_eta_real[2], linear_mode_eta_imag[2];

static void linear_mode_load (void)
{
  FILE * fp = fopen("base_state/eigenmode_metadata.dat", "r");
  if (!fp || fscanf(fp, "%lf %lf %lf %lf %lf %lf", &linear_mode_k,
      &linear_mode_amplitude, &linear_mode_eta_real[0], &linear_mode_eta_imag[0],
      &linear_mode_eta_real[1], &linear_mode_eta_imag[1]) != 6)
    case_fatal("missing or invalid eigenmode_metadata.dat");
  fclose(fp);
  if (!isfinite(linear_mode_k) || !isfinite(linear_mode_amplitude) ||
      !isfinite(linear_mode_eta_real[0]) || !isfinite(linear_mode_eta_imag[0]) ||
      !isfinite(linear_mode_eta_real[1]) || !isfinite(linear_mode_eta_imag[1]))
    case_fatal("non-finite eigenmode metadata");
  if (linear_mode_amplitude*hypot(linear_mode_eta_real[0],linear_mode_eta_imag[0]) >= H1 ||
      linear_mode_amplitude*hypot(linear_mode_eta_real[1]-linear_mode_eta_real[0],linear_mode_eta_imag[1]-linear_mode_eta_imag[0]) >= H2 ||
      linear_mode_amplitude*hypot(linear_mode_eta_real[1],linear_mode_eta_imag[1]) >= LX-LIQUID_DEPTH)
    case_fatal("eigenmode displaced phase thickness is non-positive");
  if (fabs(linear_mode_k-2.*pi*PERTURBATION_MODE_NUMBER*1./LX) > 1e-10 ||
      fabs(linear_mode_k-2.*pi*WAVE_MODE_NUMBER*1./LX) > 1e-10 ||
      linear_mode_amplitude <= 0. || linear_mode_amplitude > .01)
    case_fatal("eigenmode wavenumber/amplitude mismatch");
  const char * names[3] = {"lower", "upper", "air"};
  for (int a=0; a<3; a++) {
    char path[256]; snprintf(path,sizeof(path),"base_state/eigenmode_%s.dat",names[a]);
    fp=fopen(path,"r");
    LinearModeTable * tb=&linear_mode_table[a];
    if (!fp || fscanf(fp,"%d",&tb->n)!=1 || tb->n<2)
      case_fatal("missing/invalid eigenmode phase table");
    tb->row=malloc((size_t)tb->n*sizeof(*tb->row));
    if (!tb->row) case_fatal("eigenmode allocation failed");
    for(int j=0;j<tb->n;j++) {
      for(int c=0;c<11;c++)
        if(fscanf(fp,"%lf",&tb->row[j][c])!=1 || !isfinite(tb->row[j][c]))
          case_fatal("invalid value in eigenmode table");
      if(j && tb->row[j][0]<=tb->row[j-1][0]) case_fatal("eigenmode z not ascending");
    }
    fclose(fp);
    const double lower_expected=a==0?0.:(a==1?H1:LIQUID_DEPTH);
    const double upper_expected=a==0?H1:(a==1?LIQUID_DEPTH:LX);
    if(fabs(tb->row[0][0]-lower_expected)>1e-10 ||
       fabs(tb->row[tb->n-1][0]-upper_expected)>1e-10)
      case_fatal("eigenmode phase table does not cover expected domain");
  }
  fprintf(ferr,"# eigenmode initialization: A_surface=%g k_star=%g etaI=(%g,%g) etaS=(%g,%g); Eulerian mapping correction enabled\n",
    linear_mode_amplitude,linear_mode_k,linear_mode_eta_real[0],linear_mode_eta_imag[0],
    linear_mode_eta_real[1],linear_mode_eta_imag[1]);
}

static inline double linear_mode_interface (int a, double xx)
{
  return (a==0?H1:LIQUID_DEPTH) + linear_mode_amplitude*(
    linear_mode_eta_real[a]*cos(linear_mode_k*xx)-linear_mode_eta_imag[a]*sin(linear_mode_k*xx));
}

static void linear_mode_fields (double xx, double yy, double * ux,
                                double * uy, double * pressure)
{
  const double hI=linear_mode_interface(0,xx), hS=linear_mode_interface(1,xx);
  const int phase=yy<hI?0:(yy<hS?1:2);
  const double lower0=phase==0?0.:(phase==1?H1:LIQUID_DEPTH);
  const double upper0=phase==0?H1:(phase==1?LIQUID_DEPTH:LX);
  const double lower=phase==0?0.:(phase==1?hI:hS);
  const double upper=phase==0?hI:(phase==1?hS:LX);
  const double z=lower0+(yy-lower)*(upper0-lower0)/(upper-lower);
  LinearModeTable * tb=&linear_mode_table[phase];
  int lo=0, hi=tb->n-1;
  while(hi-lo>1) {int mid=(lo+hi)/2; if(tb->row[mid][0]>z)hi=mid;else lo=mid;}
  const double a=(z-tb->row[lo][0])/(tb->row[hi][0]-tb->row[lo][0]);
  double val[11];
  for(int c=0;c<11;c++) val[c]=(1.-a)*tb->row[lo][c]+a*tb->row[hi][c];
  const double co=linear_mode_amplitude*cos(linear_mode_k*xx);
  const double si=linear_mode_amplitude*sin(linear_mode_k*xx);
  *ux=val[1]+(yy-z)*val[2]+co*val[5]-si*val[6];
  *uy=co*val[7]-si*val[8];
  *pressure=val[3]+(yy-z)*val[4]+co*val[9]-si*val[10];
}

static void linear_mode_free (void)
{
  for(int a=0;a<3;a++) {free(linear_mode_table[a].row);linear_mode_table[a].row=NULL;}
}
#endif
#endif
