/** Runtime diagnostics only: no physical-field correction or extra projection.
 * Include after the case macros and global fields in the main source.
 */
#ifndef FRONT_RUNNER_MPI_CASE_DIAGNOSTICS_H
#define FRONT_RUNNER_MPI_CASE_DIAGNOSTICS_H

static void case_fatal (const char * message)
{
  fprintf (stderr, "FATAL [rank %d]: %s\n", pid(), message);
  fflush (stderr);
#if _MPI
  MPI_Abort (MPI_COMM_WORLD, 2);
#endif
  exit (2);
}

static void case_parallel_banner (void)
{
  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif
  fprintf (ferr, "# parallel: MPI=%d ranks=%d OpenMP_max_threads=%d "
                 "auditSteps=%d stopAfterSteps=%d fieldOutputs=%d freezeAMR=%d "
                 "initialAudit=%d volumeAudit=%d diagnosticsEvery=%d "
                 "version=v2p7 embedded=0 viscosity=standard-symmetric-stress\n",
#if _MPI
           1,
#else
           0,
#endif
           npe(), threads, CASE_AUDIT_STEPS, CASE_STOP_AFTER_STEPS,
           !CASE_DISABLE_FIELD_OUTPUTS, CASE_FREEZE_AMR_AFTER_INIT,
           CASE_INITIAL_AUDIT, STRATIFIED_VOLUME_AUDIT, CASE_DIAGNOSTICS_EVERY);
  const char * expected = getenv ("EXPECTED_MPI_RANKS");
  if (expected && *expected) {
    char * end = NULL;
    const long nr = strtol (expected, &end, 10);
    if (!end || *end || nr < 1 || nr != npe())
      case_fatal ("actual MPI size differs from EXPECTED_MPI_RANKS");
  }
  fflush (ferr);
}

static long case_global_leaf_count (void)
{
  long cells = 0;
  foreach (reduction(+:cells))
    cells++;
  return cells;
}

static void case_check_multigrid (const char * stage, mgstats mg, double tolerance)
{
#if CASE_ABORT_ON_MG_FAILURE
  /* i==0 means this solver has not yet been called, not a failure. */
  if (mg.i > 0 && (!isfinite(mg.resa) || mg.resa > tolerance)) {
    fprintf (ferr, "# FIRST FAILED SOLVE: %s i=%d t=%g dt=%g "
                   "cycles=%d resBefore=%g resAfter=%g tolerance=%g nrelax=%d\n",
             stage, iter, t, dt, mg.i, mg.resb, mg.resa, tolerance, mg.nrelax);
    fflush (ferr);
    case_fatal ("multigrid solve failed; stopping before a downstream blow-up");
  }
#endif
}

#if CASE_AUDIT_STEPS > 0
static void case_check_state (const char * stage, bool check_coefficients)
{
  long bad = 0;
  double umax = 0., rhomin = HUGE, rhomax = -HUGE;
  foreach (reduction(+:bad) reduction(max:umax)
           reduction(min:rhomin) reduction(max:rhomax)) {
    if (!isfinite(cm[]) ||
        !isfinite(f_lower[]) || !isfinite(f[]))
      bad++;
    if (cm[] > 0.) {
      if (!isfinite(u.x[]) || !isfinite(u.y[]) || !isfinite(p[]))
        bad++;
      else
        umax = max(umax, hypot(u.x[], u.y[]));
      if (check_coefficients) {
        if (!isfinite(rhov[]) || rhov[] <= 0.)
          bad++;
        else {
          const double r = rhov[]/cm[];
          rhomin = min(rhomin, r);
          rhomax = max(rhomax, r);
        }
      }
    }
  }
  if (bad)
    case_fatal ("non-finite field or non-positive fluid density in state audit");
  if (check_coefficients) {
    face vector muv = mu;
    double mumin = HUGE, mumax = -HUGE;
    double amin = HUGE, amax = -HUGE;
    long badfaces = 0;
    /* Coefficients are synchronized by the pre-viscous event in all modes. */
    foreach_face (reduction(+:badfaces) reduction(min:mumin) reduction(max:mumax)
                  reduction(min:amin) reduction(max:amax)) {
      if (!isfinite(fm.x[]) || !isfinite(muv.x[]) || !isfinite(alphav.x[]))
        badfaces++;
      else if (fm.x[] > 0.) {
        if (muv.x[] <= 0. || alphav.x[] <= 0.)
          badfaces++;
        else {
          const double viscosity = muv.x[]/fm.x[];
          const double inverse_density = alphav.x[]/fm.x[];
          mumin = min(mumin, viscosity); mumax = max(mumax, viscosity);
          amin = min(amin, inverse_density); amax = max(amax, inverse_density);
        }
      }
    }
    if (badfaces)
      case_fatal ("non-finite or non-positive open-face coefficient in audit");
    fprintf (ferr, "# state-audit stage=%s i=%d t=%g Umax=%g "
                   "rho=[%g,%g] etaFace=[%g,%g] invRhoFace=[%g,%g]\n",
             stage, iter, t, umax, rhomin, rhomax, mumin, mumax, amin, amax);
  }
  else
    fprintf (ferr, "# state-audit stage=%s i=%d t=%g Umax=%g\n",
             stage, iter, t, umax);
  fflush (ferr);
}

#endif // CASE_AUDIT_STEPS > 0

static void case_write_summary (const char * reason)
{
  long cells = 0;
  double vl = 0., vliquid = 0., momentum = 0., energy = 0., umax = 0.;
  foreach (reduction(+:cells) reduction(+:vl) reduction(+:vliquid)
           reduction(+:momentum) reduction(+:energy) reduction(max:umax)) {
    cells++;
    if (cm[] > 0.) {
      const double vol = dv();
      const double r = stratified_density_cell(f_lower[], f[]);
      const double speed2 = sq(u.x[]) + sq(u.y[]);
      vl += f_lower[]*vol;
      vliquid += f[]*vol;
      momentum += r*u.x[]*vol;
      energy += 0.5*r*speed2*vol;
      umax = max(umax, sqrt(speed2));
    }
  }
  if (pid() == 0) {
    FILE * fp = fopen ("run_summary.tsv", "w");
    if (!fp) case_fatal ("cannot open run_summary.tsv");
    fprintf (fp, "reason\tranks\ti\tt\tleaves\tVlower\tVliquid\tPx\tKE\tUmax\n");
    fprintf (fp, "%s\t%d\t%d\t%.17g\t%ld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
             reason, npe(), iter, t, cells, vl, vliquid, momentum, energy, umax);
    if (fclose(fp) != 0) case_fatal ("cannot close run_summary.tsv");
  }
}
#endif
