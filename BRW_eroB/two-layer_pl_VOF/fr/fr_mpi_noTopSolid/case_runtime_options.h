/** Production monitoring controls. All switches can be overridden with -D.
 * These change diagnostics only: not the PDE, tolerances, AMR or MPI halos.
 * Re-enable the detailed diagnostic mode with -DCASE_DEBUG_AUDITS=1.
 */
#ifndef FRONT_RUNNER_CASE_RUNTIME_OPTIONS_H
#define FRONT_RUNNER_CASE_RUNTIME_OPTIONS_H

#ifndef CASE_DEBUG_AUDITS
# define CASE_DEBUG_AUDITS 0
#endif
#ifndef CASE_AUDIT_STEPS
# if CASE_DEBUG_AUDITS
#  define CASE_AUDIT_STEPS 10
# else
#  define CASE_AUDIT_STEPS 0
# endif
#endif
#ifndef CASE_INITIAL_AUDIT
# define CASE_INITIAL_AUDIT CASE_DEBUG_AUDITS
#endif
#ifndef CASE_VERBOSE_INIT
# define CASE_VERBOSE_INIT CASE_DEBUG_AUDITS
#endif
#ifndef STRATIFIED_VOLUME_AUDIT
# define STRATIFIED_VOLUME_AUDIT CASE_DEBUG_AUDITS
#endif
#ifndef CASE_WRITE_BOUNDARY_GUARD
# define CASE_WRITE_BOUNDARY_GUARD CASE_DEBUG_AUDITS
#endif
#ifndef CASE_DIAGNOSTICS_EVERY
# if CASE_DEBUG_AUDITS
#  define CASE_DIAGNOSTICS_EVERY 100
# else
#  define CASE_DIAGNOSTICS_EVERY 1000
# endif
#endif
#ifndef CASE_STOP_AFTER_STEPS
# define CASE_STOP_AFTER_STEPS 0
#endif
#ifndef CASE_DISABLE_FIELD_OUTPUTS
# define CASE_DISABLE_FIELD_OUTPUTS 0
#endif
#ifndef CASE_FREEZE_AMR_AFTER_INIT
# define CASE_FREEZE_AMR_AFTER_INIT 0
#endif
#ifndef CASE_ABORT_ON_MG_FAILURE
# define CASE_ABORT_ON_MG_FAILURE 1
#endif

#if CASE_AUDIT_STEPS < 0 || CASE_DIAGNOSTICS_EVERY < 1
# error "CASE_AUDIT_STEPS must be non-negative; CASE_DIAGNOSTICS_EVERY must be positive"
#endif
#endif
