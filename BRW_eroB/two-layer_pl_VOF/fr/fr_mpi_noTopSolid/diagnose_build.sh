#!/usr/bin/env bash
# Optional compiler-stage isolation. Does NOT run a simulation.
# Uses a temporary copy; does not edit production inputs, outputs or executable.
# This diagnostic is deliberately fixed to pure MPI and default production
# switches. It does not accept EXTRA_QCC_FLAGS (use build_case.sh for those).
set -euo pipefail
root=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$root"
: "${QCC:=qcc}"
: "${MPICC:=mpicc}"
: "${PYTHON:=python3}"
: "${CASE_CONFIG:=case_parameters.ini}"
: "${DIAGNOSE_DIMENSIONS:=0}"
[[ $DIAGNOSE_DIMENSIONS == 0 || $DIAGNOSE_DIMENSIONS == 1 ]] || {
  echo 'DIAGNOSE_DIMENSIONS must be 0 or 1.' >&2; exit 2;
}
if [[ -n ${BASILISK:-} ]]; then export PATH="$BASILISK:$PATH"; fi
for command in "$QCC" "$MPICC" "$PYTHON"; do
  command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 127; }
done
QCC=$(command -v "$QCC")
MPICC=$(command -v "$MPICC")
PYTHON=$(command -v "$PYTHON")
# Short working paths also reduce dependence on long cluster directory names.
work=$(mktemp -d /tmp/frqcc.XXXXXX)
result=$(mktemp -d "$root/compiler_diagnosis_XXXXXX")
cleanup () { cp -a "$work"/. "$result"/; rm -rf -- "$work"; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
cp -- "$root"/*.c "$root"/*.h "$work"/
cp -- "$root"/generate_dimensionless_base_state.py "$root"/front_runner_initialization.py "$work"/
cp -- "$CASE_CONFIG" "$work/case_parameters.ini"
cd "$work"
export CC99="$MPICC -std=gnu99"
"$PYTHON" generate_dimensionless_base_state.py --config case_parameters.ini --output base_state > generation.log 2>&1
{
  date -Is
  printf 'QCC=%s\nCC99=%s\nBASILISK=%s\nstack_limit_kib=%s\n' \
    "$QCC" "$CC99" "${BASILISK:-unset}" "$(ulimit -s)"
  "$MPICC" --version
} > environment.txt 2>&1
common=(-Wall -O2 -g -D_GNU_SOURCE -D_DEFAULT_SOURCE -D_MPI=1)
echo "Compiler evidence will be saved in: $result"
echo 'Stage 1: qcc -> generated C; dimensional interpreter DISABLED.'
set +e
"$QCC" -source -disable-dimensions "${common[@]}" \
  roll_wave_amr_dimensionless.c -grid=quadtree > translation.log 2>&1
translation_status=$?
set -e
printf 'translation_exit_status=%s\n' "$translation_status" | tee status.txt
if (( translation_status != 0 )); then
  echo 'Translation failed. See translation.log; no MPI solver was run.' >&2
  exit "$translation_status"
fi
[[ -s _roll_wave_amr_dimensionless.c ]] || {
  echo 'qcc returned success without generated C.' >&2; exit 1;
}
if grep -Fq 'Basilisk C parse error' translation.log; then
  echo 'Translation returned zero but parser errors remain; no native compile/run.' >&2
  exit 65
fi
echo 'Stage 2: mpicc compiles/links the generated C and POSIX helper; qcc is NOT involved.'
includes=(-I.)
if [[ -n ${BASILISK:-} ]]; then includes+=("-I$BASILISK"); fi
set +e
"$MPICC" -std=gnu99 "${common[@]}" "${includes[@]}" \
  _roll_wave_amr_dimensionless.c gfs_tail_truncate.c -o diagnostic_solver -lm > native_compile.log 2>&1
native_status=$?
set -e
printf 'native_compile_exit_status=%s\n' "$native_status" | tee -a status.txt
if (( native_status != 0 )); then
  echo 'Native C compile/link failed. See native_compile.log.' >&2
  exit "$native_status"
fi
[[ -s diagnostic_solver && -x diagnostic_solver ]] || {
  echo 'Native compiler returned success without an executable.' >&2; exit 1;
}
# Keep relevant generated MPI code as evidence without modifying it.
grep -n -E 'mpi_sum_reduce_(init|array)|MPI_INT|int sh' \
  _roll_wave_amr_dimensionless.c > reduction_excerpts.txt || true
if [[ $DIAGNOSE_DIMENSIONS == 1 ]]; then
  # This optional A/B test may reproduce the qcc crash. Never use its output
  # as the production executable. -debug retains qcc's intermediate files.
  echo 'Stage 3: opt-in qcc dimensional-check comparison (may crash or hang).'
  set +e
  "$QCC" -source -debug -Wdimensions "${common[@]}" \
    roll_wave_amr_dimensionless.c -grid=quadtree > dimensions_enabled.log 2>&1
  dimensions_status=$?
  set -e
  printf 'dimensions_enabled_exit_status=%s\n' "$dimensions_status" | tee -a status.txt
fi
echo 'Compiler-stage checks completed; diagnostic_solver was NOT executed.'
