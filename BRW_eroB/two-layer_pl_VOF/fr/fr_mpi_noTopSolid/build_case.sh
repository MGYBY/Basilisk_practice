#!/usr/bin/env bash
# Shared build fix: top-solid v2p6a / no-top-solid v2p7a.
# POSIX GFS helper is compiled natively; qcc receives only its object file.
# Retain the earlier workaround: QCC_DIMENSIONS=off by default.
# 'warn' still runs qcc's dimensional interpreter and is NOT a bypass.
# CPU_TARGET=native assumes compilation on the same CPU type as execution.
# EXTRA_QCC_FLAGS retains the previous whitespace-separated flag interface.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
mode=${1:-mpi}
exe=${2:-rollwave}
: "${PYTHON:=python3}"
: "${QCC:=qcc}"
: "${CASE_CONFIG:=case_parameters.ini}"
if [[ -n ${BASILISK:-} ]]; then export PATH="$BASILISK:$PATH"; fi
command -v "$QCC" >/dev/null || { echo "qcc not found; set BASILISK or QCC." >&2; exit 127; }
flags=(-Wall -D_GNU_SOURCE -D_DEFAULT_SOURCE)
: "${BUILD_PROFILE:=release}"
case "$BUILD_PROFILE" in
  release) flags+=(-O3 -DNDEBUG) ;;
  debug)   flags+=(-O2 -g) ;;
  *) echo 'BUILD_PROFILE must be release or debug.' >&2; exit 2 ;;
esac
: "${CPU_TARGET:=native}"
case "$CPU_TARGET" in
  native)   flags+=(-march=native) ;;
  portable) ;; # for execution on a different CPU type
  *) echo 'CPU_TARGET must be native or portable.' >&2; exit 2 ;;
esac
# Copy only native compiler flags, before appending qcc-specific switches.
# Do not globally disable libc fortification or strip system attributes.
native_flags=("${flags[@]}")
: "${QCC_DIMENSIONS:=off}"
case "$QCC_DIMENSIONS" in
  off)    flags+=(-disable-dimensions) ;;
  warn)   flags+=(-Wdimensions) ;;
  strict) ;;  # qcc's default dimensional-consistency check
  *) echo 'QCC_DIMENSIONS must be off, warn, or strict.' >&2; exit 2 ;;
esac
# A plain output filename avoids qcc's internal unquoted command/path handling.
[[ $exe =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ &&
   ! $exe =~ \.(c|h|py|sh|ini|md|txt|log)$ ]] || {
  echo 'Use a plain executable filename, e.g. rollwave (not a source filename).' >&2
  exit 2
}
case "$mode" in
  mpi)    : "${MPICC:=mpicc}"; export CC99="$MPICC -std=gnu99"; flags+=(-D_MPI=1); compiler=$MPICC ;;
  omp)    : "${CC:=gcc}"; export CC99="$CC -std=gnu99"; flags+=(-fopenmp); compiler=$CC ;;
  serial) : "${CC:=gcc}"; export CC99="$CC -std=gnu99"; compiler=$CC ;;
  *) echo "Usage: $0 {mpi|omp|serial} [executable]" >&2; exit 2 ;;
esac
command -v "$compiler" >/dev/null || { echo "$compiler not found" >&2; exit 127; }
extra=()
if [[ -n ${EXTRA_QCC_FLAGS:-} ]]; then read -r -a extra <<< "$EXTRA_QCC_FLAGS"; fi
for required in gfs_tail_truncate.c gfs_tail_truncate.h; do
  [[ -f $required ]] || {
    echo "Missing $required: copy the complete common build update, not just build_case.sh." >&2
    exit 2
  }
done
build_revision=v2p6a
branch=top-solid
if ! grep -Eq '^[[:space:]]*#[[:space:]]*include[[:space:]]+"embed.h"' roll_wave_amr_dimensionless.c; then
  build_revision=v2p7a
  branch=no-top-solid
fi
staged="${exe}.building"
helper_object="gfs_tail_truncate.${mode}.building.o"
"$PYTHON" generate_dimensionless_base_state.py --config "$CASE_CONFIG" --output base_state
{
  date -Is
  printf 'build_revision=%s\nbranch=%s\nmode=%s\nCC99=%s\nBASILISK=%s\nQCC_DIMENSIONS=%s\n' \
    "$build_revision" "$branch" "$mode" "$CC99" "${BASILISK:-unset}" "$QCC_DIMENSIONS"
  printf 'BUILD_PROFILE=%s\nCPU_TARGET=%s\n' "$BUILD_PROFILE" "$CPU_TARGET"
  printf 'stack_limit_kib=%s\n' "$(ulimit -s)"
  printf 'native_io_command:'
  printf ' %q' "$compiler" -std=gnu99 "${native_flags[@]}" -c gfs_tail_truncate.c -o "$helper_object"
  printf '\ncommand:'
  printf ' %q' "$QCC" "${flags[@]}" "${extra[@]}" roll_wave_amr_dimensionless.c -o "$staged" "$helper_object" -lm -grid=quadtree
  printf '\n'
  "$compiler" --version
  command -v "$QCC"
  sha256sum "$(command -v "$QCC")" ./*.c ./*.h "$CASE_CONFIG" base_state/generated_case.h
  if [[ -n ${BASILISK:-} ]]; then
    for h in navier-stokes/centered.h navier-stokes/conserving.h viscosity-embed.h \
             viscosity.h poisson.h embed.h vof.h output.h grid/tree-common.h \
             grid/tree-mpi.h grid/cartesian-common.h grid/config.h curvature.h \
             qcc.c ast/translate.c ast/interpreter/dimension.c config; do
      [[ ! -f $BASILISK/$h ]] || sha256sum "$BASILISK/$h"
    done
  fi
  if [[ $mode == mpi ]]; then mpirun --version; fi
} > "build_manifest.${mode}.txt" 2>&1
# A crashed qcc can leave an executable even though the overall build failed.
# Publish only after a zero exit status; never run such an unaccepted binary.
trap 'rm -f -- "$staged" "$helper_object"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
rm -f -- "$exe" "$staged" "$helper_object"
echo "Build fix $build_revision ($branch): compiling POSIX GFS helper with $compiler, not qcc."
set +e
"$compiler" -std=gnu99 "${native_flags[@]}" -c gfs_tail_truncate.c \
  -o "$helper_object" 2>&1 | tee "build.native-io.${mode}.log"
io_status=("${PIPESTATUS[@]}")
set -e
if (( io_status[0] != 0 )); then
  echo "BUILD FAILED: native GFS helper exit status=${io_status[0]}; qcc and simulation were not started." >&2
  exit "${io_status[0]}"
fi
if (( io_status[1] != 0 )); then
  echo 'BUILD FAILED: native-helper log could not be written.' >&2
  exit "${io_status[1]}"
fi
[[ -s $helper_object ]] || {
  echo 'BUILD FAILED: native compiler returned success without the GFS helper object.' >&2
  exit 1
}
sha256sum "$helper_object" >> "build_manifest.${mode}.txt"
echo "Building ($mode), profile=$BUILD_PROFILE, CPU=$CPU_TARGET, qcc dimensional checker=$QCC_DIMENSIONS; simulation not started."
set +e
"$QCC" "${flags[@]}" "${extra[@]}" roll_wave_amr_dimensionless.c \
  -o "$staged" "$helper_object" -lm -grid=quadtree 2>&1 | tee "build.${mode}.log"
build_status=("${PIPESTATUS[@]}")
set -e
if (( build_status[0] != 0 )); then
  echo "BUILD FAILED: qcc exit status=${build_status[0]}; simulation was not launched." >&2
  echo "See build.${mode}.log and build_manifest.${mode}.txt." >&2
  if (( build_status[0] == 139 )); then
    echo 'qcc received/reported SIGSEGV. This is a build failure, not a flow divergence.' >&2
    echo 'Preserve the build logs. Run bash diagnose_build.sh for stage-separated evidence.' >&2
  fi
  exit "${build_status[0]}"
fi
if (( build_status[1] != 0 )); then
  echo 'BUILD FAILED: compiler log could not be written.' >&2
  exit "${build_status[1]}"
fi
# A parser warning means the translated program cannot be trusted, even if
# this qcc revision returns zero. Do not publish such a binary.
if grep -Fq 'Basilisk C parse error' "build.${mode}.log"; then
  echo 'BUILD FAILED: Basilisk C parser errors remain; preserve build.mpi.log and run bash diagnose_build.sh.' >&2
  exit 65
fi
[[ -s $staged && -x $staged ]] || {
  echo 'BUILD FAILED: qcc returned success but no executable was produced.' >&2
  exit 1
}
mv -f -- "$staged" "$exe"
echo "Built $exe ($mode); see build_manifest.${mode}.txt"
