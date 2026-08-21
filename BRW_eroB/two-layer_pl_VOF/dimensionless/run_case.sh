#!/usr/bin/env bash
# Generate, compile, and execute the dimensionless three-phase Basilisk case.
set -euo pipefail

case_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$case_dir"

: "${PYTHON:=python3}"
: "${CASE_CONFIG:=$case_dir/case_parameters.ini}"
: "${RUN_MODE:=openmp}"          # openmp or serial
: "${OMP_NUM_THREADS:=32}"
: "${CC_OPT:=-O3}"

if ! command -v qcc >/dev/null 2>&1; then
  cat >&2 <<'MESSAGE'
ERROR: qcc was not found in PATH.
Source your Basilisk installation first, for example:

  export BASILISK=$HOME/basilisk/src
  export PATH="$PATH:$BASILISK"

Then rerun ./run_case.sh.
MESSAGE
  exit 127
fi

"$PYTHON" generate_dimensionless_base_state.py \
  --config "$CASE_CONFIG" \
  --output "$case_dir/base_state"

common=(
  -std=c99
  -Wall
  "$CC_OPT"
  -Wdimensions
  -grid=quadtree
  roll_wave_amr_dimensionless.c
  -o rollwave_dimensionless
  -lm
)

case "$RUN_MODE" in
  serial)
    qcc "${common[@]}"
    ;;
  openmp)
    export OMP_NUM_THREADS
    export OMP_DYNAMIC="${OMP_DYNAMIC:-FALSE}"
    export OMP_PLACES="${OMP_PLACES:-cores}"
    export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
    qcc -fopenmp "${common[@]}"
    ;;
  *)
    echo "ERROR: RUN_MODE must be 'serial' or 'openmp' (got '$RUN_MODE')." >&2
    exit 2
    ;;
esac

./rollwave_dimensionless 2> log.dimensionless
