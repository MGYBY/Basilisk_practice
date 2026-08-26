#!/usr/bin/env bash
set -euo pipefail

CONFIG=${CONFIG:-steady_uniform_case_parameters.ini}
BASILISK=${BASILISK:-$HOME/basilisk/src}
QCC=${QCC:-$BASILISK/qcc}

mkdir -p verification

python3 generate_steady_uniform_case.py \
  --config "$CONFIG" \
  --output base_state \
  > verification/generation.log 2>&1

rm -rf profiles history analysis
mkdir -p profiles history analysis verification

# -------------------------------------------------------------------
# OpenMP settings
# -------------------------------------------------------------------
export OMP_NUM_THREADS=6
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores

# For 32 threads, spread often gives better memory bandwidth on large AMD EPYC nodes.
# If performance is worse than before, test OMP_PROC_BIND=close.
export OMP_PROC_BIND=close

# Reduce possible glibc malloc memory overhead for many OpenMP threads.
export MALLOC_ARENA_MAX=4

"$QCC" -std=c99 -Wall -O2 -Wdimensions -grid=quadtree \
  steady_uniform_relaxation.c \
  -o steady_uniform_relaxation -lm \
  > verification/compile.stdout \
  2> verification/compile.stderr

./steady_uniform_relaxation

# python3 postprocess_comparison.py --root . \
#   > verification/postprocess.log 2>&1
#
# python3 verify_steady_uniform_case.py \
#   > verification/automated_checks.txt 2>&1
#
# cat analysis/comparison_report.txt
# cat verification/automated_checks.txt
