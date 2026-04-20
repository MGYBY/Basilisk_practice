#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Example run script for the periodic thixotropic VOF roll-wave case.
#
# What this script does:
#   1. build the steady base-state profile used for initialisation,
#   2. compile the Basilisk case with user-adjustable dimensionless parameters,
#   3. run the simulation in the case directory so that all outputs are written there.
#
# Usage pattern:
#   export BASILISK_SRC=/path/to/basilisk/src
#   ./runCase.sh
#
# You can override any of the parameters below from the shell, e.g.
#   THIXO_GAMMA=8 THIXO_KAPPA=1e-6 ./runCase.sh
# -----------------------------------------------------------------------------

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
BASILISK_SRC="${BASILISK_SRC:-}"

THIXO_T="${THIXO_T:-1.0}"
THIXO_GAMMA="${THIXO_GAMMA:-5.0}"
THIXO_KAPPA="${THIXO_KAPPA:-2e-6}"
THIXO_A="${THIXO_A:-0.2}"
THIXO_FRV="${THIXO_FRV:-0.75}"
THIXO_SO="${THIXO_SO:-0.06}"
THIXO_RHOR="${THIXO_RHOR:-0.01}"
THIXO_MUR="${THIXO_MUR:-0.02}"

BASE_STATE_SCRIPT="$CASE_DIR/thixo_base_state_explicitFrV_fixed.py"
CASE_SOURCE="$CASE_DIR/thixo_periodic_rollwave.c"
CASE_EXE="$CASE_DIR/thixo_periodic_rollwave"
PROFILE_TXT="$CASE_DIR/base_state_profile.txt"

if [[ ! -f "$BASE_STATE_SCRIPT" ]]; then
  echo "Missing base-state generator: $BASE_STATE_SCRIPT"
  exit 1
fi

if [[ ! -f "$CASE_SOURCE" ]]; then
  echo "Missing case source: $CASE_SOURCE"
  exit 1
fi

if [[ -z "$BASILISK_SRC" ]]; then
  echo "Please set BASILISK_SRC to the Basilisk src directory containing qcc."
  echo "Example: export BASILISK_SRC=/path/to/basilisk/src"
  exit 1
fi

if [[ ! -x "$BASILISK_SRC/qcc" ]]; then
  echo "Cannot find an executable qcc at: $BASILISK_SRC/qcc"
  exit 1
fi

cd "$CASE_DIR"

echo "[1/3] Building steady base-state profile..."
python3 "$BASE_STATE_SCRIPT" \
  --kappa "$THIXO_KAPPA" \
  --a "$THIXO_A" \
  --Gamma "$THIXO_GAMMA" \
  --FrV "$THIXO_FRV" \
  --output "$PROFILE_TXT"

echo "[2/3] Compiling Basilisk case..."
"$BASILISK_SRC/qcc" -std=c99 -Wall -O2 -grid=quadtree \
  -DTHIXO_T="$THIXO_T" \
  -DTHIXO_GAMMA="$THIXO_GAMMA" \
  -DTHIXO_KAPPA="$THIXO_KAPPA" \
  -DTHIXO_A="$THIXO_A" \
  -DTHIXO_FRV="$THIXO_FRV" \
  -DTHIXO_SO="$THIXO_SO" \
  -DTHIXO_RHOR="$THIXO_RHOR" \
  -DTHIXO_MUR="$THIXO_MUR" \
  "$CASE_SOURCE" -o "$CASE_EXE" -lm

echo "[3/3] Running simulation..."
"$CASE_EXE"
