#!/usr/bin/env bash
# Remove this case's generated run outputs. Keep inputs and base profiles unless
# --base-state is explicitly passed. No broad user-case filename glob is used.
set -euo pipefail
case_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$case_dir"
if [[ $# -gt 1 || ($# -eq 1 && "$1" != "--base-state") ]]; then
  echo "Usage: ./Allclean.sh [--base-state]" >&2
  exit 2
fi
rm -rf gfs dumps interfaces slices wave_amplitude_plots
rm -f rollwave_dimensionless roll_wave_amr roll_wave_amr_smoke \
  roll_wave_amr_steady rollwave a.out _roll_wave_amr.c
rm -f out-*.gfs log log.* perfs run_configuration.txt
rm -f wave_amplitude.tsv initial_condition_audit.tsv boundary_guard.tsv
rm -rf .qcc*
if [[ "${1:-}" == "--base-state" ]]; then rm -rf base_state; fi
