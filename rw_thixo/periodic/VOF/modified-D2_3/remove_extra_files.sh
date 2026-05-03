#!/bin/bash

# keep_20s_outputs.sh
# Keep only 20-second-interval files:
#   out-0.gfs, out-20.gfs, out-40.gfs, ...
#   interface-all-0.dat, interface-all-20.dat, interface-all-40.dat, ...
# Remove the rest.

set -euo pipefail

# First run with DRY_RUN=1 to check what would be removed.
# Change to DRY_RUN=0 when you are sure.
DRY_RUN=0

INTERVAL=20

shopt -s nullglob

process_files () {
    local prefix="$1"
    local suffix="$2"

    for file in ${prefix}*${suffix}; do
        base=$(basename "$file")

        # Extract number between prefix and suffix.
        num=${base#"$prefix"}
        num=${num%"$suffix"}

        # Skip unexpected filenames.
        if ! [[ "$num" =~ ^[0-9]+$ ]]; then
            echo "Skipping unexpected file: $file"
            continue
        fi

        # Keep files whose number is divisible by INTERVAL.
        if (( num % INTERVAL != 0 )); then
            if (( DRY_RUN )); then
                echo "Would remove: $file"
            else
                rm -- "$file"
                echo "Removed: $file"
            fi
        fi
    done
}

process_files "out-" ".gfs"
process_files "interface-all-" ".dat"
