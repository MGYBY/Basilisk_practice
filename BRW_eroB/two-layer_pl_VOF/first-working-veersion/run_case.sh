#!/usr/bin/env bash

set -euo pipefail

# module purge
# module load StdEnv/2023
# module load gcc
# module load python

# Edit this path for the cluster installation.
# export BASILISK=/home/yboyuan/links/projects/rrg-sushama-ab/yboyuan/basilisk/src
# export PATH="$BASILISK:$PATH"

# case_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
# cd "$case_dir"

# ./Allgenerate
python generate_base_state.py

# Override at run time, e.g. OMP_NUM_THREADS=12 ./AllrunOpenMP
export OMP_NUM_THREADS=32
export OMP_DYNAMIC="FALSE"
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

qcc  -std=c99 -fopenmp -Wall -O3 -Wdimensions roll_wave_amr.c -o rollwave -lm -grid=quadtree -L$BASILISK/gl -lglutils -lfb_tiny
# serial
# qcc -O2 -Wall -Wdimensions -o rollwave thixo_periodic_rollwave_mu_capped_noview_restartfixed5.c -lm
./rollwave
