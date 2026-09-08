#!/bin/bash
#SBATCH --account=def-sushama-ab_cpu
#SBATCH --job-name=Frl-0p6_So-0p06_nl-0p4_nu-0p8_nl-0p4_rhoR-0p8_hR-1_REtaI1p54_OMP
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --ntasks-per-node=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=0-23:30
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=boyuan.yu@mail.mcgill.ca
#SBATCH --hint=nomultithread

set -euo pipefail

echo "Job started at: $(date)"
echo "Node list: $SLURM_JOB_NODELIST"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-1}"

if [[ "${SLURM_NTASKS}" -ne 64 ]]; then
  echo "ERROR: this script is intended for exactly 64 MPI ranks."
  exit 2
fi

module --force purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module list
module load python

export BASILISK=/home/yboyuan/links/projects/rrg-sushama-ab/yboyuan/basilisk/src
export PATH="$BASILISK:$PATH"

# Pure MPI: no hidden OpenMP oversubscription inside each rank.
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MALLOC_ARENA_MAX=2
export SLURM_EXPORT_ENV=ALL

python generate_dimensionless_base_state.py --config "./case_parameters.ini" --output "./base_state"

SRC=roll_wave_amr_dimensionless.c
EXE=rollwave

rm -f "$EXE"

echo "Compiling $SRC with Basilisk MPI..."
CC99='mpicc -std=gnu99' qcc \
  -Wall -O3 -march=native -DNDEBUG \
  -D_MPI=1 -D_GNU_SOURCE -D_DEFAULT_SOURCE \
  -Wdimensions \
  "$SRC" -o "$EXE" \
  -lm -grid=quadtree

echo "Compilation finished at: $(date)"
echo "Running with $SLURM_NTASKS MPI ranks..."

mpirun --bind-to core --map-by core \
  ./"$EXE" -parallel > "log.${SLURM_JOB_ID}" 2>&1

echo "Job finished at: $(date)"
