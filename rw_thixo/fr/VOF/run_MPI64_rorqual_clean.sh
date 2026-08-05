#!/bin/bash
#SBATCH --account=def-sushama-ab_cpu
#SBATCH --job-name=FrReal-1p75_a-0.2_T-100_Gamma-8_kappa-1e-4_long_MPI_v3_fixedMPI-1_topExtended
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --time=0-23:55
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



module --force purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module list

export BASILISK=/home/yboyuan/links/projects/rrg-sushama-ab/yboyuan/basilisk/src
export PATH="$BASILISK:$PATH"

# Pure MPI: no hidden OpenMP oversubscription inside each rank.
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MALLOC_ARENA_MAX=2
export SLURM_EXPORT_ENV=ALL

SRC=corrected_thixo_periodic_rollwave_mpi_clean.c
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
