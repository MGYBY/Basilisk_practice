#!/bin/bash
#SBATCH --account=def-sushama-ab_cpu
#SBATCH --job-name=FrReal-1_a-0.2_T-100_Gamma-4_kappa-1e-4_lambda-12_v3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=111M           # memory for the entire job across all cores 3900M
#SBATCH --time=0-12:30
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=boyuan.yu@mail.mcgill.ca
#SBATCH --hint=nomultithread

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK}"

# -------------------------------------------------------------------
# Modules
# -------------------------------------------------------------------
module purge
module load StdEnv/2023
module load gcc/12.3
module load gnuplot
module load python

# -------------------------------------------------------------------
# Basilisk environment
# -------------------------------------------------------------------
export BASILISK=/home/yboyuan/links/projects/rrg-sushama-ab/yboyuan/basilisk/src
export PATH="$BASILISK:$PATH"

# -------------------------------------------------------------------
# OpenMP settings
# -------------------------------------------------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores

# For 32 threads, spread often gives better memory bandwidth on large AMD EPYC nodes.
# If performance is worse than before, test OMP_PROC_BIND=close.
export OMP_PROC_BIND=close

# Reduce possible glibc malloc memory overhead for many OpenMP threads.
export MALLOC_ARENA_MAX=4

# Useful diagnostic information in the output file.
echo "OMP_NUM_THREADS  = $OMP_NUM_THREADS"
echo "OMP_PROC_BIND    = $OMP_PROC_BIND"
echo "OMP_PLACES       = $OMP_PLACES"
echo "OMP_DYNAMIC      = $OMP_DYNAMIC"
echo "MALLOC_ARENA_MAX = $MALLOC_ARENA_MAX"

echo "CPU information:"
lscpu | egrep 'Model name|Socket|Core|Thread|NUMA|CPU\(s\)' || true

# -------------------------------------------------------------------
# Compile Basilisk VOF OpenMP code
# -------------------------------------------------------------------
rm -f rollwave

python thixo_base_state_explicitFrV_fixed.py --kappa 1e-4 --a 0.2 --Gamma 4.0 --FrV 4.00 --n-rows 128 --output profile.txt

qcc -std=c99 -fopenmp -Wall -O3 -march=native -DNDEBUG -Wdimensions \
    corrected_thixo_periodic_rollwave.c -o rollwave \
    -lm -grid=quadtree

echo "Compilation finished at: $(date)"

# -------------------------------------------------------------------
# Run
# -------------------------------------------------------------------
srun --ntasks=1 \
     --cpus-per-task=$SLURM_CPUS_PER_TASK \
     --cpu-bind=cores \
     ./rollwave

echo "Job finished at: $(date)"
