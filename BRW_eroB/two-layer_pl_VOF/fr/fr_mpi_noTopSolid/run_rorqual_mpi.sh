#!/usr/bin/env bash
#SBATCH --account=def-sushama-ab_cpu
#SBATCH --job-name=front_runner_MPI_noSolid_v2p7
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --mem=38368M
#SBATCH --time=0-23:56
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=boyuan.yu@mail.mcgill.ca
#SBATCH --hint=nomultithread

# Pure MPI: 80 ranks, matching log.21121849; one CPU/thread per rank.
# Memory is a total node request, equal to the original 88 * 436 MiB.
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]}")}"
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 gnuplot python
export BASILISK=${BASILISK:-/home/yboyuan/links/projects/rrg-sushama-ab/yboyuan/basilisk/src}
export PATH="$BASILISK:$PATH"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MALLOC_ARENA_MAX=2
: "${SLURM_NTASKS:?Submit this script using sbatch}"
[[ ${SLURM_CPUS_PER_TASK:-1} == 1 ]] || { echo 'Pure MPI requires cpus-per-task=1.' >&2; exit 2; }
export EXPECTED_MPI_RANKS="$SLURM_NTASKS"
echo "Job started $(date -Is); requested MPI ranks=$SLURM_NTASKS; OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "Compiling v2p7 no-top-solid case; the simulation has not started yet."
echo "Host=$(hostname); cwd=$PWD; BASILISK=$BASILISK"
bash ./build_case.sh mpi rollwave
# No -parallel argument: this case's main(void) does not parse it.
# Binding itself is unchanged; suppress one verbose startup line per rank.
binding_report=()
if [[ ${REPORT_BINDINGS:-0} == 1 ]]; then binding_report=(--report-bindings); fi
echo "Build succeeded; launching MPI simulation $(date -Is)."
mpirun -np "$SLURM_NTASKS" --bind-to core --map-by core "${binding_report[@]}" \
  ./rollwave > "log.${SLURM_JOB_ID:-manual}" 2>&1
echo "Completed $(date -Is)"
