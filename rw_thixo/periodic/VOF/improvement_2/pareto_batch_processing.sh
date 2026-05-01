#!/bin/bash
#SBATCH --account=def-sushama-ab_cpu
#SBATCH --ntasks-per-node=34      # number of MPI processes
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=145M           # memory for the entire job across all cores 3900M
#SBATCH --time=0-08:10      # time (DD-HH:MM)
#SBATCH --output=%N-%j.out  # %N for node name, %j for jobID
#SBATCH --mail-type=ALL               # Type of email notification- BEGIN,END,F$
#SBATCH --mail-user=boyuan.yu@mail.mcgill.ca   # Email to which notifications will be $

module purge
module load openmpi
module load StdEnv/2023
module load gcc
module load gnuplot

export BASILISK=/home/yboyuan/links/projects/rrg-sushama-ab/yboyuan/basilisk/src
export PATH=$PATH:$BASILISK

export OMP_NUM_THREADS=32
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_DYNAMIC=FALSE
qcc  -std=c99 -fopenmp -Wall -O2 -Wdimensions corrected_thixo_periodic_rollwave.c -o rollwave -lm -grid=quadtree -L$BASILISK/gl -lglutils -lfb_tiny
# serial
# qcc -O2 -Wall -Wdimensions -o rollwave thixo_periodic_rollwave_mu_capped_noview_restartfixed5.c -lm
./rollwave

