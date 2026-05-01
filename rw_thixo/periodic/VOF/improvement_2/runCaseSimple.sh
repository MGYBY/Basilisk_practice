# python3 thixo_base_state_explicitFrV_fixed.py --kappa 1e-4 --a 0.2 --Gamma 8 --FrV 4.00 --n-rows 128 --output profile.txt

export OMP_NUM_THREADS=10
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_DYNAMIC=FALSE
qcc  -std=c99 -fopenmp -Wall -O2 -Wdimensions corrected_thixo_periodic_rollwave.c -o rollwave -lm -grid=quadtree -L$BASILISK/gl -lglutils -lfb_tiny
# serial
# qcc -O2 -Wall -Wdimensions -o rollwave thixo_periodic_rollwave_mu_capped_noview_restartfixed5.c -lm
./rollwave
