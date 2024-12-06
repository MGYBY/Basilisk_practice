# export OMP_NUM_THREADS=3
# qcc  -std=c99 -fopenmp -Wall -Wdimensions -O2 test_mpi.c -o test_mpi -lm  -L$BASILISK/gl -lglutils -lfb_tiny
# ./test_mpi

CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 test_mpi.c -o test_mpi -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
mpirun --oversubscribe -np 3 ./test_mpi -parallel > log
