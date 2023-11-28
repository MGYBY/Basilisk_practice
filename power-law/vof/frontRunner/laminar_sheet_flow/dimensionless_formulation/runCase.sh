# CC99='mpicc -std=c99' qcc -Wall -O2 -Wdimensions -D_MPI=1 rw.c -o rw -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
# mpirun --oversubscribe -np 14 ./rw -parallel > log

# make rw.tst

export OMP_NUM_THREADS=15
qcc  -std=c99 -fopenmp -Wall -O2 rw.c -o rollwave -lm -grid=octree -L$BASILISK/gl -lglutils -lfb_tiny
./rollwave
