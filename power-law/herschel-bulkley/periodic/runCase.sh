CC99='mpicc -std=c99' qcc -Wall -O2 -Wdimensions -D_MPI=1 rw.c -o rw -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
mpirun --oversubscribe -np 8 ./rw -parallel > log

# make rw.tst
