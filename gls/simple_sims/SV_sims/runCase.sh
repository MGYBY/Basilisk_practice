CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 lo_ss.c -o lo_ss -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm  -L$BASILISK/ppr -lppr -lgfortran
mpirun --oversubscribe -np 8 ./lo_ss -parallel > log
