2D VOF simulations.

Some key-points:
1. htg output applicable in MPI:
http://basilisk.fr/sandbox/sander/output_htg.h
2. Hyperthread parallelism:
```
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 example.c -o example -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
mpirun --oversubscribe -np 12 ./example -parallel > log
```
3. Interface output:
* [$wmostert/shallow.c](http://basilisk.fr/sandbox/wmostert/shallow.c)
* [$vatsal/SurfaceWavesBingham.c](http://basilisk.fr/sandbox/vatsal/GenaralizedNewtonian/SurfaceWavesBingham.c)
4. Good bview examples:

5. Other visualization & PP examples:
