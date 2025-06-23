# gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt

# cat boxes.txt >> pl_ext.gfs

#!/bin/sh
gerris2D -m -s 3  pl_ext.gfs > splitsimulation.gfs
gerris2D -b 8 splitsimulation.gfs > parallelsimulation.gfs
mpirun -np 8 gerris2D parallelsimulation.gfs
#
# gerris2D -m -b 3 pl_ext.gfs > parallelsimulation.gfs
# mpirun -np 3 gerris2D parallelsimulation.gfs

# gerris2D -m pl_ext.gfs

# mpirun -np 3 gerris2D -m pl_ext.gfs
