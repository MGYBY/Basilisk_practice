#include "grid/multigrid3D.h"
#include "embed.h"
#include "navier-stokes/centered.h"
// #define FILTERED
// #include "two-phasePL.h"
#include "two-phase.h"
// #include "./myTension.h"
#include "tension.h"
// #include "vof.h"
// alternatively, use momentum-conserving scheme
#include "navier-stokes/conserving.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"

int main() {
    dimensions (nx = 3);
    periodic (right);
    init_grid (24);

    char names[36];
    sprintf( names, "field-%d.dat", pid() );
    FILE * fp2 = fopen (names, "w");
    foreach ()
    {
        fprintf(fp2, "%g %g %g\n", x, y, z);
    }
    fclose(fp2);
    char command[80];
    sprintf(command, "LC_ALL=C  cat field* > ALLFIELD-%g.dat",t);
    system(command);// allow to use linux command in the c code to concatenate our files

}

