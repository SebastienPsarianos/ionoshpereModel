/* Two-Dimensional Array Copy */

/* This subroutine copies one two-dimensional array
   ain(1:nx,1:ny) to another two-dimensional array
   aout(1:nx,1:ny). */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_matcopy(aout, ain, nx, ny)
double **ain, **aout;
int nx, ny;
{
        int i,j;

        for (j=1; j<=ny; j++)
          for (i=1; i<=nx; i++)
            aout[i][j]=ain[i][j];

}
