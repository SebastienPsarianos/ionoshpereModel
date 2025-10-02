/* Two-Dimensional Array Zero */

/* This subroutine zeroes the entries of the
   two-dimensional array u(1:nx,1:ny). */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_matzero(u, nx, ny)
double **u;
int nx, ny;
{
        int i,j;

        for (j=1; j<=ny; j++)
          for (i=1; i<=nx; i++)
            u[i][j]=0.0;
}
