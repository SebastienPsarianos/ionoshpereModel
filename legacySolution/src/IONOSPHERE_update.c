/* Solution Update Operator */

/* This subroutine adds the prolonged solution correction
   from the coarse grid to the fine grid solution and
   thereby updates the next approximation to the solution
   on the fine grid.

   Given the solution PHI(1:nTheta,1:nPsi) and
   prolonged solution corrections RES(1:nTheta,1:nPsi),
   adds RES to PHI.
                                            
   Parameters:

   PHI     Solution array
   RES     Solution correction array
   OMEGA   Relaxation parameter
   nThetaF Number of grid points in Theta-direction
           on fine grid
   nPsiF   Number of grid points in Psi-direction
           on fine grid */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_update(PHI, RES, OMEGA, nThetaF, nPsiF)
double **PHI, **RES, OMEGA;
int nThetaF, nPsiF;
{
        int i,j;

        /* Update solution on fine grid. */
        
        for (j=1; j<=nPsiF; j++) {
	  for (i=1; i<=nThetaF; i++) {
            PHI[i][j] += OMEGA*RES[i][j];
          } /* endfor. */
        } /* endfor. */
}
