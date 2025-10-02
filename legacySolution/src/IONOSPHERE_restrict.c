/* Half-Weighting Restriction Operator */

/* This subroutine uses the half-weighting restriction operator
   to restrict the fine grid solution to the coarse grid.

   Parameters:

   UF        Fine-grid solution array
   UC        Coarse-grid solution array
   nThetaC   Number of grid points in Theta-direction for
             coarse grid
   nPsiC     Number of grid points in Psi-direction for
             coarse grid */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_restrict(UF, UC, nThetaC, nPsiC)
double **UF, **UC;
int nThetaC, nPsiC;
{
        int ic,iif,jc,jf;

        /* Interior points */
        
        for (jf=3,jc=2; jc < nPsiC; jc++, jf+=2) {
          for (iif=3,ic=2; ic < nThetaC; ic++, iif+=2) {
            UC[ic][jc]=0.5*UF[iif][jf]+
                       0.125*(UF[iif+1][jf]+UF[iif-1][jf]
                             +UF[iif][jf+1]+UF[iif][jf-1]);
                }
        }
        
        /* Boundary points */
        
        for (jc=1,ic=1; ic <= nThetaC; ic++, jc+=2) {
          UC[ic][1]=UF[jc][1];
          UC[ic][nPsiC]=UF[jc][2*nPsiC-1];
        }
        
        for (jc=1,ic=1; ic <= nPsiC; ic++, jc+=2) {
          UC[1][ic]=UF[1][jc];
          UC[nThetaC][ic]=UF[2*nThetaC-1][jc];
        }
}
