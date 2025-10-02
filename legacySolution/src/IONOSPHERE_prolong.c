/* Bilinear Interpolation Prolongation Operator */

/* This subroutine uses bilinear interpolation
   to prolong the coarse grid solution to the fine grid.

   Parameters:

   UF        Fine-grid solution array
   UC        Coarse-grid solution array
   nThetaF   Number of grid points in Theta-direction for
             fine grid
   nPsiF     Number of grid points in Psi-direction for
             fine grid */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_prolong(UF, UC, nThetaF, nPsiF)
double **UF, **UC;
int nThetaF, nPsiF;
{
        int ic,iif,jc,jf,nThetaC,nPsiC;
        nThetaC=nThetaF/2+1;
        nPsiC=nPsiF/2+1;

        /* Do elements that are direct copies. */
        
        for (jc=1, jf=1; jc<=nPsiC; jc++,jf+=2)
          for (ic=1; ic<=nThetaC; ic++) UF[2*ic-1][jf]=UC[ic][jc];

        /* Do odd-numbered columns, interpolating vertically. */
        
        for (jf=1; jf<=nPsiF; jf+=2)
          for (iif=2; iif<nThetaF; iif+=2)
            UF[iif][jf]=0.5*(UF[iif+1][jf]+UF[iif-1][jf]);

        /* Do even-numbered columns, interpolating horizontally. */

        for (jf=2; jf<nPsiF; jf+=2)
          for (iif=1; iif<=nThetaF; iif++)
            UF[iif][jf]=0.5*(UF[iif][jf+1]+UF[iif][jf-1]);
}
