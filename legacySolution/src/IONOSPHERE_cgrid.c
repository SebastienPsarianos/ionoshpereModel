/* Coarse Grid Generator */

/* This subroutine determines a coarse grid
   given the fine grid.

   Parameters:

   ThetaF    Fine grid Theta-coordinate
   PsiF      Fine grid Psi-coordinate
   ThetaC    Coarse grid Theta-coordinate
   PsiC      Coarse grid Psi-coordinate
   nThetaC   Number of grid points in Theta-direction for
             coarse grid
   nPsiC     Number of grid points in Psi-direction for
             coarse grid */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_cgrid(ThetaF, PsiF,
                      ThetaC, PsiC, nThetaC, nPsiC)
double **ThetaF, **PsiF, **ThetaC, **PsiC;
int nThetaC, nPsiC;
{
        int ic,iif,jc,jf;

        /* Interior points */
        
        for (jf=3,jc=2; jc < nPsiC; jc++, jf+=2) {
          for (iif=3,ic=2; ic < nThetaC; ic++, iif+=2) {
            ThetaC[ic][jc]=ThetaF[iif][jf];
            PsiC[ic][jc]=PsiF[iif][jf];
                }
        }
        
        /* Boundary points */
        
        for (jc=1,ic=1; ic <= nThetaC; ic++, jc+=2) {
          ThetaC[ic][1]=ThetaF[jc][1];
          PsiC[ic][1]=PsiF[jc][1];
          ThetaC[ic][nPsiC]=ThetaF[jc][2*nPsiC-1];
          PsiC[ic][nPsiC]=PsiF[jc][2*nPsiC-1];
        }
        
        for (jc=1,ic=1; ic <= nPsiC; ic++, jc+=2) {
          ThetaC[1][ic]=ThetaF[1][jc];
          PsiC[1][ic]=PsiF[1][jc];
          ThetaC[nThetaC][ic]=ThetaF[2*nThetaC-1][jc];
          PsiC[nThetaC][ic]=PsiF[2*nThetaC-1][jc];
        }
}
