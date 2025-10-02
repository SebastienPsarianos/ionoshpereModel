/* Determine Field-Aligned and Hall and Pedersen Conductivities */

/* This subroutine computes the height-integrated field-aligned and
   Hall and Pedersen conductances for the ionosphere at each
   location of the discretized solution domain.  The gradients of
   these quantities are also computed.

   Parameters:

   Sigma0  Field-aligned conductance (height-integrated conductivity)
   SigmaH  Hall Conductance (height-integrated conductivity)
   SigmaP  Petersen Conductance (height-integrated conductivity)
   PHI     Potential solution array
   Theta   Discrete values of Theta coordinate for fine grid
   Psi     Discrete values of Phi coordinate for fine grid
   Radius  Radius of sphere, Rp
   nTheta  Number of grid points in Theta-direction
           on fine grid
   nPsi    Number of grid points in Psi-direction
           on fine grid */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_conductance(Sigma0, SigmaH, SigmaP,
                            SigmaThTh, SigmaThPs, SigmaPsPs,
                            dSigmaThTh_dTheta, dSigmaThPs_dTheta,
                            dSigmaPsPs_dTheta,
                            dSigmaThTh_dPsi, dSigmaThPs_dPsi,
                            dSigmaPsPs_dPsi,
                            PHI, Theta, Psi, Radius, nTheta, nPsi)
double **Sigma0, **SigmaH, **SigmaP,
       **SigmaThTh, **SigmaThPs, **SigmaPsPs,
       **dSigmaThTh_dTheta, **dSigmaThPs_dTheta, **dSigmaPsPs_dTheta,
       **dSigmaThTh_dPsi, **dSigmaThPs_dPsi, **dSigmaPsPs_dPsi,
       **PHI, **Theta, **Psi, Radius;
int nTheta, nPsi;
{
        int i,j;
        double dTheta, dPsi, dTheta2, dPsi2;
        double dd, dd2, sn, cs, sn2, cs2, cs3, cs4, C;

        dTheta=(Theta[nTheta][1]-Theta[1][1])/(double)(nTheta-1);
        dPsi=(Psi[1][nPsi]-Psi[1][1])/(double)(nPsi-1);
        dTheta2=dTheta*dTheta;
        dPsi2=dPsi*dPsi;
        dd=dTheta;
        dd2=dTheta2;

        for (j=1; j<=nPsi; j++) {
          for (i=1; i<=nTheta; i++) {
            Sigma0[i][j]=20.00;
            SigmaH[i][j]=2.00;
            SigmaP[i][j]=2.00;

            sn = sin(Theta[i][j]);
            cs = cos(Theta[i][j]);
            sn2= sn*sn;
            cs2 = cs*cs;
            cs3 = 1.00 + 3.00*cs2;
            cs4 = sqrt(cs3);
            C = 4.00*Sigma0[i][j]*cs2 +
                SigmaP[i][j]*sn2;

	    SigmaThTh[i][j]=Sigma0[i][j]*SigmaP[i][j]*cs3/C;
            SigmaThPs[i][j]=2.00*Sigma0[i][j]*SigmaH[i][j]*
                            cs*cs4/C;
            SigmaPsPs[i][j]=SigmaP[i][j]+
                            SigmaH[i][j]*SigmaH[i][j]*
                            sn2/C;

            /* SigmaThTh[i][j]=SigmaP[i][j]*cs3/(4.00*cs2);
            SigmaThPs[i][j]=SigmaH[i][j]*cs4/(2.00*cs);
            SigmaPsPs[i][j]=SigmaP[i][j]; */

            /* SigmaThTh[i][j]=1.00;
            SigmaThPs[i][j]=2.00;
            SigmaPsPs[i][j]=1.00; */

            dSigmaThTh_dTheta[i][j]=0.00;
            dSigmaThTh_dPsi[i][j]=0.00;
            dSigmaThPs_dTheta[i][j]=0.00;
            dSigmaThPs_dPsi[i][j]=0.00;
            dSigmaPsPs_dTheta[i][j]=0.00;
            dSigmaPsPs_dPsi[i][j]=0.00;
          } /* endfor */
        } /* endfor */

        for (j=1; j<=nPsi; j++) {
          if (j > 1 && j < nPsi ) {
            for (i=2; i<nTheta; i++) {
              dSigmaThTh_dTheta[i][j]=(SigmaThTh[i+1][j]-SigmaThTh[i-1][j])/
                                      (2.00*dd);
              dSigmaThTh_dPsi[i][j]=(SigmaThTh[i][j+1]-SigmaThTh[i][j-1])/
                                    (2.00*dd);
 
              dSigmaThPs_dTheta[i][j]=(SigmaThPs[i+1][j]-SigmaThPs[i-1][j])/
                                      (2.00*dd);
              dSigmaThPs_dPsi[i][j]=(SigmaThPs[i][j+1]-SigmaThPs[i][j-1])/
                                    (2.00*dd);
 
              dSigmaPsPs_dTheta[i][j]=(SigmaPsPs[i+1][j]-SigmaPsPs[i-1][j])/
                                      (2.00*dd);
              dSigmaPsPs_dPsi[i][j]=(SigmaPsPs[i][j+1]-SigmaPsPs[i][j-1])/
                                    (2.00*dd);
            } /* endfor */
          } else if (j == 1) {
            for (i=2; i<nTheta; i++) {
              dSigmaThTh_dTheta[i][j]=(SigmaThTh[i+1][j]-SigmaThTh[i-1][j])/
                                      (2.00*dd);
              dSigmaThTh_dPsi[i][j]=(SigmaThTh[i][j+1]-SigmaThTh[i][nPsi-1])/
                                    (2.00*dd);
 
              dSigmaThPs_dTheta[i][j]=(SigmaThPs[i+1][j]-SigmaThPs[i-1][j])/
                                      (2.00*dd);
              dSigmaThPs_dPsi[i][j]=(SigmaThPs[i][j+1]-SigmaThPs[i][nPsi-1])/
                                    (2.00*dd);
 
              dSigmaPsPs_dTheta[i][j]=(SigmaPsPs[i+1][j]-SigmaPsPs[i-1][j])/
                                      (2.00*dd);
              dSigmaPsPs_dPsi[i][j]=(SigmaPsPs[i][j+1]-SigmaPsPs[i][nPsi-1])/
                                    (2.00*dd);
            } /* endfor */
          } else {
            for (i=2; i<nTheta; i++) {
              dSigmaThTh_dTheta[i][j]=(SigmaThTh[i+1][j]-SigmaThTh[i-1][j])/
                                      (2.00*dd);
              dSigmaThTh_dPsi[i][j]=(SigmaThTh[i][2]-SigmaThTh[i][j-1])/
                                    (2.00*dd);
 
              dSigmaThPs_dTheta[i][j]=(SigmaThPs[i+1][j]-SigmaThPs[i-1][j])/
                                      (2.00*dd);
              dSigmaThPs_dPsi[i][j]=(SigmaThPs[i][2]-SigmaThPs[i][j-1])/
                                    (2.00*dd);
 
              dSigmaPsPs_dTheta[i][j]=(SigmaPsPs[i+1][j]-SigmaPsPs[i-1][j])/
                                      (2.00*dd);
              dSigmaPsPs_dPsi[i][j]=(SigmaPsPs[i][2]-SigmaPsPs[i][j-1])/
                                    (2.00*dd);
            } /* endfor */
          } /* endif */
        } /* endfor */
}
