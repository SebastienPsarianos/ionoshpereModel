/* Linear Elliptic Equation Residual Operator Res */

/* This subroutine uses evaluates the residual operator
   for the linear elliptic equation governing the
   electric field potential PHI
   using a standard centred-difference discretization.
   
   Given the solution PHI(1:nTheta,1:nPsi) and
   right-hand-side term RHO(1:nTheta,1:nPsi),
   the routine returns Res(u) in RES(1:nTheta,1:nPsi).
                                            
   Parameters:

   RES        Residual array
   PHI        Solution array
   RHO        Source array
   SigmaThTh  Conductance for Theta-coordinate direction
   SigmaThPs  Transverse conductance
   SigmaPsPs  Conductance for Psi-coordinate direction
   Theta      Discrete values of Theta coordinate
   Psi        Discrete values of Phi coordinate
   Radius     Radius of sphere, Rp
   nTheta     Number of grid points in Theta-direction
   nPsi       Number of grid points in Psi-direction
   RESIDUAL   Average residual
   iout       Residual output indicator */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_residual(RES, PHI, RHO, 
                         SigmaThTh, SigmaThPs, SigmaPsPs,
                         dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta,
                         dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi,
                         Theta, Psi, Radius, nTheta, nPsi, RESIDUAL, iout)
double **RES, **PHI, **RHO, 
       **SigmaThTh, **SigmaThPs, **SigmaPsPs,
       **dSigmaThTh_dTheta, **dSigmaThPs_dTheta, **dSigmaPsPs_dTheta,
       **dSigmaThTh_dPsi, **dSigmaThPs_dPsi, **dSigmaPsPs_dPsi,
       **Theta, **Psi, Radius, *RESIDUAL;
int nTheta, nPsi, iout;
{
        int i,j,jshift,iboundary;
        double dTheta, dPsi, dTheta2, dPsi2;
        double dd, dd2, sn, cs, sn2, cs2;
        double kappa_Theta2, kappa_Theta1,
               kappa_Psi2, kappa_Psi1, kappa_ratio=2.00;
        double RHOmax, RES0, RESPi;

        /* Determine mesh spacings. */

        dTheta=(Theta[nTheta][1]-Theta[1][1])/(double)(nTheta-1);
        dPsi=(Psi[1][nPsi]-Psi[1][1])/(double)(nPsi-1);
        dTheta2=dTheta*dTheta;
        dPsi2=dPsi*dPsi;
        dd=dTheta;
        dd2=dTheta2;

        /* Interior points */

        *RESIDUAL=0.0;
        RHOmax=0.00;
        for (j = 1; j <= nPsi; j++) {
          for (i = 2; i < nTheta; i++) {

            /* Evaluate effective diffusion coefficients. */

            sn = sin(Theta[i][j]);
            cs = cos(Theta[i][j]);
            sn2= sn*sn;
            cs2 = cs*cs;

            kappa_Theta2 =   SigmaThTh[i][j]*sn2;
            kappa_Theta1 =   SigmaThTh[i][j]*sn*cs
                           + dSigmaThTh_dTheta[i][j]*sn2
                           - dSigmaThPs_dPsi[i][j]*sn;
            kappa_Psi2 =   SigmaPsPs[i][j];
            kappa_Psi1 =   dSigmaThPs_dTheta[i][j]*sn
                         + dSigmaPsPs_dPsi[i][j];

            if (fabs(kappa_Theta1) >= kappa_ratio*fabs(kappa_Theta2)/dd) {
              if (kappa_Theta1 >= 0.00) {
                kappa_Theta1 = kappa_ratio*fabs(kappa_Theta2)/dd;
              } else {
                kappa_Theta1 = -kappa_ratio*fabs(kappa_Theta2)/dd;
              } /* endif. */
            } /* endif. */
            if (fabs(kappa_Psi1) >= kappa_ratio*fabs(kappa_Psi2)/dd) {
              if (kappa_Psi1 >= 0.00) {
                kappa_Psi1 = kappa_ratio*fabs(kappa_Psi2)/dd;
              } else {
                kappa_Psi1 = -kappa_ratio*fabs(kappa_Psi2)/dd;
              } /* endif. */
            } /* endif. */
       
            /* Evaluate residual. */

            if (j > 1 && j < nPsi) {
              RES[i][j]=RHO[i][j]-
                        kappa_Theta2*(PHI[i+1][j]-2.0*PHI[i][j]+PHI[i-1][j])/dd2-
                        0.50*kappa_Theta1*(PHI[i+1][j]-PHI[i-1][j])/dd-
                        kappa_Psi2*(PHI[i][j+1]-2.0*PHI[i][j]+PHI[i][j-1])/dd2-
                        0.50*kappa_Psi1*(PHI[i][j+1]-PHI[i][j-1])/dd;
            } else if (j == 1) {
              RES[i][j]=RHO[i][j]-
                        kappa_Theta2*(PHI[i+1][j]-2.0*PHI[i][j]+PHI[i-1][j])/dd2-
                        0.50*kappa_Theta1*(PHI[i+1][j]-PHI[i-1][j])/dd-
                        kappa_Psi2*(PHI[i][j+1]-2.0*PHI[i][j]+PHI[i][nPsi-1])/dd2-
                        0.50*kappa_Psi1*(PHI[i][j+1]-PHI[i][nPsi-1])/dd;
            } else {
              RES[i][j]=RHO[i][j]-
                        kappa_Theta2*(PHI[i+1][j]-2.0*PHI[i][j]+PHI[i-1][j])/dd2-
                        0.50*kappa_Theta1*(PHI[i+1][j]-PHI[i-1][j])/dd-
                        kappa_Psi2*(PHI[i][2]-2.0*PHI[i][j]+PHI[i][j-1])/dd2-
                        0.50*kappa_Psi1*(PHI[i][2]-PHI[i][j-1])/dd;
            } /* endif */

            *RESIDUAL+=fabs(RES[i][j]);
            if (fabs(RHO[i][j]) > RHOmax) {
              RHOmax=fabs(RHO[i][j]);
            } /* endif */

          } /* endfor */

        } /* endfor */
        
        /* Theta-coordinate Boundary points. */

        iboundary = 0;
        switch (iboundary)
          {
            case 0:
              /* Zero at poles. */
              for (j=1; j<=nPsi; j++) {
                 RES[1][j]=0.00;
                 *RESIDUAL+=RES[1][j];
                 if (fabs(RHO[1][j]) > RHOmax) {
                    RHOmax=fabs(RHO[1][j]);
                 } /* endif */

	         /* RES[nTheta][j]=0.00; */
                 RES[nTheta][j]=RHO[nTheta][j]/SigmaThTh[nTheta][j]-
                                   (PHI[nTheta][j]-PHI[nTheta-1][j])/dd;
                 *RESIDUAL+=fabs(RES[nTheta][j]);
                 if (fabs(RHO[nTheta][j]) > RHOmax) {
                    RHOmax=fabs(RHO[nTheta][j]);
                 } /* endif */
              } /* endfor */
              break;

            case 1:
              /* Simplified equation update. */
              for (j=1; j<=nPsi; j++) {
                 RES[1][j]=RHO[1][j]/SigmaThTh[1][j]-
                           (PHI[2][j]-PHI[1][j])/dd;
                 *RESIDUAL+=RES[1][j];
                 if (fabs(RHO[1][j]) > RHOmax) {
                    RHOmax=fabs(RHO[1][j]);
                 } /* endif */

                 RES[nTheta][j]=RHO[nTheta][j]/SigmaThTh[nTheta][j]-
                                (PHI[nTheta][j]-PHI[nTheta-1][j])/dd;
                 *RESIDUAL+=fabs(RES[nTheta][j]);
                 if (fabs(RHO[nTheta][j]) > RHOmax) {
                    RHOmax=fabs(RHO[nTheta][j]);
                 } /* endif */
              } /* endfor */
              break;

            case 2:
              /* Full equation update. */
              for (j=1; j<=nPsi; j++) {

                 sn = sin(Theta[1][j]);
                 cs = cos(Theta[1][j]);
                 sn2= sn*sn;
                 cs2 = cs*cs;

                 kappa_Theta2 =   SigmaThTh[1][j]*sn2;
                 kappa_Theta1 =   SigmaThTh[1][j]*sn*cs
                                + dSigmaThTh_dTheta[1][j]*sn2
                                - dSigmaThPs_dPsi[1][j]*sn;
                 kappa_Psi2 =   SigmaPsPs[1][j];
                 kappa_Psi1 =   dSigmaThPs_dTheta[1][j]*sn
                              + dSigmaPsPs_dPsi[1][j];

                 if ( j >= 1 + (nPsi - 1)/2 ) {
                    jshift = j + (nPsi - 1)/2;
                 } else {
                    jshift = j - (nPsi - 1)/2;
                 }

                 if (j > 1 && j < nPsi) {
                   RES[1][j]=RHO[1][j]-
                             kappa_Theta2*(PHI[2][j]-2.0*PHI[1][j]+PHI[2][jshift])/dd2-
                             0.50*kappa_Theta1*(PHI[2][j]-PHI[2][jshift])/dd-
                             kappa_Psi2*(PHI[1][j+1]-2.0*PHI[1][j]+PHI[1][j-1])/dd2-
                             0.50*kappa_Psi1*(PHI[1][j+1]-PHI[1][j-1])/dd;
                 } else if (j == 1) {
                   RES[1][j]=RHO[1][j]-
                             kappa_Theta2*(PHI[2][j]-2.0*PHI[1][j]+PHI[2][jshift])/dd2-
                             0.50*kappa_Theta1*(PHI[2][j]-PHI[2][jshift])/dd-
                             kappa_Psi2*(PHI[1][j+1]-2.0*PHI[1][j]+PHI[1][nPsi-1])/dd2-
                             0.50*kappa_Psi1*(PHI[1][j+1]-PHI[1][nPsi-1])/dd;
                 } else {
                   RES[1][j]=RHO[1][j]-
                             kappa_Theta2*(PHI[2][j]-2.0*PHI[1][j]+PHI[2][jshift])/dd2-
                             0.50*kappa_Theta1*(PHI[2][j]-PHI[2][jshift])/dd-
                             kappa_Psi2*(PHI[1][2]-2.0*PHI[1][j]+PHI[1][j-1])/dd2-
                             0.50*kappa_Psi1*(PHI[1][2]-PHI[1][j-1])/dd;
                 } /* endif */

                 *RESIDUAL+=RES[1][j];
                 if (fabs(RHO[1][j]) > RHOmax) {
                   RHOmax=fabs(RHO[1][j]);
                 } /* endif */

                 sn = sin(Theta[nTheta][j]);
                 cs = cos(Theta[nTheta][j]);
                 sn2= sn*sn;
                 cs2 = cs*cs;

                 kappa_Theta2 =   SigmaThTh[nTheta][j]*sn2;
                 kappa_Theta1 =   SigmaThTh[nTheta][j]*sn*cs
                                + dSigmaThTh_dTheta[nTheta][j]*sn2
                                - dSigmaThPs_dPsi[nTheta][j]*sn;
                 kappa_Psi2 =   SigmaPsPs[nTheta][j];
                 kappa_Psi1 =   dSigmaThPs_dTheta[nTheta][j]*sn
                              + dSigmaPsPs_dPsi[nTheta][j];

                 if ( j >= 1 + (nPsi - 1)/2 ) {
                    jshift = j + (nPsi - 1)/2;
                 } else {
                    jshift = j - (nPsi - 1)/2;
                 }

                 if (j > 1 && j < nPsi) {
                   RES[nTheta][j]=RHO[nTheta][j]-
                                  kappa_Theta2*(PHI[nTheta-1][jshift]-2.0*PHI[nTheta][j]+PHI[nTheta-1][j])/dd2-
                                  0.50*kappa_Theta1*(PHI[nTheta-1][jshift]-PHI[nTheta-1][j])/dd-
                                  kappa_Psi2*(PHI[nTheta][j+1]-2.0*PHI[nTheta][j]+PHI[nTheta][j-1])/dd2-
                                  0.50*kappa_Psi1*(PHI[nTheta][j+1]-PHI[nTheta][j-1])/dd;
                 } else if (j == 1) {
                   RES[nTheta][j]=RHO[nTheta][j]-
                                  kappa_Theta2*(PHI[nTheta-1][jshift]-2.0*PHI[nTheta][j]+PHI[nTheta-1][j])/dd2-
                                  0.50*kappa_Theta1*(PHI[nTheta-1][jshift]-PHI[nTheta-1][j])/dd-
                                  kappa_Psi2*(PHI[nTheta][j+1]-2.0*PHI[nTheta][j]+PHI[nTheta][nPsi-1])/dd2-
                                  0.50*kappa_Psi1*(PHI[nTheta][j+1]-PHI[nTheta][nPsi-1])/dd;
                 } else {
                   RES[nTheta][j]=RHO[nTheta][j]-
                                  kappa_Theta2*(PHI[nTheta-1][jshift]-2.0*PHI[nTheta][j]+PHI[nTheta-1][j])/dd2-
                                  0.50*kappa_Theta1*(PHI[nTheta-1][jshift]-PHI[nTheta-1][j])/dd-
                                  kappa_Psi2*(PHI[nTheta][2]-2.0*PHI[nTheta][j]+PHI[nTheta][j-1])/dd2-
                                  0.50*kappa_Psi1*(PHI[nTheta][2]-PHI[nTheta][j-1])/dd;
                 } /* endif */

                 *RESIDUAL+=fabs(RES[nTheta][j]);
                 if (fabs(RHO[nTheta][j]) > RHOmax) {
                   RHOmax=fabs(RHO[nTheta][j]);
                 } /* endif */

              } /* endfor */
              break;

            case 3:
              /* Average update. */
              RES0=0.00;
              RESPi=0.00;

              for (j=1; j<=nPsi; j++) {
                 RES0=RES0+RES[2][j];
                 RESPi=RESPi+RES[nTheta-1][j];
              } /* endfor */

              RES0=RES0/(double)nPsi;
              RESPi=RESPi/(double)nPsi;
 
              for (j=1; j<=nPsi; j++) {
                 RES[1][j]=RES0;
                 *RESIDUAL+=RES[1][j];
                 if (fabs(RHO[1][j]) > RHOmax) {
                    RHOmax=fabs(RHO[1][j]);
                 } /* endif */

	         RES[nTheta][j]=RESPi;
                 *RESIDUAL+=fabs(RES[nTheta][j]);
                 if (fabs(RHO[nTheta][j]) > RHOmax) {
                    RHOmax=fabs(RHO[nTheta][j]);
                 } /* endif */
              } /* endfor */
              break;

            default:
        } /* endswitch */

       /* Determine the average normalized residual. */

        if (RHOmax <= 0.00) {
          RHOmax=1.00;
        } /* endif */
        *RESIDUAL=*RESIDUAL/(nTheta*nPsi*RHOmax);

}
