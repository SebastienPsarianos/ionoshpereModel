/* Fine Grid Generator */

/* This subroutine creates a uniform mesh in the spherical
   coordinate system for the upper hemisphere of a spherical
   shell which is used as the fine grid in the ionosphere
   current calculations.

   Parameters:

   Theta   Discrete values of Theta coordinate for fine grid
   Psi     Discrete values of Phi coordinate for fine grid
   X,Y,Z   Cartesian coordinates (x,y,z) of grid points
   Radius  Radius of sphere, Rp
   nTheta  Number of grid points in Theta-direction
           on fine grid
   nPsi    Number of grid points in Psi-direction
           on fine grid
   fhem    Fraction of hemisphere over which the solution
           is defined. */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_fgrid(Theta, Psi, X, Y, Z,
                      Radius,
                      nTheta, nPsi, fhem)
double **Theta, **Psi, **X, **Y, **Z, Radius, fhem;
int nTheta, nPsi;
{
        int i,j;
        double dTheta, dPsi;
        
        dPsi=2.00*IONOSPHERE_PI/(double)(nPsi-1);
        dTheta=dPsi;

        for (j=1;j<=nPsi;j++) {
          Theta[1][j] = IONOSPHERE_Theta_0;
          Psi[1][j] = (double)(j-1)*dPsi;
          X[1][j] = sin(Theta[1][j])*cos(Psi[1][j]);
          Y[1][j] = sin(Theta[1][j])*sin(Psi[1][j]);
          Z[1][j] = cos(Theta[1][j]);
  
          for (i=2;i<nTheta;i++) {
            Theta[i][j] = (double)(i-1)*dTheta;
            Psi[i][j] = (double)(j-1)*dPsi;
            X[i][j] = sin(Theta[i][j])*cos(Psi[i][j]);
            Y[i][j] = sin(Theta[i][j])*sin(Psi[i][j]);
            Z[i][j] = cos(Theta[i][j]);
          }

          Theta[nTheta][j] = fhem*IONOSPHERE_PI - IONOSPHERE_Theta_0;
          Psi[nTheta][j] = (double)(j-1)*dPsi;
          X[nTheta][j] = sin(Theta[nTheta][j])*cos(Psi[nTheta][j]);
          Y[nTheta][j] = sin(Theta[nTheta][j])*sin(Psi[nTheta][j]);
          Z[nTheta][j] = cos(Theta[nTheta][j]);
        }
}
