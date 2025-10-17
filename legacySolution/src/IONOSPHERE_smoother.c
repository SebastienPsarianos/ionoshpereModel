/* Red-Black Gauss-Seidel Smoothing Operator */

/* This subroutine uses a red-black Gauss-Seidel scheme
   to smooth the solution to the linear elliptic equation
   governing the electric field potential PHI
   as defined on the domain 0 < theta < PI, 0 < psi < 2 PI,
   and subject to the Dirichlet boundary data

        PHI(0,Psi) = 0 ,

        PHI(PI,Psi) = 0 ,

        PHI(Theta,0) = PHI(Theta,2*PI).

   The current value of the solution, PHI, is updated using
   Gauss-Seidel centred-difference discretization procedure
   and the right-hand side function, RHO.

   Parameters:

   PHI        Solution array
   RHO        Source array
   SigmaThTh  Conductance for Theta-coordinate direction
   SigmaThPs  Transverse conductance
   SigmaPsPs  Conductance for Psi-coordinate direction
   Theta      Discrete values of Theta coordinate
   Psi        Discrete values of Phi coordinate
   Radius     Radius of sphere, Rp
   nTheta     Number of grid points in Theta-direction
   nPsi       Number of grid points in Psi-direction */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_smoother(PHI, RHO, SigmaThTh, SigmaThPs, SigmaPsPs,
                         dSigmaThTh_dTheta, dSigmaThPs_dTheta,
                         dSigmaPsPs_dTheta, dSigmaThTh_dPsi, dSigmaThPs_dPsi,
                         dSigmaPsPs_dPsi, Theta, Psi, Radius, nTheta,
                         nPsi) double **PHI,
    **RHO, **SigmaThTh, **SigmaThPs, **SigmaPsPs, **dSigmaThTh_dTheta,
    **dSigmaThPs_dTheta, **dSigmaPsPs_dTheta, **dSigmaThTh_dPsi,
    **dSigmaThPs_dPsi, **dSigmaPsPs_dPsi, **Theta, **Psi, Radius;
int nTheta, nPsi;
{
  int i, ipass, isw, j, jsw, jshift, iboundary;
  double dTheta, dPsi, dTheta2, dPsi2;
  double dd, dd2, sn, cs, sn2, cs2;
  double kappa_Theta2, kappa_Theta1, kappa_Psi2, kappa_Psi1, kappa_ratio = 2.00;
  double Phi0, PhiPi;

  /* Determine mesh spacings. */

  dTheta = (Theta[nTheta][1] - Theta[1][1]) / (double)(nTheta - 1);
  dPsi = (Psi[1][nPsi] - Psi[1][1]) / (double)(nPsi - 1);
  dTheta2 = dTheta * dTheta;
  dPsi2 = dPsi * dPsi;
  dd = dTheta;
  dd2 = dTheta2;

  /* Red and black sweeps. */

  jsw = 1;
  for (ipass = 1; ipass <= 2; ipass++, jsw = 3 - jsw) {
    isw = jsw;
    for (j = 1; j <= nPsi; j++, isw = 3 - isw) {
      for (i = isw + 1; i < nTheta; i += 2) {

        /* Evaluate effective diffusion coefficients. */

        sn = sin(Theta[i][j]);
        cs = cos(Theta[i][j]);
        sn2 = sn * sn;
        cs2 = cs * cs;

        kappa_Theta2 = SigmaThTh[i][j] * sn2;
        kappa_Theta1 = SigmaThTh[i][j] * sn * cs +
                       dSigmaThTh_dTheta[i][j] * sn2 -
                       dSigmaThPs_dPsi[i][j] * sn;
        kappa_Psi2 = SigmaPsPs[i][j];
        kappa_Psi1 = dSigmaThPs_dTheta[i][j] * sn + dSigmaPsPs_dPsi[i][j];

        if (fabs(kappa_Theta1) >= kappa_ratio * fabs(kappa_Theta2) / dd) {
          if (kappa_Theta1 >= 0.00) {
            kappa_Theta1 = kappa_ratio * fabs(kappa_Theta2) / dd;
          } else {
            kappa_Theta1 = -kappa_ratio * fabs(kappa_Theta2) / dd;
          } /* endif. */
        } /* endif. */
        if (fabs(kappa_Psi1) >= kappa_ratio * fabs(kappa_Psi2) / dd) {
          if (kappa_Psi1 >= 0.00) {
            kappa_Psi1 = kappa_ratio * fabs(kappa_Psi2) / dd;
          } else {
            kappa_Psi1 = -kappa_ratio * fabs(kappa_Psi2) / dd;
          } /* endif. */
        } /* endif. */

        /* Gauss-Seidel formula. */

        if (j > 1 && j < nPsi) {
          PHI[i][j] =
              (kappa_Theta2 * (PHI[i + 1][j] + PHI[i - 1][j]) +
               0.50 * dd * kappa_Theta1 * (PHI[i + 1][j] - PHI[i - 1][j]) +
               kappa_Psi2 * (PHI[i][j + 1] + PHI[i][j - 1]) +
               0.50 * dd * kappa_Psi1 * (PHI[i][j + 1] - PHI[i][j - 1]) -
               dd2 * RHO[i][j]) /
              (2.00 * (kappa_Theta2 + kappa_Psi2));
        } else if (j == 1) {
          PHI[i][j] =
              (kappa_Theta2 * (PHI[i + 1][j] + PHI[i - 1][j]) +
               0.50 * dd * kappa_Theta1 * (PHI[i + 1][j] - PHI[i - 1][j]) +
               kappa_Psi2 * (PHI[i][j + 1] + PHI[i][nPsi - 1]) +
               0.50 * dd * kappa_Psi1 * (PHI[i][j + 1] - PHI[i][nPsi - 1]) -
               dd2 * RHO[i][j]) /
              (2.00 * (kappa_Theta2 + kappa_Psi2));
        } else {
          PHI[i][j] =
              (kappa_Theta2 * (PHI[i + 1][j] + PHI[i - 1][j]) +
               0.50 * dd * kappa_Theta1 * (PHI[i + 1][j] - PHI[i - 1][j]) +
               kappa_Psi2 * (PHI[i][2] + PHI[i][j - 1]) +
               0.50 * dd * kappa_Psi1 * (PHI[i][2] - PHI[i][j - 1]) -
               dd2 * RHO[i][j]) /
              (2.00 * (kappa_Theta2 + kappa_Psi2));
        } /* endif */
      } /* endfor */
    } /* endfor */
  } /* endfor */

  /* Theta-coordinate boundary points. */

  iboundary = 0;
  switch (iboundary) {
  case 0:
    /* Zero at poles. */
    for (j = 1; j <= nPsi; j++) {
      PHI[1][j] = 0.00;
      /* PHI[nTheta][j]=0.00; */
      PHI[nTheta][j] =
          (PHI[nTheta - 1][j] + dd * RHO[nTheta][j] / SigmaThTh[nTheta][j]);
    } /* endfor */
    break;

  case 1:
    /* Simplified equation update. */
    Phi0 = 0.00;
    PhiPi = 0.00;

    for (j = 1; j <= nPsi; j++) {
      PHI[1][j] = (PHI[2][j] - dd * RHO[1][j] / SigmaThTh[1][j]);
      Phi0 = Phi0 + PHI[1][j];

      PHI[nTheta][j] =
          (PHI[nTheta - 1][j] + dd * RHO[nTheta][j] / SigmaThTh[nTheta][j]);
      PhiPi = PhiPi + PHI[nTheta][j];
    } /* endfor */

    Phi0 = Phi0 / (double)nPsi;
    PhiPi = PhiPi / (double)nPsi;

    for (j = 1; j <= nPsi; j++) {
      PHI[1][j] = Phi0;
      PHI[nTheta][j] = PhiPi;
    } /* endfor */
    break;

  case 2:
    /* Full equation update. */
    Phi0 = 0.00;
    PhiPi = 0.00;

    for (j = 1; j <= nPsi; j++) {

      sn = sin(Theta[1][j]);
      cs = cos(Theta[1][j]);
      sn2 = sn * sn;
      cs2 = cs * cs;

      kappa_Theta2 = SigmaThTh[1][j] * sn2;
      kappa_Theta1 = SigmaThTh[1][j] * sn * cs + dSigmaThTh_dTheta[1][j] * sn2 -
                     dSigmaThPs_dPsi[1][j] * sn;
      kappa_Psi2 = SigmaPsPs[1][j];
      kappa_Psi1 = dSigmaThPs_dTheta[1][j] * sn + dSigmaPsPs_dPsi[1][j];

      if (j >= 1 + (nPsi - 1) / 2) {
        jshift = j + (nPsi - 1) / 2;
      } else {
        jshift = j - (nPsi - 1) / 2;
      }

      if (j > 1 && j < nPsi) {
        PHI[1][j] = (kappa_Theta2 * (PHI[2][j] + PHI[2][jshift]) +
                     0.50 * dd * kappa_Theta1 * (PHI[2][j] - PHI[2][jshift]) +
                     kappa_Psi2 * (PHI[1][j + 1] + PHI[1][j - 1]) +
                     0.50 * dd * kappa_Psi1 * (PHI[1][j + 1] - PHI[1][j - 1]) -
                     dd2 * RHO[1][j]) /
                    (2.00 * (kappa_Theta2 + kappa_Psi2));
      } else if (j == 1) {
        PHI[1][j] =
            (kappa_Theta2 * (PHI[2][j] + PHI[2][jshift]) +
             0.50 * dd * kappa_Theta1 * (PHI[2][j] - PHI[2][jshift]) +
             kappa_Psi2 * (PHI[1][j + 1] + PHI[1][nPsi - 1]) +
             0.50 * dd * kappa_Psi1 * (PHI[1][j + 1] - PHI[1][nPsi - 1]) -
             dd2 * RHO[1][j]) /
            (2.00 * (kappa_Theta2 + kappa_Psi2));
      } else {
        PHI[1][j] = (kappa_Theta2 * (PHI[2][j] + PHI[2][jshift]) +
                     0.50 * dd * kappa_Theta1 * (PHI[2][j] - PHI[2][jshift]) +
                     kappa_Psi2 * (PHI[1][2] + PHI[1][j - 1]) +
                     0.50 * dd * kappa_Psi1 * (PHI[1][2] - PHI[1][j - 1]) -
                     dd2 * RHO[1][j]) /
                    (2.00 * (kappa_Theta2 + kappa_Psi2));
      } /* endif */

      Phi0 = Phi0 + PHI[1][j];

      sn = sin(Theta[nTheta][j]);
      cs = cos(Theta[nTheta][j]);
      sn2 = sn * sn;
      cs2 = cs * cs;

      kappa_Theta2 = SigmaThTh[nTheta][j] * sn2;
      kappa_Theta1 = SigmaThTh[nTheta][j] * sn * cs +
                     dSigmaThTh_dTheta[nTheta][j] * sn2 -
                     dSigmaThPs_dPsi[nTheta][j] * sn;
      kappa_Psi2 = SigmaPsPs[nTheta][j];
      kappa_Psi1 =
          dSigmaThPs_dTheta[nTheta][j] * sn + dSigmaPsPs_dPsi[nTheta][j];

      if (j >= 1 + (nPsi - 1) / 2) {
        jshift = j + (nPsi - 1) / 2;
      } else {
        jshift = j - (nPsi - 1) / 2;
      }

      if (j > 1 && j < nPsi) {
        PHI[nTheta][j] =
            (kappa_Theta2 * (PHI[nTheta - 1][jshift] + PHI[nTheta - 1][j]) +
             0.50 * dd * kappa_Theta1 *
                 (PHI[nTheta - 1][jshift] - PHI[nTheta - 1][j]) +
             kappa_Psi2 * (PHI[nTheta][j + 1] + PHI[nTheta][j - 1]) +
             0.50 * dd * kappa_Psi1 *
                 (PHI[nTheta][j + 1] - PHI[nTheta][j - 1]) -
             dd2 * RHO[nTheta][j]) /
            (2.00 * (kappa_Theta2 + kappa_Psi2));
      } else if (j == 1) {
        PHI[nTheta][j] =
            (kappa_Theta2 * (PHI[nTheta - 1][jshift] + PHI[nTheta - 1][j]) +
             0.50 * dd * kappa_Theta1 *
                 (PHI[nTheta - 1][j] - PHI[nTheta - 1][jshift]) +
             kappa_Psi2 * (PHI[nTheta][j + 1] + PHI[nTheta][nPsi - 1]) +
             0.50 * dd * kappa_Psi1 *
                 (PHI[nTheta][j + 1] - PHI[nTheta][nPsi - 1]) -
             dd2 * RHO[nTheta][j]) /
            (2.00 * (kappa_Theta2 + kappa_Psi2));
      } else {
        PHI[nTheta][j] =
            (kappa_Theta2 * (PHI[nTheta - 1][jshift] + PHI[nTheta - 1][j]) +
             0.50 * dd * kappa_Theta1 *
                 (PHI[nTheta - 1][jshift] - PHI[nTheta - 1][j]) +
             kappa_Psi2 * (PHI[nTheta][2] + PHI[nTheta][j - 1]) +
             0.50 * dd * kappa_Psi1 * (PHI[nTheta][2] - PHI[nTheta][j - 1]) -
             dd2 * RHO[nTheta][j]) /
            (2.00 * (kappa_Theta2 + kappa_Psi2));
      } /* endif */

      PhiPi = PhiPi + PHI[nTheta][j];
    } /* endfor */

    Phi0 = Phi0 / (double)nPsi;
    PhiPi = PhiPi / (double)nPsi;

    for (j = 1; j <= nPsi; j++) {
      PHI[1][j] = Phi0;
      PHI[nTheta][j] = PhiPi;
    } /* endfor */
    break;

  case 3:
    /* Average update. */
    Phi0 = 0.00;
    PhiPi = 0.00;

    for (j = 1; j <= nPsi; j++) {
      Phi0 = Phi0 + PHI[2][j];
      PhiPi = PhiPi + PHI[nTheta - 1][j];
    } /* endfor */

    Phi0 = Phi0 / (double)nPsi;
    PhiPi = PhiPi / (double)nPsi;

    for (j = 1; j <= nPsi; j++) {
      PHI[1][j] = Phi0;
      PHI[nTheta][j] = PhiPi;
    } /* endfor */
    break;

  default:
  } /* endswitch */
}
