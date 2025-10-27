/* Dimensionalization Procedure for
   Magnetospheric Variables */

/* This subroutine dimensionalizes the
   appropriate magnetospheric variables.

   Parameters:

   Thetamag   Discrete values of Theta coordinate
              (magnetosphere calculation)
   Psimag     Discrete values of Phi coordinate
              (magnetosphere calculation)
   Xmag,Ymag, Cartesian coordinates (x,y,z) where
   Zmag,      FAC is defined (magnetosphere calculation)
   Jmag_x,    Components of the field-aligned current (FAC)
   Jmag_y,    (magnetosphere calculation)
   Jmag_z,
   Jmag_R
   Bmag_x,    Components of the magnetic field
   Bmag_y,    (magnetosphere calculation)
   Bmag_z
   Radius     Radius of sphere, Rp
   nmag       Number of magnetospheric solution points
              on ionospheric surface
   TriMagNode1,TriMagNode2,TriMagNode3
              Nodes of triangular elements used to
              defined magnetospheric solution
   ntrimag    Number of triangular elements used in
              defining the magnetospheric solution
              on the ionospheric surface
   ireinit    Integer parameter indicating whether
              or not to re-initialize the coordinate
              values for the nodes of the
              triangular elements of the magnetospheric
              solution */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_magdim(Thetamag, Psimag, Xmag, Ymag, Zmag, Jmag_x, Jmag_y,
                       Jmag_z, Jmag_R, Bmag_x, Bmag_y, Bmag_z, nmag,
                       MagTriNode1, MagTriNode2, MagTriNode3, ntrimag,
                       ireinit) double *Thetamag,
    *Psimag, *Xmag, *Ymag, *Zmag, *Jmag_x, *Jmag_y, *Jmag_z, *Jmag_R, *Bmag_x,
    *Bmag_y, *Bmag_z;
int *MagTriNode1, *MagTriNode2, *MagTriNode3;
int nmag, ntrimag, ireinit;
{
  int i;
  double Radius, ConstantB, ConstantJ;

  ConstantB = sqrt(IONOSPHERE_MU * IONOSPHERE_Density_Earth *
                   IONOSPHERE_SoundSpeed_Earth * IONOSPHERE_SoundSpeed_Earth);
  ConstantJ = ConstantB / (IONOSPHERE_MU * IONOSPHERE_Radius_Earth);

  for (i = 1; i <= nmag; ++i) {
    if (ireinit != 0) {
      Radius = sqrt(Xmag[i] * Xmag[i] + Ymag[i] * Ymag[i] + Zmag[i] * Zmag[i]);
      Xmag[i] = Xmag[i] / Radius;
      Ymag[i] = Ymag[i] / Radius;
      Zmag[i] = Zmag[i] / Radius;
      if (Ymag[i] == 0.00 && Xmag[i] == 0.00) {
        Psimag[i] = 0.00;
      } else {
        Psimag[i] = atan2(Ymag[i], Xmag[i]);
      } /* end if */
      if (Psimag[i] < 0.00)
        Psimag[i] = Psimag[i] + 2.00 * IONOSPHERE_PI;
      Radius = sqrt(Xmag[i] * Xmag[i] + Ymag[i] * Ymag[i]);
      Thetamag[i] = atan2(Radius, Zmag[i]);
    }

    Jmag_R[i] = Jmag_x[i] * sin(Thetamag[i]) * cos(Psimag[i]) +
                Jmag_y[i] * sin(Thetamag[i]) * sin(Psimag[i]) +
                Jmag_z[i] * cos(Thetamag[i]);
    Jmag_R[i] = Jmag_R[i] * ConstantJ;
    Bmag_x[i] = Bmag_x[i] * ConstantB;
    Bmag_y[i] = Bmag_y[i] * ConstantB;
    Bmag_z[i] = Bmag_z[i] * ConstantB;

  } /* endfor */
}
