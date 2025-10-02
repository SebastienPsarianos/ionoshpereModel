/* Bilinear Interpolation Procedure for
   Magnetospheric Convection Velocities */

/* This subroutine determines the components
   of the magnetospheric convection velocities
   as determined by the ionospheric solution for
   the convection velocities (based on electric field
   and current systems and enforcement of zero
   normal velocity).  Bilinear interpolation is
   used to 

   Parameters:

   PHI        Potential solution array
   Ux,Uy,Uz   Convection velocity components
              (ionosphere calculation)
   Theta      Discrete values of Theta coordinate
              (ionosphere calculation)
   Psi        Discrete values of Phi coordinate
              (ionosphere calculation)
   X,Y,Z      Cartesian coordinates (x,y,z)  
              (ionosphere calculation)
   Thetamag   Discrete values of Theta coordinate
              (magnetosphere calculation)
   Psimag     Discrete values of Phi coordinate
              (magnetosphere calculation)
   Xmag,Ymag, Cartesian coordinates (x,y,z) where
   Zmag,      FAC is defined (magnetosphere calculation)
   Umag_x,    Components of the velocity field
   Umag_y,    (magnetosphere calculation)
   Umag_z
   Radius     Radius of sphere, Rp
   nTheta     Number of grid points in Theta-direction
              on finest grid
   nPsi       Number of grid points in Psi-direction
              on finest grid
   nmag       Number of magnetospheric solution points
              on ionospheric surface
   InMagTri   Pointer to the triangular element of the
              magnetospheric solution in which the
	      ionospheric grid point is found
   TriMagNode1,TriMagNode2,TriMagNode3
              Nodes of triangular elements used to
	      defined magnetospheric solution
   ntrimag    Number of triangular elements used in
              defining the magnetospheric solution
	      on the ionospheric surface */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_magvel(PHI, Ux, Uy, Uz,
                       Theta, Psi, X, Y, Z, Radius,
                       nTheta, nPsi,
                       Umag_x, Umag_y, Umag_z,
                       Thetamag, Psimag,
                       Xmag, Ymag, Zmag, nmag,
		       InMagTri, MagTriNode1, MagTriNode2,
		       MagTriNode3, ntrimag)
double **PHI, **Ux, **Uy, **Uz,
       **Theta, **Psi, **X, **Y, **Z, Radius;
double *Thetamag, *Psimag,
       *Xmag, *Ymag, *Zmag,
       *Umag_x, *Umag_y, *Umag_z;
int **InMagTri, *MagTriNode1, *MagTriNode2, *MagTriNode3;
int nTheta, nPsi, nmag, ntrimag;
{
        int n,i,j;
        double dTheta, dPsi;
	double x,y,x1,x2,x3,y1,y2,y3,x4,y4,
	       f1,f2,f3,f4,a,b,c,d;
	
        dPsi=2.00*IONOSPHERE_PI/(double)(nPsi-1);
        dTheta=dPsi;

	for (n = 1; n <= nmag; ++n) {

	   /* Magnetospheric point. */

	   x=Thetamag[n];
	   if (x <= 0.50*IONOSPHERE_PI) {
	     x = x + IONOSPHERE_Theta_0;
	   } else {
             x = x - IONOSPHERE_Theta_0;
	   }
	   y=Psimag[n];
	   
	   i=x/dTheta+1;
	   j=y/dPsi+1;

	   /* Four corners of ionospheric cell. */

           x1=Theta[i][j];
	   y1=Psi[i][j];
           x2=Theta[i][j+1];
	   y2=Psi[i][j+1];
           x3=Theta[i+1][j+1];
	   y3=Psi[i+1][j+1];
           x4=Theta[i+1][j];
	   y4=Psi[i+1][j];

	   /* Interpolate each velocity component. */

           f1=Ux[i][j];
           f2=Ux[i][j+1];
           f3=Ux[i+1][j+1];
           f4=Ux[i+1][j];
           a=f1;
	   b=(f4-f1)/dTheta;
	   c=(f2-f1)/dPsi;
	   d=(f3+f1-f4-f2)/(dTheta*dPsi);
	   Umag_x[n] = a+b*(x-x1)+c*(y-y1)+d*(x-x1)*(y-y1);

           f1=Uy[i][j];
           f2=Uy[i][j+1];
           f3=Uy[i+1][j+1];
           f4=Uy[i+1][j];
           a=f1;
	   b=(f4-f1)/dTheta;
	   c=(f2-f1)/dPsi;
	   d=(f3+f1-f4-f2)/(dTheta*dPsi);
	   Umag_y[n] = a+b*(x-x1)+c*(y-y1)+d*(x-x1)*(y-y1);

           f1=Uz[i][j];
           f2=Uz[i][j+1];
           f3=Uz[i+1][j+1];
           f4=Uz[i+1][j];
           a=f1;
	   b=(f4-f1)/dTheta;
	   c=(f2-f1)/dPsi;
	   d=(f3+f1-f4-f2)/(dTheta*dPsi);
	   Umag_z[n] = a+b*(x-x1)+c*(y-y1)+d*(x-x1)*(y-y1);

	   /* Non-dimensionalize velocity components. */

	   Umag_x[n] = Umag_x[n]/IONOSPHERE_SoundSpeed_Earth;
	   Umag_y[n] = Umag_y[n]/IONOSPHERE_SoundSpeed_Earth;
	   Umag_z[n] = Umag_z[n]/IONOSPHERE_SoundSpeed_Earth;

        } /* endfor */

}
