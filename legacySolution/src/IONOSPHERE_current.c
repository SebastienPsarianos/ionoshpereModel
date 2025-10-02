/* Determine Ionospheric Electric Field, Current Systems,
   and Convection Velocity for Ionospheric Solution */

/* This subroutine computes the ionospheric electric field,
   current, and convection velocity at each location of the 
   discretized solution domain given the height-integrated
   conductivities, electric field potential, and magnetic
   field.

   Parameters:

   ETh,EPs        Electric field components (spherical)
   Ex, Ey, Ez     Electric field components (Cartesian)
   Jx, Jy, Jz     Ionospheric current components (Cartesian)
   Ux, Uy, Uz     Ionospheric convection velocity components (Cartesian)
   PHI            Potential solution array
   SigmaThTh, SigmaThPs, SigmaPsPs
                  Conductances (height-integrated conductivities)
   Bx, By, Bz     Magnetic field components (Cartesian)
   Theta          Discrete values of Theta coordinate for fine grid
   Psi            Discrete values of Phi coordinate for fine grid
   X,Y,Z          Cartesian coordinates for fine grid
   Radius         Radius of sphere, Rp
   nTheta         Number of grid points in Theta-direction
                  on fine grid
   nPsi           Number of grid points in Psi-direction
                  on fine grid */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_current(ETh, EPs, Ex, Ey, Ez, Jx, Jy, Jz,
                        Ux, Uy, Uz,
                        PHI, SigmaThTh, SigmaThPs, SigmaPsPs,
                        Bx, By, Bz,
			Theta, Psi, X, Y, Z,
                        Radius, nTheta, nPsi)
double **ETh, **EPs, **Ex, **Ey, **Ez,
       **Jx, **Jy, **Jz, **Ux, **Uy, **Uz,
       **PHI, **SigmaThTh, **SigmaThPs, **SigmaPsPs,
       **Bx, **By, **Bz,
       **Theta, **Psi, **X, **Y, **Z, Radius;
int nTheta, nPsi;
{
        int i,j, nhemi;
        double dTheta, dPsi, dTheta2, dPsi2;
        double dd, dd2;
	double cosTheta, sinTheta, cosPhi, sinPhi, 
               JR, JTh, JPs, ER, 
	       xx, yy, zz, RR,
	       bx, by, bz, BB,
	       Vll, Vp_x, Vp_y, Vp_z, VR;

        dTheta=(Theta[nTheta][1]-Theta[1][1])/(double)(nTheta-1);
        dPsi=(Psi[1][nPsi]-Psi[1][1])/(double)(nPsi-1);
        dTheta2=dTheta*dTheta;
        dPsi2=dPsi*dPsi;
        dd=dTheta;
        dd2=dTheta2;

        for (j=1; j<=nPsi; j++) {
          if (j > 1 && j < nPsi ) {

            for (i=2; i<nTheta; i++) {
              sinTheta = sin(Theta[i][j]);
              ETh[i][j]=-(PHI[i+1][j]-PHI[i-1][j])/(2.00*dd*Radius);
              EPs[i][j]=-(PHI[i][j+1]-PHI[i][j-1])/(2.00*dd*Radius*sinTheta);
            } /* endfor */
            ETh[1][j]=-(PHI[2][j]-PHI[1][j])/(dd*Radius);
            EPs[1][j]=EPs[2][j];
            ETh[nTheta][j]=-(PHI[nTheta][j]-PHI[nTheta-1][j])/(dd*Radius);
            EPs[nTheta][j]=EPs[nTheta-1][j];

          } else if (j == 1) {

            for (i=2; i<nTheta; i++) {
              sinTheta  = sin(Theta[i][j]);
              ETh[i][j]=-(PHI[i+1][j]-PHI[i-1][j])/(2.00*dd*Radius);
              EPs[i][j]=-(PHI[i][j+1]-PHI[i][nPsi-1])/(2.00*dd*Radius*sinTheta);
             } /* endfor */
            ETh[1][j]=-(PHI[2][j]-PHI[1][j])/(dd*Radius);
            EPs[1][j]=EPs[2][j];
            ETh[nTheta][j]=-(PHI[nTheta][j]-PHI[nTheta-1][j])/(dd*Radius);
            EPs[nTheta][j]=EPs[nTheta-1][j];

          } else {

            for (i=2; i<nTheta; i++) {
              sinTheta = sin(Theta[i][j]);
              ETh[i][j]=-(PHI[i+1][j]-PHI[i-1][j])/(2.00*dd*Radius);
              EPs[i][j]=-(PHI[i][2]-PHI[i][j-1])/(2.00*dd*Radius*sinTheta);
            } /* endfor */
            ETh[1][j]=-(PHI[2][j]-PHI[1][j])/(dd*Radius);
            EPs[1][j]=EPs[2][j];
            ETh[nTheta][j]=-(PHI[nTheta][j]-PHI[nTheta-1][j])/(dd*Radius);
            EPs[nTheta][j]=EPs[nTheta-1][j];

          } /* endif */
        } /* endfor */

        nhemi = (nTheta-1)/2 + 1;

        for (j=1; j<=nPsi; j++) {
          for (i=1; i<=nTheta; i++) {

            cosTheta = cos(Theta[i][j]);
            sinTheta = sin(Theta[i][j]);
            cosPhi = cos(Psi[i][j]);
            sinPhi = sin(Psi[i][j]);

 	    RR = sqrt(X[i][j]*X[i][j]+
	              Y[i][j]*Y[i][j]+
	              Z[i][j]*Z[i][j]);
	    xx = X[i][j]/RR;
	    yy = Y[i][j]/RR;
	    zz = Z[i][j]/RR;

            ER=0.00;
            Ex[i][j] = ER*sinTheta*cosPhi + ETh[i][j]*cosTheta*cosPhi -
                       EPs[i][j]*sinPhi;
            Ey[i][j] = ER*sinTheta*sinPhi + ETh[i][j]*cosTheta*sinPhi +
                       EPs[i][j]*cosPhi;
            Ez[i][j] = ER*cosTheta - ETh[i][j]*sinTheta;

            JR=0.00;
            JTh=SigmaThTh[i][j]*ETh[i][j]+SigmaThPs[i][j]*EPs[i][j];
            JPs=-SigmaThPs[i][j]*ETh[i][j]+SigmaPsPs[i][j]*EPs[i][j];

            Jx[i][j] = JR*sinTheta*cosPhi + JTh*cosTheta*cosPhi -
                       JPs*sinPhi;
            Jy[i][j] = JR*sinTheta*sinPhi + JTh*cosTheta*sinPhi +
                       JPs*cosPhi;
            Jz[i][j] = JR*cosTheta - JTh*sinTheta;

            if (i != nhemi) {
     	       BB = sqrt(Bx[i][j]*Bx[i][j]+
		         By[i][j]*By[i][j]+
		         Bz[i][j]*Bz[i][j]);
	       bx = Bx[i][j]/BB;
	       by = By[i][j]/BB;
	       bz = Bz[i][j]/BB;

               Vp_x = (Ey[i][j]*bz - Ez[i][j]*by)/BB;
	       Vp_y = (Ez[i][j]*bx - Ex[i][j]*bz)/BB;
               Vp_z = (Ex[i][j]*by - Ey[i][j]*bx)/BB;

	       VR = Vp_x*xx + Vp_y*yy + Vp_z*zz;
	       Vll = -VR/(bx*xx+by*yy+bz*zz);

	       Ux[i][j] = Vp_x + Vll*bx;
	       Uy[i][j] = Vp_y + Vll*by;
	       Uz[i][j] = Vp_z + Vll*bz;
            }

          } /* endfor */
        } /* endfor */

        for (j=1; j<=nPsi; j++) {
          i=nhemi;
	  Ux[i][j] = 0.50*(Ux[i+1][j]+Ux[i-1][j]);
	  Uy[i][j] = 0.50*(Uy[i+1][j]+Uy[i-1][j]);
	  Uz[i][j] = 0.50*(Uz[i+1][j]+Uz[i-1][j]);
        } /* endfor */

}

