/* Bilinear Interpolation Procedure for
   Magnetospheric Magnetic Field */

/* This subroutine obtains the components of the magnetic
   field of the magnetospheric solution used in ionosphere
   model.  Bilinear interpolation is used to map the 
   various components from the magnetosphere calculation to
   the two-dimensional uniform mesh in the co-latitude and
   longitude coordinate system (Theta, Psi) used in the
   ionosphere model.

   NOTE: THIS INTERPOLATION PROCEDURE ASSUMES 
         SYMMETRY OF MAGNETOSPHERE GRID ABOUT
         X=0, Y=0, AND Z=0 PLANES.

   Parameters:

   Bx,By,Bz   Components of the magnetic field
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
   Bmag_x,    Components of the magnetic field
   Bmag_y,    (magnetosphere calculation)
   Bmag_z
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
	      on the ionospheric surface
   ireinit    Integer parameter indicating whether
              or not to re-initialize the pointers to
	      the triangular elements of the magnetospheric
	      solution */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_magbfield(Bx, By, Bz, Theta, Psi, X, Y, Z, Radius,
                          nTheta, nPsi,
                          Bmag_x, Bmag_y, Bmag_z,
                          Thetamag, Psimag,
                          Xmag, Ymag, Zmag, nmag,
		          InMagTri, MagTriNode1, MagTriNode2,
		          MagTriNode3, ntrimag, ireinit)
double **Bx, **By, **Bz, **Theta, **Psi, **X, **Y, **Z, Radius;
double *Thetamag, *Psimag,
       *Xmag, *Ymag, *Zmag,
       *Bmag_x, *Bmag_y, *Bmag_z;
int **InMagTri, *MagTriNode1, *MagTriNode2, *MagTriNode3;
int nTheta, nPsi, nmag, ntrimag, ireinit;
{
        int i,j,k,n,itrimag;
	double x1,x2,x3,y1,y2,y3,det0,det1,det2,det3,
	       f1,f2,f3,a,b,c,x,y;
	double z1,z2,z3,z;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	
        for (j = 1; j <= nPsi; j++) {
          for (i = 1; i <= nTheta; i++) {

            x=X[i][j];
	    y=Y[i][j];
	    z=Z[i][j];

	    if (ireinit == 0  && 
                (Theta[i][j] <= 0.45*IONOSPHERE_PI ||
                 Theta[i][j] >= 0.55*IONOSPHERE_PI) ) {

               itrimag=InMagTri[i][j];
	       
               x1=Xmag[MagTriNode1[itrimag]];
	       y1=Ymag[MagTriNode1[itrimag]];
	       x2=Xmag[MagTriNode2[itrimag]];
	       y2=Ymag[MagTriNode2[itrimag]];
	       x3=Xmag[MagTriNode3[itrimag]];
	       y3=Ymag[MagTriNode3[itrimag]];

  	       f1=Bmag_x[MagTriNode1[itrimag]];
	       f2=Bmag_x[MagTriNode2[itrimag]];
	       f3=Bmag_x[MagTriNode3[itrimag]];

	       det0 = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
               a =   (x2*y3-y2*x3)*f1/det0
	           + (y1*x3-x1*y3)*f2/det0
	           + (x1*y2-y1*x2)*f3/det0;
               b =   (y2-y3)*f1/det0
	           + (y3-y1)*f2/det0
                   + (y1-y2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bx[i][j]=a+b*x+c*y;

  	       f1=Bmag_y[MagTriNode1[itrimag]];
	       f2=Bmag_y[MagTriNode2[itrimag]];
	       f3=Bmag_y[MagTriNode3[itrimag]];

	       det0 = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
               a =   (x2*y3-y2*x3)*f1/det0
	           + (y1*x3-x1*y3)*f2/det0
	           + (x1*y2-y1*x2)*f3/det0;
               b =   (y2-y3)*f1/det0
	           + (y3-y1)*f2/det0
                   + (y1-y2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               By[i][j]=a+b*x+c*y;

  	       f1=Bmag_z[MagTriNode1[itrimag]];
	       f2=Bmag_z[MagTriNode2[itrimag]];
	       f3=Bmag_z[MagTriNode3[itrimag]];

	       det0 = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
               a =   (x2*y3-y2*x3)*f1/det0
	           + (y1*x3-x1*y3)*f2/det0
	           + (x1*y2-y1*x2)*f3/det0;
               b =   (y2-y3)*f1/det0
	           + (y3-y1)*f2/det0
                   + (y1-y2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bz[i][j]=a+b*x+c*y;

            } else if (ireinit == 0 && 
                       ((Psi[i][j] >= 0.05*IONOSPHERE_PI &&
                         Psi[i][j] <= 0.95*IONOSPHERE_PI) ||
                        (Psi[i][j] >= 1.05*IONOSPHERE_PI &&
                         Psi[i][j] <= 1.95*IONOSPHERE_PI) ) ) {

               itrimag=InMagTri[i][j];
	       
               x1=Xmag[MagTriNode1[itrimag]];
	       z1=Zmag[MagTriNode1[itrimag]];
	       x2=Xmag[MagTriNode2[itrimag]];
	       z2=Zmag[MagTriNode2[itrimag]];
	       x3=Xmag[MagTriNode3[itrimag]];
	       z3=Zmag[MagTriNode3[itrimag]];

  	       f1=Bmag_x[MagTriNode1[itrimag]];
	       f2=Bmag_x[MagTriNode2[itrimag]];
	       f3=Bmag_x[MagTriNode3[itrimag]];

	       det0 = x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2;
               a =   (x2*z3-z2*x3)*f1/det0
	           + (z1*x3-x1*z3)*f2/det0
	           + (x1*z2-z1*x2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bx[i][j]=a+b*x+c*z;

  	       f1=Bmag_y[MagTriNode1[itrimag]];
	       f2=Bmag_y[MagTriNode2[itrimag]];
	       f3=Bmag_y[MagTriNode3[itrimag]];

	       det0 = x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2;
               a =   (x2*z3-z2*x3)*f1/det0
	           + (z1*x3-x1*z3)*f2/det0
	           + (x1*z2-z1*x2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               By[i][j]=a+b*x+c*z;

  	       f1=Bmag_z[MagTriNode1[itrimag]];
	       f2=Bmag_z[MagTriNode2[itrimag]];
	       f3=Bmag_z[MagTriNode3[itrimag]];

	       det0 = x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2;
               a =   (x2*z3-z2*x3)*f1/det0
	           + (z1*x3-x1*z3)*f2/det0
	           + (x1*z2-z1*x2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bz[i][j]=a+b*x+c*z;

            } else if (ireinit == 0) {

               itrimag=InMagTri[i][j];
	       
               y1=Ymag[MagTriNode1[itrimag]];
	       z1=Zmag[MagTriNode1[itrimag]];
	       y2=Ymag[MagTriNode2[itrimag]];
	       z2=Zmag[MagTriNode2[itrimag]];
	       y3=Ymag[MagTriNode3[itrimag]];
	       z3=Zmag[MagTriNode3[itrimag]];

  	       f1=Bmag_x[MagTriNode1[itrimag]];
	       f2=Bmag_x[MagTriNode2[itrimag]];
	       f3=Bmag_x[MagTriNode3[itrimag]];

	       det0 = y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2;
               a =   (y2*z3-z2*y3)*f1/det0
	           + (z1*y3-y1*z3)*f2/det0
	           + (y1*z2-z1*y2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (y3-y2)*f1/det0
                   + (y1-y3)*f2/det0
                   + (y2-y1)*f3/det0;
               Bx[i][j]=a+b*y+c*z;

  	       f1=Bmag_y[MagTriNode1[itrimag]];
	       f2=Bmag_y[MagTriNode2[itrimag]];
	       f3=Bmag_y[MagTriNode3[itrimag]];

	       det0 = y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2;
               a =   (y2*z3-z2*y3)*f1/det0
	           + (z1*y3-y1*z3)*f2/det0
	           + (y1*z2-z1*y2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (y3-y2)*f1/det0
                   + (y1-y3)*f2/det0
                   + (y2-y1)*f3/det0;
               By[i][j]=a+b*y+c*z;

  	       f1=Bmag_z[MagTriNode1[itrimag]];
	       f2=Bmag_z[MagTriNode2[itrimag]];
	       f3=Bmag_z[MagTriNode3[itrimag]];

	       det0 = y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2;
               a =   (y2*z3-z2*y3)*f1/det0
	           + (z1*y3-y1*z3)*f2/det0
	           + (y1*z2-z1*y2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (y3-y2)*f1/det0
                   + (y1-y3)*f2/det0
                   + (y2-y1)*f3/det0;
               Bz[i][j]=a+b*y+c*z;

	    } else if (Theta[i][j] <= 0.45*IONOSPHERE_PI ||
                       Theta[i][j] >= 0.55*IONOSPHERE_PI) {
	      
               for (n = 1, itrimag=0; n <= ntrimag && itrimag < 1; ++n) {

                 x1=Xmag[MagTriNode1[n]];
	         y1=Ymag[MagTriNode1[n]];
	         z1=Zmag[MagTriNode1[n]];
	         x2=Xmag[MagTriNode2[n]];
	         y2=Ymag[MagTriNode2[n]];
	         z2=Zmag[MagTriNode2[n]];
	         x3=Xmag[MagTriNode3[n]];
	         y3=Ymag[MagTriNode3[n]];
       	         z3=Zmag[MagTriNode3[n]];

	         xmin=IONOSPHERE_DMIN(x1,x2);
	         xmin=IONOSPHERE_DMIN(xmin,x3);
	         ymin=IONOSPHERE_DMIN(y1,y2);
	         ymin=IONOSPHERE_DMIN(ymin,y3);

	         xmax=IONOSPHERE_DMAX(x1,x2);
	         xmax=IONOSPHERE_DMAX(xmax,x3);
	         ymax=IONOSPHERE_DMAX(y1,y2);
	         ymax=IONOSPHERE_DMAX(ymax,y3);

	         if (x >= xmin && x <= xmax &&
		     y >= ymin && y <= ymax &&
		     z*z1 >= 0 && z*z2 >= 0 &&
		     z*z3 >= 0) {
		   det1=(x-x1)*(y2-y1)-(y-y1)*(x2-x1);
		   det2=(x-x2)*(y3-y2)-(y-y2)*(x3-x2);
		   det3=(x-x3)*(y1-y3)-(y-y3)*(x1-x3);
		   if (det1 <= 0) {
		     det1=-det1;
		     det2=-det2;
		     det3=-det3;
   		   } /* endif */
	  	   if (det1 >= 0 && det2 >=0 && det3 >= 0)
	             itrimag = n;
	         } else if ( n == ntrimag) {
		   itrimag = ntrimag + 1;
	         } /* endif */
	      
               } /* endfor */

               InMagTri[i][j]=itrimag;

  	       f1=Bmag_x[MagTriNode1[itrimag]];
	       f2=Bmag_x[MagTriNode2[itrimag]];
	       f3=Bmag_x[MagTriNode3[itrimag]];

	       det0 = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
               a =   (x2*y3-y2*x3)*f1/det0
	           + (y1*x3-x1*y3)*f2/det0
	           + (x1*y2-y1*x2)*f3/det0;
               b =   (y2-y3)*f1/det0
	           + (y3-y1)*f2/det0
                   + (y1-y2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bx[i][j]=a+b*x+c*y;

  	       f1=Bmag_y[MagTriNode1[itrimag]];
	       f2=Bmag_y[MagTriNode2[itrimag]];
	       f3=Bmag_y[MagTriNode3[itrimag]];

	       det0 = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
               a =   (x2*y3-y2*x3)*f1/det0
	           + (y1*x3-x1*y3)*f2/det0
	           + (x1*y2-y1*x2)*f3/det0;
               b =   (y2-y3)*f1/det0
	           + (y3-y1)*f2/det0
                   + (y1-y2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               By[i][j]=a+b*x+c*y;

  	       f1=Bmag_z[MagTriNode1[itrimag]];
	       f2=Bmag_z[MagTriNode2[itrimag]];
	       f3=Bmag_z[MagTriNode3[itrimag]];

	       det0 = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
               a =   (x2*y3-y2*x3)*f1/det0
	           + (y1*x3-x1*y3)*f2/det0
	           + (x1*y2-y1*x2)*f3/det0;
               b =   (y2-y3)*f1/det0
	           + (y3-y1)*f2/det0
                   + (y1-y2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bz[i][j]=a+b*x+c*y;

	    } else if ((Psi[i][j] >= 0.05*IONOSPHERE_PI &&
                        Psi[i][j] <= 0.95*IONOSPHERE_PI) ||
                       (Psi[i][j] >= 1.05*IONOSPHERE_PI &&
                        Psi[i][j] <= 1.95*IONOSPHERE_PI) ) {
	      
               for (n = 1, itrimag=0; n <= ntrimag && itrimag < 1; ++n) {

                 x1=Xmag[MagTriNode1[n]];
	         y1=Ymag[MagTriNode1[n]];
	         z1=Zmag[MagTriNode1[n]];
	         x2=Xmag[MagTriNode2[n]];
	         y2=Ymag[MagTriNode2[n]];
	         z2=Zmag[MagTriNode2[n]];
	         x3=Xmag[MagTriNode3[n]];
	         y3=Ymag[MagTriNode3[n]];
       	         z3=Zmag[MagTriNode3[n]];

	         xmin=IONOSPHERE_DMIN(x1,x2);
	         xmin=IONOSPHERE_DMIN(xmin,x3);
	         zmin=IONOSPHERE_DMIN(z1,z2);
	         zmin=IONOSPHERE_DMIN(zmin,z3);

	         xmax=IONOSPHERE_DMAX(x1,x2);
	         xmax=IONOSPHERE_DMAX(xmax,x3);
	         zmax=IONOSPHERE_DMAX(z1,z2);
	         zmax=IONOSPHERE_DMAX(zmax,z3);

	         if (x >= xmin && x <= xmax &&
		     z >= zmin && z <= zmax &&
		     y*y1 >= 0 && y*y2 >= 0 &&
		     y*y3 >= 0) {
		   det1=(x-x1)*(z2-z1)-(z-z1)*(x2-x1);
		   det2=(x-x2)*(z3-z2)-(z-z2)*(x3-x2);
		   det3=(x-x3)*(z1-z3)-(z-z3)*(x1-x3);
		   if (det1 <= 0) {
		     det1=-det1;
		     det2=-det2;
		     det3=-det3;
   		   } /* endif */
	  	   if (det1 >= 0 && det2 >=0 && det3 >= 0)
	             itrimag = n;
	         } else if ( n == ntrimag) {
		   itrimag = ntrimag + 1;
	         } /* endif */
	      
               } /* endfor */

               InMagTri[i][j]=itrimag;

  	       f1=Bmag_x[MagTriNode1[itrimag]];
	       f2=Bmag_x[MagTriNode2[itrimag]];
	       f3=Bmag_x[MagTriNode3[itrimag]];

	       det0 = x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2;
               a =   (x2*z3-z2*x3)*f1/det0
	           + (z1*x3-x1*z3)*f2/det0
	           + (x1*z2-z1*x2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bx[i][j]=a+b*x+c*z;

  	       f1=Bmag_y[MagTriNode1[itrimag]];
	       f2=Bmag_y[MagTriNode2[itrimag]];
	       f3=Bmag_y[MagTriNode3[itrimag]];

	       det0 = x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2;
               a =   (x2*z3-z2*x3)*f1/det0
	           + (z1*x3-x1*z3)*f2/det0
	           + (x1*z2-z1*x2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               By[i][j]=a+b*x+c*z;

  	       f1=Bmag_z[MagTriNode1[itrimag]];
	       f2=Bmag_z[MagTriNode2[itrimag]];
	       f3=Bmag_z[MagTriNode3[itrimag]];

	       det0 = x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2;
               a =   (x2*z3-z2*x3)*f1/det0
	           + (z1*x3-x1*z3)*f2/det0
	           + (x1*z2-z1*x2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (x3-x2)*f1/det0
                   + (x1-x3)*f2/det0
                   + (x2-x1)*f3/det0;
               Bz[i][j]=a+b*x+c*z;

	    } else {
	      
               for (n = 1, itrimag=0; n <= ntrimag && itrimag < 1; ++n) {

                 x1=Xmag[MagTriNode1[n]];
	         y1=Ymag[MagTriNode1[n]];
	         z1=Zmag[MagTriNode1[n]];
	         x2=Xmag[MagTriNode2[n]];
	         y2=Ymag[MagTriNode2[n]];
	         z2=Zmag[MagTriNode2[n]];
	         x3=Xmag[MagTriNode3[n]];
	         y3=Ymag[MagTriNode3[n]];
       	         z3=Zmag[MagTriNode3[n]];

	         ymin=IONOSPHERE_DMIN(y1,y2);
	         ymin=IONOSPHERE_DMIN(ymin,y3);
	         zmin=IONOSPHERE_DMIN(z1,z2);
	         zmin=IONOSPHERE_DMIN(zmin,z3);

	         ymax=IONOSPHERE_DMAX(y1,y2);
	         ymax=IONOSPHERE_DMAX(ymax,y3);
	         zmax=IONOSPHERE_DMAX(z1,z2);
	         zmax=IONOSPHERE_DMAX(zmax,z3);

	         if (y >= ymin && y <= ymax &&
		     z >= zmin && z <= zmax &&
		     x*x1 >= 0 && x*x2 >= 0 &&
		     x*x3 >= 0) {
		   det1=(y-y1)*(z2-z1)-(z-z1)*(y2-y1);
		   det2=(y-y2)*(z3-z2)-(z-z2)*(y3-y2);
		   det3=(y-y3)*(z1-z3)-(z-z3)*(y1-y3);
		   if (det1 <= 0) {
		     det1=-det1;
		     det2=-det2;
		     det3=-det3;
   		   } /* endif */
	  	   if (det1 >= 0 && det2 >=0 && det3 >= 0)
	             itrimag = n;
	         } else if ( n == ntrimag) {
		   itrimag = ntrimag + 1;
	         } /* endif */
	      
               } /* endfor */

               InMagTri[i][j]=itrimag;

  	       f1=Bmag_x[MagTriNode1[itrimag]];
	       f2=Bmag_x[MagTriNode2[itrimag]];
	       f3=Bmag_x[MagTriNode3[itrimag]];

	       det0 = y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2;
               a =   (y2*z3-z2*y3)*f1/det0
	           + (z1*y3-y1*z3)*f2/det0
	           + (y1*z2-z1*y2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (y3-y2)*f1/det0
                   + (y1-y3)*f2/det0
                   + (y2-y1)*f3/det0;
               Bx[i][j]=a+b*y+c*z;

  	       f1=Bmag_y[MagTriNode1[itrimag]];
	       f2=Bmag_y[MagTriNode2[itrimag]];
	       f3=Bmag_y[MagTriNode3[itrimag]];

	       det0 = y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2;
               a =   (y2*z3-z2*y3)*f1/det0
	           + (z1*y3-y1*z3)*f2/det0
	           + (y1*z2-z1*y2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (y3-y2)*f1/det0
                   + (y1-y3)*f2/det0
                   + (y2-y1)*f3/det0;
               By[i][j]=a+b*y+c*z;

  	       f1=Bmag_z[MagTriNode1[itrimag]];
	       f2=Bmag_z[MagTriNode2[itrimag]];
	       f3=Bmag_z[MagTriNode3[itrimag]];

	       det0 = y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2;
               a =   (y2*z3-z2*y3)*f1/det0
	           + (z1*y3-y1*z3)*f2/det0
	           + (y1*z2-z1*y2)*f3/det0;
               b =   (z2-z3)*f1/det0
	           + (z3-z1)*f2/det0
                   + (z1-z2)*f3/det0;
               c =   (y3-y2)*f1/det0
                   + (y1-y3)*f2/det0
                   + (y2-y1)*f3/det0;
               Bz[i][j]=a+b*y+c*z;
	       
	    } /* endif */
	    
          } /* endfor */

        } /* endfor */

}
