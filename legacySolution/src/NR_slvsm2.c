/* Solution of Model Problem on Coarsest Grid */

/* This subroutine provides an exact solution to the
   model problem on the coarsest grid where h=1/2 and
   there are only 3 grid points in the x- and y-directions.
   The right-hand-side array is input in rhs(1:3,1:3) and
   the solution is returned in u(1:3,1:3).
   
   Parameters:

   u      Solution array on coarsest grid
   rhs    rhs array

   Additional routines:

   NR_fill0 */

#include <math.h>

void NR_slvsm2(u,rhs)
double **rhs,**u;
{
	void NR_fill0();
	double disc,fact,h=0.5;

	NR_fill0(u,3);
	fact=2.0/(h*h);
	disc=sqrt(fact*fact+rhs[2][2]);
	u[2][2] = -rhs[2][2]/(fact+disc);
}
