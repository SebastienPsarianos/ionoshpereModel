/* Full Linear Multigrid Solution Algorithm */

/* This subroutine applies a full linear multigrid solution
   algorithm to the problem of determining the solution of the
   boundary value problem (BVP) for the linear elliptic equation

          d^2 u      d^2 u
	  -----  +   -----   =   r   ,
	  dx^2       dx^2

   defined on the square domain 0 < x < 1, 0 < y < 1, and
   subject to the Dirichlet boundary data

          u(0,y) = u(1,y) = u(x,0) = u(x,1) = 0.

   A red-black Gauss-Seidel scheme is used as the smoothing
   operator, bilinear interpolation is used for the prolongation
   operator, and half-weighting is used for the restriction
   operator.
	  
   On input the two-diminsional solution array u(1:n,1:n)
   contains the right-hand-side source term, r, while on
   output it returns the solution to the elliptic PDE.  The
   number of grid points n should be related to the number of grid
   levels used in the solution, ng, by n = 2^ng+1.

   Parameters:

   u      Solution array
   n      Number of grid points in x- and y-directions on
          finest grid
   ncycle Number of V-cyles to be used at each level
   ng     Number of grid levels used
   NPRE   Number of relaxation sweeps before coarse-grid
          correction is computed
   NPOST  Number of relaxation sweeps
	
   Additional routines:

   NR_addint, NR_copy, NR_fill0, NR_interp, NR_relax, NR_resid,
   NR_rstrct, NR_slvsml. */

#include "NR_util.h"
#define NPRE 1
#define NPOST 1
#define NGMAX 15

void NR_mglin(u,n,ncycle)
double **u;
int n,ncycle;
{
	void NR_addint(),NR_copy(),NR_fill0(),NR_interp(),NR_relax(),
	     NR_resid(),NR_rstrct(),NR_slvsml();
	unsigned int j,jcycle,jj,jpost,jpre,nf,ng=0,ngrid,nn;
	double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];

	nn=n;
	while (nn >>= 1) ng++;
	if (n != 1+(1L << ng)) NR_error("n-1 must be a power of 2 in mglin.");
	if (ng > NGMAX) NR_error("increase NGMAX in mglin.");
	
        /* Allocate storage for rhs on grid ng-1 and fill it by
	   restricting rhs from the fine grid. */
	
	nn=n/2+1;
	ngrid=ng-1;
	irho[ngrid]=NR_dmatrix(1,nn,1,nn);
	NR_rstrct(irho[ngrid],u,nn);

	/* Similarily allocate storage and fill rhs on all
	   coarse grids. */
	
	while (nn > 3) {
		nn=nn/2+1;
		irho[--ngrid]=NR_dmatrix(1,nn,1,nn);
		NR_rstrct(irho[ngrid],irho[ngrid+1],nn);
	}
	nn=3;
	iu[1]=NR_dmatrix(1,nn,1,nn);
	irhs[1]=NR_dmatrix(1,nn,1,nn);
	
	/* Initial solution on coarsest grid. */
	
	NR_slvsml(iu[1],irho[1]);
	NR_free_dmatrix(irho[1],1,nn,1,nn);

	/* Nested iteration loop. */
	
	ngrid=ng;
	for (j=2;j<=ngrid;j++) {
		nn=2*nn-1;
		iu[j]=NR_dmatrix(1,nn,1,nn);
		irhs[j]=NR_dmatrix(1,nn,1,nn);
		ires[j]=NR_dmatrix(1,nn,1,nn);

		/* Interpolate from coarse grid to next
		   finer grid. */
		
		NR_interp(iu[j],iu[j-1],nn);

		/* Set up rhs. */
		
		NR_copy(irhs[j],(j != ngrid ? irho[j] : u),nn);

		/* V-cycle loop. */
		
		for (jcycle=1;jcycle<=ncycle;jcycle++) {
			nf=nn;
			for (jj=j;jj>=2;jj--) {

			  /* Pre-smoothing on downward stroke of
			     V-cycle */

		          for (jpre=1;jpre<=NPRE;jpre++)
		               NR_relax(iu[jj],irhs[jj],nf);

			  /* Restricted residual is the next rhs. */
			
			  NR_resid(ires[jj],iu[jj],irhs[jj],nf);
			  nf=nf/2+1;
			  NR_rstrct(irhs[jj-1],ires[jj],nf);

		    	  /* Zero for initial guess in next relaxation. */
			
			  NR_fill0(iu[jj-1],nf);
			}

			/* Bottom of V-cycle; solve on coarsest grid. */
			
			NR_slvsml(iu[1],irhs[1]);
			nf=3;
			for (jj=2;jj<=j;jj++) {
 			  nf=2*nf-1;
			  NR_addint(iu[jj],iu[jj-1],ires[jj],nf);

			  /* Post-smoothing on upward stroke of
			     V-cycle. */
			  
			  for (jpost=1;jpost<=NPOST;jpost++)
				  NR_relax(iu[jj],irhs[jj],nf);
			}
		}
	}

	/* Return solution in u. */
	
	NR_copy(u,iu[ngrid],n);

	/* De-allocate storage. */
	
	for (nn=n,j=ng;j>=2;j--,nn=nn/2+1) {
		NR_free_dmatrix(ires[j],1,nn,1,nn);
		NR_free_dmatrix(irhs[j],1,nn,1,nn);
		NR_free_dmatrix(iu[j],1,nn,1,nn);
		if (j != ng) NR_free_dmatrix(irho[j],1,nn,1,nn);
	}
	NR_free_dmatrix(irhs[1],1,3,1,3);
	NR_free_dmatrix(iu[1],1,3,1,3);
}
#undef NPRE
#undef NPOST
#undef NGMAX
