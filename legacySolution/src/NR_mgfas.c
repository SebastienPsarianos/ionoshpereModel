/* Full Multigrid Full Approximation Storage (FAS) Solution Algorithm */

/* This subroutine applies a full multigrid FAS solution
   algorithm to the problem of determining the solution of the
   boundary value problem (BVP) for the non-linear elliptic equation

          d^2 u      d^2 u
	  -----  +   -----   + u^2  =  r  ,
	  dx^2       dx^2

   defined on the square domain 0 < x < 1, 0 < y < 1, and
   subject to the Dirichlet boundary data

          u(0,y) = u(1,y) = u(x,0) = u(x,1) = 0.

   A red-black Gauss-Seidel scheme with Newton acceleration
   is used as the smoothing operator, bilinear interpolation is
   used for the prolongation operator, and half-weighting is
   used for the restriction operator.  Both residuals (defect
   corrections) and solutions are restricted to the coarser
   grids.
	  
   On input the two-diminsional solution array u(1:n,1:n)
   contains the right-hand-side source term, r, while on
   output it returns the solution to the elliptic PDE.  The
   number of grid points n should be related to the number of grid
   levels used in the solution, ng, by n = 2^ng+1.

   Parameters:

   u      Solution array
   n      Number of grid points in x- and y-directions on
          finest grid
   maxcyc Number of V-cyles to be used at each level
   ng     Number of grid levels used
   NPRE   Number of relaxation sweeps before coarse-grid
          correction is computed
   NPOST  Number of relaxation sweeps
   ALPHA  Relates estimated truncation error to the norm of
          the residual.
	
   Additional routines:

   NR_anorm2, NR_copy, NR_interp, NR_lop, NR_matadd, NR_matsub,
   NR_relax2, NR_rstrct, NR_slvsml2. */

#include "NR_util.h"
#define NPRE 1
#define NPOST 1
#define ALPHA 0.33
#define NGMAX 15

void NR_mgfas(u,n,maxcyc)
double **u;
int maxcyc,n;
{
	double NR_anorm2();
	void NR_copy(),NR_interp(),NR_lop(),NR_matadd(),NR_matsub(),
	     NR_relax2(),NR_rstrct(),NR_slvsm2();
	unsigned int j,jcycle,jj,jm1,jpost,jpre,nf,ng=0,ngrid,nn;
	double **irho[NGMAX+1],**irhs[NGMAX+1],**itau[NGMAX+1],
		**itemp[NGMAX+1],**iu[NGMAX+1];
	double res,trerr;

	nn=n;
	while (nn >>= 1) ng++;
	if (n != 1+(1L << ng)) NR_error("n-1 must be a power of 2 in mgfas.");
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
	itau[1]=NR_dmatrix(1,nn,1,nn);
	itemp[1]=NR_dmatrix(1,nn,1,nn);

	/* Initial solution on coarsest grid. */

	NR_slvsm2(iu[1],irho[1]);
	NR_free_dmatrix(irho[1],1,nn,1,nn);

	/* Nested iteration loop. */

	ngrid=ng;
	for (j=2;j<=ngrid;j++) {
		nn=2*nn-1;
		iu[j]=NR_dmatrix(1,nn,1,nn);
		irhs[j]=NR_dmatrix(1,nn,1,nn);
		itau[j]=NR_dmatrix(1,nn,1,nn);
		itemp[j]=NR_dmatrix(1,nn,1,nn);

		/* Interpolate from coarse grid to next
		   finer grid. */

		NR_interp(iu[j],iu[j-1],nn);

		/* Set up rhs. */

		NR_copy(irhs[j],(j != ngrid ? irho[j] : u),nn);
		
		/* V-cycle loop. */

		for (jcycle=1;jcycle<=maxcyc;jcycle++) {
		nf=nn;
			for (jj=j;jj>=2;jj--) {
			  
			  /* Pre-smoothing on downward stroke of
			     V-cycle */
			  
			  for (jpre=1;jpre<=NPRE;jpre++)
		               NR_relax2(iu[jj],irhs[jj],nf);

			  /* Compute the residual Lh(uh) . */
			  
		          NR_lop(itemp[jj],iu[jj],nf);

			  /* Restrict the residual R(Lh(uh)) and
			     solution R(uh). */

			  nf=nf/2+1;
			  jm1=jj-1;
			  NR_rstrct(itemp[jm1],itemp[jj],nf);
			  NR_rstrct(iu[jm1],iu[jj],nf);

			  /* Determine relative truncation error for
			     the coarse grid relative to the fine grid
			     tH = LH(R(uh)) - R(Lh(uh)). */
			  
			  NR_lop(itau[jm1],iu[jm1],nf);
			  NR_matsub(itau[jm1],itemp[jm1],itau[jm1],nf);
			  
			  /* Estimate truncation error of
			     multigrid solution. */

			  if (jj == j)
		              trerr=ALPHA*NR_anorm2(itau[jm1],nf);

			  /* Restrict the rhs fh and determine coarse
			     grid forcing function 
			     rhs = fH = R(fh) + tH. */
			  
			  NR_rstrct(irhs[jm1],irhs[jj],nf);
			  NR_matadd(irhs[jm1],itau[jm1],irhs[jm1],nf);
			}

			/* Bottom of V-cycle; solve on coarsest grid. */

			NR_slvsm2(iu[1],irhs[1]);
			nf=3;
			for (jj=2;jj<=j;jj++) {
			  jm1=jj-1;

			  /* Determine R(uh). */
			  
		          NR_rstrct(itemp[jm1],iu[jj],nf);

			  /* Evaluate uH - R(uh). */
			  
			  NR_matsub(iu[jm1],itemp[jm1],itemp[jm1],nf);

			  /* Prolong the solution corrections to the
			     finer grid, th = P(uH - R(uh)). */
			  
			  nf=2*nf-1;
			  NR_interp(itau[jj],itemp[jm1],nf);
			  
			  /* Update the solution on the fine grid
			     uh = uh + th. */
			  
			  NR_matadd(iu[jj],itau[jj],iu[jj],nf);
			
		          /* Post-smoothing on upward stroke of
		             V-cycle. */

			  for (jpost=1;jpost<=NPOST;jpost++)
			       NR_relax2(iu[jj],irhs[jj],nf);
			}

			/* Form solution residual on finest grid
			   and check for convergence. */
			
			NR_lop(itemp[j],iu[j],nf);
			NR_matsub(itemp[j],irhs[j],itemp[j],nf);
			res=NR_anorm2(itemp[j],nf);
			if (res < trerr) break;
		}
	}

	/* Return solution in u. */

	NR_copy(u,iu[ngrid],n);

	/* De-allocate storage. */

	for (nn=n,j=ng;j>=1;j--,nn=nn/2+1) {
		NR_free_dmatrix(itemp[j],1,nn,1,nn);
		NR_free_dmatrix(itau[j],1,nn,1,nn);
		NR_free_dmatrix(irhs[j],1,nn,1,nn);
		NR_free_dmatrix(iu[j],1,nn,1,nn);
		if (j != ng && j != 1) NR_free_dmatrix(irho[j],1,nn,1,nn);
	}
}
#undef NGMAX
#undef NPRE
#undef NPOST
#undef ALPHA
