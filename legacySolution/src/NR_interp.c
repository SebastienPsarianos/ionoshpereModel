/* Bilinear Interpolation Prolongation Operator */

/* This subroutine uses bilinear interpolation
   to prolong the coarse grid solution to the fine grid.

   Parameters:

   uc     Coarse-grid solution array
   uf     Fine-grid solution array
   nf     Number of grid points in x- and y-directions on
          fine grid */

void NR_interp(uf,uc,nf)
double **uc,**uf;
int nf;
{
	int ic,iif,jc,jf,nc;
	nc=nf/2+1;

	/* Do elements that are direct copies. */
	
	for (jc=1,jf=1;jc<=nc;jc++,jf+=2)
		for (ic=1;ic<=nc;ic++) uf[2*ic-1][jf]=uc[ic][jc];

	/* Do odd-numbered columns, interpolating vertically. */
	
	for (jf=1;jf<=nf;jf+=2)
		for (iif=2;iif<nf;iif+=2)
			uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);

        /* Do even-numbered columns, interpolating horizontally. */

	for (jf=2;jf<nf;jf+=2)
		for (iif=1;iif <= nf;iif++)
			uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);
}
