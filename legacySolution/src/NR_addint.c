/* Coarse Grid Correction (CGC) Operator */

/* This subroutine adds the CGC prolonged from the coarse
   grid to the fine grid and computes (updates) the next
   approximation to the solution on the fine grid.

   Parameters:

   uf     Fine-grid solution array
   uc     Coarse-grid solution correction array
   nf     Number of grid points in x- and y-directions on
          fine grid
   res    Temporary storage containing prolonged corrections to
          solution for fine grid

   Additional routines:

   NR_interp  */

void NR_addint(uf,uc,res,nf)
double **res,**uc,**uf;
int nf;
{
	void NR_interp();
	int i,j;

	/* Prolong solution corrections from coarse to fine
	   grid. */

	NR_interp(res,uc,nf);

	/* Update solution on fine grid. */
	
	for (j=1;j<=nf;j++)
		for (i=1;i<=nf;i++)
			uf[i][j] += res[i][j];
}
