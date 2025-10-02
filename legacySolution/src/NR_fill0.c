/* Two-Dimensional Array Zero */

/* This subroutine zeroes the entries of the
   two-dimensional array u(1:n,1:n). */

void NR_fill0(u,n)
double **u;
int n;
{
	int i,j;
	for (j=1;j<=n;j++)
		for (i=1;i<=n;i++)
			u[i][j]=0.0;
}
