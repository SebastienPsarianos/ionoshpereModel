/* Two-Dimensional Array Add */

/* This subroutine adds one two-dimensional matrix
   a(1:n,1:n) to another two-dimensional matrix
   b(1:n,1:n) and returns result in c(1:n,1:n). */

void NR_matadd(a,b,c,n)
double **a,**b,**c;
int n;
{
	int i,j;

	for (j=1;j<=n;j++)
		for (i=1;i<=n;i++)
			c[i][j]=a[i][j]+b[i][j];
}
