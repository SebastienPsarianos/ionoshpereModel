/* Two-Dimensional Array Copy */

/* This subroutine copies one two-dimensional array
   ain(1:n,1:n) to another two-dimensional array
   aout(1:n,1:n). */

void NR_copy(aout,ain,n)
double **ain,**aout;
int n;
{
	int i,j;
	for (i=1;i<=n;i++)
		for (j=1;j<=n;j++)
			aout[j][i]=ain[j][i];

}
