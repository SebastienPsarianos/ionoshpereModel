/* Two-Dimensional Array Norm */

/* This subroutine returns the Euclidean norm of
   the two-dimensional array a(1:n,1:n). */

#include <math.h>

double NR_anorm2(a,n)
double **a;
int n;
{
	int i,j;
	double sum=0.0;

	for (j=1;j<=n;j++)
		for (i=1;i<=n;i++)
			sum += a[i][j]*a[i][j];
	return sqrt(sum)/n;
}
