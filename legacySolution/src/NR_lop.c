/* Nonlinear Elliptic Equation Solution Operator L */

/* This subroutine uses evaluates the nonlinear
   elliptic equation solution operator

          d^2 u      d^2 u
	  -----  +   -----   + u^2  =  L(u) ,
	  dx^2       dx^2

   using the discretized formulation

                      n        n+1      n       n+1        n
        L(u   )  =  (u      + u      + u     + u      - 4 u    ) /
	   i,j        i+1,j    i-1,j    i,j+1   i,j-1      i,j

                       2     n   2    
	             dx  + (u   )
		             i,j

   Given u(1:n,1:n), the routine returns L(u) in
   out(1:n,1:n).
	                                    
   Parameters:

   u      Solution array 
   out    L(u)
   n      Number of grid points in x- and y-directions */

void NR_lop(out,u,n)
double **out,**u;
int n;
{
	int i,j;
	double h,h2i;

	h=1.0/(n-1);
	h2i=1.0/(h*h);

	/* Interior points */

	for (j=2;j<n;j++)
	  for (i=2;i<n;i++)
	    out[i][j]=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
		      4.0*u[i][j])+u[i][j]*u[i][j];

	/* Boundary points */
	
	for (i=1;i<=n;i++)
	  out[i][1]=out[i][n]=out[1][i]=out[n][i]=0.0;
}
