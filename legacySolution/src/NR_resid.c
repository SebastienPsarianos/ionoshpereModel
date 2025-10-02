/* Linear Elliptic Equation Residual Operator Res */

/* This subroutine uses evaluates the linear
   elliptic equation residual operator

           d^2 u      d^2 u
      r -  -----  +   -----  =  Res(u) ,
	   dx^2       dx^2

   using the discretized formulation

                n        n+1      n       n+1        n
      Res = - (u      + u      + u     + u      - 4 u    ) /
                i,j       i+1,j    i-1,j    i,j+1   i,j-1      i,j

              2
            dx  + r
       	      i,j

   Given the solution u(1:n,1:n) and right-hand-side term rhs(1:n,1:n),
   the routine returns Res(u) in res(1:n,1:n).
	                                    
   Parameters:

   res    Residual array
   u      Solution array 
   rhs    rhs
   n      Number of grid points in x- and y-directions */

void NR_resid(res,u,rhs,n)
double **res,**rhs,**u;
int n;
{
	int i,j;
	double h,h2i;

	h=1.0/(n-1);
	h2i=1.0/(h*h);
	

	/* Interior points */

	for (j=2;j<n;j++)
	  for (i=2;i<n;i++)
	    res[i][j]=-h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
		      4.0*u[i][j])+rhs[i][j];
	
	/* Boundary points */

	for (i=1;i<=n;i++)
		res[i][1]=res[i][n]=res[1][i]=res[n][i]=0.0;
}
