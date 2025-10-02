/* Red-Black Gauss-Seidel Smoothing Operator */

/* This subroutine uses a red-black Gauss-Seidel scheme
   to smooth the solution to the linear elliptic equation

          d^2 u      d^2 u
	  -----  +   -----   =   r   ,
	  dx^2       dx^2

   defined on the square domain 0 < x < 1, 0 y < 1, and
   subject to the Dirichlet boundary data

          u(0,y) = u(1,y) = u(x,0) = u(x,1) = 0.

   The current value of the solution , u, is updated using
   the right-hand side function, rhs, as follows:

           n+1    1   n        n+1      n       n+1
	  u    =  - (u      + u      + u     + u     ) -
	   i,j    4   i+1,j    i-1,j    i,j+1   i,j-1

                    dx^2
	          - ---- r
		      4   i,j

   Parameters:

   u      Solution array 
   rhs    rhs array
   n      Number of grid points in x- and y-directions */

void NR_relax(u,rhs,n)
double **rhs,**u;
int n;
{
	int i,ipass,isw,j,jsw=1;
	double h,h2;

	h=1.0/(n-1);
	h2=h*h;

	/* Red and black sweeps. */
	
	for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
		isw=jsw;
		for (j=2;j<n;j++,isw=3-isw)
			for (i=isw+1;i<n;i+=2)
			  
			  /* Gauss-Seidel formula. */
			  
			  u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]
			          +u[i][j-1]-h2*rhs[i][j]);
	}
}
