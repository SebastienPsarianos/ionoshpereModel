/* Newton Red-Black Gauss-Seidel Smoothing Operator */

/* This subroutine uses a red-black Gauss-Seidel scheme
   with Newton accleration to smooth the solution to the
   linear elliptic equation

          d^2 u      d^2 u
	  -----  +   -----   + u^2  =  r  ,
	  dx^2       dx^2

   defined on the square domain 0 < x < 1, 0 y < 1, and
   subject to the Dirichlet boundary data

          u(0,y) = u(1,y) = u(x,0) = u(x,1) = 0.

   The current value of the solution , u, is updated using
   the right-hand side function, rhs, as follows:

                     n        n+1      n       n+1        n
	  du    =  (u      + u      + u     + u      - 4 u    ) /
	    i,j      i+1,j    i-1,j    i,j+1   i,j-1      i,j

                     2     n   2    
	           dx  + (u   )   -  r
		           i,j        i,j

			   
           n+1    n                  n       4
	  u    = u    - du    / ( 2 u    -  ---)
	   i,j    i,j     i,j        i,j      2
                                            dx
	                                    
   Parameters:

   u      Solution array 
   rhs    rhs array, r(i,j)
   n      Number of grid points in x- and y-directions */

void NR_relax2(u,rhs,n)
double **rhs,**u;
int n;
{
	int i,ipass,isw,j,jsw=1;
	double foh2,h,h2i,res;

	h=1.0/(n-1);
	h2i=1.0/(h*h);
	foh2 = -4.0*h2i;

	/* Red and black sweeps. */
	
	for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
		isw=jsw;
		for (j=2;j<n;j++,isw=3-isw) {
		  for (i=isw+1;i<n;i+=2) {

		    /* Newton Gauss-Seidel formula. */
			  
		    res=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
			4.0*u[i][j])+u[i][j]*u[i][j]-rhs[i][j];
		        u[i][j] -= res/(foh2+2.0*u[i][j]);
		  }
		}
	}
}
