/* Full Linear Multigrid Solution Algorithm for the
   Linear Elliptic PDE Governing Ionospheric Currents */

/* This subroutine applies a full linear multigrid solution
   algorithm to the problem of determining the solution of the
   boundary value problem (BVP) for the linear elliptic equation
   that governs planetary ionospheric currents on a spherical
   surface of radius Rp.  The linear elliptic PDE for the
   electric field potential PHI is defined on the domain 
   0 < theta < PI, 0 < psi < 2 PI, and subject to the Dirichlet
   boundary data

        PHI(0,Psi) = 0 ,

        PHI(PI,Psi) = 0 ,

        PHI(Theta,0) = PHI(Theta,2*PI).

   A red-black Gauss-Seidel scheme is used as the smoothing
   operator, bilinear interpolation is used for the prolongation
   operator, and half-weighting is used for the restriction
   operator.
  
   On input the two-diminsional solution array PHI(1:nTheta,1:nPsi)
   contains the right-hand-side source term, RHO, (this is the
   radial component of the magnetospheric current) while on
   output it returns the solution to the elliptic PDE.  The
   number of grid points nTheta should be related to the number of
   grid levels used in the solution, ng, by nTheta = 2^ng+1.

   Parameters:

   PHI     Solution array
   Sigma0  Field aligned conductance (height-integrated conductivity)
   SigmaH  Hall Conductance (height-integrated conductivity)
   SigmaP  Petersen Conductance (height-integrated conductivity)
   Theta   Discrete values of Theta coordinate
   Psi     Discrete values of Phi coordinate
   Radius  Radius of sphere, Rp
   nTheta  Number of grid points in Theta-direction
           on finest grid
   nPsi    Number of grid points in Psi-direction
           on finest grid
   ncycle  Number of V-cyles to be used at each level
   ng      Number of grid levels used
   
   IONOSPHERE_NPRE    Number of relaxation sweeps before coarse-grid
                      correction is computed
   IONOSPHERE_NPOST   Number of relaxation sweeps

   Additional routines:

   IONOSPHERE_cgrid, IONOSPHERE_restrict, IONOSPHERE_prolong,
   IONOSPHERE_smoother, IONOSPHERE_residual, IONOSPHERE_update,
   IONOSPHERE_matcopy, IONOSPHERE_matzero */

#include "IONOSPHERE.h"
#include "IONOSPHERE_util.h"

void IONOSPHERE_multigrid(PHI, Sigma0, SigmaH, SigmaP,
                          SigmaThTh, SigmaThPs, SigmaPsPs,
                          dSigmaThTh_dTheta, dSigmaThPs_dTheta,
                          dSigmaPsPs_dTheta,
                          dSigmaThTh_dPsi, dSigmaThPs_dPsi,
                          dSigmaPsPs_dPsi,
                          Theta, Psi, Radius,
                          nTheta, nPsi, ncycle)
double **PHI, **Sigma0, **SigmaH, **SigmaP,
       **SigmaThTh, **SigmaThPs, **SigmaPsPs,
       **dSigmaThTh_dTheta, **dSigmaThPs_dTheta, **dSigmaPsPs_dTheta,
       **dSigmaThTh_dPsi, **dSigmaThPs_dPsi, **dSigmaPsPs_dPsi,
       **Theta, **Psi, Radius;
int nTheta, nPsi, ncycle;
{
        unsigned int i, j, jcycle, jj, jpost, jpre,
                     ng=0, ngrid,
                     nThetaC, nPsiC, nThetaF, nPsiF,
	             irepeat;
        double **ires[IONOSPHERE_NGMAX+1], **irho[IONOSPHERE_NGMAX+1],
               **irhs[IONOSPHERE_NGMAX+1], **iu[IONOSPHERE_NGMAX+1],
               **ix[IONOSPHERE_NGMAX+1], **iy[IONOSPHERE_NGMAX+1],
               **iS0[IONOSPHERE_NGMAX+1], **iSH[IONOSPHERE_NGMAX+1],
               **iSP[IONOSPHERE_NGMAX+1], **iSxx[IONOSPHERE_NGMAX+1],
               **iSxy[IONOSPHERE_NGMAX+1], **iSyy[IONOSPHERE_NGMAX+1],
               **idSxxdx[IONOSPHERE_NGMAX+1], **idSxxdy[IONOSPHERE_NGMAX+1],
               **idSxydx[IONOSPHERE_NGMAX+1], **idSxydy[IONOSPHERE_NGMAX+1],
               **idSyydx[IONOSPHERE_NGMAX+1], **idSyydy[IONOSPHERE_NGMAX+1];
	double RESIDUAL;

        nThetaC=nTheta;
        while (nThetaC >>= 1) ng++;
        if (nTheta != 1+(1L << ng))
          IONOSPHERE_error("nTheta-1 must be a power of 2 in IONOSPHERE_multigrid.");
        if (ng > IONOSPHERE_NGMAX)
          IONOSPHERE_error("increase IONOSPHERE_NGMAX in IONOSPHERE.h.");

        /* Multiply the right-hand-side source terms
           (in effect the field-aligned current) by
           Radius^2 sin^2(Theta). */

        for (j=1; j<=nPsi; j++) {
          for (i=2; i<nTheta; i++) {
            PHI[i][j]=PHI[i][j]*Radius*Radius*
                      sin(Theta[i][j])*sin(Theta[i][j]);
          } /* endfor */
          PHI[1][j]=PHI[nTheta][j]=0.00;
        } /* endfor */

        /* Allocate storage for the coarse grid and
           rhs on grids ng and ng-1.  Create coarse grid and fill
           rhs by restricting rhs from the fine grid. 
           Also determine the field aligned and Hall and
           Pedersen conductances.  */

        ix[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        iy[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        IONOSPHERE_matcopy(ix[ng], Theta, nTheta, nPsi);
        IONOSPHERE_matcopy(iy[ng], Psi, nTheta, nPsi);

        iS0[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        iSH[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        iSP[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        iSxx[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        iSxy[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        iSyy[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        idSxxdx[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        idSxydx[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        idSyydx[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        idSxxdy[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        idSxydy[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        idSyydy[ng]=IONOSPHERE_dmatrix(1,nTheta,1,nPsi);
        IONOSPHERE_conductance(iS0[ng], iSH[ng], iSP[ng],
                               iSxx[ng], iSxy[ng], iSyy[ng],
                               idSxxdx[ng], idSxydx[ng], idSyydx[ng],
                               idSxxdy[ng], idSxydy[ng], idSyydy[ng],
                               PHI, ix[ng], iy[ng],
                               Radius, nTheta, nPsi);
        
        nThetaC=nTheta/2+1;
        nPsiC=nPsi/2+1;
        ngrid=ng-1;
        irho[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        ix[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        iy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        IONOSPHERE_cgrid(Theta, Psi, ix[ngrid], iy[ngrid], nThetaC, nPsiC);
        IONOSPHERE_restrict(PHI, irho[ngrid], nThetaC, nPsiC);

        iS0[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        iSH[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        iSP[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        iSxx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        iSxy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        iSyy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        idSxxdx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        idSxydx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        idSyydx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        idSxxdy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        idSxydy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        idSyydy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
        IONOSPHERE_conductance(iS0[ngrid], iSH[ngrid], iSP[ngrid],
                               iSxx[ngrid], iSxy[ngrid], iSyy[ngrid],
                               idSxxdx[ngrid], idSxydx[ngrid], idSyydx[ngrid],
                               idSxxdy[ngrid], idSxydy[ngrid], idSyydy[ngrid],
                               irho[ngrid], ix[ngrid], iy[ngrid],
                               Radius, nThetaC, nPsiC);

        /* Similarily allocate storage for other coarse
          grids and fill rhs on all coarse grids. */
       
        while (nThetaC > 3) {
          nThetaC=nThetaC/2+1;
          nPsiC=nPsiC/2+1;
          irho[--ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          ix[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          iy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          IONOSPHERE_cgrid(ix[ngrid+1], iy[ngrid+1],
                           ix[ngrid], iy[ngrid], nThetaC, nPsiC);
          IONOSPHERE_restrict(irho[ngrid+1], irho[ngrid], nThetaC, nPsiC);

          iS0[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          iSH[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          iSP[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          iSxx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          iSxy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          iSyy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          idSxxdx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          idSxydx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          idSyydx[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          idSxxdy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          idSxydy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          idSyydy[ngrid]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
          IONOSPHERE_conductance(iS0[ngrid], iSH[ngrid], iSP[ngrid],
                                 iSxx[ngrid], iSxy[ngrid], iSyy[ngrid],
                                 idSxxdx[ngrid], idSxydx[ngrid], idSyydx[ngrid],
                                 idSxxdy[ngrid], idSxydy[ngrid], idSyydy[ngrid],
                                 irho[ngrid], ix[ngrid], iy[ngrid],
                                 Radius, nThetaC, nPsiC);
       }
       
       /* Initial solution on coarsest grid.  Solve this
          problem using IONOSPHERE_NEXACT Gauss-Seidel interations. */

       iu[1]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
       irhs[1]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
       IONOSPHERE_matzero(iu[1], nThetaC, nPsiC);
       for (j=1; j <= IONOSPHERE_NEXACT; j++)
         IONOSPHERE_smoother(iu[1], irho[1],
                             iSxx[1], iSxy[1], iSyy[1],
                             idSxxdx[1], idSxydx[1], idSyydx[1],
                             idSxxdy[1], idSxydy[1], idSyydy[1],
                             ix[1], iy[1], Radius, nThetaC, nPsiC);
	
       IONOSPHERE_free_dmatrix(irho[1],1,nThetaC,1,nPsiC);

       /* Nested iteration loop. */
       
       ngrid=ng;
       for (j=2; j<=ngrid; j++) {

         nThetaC=2*nThetaC-1;
         nPsiC=2*nPsiC-1;
         iu[j]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
         irhs[j]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);
         ires[j]=IONOSPHERE_dmatrix(1,nThetaC,1,nPsiC);

         /* Prolong from coarse grid to next
            finer grid using bilinear interpolation. */
              
         IONOSPHERE_prolong(iu[j], iu[j-1], nThetaC, nPsiC);
	     
         /* Set up rhs. */
              
         IONOSPHERE_matcopy(irhs[j],(j != ngrid ? irho[j] : PHI),
                            nThetaC, nPsiC);
         
         /* V-cycle loop. */

	 irepeat = 0;
         for (jcycle=1;jcycle<=ncycle;jcycle++) {

           begin_V_cycle: ;

           nThetaF=nThetaC;
           nPsiF=nPsiC;
           for (jj=j; jj>=2; jj--) {

             /* Pre-smoothing on downward stroke of
                V-cycle */

             for (jpre=1; jpre<=IONOSPHERE_NPRE; jpre++)
               IONOSPHERE_smoother(iu[jj], irhs[jj],
                                   iSxx[jj], iSxy[jj], iSyy[jj],
                                   idSxxdx[jj], idSxydx[jj], idSyydx[jj],
                                   idSxxdy[jj], idSxydy[jj], idSyydy[jj],
                                   ix[jj], iy[jj],
                                   Radius, nThetaF, nPsiF);
             
             /* Restricted residual is the next rhs. */
                    
             IONOSPHERE_residual(ires[jj], iu[jj], irhs[jj],
                                 iSxx[jj], iSxy[jj], iSyy[jj],
                                 idSxxdx[jj], idSxydx[jj], idSyydx[jj],
                                 idSxxdy[jj], idSxydy[jj], idSyydy[jj],
                                 ix[jj], iy[jj],
                                 Radius, nThetaF, nPsiF, &RESIDUAL, 0);

	     nThetaF=nThetaF/2+1;
             nPsiF=nPsiF/2+1;

             IONOSPHERE_restrict(ires[jj], irhs[jj-1], 
                                 nThetaF, nPsiF);

             /* Zero for initial guess in next relaxation. */
                     
             IONOSPHERE_matzero(iu[jj-1], nThetaF, nPsiF);
           }

           /* Bottom of V-cycle; solve on coarsest grid. */
                     
           for (jj=1; jj <= IONOSPHERE_NEXACT; jj++)
             IONOSPHERE_smoother(iu[1], irhs[1],
                                 iSxx[1], iSxy[1], iSyy[1],
                                 idSxxdx[1], idSxydx[1], idSyydx[1],
                                 idSxxdy[1], idSxydy[1], idSyydy[1],
                                 ix[1], iy[1],
                                 Radius, nThetaF, nPsiF);

           /* Go back up. */
	   
           for (jj=2; jj<=j; jj++) {
             nThetaF=2*nThetaF-1;
             nPsiF=2*nPsiF-1;

             /* Prolong solution corrections from coarse to finer
                grid. */
             
             IONOSPHERE_prolong(ires[jj], iu[jj-1], nThetaF, nPsiF);
	     
             /* Update the solution on the finer grid. */
             
             IONOSPHERE_update(iu[jj], ires[jj], IONOSPHERE_OMEGA, nThetaF, nPsiF);
             
             /* Post-smoothing on upward stroke of
                V-cycle. */
                       
             for (jpost=1; jpost<=IONOSPHERE_NPOST; jpost++)
               IONOSPHERE_smoother(iu[jj], irhs[jj],
                                   iSxx[jj], iSxy[jj], iSyy[jj],
                                   idSxxdx[jj], idSxydx[jj], idSyydx[jj],
                                   idSxxdy[jj], idSxydy[jj], idSyydy[jj],
                                   ix[jj], iy[jj],
                                   Radius, nThetaF, nPsiF);

             /* Determine the residuals on the finer grid. */

             IONOSPHERE_residual(ires[jj], iu[jj], irhs[jj],
                                 iSxx[jj], iSxy[jj], iSyy[jj],
                                 idSxxdx[jj], idSxydx[jj], idSyydx[jj],
                                 idSxxdy[jj], idSxydy[jj], idSyydy[jj],
                                 ix[jj], iy[jj],
                                 Radius, nThetaF, nPsiF, &RESIDUAL, 1);
             printf("\n Grid: %d x %d, Residual = %g",
                    nThetaF-1, nPsiF-1, RESIDUAL);
	     
           } /* endfor. */

	   /* Check for convergence at this grid level
              and continue futher V-cycles as necessary. */

           /* if (j == ngrid && jcycle == ncycle &&
	       RESIDUAL > IONOSPHERE_TOLER && irepeat <= 2) {
	     irepeat = irepeat + 1;
             goto begin_V_cycle;
           } */

         } /* endfor (V-cyle). */

       } /* endfor (FMG). */ 

       /* Return solution in u. */
       
       IONOSPHERE_matcopy(PHI, iu[ngrid], nTheta, nPsi);
       
       /* De-allocate storage. */

       for (nThetaC=nTheta, nPsiC=nPsi, j=ng;
            j>=2;
            j--, nThetaC=nThetaC/2+1, nPsiC=nPsiC/2+1) {
         IONOSPHERE_free_dmatrix(ires[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(irhs[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iu[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(ix[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iy[j],1,nThetaC,1,nPsiC);
         if (j != ng) IONOSPHERE_free_dmatrix(irho[j],1,nThetaC,1,nPsiC);

         IONOSPHERE_free_dmatrix(iS0[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iSH[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iSP[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iSxx[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iSxy[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(iSyy[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(idSxxdx[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(idSxydx[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(idSyydx[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(idSxxdy[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(idSxydy[j],1,nThetaC,1,nPsiC);
         IONOSPHERE_free_dmatrix(idSyydy[j],1,nThetaC,1,nPsiC);
       }
       IONOSPHERE_free_dmatrix(irhs[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iu[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(ix[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iy[1],1,nThetaC,1,nPsiC);

       IONOSPHERE_free_dmatrix(iS0[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iSH[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iSP[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iSxx[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iSxy[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(iSyy[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(idSxxdx[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(idSxydx[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(idSyydx[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(idSxxdy[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(idSxydy[1],1,nThetaC,1,nPsiC);
       IONOSPHERE_free_dmatrix(idSyydy[1],1,nThetaC,1,nPsiC);
}

