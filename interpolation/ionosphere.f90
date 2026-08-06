!====================================
!                                   |
!         Ionosphere Model          |
!                                   |
!====================================

subroutine ionosphere(iter,iAction)
  !\
  ! Ionosphere driver (main or controlling) routine.
  !/
  use ModIonosphere
  implicit none

  integer, intent(in) :: iter,iAction

  real :: Radius

  select case (iAction)
     !\
     ! Create fine grids for north and south 
     ! hemisphere ionospheric solutions and evaluate the
     ! appropriate conductances.
     !/
     case (1)
        IONO_PI = 2.00*asin(1.00)
        call ionosphere_fine_grid
        call ionosphere_init
        call ionosphere_conductance(IONO_NORTH_Sigma0, IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
                                    IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                    IONO_NORTH_SigmaPsPs, &
                                    IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
                                    IONO_NORTH_dSigmaPsPs_dTheta, &
                                    IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
                                    IONO_NORTH_dSigmaPsPs_dPsi, &
                                    IONO_NORTH_Theta, IONO_NORTH_Psi, IONO_nTheta, IONO_nPsi)
        call ionosphere_conductance(IONO_SOUTH_Sigma0, IONO_SOUTH_SigmaH, IONO_SOUTH_SigmaP, &
                                    IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                    IONO_SOUTH_SigmaPsPs, &
                                    IONO_SOUTH_dSigmaThTh_dTheta, IONO_SOUTH_dSigmaThPs_dTheta, &
                                    IONO_SOUTH_dSigmaPsPs_dTheta, &
                                    IONO_SOUTH_dSigmaThTh_dPsi, IONO_SOUTH_dSigmaThPs_dPsi, &
                                    IONO_SOUTH_dSigmaPsPs_dPsi, &
                                    IONO_SOUTH_Theta, IONO_SOUTH_Psi, IONO_nTheta, IONO_nPsi)
        if (IONO_NORTH_nMagBndPts >= 32 .and. &
            IONO_SOUTH_nMagBndPts >= 32) call ionosphere_fac

     !\
     ! Perform ionophere calculations for north and south 
     ! hemispheres using magnetospheric FAC solution.
     ! A multigrid technique is used to solve the elliptic
     ! equations for the ionospheric potential on each of
     ! the hemispheres.
     !/
     case (2)
        if (IONO_NORTH_nMagBndPts >= 32 .and. &
            IONO_SOUTH_nMagBndPts >= 32) call ionosphere_fac
        Radius = IONO_Radius + IONO_Height

        IONO_NORTH_PHI = IONO_NORTH_JR
        call ionosphere_multigrid(IONO_NORTH_PHI, &
                                  IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                  IONO_NORTH_SigmaPsPs, &
                                  IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
                                  IONO_NORTH_dSigmaPsPs_dTheta, &
                                  IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
                                  IONO_NORTH_dSigmaPsPs_dPsi, &
                                  IONO_NORTH_Theta, IONO_NORTH_Psi, Radius, IONO_nTheta, IONO_nPsi, &
                                  4, 0)
        call ionosphere_currents(IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz, &
                                 IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
                                 IONO_NORTH_ETh,IONO_NORTH_EPs, &
                                 IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
                                 IONO_NORTH_PHI, &
                                 IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                 IONO_NORTH_SigmaPsPs, &
                                 IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                                 IONO_NORTH_Theta, IONO_NORTH_Psi, &
                                 Radius, IONO_nTheta, IONO_nPsi)
        call ionosphere_magBCs(IONO_NORTH_PHI_BC, IONO_NORTH_ETh_BC, IONO_NORTH_EPs_BC, &
                               IONO_NORTH_UR_BC, IONO_NORTH_UTh_BC, IONO_NORTH_UPs_BC, &
                               IONO_Radius_Mag_Boundary, IONO_NORTH_PHI, &
                               IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                               IONO_NORTH_Theta, IONO_NORTH_Psi, &
                               Radius, IONO_nTheta, IONO_nPsi)

        IONO_SOUTH_PHI = IONO_SOUTH_JR
        call ionosphere_multigrid(IONO_SOUTH_PHI, &
                                  IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                  IONO_SOUTH_SigmaPsPs, &
                                  IONO_SOUTH_dSigmaThTh_dTheta, IONO_SOUTH_dSigmaThPs_dTheta, &
                                  IONO_SOUTH_dSigmaPsPs_dTheta, &
                                  IONO_SOUTH_dSigmaThTh_dPsi, IONO_SOUTH_dSigmaThPs_dPsi, &
                                  IONO_SOUTH_dSigmaPsPs_dPsi, &
                                  IONO_SOUTH_Theta, IONO_SOUTH_Psi, Radius, IONO_nTheta, IONO_nPsi, &
                                  4, 1)
        call ionosphere_currents(IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz, &
                                 IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
                                 IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
                                 IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
                                 IONO_SOUTH_PHI, &
                                 IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                 IONO_SOUTH_SigmaPsPs, &
                                 IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
                                 IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                                 Radius, IONO_nTheta, IONO_nPsi)
        call ionosphere_magBCs(IONO_SOUTH_PHI_BC, IONO_SOUTH_ETh_BC, IONO_SOUTH_EPs_BC, &
                               IONO_SOUTH_UR_BC, IONO_SOUTH_UTh_BC, IONO_SOUTH_UPs_BC, &
                               IONO_Radius_Mag_Boundary, IONO_SOUTH_PHI, &
                               IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
                               IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                               Radius, IONO_nTheta, IONO_nPsi)

     !\
     ! Write ionospheric solution for north and south 
     ! hemispheres to output data file.
     !/
     case (3)
        call ionosphere_write_output(iter)

     !\
     ! Create a restart solution file containing the
     ! FAC driving the ionospheric solution.
     !/
     case (4)
        call ionosphere_write_restart_file

     !\
     ! Create fine grids for north and south 
     ! hemisphere ionospheric solutions and evaluate the
     ! appropriate conductances.  Then perform ionophere
     ! calculations for north and south hemispheres using 
     ! FAC solution from restart file.
     !/
     case (5)
        IONO_PI = 2.00*asin(1.00)
        call ionosphere_fine_grid
        call ionosphere_init
        call ionosphere_conductance(IONO_NORTH_Sigma0, IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
                                    IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                    IONO_NORTH_SigmaPsPs, &
                                    IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
                                    IONO_NORTH_dSigmaPsPs_dTheta, &
                                    IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
                                    IONO_NORTH_dSigmaPsPs_dPsi, &
                                    IONO_NORTH_Theta, IONO_NORTH_Psi, IONO_nTheta, IONO_nPsi)
        call ionosphere_conductance(IONO_SOUTH_Sigma0, IONO_SOUTH_SigmaH, IONO_SOUTH_SigmaP, &
                                    IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                    IONO_SOUTH_SigmaPsPs, &
                                    IONO_SOUTH_dSigmaThTh_dTheta, IONO_SOUTH_dSigmaThPs_dTheta, &
                                    IONO_SOUTH_dSigmaPsPs_dTheta, &
                                    IONO_SOUTH_dSigmaThTh_dPsi, IONO_SOUTH_dSigmaThPs_dPsi, &
                                    IONO_SOUTH_dSigmaPsPs_dPsi, &
                                    IONO_SOUTH_Theta, IONO_SOUTH_Psi, IONO_nTheta, IONO_nPsi)

        call ionosphere_read_restart_file
        Radius = IONO_Radius + IONO_Height

        IONO_NORTH_PHI = IONO_NORTH_JR
        call ionosphere_multigrid(IONO_NORTH_PHI, &
                                  IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                  IONO_NORTH_SigmaPsPs, &
                                  IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
                                  IONO_NORTH_dSigmaPsPs_dTheta, &
                                  IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
                                  IONO_NORTH_dSigmaPsPs_dPsi, &
                                  IONO_NORTH_Theta, IONO_NORTH_Psi, Radius, IONO_nTheta, IONO_nPsi, &
                                  4, 0)
!!$        call ionosphere_hall_conductance(IONO_NORTH_Sigma0, IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
!!$                                         IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
!!$                                         IONO_NORTH_SigmaPsPs, &
!!$                                         IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
!!$                                         IONO_NORTH_dSigmaPsPs_dTheta, &
!!$                                         IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
!!$                                         IONO_NORTH_dSigmaPsPs_dPsi, &
!!$                                         IONO_NORTH_Theta, IONO_NORTH_Psi, IONO_nTheta, IONO_nPsi)
!!$        call ionosphere_ped_conductance(IONO_NORTH_Sigma0, IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
!!$                                        IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
!!$                                        IONO_NORTH_SigmaPsPs, &
!!$                                        IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
!!$                                        IONO_NORTH_dSigmaPsPs_dTheta, &
!!$                                        IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
!!$                                        IONO_NORTH_dSigmaPsPs_dPsi, &
!!$                                        IONO_NORTH_Theta, IONO_NORTH_Psi, IONO_nTheta, IONO_nPsi)
        call ionosphere_currents(IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz, &
                                 IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
                                 IONO_NORTH_ETh,IONO_NORTH_EPs, &
                                 IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
                                 IONO_NORTH_PHI, &
                                 IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                 IONO_NORTH_SigmaPsPs, &
                                 IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                                 IONO_NORTH_Theta, IONO_NORTH_Psi, &
                                 Radius, IONO_nTheta, IONO_nPsi)

        IONO_SOUTH_PHI = IONO_SOUTH_JR
        call ionosphere_multigrid(IONO_SOUTH_PHI, &
                                  IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                  IONO_SOUTH_SigmaPsPs, &
                                  IONO_SOUTH_dSigmaThTh_dTheta, IONO_SOUTH_dSigmaThPs_dTheta, &
                                  IONO_SOUTH_dSigmaPsPs_dTheta, &
                                  IONO_SOUTH_dSigmaThTh_dPsi, IONO_SOUTH_dSigmaThPs_dPsi, &
                                  IONO_SOUTH_dSigmaPsPs_dPsi, &
                                  IONO_SOUTH_Theta, IONO_SOUTH_Psi, Radius, IONO_nTheta, IONO_nPsi, &
                                  4, 1)
        call ionosphere_currents(IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz, &
                                 IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
                                 IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
                                 IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
                                 IONO_SOUTH_PHI, &
                                 IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                 IONO_SOUTH_SigmaPsPs, &
                                 IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
                                 IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                                 Radius, IONO_nTheta, IONO_nPsi)

  end select  

end subroutine ionosphere

subroutine ionosphere_fine_grid
  !\
  ! This routine sets the fine grid meshes for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  implicit none

  integer :: i,j
  real :: dTheta, dPsi
        
  dTheta = 0.50*IONO_PI/real(IONO_nTheta-1)
  dPsi = dTheta

  do j = 1, IONO_nPsi
     IONO_NORTH_Theta(1,j) = IONO_Theta_0
     IONO_NORTH_Psi(1,j) = real(j-1)*dPsi
     IONO_NORTH_X(1,j) = sin(IONO_NORTH_Theta(1,j))* &
                         cos(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Y(1,j) = sin(IONO_NORTH_Theta(1,j))* &
                         sin(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Z(1,j) = cos(IONO_NORTH_Theta(1,j))

     do i = 1, IONO_nTheta
        IONO_NORTH_Theta(i,j) = real(i-1)*dTheta
        IONO_NORTH_Psi(i,j) = real(j-1)*dPsi
        IONO_NORTH_X(i,j) = sin(IONO_NORTH_Theta(i,j))* &
                            cos(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Y(i,j) = sin(IONO_NORTH_Theta(i,j))* &
                            sin(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Z(i,j) = cos(IONO_NORTH_Theta(i,j))
     end do  

     IONO_NORTH_Theta(IONO_nTheta,j) = 0.50*IONO_PI
     IONO_NORTH_Psi(IONO_nTheta,j) = real(j-1)*dPsi
     IONO_NORTH_X(IONO_nTheta,j) = sin(IONO_NORTH_Theta(IONO_nTheta,j))* &
                                   cos(IONO_NORTH_Psi(IONO_nTheta,j))
     IONO_NORTH_Y(IONO_nTheta,j) = sin(IONO_NORTH_Theta(IONO_nTheta,j))* &
                                   sin(IONO_NORTH_Psi(IONO_nTheta,j))
     IONO_NORTH_Z(IONO_nTheta,j) = cos(IONO_NORTH_Theta(IONO_nTheta,j))
  end do

  do j = 1, IONO_nPsi
     IONO_SOUTH_Theta(1,j) = IONO_PI-0.50*IONO_PI
     IONO_SOUTH_Psi(1,j) = real(j-1)*dPsi
     IONO_SOUTH_X(1,j) = sin(IONO_SOUTH_Theta(1,j))* &
                         cos(IONO_SOUTH_Psi(1,j))
     IONO_SOUTH_Y(1,j) = sin(IONO_SOUTH_Theta(1,j))* &
                         sin(IONO_SOUTH_Psi(1,j))
     IONO_SOUTH_Z(1,j) = cos(IONO_SOUTH_Theta(1,j))

     do i = 1, IONO_nTheta
        IONO_SOUTH_Theta(i,j) = IONO_PI-0.50*IONO_PI+real(i-1)*dTheta
        IONO_SOUTH_Psi(i,j) = real(j-1)*dPsi
        IONO_SOUTH_X(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
                            cos(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Y(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
                            sin(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Z(i,j) = cos(IONO_SOUTH_Theta(i,j))
     end do  

     IONO_SOUTH_Theta(IONO_nTheta,j) = IONO_PI-IONO_Theta_0
     IONO_SOUTH_Psi(IONO_nTheta,j) = real(j-1)*dPsi
     IONO_SOUTH_X(IONO_nTheta,j) = sin(IONO_SOUTH_Theta(IONO_nTheta,j))* &
                                   cos(IONO_SOUTH_Psi(IONO_nTheta,j))
     IONO_SOUTH_Y(IONO_nTheta,j) = sin(IONO_SOUTH_Theta(IONO_nTheta,j))* &
                                   sin(IONO_SOUTH_Psi(IONO_nTheta,j))
     IONO_SOUTH_Z(IONO_nTheta,j) = cos(IONO_SOUTH_Theta(IONO_nTheta,j))
  end do

end subroutine ionosphere_fine_grid

subroutine ionosphere_init
  !\
  ! This routine initializes the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  implicit none

  IONO_NORTH_PHI = 0.00
  IONO_NORTH_JR = 0.00 
  IONO_NORTH_Jx = 0.00
  IONO_NORTH_Jy = 0.00
  IONO_NORTH_Jz = 0.00
  IONO_NORTH_Ex = 0.00
  IONO_NORTH_Ey = 0.00
  IONO_NORTH_Ez = 0.00
  IONO_NORTH_ETh = 0.00
  IONO_NORTH_EPs = 0.00
  IONO_NORTH_Ux = 0.00
  IONO_NORTH_Uy = 0.00
  IONO_NORTH_Uz = 0.00

  IONO_SOUTH_PHI = 0.00
  IONO_SOUTH_JR = 0.00
  IONO_SOUTH_Jx = 0.00
  IONO_SOUTH_Jy = 0.00
  IONO_SOUTH_Jz = 0.00
  IONO_SOUTH_Ex = 0.00
  IONO_SOUTH_Ey = 0.00
  IONO_SOUTH_Ez = 0.00
  IONO_SOUTH_ETh = 0.00
  IONO_SOUTH_EPs = 0.00
  IONO_SOUTH_Ux = 0.00
  IONO_SOUTH_Uy = 0.00
  IONO_SOUTH_Uz = 0.00  

  IONO_NORTH_PHI_BC = 0.00
  IONO_NORTH_ETh_BC = 0.00
  IONO_NORTH_EPs_BC = 0.00
  IONO_NORTH_UR_BC = 0.00 
  IONO_NORTH_UTh_BC = 0.00 
  IONO_NORTH_UPs_BC = 0.00

  IONO_SOUTH_PHI_BC = 0.00
  IONO_SOUTH_ETh_BC = 0.00
  IONO_SOUTH_EPs_BC = 0.00
  IONO_SOUTH_UR_BC = 0.00
  IONO_SOUTH_UTh_BC = 0.00 
  IONO_SOUTH_UPs_BC = 0.00

end subroutine ionosphere_init

subroutine ionosphere_write_output(iter)
  !\
  ! This routine writes out the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres to
  ! a data file.
  !/
  use ModIonosphere
  implicit none

  integer, intent(in) :: iter

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=7) :: NUMiter

  write(*,*) '=> Writing output datafiles for ionosphere.'

  if (iter < 10) then
     write (NUMiter,'(a6,i1)') 'n00000',iter
  elseif (iter < 100) then
     write (NUMiter,'(a5,i2)') 'n0000',iter
  elseif (iter < 1000) then
     write (NUMiter,'(a4,i3)') 'n000',iter
  elseif (iter < 10000) then
     write (NUMiter,'(a3,i4)') 'n00',iter
  elseif (iter < 100000) then
     write (NUMiter,'(a2,i5)') 'n0',iter
  elseif (iter < 1000000) then
     write (NUMiter,'(a1,i6)') 'n',iter
  else
     write(*,*) "ionosphere_output: iter = ",iter, &
                " Error, too many iterations for filename"
     stop
  end if

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere"//"_"//NUMiter//IO_ext, &
       status="unknown")

  write(Iunit, *)  'TITLE="BATSRUS: Ionospheric Potential Solution"'
  write(Iunit, *)  'VARIABLES= "X [R]","Y [R]","Z [R]"'
  write(Iunit, *)  ' "Theta [deg]","Psi [deg]"'
  write(Iunit, *)  ' "SigmaH [S]","SigmaP [S]"'
  write(Iunit, *)  ' "JR [`mA/m^2]","PHI [kV]"'
  write(Iunit, *)  ' "Ex [mV/m]","Ey [mV/m]","Ez [mV/m]"'
  write(Iunit, *)  ' "Jx [`mA/m^2]","Jy [`mA/m^2]","Jz [`mA/m^2]"'
  write(Iunit, *)  ' "Ux [km/s]","Uy [km/s]","Uz [km/s]"'
  write(Iunit, *)  'ZONE T="Northern Hemisphere" ','I= ',IONO_nTheta, &
                   ' J= ',IONO_nPsi,' F=POINT'
  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,fmt="(18(E13.5))")  &
              IONO_NORTH_X(i,j),IONO_NORTH_Y(i,j),IONO_NORTH_Z(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_NORTH_Psi(i,j), &
              IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
              1.0e06*IONO_NORTH_JR(i,j),1.0e-03*IONO_NORTH_PHI(i,j), &
              1.0e03*IONO_NORTH_Ex(i,j),1.0e03*IONO_NORTH_Ey(i,j), &
              1.0e03*IONO_NORTH_Ez(i,j), &
              1.0e06*IONO_NORTH_Jx(i,j),1.0e06*IONO_NORTH_Jy(i,j), &
              1.0e06*IONO_NORTH_Jz(i,j), &
              1.0e-03*IONO_NORTH_Ux(i,j),1.0e-03*IONO_NORTH_Uy(i,j), &
              1.0e-03*IONO_NORTH_Uz(i,j)
     end do
  end do
  write(Iunit, *)  'ZONE T="Southern Hemisphere" ','I= ',IONO_nTheta, &
                   ' J= ',IONO_nPsi,' F=POINT'
  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,fmt="(18(E13.5))")  &
              IONO_SOUTH_X(i,j),IONO_SOUTH_Y(i,j),IONO_SOUTH_Z(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j), &
              IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
              1.0e06*IONO_SOUTH_JR(i,j),1.0e-03*IONO_SOUTH_PHI(i,j), &
              1.0e03*IONO_SOUTH_Ex(i,j),1.0e03*IONO_SOUTH_Ey(i,j), &
              1.0e03*IONO_SOUTH_Ez(i,j), &
              1.0e06*IONO_SOUTH_Jx(i,j),1.0e06*IONO_SOUTH_Jy(i,j), &
              1.0e06*IONO_SOUTH_Jz(i,j), &
              1.0e-03*IONO_SOUTH_Ux(i,j),1.0e-03*IONO_SOUTH_Uy(i,j), &
              1.0e-03*IONO_SOUTH_Uz(i,j)
     end do
  end do

  close(UNIT=Iunit)

end subroutine ionosphere_write_output

subroutine ionosphere_read_restart_file
  !\
  ! This routine creates an ionospheric 
  ! restart solution file.
  !/
  use ModIonosphere
  implicit none

  integer :: Iunit, i, j, i2, j2, nhemi
  real :: xx,yy
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".rst"

  write(*,*) '=> Reading restart file for ionosphere.'

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere"//IO_ext, &
       status="old")

  read(Iunit,*) IONO_Radius, IONO_Height, &
                IONO_Ref_Density, IONO_Ref_SoundSpeed, &
                IONO_Bdp, IONO_Radius_Mag_Boundary, &
                IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        read(Iunit,*) i2,j2,xx,yy,IONO_NORTH_JR(i,j)
     end do
  end do

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        read(Iunit,*) i2,j2,xx,yy,IONO_SOUTH_JR(i,j) 
     end do
  end do

!!$  nhemi = (IONO_nPsi-1)/2 + 1
!!$  do j = nhemi+1, IONO_nPsi
!!$     do i = 1, IONO_nTheta
!!$        IONO_NORTH_JR(i,j) = - IONO_NORTH_JR(i,nhemi-(j-nhemi))
!!$     end do
!!$  end do
!!$

  close(UNIT=Iunit)

end subroutine ionosphere_read_restart_file

subroutine ionosphere_write_restart_file
  !\
  ! This routine reads in an ionospheric 
  ! restart solution file.
  !/
  use ModIonosphere
  implicit none

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".rst"

  write(*,*) '=> Writing restart file for ionosphere.'

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere"//IO_ext, &
       status="unknown")

  write(Iunit,*) IONO_Radius, IONO_Height, &
                 IONO_Ref_Density, IONO_Ref_SoundSpeed, &
                 IONO_Bdp, IONO_Radius_Mag_Boundary, &
                 IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,*) i,j, &
                       IONO_NORTH_Theta(i,j),IONO_NORTH_Psi(i,j), &
                       IONO_NORTH_JR(i,j)
     end do
  end do

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,*) i,j, &
                       IONO_SOUTH_Theta(i,j),IONO_SOUTH_Psi(i,j), &
                       IONO_SOUTH_JR(i,j) 
     end do
  end do

  close(UNIT=Iunit)

end subroutine ionosphere_write_restart_file

subroutine ionosphere_fac
  !\
  ! Given the solution for the current from the magnetosphere,
  ! this routine maps the solution for the current down to the
  ! ionospheric boundary and evaluates the incoming and outgoing
  ! field-aligned currents (FAC).
  !/
  use ModIonosphere
  implicit none

  integer :: i, j, n, n_near
  real :: Radius, Radius2, Rmax_N, Rmax_S, &
          ConstantB, ConstantJ, temp
  real :: x_test, y_test, z_test, r_test, jr_test
  real, dimension(8) :: x_near, y_near, z_near, r_near, jr_near

  ! Map the magnetospheric field-aligned current (FAC) at the 
  ! inner magnetospheric boundary down to the ionospheric boundary
  ! and then dimensionalize the FAC.

  ConstantB = sqrt(IONO_MU*IONO_Ref_Density* &
                   IONO_Ref_SoundSpeed*IONO_Ref_SoundSpeed)
  ConstantJ = ConstantB/(IONO_MU*IONO_Radius)

  Radius2 = (IONO_Radius + IONO_Height)/IONO_Radius

  ! Northern hemisphere

  Rmax_N = 0.00
  do i = 1, IONO_NORTH_nMagBndPts
     ! Determine Theta, Psi, and R at magnetospheric boundary.

     if (MAG_NORTH_X(i) == 0.00 .and. MAG_NORTH_Y(i) == 0.00) then
        MAG_NORTH_Psi(i) = 0.00
     else
        MAG_NORTH_Psi(i) = atan2(MAG_NORTH_Y(i), MAG_NORTH_X(i))
     end if
     if (MAG_NORTH_Psi(i) < 0.00) MAG_NORTH_Psi(i) = MAG_NORTH_Psi(i) + &
                                                     2.00*IONO_PI
     Radius = sqrt(MAG_NORTH_X(i)*MAG_NORTH_X(i) + &
                   MAG_NORTH_Y(i)*MAG_NORTH_Y(i))
     MAG_NORTH_Theta(i) = atan2(Radius, MAG_NORTH_Z(i))
     Radius = sqrt(MAG_NORTH_X(i)*MAG_NORTH_X(i) + &
                   MAG_NORTH_Y(i)*MAG_NORTH_Y(i) + &
                   MAG_NORTH_Z(i)*MAG_NORTH_Z(i))
     Rmax_N = max(Rmax_N, Radius)
     MAG_NORTH_JR(i) = MAG_NORTH_Jx(i)*sin(MAG_NORTH_Theta(i))*cos(MAG_NORTH_Psi(i)) + &
                       MAG_NORTH_Jy(i)*sin(MAG_NORTH_Theta(i))*sin(MAG_NORTH_Psi(i)) + &
                       MAG_NORTH_Jz(i)*cos(MAG_NORTH_Theta(i))

     ! Map solution to ionospheric boundary.

     temp = sqrt(1.00+3.00*(cos(MAG_NORTH_Theta(i))**2))
     MAG_NORTH_JR(i) = -MAG_NORTH_JR(i)*(temp/(2.00*cos(MAG_NORTH_Theta(i))))

!--- Map total current, not just field aligned component!
!!$     MAG_NORTH_JR(i) = sign(sqrt(MAG_NORTH_Jx(i)*MAG_NORTH_Jx(i)+&
!!$                                 MAG_NORTH_Jy(i)*MAG_NORTH_Jy(i)+&
!!$                                 MAG_NORTH_Jz(i)*MAG_NORTH_Jz(i)), &
!!$                            MAG_NORTH_JR(i))
!---

     MAG_NORTH_Psi(i) = MAG_NORTH_Psi(i)
     MAG_NORTH_Theta(i) = asin(sqrt((Radius2/Radius)* &
                                    (sin(MAG_NORTH_Theta(i))**2)))
     MAG_NORTH_JR(i) = MAG_NORTH_JR(i)*(Radius*Radius*Radius/ &
                                       (Radius2*Radius2*Radius2))/temp
     temp = sqrt(1.00+3.00*(cos(MAG_NORTH_Theta(i))**2))
     MAG_NORTH_JR(i) = MAG_NORTH_JR(i)*temp
     MAG_NORTH_JR(i) = - (2.00*cos(MAG_NORTH_Theta(i))/temp)*MAG_NORTH_JR(i)
  
     ! Dimensionalize the incoming ionospheric FAC.

     MAG_NORTH_JR(i) = MAG_NORTH_JR(i)*ConstantJ

  end do

  ! Southern hemisphere

  Rmax_S = 0.00
  do i = 1, IONO_SOUTH_nMagBndPts
     ! Determine Theta, Psi, and R at magnetospheric boundary.

     if (MAG_SOUTH_X(i) == 0.00 .and. MAG_SOUTH_Y(i) == 0.00) then
        MAG_SOUTH_Psi(i) = 0.00
     else
        MAG_SOUTH_Psi(i) = atan2(MAG_SOUTH_Y(i), MAG_SOUTH_X(i))
     end if
     if (MAG_SOUTH_Psi(i) < 0.00) MAG_SOUTH_Psi(i) = MAG_SOUTH_Psi(i) + &
                                                     2.00*IONO_PI
     Radius = sqrt(MAG_SOUTH_X(i)*MAG_SOUTH_X(i) + &
                   MAG_SOUTH_Y(i)*MAG_SOUTH_Y(i))
     MAG_SOUTH_Theta(i) = atan2(Radius, MAG_SOUTH_Z(i))
     Radius = sqrt(MAG_SOUTH_X(i)*MAG_SOUTH_X(i) + &
                   MAG_SOUTH_Y(i)*MAG_SOUTH_Y(i) + &
                   MAG_SOUTH_Z(i)*MAG_SOUTH_Z(i))
     Rmax_S = max(Rmax_S, Radius)
     MAG_SOUTH_JR(i) = MAG_SOUTH_Jx(i)*sin(MAG_SOUTH_Theta(i))*cos(MAG_SOUTH_Psi(i)) + &
                       MAG_SOUTH_Jy(i)*sin(MAG_SOUTH_Theta(i))*sin(MAG_SOUTH_Psi(i)) + &
                       MAG_SOUTH_Jz(i)*cos(MAG_SOUTH_Theta(i))

     ! Map solution to ionospheric boundary.

     temp = sqrt(1.00+3.00*(cos(MAG_SOUTH_Theta(i))**2))
     MAG_SOUTH_JR(i) = -MAG_SOUTH_JR(i)*(temp/(2.00*cos(MAG_SOUTH_Theta(i))))

!--- Map total current, not just field aligned component!
!!$     MAG_SOUTH_JR(i) = sign(sqrt(MAG_SOUTH_Jx(i)*MAG_SOUTH_Jx(i)+&
!!$                                 MAG_SOUTH_Jy(i)*MAG_SOUTH_Jy(i)+&
!!$                                 MAG_SOUTH_Jz(i)*MAG_SOUTH_Jz(i)), &
!!$                            MAG_SOUTH_JR(i))
!---

     MAG_SOUTH_Psi(i) = MAG_SOUTH_Psi(i)
     MAG_SOUTH_Theta(i) = IONO_PI - asin(sqrt((Radius2/Radius)* &
                                    (sin(MAG_SOUTH_Theta(i))**2)))
     MAG_SOUTH_JR(i) = MAG_SOUTH_JR(i)*(Radius*Radius*Radius/ &
                                       (Radius2*Radius2*Radius2))/temp
     temp = sqrt(1.00+3.00*(cos(MAG_SOUTH_Theta(i))**2))
     MAG_SOUTH_JR(i) = MAG_SOUTH_JR(i)*temp
     MAG_SOUTH_JR(i) = - (2.00*cos(MAG_SOUTH_Theta(i))/temp)*MAG_SOUTH_JR(i)
  
     ! Dimensionalize the incoming ionospheric FAC.

     MAG_SOUTH_JR(i) = MAG_SOUTH_JR(i)*ConstantJ

  end do

  ! For each point in the ionospheric solution, find the five
  ! nearest magnetospheric FAC solution points and interpolate
  ! this solution to the desired ionospheric solution point. 
  ! For points outside of the mapping region, zero the FAC.

  ! Northern hemisphere

  IONO_NORTH_JR = 0.00
  IONO_NORTH_Theta_Max = asin(sqrt(Radius2/Rmax_N))
  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        if (IONO_NORTH_Theta(i,j) <= IONO_NORTH_Theta_Max ) then
           x_near = sin(MAG_NORTH_Theta(1:8))*cos(MAG_NORTH_Psi(1:8))
           y_near = sin(MAG_NORTH_Theta(1:8))*sin(MAG_NORTH_Psi(1:8))
           z_near = cos(MAG_NORTH_Theta(1:8))
           r_near = sqrt((x_near-IONO_NORTH_X(i,j))**2 + &
                         (y_near-IONO_NORTH_Y(i,j))**2 + &
                         (z_near-IONO_NORTH_Z(i,j))**2)
           jr_near = MAG_NORTH_JR(1:8)
           do n = 9, IONO_NORTH_nMagBndPts
              x_test = sin(MAG_NORTH_Theta(n))*cos(MAG_NORTH_Psi(n))
              y_test = sin(MAG_NORTH_Theta(n))*sin(MAG_NORTH_Psi(n))
              z_test = cos(MAG_NORTH_Theta(n))
              r_test = sqrt((x_test-IONO_NORTH_X(i,j))**2 + &
                            (y_test-IONO_NORTH_Y(i,j))**2 + &
                            (z_test-IONO_NORTH_Z(i,j))**2)
              jr_test = MAG_NORTH_JR(n)
              do n_near = 1, 8
                 if (r_test < r_near(n_near)) then
                    temp = r_near(n_near)
                    r_near(n_near) = r_test
                    r_test = temp
                    temp = jr_near(n_near)
                    jr_near(n_near) = jr_test
                    jr_test = temp                 
                 end if
              end do
           end do
           r_near = 1.00/(r_near+1.0e-10)**2
           IONO_NORTH_JR(i,j) = (r_near(1)*jr_near(1) + &
                                 r_near(2)*jr_near(2) + &
                                 r_near(3)*jr_near(3) + &
                                 r_near(4)*jr_near(4) + &
                                 r_near(5)*jr_near(5) + &
                                 r_near(6)*jr_near(6) + &
                                 r_near(7)*jr_near(7) + &
                                 r_near(8)*jr_near(8)) / &
                                (r_near(1) + &
                                 r_near(2) + &
                                 r_near(3) + &
                                 r_near(4) + &
                                 r_near(5) + &
                                 r_near(6) + &
                                 r_near(7) + &
                                 r_near(8))
        end if
     end do
  end do

  IONO_NORTH_Theta_Max = min(IONO_NORTH_Theta_Max+18.00*IONO_PI/180.00, &
                             85.00*IONO_PI/180.00)

  ! Southern hemisphere

  IONO_SOUTH_JR = 0.00
  IONO_SOUTH_Theta_Min = IONO_PI - asin(sqrt(Radius2/Rmax_S))
  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        if (IONO_SOUTH_Theta(i,j) >= IONO_SOUTH_Theta_Min ) then
           x_near = sin(MAG_SOUTH_Theta(1:8))*cos(MAG_SOUTH_Psi(1:8))
           y_near = sin(MAG_SOUTH_Theta(1:8))*sin(MAG_SOUTH_Psi(1:8))
           z_near = cos(MAG_SOUTH_Theta(1:8))
           r_near = sqrt((x_near-IONO_SOUTH_X(i,j))**2 + &
                         (y_near-IONO_SOUTH_Y(i,j))**2 + &
                         (z_near-IONO_SOUTH_Z(i,j))**2)
           jr_near = MAG_SOUTH_JR(1:8)
           do n = 9, IONO_SOUTH_nMagBndPts
              x_test = sin(MAG_SOUTH_Theta(n))*cos(MAG_SOUTH_Psi(n))
              y_test = sin(MAG_SOUTH_Theta(n))*sin(MAG_SOUTH_Psi(n))
              z_test = cos(MAG_SOUTH_Theta(n))
              r_test = sqrt((x_test-IONO_SOUTH_X(i,j))**2 + &
                            (y_test-IONO_SOUTH_Y(i,j))**2 + &
                            (z_test-IONO_SOUTH_Z(i,j))**2)
              jr_test = MAG_SOUTH_JR(n)
              do n_near = 1, 8
                 if (r_test < r_near(n_near)) then
                    temp = r_near(n_near)
                    r_near(n_near) = r_test
                    r_test = temp
                    temp = jr_near(n_near)
                    jr_near(n_near) = jr_test
                    jr_test = temp                 
                 end if
              end do
           end do
           r_near = 1.00/(r_near+1.0e-10)**2
           IONO_SOUTH_JR(i,j) = (r_near(1)*jr_near(1) + &
                                 r_near(2)*jr_near(2) + &
                                 r_near(3)*jr_near(3) + &
                                 r_near(4)*jr_near(4) + &
                                 r_near(5)*jr_near(5) + &
                                 r_near(6)*jr_near(6) + &
                                 r_near(7)*jr_near(7) + &
                                 r_near(8)*jr_near(8)) / &
                                (r_near(1) + &
                                 r_near(2) + &
                                 r_near(3) + &
                                 r_near(4) + &
                                 r_near(5) + &
                                 r_near(6) + &
                                 r_near(7) + &
                                 r_near(8))
        end if
     end do
  end do

  IONO_SOUTH_Theta_Min = max(IONO_SOUTH_Theta_Min-18.00*IONO_PI/180.00, &
                             95.00*IONO_PI/180.00)

  ! Ensure that the FAC is periodic.

  IONO_NORTH_JR(1:IONO_nTheta,IONO_nPsi) = &
     0.50*(IONO_NORTH_JR(1:IONO_nTheta,IONO_nPsi)+ &
           IONO_NORTH_JR(1:IONO_nTheta,1)) 
  IONO_NORTH_JR(1:IONO_nTheta,1) = &
     IONO_NORTH_JR(1:IONO_nTheta,IONO_nPsi)

  IONO_SOUTH_JR(1:IONO_nTheta,IONO_nPsi) = &
     0.50*(IONO_SOUTH_JR(1:IONO_nTheta,IONO_nPsi)+ &
           IONO_SOUTH_JR(1:IONO_nTheta,1)) 
  IONO_SOUTH_JR(1:IONO_nTheta,1) = &
     IONO_SOUTH_JR(1:IONO_nTheta,IONO_nPsi)

end subroutine ionosphere_fac

subroutine ionosphere_currents(Jx,Jy,Jz, &
                               Ex,Ey,Ez,ETh,EPs, &
                               Ux,Uy,Uz, &
                               PHI, SigmaThTh, SigmaThPs, SigmaPsPs, &
                               X, Y, Z, &
                               Theta, Psi, Radius, nTheta, nPsi)
  !\
  ! For the calculated ionospheric potential solution,
  ! this routine determines the ionospheric currents and
  ! electric fields, as well as convection velocities.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI, SigmaThTh, SigmaThPs, SigmaPsPs, &
                  Jx,Jy,Jz, &
                  Ex,Ey,Ez,ETh,EPs, &
                  Ux,Uy,Uz, &
                  X, Y, Z, &
                  Theta, Psi

  logical :: north
  integer :: i, j
  real :: dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, cosTheta, sinTheta, cosPhi, sinPhi, &
          cs2, cs3, cs4, &
          ER, JR, JTh, JPs, &
          xx, yy, zz, &
          bx, by, bz, bR, bTh, bPs, BB, BBabs, &
          Vll, Vp_x, Vp_y, Vp_z, VR

  if (Theta(1,1) < 2.00*IONO_Theta_0) then
     north = .true.
  else
     north = .false.
  end if

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  ! Compute the ionospheric electric field.

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/ &
                      (2.00*dd*Radius)
           EPs(i,j) = -(PHI(i,j+1)-PHI(i,j-1))/ &
                      (2.00*dd*Radius*sinTheta)
        end do
        ETh(1,j) = -(PHI(2,j)-PHI(1,j))/ &
                   (dd*Radius)
        EPs(1,j) = EPs(2,j)
        ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/ &
                        (dd*Radius)
        EPs(nTheta,j) = EPs(nTheta-1,j)
     else if (j == 1) then
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/ &
                      (2.00*dd*Radius)
           EPs(i,j) = -(PHI(i,j+1)-PHI(i,nPsi-1))/ &
                      (2.00*dd*Radius*sinTheta)
        end do
        ETh(1,j) = -(PHI(2,j)-PHI(1,j))/ &
                   (dd*Radius)
        EPs(1,j) = EPs(2,j)
        ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/ &
                        (dd*Radius)
        EPs(nTheta,j) = EPs(nTheta-1,j)
     else
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/ &
                      (2.00*dd*Radius)
           EPs(i,j) = -(PHI(i,2)-PHI(i,j-1))/ &
                      (2.00*dd*Radius*sinTheta)
        end do
        ETh(1,j) = -(PHI(2,j)-PHI(1,j))/ &
                   (dd*Radius)
        EPs(1,j) = EPs(2,j)
        ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/ &
                        (dd*Radius)
        EPs(nTheta,j) = EPs(nTheta-1,j)
     end if
  end do

  ! Compute the ionospheric currents convection velocities.

  do j = 1, nPsi
     do i = 1, nTheta
        cosTheta = cos(Theta(i,j))
        sinTheta = sin(Theta(i,j))
        cosPhi = cos(Psi(i,j))
        sinPhi = sin(Psi(i,j))

        if (north .and. i == nTheta) then
           ER = 0.00
        else if (.not.north .and. i == 1) then
           ER = 0.00
        else
           ER = -0.50*(sinTheta/(cosTheta+IONO_Toler**2))*ETh(i,j)
        end if

        Ex(i,j) = ER*sinTheta*cosPhi + ETh(i,j)*cosTheta*cosPhi - &
                  EPs(i,j)*sinPhi
        Ey(i,j) = ER*sinTheta*sinPhi + ETh(i,j)*cosTheta*sinPhi + &
                  EPs(i,j)*cosPhi
        Ez(i,j) = ER*cosTheta - ETh(i,j)*sinTheta
        
        JR = 0.00
        JTh =  SigmaThTh(i,j)*ETh(i,j) + SigmaThPs(i,j)*EPs(i,j)
        JPs = -SigmaThPs(i,j)*ETh(i,j) + SigmaPsPs(i,j)*EPs(i,j)
        
        Jx(i,j) = JR*sinTheta*cosPhi + JTh*cosTheta*cosPhi - &
                  JPs*sinPhi
        Jy(i,j) = JR*sinTheta*sinPhi + JTh*cosTheta*sinPhi + &
                  JPs*cosPhi
        Jz(i,j) = JR*cosTheta - JTh*sinTheta
        
        cs2 = cosTheta*cosTheta
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)

        BB = IONO_Bdp*((IONO_Radius/Radius)**3)
        BBabs = abs(BB)
        bR = (BB/BBabs)*2.00*cosTheta/cs4
        bTh = (BB/BBabs)*sinTheta/cs4 
        bPs = 0.00  
        bx = bR*sinTheta*cosPhi + bTh*cosTheta*cosPhi - &
             bPs*sinPhi
        by = bR*sinTheta*sinPhi + bTh*cosTheta*sinPhi + &
             bPs*cosPhi
        bz = bR*cosTheta - bTh*sinTheta
        
        xx = sinTheta*cosPhi
        yy = sinTheta*sinPhi
        zz = cosTheta

        Vp_x = (Ey(i,j)*bz - Ez(i,j)*by)/BBabs
        Vp_y = (Ez(i,j)*bx - Ex(i,j)*bz)/BBabs
        Vp_z = (Ex(i,j)*by - Ey(i,j)*bx)/BBabs

!!$        VR = Vp_x*xx + Vp_y*yy + Vp_z*zz
!!$        if (north .and. i == nTheta) then
!!$           Vll = 0.00
!!$        else if (.not.north .and. i == 1) then
!!$           Vll = 0.00
!!$        else
!!$           Vll = -VR/(bx*xx+by*yy+bz*zz)
!!$        end if

        Vll = 0.00
        
        Ux(i,j) = Vp_x + Vll*bx
        Uy(i,j) = Vp_y + Vll*by
        Uz(i,j) = Vp_z + Vll*bz
     end do  
  end do

end subroutine ionosphere_currents

subroutine ionosphere_magBCs(PHI_BC, ETh_BC, EPs_BC, &
                             UR_BC, UTh_BC, UPs_BC, Radius_BC, &
                             PHI, X, Y, Z, &
                             Theta, Psi, Radius, nTheta, nPsi)
  !\
  ! For the calculated ionospheric potential solution,
  ! map the potential solution out to the magnetospheric
  ! inner boundary at R=Radius_BC and determine the convection
  ! velocities to be used in the magnetospheric boundary
  ! conditions.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi
  real :: Radius, Radius_BC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI_BC, ETh_BC, EPs_BC, &
                  UR_BC, UTh_BC, UPs_BC, &
                  PHI, X, Y, Z, Theta, Psi

  logical :: north
  integer :: i, j, i0
  real :: Theta0, dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, cosTheta, sinTheta, cosPhi, sinPhi, &
          cs2, cs3, cs4, xx, yy, zz, &
          ER, Ex, Ey, Ez, Ux, Uy, Uz, &
          bx, by, bz, bR, bTh, bPs, BB, BBabs, &
          Vll, Vp_x, Vp_y, Vp_z, VR

  if (Theta(1,1) < 2.00*IONO_Theta_0) then
     north = .true.
  else
     north = .false.
  end if

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  ! Determine the potential at the magnetospheric inner boundary.

  if (Radius_BC/IONO_Radius < 1.05) then
     PHI_BC = PHI
  else if (north) then
     do j = 1, nPsi
        PHI_BC(1,j) = PHI(1,j)
        do i = 2, nTheta
           Theta0 = asin(sqrt((Radius/Radius_BC)* &
                              (sin(Theta(i,j))**2)))
	   i0 = Theta0/dd + 1
	   PHI_BC(i,j) = PHI(i0,j) + &
                         (Theta0-Theta(i0,j))* &
                         (PHI(i0+1,j)-PHI(i0,j))/ &
                         dd
        end do
     end do
  else
     do j = 1, nPsi
        PHI_BC(nTheta,j) = PHI(nTheta,j)
        do i = 1, nTheta-1
           Theta0 = IONO_PI - asin(sqrt((Radius/Radius_BC)* &
                                        (sin(Theta(i,j))**2)))
	   i0 = (Theta0-0.50*IONO_PI)/dd + 1
	   PHI_BC(i,j) = PHI(i0,j) + &
                         (Theta0-Theta(i0,j))* &
                         (PHI(i0+1,j)-PHI(i0,j))/ &
                         dd
        end do
     end do
  end if

  ! Compute the electric field at the magnetospheric inner boundary.

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,j-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     else if (j == 1) then
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,nPsi-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     else
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,2)-PHI_BC(i,j-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     end if
  end do

  ! Compute the convection velocities at the magnetospheric inner boundary.

  do j = 1, nPsi
     do i = 1, nTheta
        cosTheta = cos(Theta(i,j))
        sinTheta = sin(Theta(i,j))
        cosPhi = cos(Psi(i,j))
        sinPhi = sin(Psi(i,j))

        if (north .and. i == nTheta) then
           ER = 0.00
        else if (.not.north .and. i == 1) then
           ER = 0.00
        else
           ER = -0.50*(sinTheta/(cosTheta+IONO_Toler**2))*ETh_BC(i,j)
        end if
      
        Ex = ER*sinTheta*cosPhi + ETh_BC(i,j)*cosTheta*cosPhi - &
             EPs_BC(i,j)*sinPhi
        Ey = ER*sinTheta*sinPhi + ETh_BC(i,j)*cosTheta*sinPhi + &
             EPs_BC(i,j)*cosPhi
        Ez = ER*cosTheta - ETh_BC(i,j)*sinTheta
        
        cs2 = cosTheta*cosTheta
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)

        BB = IONO_Bdp*((IONO_Radius/Radius_BC)**3)
        BBabs = abs(BB)
        bR = (BB/BBabs)*2.00*cosTheta/cs4
        bTh = (BB/BBabs)*sinTheta/cs4 
        bPs = 0.00  
        bx = bR*sinTheta*cosPhi + bTh*cosTheta*cosPhi - &
             bPs*sinPhi
        by = bR*sinTheta*sinPhi + bTh*cosTheta*sinPhi + &
             bPs*cosPhi
        bz = bR*cosTheta - bTh*sinTheta
        
        Vp_x = (Ey*bz - Ez*by)/BBabs
        Vp_y = (Ez*bx - Ex*bz)/BBabs
        Vp_z = (Ex*by - Ey*bx)/BBabs
        
        xx = sinTheta*cosPhi
        yy = sinTheta*sinPhi
        zz = cosTheta

!!$        VR = Vp_x*xx + Vp_y*yy + Vp_z*zz
!!$        if (north .and. i == nTheta) then
!!$           Vll = 0.00
!!$        else if (.not.north .and. i == 1) then
!!$           Vll = 0.00
!!$        else
!!$           Vll = -VR/(bx*xx+by*yy+bz*zz)
!!$        end if

        Vll = 0.00
        
        Ux = Vp_x + Vll*bx
        Uy = Vp_y + Vll*by
        Uz = Vp_z + Vll*bz

        UR_BC(i,j) = Ux*xx + Uy*yy + Uz*zz
        UTh_BC(i,j) = ((Ux*xx + Uy*yy)*zz - &
                       Uz*(xx**2+yy**2)) / &
                      sqrt(xx**2+yy**2+IONO_TOLER**2)
        UPs_BC(i,j) = (Uy*xx - Ux*yy)*sinTheta / &
                      (xx**2+yy**2+IONO_TOLER**2)

        UR_BC(i,j) = UR_BC(i,j) / IONO_Ref_SoundSpeed
        UTh_BC(i,j) = UTh_BC(i,j) / IONO_Ref_SoundSpeed
        UPs_BC(i,j) = UPs_BC(i,j) / IONO_Ref_SoundSpeed
     end do  
  end do

end subroutine ionosphere_magBCs

subroutine ionosphere_multigrid(PHI, &
                                SigmaThTh, SigmaThPs, SigmaPsPs, &
                                dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                                dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                                Theta, Psi, Radius, nTheta, nPsi, &
                                ncycle, iboundary)
  !\
  ! This subroutine applies a full linear multigrid solution
  ! algorithm to the problem of determining the solution of the
  ! boundary value problem (BVP) for the linear elliptic equation
  ! that governs planetary ionospheric currents on a spherical
  ! surface of radius Radius.  The linear elliptic PDE for the
  ! electric field potential PHI is defined on the domain 
  ! 0 < Theta < PI/2 (this is for the northern hemisphere, 
  ! PI/2 < Theta < PI for southern hemisphere), 0 < Psi < 2 PI, 
  ! and subject to the Dirichlet boundary data
  !
  !      PHI(0,Psi) = PHI(PI,Psi) = 0,
  !
  !      PHI(PI/2,Psi) = 0,
  !
  !      PHI(Theta,0) = PHI(Theta,2*PI).
  ! 
  ! A red-black Gauss-Seidel scheme is used as the smoothing
  ! operator, bilinear interpolation is used for the prolongation
  ! operator, and half-weighting is used for the restriction
  ! operator.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi, ncycle, iboundary
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI, &
                  SigmaThTh, SigmaThPs, SigmaPsPs, &
                  dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                  dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                  Theta, Psi

  integer :: nThetaC, nPsiC, nThetaF, nPsiF
  integer :: i, j, jj, jcycle, jpre, jpost, irepeat, ngrid
  integer :: i2, j2
  real :: RESIDUAL

  ! Multiply the right-hand-side source terms
  ! (in effect the field-aligned current) by
  ! Radius^2 sin^2(Theta).
  
  do j = 1, nPsi
     do i = 2, nTheta-1
        PHI(i,j) = PHI(i,j)*Radius*Radius* &
                   sin(Theta(i,j))*sin(Theta(i,j))
     end do
     PHI(1,j) = 0.00
     PHI(nTheta,j) = 0.00
  end do

  ! Create coarse grid and fill
  ! rhs by restricting rhs from the fine grid. 
  ! Also determine the field aligned and Hall and
  ! Pedersen conductances.
  
  ngrid = IONO_NGMAX

  nThetaC=nTheta
  nPsiC=nPsi
  IONO_X(:,:,ngrid) = Theta
  IONO_Y(:,:,ngrid) = Psi
  IONO_RHO(:,:,ngrid) = PHI
  call ionosphere_conductance(IONO_S0(:,:,ngrid),IONO_SH(:,:,ngrid),IONO_SP(:,:,ngrid), &
                              IONO_Sxx(:,:,ngrid),IONO_Sxy(:,:,ngrid),IONO_Syy(:,:,ngrid), &
                              IONO_dSxxdx(:,:,ngrid),IONO_dSxydx(:,:,ngrid), &
                              IONO_dSyydx(:,:,ngrid), &
                              IONO_dSxxdy(:,:,ngrid),IONO_dSxydy(:,:,ngrid), &
                              IONO_dSyydy(:,:,ngrid), &
                              IONO_X(:,:,ngrid), IONO_Y(:,:,ngrid), &
                              nThetaC, nPsiC)
  do j = ngrid-1, 1, -1
     nThetaF = nThetaC
     nPsiF = nPsiC
     nThetaC = nThetaF/2+1
     nPsiC = nPsiF/2+1
     call ionosphere_coarse_grid(IONO_X(:,:,j+1), IONO_Y(:,:,j+1), nThetaF, nPsiF, &
                                 IONO_X(:,:,j), IONO_Y(:,:,j), nThetaC, nPsiC)
     call ionosphere_restrict(IONO_RHO(:,:,j+1), nThetaF, nPsiF, &
                              IONO_RHO(:,:,j), nThetaC, nPsiC)
     call ionosphere_conductance(IONO_S0(:,:,j),IONO_SH(:,:,j),IONO_SP(:,:,j), &
                                 IONO_Sxx(:,:,j),IONO_Sxy(:,:,j),IONO_Syy(:,:,j), &
                                 IONO_dSxxdx(:,:,j),IONO_dSxydx(:,:,j), &
                                 IONO_dSyydx(:,:,j), &
                                 IONO_dSxxdy(:,:,j),IONO_dSxydy(:,:,j), &
                                 IONO_dSyydy(:,:,j), &
                                 IONO_X(:,:,j), IONO_Y(:,:,j), &
                                 nThetaC, nPsiC)
  end do

  ! Initial solution on coarsest grid.  Solve this
  ! problem using IONO_NEXACT Gauss-Seidel interations. 

  IONO_U(:,:,1) = 0.00
  do j = 1, IONO_NEXACT
     call ionosphere_smoother(IONO_U(:,:,1), IONO_RHO(:,:,1), &
                              IONO_Sxx(:,:,1),IONO_Sxy(:,:,1),IONO_Syy(:,:,1), &
                              IONO_dSxxdx(:,:,1),IONO_dSxydx(:,:,1), &
                              IONO_dSyydx(:,:,1), &
                              IONO_dSxxdy(:,:,1),IONO_dSxydy(:,:,1), &
                              IONO_dSyydy(:,:,1), &
                              IONO_X(:,:,1), IONO_Y(:,:,1), &
                              Radius, nThetaC, nPsiC, &
                              iboundary)
  end do

  ! Nested FMG iteration loop.

  nThetaF = nThetaC
  nPsiF = nPsiC
  do j = 2, ngrid
     nThetaC = nThetaF
     nPsiC = nPsiF  
     nThetaF = 2*nThetaC-1
     nPsiF = 2*nPsiC-1

     ! Prolong from coarse grid to next
     ! finer grid using bilinear interpolation.
         
     call ionosphere_prolong(IONO_U(:,:,j), nThetaF, nPsiF, &
                             IONO_U(:,:,j-1), nThetaC, nPsiC)
  
     ! Set up rhs.
         
     IONO_RHS(:,:,j) = IONO_RHO(:,:,j)  
    
     ! V-cycle loop.
  
     irepeat = 0
     do jcycle = 1, ncycle
        ! Begin V-cycle.
1000    continue
  
        nThetaC=nThetaF
        nPsiC=nPsiF

        do jj = j, 2, -1
           ! Pre-smoothing on downward stroke of
           ! V-cycle.

           do  jpre = 1, IONO_NPRE
              call ionosphere_smoother(IONO_U(:,:,jj), IONO_RHS(:,:,jj), &
                                       IONO_Sxx(:,:,jj),IONO_Sxy(:,:,jj),IONO_Syy(:,:,jj), &
                                       IONO_dSxxdx(:,:,jj),IONO_dSxydx(:,:,jj), &
                                       IONO_dSyydx(:,:,jj), &
                                       IONO_dSxxdy(:,:,jj),IONO_dSxydy(:,:,jj), &
                                       IONO_dSyydy(:,:,jj), &
                                       IONO_X(:,:,jj), IONO_Y(:,:,jj), &
                                       Radius, nThetaC, nPsiC, &
                                       iboundary)
           end do
        
           ! Restricted residual is the next rhs.
               
           call ionosphere_residual(IONO_RES(:,:,jj), IONO_U(:,:,jj), IONO_RHS(:,:,jj), &
                                    IONO_Sxx(:,:,jj),IONO_Sxy(:,:,jj),IONO_Syy(:,:,jj), &
                                    IONO_dSxxdx(:,:,jj),IONO_dSxydx(:,:,jj), &
                                    IONO_dSyydx(:,:,jj), &
                                    IONO_dSxxdy(:,:,jj),IONO_dSxydy(:,:,jj), &
                                    IONO_dSyydy(:,:,jj), &
                                    IONO_X(:,:,jj), IONO_Y(:,:,jj), &
                                    Radius, nThetaC, nPsiC, &
                                    iboundary, RESIDUAL)

           call ionosphere_restrict(IONO_RES(:,:,jj), nThetaC, nPsiC, &
                                    IONO_RHS(:,:,jj-1), nThetaC/2+1, nPsiC/2+1)

           nThetaC = nThetaC/2+1
           nPsiC = nPsiC/2+1

           ! Zero U for initial guess in next relaxation.
                
           IONO_U(:,:,jj-1) = 0.00

        end do ! Downward V-cycle.
  
        ! Bottom of V-cycle; solve on coarsest grid.
                
        do jj = 1, IONO_NEXACT
           call ionosphere_smoother(IONO_U(:,:,1), IONO_RHS(:,:,1), &
                                    IONO_Sxx(:,:,1),IONO_Sxy(:,:,1),IONO_Syy(:,:,1), &
                                    IONO_dSxxdx(:,:,1),IONO_dSxydx(:,:,1), &
                                    IONO_dSyydx(:,:,1), &
                                    IONO_dSxxdy(:,:,1),IONO_dSxydy(:,:,1), &
                                    IONO_dSyydy(:,:,1), &
                                    IONO_X(:,:,1), IONO_Y(:,:,1), &
                                    Radius, nThetaC, nPsiC, &
                                    iboundary)
        end do
      
        ! Go back up.
  
        do jj = 2, j
           ! Prolong solution corrections from coarse to finer
           ! grid.
        
           call ionosphere_prolong(IONO_RES(:,:,jj), 2*nThetaC-1, 2*nPsiC-1, &
                                   IONO_U(:,:,jj-1), nThetaC, nPsiC)

           nThetaC = 2*nThetaC-1
           nPsiC = 2*nPsiC-1
     
           ! Update the solution on the finer grid.
        
           IONO_U(:,:,jj) = IONO_U(:,:,jj)+IONO_OMEGA*IONO_RES(:,:,jj)
        
           ! Post-smoothing on upward stroke of V-cycle.
                  
           do  jpost = 1, IONO_NPOST
              call ionosphere_smoother(IONO_U(:,:,jj), IONO_RHS(:,:,jj), &
                                       IONO_Sxx(:,:,jj),IONO_Sxy(:,:,jj),IONO_Syy(:,:,jj), &
                                       IONO_dSxxdx(:,:,jj),IONO_dSxydx(:,:,jj), &
                                       IONO_dSyydx(:,:,jj), &
                                       IONO_dSxxdy(:,:,jj),IONO_dSxydy(:,:,jj), &
                                       IONO_dSyydy(:,:,jj), &
                                       IONO_X(:,:,jj), IONO_Y(:,:,jj), &
                                       Radius, nThetaC, nPsiC, &
                                       iboundary)
           end do
  
           ! Determine the residuals on the finer grid.
  
           call ionosphere_residual(IONO_RES(:,:,jj), IONO_U(:,:,jj), IONO_RHS(:,:,jj), &
                                    IONO_Sxx(:,:,jj),IONO_Sxy(:,:,jj),IONO_Syy(:,:,jj), &
                                    IONO_dSxxdx(:,:,jj),IONO_dSxydx(:,:,jj), &
                                    IONO_dSyydx(:,:,jj), &
                                    IONO_dSxxdy(:,:,jj),IONO_dSxydy(:,:,jj), &
                                    IONO_dSyydy(:,:,jj), &
                                    IONO_X(:,:,jj), IONO_Y(:,:,jj), &
                                    Radius, nThetaC, nPsiC, &
                                    iboundary, RESIDUAL)

!!$           write(*,*) "IONOSPHERE_MULTIGRID: Grid = ",nThetaC-1," x ",nPsiC-1,RESIDUAL
  
        end do ! Upward V-cycle.

	! Check for convergence at this grid level
        ! and continue futher V-cycles as necessary.
        
        if (j == ngrid .and. jcycle == ncycle .and. &
            RESIDUAL > IONO_TOLER .and. irepeat <= 4) then
           irepeat = irepeat + 1
           goto 1000
        end if
          
     end do ! (V-cyle).
  
  end do ! (FMG).
  
  ! Return solution in PHI.
  
  PHI = IONO_U(:,:,ngrid)
  
end subroutine ionosphere_multigrid

subroutine ionosphere_conductance(Sigma0, SigmaH, SigmaP, &
                                  SigmaThTh, SigmaThPs, SigmaPsPs, &
                                  dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                                  dSigmaPsPs_dTheta, &
                                  dSigmaThTh_dPsi, dSigmaThPs_dPsi, &
                                  dSigmaPsPs_dPsi, &
                                  Theta, Psi, nTheta, nPsi)
  !\
  ! This subroutine computes the height-integrated field-aligned and
  ! Hall and Pedersen conductances for the ionosphere at each
  ! location of the discretized solution domain.  The gradients of
  ! these quantities are also computed.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta,nPsi
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                Sigma0, SigmaH, SigmaP, &
                SigmaThTh, SigmaThPs, SigmaPsPs, &
                dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                Theta, Psi

  integer :: i,j
  real :: dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, sn, cs, sn2, cs2, cs3, cs4, C
  real :: cos_SZA, SigmaH_EUV, SigmaP_EUV, SigmaH_SCAT, SigmaP_SCAT, &
          SigmaH_STAR, SigmaP_STAR
  

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  do j = 1, nPsi
     do i = 1, nTheta
!!$        !\
!!$        ! Southward IMF Values.
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 8.00
!!$        SigmaP(i,j) = 6.00

!!$        !\
!!$        ! Northward IMF Values.
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 1.25
!!$        SigmaP(i,j) = 1.00

!!$        !\
!!$        ! Saturn Values.
!!$        !/
!!$        Sigma0(i,j) = 10000.00
!!$        SigmaH(i,j) = 0.00
!!$        SigmaP(i,j) = 10.00
        
!!$        !\
!!$        ! GGCM Values.
!!$        !/
!!$        Sigma0(i,j) = 5000.00
!!$        SigmaH(i,j) = 0.00
!!$        SigmaP(i,j) = 5.00

        !\
        ! Test Case.
        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 0.00
!!$        SigmaP(i,j) = 1.00

        Sigma0(i,j) = 1000.00
        cos_SZA = sin(Theta(i,j))*cos(Psi(i,j))
        SigmaH_EUV = 3.20*sqrt(max(0.00,cos_SZA))
        SigmaP_EUV = 1.60*sqrt(max(0.00,cos_SZA))
        if (cos_SZA > 0) then
           SigmaH_SCAT = 1.00
           SigmaP_SCAT = 0.50
        else
           SigmaH_SCAT = 1.00*(10.00**cos_SZA)
           SigmaP_SCAT = 0.50*(10.00**cos_SZA)
        end if
        SigmaH_STAR = 0.40
        SigmaP_STAR = 0.20
        SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + &
                           SigmaH_SCAT*SigmaH_SCAT + &
                           SigmaH_STAR*SigmaH_STAR)
        SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV + &
                           SigmaP_SCAT*SigmaP_SCAT + &
                           SigmaP_STAR*SigmaP_STAR)

!!$        !\
!!$        ! Rasmussen and Schunk model B (JGR, Vol. 92, pp. 4491-4504, 1987).
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        cos_SZA = sin(Theta(i,j))*cos(Psi(i,j))
!!$        SigmaH_EUV = 10.00*sqrt(max(0.00,cos_SZA))
!!$        SigmaP_EUV = 5.00*sqrt(max(0.00,cos_SZA))
!!$        if (cos_SZA > 0) then
!!$           SigmaH_SCAT = 1.00
!!$           SigmaP_SCAT = 0.50
!!$        else
!!$           SigmaH_SCAT = 1.00*(10.00**cos_SZA)
!!$           SigmaP_SCAT = 0.50*(10.00**cos_SZA)
!!$        end if
!!$        SigmaH_STAR = 0.10
!!$        SigmaP_STAR = 0.05
!!$        SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + &
!!$                           SigmaH_SCAT*SigmaH_SCAT + &
!!$                           SigmaH_STAR*SigmaH_STAR)
!!$        SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV + &
!!$                           SigmaP_SCAT*SigmaP_SCAT + &
!!$                           SigmaP_STAR*SigmaP_STAR)
        
        sn = sin(Theta(i,j))
        cs = cos(Theta(i,j))
        sn2= sn*sn
        cs2 = cs*cs
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)
        C = 4.00*Sigma0(i,j)*cs2 + &
            SigmaP(i,j)*sn2
        
        SigmaThTh(i,j) = Sigma0(i,j)*SigmaP(i,j)*cs3/C
        SigmaThPs(i,j) = 2.00*Sigma0(i,j)*SigmaH(i,j)* &
                         cs*cs4/C
        SigmaPsPs(i,j) = SigmaP(i,j)+ &
                         SigmaH(i,j)*SigmaH(i,j)* &
                         sn2/C
        
!!$        SigmaThTh(i,j) = SigmaP(i,j)*cs3/(4.00*cs2)
!!$        SigmaThPs(i,j) = SigmaH(i,j)*cs4/(2.00*cs)
!!$        SigmaPsPs(i,j) = SigmaP(i,j)
        
        dSigmaThTh_dTheta(i,j) = 0.00
        dSigmaThTh_dPsi(i,j) = 0.00
        dSigmaThPs_dTheta(i,j) = 0.00
        dSigmaThPs_dPsi(i,j) = 0.00
        dSigmaPsPs_dTheta(i,j) = 0.00
        dSigmaPsPs_dPsi(i,j) = 0.00
     end do  
  end do

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/ &
                                  (2.00*dd)
        end do
     else if (j == 1) then
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nPsi-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nPsi-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nPsi-1))/ &
                                  (2.00*dd)
        end do
     else
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,2)-SigmaThTh(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,2)-SigmaThPs(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,2)-SigmaPsPs(i,j-1))/ &
                                  (2.00*dd)
        end do
     end if
  end do

end subroutine ionosphere_conductance

subroutine ionosphere_hall_conductance(Sigma0, SigmaH, SigmaP, &
                                       SigmaThTh, SigmaThPs, SigmaPsPs, &
                                       dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                                       dSigmaPsPs_dTheta, &
                                       dSigmaThTh_dPsi, dSigmaThPs_dPsi, &
                                       dSigmaPsPs_dPsi, &
                                       Theta, Psi, nTheta, nPsi)
  !\
  ! This subroutine computes the height-integrated field-aligned and
  ! Hall conductances for the ionosphere at each location of the 
  ! discretized solution domain.  The Pedersen conductances are taken
  ! to be zero.  The gradients of these quantities are also computed.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta,nPsi
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                Sigma0, SigmaH, SigmaP, &
                SigmaThTh, SigmaThPs, SigmaPsPs, &
                dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                Theta, Psi

  integer :: i,j
  real :: dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, sn, cs, sn2, cs2, cs3, cs4, C

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  do j = 1, nPsi
     do i = 1, nTheta
!!$        !\
!!$        ! Southward IMF Values.
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 8.00
!!$        SigmaP(i,j) = 0.00

!!$        !\
!!$        ! Northward IMF Values.
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 1.25
!!$        SigmaP(i,j) = 0.00
        
        !\
        ! GGCM Values.
        !/
        Sigma0(i,j) = 5000.00
        SigmaH(i,j) = 0.00
        SigmaP(i,j) = 0.00

        sn = sin(Theta(i,j))
        cs = cos(Theta(i,j))
        sn2= sn*sn
        cs2 = cs*cs
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)
        C = 4.00*Sigma0(i,j)*cs2 + &
            SigmaP(i,j)*sn2
        
        SigmaThTh(i,j) = Sigma0(i,j)*SigmaP(i,j)*cs3/C
        SigmaThPs(i,j) = 2.00*Sigma0(i,j)*SigmaH(i,j)* &
                         cs*cs4/C
        SigmaPsPs(i,j) = SigmaP(i,j)+ &
                         SigmaH(i,j)*SigmaH(i,j)* &
                         sn2/C
        
        dSigmaThTh_dTheta(i,j) = 0.00
        dSigmaThTh_dPsi(i,j) = 0.00
        dSigmaThPs_dTheta(i,j) = 0.00
        dSigmaThPs_dPsi(i,j) = 0.00
        dSigmaPsPs_dTheta(i,j) = 0.00
        dSigmaPsPs_dPsi(i,j) = 0.00
     end do  
  end do

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/ &
                                  (2.00*dd)
        end do
     else if (j == 1) then
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nPsi-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nPsi-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nPsi-1))/ &
                                  (2.00*dd)
        end do
     else
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,2)-SigmaThTh(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,2)-SigmaThPs(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,2)-SigmaPsPs(i,j-1))/ &
                                  (2.00*dd)
        end do
     end if
  end do

end subroutine ionosphere_hall_conductance

subroutine ionosphere_ped_conductance(Sigma0, SigmaH, SigmaP, &
                                      SigmaThTh, SigmaThPs, SigmaPsPs, &
                                      dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                                      dSigmaPsPs_dTheta, &
                                      dSigmaThTh_dPsi, dSigmaThPs_dPsi, &
                                      dSigmaPsPs_dPsi, &
                                      Theta, Psi, nTheta, nPsi)
  !\
  ! This subroutine computes the height-integrated field-aligned and
  ! Pedersen conductances for the ionosphere at each location of the 
  ! discretized solution domain.  The Hall conductances are taken
  ! to be zero.  The gradients of these quantities are also computed.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta,nPsi
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                Sigma0, SigmaH, SigmaP, &
                SigmaThTh, SigmaThPs, SigmaPsPs, &
                dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                Theta, Psi

  integer :: i,j
  real :: dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, sn, cs, sn2, cs2, cs3, cs4, C

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  do j = 1, nPsi
     do i = 1, nTheta
!!$        !\
!!$        ! Southward IMF Values.
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 0.00
!!$        SigmaP(i,j) = 6.00

!!$        !\
!!$        ! Northward IMF Values.
!!$        !/
!!$        Sigma0(i,j) = 1000.00
!!$        SigmaH(i,j) = 0.00
!!$        SigmaP(i,j) = 1.00

        !\
        ! GGCM Values.
        !/
        Sigma0(i,j) = 5000.00
        SigmaH(i,j) = 0.00
        SigmaP(i,j) = 5.00

        sn = sin(Theta(i,j))
        cs = cos(Theta(i,j))
        sn2= sn*sn
        cs2 = cs*cs
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)
        C = 4.00*Sigma0(i,j)*cs2 + &
            SigmaP(i,j)*sn2
        
        SigmaThTh(i,j) = Sigma0(i,j)*SigmaP(i,j)*cs3/C
        SigmaThPs(i,j) = 2.00*Sigma0(i,j)*SigmaH(i,j)* &
                         cs*cs4/C
        SigmaPsPs(i,j) = SigmaP(i,j)+ &
                         SigmaH(i,j)*SigmaH(i,j)* &
                         sn2/C
        
        dSigmaThTh_dTheta(i,j) = 0.00
        dSigmaThTh_dPsi(i,j) = 0.00
        dSigmaThPs_dTheta(i,j) = 0.00
        dSigmaThPs_dPsi(i,j) = 0.00
        dSigmaPsPs_dTheta(i,j) = 0.00
        dSigmaPsPs_dPsi(i,j) = 0.00
     end do  
  end do

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/ &
                                  (2.00*dd)
        end do
     else if (j == 1) then
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nPsi-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nPsi-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nPsi-1))/ &
                                  (2.00*dd)
        end do
     else
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,2)-SigmaThTh(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,2)-SigmaThPs(i,j-1))/ &
                                  (2.00*dd)
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (2.00*dd)
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,2)-SigmaPsPs(i,j-1))/ &
                                  (2.00*dd)
        end do
     end if
  end do

end subroutine ionosphere_ped_conductance

subroutine ionosphere_coarse_grid(ThetaF, PsiF, nThetaF, nPsiF, &
                                  ThetaC, PsiC, nThetaC, nPsiC)
  !\
  ! Given a fine grid, this routine determines the corresponding
  ! coarse grid.
  !/
  use ModIonosphere
  implicit none

  integer :: nThetaF,nPsiF, &
             nThetaC,nPsiC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                ThetaF, PsiF
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                ThetaC, PsiC

  integer :: ic,iif,jc,jf

  ! Interior points

  jf = 1
  do jc = 2, nPsiC-1
     jf = jf + 2
     iif = 1
     do ic = 2, nThetaC-1  
        iif = iif + 2
        ThetaC(ic,jc)=ThetaF(iif,jf)
        PsiC(ic,jc)=PsiF(iif,jf)
     end do
  end do
        
  ! Boundary points

  iif = -1
  do ic = 1, nThetaC
    iif = iif + 2
    ThetaC(ic,1) = ThetaF(iif,1)
    PsiC(ic,1) = PsiF(iif,1)
    ThetaC(ic,nPsiC) = ThetaF(iif,2*nPsiC-1)
    PsiC(ic,nPsiC) = PsiF(iif,2*nPsiC-1)
  end do
        
  jf = -1
  do jc = 1, nPsiC
    jf = jf + 2
    ThetaC(1,jc) = ThetaF(1,jf)
    PsiC(1,jc) = PsiF(1,jf)
    ThetaC(nThetaC,jc) = ThetaF(2*nThetaC-1,jf)
    PsiC(nThetaC,jc) = PsiF(2*nThetaC-1,jf)
  end do

end subroutine ionosphere_coarse_grid

subroutine ionosphere_restrict(UF, nThetaF, nPsiF, &
                               UC, nThetaC, nPsiC)
  !\
  ! This subroutine uses the half-weighting restriction operator
  ! to restrict the fine grid solution to the coarse grid.
  !/
  use ModIonosphere
  implicit none

  integer :: nThetaF,nPsiF, &
             nThetaC,nPsiC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) :: UF
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) :: UC

  integer :: ic,iif,jc,jf

  ! Interior points.

  jf = 3
  do jc = 2, nPsiC-1
     iif = 3
     do ic = 2, nThetaC-1  
        UC(ic,jc) =   0.500*UF(iif,jf) &
                    + 0.125*( UF(iif+1,jf)+UF(iif-1,jf) &
                    +         UF(iif,jf+1)+UF(iif,jf-1) )
        iif = iif + 2
     end do
     jf = jf + 2
  end do
        
  ! Boundary points.

  iif = 3
  do ic = 2, nThetaC-1
    UC(ic,1) =   0.500*UF(iif,1) &
               + 0.125*( UF(iif+1,1)+UF(iif-1,1) &
               +         UF(iif,2)+UF(iif,nPsiF-1) )
    UC(ic,nPsiC) =   0.500*UF(iif,nPsiF) &
                   + 0.125*( UF(iif+1,nPsiF)+UF(iif-1,nPsiF) &
                   +         UF(iif,2)+UF(iif,nPsiF-1) )
    iif = iif + 2
  end do

  jf = 1
  do jc = 1, nPsiC
    UC(1,jc) = UF(1,jf)
    UC(nThetaC,jc) = UF(nThetaF,jf)
    jf = jf + 2
  end do

  ! Ensure restricted solution is periodic.

  UC(1:nThetaC,nPsiC) = &
     0.50*(UC(1:nThetaC,1)+UC(1:nThetaC,nPsiC))
  UC(1:nThetaC,1) = UC(1:nThetaC,nPsiC)

end subroutine ionosphere_restrict

subroutine ionosphere_prolong(UF, nThetaF, nPsiF, &
                              UC, nThetaC, nPsiC)
  !\
  ! This subroutine uses bilinear interpolation
  ! to prolong the coarse grid solution to the fine grid.
  !/
  use ModIonosphere
  implicit none

  integer :: nThetaF,nPsiF, &
             nThetaC,nPsiC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) :: UF
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) :: UC

  integer :: ic,iif,jc,jf

  ! Do elements that are direct copies.

  jf = 1
  do jc = 1, nPsiC
     do ic = 1, nThetaC  
        UF(2*ic-1,jf) = UC(ic,jc)
     end do
     jf = jf + 2
  end do
        
  ! Do odd-numbered columns, interpolating vertically.

  do jf = 1, nPsiF, 2
     do iif = 2, nThetaF-1, 2
        UF(iif,jf) = 0.50*(UF(iif+1,jf)+UF(iif-1,jf))
     end do
  end do
        
  ! Do even-numbered columns, interpolating horizontally.

  do jf = 2, nPsiF-1, 2
     do iif = 1, nThetaF
        UF(iif,jf) = 0.50*(UF(iif,jf+1)+UF(iif,jf-1))
     end do
  end do

  ! Ensure prolonged solution is periodic.

  UF(1:nThetaF,nPsiF) = &
     0.50*(UF(1:nThetaF,1)+UF(1:nThetaF,nPsiF))
  UF(1:nThetaF,1) = UF(1:nThetaF,nPsiF)

end subroutine ionosphere_prolong

subroutine ionosphere_residual(RES, PHI, RHO, &
                               SigmaThTh, SigmaThPs, SigmaPsPs, &
                               dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                               dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                               Theta, Psi, Radius, nTheta, nPsi, &
                               iboundary, RESIDUAL)
  !\
  ! This subroutine uses evaluates the residual operator
  ! for the linear elliptic equation governing the
  ! electric field potential PHI
  ! using a standard centred-difference discretization.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi, iboundary
  real :: Radius, RESIDUAL
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  RES, PHI, RHO, &
                  SigmaThTh, SigmaThPs, SigmaPsPs, &
                  dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                  dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                  Theta, Psi

  integer :: i, j, jshift
  real :: dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, sn, cs, sn2, cs2
  real :: kappa_Theta2, kappa_Theta1, &
          kappa_Psi2, kappa_Psi1, kappa_ratio
  real :: RHOmax, PHIpole

  ! Determine mesh spacings.

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  kappa_ratio = 2.00

  ! Interior points.

  RESIDUAL=0.0
  RHOmax=0.00

  do j = 1, nPsi
     do i = 2, nTheta-1 
        if ((iboundary == 0 .and. &
             Theta(i,j) > IONO_NORTH_Theta_Max) .or. &
            (iboundary == 1 .and. &
             Theta(i,j) < IONO_SOUTH_Theta_Min)) then
           RES(i,j) = 0.00
        else
           ! Evaluate effective diffusion coefficients.

           sn = sin(Theta(i,j))
           cs = cos(Theta(i,j)) 
           sn2= sn*sn 
           cs2 = cs*cs 

           kappa_Theta2 =   SigmaThTh(i,j)*sn2 
           kappa_Theta1 =   SigmaThTh(i,j)*sn*cs &
                          + dSigmaThTh_dTheta(i,j)*sn2 &
                          - dSigmaThPs_dPsi(i,j)*sn 
           kappa_Psi2 =   SigmaPsPs(i,j) 
           kappa_Psi1 =   dSigmaThPs_dTheta(i,j)*sn &
                        + dSigmaPsPs_dPsi(i,j) 

!!$           if (abs(kappa_Theta1) >= kappa_ratio*abs(kappa_Theta2)/dd) then
!!$              if (kappa_Theta1 >= 0.00) then
!!$                 kappa_Theta1 = kappa_ratio*abs(kappa_Theta2)/dd 
!!$              else
!!$                 kappa_Theta1 = -kappa_ratio*abs(kappa_Theta2)/dd 
!!$              end if
!!$           end if
!!$           if (abs(kappa_Psi1) >= kappa_ratio*abs(kappa_Psi2)/dd) then
!!$              if (kappa_Psi1 >= 0.00) then
!!$                 kappa_Psi1 = kappa_ratio*abs(kappa_Psi2)/dd 
!!$              else
!!$                 kappa_Psi1 = -kappa_ratio*abs(kappa_Psi2)/dd 
!!$              end if
!!$           end if
       
           ! Evaluate residual.

           if (j > 1 .and. j < nPsi) then
              RES(i,j)=RHO(i,j)- &
                       kappa_Theta2*(PHI(i+1,j)-2.0*PHI(i,j)+PHI(i-1,j))/dd2- &
                       0.50*kappa_Theta1*(PHI(i+1,j)-PHI(i-1,j))/dd- &
                       kappa_Psi2*(PHI(i,j+1)-2.0*PHI(i,j)+PHI(i,j-1))/dd2- &
                       0.50*kappa_Psi1*(PHI(i,j+1)-PHI(i,j-1))/dd 
           else if (j == 1) then
              RES(i,j)=RHO(i,j)- &
                       kappa_Theta2*(PHI(i+1,j)-2.0*PHI(i,j)+PHI(i-1,j))/dd2- &
                       0.50*kappa_Theta1*(PHI(i+1,j)-PHI(i-1,j))/dd- &
                       kappa_Psi2*(PHI(i,j+1)-2.0*PHI(i,j)+PHI(i,nPsi-1))/dd2- &
                       0.50*kappa_Psi1*(PHI(i,j+1)-PHI(i,nPsi-1))/dd
           else
              RES(i,j)=RHO(i,j)- &
                       kappa_Theta2*(PHI(i+1,j)-2.0*PHI(i,j)+PHI(i-1,j))/dd2- &
                       0.50*kappa_Theta1*(PHI(i+1,j)-PHI(i-1,j))/dd- &
                       kappa_Psi2*(PHI(i,2)-2.0*PHI(i,j)+PHI(i,j-1))/dd2- &
                       0.50*kappa_Psi1*(PHI(i,2)-PHI(i,j-1))/dd
           end if
        end if

        RESIDUAL = RESIDUAL + abs(RES(i,j))
        if (abs(RHO(i,j)) > RHOmax) RHOmax = abs(RHO(i,j))

     end do
  end do

  ! Theta-coordinate Boundary points.

  select case (iboundary)
     case (0)
!!$        ! Zero at north pole, zero at equator.
!!$        do j = 1, nPsi
!!$           RES(1,j) = 0.00
!!$           RESIDUAL = RESIDUAL + abs(RES(1,j))
!!$           if (abs(RHO(1,j)) > RHOmax) RHOmax = abs(RHO(1,j))
!!$
!!$           RES(nTheta,j) = 0.00
!!$           RESIDUAL = RESIDUAL + abs(RES(nTheta,j))
!!$           if (abs(RHO(nTheta,j)) > RHOmax) RHOmax = abs(RHO(nTheta,j))
!!$        end do

        ! Average at north pole, zero at equator.
        PHIpole = 0.00
        if (nTheta > 3) then
           do j = 1, nPsi
              PHIpole = PHIpole + &
                        (PHI(2,j)-&
                        0.00*(-3.00*PHI(2,j)+4.00*PHI(3,j)-PHI(4,j)))
           end do
           PHIpole = PHIpole/real(nPsi)
        else
        end if

        do j = 1, nPsi
           RES(1,j) = PHIpole - PHI(1,j)
           RESIDUAL = RESIDUAL + abs(RES(1,j))
           if (abs(RHO(1,j)) > RHOmax) RHOmax = abs(RHO(1,j))

           RES(nTheta,j) = 0.00
           RESIDUAL = RESIDUAL + abs(RES(nTheta,j))
           if (abs(RHO(nTheta,j)) > RHOmax) RHOmax = abs(RHO(nTheta,j))
        end do

     case (1)
!!$        ! Zero at south pole, zero at equator.
!!$        do j = 1, nPsi
!!$           RES(1,j) = 0.00
!!$           RESIDUAL = RESIDUAL + abs(RES(1,j))
!!$           if (abs(RHO(1,j)) > RHOmax) RHOmax = abs(RHO(1,j))
!!$
!!$           RES(nTheta,j)= 0.00
!!$           RESIDUAL = RESIDUAL + abs(RES(nTheta,j))
!!$           if (abs(RHO(nTheta,j)) > RHOmax) RHOmax = abs(RHO(nTheta,j))
!!$        end do

        ! Average at south pole, zero at equator.

        PHIpole = 0.00
        if (nTheta > 3) then
           do j = 1, nPsi
              PHIpole = PHIpole + &
                        (PHI(nTheta-1,j)+ &
                        0.00*(3.00*PHI(nTheta-1,j)-4.00*PHI(nTheta-2,j)+PHI(nTheta-3,j)))
           end do
           PHIpole = PHIpole/real(nPsi)
        else
        end if

        do j = 1, nPsi
           RES(1,j) = 0.00
           RESIDUAL = RESIDUAL + abs(RES(1,j))
           if (abs(RHO(1,j)) > RHOmax) RHOmax = abs(RHO(1,j))

           RES(nTheta,j)= PHIpole - PHI(nTheta,j)
           RESIDUAL = RESIDUAL + abs(RES(nTheta,j))
           if (abs(RHO(nTheta,j)) > RHOmax) RHOmax = abs(RHO(nTheta,j))
        end do

  end select 

  ! Determine the average normalized residual.

  if (RHOmax <= 0.00) RHOmax = 1.00
  RESIDUAL = RESIDUAL/(nTheta*nPsi*RHOmax)

end subroutine ionosphere_residual

subroutine ionosphere_smoother(PHI, RHO, &
                               SigmaThTh, SigmaThPs, SigmaPsPs, &
                               dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                               dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                               Theta, Psi, Radius, nTheta, nPsi, &
                               iboundary)
  !\
  ! This subroutine uses a red-black Gauss-Seidel scheme
  ! to smooth the solution to the linear elliptic equation
  ! governing the electric field potential PHI
  ! as defined on the domain 0 < Theta < PI/2 (this is for
  ! the northern hemisphere, PI/2 < Theta < PI for southern
  ! hemisphere), , 0 < Psi < 2 PI, and subject to the 
  ! Dirichlet boundary data
  !
  !      PHI(0,Psi) = PHI(PI,Psi) = 0,
  !
  !      PHI(PI/2,Psi) = 0,
  !
  !      PHI(Theta,0) = PHI(Theta,2*PI).
  !
  ! The current value of the solution, PHI, is updated using
  ! Gauss-Seidel centred-difference discretization procedure 
  ! and the right-hand side function, RHO.
  !/
  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi, iboundary
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI, RHO, &
                  SigmaThTh, SigmaThPs, SigmaPsPs, &
                  dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                  dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                  Theta, Psi

  integer :: i, j, ipass, jsw
  real :: dTheta, dPsi, dTheta2, dPsi2, &
          dd, dd2, sn, cs, sn2, cs2
  real :: kappa_Theta2, kappa_Theta1, &
          kappa_Psi2, kappa_Psi1, kappa_ratio
  real :: PHIpole

  ! Determine mesh spacings.

  dTheta=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
  dPsi=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
  dTheta2=dTheta*dTheta
  dPsi2=dPsi*dPsi
  dd=dPsi
  dd2=dPsi2

  kappa_ratio = 2.00

  ! Red and black sweeps.
  
  do ipass = 1,2
     if (ipass == 1) then
        jsw = 1
     else
        jsw = 0
     end if

     do i = 2, nTheta-1
        do j = jsw+1, nPsi, 2
           if ((iboundary == 0 .and. &
                Theta(i,j) > IONO_NORTH_Theta_Max) .or. &
               (iboundary == 1 .and. &
                Theta(i,j) < IONO_SOUTH_Theta_Min)) then
              PHI(i,j) = 0.00
           else
              ! Evaluate effective diffusion coefficients.

              sn = sin(Theta(i,j))
              cs = cos(Theta(i,j)) 
              sn2= sn*sn 
              cs2 = cs*cs 
  
              kappa_Theta2 =   SigmaThTh(i,j)*sn2 
              kappa_Theta1 =   SigmaThTh(i,j)*sn*cs &
                             + dSigmaThTh_dTheta(i,j)*sn2 &
                             - dSigmaThPs_dPsi(i,j)*sn 
              kappa_Psi2 =   SigmaPsPs(i,j) 
              kappa_Psi1 =   dSigmaThPs_dTheta(i,j)*sn &
                           + dSigmaPsPs_dPsi(i,j)

!!$              if (abs(kappa_Theta1) >= kappa_ratio*abs(kappa_Theta2)/dd) then
!!$                 if (kappa_Theta1 >= 0.00) then
!!$                    kappa_Theta1 = kappa_ratio*abs(kappa_Theta2)/dd 
!!$                 else
!!$                    kappa_Theta1 = -kappa_ratio*abs(kappa_Theta2)/dd 
!!$                 end if
!!$              end if
!!$              if (abs(kappa_Psi1) >= kappa_ratio*abs(kappa_Psi2)/dd) then
!!$                 if (kappa_Psi1 >= 0.00) then
!!$                    kappa_Psi1 = kappa_ratio*abs(kappa_Psi2)/dd 
!!$                 else
!!$                    kappa_Psi1 = -kappa_ratio*abs(kappa_Psi2)/dd 
!!$                 end if
!!$              end if
  
              ! Gauss-Seidel formula.
          
              if (j > 1 .and. j < nPsi) then
                 PHI(i,j) = (kappa_Theta2*(PHI(i+1,j)+PHI(i-1,j)) &
                           +0.50*dd*kappa_Theta1*(PHI(i+1,j)-PHI(i-1,j)) &
                           +kappa_Psi2*(PHI(i,j+1)+PHI(i,j-1)) &
                           +0.50*dd*kappa_Psi1*(PHI(i,j+1)-PHI(i,j-1)) &
                           -dd2*RHO(i,j))/ &
                           (2.00*(kappa_Theta2+kappa_Psi2))
              else if (j == 1) then
                 PHI(i,j) = (kappa_Theta2*(PHI(i+1,j)+PHI(i-1,j)) &
                           +0.50*dd*kappa_Theta1*(PHI(i+1,j)-PHI(i-1,j)) &
                           +kappa_Psi2*(PHI(i,j+1)+PHI(i,nPsi-1)) &
                           +0.50*dd*kappa_Psi1*(PHI(i,j+1)-PHI(i,nPsi-1)) &
                           -dd2*RHO(i,j))/ &
                           (2.00*(kappa_Theta2+kappa_Psi2))
              else
                 PHI(i,j) = (kappa_Theta2*(PHI(i+1,j)+PHI(i-1,j)) &
                           +0.50*dd*kappa_Theta1*(PHI(i+1,j)-PHI(i-1,j)) &
                           +kappa_Psi2*(PHI(i,2)+PHI(i,j-1)) &
                           +0.50*dd*kappa_Psi1*(PHI(i,2)-PHI(i,j-1)) &
                           -dd2*RHO(i,j))/ &
                           (2.00*(kappa_Theta2+kappa_Psi2))
              end if
           end if

        end do
     end do
  end do
  
  ! Theta-coordinate boundary points.

  select case (iboundary)
     case (0)
!!$        ! Zero at north pole, zero at equator.
!!$        do j = 1, nPsi
!!$           PHI(1,j) = 0.00
!!$           PHI(nTheta,j) = 0.00
!!$        end do

        ! Average at north pole, zero at equator.
        PHIpole = 0.00
        if (nTheta > 3) then
           do j = 1, nPsi
              PHIpole = PHIpole + &
                        (PHI(2,j)-0.00*(-3.00*PHI(2,j)+4.00*PHI(3,j)-PHI(4,j)))
           end do
           PHIpole = PHIpole/real(nPsi)
        else
        end if

        do j = 1, nPsi
           PHI(1,j) = PHI(1,j)+0.50*(PHIpole-PHI(1,j))
           PHI(nTheta,j) = 0.00
        end do

     case (1)
!!$        ! Zero at south pole, zero at equator.
!!$        do j = 1, nPsi
!!$           PHI(1,j) = 0.00
!!$           PHI(nTheta,j) = 0.00
!!$        end do

        ! Average at south pole, zero at equator.
        PHIpole = 0.00
        if (nTheta > 3) then
           do j = 1, nPsi
              PHIpole = PHIpole + &
                        (PHI(nTheta-1,j)+ &
                        0.00*(3.00*PHI(nTheta-1,j)-4.00*PHI(nTheta-2,j)+PHI(nTheta-3,j)))
           end do
           PHIpole = PHIpole/real(nPsi)
        else
        end if

        do j = 1, nPsi
           PHI(1,j) = 0.00
           PHI(nTheta,j) = PHI(nTheta,j)+0.50*(PHIpole-PHI(nTheta,j))
        end do

  end select

end subroutine ionosphere_smoother
