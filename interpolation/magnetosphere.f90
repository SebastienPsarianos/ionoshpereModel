!============================================
!                                           |
!     Magnetosphere/Ionosphere Coupling     |
!                                           |
!============================================

subroutine magnetosphere_fac(Init)
  use ModMain
  use ModPhysics
  use ModParallelAMR
  use ModIonosphere
  implicit none
  include "mpif90.h"

  integer, intent(in) :: Init

  integer :: i, j, k, nn, ns, iBLK, iPE, ierror, &
             itag = 99, status(MPI_STATUS_SIZE)

  select case (Init)
     !\
     ! Determine the mapping points in the magnetosphere for the 
     ! incoming and outgoing field-aligned currents (FAC) to the
     ! ionosphere.  The mapping points are located at a radius of
     ! Rcurrents and storage for this solution information must be
     ! allocated.
     !/
     case (1)
        IONO_PI = 2.00*asin(1.00)

        select case (problem_type)
        case (11)
           !\
           ! Earth magnetosphere problem.
           !/
           IONO_Radius = Rearth
           IONO_Height = 110000.00
           IONO_Ref_Density = 1.0e06*CON_mp*SW_rho_dim
           IONO_Ref_SoundSpeed = SSPearth 

        case (12)
           !\
           ! Saturn magnetosphere problem.
           !/
           IONO_Radius = Rsaturn
           IONO_Height = 1000000.00
           IONO_Ref_Density = 1.0e06*CON_mp*SW_rho_dim
           IONO_Ref_SoundSpeed = SSPsaturn

        case (13)
           !\
           ! Jupiter magnetosphere problem.
           !/
           IONO_Radius = Rjupiter
           IONO_Height = 1000000.00
           IONO_Ref_Density = 1.0e06*CON_mp*SW_rho_dim
           IONO_Ref_SoundSpeed = SSPjupiter

        case default
           !\
           ! Defaults.
           !/
           IONO_Radius = 6378000.00
           IONO_Height = 110000.00
           IONO_Ref_Density = 0.8363000e-20
           IONO_Ref_SoundSpeed = 50000.00

        end select

        IONO_Bdp = Bdp* &
                   sqrt(IONO_MU*IONO_Ref_Density* &
                        IONO_Ref_SoundSpeed*IONO_Ref_SoundSpeed)
        IONO_Radius_Mag_Boundary = Rbody*IONO_Radius
    
        allocate(nMagBndPts_North_PE(numprocs), stat = ierror)
        if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                   " allocation error for nMagBndPts_North_PE", &
                                 & " numprocs = ",numprocs
        nMagBndPts_North_PE = 0

        allocate(nMagBndPts_South_PE(numprocs), stat = ierror)
        if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                   " allocation error for nMagBndPts_South_PE", &
                                 & " numprocs = ",numprocs
        nMagBndPts_South_PE = 0

        nMagBndPts_North = 0
        nMagBndPts_South = 0
        do iBLK = 1,nMultiBlkLevels
           if (.not.unusedBLK(iBLK) .and. Rmin_BLK(iBLK) < Rcurrents) then
              do k = 1, nCellsK
                 do j = 1, nCellsJ
                    do i = 1, nCellsI
  
                       if (R_BLK(i,j,k,iBLK) > Rcurrents .and. &
                           ((R_BLK(i+1,j,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i-1,j,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j+1,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j-1,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j,k+1,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j,k-1,iBLK) <= Rcurrents)) ) then
                          if (z_BLK(i,j,k,iBLK) >= 0.00) then
                             nMagBndPts_North = nMagBndPts_North + 1 
                          else
                             nMagBndPts_South = nMagBndPts_South + 1 
                          end if
                       end if
  
                    end do !end i loop
                 end do !end j loop
              end do !end k loop

           end if
        end do

        if (nMagBndPts_North > 0) then
           allocate(Xmag_North(nMagBndPts_North), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for Xmag_North", &
                                    & " nMagBndPts_North = ",nMagBndPts_North 
           Xmag_North = 0.00
           allocate(Ymag_North(nMagBndPts_North), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for Ymag_North", &
                                    & " nMagBndPts_North = ",nMagBndPts_North 
           Ymag_North = 0.00
           allocate(Zmag_North(nMagBndPts_North), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for Zmag_North", &
                                    & " nMagBndPts_North = ",nMagBndPts_North
           Zmag_North = 0.00
           allocate(JXmag_North(nMagBndPts_North), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for JXmag_North", &
                                    & " nMagBndPts_North = ",nMagBndPts_North
           JXmag_North = 0.00
           allocate(JYmag_North(nMagBndPts_North), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for JYmag_North", &
                                    & " nMagBndPts_North = ",nMagBndPts_North
           JYmag_North = 0.00
           allocate(JZmag_North(nMagBndPts_North), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for JZmag_North", &
                                    & " nMagBndPts_North = ",nMagBndPts_North
           JZmag_North = 0.00
        end if

        if (nMagBndPts_South > 0) then
           allocate(Xmag_South(nMagBndPts_South), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for Xmag_South", &
                                    & " nMagBndPts_South = ",nMagBndPts_South
           Xmag_South = 0.00
           allocate(Ymag_South(nMagBndPts_South), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for Ymag_South", &
                                    & " nMagBndPts_South = ",nMagBndPts_South
           Ymag_South = 0.00
           allocate(Zmag_South(nMagBndPts_South), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for Zmag_South", &
                                    & " nMagBndPts_South = ",nMagBndPts_South
           Zmag_South = 0.00
           allocate(JXmag_South(nMagBndPts_South), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for JXmag_South", &
                                    & " nMagBndPts_South = ",nMagBndPts_South
           JXmag_South = 0.00
           allocate(JYmag_South(nMagBndPts_South), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for JYmag_South", &
                                    & " nMagBndPts_South = ",nMagBndPts_South
           JYmag_South = 0.00
           allocate(JZmag_South(nMagBndPts_South), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for JZmag_South", &
                                    & " nMagBndPts_South = ",nMagBndPts_South
           JZmag_South = 0.00
        end if

        nn = 0
        ns = 0
        do iBLK = 1,nMultiBlkLevels
           if (.not.unusedBLK(iBLK) .and. Rmin_BLK(iBLK) < Rcurrents) then
              do k = 1, nCellsK
                 do j = 1, nCellsJ
                    do i = 1, nCellsI
  
                       if (R_BLK(i,j,k,iBLK) > Rcurrents .and. &
                           ((R_BLK(i+1,j,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i-1,j,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j+1,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j-1,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j,k+1,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j,k-1,iBLK) <= Rcurrents)) ) then
                          if (z_BLK(i,j,k,iBLK) >= 0.00) then
                             nn = nn + 1
                             Xmag_North(nn) = x_BLK(i,j,k,iBLK)
                             Ymag_North(nn) = y_BLK(i,j,k,iBLK)
                             Zmag_North(nn) = z_BLK(i,j,k,iBLK)
                          else
                             ns = ns + 1
                             Xmag_South(ns) = x_BLK(i,j,k,iBLK)
                             Ymag_South(ns) = y_BLK(i,j,k,iBLK)
                             Zmag_South(ns) = z_BLK(i,j,k,iBLK)
                          end if
                       end if
  
                    end do !end i loop
                 end do !end j loop
              end do !end k loop

           end if
        end do

        call MPI_ALLGATHER(nMagBndPts_North, 1, MPI_INTEGER, &
                           nMagBndPts_North_PE(1), 1, MPI_INTEGER, &
                           MPI_COMM_WORLD, ierror)

        call MPI_ALLGATHER(nMagBndPts_South, 1, MPI_INTEGER, &
                           nMagBndPts_South_PE(1), 1, MPI_INTEGER, &
                           MPI_COMM_WORLD, ierror)

        if (me_world == 0) then
           IONO_NORTH_nMagBndPts = 0
           do iPE = 1, numprocs 
              IONO_NORTH_nMagBndPts = IONO_NORTH_nMagBndPts + nMagBndPts_North_PE(iPE)
           end do

           IONO_South_nMagBndPts = 0
           do iPE = 1, numprocs 
              IONO_South_nMagBndPts = IONO_South_nMagBndPts + nMagBndPts_South_PE(iPE)
           end do

           allocate(MAG_NORTH_X(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_X"
           MAG_NORTH_X = 0.00
           allocate(MAG_NORTH_Y(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Y"
           MAG_NORTH_Y = 0.00
           allocate(MAG_NORTH_Z(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Z"
           MAG_NORTH_Z = 0.00
           allocate(MAG_NORTH_Theta(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Theta"
           MAG_NORTH_Theta = 0.00
           allocate(MAG_NORTH_Psi(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Psi"
           MAG_NORTH_Psi = 0.00
           allocate(MAG_NORTH_Jx(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Jx"
           MAG_NORTH_Jx = 0.00
           allocate(MAG_NORTH_Jy(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Jy"
           MAG_NORTH_Jy = 0.00
           allocate(MAG_NORTH_Jz(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_Jz"
           MAG_NORTH_Jz = 0.00
           allocate(MAG_NORTH_JR(IONO_NORTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_NORTH_JR"
           MAG_NORTH_JR = 0.00

           allocate(MAG_SOUTH_X(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_X"
           MAG_SOUTH_X = 0.00
           allocate(MAG_SOUTH_Y(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Y"
           MAG_SOUTH_Y = 0.00
           allocate(MAG_SOUTH_Z(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Z"
           MAG_SOUTH_Z = 0.00
           allocate(MAG_SOUTH_Theta(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Theta"
           MAG_SOUTH_Theta = 0.00
           allocate(MAG_SOUTH_Psi(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Psi"
           MAG_SOUTH_Psi = 0.00
           allocate(MAG_SOUTH_Jx(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Jx"
           MAG_SOUTH_Jx = 0.00
           allocate(MAG_SOUTH_Jy(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Jy"
           MAG_SOUTH_Jy = 0.00
           allocate(MAG_SOUTH_Jz(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_Jz"
           MAG_SOUTH_Jz = 0.00
           allocate(MAG_SOUTH_JR(IONO_SOUTH_nMagBndPts), stat = ierror)
           if (ierror > 0) write(*,*) "magnetosphere_fac: PE = ",me_world, &
                                      " allocation error for MAG_SOUTH_JR"
           MAG_SOUTH_JR = 0.00

           nn = 0
           if (nMagBndPts_North_PE(1) > 0) then
              MAG_NORTH_X(nn+1:nn+nMagBndPts_North_PE(1)) = Xmag_North
              MAG_NORTH_Y(nn+1:nn+nMagBndPts_North_PE(1)) = Ymag_North
              MAG_NORTH_Z(nn+1:nn+nMagBndPts_North_PE(1)) = Zmag_North
              nn = nn + nMagBndPts_North_PE(1)
           end if
           do iPE = 2, numprocs 
              if (nMagBndPts_North_PE(iPE) > 0) then
                 call MPI_RECV(MAG_NORTH_X(nn+1), nMagBndPts_North_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_NORTH_Y(nn+1), nMagBndPts_North_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_NORTH_Z(nn+1), nMagBndPts_North_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 nn = nn + nMagBndPts_North_PE(iPE)
              end if
           end do

           ns = 0
           if (nMagBndPts_South_PE(1) > 0) then
              MAG_SOUTH_X(ns+1:ns+nMagBndPts_South_PE(1)) = Xmag_South
              MAG_SOUTH_Y(ns+1:ns+nMagBndPts_South_PE(1)) = Ymag_South
              MAG_SOUTH_Z(ns+1:ns+nMagBndPts_South_PE(1)) = Zmag_South
              ns = ns + nMagBndPts_South_PE(1)
           end if
           do iPE = 2, numprocs 
              if (nMagBndPts_South_PE(iPE) > 0) then
                 call MPI_RECV(MAG_SOUTH_X(ns+1), nMagBndPts_South_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_SOUTH_Y(ns+1), nMagBndPts_South_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_SOUTH_Z(ns+1), nMagBndPts_South_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 ns = ns + nMagBndPts_South_PE(iPE)
              end if
           end do
        else
           if (nMagBndPts_North > 0) then
              call MPI_SEND(Xmag_North(1), nMagBndPts_North, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(Ymag_North(1), nMagBndPts_North, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(Zmag_North(1), nMagBndPts_North, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
           end if

           if (nMagBndPts_South > 0) then
              call MPI_SEND(Xmag_South(1), nMagBndPts_South, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(Ymag_South(1), nMagBndPts_South, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(Zmag_South(1), nMagBndPts_South, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
           end if
        end if

     !\
     ! Deallocate storage for FAC currents.
     !/
     case (-1)
        deallocate(nMagBndPts_North_PE)
        deallocate(nMagBndPts_South_PE)

        if (nMagBndPts_North > 0) then
           deallocate(Xmag_North)
           deallocate(Ymag_North)
           deallocate(Zmag_North)
           deallocate(JXmag_North)
           deallocate(JYmag_North)
           deallocate(JZmag_North)
        end if

        if (nMagBndPts_South > 0) then
           deallocate(Xmag_South)
           deallocate(Ymag_South)
           deallocate(Zmag_South)
           deallocate(JXmag_South)
           deallocate(JYmag_South)
           deallocate(JZmag_South)
        end if

        if (me_world == 0) then
           deallocate(MAG_NORTH_X)
           deallocate(MAG_NORTH_Y)
           deallocate(MAG_NORTH_Z)
           deallocate(MAG_NORTH_Theta)
           deallocate(MAG_NORTH_Psi)
           deallocate(MAG_NORTH_Jx)
           deallocate(MAG_NORTH_Jy)
           deallocate(MAG_NORTH_Jz)
           deallocate(MAG_NORTH_JR)
          
           deallocate(MAG_SOUTH_X)
           deallocate(MAG_SOUTH_Y)
           deallocate(MAG_SOUTH_Z)
           deallocate(MAG_SOUTH_Theta)
           deallocate(MAG_SOUTH_Psi)
           deallocate(MAG_SOUTH_Jx)
           deallocate(MAG_SOUTH_Jy)
           deallocate(MAG_SOUTH_Jz)
           deallocate(MAG_SOUTH_JR)
        end if

     !\
     ! Calculate the incoming and outgoing field-aligned currents
     ! (FAC) to the ionosphere at the mapping points in the magnetosphere
     ! defined at a radius of Rcurrents and send this solution 
     ! information to the PE that is performing the ionosphere
     ! potential solution calculation.
     !/
     case default
        IONO_Bdp = Bdp* &
                   sqrt(IONO_MU*IONO_Ref_Density* &
                        IONO_Ref_SoundSpeed*IONO_Ref_SoundSpeed)
        IONO_Radius_Mag_Boundary = Rbody*IONO_Radius

        nn = 0
        ns = 0
        do iBLK = 1,nMultiBlkLevels
           if (.not.unusedBLK(iBLK) .and. Rmin_BLK(iBLK) < Rcurrents) then
              do k = 1, nCellsK
                 do j = 1, nCellsJ
                    do i = 1, nCellsI
  
                       if (R_BLK(i,j,k,iBLK) > Rcurrents .and. &
                           ((R_BLK(i+1,j,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i-1,j,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j+1,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j-1,k,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j,k+1,iBLK) <= Rcurrents) .or. &
                            (R_BLK(i,j,k-1,iBLK) <= Rcurrents)) ) then
                          if (z_BLK(i,j,k,iBLK) >= 0.00) then
                             nn = nn + 1
                             JXmag_North(nn) = 0.50* &
                                               ((Bz_BLK(i,j+1,k,iBLK) - &
                                                 Bz_BLK(i,j-1,k,iBLK))/dy_BLK(iBLK) - &
                                                (By_BLK(i,j,k+1,iBLK) - &
                                                 By_BLK(i,j,k-1,iBLK))/dz_BLK(iBLK))
                             JYmag_North(nn) = 0.50* &
                                               ((Bx_BLK(i,j,k+1,iBLK) - &
                                                 Bx_BLK(i,j,k-1,iBLK))/dz_BLK(iBLK) - &
                                                (Bz_BLK(i+1,j,k,iBLK) - &
                                                 Bz_BLK(i-1,j,k,iBLK))/dx_BLK(iBLK))
                             JZmag_North(nn) = 0.50* &
                                               ((By_BLK(i+1,j,k,iBLK) - &
                                                 By_BLK(i-1,j,k,iBLK))/dx_BLK(iBLK) - &
                                                (Bx_BLK(i,j+1,k,iBLK) - &
                                                 Bx_BLK(i,j-1,k,iBLK))/dy_BLK(iBLK))
                          else
                             ns = ns + 1
                             JXmag_South(ns) = 0.50* &
                                               ((Bz_BLK(i,j+1,k,iBLK) - &
                                                 Bz_BLK(i,j-1,k,iBLK))/dy_BLK(iBLK) - &
                                                (By_BLK(i,j,k+1,iBLK) - &
                                                 By_BLK(i,j,k-1,iBLK))/dz_BLK(iBLK))
                             JYmag_South(ns) = 0.50* &
                                               ((Bx_BLK(i,j,k+1,iBLK) - &
                                                 Bx_BLK(i,j,k-1,iBLK))/dz_BLK(iBLK) - &
                                                (Bz_BLK(i+1,j,k,iBLK) - &
                                                 Bz_BLK(i-1,j,k,iBLK))/dx_BLK(iBLK))
                             JZmag_South(ns) = 0.50* &
                                               ((By_BLK(i+1,j,k,iBLK) - &
                                                 By_BLK(i-1,j,k,iBLK))/dx_BLK(iBLK) - &
                                                (Bx_BLK(i,j+1,k,iBLK) - &
                                                 Bx_BLK(i,j-1,k,iBLK))/dy_BLK(iBLK))
                          end if
                       end if
  
                    end do !end i loop
                 end do !end j loop
              end do !end k loop

           end if
        end do

        if (me_world == 0) then
           nn = 0
           if (nMagBndPts_North_PE(1) > 0) then
              MAG_NORTH_Jx(nn+1:nn+nMagBndPts_North_PE(1)) = JXmag_North
              MAG_NORTH_Jy(nn+1:nn+nMagBndPts_North_PE(1)) = JYmag_North
              MAG_NORTH_Jz(nn+1:nn+nMagBndPts_North_PE(1)) = JZmag_North
              nn = nn + nMagBndPts_North_PE(1)
           end if
           do iPE = 2, numprocs 
              if (nMagBndPts_North_PE(iPE) > 0) then
                 call MPI_RECV(MAG_NORTH_Jx(nn+1), nMagBndPts_North_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_NORTH_Jy(nn+1), nMagBndPts_North_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_NORTH_Jz(nn+1), nMagBndPts_North_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 nn = nn + nMagBndPts_North_PE(iPE)
              end if
           end do

           ns = 0
           if (nMagBndPts_South_PE(1) > 0) then
              MAG_SOUTH_Jx(ns+1:ns+nMagBndPts_South_PE(1)) = JXmag_South
              MAG_SOUTH_Jy(ns+1:ns+nMagBndPts_South_PE(1)) = JYmag_South
              MAG_SOUTH_Jz(ns+1:ns+nMagBndPts_South_PE(1)) = JZmag_South
              ns = ns + nMagBndPts_South_PE(1)
           end if
           do iPE = 2, numprocs 
              if (nMagBndPts_South_PE(iPE) > 0) then
                 call MPI_RECV(MAG_SOUTH_Jx(ns+1), nMagBndPts_South_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_SOUTH_Jy(ns+1), nMagBndPts_South_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 call MPI_RECV(MAG_SOUTH_Jz(ns+1), nMagBndPts_South_PE(iPE), &
                               MPI_REAL, iPE-1, itag, MPI_COMM_WORLD, status, ierror)
                 ns = ns + nMagBndPts_South_PE(iPE)
              end if
           end do
        else
           if (nMagBndPts_North > 0) then
              call MPI_SEND(JXmag_North(1), nMagBndPts_North, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(JYmag_North(1), nMagBndPts_North, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(JZmag_North(1), nMagBndPts_North, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
           end if

           if (nMagBndPts_South > 0) then
              call MPI_SEND(JXmag_South(1), nMagBndPts_South, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(JYmag_South(1), nMagBndPts_South, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
              call MPI_SEND(JZmag_South(1), nMagBndPts_South, &
                            MPI_REAL, 0, itag, MPI_COMM_WORLD, ierror)
           end if
        end if

  end select

end subroutine magnetosphere_fac

subroutine magnetosphere_BCs(Init)
  use ModMain
  use ModPhysics
  use ModParallelAMR
  use ModIonosphere
  implicit none
  include "mpif90.h"

  integer, intent(in) :: Init

  integer :: isize, ierror

  !\
  ! Send the predicted ionospheric convection velocities to
  ! all PEs for use in the application of the ionospheric boundary
  ! conditions of magnetosphere calculations.
  !/
  isize = IONO_nTheta*IONO_nPsi

  if (Init.ne.0) then
     call MPI_Bcast(IONO_NORTH_Theta(1,1), &
                    isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(IONO_NORTH_Psi(1,1), &
                    isize,MPI_REAL,0,MPI_COMM_WORLD,ierror) 
     call MPI_Bcast(IONO_SOUTH_Theta(1,1), &
                    isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(IONO_SOUTH_Psi(1,1), &
                    isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  end if

  call MPI_Bcast(IONO_NORTH_UR_BC(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(IONO_NORTH_UTh_BC(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(IONO_NORTH_UPs_BC(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(IONO_SOUTH_UR_BC(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(IONO_SOUTH_UTh_BC(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(IONO_SOUTH_UPs_BC(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)

end subroutine magnetosphere_BCs

subroutine magnetosphere_velB0
  use ModMain
  use ModPhysics
  use ModParallelAMR
  use ModIonosphere
  implicit none
  include "mpif90.h"

  integer :: i, j, k, i0, j0, k0, iBLK, isize, ierror

  real :: Radius, dRadial
  real :: VrCell, VthetaCell, VphiCell, &
          ThetaCell, PhiCell, &
          cosThetaCell, sinThetaCell, cosPhiCell, sinPhiCell, &
          dTheta_N, dPsi_N, dd_N, dTheta_S, dPsi_S, dd_S
  real :: xC1, yC1, zC1, xC2, yC2, zC2, &
          xC3, yC3, zC3, xC4, yC4, zC4, &
          xC5, yC5, zC5, xC6, yC6, zC6, &
          xC7, yC7, zC7, xC8, yC8, zC8, &
          fC1, fC2, fC3, fC4, fC5, fC6, fC7, fC8, &
          a_trilin, b_trilin, c_trilin, d_trilin, &
          e_trilin, f_trilin, g_trilin, h_trilin

  !\
  ! Send the predicted ionospheric potential solution to
  ! all PEs for use in the calculation of the B0 convection
  ! velocity field.
  !/
  isize = IONO_nTheta*IONO_nPsi
  call MPI_Bcast(IONO_NORTH_PHI(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(IONO_SOUTH_PHI(1,1), &
                 isize,MPI_REAL,0,MPI_COMM_WORLD,ierror)

  !\
  ! Using this ionospheric potential solution, compute
  ! the B0 convection velocity field on a spherical grid.
  !/
  Radius = IONO_Radius + IONO_Height
  call magnetosphere_calc_velB0(IONO_NORTH_PHI_BC, &
                                IONO_NORTH_ETh_BC, IONO_NORTH_EPs_BC, &
                                IONO_NORTH_UR_velB0, &
                                IONO_NORTH_UTh_velB0, IONO_NORTH_UPs_velB0, &
                                IONO_NORTH_PHI, &
                                IONO_NORTH_Theta, IONO_NORTH_Psi, &
                                Radius, IONO_nTheta, IONO_nPsi, IONO_nRADIAL)
  call magnetosphere_calc_velB0(IONO_SOUTH_PHI_BC, &
                                IONO_SOUTH_ETh_BC, IONO_SOUTH_EPs_BC, &
                                IONO_SOUTH_UR_velB0, &
                                IONO_SOUTH_UTh_velB0, IONO_SOUTH_UPs_velB0, &
                                IONO_SOUTH_PHI, &
                                IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                                Radius, IONO_nTheta, IONO_nPsi, IONO_nRADIAL)

  !\
  ! Interpolate the B0 convection velocity field from the
  ! spherical grid to the unstructured block-based Cartesian
  ! mesh.
  !/
  dTheta_N = (IONO_NORTH_Theta(IONO_nTheta,1)- &
              IONO_NORTH_Theta(1,1))/real(IONO_nTheta-1)
  dPsi_N = (IONO_NORTH_Psi(1,IONO_nPsi)- &
            IONO_NORTH_Psi(1,1))/real(IONO_nPsi-1)
  dd_N = dPsi_N

  dTheta_S = (IONO_SOUTH_Theta(IONO_nTheta,1)- &
              IONO_SOUTH_Theta(1,1))/real(IONO_nTheta-1)
  dPsi_S = (IONO_SOUTH_Psi(1,IONO_nPsi)- &
            IONO_SOUTH_Psi(1,1))/real(IONO_nPsi-1)
  dd_S = dPsi_S

  Radius = (IONO_Radius + IONO_Height)/IONO_Radius
  dRadial = (20.00*Radius-Radius)/real(IONO_nRadial-1)

  do iBLK = 1,nMultiBlkLevels
     if (.not.unusedBLK(iBLK)) then
        do k = 0, nCellsK+1
           do j = 0, nCellsJ+1
              do i = 0, nCellsI+1
  
                 if (R_BLK(i,j,k,iBLK) > Radius .and. &
                     R_BLK(i,j,k,iBLK) < 20.00*Radius) then

                    if (z_BLK(i,j,k,iBLK) > 0.00) then
                       !\
                       ! Northern Hemisphere.
                       !/
                       if (x_BLK(i,j,k,iBLK) == 0.00 .and. &
                           y_BLK(i,j,k,iBLK) == 0.00) then
                          PhiCell = 0.00
                       else
                          PhiCell = atan2(y_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK), &
                                          x_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK))
                       end if
                       if (PhiCell < 0.00) PhiCell = PhiCell + 2.00*IONO_PI
                       ThetaCell = atan2(sqrt(x_BLK(i,j,k,iBLK)*x_BLK(i,j,k,iBLK)+ &
                                              y_BLK(i,j,k,iBLK)*y_BLK(i,j,k,iBLK))/ &
                                         R_BLK(i,j,k,iBLK), &
                                         z_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK))

                       i0 = ThetaCell/dd_N + 1
                       j0 = PhiCell/dd_N + 1
                       k0 = R_BLK(i,j,k,iBLK)/dRadial + 1

                       xC1 = IONO_NORTH_Theta(i0,j0)
                       yC1 = IONO_NORTH_Psi(i0,j0)
                       zC1 = Radius+real(k0-1)*dRadial
                       xC2 = IONO_NORTH_Theta(i0,j0+1)
                       yC2 = IONO_NORTH_Psi(i0,j0+1)
                       zC2 = Radius+real(k0-1)*dRadial
                       xC3 = IONO_NORTH_Theta(i0+1,j0+1)
                       yC3 = IONO_NORTH_Psi(i0+1,j0+1)
                       zC3 = Radius+real(k0-1)*dRadial
                       xC4 = IONO_NORTH_Theta(i0+1,j0)
                       yC4 = IONO_NORTH_Psi(i0+1,j0)
                       zC4 = Radius+real(k0-1)*dRadial
                       xC5 = IONO_NORTH_Theta(i0,j0)
                       yC5 = IONO_NORTH_Psi(i0,j0)
                       zC5 = Radius+real(k0)*dRadial
                       xC6 = IONO_NORTH_Theta(i0,j0+1)
                       yC6 = IONO_NORTH_Psi(i0,j0+1)
                       zC6 = Radius+real(k0)*dRadial
                       xC7 = IONO_NORTH_Theta(i0+1,j0+1)
                       yC7 = IONO_NORTH_Psi(i0+1,j0+1)
                       zC7 = Radius+real(k0)*dRadial
                       xC8 = IONO_NORTH_Theta(i0+1,j0)
                       yC8 = IONO_NORTH_Psi(i0+1,j0)
                       zC8 = Radius+real(k0)*dRadial
       
       
                       fC1 = IONO_NORTH_UR_velB0(i0,j0,k0)
                       fC2 = IONO_NORTH_UR_velB0(i0,j0+1,k0)
                       fC3 = IONO_NORTH_UR_velB0(i0+1,j0+1,k0)
                       fC4 = IONO_NORTH_UR_velB0(i0+1,j0,k0)
                       fC5 = IONO_NORTH_UR_velB0(i0,j0,k0+1)
                       fC6 = IONO_NORTH_UR_velB0(i0,j0+1,k0+1)
                       fC7 = IONO_NORTH_UR_velB0(i0+1,j0+1,k0+1)
                       fC8 = IONO_NORTH_UR_velB0(i0+1,j0,k0+1)
                       a_trilin = fC1
                       b_trilin = (fC4-fC1)/dd_N
                       c_trilin = (fC2-fC1)/dd_N
                       d_trilin = (fC3+fC1-fC4-fC2)/(dd_N*dd_N)
                       e_trilin = (fC5-fC1)/dRadial
                       f_trilin = (fC8+fC1-fC5-fC4)/(dd_N*dRadial)
                       g_trilin = (fC6+fC1-fC2-fC5)/(dd_N*dRadial)
                       h_trilin = (fC7+fC5+fC4+fC2-fC8-fC6-fC3-fC1)/ &
                                  (dd_N*dd_N*dRadial)
                       VrCell = a_trilin + &
                                b_trilin*(ThetaCell-xC1) + &
                                c_trilin*(PhiCell-yC1) + &
                                d_trilin*(ThetaCell-xC1)*(PhiCell-yC1) + &
                                e_trilin*(R_BLK(i,j,k,iBLK)-zC1) + &
                                f_trilin*(ThetaCell-xC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                f_trilin*(PhiCell-yC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                h_trilin*(ThetaCell-xC1)*(PhiCell-yC1)* &
                                (R_BLK(i,j,k,iBLK)-zC1)

                       fC1 = IONO_NORTH_UTh_velB0(i0,j0,k0)
                       fC2 = IONO_NORTH_UTh_velB0(i0,j0+1,k0)
                       fC3 = IONO_NORTH_UTh_velB0(i0+1,j0+1,k0)
                       fC4 = IONO_NORTH_UTh_velB0(i0+1,j0,k0)
                       fC5 = IONO_NORTH_UTh_velB0(i0,j0,k0+1)
                       fC6 = IONO_NORTH_UTh_velB0(i0,j0+1,k0+1)
                       fC7 = IONO_NORTH_UTh_velB0(i0+1,j0+1,k0+1)
                       fC8 = IONO_NORTH_UTh_velB0(i0+1,j0,k0+1)
                       a_trilin = fC1
                       b_trilin = (fC4-fC1)/dd_N
                       c_trilin = (fC2-fC1)/dd_N
                       d_trilin = (fC3+fC1-fC4-fC2)/(dd_N*dd_N)
                       e_trilin = (fC5-fC1)/dRadial
                       f_trilin = (fC8+fC1-fC5-fC4)/(dd_N*dRadial)
                       g_trilin = (fC6+fC1-fC2-fC5)/(dd_N*dRadial)
                       h_trilin = (fC7+fC5+fC4+fC2-fC8-fC6-fC3-fC1)/ &
                                  (dd_N*dd_N*dRadial)
                       VthetaCell = a_trilin + &
                                    b_trilin*(ThetaCell-xC1) + &
                                    c_trilin*(PhiCell-yC1) + &
                                    d_trilin*(ThetaCell-xC1)*(PhiCell-yC1) + &
                                    e_trilin*(R_BLK(i,j,k,iBLK)-zC1) + &
                                    f_trilin*(ThetaCell-xC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                    f_trilin*(PhiCell-yC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                    h_trilin*(ThetaCell-xC1)*(PhiCell-yC1)* &
                                    (R_BLK(i,j,k,iBLK)-zC1)
       
                       fC1 = IONO_NORTH_UPs_velB0(i0,j0,k0)
                       fC2 = IONO_NORTH_UPs_velB0(i0,j0+1,k0)
                       fC3 = IONO_NORTH_UPs_velB0(i0+1,j0+1,k0)
                       fC4 = IONO_NORTH_UPs_velB0(i0+1,j0,k0)
                       fC5 = IONO_NORTH_UPs_velB0(i0,j0,k0+1)
                       fC6 = IONO_NORTH_UPs_velB0(i0,j0+1,k0+1)
                       fC7 = IONO_NORTH_UPs_velB0(i0+1,j0+1,k0+1)
                       fC8 = IONO_NORTH_UPs_velB0(i0+1,j0,k0+1)
                       a_trilin = fC1
                       b_trilin = (fC4-fC1)/dd_N
                       c_trilin = (fC2-fC1)/dd_N
                       d_trilin = (fC3+fC1-fC4-fC2)/(dd_N*dd_N)
                       e_trilin = (fC5-fC1)/dRadial
                       f_trilin = (fC8+fC1-fC5-fC4)/(dd_N*dRadial)
                       g_trilin = (fC6+fC1-fC2-fC5)/(dd_N*dRadial)
                       h_trilin = (fC7+fC5+fC4+fC2-fC8-fC6-fC3-fC1)/ &
                                  (dd_N*dd_N*dRadial)
                       VphiCell = a_trilin + &
                                  b_trilin*(ThetaCell-xC1) + &
                                  c_trilin*(PhiCell-yC1) + &
                                  d_trilin*(ThetaCell-xC1)*(PhiCell-yC1) + &
                                  e_trilin*(R_BLK(i,j,k,iBLK)-zC1) + &
                                  f_trilin*(ThetaCell-xC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                  f_trilin*(PhiCell-yC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                  h_trilin*(ThetaCell-xC1)*(PhiCell-yC1)* &
                                  (R_BLK(i,j,k,iBLK)-zC1)
                    else
                       !\
                       ! Southern Hemisphere.
                       !/
                       if (x_BLK(i,j,k,iBLK) == 0.00 .and. &
                           y_BLK(i,j,k,iBLK) == 0.00) then
                          PhiCell = 0.00
                       else
                          PhiCell = atan2(y_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK), &
                                          x_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK))
                       end if           
                       if (PhiCell < 0.00) PhiCell = PhiCell + 2.00*IONO_PI
                       ThetaCell = IONO_PI - &
                                   atan2(sqrt(x_BLK(i,j,k,iBLK)*x_BLK(i,j,k,iBLK)+ &
                                              y_BLK(i,j,k,iBLK)*y_BLK(i,j,k,iBLK))/ &
                                         R_BLK(i,j,k,iBLK), &
                                         -z_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK))

                       i0 = (ThetaCell-0.50*IONO_PI)/dd_S + 1
                       j0 = PhiCell/dd_S + 1
                       k0 = R_BLK(i,j,k,iBLK)/dRadial + 1

                       xC1 = IONO_SOUTH_Theta(i0,j0)
                       yC1 = IONO_SOUTH_Psi(i0,j0)
                       zC1 = Radius+real(k0-1)*dRadial
                       xC2 = IONO_SOUTH_Theta(i0,j0+1)
                       yC2 = IONO_SOUTH_Psi(i0,j0+1)
                       zC2 = Radius+real(k0-1)*dRadial
                       xC3 = IONO_SOUTH_Theta(i0+1,j0+1)
                       yC3 = IONO_SOUTH_Psi(i0+1,j0+1)
                       zC3 = Radius+real(k0-1)*dRadial
                       xC4 = IONO_SOUTH_Theta(i0+1,j0)
                       yC4 = IONO_SOUTH_Psi(i0+1,j0)
                       zC4 = Radius+real(k0-1)*dRadial
                       xC5 = IONO_SOUTH_Theta(i0,j0)
                       yC5 = IONO_SOUTH_Psi(i0,j0)
                       zC5 = Radius+real(k0)*dRadial
                       xC6 = IONO_SOUTH_Theta(i0,j0+1)
                       yC6 = IONO_SOUTH_Psi(i0,j0+1)
                       zC6 = Radius+real(k0)*dRadial
                       xC7 = IONO_SOUTH_Theta(i0+1,j0+1)
                       yC7 = IONO_SOUTH_Psi(i0+1,j0+1)
                       zC7 = Radius+real(k0)*dRadial
                       xC8 = IONO_SOUTH_Theta(i0+1,j0)
                       yC8 = IONO_SOUTH_Psi(i0+1,j0)
                       zC8 = Radius+real(k0)*dRadial

                       fC1 = IONO_SOUTH_UR_velB0(i0,j0,k0)
                       fC2 = IONO_SOUTH_UR_velB0(i0,j0+1,k0)
                       fC3 = IONO_SOUTH_UR_velB0(i0+1,j0+1,k0)
                       fC4 = IONO_SOUTH_UR_velB0(i0+1,j0,k0)
                       fC5 = IONO_SOUTH_UR_velB0(i0,j0,k0+1)
                       fC6 = IONO_SOUTH_UR_velB0(i0,j0+1,k0+1)
                       fC7 = IONO_SOUTH_UR_velB0(i0+1,j0+1,k0+1)
                       fC8 = IONO_SOUTH_UR_velB0(i0+1,j0,k0+1)
                       a_trilin = fC1
                       b_trilin = (fC4-fC1)/dd_S
                       c_trilin = (fC2-fC1)/dd_S
                       d_trilin = (fC3+fC1-fC4-fC2)/(dd_S*dd_S)
                       e_trilin = (fC5-fC1)/dRadial
                       f_trilin = (fC8+fC1-fC5-fC4)/(dd_S*dRadial)
                       g_trilin = (fC6+fC1-fC2-fC5)/(dd_S*dRadial)
                       h_trilin = (fC7+fC5+fC4+fC2-fC8-fC6-fC3-fC1)/ &
                                  (dd_S*dd_S*dRadial)
                       VrCell = a_trilin + &
                                b_trilin*(ThetaCell-xC1) + &
                                c_trilin*(PhiCell-yC1) + &
                                d_trilin*(ThetaCell-xC1)*(PhiCell-yC1) + &
                                e_trilin*(R_BLK(i,j,k,iBLK)-zC1) + &
                                f_trilin*(ThetaCell-xC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                f_trilin*(PhiCell-yC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                h_trilin*(ThetaCell-xC1)*(PhiCell-yC1)* &
                                (R_BLK(i,j,k,iBLK)-zC1)

                       fC1 = IONO_SOUTH_UTh_velB0(i0,j0,k0)
                       fC2 = IONO_SOUTH_UTh_velB0(i0,j0+1,k0)
                       fC3 = IONO_SOUTH_UTh_velB0(i0+1,j0+1,k0)
                       fC4 = IONO_SOUTH_UTh_velB0(i0+1,j0,k0)
                       fC5 = IONO_SOUTH_UTh_velB0(i0,j0,k0+1)
                       fC6 = IONO_SOUTH_UTh_velB0(i0,j0+1,k0+1)
                       fC7 = IONO_SOUTH_UTh_velB0(i0+1,j0+1,k0+1)
                       fC8 = IONO_SOUTH_UTh_velB0(i0+1,j0,k0+1)
                       a_trilin = fC1
                       b_trilin = (fC4-fC1)/dd_S
                       c_trilin = (fC2-fC1)/dd_S
                       d_trilin = (fC3+fC1-fC4-fC2)/(dd_S*dd_S)
                       e_trilin = (fC5-fC1)/dRadial
                       f_trilin = (fC8+fC1-fC5-fC4)/(dd_S*dRadial)
                       g_trilin = (fC6+fC1-fC2-fC5)/(dd_S*dRadial)
                       h_trilin = (fC7+fC5+fC4+fC2-fC8-fC6-fC3-fC1)/ &
                                  (dd_S*dd_S*dRadial)
                       VthetaCell = a_trilin + &
                                    b_trilin*(ThetaCell-xC1) + &
                                    c_trilin*(PhiCell-yC1) + &
                                    d_trilin*(ThetaCell-xC1)*(PhiCell-yC1) + &
                                    e_trilin*(R_BLK(i,j,k,iBLK)-zC1) + &
                                    f_trilin*(ThetaCell-xC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                    f_trilin*(PhiCell-yC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                    h_trilin*(ThetaCell-xC1)*(PhiCell-yC1)* &
                                    (R_BLK(i,j,k,iBLK)-zC1)

                       fC1 = IONO_SOUTH_UPs_velB0(i0,j0,k0)
                       fC2 = IONO_SOUTH_UPs_velB0(i0,j0+1,k0)
                       fC3 = IONO_SOUTH_UPs_velB0(i0+1,j0+1,k0)
                       fC4 = IONO_SOUTH_UPs_velB0(i0+1,j0,k0)
                       fC5 = IONO_SOUTH_UPs_velB0(i0,j0,k0+1)
                       fC6 = IONO_SOUTH_UPs_velB0(i0,j0+1,k0+1)
                       fC7 = IONO_SOUTH_UPs_velB0(i0+1,j0+1,k0+1)
                       fC8 = IONO_SOUTH_UPs_velB0(i0+1,j0,k0+1)
                       a_trilin = fC1
                       b_trilin = (fC4-fC1)/dd_S
                       c_trilin = (fC2-fC1)/dd_S
                       d_trilin = (fC3+fC1-fC4-fC2)/(dd_S*dd_S)
                       e_trilin = (fC5-fC1)/dRadial
                       f_trilin = (fC8+fC1-fC5-fC4)/(dd_S*dRadial)
                       g_trilin = (fC6+fC1-fC2-fC5)/(dd_S*dRadial)
                       h_trilin = (fC7+fC5+fC4+fC2-fC8-fC6-fC3-fC1)/ &
                                  (dd_S*dd_S*dRadial)
                       VphiCell = a_trilin + &
                                  b_trilin*(ThetaCell-xC1) + &
                                  c_trilin*(PhiCell-yC1) + &
                                  d_trilin*(ThetaCell-xC1)*(PhiCell-yC1) + &
                                  e_trilin*(R_BLK(i,j,k,iBLK)-zC1) + &
                                  f_trilin*(ThetaCell-xC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                  f_trilin*(PhiCell-yC1)*(R_BLK(i,j,k,iBLK)-zC1) + &
                                  h_trilin*(ThetaCell-xC1)*(PhiCell-yC1)* &
                                  (R_BLK(i,j,k,iBLK)-zC1)
                    end if
                 else
                    !\
                    ! Outside of spherical grid.
                    !/
                    VrCell = 0.00
                    VthetaCell = 0.00
                    VphiCell = 0.00
                 end if

                 cosThetaCell = z_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK)
                 sinThetaCell = sqrt(x_BLK(i,j,k,iBLK)**2+y_BLK(i,j,k,iBLK)**2)/ &
                                R_BLK(i,j,k,iBLK)
                 cosPhiCell = x_BLK(i,j,k,iBLK)/ &
                              sqrt(x_BLK(i,j,k,iBLK)**2+y_BLK(i,j,k,iBLK)**2+rTOL**2)
                 sinPhiCell = y_BLK(i,j,k,iBLK)/ &
                              sqrt(x_BLK(i,j,k,iBLK)**2+y_BLK(i,j,k,iBLK)**2+rTOL**2)

                 U0x_BLK(i,j,k,iBLK)  = VrCell*x_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK) + &
                                        VthetaCell*cosThetaCell*cosPhiCell - &
                                        VphiCell*sinPhiCell
                 U0y_BLK(i,j,k,iBLK)  = VrCell*y_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK) + &
                                        VthetaCell*cosThetaCell*sinPhiCell + &
                                        VphiCell*cosPhiCell
                 U0z_BLK(i,j,k,iBLK)  = VrCell*z_BLK(i,j,k,iBLK)/R_BLK(i,j,k,iBLK) - &
                                        VthetaCell*sinThetaCell
  
              end do !end i loop
           end do !end j loop
        end do !end k loop
  
     end if
  end do

end subroutine magnetosphere_velB0

subroutine magnetosphere_calc_velB0(PHI_BC, ETh_BC, EPs_BC, &
                                    UR_velB0,UTh_velB0, UPs_velB0, &
                                    PHI, Theta, Psi, Radius, &
                                    nTheta, nPsi, nRadial)
  !\
  ! For the calculated ionospheric potential solution,
  ! map the potential solution out to the magnetospheric
  ! inner boundary at various radii and determine the convection
  ! velocities of the intrinsic magnetic field to be used
  ! in the magnetosphere calculation.
  !/
  use ModMain
  use ModPhysics
  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi, nRadial
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI_BC, ETh_BC, EPs_BC, &
                  PHI, Theta, Psi
  real, dimension(1:IONO_nTheta,1:IONO_nPsi,1:IONO_nRadial) ::  &
                  UR_velB0,UTh_velB0, UPs_velB0

  logical :: north
  integer :: i, j, k, i0
  real :: Radius_velB0, dRadial, dRadial2, &
          Theta0, dTheta, dPsi, dTheta2, dPsi2, &
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

  dRadial = (20.00*Radius-Radius)/real(nRadial-1)
  dRadial2 = dRadial*dRadial

  do k= 1, nRadial
     ! Determine the radius.
     Radius_velB0 = Radius+real(k-1)*dRadial

     ! Determine the potential at this radial location.
     if (Radius_velB0/IONO_Radius < 1.05) then
        PHI_BC = PHI
     else if (north) then
        do j = 1, nPsi
           PHI_BC(1,j) = PHI(1,j)
           do i = 2, nTheta
              Theta0 = asin(sqrt((Radius/Radius_velB0)* &
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
              Theta0 = IONO_PI - asin(sqrt((Radius/Radius_velB0)* &
                                           (sin(Theta(i,j))**2)))
	      i0 = (Theta0-0.50*IONO_PI)/dd + 1
	      PHI_BC(i,j) = PHI(i0,j) + &
                               (Theta0-Theta(i0,j))* &
                               (PHI(i0+1,j)-PHI(i0,j))/ &
                               dd
           end do
        end do
     end if

     ! Compute the electric field at this radial location.
     do j = 1, nPsi
        if (j > 1 .and. j < nPsi ) then 
           do i = 2, nTheta-1
              sinTheta = sin(Theta(i,j))
              ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                            (2.00*dd*Radius_velB0)
              EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,j-1))/ &
                            (2.00*dd*Radius_velB0*sinTheta)
           end do
           ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                         (dd*Radius_velB0)
           EPs_BC(1,j) = EPs_BC(2,j)
           ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                              (dd*Radius_velB0)
           EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
        else if (j == 1) then
           do i = 2, nTheta-1
              sinTheta = sin(Theta(i,j))
              ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                            (2.00*dd*Radius_velB0)
              EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,nPsi-1))/ &
                            (2.00*dd*Radius_velB0*sinTheta)
           end do
           ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                         (dd*Radius_velB0)
           EPs_BC(1,j) = EPs_BC(2,j)
           ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                              (dd*Radius_velB0)
           EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
        else
           do i = 2, nTheta-1
              sinTheta = sin(Theta(i,j))
              ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                            (2.00*dd*Radius_velB0)
              EPs_BC(i,j) = -(PHI_BC(i,2)-PHI_BC(i,j-1))/ &
                            (2.00*dd*Radius_velB0*sinTheta)
           end do
           ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                         (dd*Radius_velB0)
           EPs_BC(1,j) = EPs_BC(2,j)
           ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                              (dd*Radius_velB0)
           EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
        end if
     end do

     ! Compute the convection velocities at this radial location.
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

           BB = IONO_Bdp*((IONO_Radius/Radius_velB0)**3)
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

!!$           VR = Vp_x*xx + Vp_y*yy + Vp_z*zz
!!$           if (north .and. i == nTheta) then
!!$              Vll = 0.00
!!$           else if (.not.north .and. i == 1) then
!!$              Vll = 0.00
!!$           else
!!$              Vll = -VR/(bx*xx+by*yy+bz*zz)
!!$           end if

           Vll = 0.00
        
           Ux = Vp_x + Vll*bx
           Uy = Vp_y + Vll*by
           Uz = Vp_z + Vll*bz

           UR_velB0(i,j,k) = Ux*xx + Uy*yy + Uz*zz
           UTh_velB0(i,j,k) = ((Ux*xx + Uy*yy)*zz - &
                              Uz*(xx**2+yy**2)) / &
                              sqrt(xx**2+yy**2+IONO_TOLER**2)
           UPs_velB0(i,j,k) = (Uy*xx - Ux*yy)*sinTheta / &
                              (xx**2+yy**2+IONO_TOLER**2)

           UR_velB0(i,j,k) = UR_velB0(i,j,k) / IONO_Ref_SoundSpeed
           UTh_velB0(i,j,k) = UTh_velB0(i,j,k) / IONO_Ref_SoundSpeed
           UPs_velB0(i,j,k) = UPs_velB0(i,j,k) / IONO_Ref_SoundSpeed
        end do  
     end do

  end do

end subroutine magnetosphere_calc_velB0
