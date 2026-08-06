!=====================================
! IONO: Stand-Alone Ionosphere Model |
!=====================================
program IONO
  use ModIonosphere
  implicit none

  write(*,'(/,1X,''IONO: Ionosphere Model.'',/)')

  !\
  ! Peform ionosphere calculations.
  !/
  call ionosphere(0,5)

  !\
  ! Save ionosphere solution in output data file.
  !/
  call ionosphere(0,3)

  write(*,'(/,1X,''IONO: Execution complete.'',/)')

end program IONO
