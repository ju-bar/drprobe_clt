! ---------------------------------------------------------------------
!
! File: 'ConsoleProgressBar.f90'
! Author: J. Barthel, juribarthel@gmx.de
! Date: 06.04.2011
! Last modified: 21.09.2015 (JB)
!
! Subroutines for display of a simple progress-bar in console apps
! - use CSPROG
! - call CSPROG_INIT to initialize a new progress bar for a certain
!   number of function steps and a certain display length
! - increase the value of FSTEP_CUR inside your function to track the
!   progress from state 0 to the max. number of function steps as
!   defined by the call of CSPROG_INIT
! - call CSPROG_UPDATE to update the display
! - call CSPROG_FINISH to finish the display. Do this also if you
!   exit the progress-controlled function before the total ammount
!   of steps is reached. 
!
! ---------------------------------------------------------------------

module CSPROG

  integer*4, public, parameter :: STEP_STDOUT = 6
  integer*4, public :: FSTEP_CUR
  integer*4, public :: FSTEP_TOT
  integer*4, public :: FSTEP_SHW
  integer*4, public :: FSTEP_SHN
  integer*4, public :: FSTEP_RUN

CONTAINS



! ---------------------------------------------------------------------
subroutine CSPROG_INIT(ntot,nshow)

  implicit none

  integer*4, intent(in) :: ntot, nshow
  integer*4 :: i
  
  if (FSTEP_RUN/=0) then
    call CSPROG_FINISH()
  end if
  
  FSTEP_CUR = 0
  FSTEP_TOT = ntot
  FSTEP_SHW = nshow
  FSTEP_SHN = 0
  
  write(unit=STEP_STDOUT,fmt='(A,$)') "0"
  do i=1, nshow-1
    write(unit=STEP_STDOUT,fmt='(A,$)') "-"
  end do
  write(unit=STEP_STDOUT,fmt='(A)') "1"
  write(unit=STEP_STDOUT,fmt='(A,$)') "*"
  FSTEP_RUN = 1
  
  return

end subroutine CSPROG_INIT
! ---------------------------------------------------------------------




! ---------------------------------------------------------------------
subroutine CSPROG_FINISH()

  implicit none
  
  integer*4 :: i
  
  if (FSTEP_RUN==0) return
  
  if (FSTEP_SHN==FSTEP_SHW) then
    FSTEP_RUN = 0
    FSTEP_CUR = FSTEP_TOT
    write(unit=STEP_STDOUT,fmt='(A)') " "
    return
  end if
  
  do i=FSTEP_SHN+1, FSTEP_SHW
    write(unit=STEP_STDOUT,fmt='(A,$)') "*"
    FSTEP_SHN = FSTEP_SHN + 1
  end do
  FSTEP_RUN = 0
  FSTEP_CUR = FSTEP_TOT
  write(unit=STEP_STDOUT,fmt='(A)') " "
  
  return

end subroutine CSPROG_FINISH
! ---------------------------------------------------------------------




! ---------------------------------------------------------------------
subroutine CSPROG_UPDATE()

  implicit none
  
  integer*4 :: ifract
  real*4 :: fract
  
  if (FSTEP_RUN==0) return
    
  fract = real(FSTEP_CUR*FSTEP_SHW)/real(FSTEP_TOT)
  ifract = int(fract)
  
  if (ifract<=FSTEP_SHN) return
  
  do while (ifract>FSTEP_SHN)
    write(unit=STEP_STDOUT,fmt='(A,$)') "*"
    FSTEP_SHN = FSTEP_SHN + 1
  end do
  
  return
  
end subroutine CSPROG_UPDATE
! ---------------------------------------------------------------------

end module CSPROG