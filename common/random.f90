!**********************************************************************!
!**********************************************************************!
!                                                                      !
!                    file   "random".f95                               !
!                                                                      !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                  !
!    Version  :  1.0.0, November 21, 2003                              !
!    partial code source: Numerical Recipies for F77 and F90           !
!                                                                      !
!                                                                      !
!**********************************************************************!
!                                                                       
!   Author: Juri Barthel                                                
!           Ernst Ruska-Centre                                          
!           Forschungszentrum Jülich GmbH, 52425 Jülich, Germany        
!           RWTH Aachen University, 52074 Aachen, Germany               
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
! This program is free software: you can redistribute it and/or modify  
! it under the terms of the GNU General Public License as published by  
! the Free Software Foundation, either version 3 of the License, or     
! (at your option) any later version.                                   
!                                                                       
! This program is distributed in the hope that it will be useful,       
! but WITHOUT ANY WARRANTY; without even the implied warranty of        
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
! GNU General Public License for more details.                          
!                                                                       
! You should have received a copy of the GNU General Public License     
! along with this program. If not, see <http://www.gnu.org/licenses/>.  
!                                                                       
!-----------------------------------------------------------------------


!**********************************************************************!
!                                                                      !
!   Purpose: implementing subroutines for random number geneartion     !
!            accessable from F90/95 AND F77                            !
!                                                                      !
!   External  sources to link:                                         !
!                                                                      !
!**********************************************************************!
!                                                                      !
!   CONTAINS:                                                          !
!      1) RANDOM Subroutines + Functions                               !
!         SUBROUTINE InitRand()                                        !
!         REAL*4 FUNCTION UniRand()                                    !
!         REAL*4 FUNCTION GaussRand()                                  !
!         integer*4 FUNCTION PoissonRand(mean)                         !
!                                                                      !
!**********************************************************************!








!**********************************************************************!
!*********************** RANDOM Subroutines ***************************!
!**********************************************************************!











!**********************************************************************!
!**********************************************************************!
SUBROUTINE InitRand()
! function: initiates the internal pseudo random number generator with
!           a processor dependent value
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer :: seed_size, newseed, values(8)
  integer,allocatable :: seed(:)
  character*8 :: cldate
  character*10 :: cltime
  character*5 :: clzone
! ------------





! ------------
  call RANDOM_SEED() ! initialize with system generated seed
	
  call RANDOM_SEED(SIZE=seed_size) ! find out size of seed
  allocate(seed(seed_size))

  call date_and_time(cldate,cltime,clzone,values)
	
  newseed=values(1)+values(2)+values(3)
  newseed = newseed+values(7)*1000+values(6)*60000+values(5)*3600000
  seed = newseed
  call RANDOM_SEED(PUT=seed) ! put new seed

  deallocate(seed)           ! safe
! ------------



! ------------
  RETURN

END SUBROUTINE InitRand
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE InitRand2()
! function: initiates the internal pseudo random number generator with
!           a processor dependent value
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer :: seed_size, newseed, values(8)
  integer,allocatable :: seed(:)
  character*8 :: cldate
  character*10 :: cltime
  character*5 :: clzone
! ------------

! ------------
  call RANDOM_SEED() ! initialize with system generated seed
  call RANDOM_SEED(SIZE=seed_size) ! find out size of seed
  allocate(seed(seed_size))
  call date_and_time(cldate,cltime,clzone,values)
	
  newseed = values(1)+values(2)+values(3)
  newseed = newseed+values(7)*1000+values(6)*60000+values(5)*3600000
  seed = newseed
  call RANDOM_SEED(PUT=seed) ! put new seed

  deallocate(seed)           ! safe
! ------------

! ------------
  RETURN

END SUBROUTINE InitRand2
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE InitRandFix()
! function: initiates the internal pseudo random number generator with
!           a fix seed of 42
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer :: seed_size, newseed
  integer,allocatable :: seed(:)
! ------------

! ------------
  call RANDOM_SEED() ! initialize with system generated seed
  call RANDOM_SEED(SIZE=seed_size) ! find out size of seed
  allocate(seed(seed_size))
  newseed = 42
  seed = newseed
  call RANDOM_SEED(PUT=seed) ! put new seed
  deallocate(seed)           ! safe
! ------------

! ------------
  RETURN

END SUBROUTINE InitRandFix
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
REAL*4 FUNCTION UniRand()
! function: returns a uniformly distributed pseudo random number in the
!           (REAL*4)-interval 0.0.le.x.lt.1.0
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !


  implicit none
  
! ------------
! declaration
	REAL*4 :: harvest
! ------------


! ------------
	CALL RANDOM_NUMBER(harvest)
	UniRand = harvest
! ------------	



! ------------	
	RETURN

END FUNCTION UniRand
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
REAL*4 FUNCTION urnd(r0,r1)
! function: returns a uniformly distributed pseudo random number in the
!           (REAL*4)-interval r0.le.x.ge.r1
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !


  implicit none
  
! ------------
! declaration
	REAL*4 :: harvest
	REAL*4, intent(in) :: r0, r1
	REAL*4 :: ir0, ir1
! ------------


! ------------
    ir0 = r0
    ir1 = r1
    if (ir0.gt.ir1) then
      ir0 = r1
      ir1 = r0
    end if
	CALL RANDOM_NUMBER(harvest)
	urnd = ir0+harvest*(ir1-ir0)
! ------------	



! ------------	
	RETURN

END FUNCTION urnd
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
REAL*4 FUNCTION GaussRand()
! function: returns a normally distributed pseudo random number in
!           (REAL*4), variance is 1.0, mean is 0.0
!           using box mueller transformation, see Numerical Recipies 7.2
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none
  
  
! ------------
! declaration
	real*4 :: v1,v2,rsq,fac,gset
	integer*4 :: iset
	real*4, external :: UniRand
	SAVE iset, gset
	DATA iset /0/
! ------------


! ------------
  if (iset==0) then

! ------------
! get two uniformly distributed PRN within unit circle
	DO
		v1 = 2.0*UniRand()-1.0
		v2 = 2.0*UniRand()-1.0
		rsq = v1**2.0+v2**2.0
		IF ((rsq<=1.0).AND.(rsq>0.0)) EXIT
	END DO

!	write(*,*) 'GaussRand: v1', v1
!	write(*,*) 'GaussRand: v2', v2
!	write(*,*) 'GaussRand: rsq', rsq
! ------------




! ------------
! Do box-mueller transformation now
! one additional Gaussian random number can be createt, but its skipped
! in this version of this function
	fac = SQRT( -2.0*LOG(rsq)/rsq )
	
!	write(*,*) 'GaussRand: fac', fac
!	write(*,*) 'GaussRand: result', v1*fac
	
	GaussRand = v1 * fac
! 2nd-GaussRan
    gset = v2 * fac	
    iset = 1
! ------------

  else ! iset > 0
    
    GaussRand = gset
    iset = 0
  
  end if




! ------------	
	RETURN

END FUNCTION GaussRand
!**********************************************************************!




!**********************************************************************!
!
! integer*4 FUNCTION PoissonRand
!
! adds poisson noise to a given mean value number
!
! parameter: real*4 :: mean
integer*4 FUNCTION PoissonRand(mean)  

  implicit none
  
  real*4, intent(in) :: mean
  real*4, external :: UniRand
  
  integer*4 :: k, k_max
  real*8 :: m, l, p, sum
  
  k = 0                             ! counter of events
  m = abs(mean)                     ! expected mean
  k_max = nint(10*m)        ! upper limit for counts
  l = exp(-m)                       ! probability
  sum = l                           ! cumulant
  p = UniRand()                     ! a uniform random number
  
  do while (sum<p .and. k<k_max)    ! loop over all events up to a random cumulant value
    k = k + 1                       ! next k
    l = l * m / real(k)             ! next probability
    sum = sum + l                   ! add to cumulant
  end do
  
  PoissonRand = k
  
  return

end FUNCTION PoissonRand
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
