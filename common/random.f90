!**********************************************************************!
!**********************************************************************!
!                                                                      !
!                    file   "random".f90                               !
!                                                                      !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                  !
!    Version  :  1.0.0, November 21, 2003                              !
!                                                                      !
!                                                                      !
!**********************************************************************!
!                                                                       
!   Author: Juri Barthel                                                
!           Ernst Ruska-Centre                                          
!           Forschungszentrum Jülich GmbH, 52425 Jülich, Germany        
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
SUBROUTINE InitRand2(newseed)
! function: initiates the internal pseudo random number generator with
!           a given value
! -------------------------------------------------------------------- !
! parameter: integer newseed (seed number for the rng)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer, intent(inout) :: newseed
  integer :: seed_size
  integer,allocatable :: seed(:)
! ------------

! ------------
  call RANDOM_SEED() ! initialize with system generated seed
  call RANDOM_SEED(SIZE=seed_size) ! find out size of seed
  allocate(seed(seed_size))
  seed = newseed
  call RANDOM_SEED(PUT=seed) ! put new seed
  deallocate(seed)
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
  
  integer*4 :: k
  real*8 :: m, l, p
  
  m = dabs(dble(mean))
  
  ! Donold Knuth's algorithm is fast for small mean values < 30
  l = dexp(-m)
  k = 0
  p = 1.D+0
  do while (p > l)
    k = k + 1
    p = p * dble( UniRand() )
  end do
  PoissonRand = k-1
  return

  end FUNCTION PoissonRand
!**********************************************************************!



  
!**********************************************************************72
!
! RNG_DISTR_REJ returns a pseudorandom R4 following a probability
!               distribution function (PDF) given by samples (x, y)
!               in an array of length n.
!
!  Discussion:
!
!    This routine implements the rejection sampling method.
!
!    It makes use of function urnd to generate uniform random numbers
!    for the range of samples x and y.
!
!    A random is x is accepted as return value if the random y is smaller
!    than or equal to the value of y at x.
!
!    Linear interpolation is used between samples of x.
!
!    x is expected to be in ascending order
!    with corresponding y = PDF(x).
!
!    It is not required to input a normalized PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2023
!
!  Author:
!
!    Juri Barthel
!
!  Reference:
!
!    rejection sampling:  https://en.wikipedia.org/wiki/Rejection_sampling
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of samples in following input.
!    N must be greater than 0. N = 1 assumes a delta ditribution.
!
!    Input, real ( kind = 4 , dimension = (2,N) ) X,Y samples
!               of Y=PDF(X) corresponding to the input X.
!
!    Output, real ( kind = 4 ) RNG_DISTR_REJ, a new pseudorandom variate,
!    strictly in the range of minval(pdf(1,:)) and maxval(pdf(1,:)).
!
!
REAL*4 FUNCTION rng_distr_rej( n, pdf )
  
	implicit none
    
  integer*4, parameter :: num_rej_max = 1048576

	integer*4, intent(in) :: n
  real*4, intent(in), dimension(2,n) :: pdf
	
  real*4 :: x0, x1, y0, y1 ! x and y ranges
  real*4 :: urn_x, urn_y ! uniform random numbers
  real*4 :: pdf_urn_x ! interpolated PDF value at urn_x
  integer*4 :: ip ! interpolation base
  integer*4 :: irej ! rejection count
  
  real*4, external :: urnd
    
  if ( n == 0 ) then
		write ( *, '(a)' ) ' '
		write ( *, '(a)' ) 'RNG_DISTR_REJ - Fatal error!'
		write ( *, '(a)' ) '  Input value of N = 0.'
		stop 1
  end if
    
  if ( n == 1 ) then ! delta distribution
      
    rng_distr_rej = pdf(1,1)
      
  end if
    
  if ( n > 1 ) then
      
    ! determine range of x samples
    x0 = minval(pdf(1,:))
    x1 = maxval(pdf(1,:))
    ! determine range of y samples
    y0 = minval(pdf(2,:))
    y1 = maxval(pdf(2,:))
    !
    if (y0 < 0.0D+0 .or. y1 < 0.0D+0 ) then
      write ( *, '(a)' ) ' '
			write ( *, '(a)' ) 'RNG_DISTR_REJ - Fatal error!'
			write ( *, '(a)' ) '  Negative probability input encountered.'
			stop 2
    end if
      
    if (x1 > x0) then ! good range in x
        
      if (y1 > y0) then ! good range in y
          
        ! implement rejection sampling
        do irej=1, num_rej_max
            
          ! get a random x, uniform sample of the x range
          urn_x = urnd(x0, x1)
          ! get a random y, uniform sample of the range 0 to max(y)
          urn_y = urnd(0.0, y1)
            
          ! get the value of the PDF a urn_x by linear interpolation
          ! on the samples (x, y)
          pdf_urn_x = 0.0 ! init some value of the PDF, use zero to catch problems
          do ip=1, n
              
            ! find first x sample that is bigger then urn_x
            if ( pdf(1,ip) > urn_x ) then
                
              ! handle strange cases first
              if ( ip == 1 ) then
                  
                pdf_urn_x = pdf(2,1)
                  
                exit ! exit loop over ip
                  
              else if ( ip > 1 .and. ip <= n ) then ! in interpolation range
                  
                ! linear interpolation
                pdf_urn_x = pdf(2,ip-1) + ( pdf(2,ip) - pdf(2,ip-1) )  &
                    / ( pdf(1,ip) - pdf(1,ip-1) ) * ( urn_x - pdf(1,ip-1) )
                  
                exit ! exit loop over ip
                  
              end if
                
            end if
              
          end do
            
          ! loopexit, check whether this was because none of the if conditions aplied
          if (ip > n) then
              
              pdf_urn_x = pdf(2,n) ! use end value
               
          end if
            
          ! acceptance ?
          if ( urn_y <= pdf_urn_x ) then ! yes, accept this random x
              
            rng_distr_rej = urn_x
              
            return
              
          end if
            
        end do
          
        ! if the code gets to this point the number of rejections
        ! is too large. Something is fishy here. Report a problem.
        write ( *, '(a)' ) ' '
			  write ( *, '(a)' ) 'RNG_DISTR_REJ - Fatal error!'
			  write ( *, '(a)' ) '  Too many rejections.'
			  stop 3
          
      else ! invalid range in y, assuming uniform distribution
          
        rng_distr_rej = urnd(x0, x1)
          
        return
          
      end if
        
    else ! invalid range in x, assuming delta distribution
        
      rng_distr_rej = pdf(1,1)
        
      return
        
    end if
  end if
  
end function
  

!**********************************************************************!
!**********************************************************************!
