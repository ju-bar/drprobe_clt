!*********************************************************************!
!
! J. Barthel, (c) 2017, RWTH Aachen
!
! ju.barthel-at-fz-juelich.de
!
!*********************************************************************!

!*********************************************************************!
!
! Numerical integration routines
!
! - implementation in double precision for real variables
!
! - trapezoid with simpsons rule for 1d integrals
! - simple serial integration for 2d integrals
!
!*********************************************************************!
!----------------------------------------------------------------------
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
!----------------------------------------------------------------------

!*********************************************************************!
!
! link dependencies:
!
!*********************************************************************!


!*********************************************************************!
!
! History
!
! JB (180529)
! - implementation of dsgrid2d
!   grid based integration (routine used by celslc)
! JB (170705)
! - implementations of dtrapez and dsimpson of trapezoidal and simpson
!   rule integrations
!
!*********************************************************************!







!*********************************************************************!
!
! 1D Integration routines
!
!*********************************************************************!



!*********************************************************************!
!
! SUBROUTINE dtrapez
!
! This routine computes the n-th iteration of an extended
! trapezoidal integration.
! https://en.wikipedia.org/wiki/Trapezoidal_rule
!
! foo(x) is the name of the function to be integrated between limits
! a and b.
! n=1 initializes and returns s = (foo(a)+foo(b))/2/(b-a)
! Subsequent calls with n=2,3,... (in this order) will
! improve the accuracy of the integral s by adding 2n-2 additional
! interior points.
! Do not modify s between calls of increasing n.
!
SUBROUTINE dtrapez(foo,a,b,s,n)
  
  implicit none
  
  integer*4, intent(in) :: n
  real*8, intent(in) :: a, b
  real*8, intent(inout) :: s
  
  real*8, external :: foo

  integer*8 :: j, i
  real*8 :: x, sum, delta, renorm
    
  if (n==1) then ! initialize and 1st trapezoid
    s   = 0.5d+0 * (b-a)*(foo(a)+foo(b))
  else ! evaluations for n>1 assuming previous call with n--
    i  = 2**(n-2) ! number of new samples of integration 1, 2, 4 ...
	renorm = 1.0d+0 / dble(i)
	delta = (b-a) * renorm ! spacing between new points of integration
	x   = a + 0.5d+0 * delta ! first new point
	sum = 0.0d+0
	do j=1,i
	  sum = sum + foo(x)
	  x = x + delta
	end do
	s   = 0.5d+0*(s+(b-a)*sum*renorm) ! better s value
  end if
  return

  END SUBROUTINE dtrapez
!*********************************************************************!



!*********************************************************************!
!
! SUBROUTINE dsimpson
!
! Returns s as the integral of the function foo from a to b applying
! Simpson's rule.
! https://en.wikipedia.org/wiki/Simpson%27s_rule
!
! foo(x) is the name of the function to be integrated between a and b.
! Set EPSIL to the desired fractional accuracy.
! Set 2**(ITMAX-1) is the maximum allowed number of steps.
! Integration is performed by Simpson’s rule.
!
SUBROUTINE dsimpson(foo,a,b,s)
  
  implicit none
  
  real*8, intent(in) :: a, b
  real*8, intent(inout) :: s
  
  real*8, external :: foo
  external :: dtrapez

  ! Change ITMAX to avoid long non-converging runs.
  ! ITMAX = 20 means 2^18 = 262144 points added
  integer*4, parameter :: ITMAX = 20
  ! Change EPSIL to set target convergence (relative)
  real*8, parameter :: EPSIL = 1.0d-06
  
  integer*4 :: i
  real*8 :: pres, prest, stra

  ! initialize old s values with some strange values  
  prest = -1.d+30
  pres  = -1.d+30
  do i=1,ITMAX
    call dtrapez(foo,a,b,stra,i) ! update trapezoidal estimate
    s=(4.d+0 * stra - prest)/3.d+0  ! update simpsons estimate
    if (i > 5) then ! loop at least 4 times, then check for convergence
      if (dabs(s-pres)<EPSIL*dabs(pres).or. &
      &  (pres==0.d+0.and.s==0.d+0)) then
        return
      end if
    end if
    pres = s ! store current estimate of s
    prest = stra ! store current trapezoidal estimate of s
  end do
  return

  END SUBROUTINE dsimpson
!*********************************************************************!








!*********************************************************************!
!
! 2D integration routines
!
!*********************************************************************!



!*********************************************************************!
!
! SUBROUTINE dtrapz2d1
!
! This routine computes the nth stage of reﬁnement of an extended
! trapezoidal rule on dimension 1 of a 2d function.
!
! foo(x,y) is the name of the function to be integrated
! between limits a and b on x for constant y.
! When called with n=1, the routine returns initializes with 1 trapez
! and returns a rough estimate s.
! Further calls with n=2,3,... improve s by adding 2n-2 additional
! interpolation samples.
! Do not modify s between calls!
!
SUBROUTINE dtrapz2d1(foo,y,a,b,s,n)
  
  implicit none
  
  integer*4, intent(in) :: n
  real*8, intent(in) :: a, b, y
  real*8, intent(inout) :: s
  
  real*8, external :: foo

  integer*8 :: j, i
  real*8 :: x, delta, renorm, sum
    
  if (n==1) then ! initial
    s   = 0.5d+0 * (b-a)*(foo(a,y)+foo(b,y))
  else ! evaluation for n>1
    i  = 2**(n-2) ! number of new samples of integration 1, 2, 4 ...
	renorm = 1.0d+0 / dble(i)
	delta = (b-a) * renorm ! spacing between new points of integration
	x   = a + 0.5d+0 * delta ! first new point
	sum = 0.0d+0
	do j=1,i
	  sum = sum + foo(x,y)
	  x = x + delta
	end do
	s   = 0.5d+0*(s+(b-a)*sum*renorm) ! better s value
  end if
  return

END SUBROUTINE dtrapz2d1
!*********************************************************************!



!*********************************************************************!
!
! SUBROUTINE dtrapz2d2
!
! This routine computes the nth stage of reﬁnement of an extended
! trapezoidal rule on dimension 2 of a 2d function.
!
! foo(x,y) is the name of the function to be integrated
! between limits a and b on x for constant y.
! When called with n=1, the routine returns initializes with 1 trapez
! and returns a rough estimate s.
! Further calls with n=2,3,... improve s by adding 2n-2 additional
! interpolation samples.
! Do not modify s between calls!
!
SUBROUTINE dtrapz2d2(foo,x,a,b,s,n)
  
  implicit none
  
  integer*4, intent(in) :: n
  real*8, intent(in) :: a, b, x
  real*8, intent(inout) :: s
  
  real*8, external :: foo

  integer*8 :: j, i
  real*8 :: y, delta, renorm, sum
    
  if (n==1) then ! initial
    s   = 0.5d+0 * (b-a)*(foo(x,a)+foo(x,b))
  else ! evaluation for n>1
    i  = 2**(n-2) ! number of new samples of integration 1, 2, 4 ...
	renorm = 1.0d+0 / dble(i)
	delta = (b-a) * renorm ! spacing between new points of integration
	y   = a + 0.5d+0 * delta ! first new point
	sum = 0.0d+0
	do j=1,i
	  sum = sum + foo(x,y)
	  y = y + delta
	end do
	s   = 0.5d+0*(s+(b-a)*sum*renorm) ! better s value
  end if
  return

END SUBROUTINE dtrapz2d2
!*********************************************************************!



!*********************************************************************!
!
! FUNCTION dsimp2d1(foo,y,a,b)
!
! Returns the integral of the function foo(x,y) over the
! range a <= x <= b for a fix y.
!
! foo is input as a function of two variables : foo(x,y)
!
! The result of integration is returned as function value.
!
FUNCTION dsimp2d1(foo,y,a,b)

  implicit none
  
  real*8, intent(in) :: y, a, b
  real*8, external :: foo
  
  external :: dtrapz2d2
  
  real*8 :: dsimp2d1
  
  ! Change ITMAX to avoid long non-converging runs.
  ! ITMAX = 20 means 2^18 = 262144 points added
  integer*4, parameter :: ITMAX = 20
  ! Change EPSIL to set target convergence (relative)
  real*8, parameter :: EPSIL = 1.0d-06
  
  integer*4 :: j
  real*8 :: s, strap, pres, prest

  ! initialize old s values with some strange values  
  prest = -1.d+30
  pres  = -1.d+30
  do j=1,ITMAX
    call dtrapz2d1(foo,y,a,b,strap,j)
    s=(4.d+0 * strap - prest)/3.d+0 
    if (j > 5) then ! Avoid spurious early convergence.
      if (  dabs(s-pres)<EPSIL*dabs(pres) .or. &
       &   (s==0.d+0 .and. pres==0.d+0)) then ! converged
        dsimp2d1 = s ! terminate
        return
      end if
    end if
    pres = s
    prest = strap
  end do
  dsimp2d1 = s ! terminate without reaching convergence
  
  return
  
END FUNCTION dsimp2d1
!*********************************************************************!



!*********************************************************************!
!
! FUNCTION dsimp2d2(foo,x,a,b)
!
! Returns the integral of the function foo(x,y) over the
! range a <= y <= b for a fix x.
!
! foo is input as a function of two variables : foo(x,y)
!
! The result of integration is returned as function value.
!
FUNCTION dsimp2d2(foo,x,a,b)

  implicit none
  
  real*8, intent(in) :: x, a, b
  real*8, external :: foo
  
  external :: dtrapz2d2
  
  real*8 :: dsimp2d2
  
  ! Change ITMAX to avoid long non-converging runs.
  ! ITMAX = 20 means 2^18 = 262144 points added
  integer*4, parameter :: ITMAX = 20
  ! Change EPSIL to set target convergence (relative)
  real*8, parameter :: EPSIL = 1.0d-06
  
  integer*4 :: j
  real*8 :: s, strap, pres, prest

  ! initialize old s values with some strange values  
  prest = -1.d+30
  pres  = -1.d+30
  do j=1,ITMAX
    call dtrapz2d2(foo,x,a,b,strap,j)
    s=(4.d+0 * strap - prest)/3.d+0 
    if (j > 5) then ! Avoid spurious early convergence.
      if (  dabs(s-pres)<EPSIL*dabs(pres) .or. &
       &   (s==0.d+0 .and. pres==0.d+0)) then ! converged
        dsimp2d2 = s ! terminate
        return
      end if
    end if
    pres = s
    prest = strap
  end do
  dsimp2d2 = s ! terminate without reaching convergence
  
  return
  
END FUNCTION dsimp2d2
!*********************************************************************!


!*********************************************************************!
!
! SUBROUTINE dtrapz2do1
!
! This routine computes the nth stage of reﬁnement of an extended
! trapezoidal rule for the outer integral of a 2d integration done
! on the first variable of a function foo(x,y).
!
! foo(x,y) is input as the name of the function to be integrated
! between limits a1 and b1 for x and a2 and b2 for y, also input.
! When called with n=1, the routine returns as s the crudest estimate
! of Int_a1^b1 [ dx Int_a2_b2 dy f(x,y) ]. Where the inner integral
! is also evaluated.
! Subsequent calls with n=2,3,... (in that sequential order) will
! improve the accuracy of s by adding 2n-2 additional interior points.
! s should not be modiﬁed between sequential calls.
!
SUBROUTINE dtrapz2do1(foo,a1,b1,a2,b2,s,n)
  
  implicit none
  
  integer*4, intent(in) :: n
  real*8, intent(in) :: a1, b1, a2, b2
  real*8, intent(inout) :: s
  
  real*8, external :: foo ! 2d function reference
  real*8, external :: dsimp2d2 ! reference of function doing inner integration over y

  integer*8 :: j, i
  real*8 :: x, delta, renorm, sum
    
  if (n==1) then ! initial and crude
    s   = 0.5d+0 * (b1-a1)*(dsimp2d2(foo,a1,a2,b2)+dsimp2d2(foo,b1,a2,b2))
  else ! evaluation for n>1
    i  = 2**(n-2) ! number of new samples of integration 1, 2, 4 ...
	renorm = 1.0d+0 / dble(i)
	delta = (b1-a1) * renorm ! spacing between new points of integration
	x   = a1 + 0.5d+0 * delta ! first new point
	sum = 0.0d+0
	do j=1,i
	  sum = sum + dsimp2d2(foo,x,a2,b2)
	  x = x + delta
	end do
	s   = 0.5d+0*(s+(b1-a1)*sum*renorm) ! refined s value
  end if
  return

END SUBROUTINE dtrapz2do1
!*********************************************************************!


!*********************************************************************!
!
! SUBROUTINE dtrapz2do2
!
! This routine computes the nth stage of reﬁnement of an extended
! trapezoidal rule for the outer integral of a 2d integration done
! on the second variable of a function foo(x,y).
!
! foo(x,y) is input as the name of the function to be integrated
! between limits a1 and b1 for x and a2 and b2 for y, also input.
! When called with n=1, the routine returns as s the crudest estimate
! of Int_a2^b2 [ dy Int_a1_b1 dx f(x,y) ]. Where the inner integral
! is also evaluated.
! Subsequent calls with n=2,3,... (in that sequential order) will
! improve the accuracy of s by adding 2n-2 additional interior points.
! s should not be modiﬁed between sequential calls.
!
SUBROUTINE dtrapz2do2(foo,a1,b1,a2,b2,s,n)
  
  implicit none
  
  integer*4, intent(in) :: n
  real*8, intent(in) :: a1, b1, a2, b2
  real*8, intent(inout) :: s
  
  real*8, external :: foo ! 2d function reference
  real*8, external :: dsimp2d1 ! reference of function doing inner integration over x

  integer*8 :: j, i
  real*8 :: y, delta, renorm, sum
    
  if (n==1) then ! initial and crude
    s   = 0.5d+0 * (b2-a2)*(dsimp2d1(foo,a2,a1,b1)+dsimp2d1(foo,b2,a1,b1))
  else ! evaluation for n>1
    i  = 2**(n-2) ! number of new samples of integration 1, 2, 4 ...
	renorm = 1.0d+0 / dble(i)
	delta = (b2-a2) * renorm ! spacing between new points of integration
	y   = a2 + 0.5d+0 * delta ! first new point
	sum = 0.0d+0
	do j=1,i
	  sum = sum + dsimp2d1(foo,y,a1,b1)
	  y = y + delta
	end do
	s   = 0.5d+0*(s+(b2-a2)*sum*renorm) ! refined s value
  end if
  return

END SUBROUTINE dtrapz2do2
!*********************************************************************!



!*********************************************************************!
!
! SUBROUTINE dqsimp2d
!
! Returns s as the integral of the function foo(x,y) over the
! simple ranges a1 <= x <= b1 and a2 <= y <= b2.
!
! The inner integral will be over x if d=1 otherwise over y.
! d == 1:
! s = Int_a2^b2[ dy * Int_a1^b1 [ dx * foo(x,y) ] ]
! d /= 1:
! s = Int_a1^b1[ dx * Int_a2^b2 [ dy * foo(x,y) ] ]
!
! foo(x) is input as the name of the function to be integrated
! between limits a and b, also input.
! The parameters EPSIL can be set to the desired fractional accuracy
! and ITMAX so that 2 to the power ITMAX-1 is the maximum allowed number
! of steps. Integration is performed by Simpson’s rule.
!
SUBROUTINE dqsimp2d(foo,a1,b1,a2,b2,d,s)
  
  implicit none
  
  integer*4, intent(in) :: d
  real*8, intent(in) :: a1, b1, a2, b2
  real*8, intent(inout) :: s
  
  real*8, external :: foo
  external :: dtrapz2do2

  ! Change ITMAX to avoid long non-converging runs.
  ! ITMAX = 20 means 2^18 = 262144 points added
  integer*4, parameter :: ITMAX = 20
  ! Change EPSIL to set target convergence (relative)
  real*8, parameter :: EPSIL = 1.0d-06
  
  integer*4 :: j
  real*8 :: pres, prest, strap

  ! initialize old s values with some strange values  
  prest = -1.d+30
  pres  = -1.d+30
  do j=1,ITMAX
    if (d==2) then ! outer integral is over x
      call dtrapz2do1(foo,a1,b1,a2,b2,strap,j)
    else ! outer integral is over y (default)
      call dtrapz2do2(foo,a1,b1,a2,b2,strap,j)
    end if
    s=(4.d+0 * strap - prest)/3.d+0 
    if (j > 5) then ! Avoid spurious early convergence.
      if (  dabs(s-pres)<EPSIL*dabs(pres) .or. &
       &   (s==0.d+0 .and. pres==0.d+0)) then
        return
      end if
    end if
    pres = s
    prest = strap
  end do
  return

END SUBROUTINE dqsimp2d
!*********************************************************************!



!*********************************************************************!
!
! FUNCTION dsgrid1d
!
! Integrates a 1d function foo on a grid applying power law sampling
! between the integration boundaries. The summation on the grid is
! done using the trapezoidal rule.
!
FUNCTION dsgrid1d(foo,a1,b1,p1,n1)

  implicit none
  
  integer*4, parameter :: nmin = 8 ! min. grid dimension
  integer*4, parameter :: nmax = 2048 ! max. grid dimension
  real*8, parameter :: pmin = 1.0d-1 ! min. grid sampling power
  real*8, parameter :: pmax = 1.0d+1 ! max. grid sampling power
  
  real*8, external :: foo  ! integrand function of one variable
                           ! foo(x1)
  
  real*8 :: dsgrid1d ! function type
                                  
  real*8, intent(in) :: a1, b1 ! integral boundaries (closed)
                           ! a1 <= x1 <= b1
  real*8, intent(in) :: p1 ! grid sampling power law
                           ! x1_i = a1 + (b1-a1)*((i-1)/(n1-1))**p1
  integer*4, intent(in) :: n1 ! grid size
  
  integer*4 :: i1 ! iterator
  integer*4 :: ni1 ! internal array dimensions
  integer*4 :: nalloc ! allocation status
  real*8 :: cy, cx1 ! temp values
  real*8 :: pi1 ! internal sampling power
  
  real*8, allocatable, dimension(:) :: tx1 ! x samples
  real*8, allocatable, dimension(:) :: y ! integrand samples
  
  ! init
  dsgrid1d = 0.0d+0
  nalloc = 0
  ni1 = max(nmin,min(nmax,n1)) ! put array dimensions into limits
  pi1 = max(pmin,min(pmax,p1)) ! put sampling powers into limits
  allocate(tx1(ni1), stat=nalloc ) ! allocate x array
  if (nalloc/=0) return
  ! fill the grid array
  cx1 = 1.d+0 / dble(ni1-1)
  do i1=1, ni1
    tx1(i1) = a1 + (b1-a1)*(dble(i1-1)*cx1)**pi1
  end do
  allocate(y(ni1), stat=nalloc) ! allocate integrand array
  if (nalloc/=0) return
  y = 0.0d+0
  ! fill the integrand array
  do i1=1, ni1
    cx1 = tx1(i1)
    y(i1) = foo(cx1)
  end do
  ! integrate by trapezoidal summation
  cy = 0.0d+0
  ! integral over x1
  do i1=1, ni1-1
    cy = cy + (y(i1+1)+y(i1))*(tx1(i1+1)-tx1(i1))
  end do
  cy = cy * 0.5d+0 ! half-sum (trapez)
  dsgrid1d = cy
  ! finish
  deallocate(tx1,y, stat=nalloc)
  return

END FUNCTION dsgrid1d



!*********************************************************************!
!
! FUNCTION dsgrid2d
!
! Integrates a 2d function foo on a grid applying power law sampling
! between the integration boundaries. The summation on the grid is
! done using the trapezoidal rule.
!
FUNCTION dsgrid2d(foo,a1,b1,a2,b2,p1,p2,n1,n2)

  implicit none
  
  integer*4, parameter :: nmin = 8 ! min. grid dimension
  integer*4, parameter :: nmax = 2048 ! max. grid dimension
  real*8, parameter :: pmin = 1.0d-1 ! min. grid sampling power
  real*8, parameter :: pmax = 1.0d+1 ! max. grid sampling power
  
  real*8, external :: foo  ! integrand function of two variables
                           ! foo(x1,x2)
  
  real*8 :: dsgrid2d ! function type
                                  
  real*8, intent(in) :: a1, b1, a2, b2 ! integral boundaries (closed)
                           ! a1 <= x1 <= b1, a2 <= x2 <= b2
  real*8, intent(in) :: p1, p2 ! grid sampling power law
                           ! xj_i = aj + (bj-aj)*((i-1)/(nj-1))**pj
  integer*4, intent(in) :: n1, n2 ! grid sizes for the two dimensions
  
  integer*4 :: i1, i2 ! iterators
  integer*4 :: ni1, ni2 ! internal array dimensions
  integer*4 :: nalloc ! allocation status
  real*8 :: cy, cx1, cx2 ! temp values
  real*8 :: pi1, pi2 ! internal sampling power
  
  real*8, allocatable, dimension(:) :: tx1, tx2 ! x_j samples
  real*8, allocatable, dimension(:) :: sx2 ! partial integral
  real*8, allocatable, dimension(:,:) :: y ! integrand samples
  
  ! init
  dsgrid2d = 0.0d+0
  nalloc = 0
  ni1 = max(nmin,min(nmax,n1)) ! put array dimensions into limits
  ni2 = max(nmin,min(nmax,n2)) ! ...
  pi1 = max(pmin,min(pmax,p1)) ! put sampling powers into limits
  pi2 = max(pmin,min(pmax,p2)) ! ...
  allocate(tx1(ni1),tx2(ni2),sx2(ni2) , stat=nalloc ) ! allocate xj arrays
  if (nalloc/=0) return
  ! fill the grid arrays
  cx1 = 1.d+0 / dble(ni1-1)
  do i1=1, ni1
    tx1(i1) = a1 + (b1-a1)*(dble(i1-1)*cx1)**pi1
  end do
  cx2 = 1.d+0 / dble(ni2-1)
  do i2=1, ni2
    tx2(i2) = a2 + (b2-a2)*(dble(i2-1)*cx2)**pi2
  end do
  allocate(y(ni1,ni2), stat=nalloc) ! allocate integrand array
  if (nalloc/=0) return
  y = 0.0d+0
  ! fill the integrand array
  do i2=1, ni2
    cx2 = tx2(i2)
    do i1=1, ni1
      cx1 = tx1(i1)
      y(i1,i2) = foo(cx1,cx2)
    end do
  end do
  ! integrate by trapezoidal summation
  cy = 0.0d+0
  sx2 = 0.0d+0
  ! inner integral over x1 / fast memory sequence
  do i2=1, ni2
    do i1=1, ni1-1
      sx2(i2) = sx2(i2) + (y(i1+1,i2)+y(i1,i2))*(tx1(i1+1)-tx1(i1))
    end do
  end do
  sx2 = sx2 * 0.5d+0 ! half-sum (trapez)
  ! outer integral over x2
  do i2=1, ni2-1
    cy = cy + (sx2(i2+1)+sx2(i2))*(tx2(i2+1)-tx2(i2))
  end do
  cy = cy * 0.5d+0 ! half-sum (trapez)
  dsgrid2d = cy
  ! finish
  deallocate(tx1,tx2,sx2,y, stat=nalloc)
  return

END FUNCTION dsgrid2d



!*********************************************************************!
!*********************************************************************!

