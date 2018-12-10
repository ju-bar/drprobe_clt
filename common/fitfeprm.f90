!**********************************************************************!
!
! file: fitfeprm.f90
! author: j.barthel (ju.barthel-at-fz-juelich.de) 2013
!
! Implementation of routines for least-squares minimization
! fitting of a parameterization to a table of given 
! electron scattering factors.
!
! Incorporate by USE fitfeprm
!
! Linked files: 
!   random.f90
! + Watch out for dependencies therein!
!
!**********************************************************************!
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


MODULE fitfeprm

  IMPLICIT NONE
  
  SAVE
  
  INTEGER*4, PUBLIC, PARAMETER :: FFE_PRM_MAX = 12 ! using 12 function parameters
  REAL*4, PUBLIC, PARAMETER :: FFE_PREFAC = 0.0239336609804 ! translation factor m0 * e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) * 10^(-10)     ->[A]
  ! process control flags
  INTEGER*4, PUBLIC :: FFE_SOUT ! standard output channel (may be redirected to files)
  DATA FFE_SOUT /6/
  INTEGER*4, PUBLIC :: FFE_SHOW_MSG ! Control flag, set unequal zero to activate console messages
  DATA FFE_SHOW_MSG /0/
  INTEGER*4, PUBLIC :: FFE_SHOW_DBG ! Control flag, set unequal zero to activate certain debug levels
  DATA FFE_SHOW_DBG /0/
  INTEGER*4, PUBLIC :: FFE_SHOW_ERR ! Control flag, set unequal zero to activate error messages
  DATA FFE_SHOW_ERR /0/
  INTEGER*4, PUBLIC :: FFE_ERR_STOP ! Control flag, set unequal zero to cause a stop in FFE_ERROR
  DATA FFE_ERR_STOP /0/
  INTEGER*4, PUBLIC :: FFE_ERR_ALLOC
  INTEGER*4, PUBLIC :: FFE_ERR_CUR
  INTEGER*4, PUBLIC :: FFE_USE_WEIGHTS
  DATA FFE_USE_WEIGHTS /0/
  INTEGER*4, PUBLIC :: FFE_USE_RANGE
  DATA FFE_USE_RANGE /0/
  ! function parameters: {A1, B1, C1, D1, A2, B2, C2, D2, A3, B3, C3, D3}
  ! function: f(x) = SUM_i=1,N [   Ai^2 / ( Bi^2 + x^2 + ceps)
  !                              + Ci^2 * Exp( - (Di^2 + ceps) * x^2 ) ]
  ! - 3 Lorentzians(x;Ai,Bi) and 3 Gaussians(x;Ci,Di)
  REAL*4, PUBLIC :: FFE_PCUR(FFE_PRM_MAX) ! current parameter set (backup)
  REAL*4, PUBLIC :: FFE_PTMP(FFE_PRM_MAX) ! temporary parameter set
  REAL*4, PUBLIC :: FFE_FTMP              ! temporary function value
  REAL*4, PUBLIC :: FFE_DFDP(FFE_PRM_MAX) ! temporary function derivatives
  REAL*4, PUBLIC :: FFE_PRNG(2,FFE_PRM_MAX)
  REAL*4, PUBLIC :: FFE_ZRED ! Z-dZ comparison value for the extra condition
                             ! MUST BE SET BEFORE FITTING!
                             !      Sum of Lorentz Ai^2 parameters divided by
                             !      FFE_PREFAC should be equal to Z-dZ
                             !      This is important for the large s asymptotic
                             !      behaviour.
                             !      See: Weickenmeier & Kohl Acta Cryst. A 47 (1991)
  REAL*4, PUBLIC :: FFE_W_EXTRA ! weight for the extra condition
                             ! MUST BE SET BEFORE FITTING!
  DATA FFE_W_EXTRA /0/
  !
  INTEGER*4, PUBLIC :: FFE_ITMAX ! Maximum number of optimization iterations
  DATA FFE_ITMAX /2000/
  !
  REAL*4, PUBLIC, ALLOCATABLE, DIMENSION(:) :: FFE_S, FFE_F, FFE_W
  
  
  ! SUPPORTED METHODS
  INTEGER*4, PUBLIC, PARAMETER :: FFE_MET_MAX = 2
  CHARACTER(LEN=256), PUBLIC :: FFE_MET_MSG(FFE_MET_MAX)
  DATA FFE_MET_MSG(1) /"Quasi Newton (BFGS)."/
  DATA FFE_MET_MSG(2) /"Simulated Annealing Simplex."/
  
  ! CONTROL PARAMETERS FOR SPECIFIC METHODS, WITH DEFAULT INITIAL VALUES
  INTEGER*4, PUBLIC :: FFE_SAS_NTRY
  DATA FFE_SAS_NTRY /100/
  REAL*4, PUBLIC :: FFE_SAS_TRED
  DATA FFE_SAS_TRED /0.1/
  
  ! ERROR MESSAGES
  INTEGER*4, PUBLIC, PARAMETER :: FFE_ERR_MAX = 10
  CHARACTER(LEN=256), PUBLIC :: FFE_ERR_MSG(FFE_ERR_MAX)
  DATA FFE_ERR_MSG(1)  /"Memory allocation failed."/
  DATA FFE_ERR_MSG(2)  /"Memory deallocation failed."/
  DATA FFE_ERR_MSG(3)  /"Number of function parameters is out of range."/
  DATA FFE_ERR_MSG(4)  /"Invalid data table size (0)."/
  DATA FFE_ERR_MSG(5)  /"Invalid minimization method."/
  DATA FFE_ERR_MSG(6)  /"Parameter arrays have inconsitent sizes."/
  DATA FFE_ERR_MSG(7)  /"Round-off problem in FFE_LNSRCH."/
  DATA FFE_ERR_MSG(8)  /"Too many Davidon-Fletcher-Powell iterations."/
  DATA FFE_ERR_MSG(9)  /"Parameter index out of range."/
  DATA FFE_ERR_MSG(10) /"Parameter variation went out of range."/
  
  PUBLIC :: FFE_swapr, FFE_swaprv
  PUBLIC :: FFE_iminloc, FFE_imaxloc
  PUBLIC :: FFE_outerprod
  PRIVATE :: FFE_MSG, FFE_DBG, FFE_ERR
  PUBLIC :: FFE_SET_S, FFE_SET_F, FFE_SET_W, FFE_SET_PRNG
  PUBLIC :: FFE_FCHI, FFE_DFCHI
  PUBLIC :: FFE_DFP
  PRIVATE :: FFE_LNSRCH
  PUBLIC :: FFE_SAS
  PRIVATE :: FFE_AMEBSA
  PUBLIC :: FFE_FINDMIN
  PUBLIC :: FFE_FUNCPRM
  PUBLIC :: FFE_FEPRM
  PUBLIC :: FFE_FEPRMY
  
  CONTAINS



  
!**********************************************************************!
!
SUBROUTINE FFE_swapr(r1,r2)
  IMPLICIT NONE
  REAL*4, INTENT(INOUT) :: r1, r2
  REAL*4 :: rtmp
  rtmp = r2
  r2 = r1
  r1 = rtmp
END SUBROUTINE FFE_swapr

!**********************************************************************!
!
SUBROUTINE FFE_swaprv(rv1,rv2)
  IMPLICIT NONE
  REAL*4, INTENT(INOUT), DIMENSION(:) :: rv1, rv2
  REAL*4, DIMENSION(max(size(rv1),size(rv2))) :: rvtmp
  INTEGER*4 :: n
  n = max(size(rv1),size(rv2))
  rvtmp(1:n) = rv2(1:n)
  rv2(1:n) = rv1(1:n)
  rv1(1:n) = rvtmp(1:n)
END SUBROUTINE FFE_swaprv

!**********************************************************************!
!
FUNCTION FFE_iminloc(rv)
  IMPLICIT NONE
  REAL*4, DIMENSION(:), INTENT(IN) :: rv
  INTEGER*4 :: FFE_iminloc
  INTEGER*4 :: iloc(1)
  FFE_iminloc = 0
  if (size(rv>0)) then
    iloc = minloc(rv)
    FFE_iminloc = iloc(1)
  end if
END FUNCTION FFE_iminloc  

!**********************************************************************!
!
FUNCTION FFE_imaxloc(rv)
  IMPLICIT NONE
  REAL*4, DIMENSION(:), INTENT(IN) :: rv
  INTEGER*4 :: FFE_imaxloc
  INTEGER*4 :: iloc(1)
  FFE_imaxloc = 0
  if (size(rv>0)) then
    iloc = maxloc(rv)
    FFE_imaxloc = iloc(1)
  end if
END FUNCTION FFE_imaxloc  

!**********************************************************************!
!
FUNCTION FFE_outerprod(a,b)
  IMPLICIT NONE
  REAL*4, DIMENSION(:), INTENT(IN) :: a,b
  REAL*4, DIMENSION(size(a),size(b)) :: FFE_outerprod
  FFE_outerprod = spread(a,dim=2,ncopies=size(b)) * &
                & spread(b,dim=1,ncopies=size(a))
END FUNCTION FFE_outerprod  
 
!**********************************************************************!
!
SUBROUTINE FFE_MSG(msg)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: msg
  IF (FFE_SHOW_MSG/=0) THEN
    write(unit=FFE_SOUT,fmt='(A)') "FFE > "//trim(msg)
  END IF
END SUBROUTINE FFE_MSG

!**********************************************************************!
!
SUBROUTINE FFE_DBG(msg,ndbg)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: msg
  INTEGER*4, INTENT(IN) :: ndbg
  IF (FFE_SHOW_DBG>=ndbg) THEN
    write(unit=FFE_SOUT,fmt='(A)') "FFE-DBG > "//trim(msg)
  END IF
END SUBROUTINE FFE_DBG

!**********************************************************************!
!
SUBROUTINE FFE_ERR(nerr)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: nerr
  CHARACTER(len=256) :: emsg
  FFE_ERR_CUR = nerr
  IF (FFE_SHOW_ERR/=0) THEN
    emsg = "Unknown."
    if (nerr>0 .and. nerr<=FFE_ERR_MAX) emsg = FFE_ERR_MSG(nerr)
    write(unit=FFE_SOUT,fmt='(A,I4.4,A)') &
         & "FFE-ERR(",nerr,") > "//trim(emsg)
  END IF
  IF (FFE_ERR_STOP/=0) STOP 'Terminated by module fitfeprm.'
END SUBROUTINE FFE_ERR
 
!**********************************************************************!
!
SUBROUTINE FFE_SET_S(n,x)
  implicit none
  integer*4, intent(in) :: n
  real*4, intent(in) :: x(n)
  integer*4 :: n1
  n1 = 0
  FFE_ERR_ALLOC = 0
  if (allocated(FFE_S)) n1 = size(FFE_S)
  if (n>0) then
    if (n/=n1) then
      if (n1>0) deallocate(FFE_S,stat=FFE_ERR_ALLOC)
      if (FFE_ERR_ALLOC/=0) goto 102
      allocate(FFE_S(n),stat=FFE_ERR_ALLOC)
      if (FFE_ERR_ALLOC/=0) goto 101
    end if
    FFE_S(1:n) = x(1:n)
    return
  else
    deallocate(FFE_S,stat=FFE_ERR_ALLOC)
    if (FFE_ERR_ALLOC/=0) goto 102
  end if
  return
101 call FFE_ERR(1)
    return
102 call FFE_ERR(2)
    return
END SUBROUTINE FFE_SET_S

!**********************************************************************!
!
SUBROUTINE FFE_SET_F(n,y)
  implicit none
  integer*4, intent(in) :: n
  real*4, intent(in) :: y(n)
  integer*4 :: n1
  n1 = 0
  FFE_ERR_ALLOC = 0
  if (allocated(FFE_F)) n1 = size(FFE_F)
  if (n>0) then
    if (n/=n1) then
      if (n1>0) deallocate(FFE_F,stat=FFE_ERR_ALLOC)
      if (FFE_ERR_ALLOC/=0) return
      allocate(FFE_F(n),stat=FFE_ERR_ALLOC)
      if (FFE_ERR_ALLOC/=0) return
    end if
    FFE_F(1:n) = y(1:n)
    return
  else
    deallocate(FFE_F,stat=FFE_ERR_ALLOC)
    if (FFE_ERR_ALLOC/=0) return
  end if
  return
101 call FFE_ERR(1)
    return
102 call FFE_ERR(2)
    return
END SUBROUTINE FFE_SET_F

!**********************************************************************!
!
SUBROUTINE FFE_SET_W(n,y)
  implicit none
  integer*4, intent(in) :: n
  real*4, intent(in) :: y(n)
  integer*4 :: n1
  n1 = 0
  FFE_ERR_ALLOC = 0
  if (allocated(FFE_W)) n1 = size(FFE_W)
  if (n>0) then
    if (n/=n1) then
      if (n1>0) deallocate(FFE_W,stat=FFE_ERR_ALLOC)
      if (FFE_ERR_ALLOC/=0) return
      allocate(FFE_W(n),stat=FFE_ERR_ALLOC)
      if (FFE_ERR_ALLOC/=0) return
    end if
    FFE_W(1:n) = y(1:n)
    return
  else
    deallocate(FFE_W,stat=FFE_ERR_ALLOC)
    if (FFE_ERR_ALLOC/=0) return
  end if
  return
101 call FFE_ERR(1)
    return
102 call FFE_ERR(2)
    return
END SUBROUTINE FFE_SET_W

!**********************************************************************!
!
SUBROUTINE FFE_SET_PRNG(prm_idx, plow, phigh)
  implicit none
  integer*4, intent(in) :: prm_idx
  real*4, intent(in) :: plow, phigh
  IF (prm_idx<0 .OR. prm_idx>FFE_PRM_MAX) goto 109
  FFE_PRNG(1,prm_idx) = plow
  FFE_PRNG(2,prm_idx) = phigh
  return
109 call FFE_ERR(9)
    return
END SUBROUTINE FFE_SET_PRNG



!**********************************************************************!
!
! FUNCTION FFE_FCHI
!
! Calculates and returns chi-square w.r.t. the
! list of parameters p(:) of the function FFE_FUNCPRM
!
! REMARKS:
! - Will modify FFE_DYDP.
!
FUNCTION FFE_FCHI(p)
!
  IMPLICIT NONE
!
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(INOUT) :: p
  REAL*4 :: FFE_FCHI
!
  INTEGER*4 :: i, n, nx, ny, np, nw
  REAL*4 :: vy, ctmp, diff, wgt
!
  FFE_FCHI = 0.0
  nx = 0
  ny = 0
  nw = 0
  np = FFE_PRM_MAX
  if (allocated(FFE_S)) nx = size(FFE_S)
  if (allocated(FFE_F)) ny = size(FFE_F)
  if (allocated(FFE_W)) nw = size(FFE_W)
  n = min( nx, ny )
  if (FFE_USE_WEIGHTS/=0) n = min(n,nw)
  if (n>0) then
    ctmp = 0.0
    if (FFE_USE_WEIGHTS==0) then
      ctmp = ( FFE_ZRED - &
             & (p(1)*p(1)+p(5)*p(5)+p(9)*p(9))/FFE_PREFAC)**2.0
      do i=1, n
        call FFE_FUNCPRM(FFE_S(i), p(1:np), vy, FFE_DFDP(1:np), np)
        diff = FFE_F(i) - vy
        ctmp = ctmp + diff*diff
      end do
    else
      ctmp = FFE_W_EXTRA*(FFE_ZRED - &
              & (p(1)*p(1)+p(5)*p(5)+p(9)*p(9))/FFE_PREFAC)**2.0
      do i=1, n
        call FFE_FUNCPRM(FFE_S(i), p(1:np), vy, FFE_DFDP(1:np), np)
        diff = FFE_F(i) - vy
        wgt = FFE_W(i)
        ctmp = ctmp + wgt*diff*diff
      end do
    end if
    FFE_FCHI = ctmp
  else
    goto 104
  end if
  return
104 FFE_ERR_CUR = 4
    return
END FUNCTION FFE_FCHI


!**********************************************************************!
!
! FUNCTION FFE_DFCHI
!
! Calculates and returns the gradient of chi-square w.r.t. the
! list of parameters p(:) of the function FFE_FUNCPRM
!
! REMARKS:
! - Will modify FFE_DFDP and FFE_PTMP.
!
FUNCTION FFE_DFCHI(p)
!
  IMPLICIT NONE
!
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(INOUT) :: p
  REAL*4, DIMENSION(FFE_PRM_MAX) :: FFE_DFCHI
!
  INTEGER*4 :: i, j, n, nx, ny, np, nw
  REAL*4 :: vy, ctmp, diff, wgt, econ
!
  FFE_DFCHI = 0.0
  nx = 0
  ny = 0
  nw = 0
  np = FFE_PRM_MAX
  if (allocated(FFE_S)) nx = size(FFE_S)
  if (allocated(FFE_F)) ny = size(FFE_F)
  if (allocated(FFE_W)) nw = size(FFE_W)
  n = min( nx, ny )
  if (FFE_USE_WEIGHTS/=0) n = min(n,nw)
  if (n>0) then
    FFE_PTMP = 0.0
    econ = FFE_ZRED - (p(1)*p(1)+p(5)*p(5)+p(9)*p(9))/FFE_PREFAC
    if (FFE_USE_WEIGHTS==0) then
      FFE_PTMP(1) = 2.0*p(1)/FFE_PREFAC*econ
      FFE_PTMP(5) = 2.0*p(5)/FFE_PREFAC*econ
      FFE_PTMP(9) = 2.0*p(9)/FFE_PREFAC*econ
      do i=1, n
        call FFE_FUNCPRM(FFE_S(i), p(1:np), vy, FFE_DFDP(1:np), np)
        diff = FFE_F(i) - vy
        do j=1, np
          FFE_PTMP(j) = FFE_PTMP(j) - 2.0*diff*FFE_DFDP(j)
        end do
      end do
    else
      FFE_PTMP(1) = FFE_W_EXTRA*2.0*p(1)/FFE_PREFAC*econ
      FFE_PTMP(5) = FFE_W_EXTRA*2.0*p(5)/FFE_PREFAC*econ
      FFE_PTMP(9) = FFE_W_EXTRA*2.0*p(9)/FFE_PREFAC*econ
      do i=1, n
        call FFE_FUNCPRM(FFE_S(i), p(1:np), vy, FFE_DFDP(1:np), np)
        diff = FFE_F(i) - vy
        wgt = FFE_W(i)
        do j=1, np
          FFE_PTMP(j) = FFE_PTMP(j) - 2.0*diff*FFE_DFDP(j)*wgt
        end do
      end do
    end if
    FFE_DFCHI(1:np) = FFE_PTMP(1:np)
  else
    goto 104
  end if
  return
104 call FFE_ERR(4)
    return
END FUNCTION FFE_DFCHI



!**********************************************************************!
!
! SUBROUTINE FFE_DFP
!
! Given a starting point p that is a vector of length N,
! the Broyden-Fletcher-Goldfarb-Shanno variant
! of Davidon-Fletcher-Powell minimization is performed
! on the implemented chi-square norm.
! The convergence requirement on zeroing the gradient is input as gtol.
! Returned quantities are
!   p (the location of the minimum),
!   iter (the number of iterations that were performed), and
!   fret (the minimum value of the function).
! The routine lnsrch is called to perform approximate line
! minimizations.
! Parameters:
!   ITMAX is the maximum allowed number of iterations;
!   STPMX is the scaled maximum step length allowed in line searches;
!   EPS is the machine precision;
!   TOLX is the convergence criterion on x values.
!
SUBROUTINE FFE_DFP(p,gtol,iter,fret)
!
  IMPLICIT NONE
!
  INTEGER*4, INTENT(OUT) :: iter
  REAL*4, INTENT(IN) :: gtol
  REAL*4, INTENT(OUT) :: fret
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(INOUT) :: p
!
  INTEGER*4 :: ITMAX
  REAL*4, PARAMETER :: STPMX=100.0,EPS=epsilon(p),TOLX=4.0*EPS
!
  INTEGER*4 :: its, i, np
  LOGICAL :: check
  REAL*4 :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
  REAL*4, DIMENSION(FFE_PRM_MAX) :: dg,g,hdg,pnew,xi
  REAL*4, DIMENSION(FFE_PRM_MAX,FFE_PRM_MAX) :: hessin
!
  ITMAX = FFE_ITMAX
  np = FFE_PRM_MAX
  fp=FFE_FCHI(p) ! Calculate starting function value
  g=FFE_DFCHI(p) ! and gradient
  hessin = 0.0        ! Init inverse Hessian to the unit matrix.
  do i=1, np
    hessin(i,i) = 1.0
  end do
  xi=-g ! Initial line direction.
  stpmax=STPMX*max(sqrt(dot_product(p(:),p(:))),real(np))
  do its=1, ITMAX ! Main loop over the iterations.
    iter=its
    ! The new function evaluation occurs in lnsrch;
    ! save the function value in fp for the next line search.
    ! It is usually safe to ignore the value of check.
    call FFE_LNSRCH(p,fp,g,xi,pnew,fret,stpmax,check)
    fp=fret
    xi=pnew-p ! Update the line direction, and the current point.
    p=pnew 
    ! Test for convergence on Dx.
    if (maxval(abs(xi)/max(abs(p),1.0)) < TOLX) RETURN
    dg=g ! Save the old gradient,
    g=FFE_DFCHI(p) ! and get the new gradient.
    den=max(fret,1.0)
    ! Test for convergence on zero gradient.
    if (maxval(abs(g)*max(abs(p),1.0)/den) < gtol) RETURN
    if (FFE_USE_RANGE/=0) then ! Test for range violation
      do i=1, np
        IF (p(i)>FFE_PRNG(2,i) .or. p(i)<FFE_PRNG(1,i)) GOTO 110
      end do
    end if
    dg=g-dg ! Compute difference of gradients,
    hdg=matmul(hessin,dg) ! and difference times current matrix.
    fac=dot_product(dg,xi) ! Calculate dot products for the denominators.
    fae=dot_product(dg,hdg)
    sumdg=dot_product(dg,dg)
    sumxi=dot_product(xi,xi)
    if (fac > sqrt(EPS*sumdg*sumxi)) then
      ! Skip update if fac not sufficiently positive.
      fac=1.0/fac
      fad=1.0/fae
      dg=fac*xi-fad*hdg ! Vector that makes BFGS different from DFP.
      hessin=hessin+fac*FFE_outerprod(xi,xi)-& !  The BFGS updating formula.
             & fad*FFE_outerprod(hdg,hdg)+fae*FFE_outerprod(dg,dg)
    end if
    xi=-matmul(hessin,g) ! Now calculate the next direction to go,
  end do ! and go back for another iteration.
  !
108 call FFE_ERR(8)
    RETURN
110 call FFE_ERR(10)
    RETURN
  !
END SUBROUTINE FFE_DFP



!**********************************************************************!
!
! SUBROUTINE FFE_LNSRCH
!
! Given an N-dimensional point xold, the value of the chi-square and
! gradient there, fold and g, and a direction p, finds a new point
! x along the direction p from xold where the function func has
! decreased "sufficiently."
! xold, g, p, and x are all arrays of length N.
! The new chi-square value is returned in f.
!   stpmax is an input quantity that limits the length of the steps
!          so that you do not try to evaluate the function in regions
!          where it is undefined or subject to overflow.
!   p is usually the Newton direction.
! The output quantity check is false on a normal exit.
! It is true when x is too close to xold.
! In a minimization algorithm, this usually signals convergence and
! can be ignored. However, in a zero-finding algorithm the calling
! program should check whether the convergence is spurious.
! Parameters:
!   ALF ensures sufficient decrease in function value;
!   TOLX is the convergence criterion on Dx.
!
SUBROUTINE FFE_LNSRCH(xold,fold,g,p,xnew,f,stpmax,check)
!
  IMPLICIT NONE
!
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(IN) :: xold,g
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(INOUT) :: p
  REAL*4, INTENT(IN) :: fold,stpmax
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(OUT) :: xnew
  REAL*4, INTENT(OUT) :: f
  LOGICAL, INTENT(OUT) :: check
!
  REAL*4, PARAMETER :: ALF=1.0e-4,TOLX=epsilon(xnew)
!
  INTEGER*4 :: ndum
  REAL*4 :: av,alam,alam2,alamin,bv,disc,f2,pabs,rhs1,rhs2,slope,tmplam
!  
  ndum=FFE_PRM_MAX
  check=.false.
  pabs=sqrt(dot_product(p(:),p(:)))
  if (pabs > stpmax) p(:)=p(:)*stpmax/pabs ! Scale if attempted step is too big.
  slope=dot_product(g,p)
  if (slope >= 0.0) goto 107
  alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0)) ! Compute lambda_min.
  alam=1.0 ! Always try full Newton step first.
  do ! Start of iteration loop.
    xnew(:)=xold(:)+alam*p(:)
    f=FFE_FCHI(xnew)
    if (alam < alamin) then ! Convergence on Dx. For zero finding,
                            ! the calling program should
                            ! verify the convergence.
      xnew(:)=xold(:)
      check=.true.
      RETURN
    else if (f <= fold+ALF*alam*slope) then ! Sufficient function decrease.
      RETURN
    else ! Backtrack.
      if (alam == 1.0) then ! First time.
        tmplam=-slope/(2.0*(f-fold-slope))
      else ! Subsequent backtracks.
        rhs1=f-fold-alam*slope
        rhs2=f2-fold-alam2*slope
        av=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        bv=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
            & (alam-alam2)
        if (av == 0.0) then
          tmplam=-slope/(2.0*bv)
        else
          disc=bv*bv-3.0*av*slope
          if (disc < 0.0) then
            tmplam=0.5*alam
          else if (bv <= 0.0) then
            tmplam=(-bv+sqrt(disc))/(3.0*av)
          else
            tmplam=-slope/(bv+sqrt(disc))
          end if
        end if
        if (tmplam > 0.5*alam) tmplam=0.5*alam ! lambda <=0.5*lambda1.
      end if
    end if
    alam2=alam
    f2=f
    alam=max(tmplam,0.1*alam) ! lambda >= 0.1*lambda1.
  end do ! Try again.
  return
106 call FFE_ERR(6)
    return
107 call FFE_ERR(7)
    return
END SUBROUTINE FFE_LNSRCH




!**********************************************************************!
!
! Simulated Annealing Simplex - Minimization
!
! SUBROUTINE FFE_SAS
!
! PARAMETERS:
!   real*4 p(:) = (INPUT) starting position
!               = (OUTPUT) position of the found minimum
!   real*4 gtol = (INPUT) fractional convergence tolerance
!   integer*4 iter = (OUTPUT) number of annealing steps performed
!   integer*4 ntry = (INPUT) number of trials asked per annealing step
!                  = (OUTPUT) number of trials used per annealing step
!   real*4 tred = (INPUT) temperature reduction with annealing
!   real*4 fret = (OUTPUT) minimum value
!
SUBROUTINE FFE_SAS(p, gtol, iter, fret)
  !
  IMPLICIT NONE
  !
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(INOUT) :: p
  REAL*4, INTENT(IN) :: gtol
  INTEGER*4, INTENT(INOUT) :: iter
  REAL*4, INTENT(OUT) :: fret
  !
  INTERFACE
    FUNCTION UniRand()
      IMPLICIT NONE
      REAL*4 :: UniRand
    END FUNCTION UniRand
  END INTERFACE
  !
  INTEGER*4, PARAMETER :: MINTRYFAC = 10 ! min number of trials wrt to problem dimension
  REAL*4, PARAMETER :: TRED_DEFAULT = 0.1 ! default temperature reduction
  REAL*4, PARAMETER :: EPS = EPSILON(gtol)
  !
  INTEGER*4 :: ITMAX ! max. number of annealing steps
  INTEGER*4 :: i, j, k ! standard iterators
  INTEGER*4 :: it ! iteration count
  INTEGER*4 :: ittry ! iteration count of trials
  INTEGER*4 :: numtry ! used number of  trials
  INTEGER*4 :: np ! dimension of the problem
  INTEGER*4 :: nstat ! status catcher
  REAL*4 :: temptr, tempred ! annealing temperatur and reduction factor
  REAL*4, DIMENSION(FFE_PRM_MAX+1,FFE_PRM_MAX) :: simp ! current simplex
  REAL*4, DIMENSION(FFE_PRM_MAX+1) :: simy ! current simplex energy
  REAL*4, DIMENSION(FFE_PRM_MAX) :: g ! initial gradient
  REAL*4 :: pchg
  !
  ! Initialize the routine
  ITMAX = FFE_ITMAX
  i = 0
  j = 0
  it = 0
  ittry = 0
  temptr = 0.0
  tempred = 1.0 - FFE_SAS_TRED
  if (tempred>0.999.or.tempred<0.001) tempred = 1.0 - TRED_DEFAULT
  np = FFE_PRM_MAX
  simp = 0.0
  simy = 0.0
  fret = FFE_FCHI(p) ! initialize with the current status
  g = FFE_DFCHI(p) ! get initial gradients
  numtry = max( FFE_SAS_NTRY, np*MINTRYFAC ) ! Set the used number of trials
                                     ! high enough to tickle alot of
                                     ! the simplex points per annealing
                                     ! step.
  !
  ! Determine the initial annealing temperature as the maximum 
  ! energy observed for the additional np simplex points
  ! - set the first np points with one dimension randomized
  !   around the current position by a random factor (0.5,1.5)
  SIMPLEXINIT: DO i=1, np
    simp(i,1:np) = p(1:np)
    pchg = 1.0
    simp(i,i) = simp(i,i)*(1.0+pchg*(UniRand()-0.5)) ! randomize one coordinate
    IF (FFE_USE_RANGE/=0) THEN ! clip made change to search range
      simp(i,i) = max(min(simp(i,i),FFE_PRNG(2,i)),FFE_PRNG(1,i))
    END IF
    simy(i) = FFE_FCHI(simp(i,1:np)) ! update the energy
  END DO SIMPLEXINIT
  ! - copy the input coordinate to the last point of the simplex
  simp(np+1,1:np) = p(1:np)
  simy(np+1) = fret
  !
  ! set the start temperature to twice the current maximum
  temptr = maxval(simy(:))
  !
  ! We have now a maximum energy (temperature estimate) temptr
  ! and a start simplex simp with pre-calculated energies simy.
  ! Run the annealing until it converges
  ANNEAL: DO it=1, ITMAX
    iter = it ! store number of iterations
    j = numtry
    call FFE_AMEBSA(simp,simy,p,fret,gtol,j,temptr)
    IF (j>0) EXIT ANNEAL ! convergence, stop here, further improvement not forseen
    IF (FFE_USE_RANGE/=0) THEN ! check for range violation
      ! throw error if any simplex point is out-of-range
      DO i=1,np
        IF (maxval(simp(:,i))>FFE_PRNG(2,i)) GOTO 110
        IF (minval(simp(:,i))<FFE_PRNG(1,i)) GOTO 110
      END DO
    END IF
    temptr = tempred*temptr
  END DO ANNEAL
  !
  RETURN
  !
101 call FFE_ERR(1)
    return
102 call FFE_ERR(2)
    return
103 call FFE_ERR(3)
    return
104 call FFE_ERR(4)
    return
110 call FFE_ERR(10)
    return
END SUBROUTINE FFE_SAS




!**********************************************************************!
!
! Single Simulated Annealing Step including a modified simplex minimizer
! using the internal chi-square function based on the input reference
! function func
!
SUBROUTINE FFE_AMEBSA(p,y,pb,yb,ftol,iter,temptr)
  IMPLICIT NONE
  REAL*4, DIMENSION(FFE_PRM_MAX+1,FFE_PRM_MAX), INTENT(INOUT) :: p
  REAL*4, DIMENSION(FFE_PRM_MAX+1), INTENT(INOUT) :: y
  REAL*4, DIMENSION(FFE_PRM_MAX), INTENT(INOUT) :: pb
  REAL*4, INTENT(INOUT) :: yb
  REAL*4, INTENT(IN) :: ftol,temptr
  INTEGER*4, INTENT(INOUT) :: iter
  INTERFACE
    FUNCTION UniRand()
      IMPLICIT NONE
      REAL*4 :: UniRand
    END FUNCTION UniRand
  END INTERFACE
  INTEGER*4, PARAMETER :: NMAX=2000
! Minimization of the N-dimensional function func by simulated annealing
! combined with the downhill simplex method of Nelder and Mead.
! The (N+1) * N matrix p is input. Its N+1 rows are N-dimensional
! vectors that are the vertices of the starting simplex. Also input is
! the vector y of length N+1, whose components must be preinitialized
! to the values of func evaluated at the N+1 vertices (rows) of p;
! ftol, the fractional convergence tolerance to be achieved in the
! function value for an early return;
! iter, and temptr. The routine makes iter function evaluations at an
! annealing temperature temptr, then returns. You should then decrease
! temptr according to your annealing schedule, reset iter, and call the
! routine again (leaving other arguments unaltered between calls).
! If iter is returned with a positive value, then early convergence and
! return occurred. If you initialize yb to a very large value on the
! first call, then yb and pb (an array of length N) will subsequently
! return the best function value and point ever encountered (even if it
! is no longer a point in the simplex).
  INTEGER*4 :: ihi,ndim ! Global variables.
  REAL*4 :: yhi
  REAL*4, DIMENSION(FFE_PRM_MAX) :: psum
  call amebsa_private
  
  CONTAINS

SUBROUTINE amebsa_private
  INTEGER*4 :: i,ilo,inhi
  REAL*4 :: rtol,ylo,ynhi,ysave,ytry
  REAL*4, DIMENSION(FFE_PRM_MAX+1) :: yt,harvest
  !ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amebsa')
  ndim = FFE_PRM_MAX
  psum(:)=sum(p(:,:),dim=1)
  do ! Iteration loop.
    do i=1, ndim+1
      harvest(i) = UniRand()
    end do
    yt(:)=y(:)-temptr*log(harvest)
! Whenever we 'look at' a vertex, it gets a random thermal fluctuation.
    ilo=FFE_iminloc(yt(:)) !  Determine which point is the highest (worst),
    ylo=yt(ilo)        ! next-highest, and lowest (best).
    ihi=FFE_imaxloc(yt(:))
    yhi=yt(ihi)
    yt(ihi)=ylo
    inhi=FFE_imaxloc(yt(:))
    ynhi=yt(inhi)
    rtol=2.0*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
! Compute the fractional range from highest to lowest and return if satisfactory.
    if (rtol < ftol .or. iter < 0) then ! If returning, put best point and value in
      call FFE_swapr(y(1),y(ilo))            ! slot 1.
      call FFE_swaprv(p(1,:),p(ilo,:))
      RETURN
    end if
! Begin a new iteration. First extrapolate by a factor -1 through the
! face of the simplex across from the high point, i.e., reflect the
! simplex from the high point.
    ytry=amotsa(-1.0)
    iter=iter-1
    if (ytry <= ylo) then ! Gives a result better than the best point, so
                          ! try an additional extrapolation by a factor of 2.
      ytry=amotsa(2.0)
      iter=iter-1
    else if (ytry >= ynhi) then ! The reflected point is worse than the secondhighest,
                                ! so look for an intermediate lower point,
                                ! i.e. do a one-dimensional contraction.
      ysave=yhi
      ytry=amotsa(0.5)
      iter=iter-1
      if (ytry >= ysave) then
! Can't seem to get rid of that high point. Better contract around the lowest
! (best) point.
        p(:,:)=0.5*(p(:,:)+spread(p(ilo,:),1,ndim+1))
        do i=1,ndim+1
          if (i /= ilo) y(i)=FFE_FCHI(p(i,:))
        end do
        iter=iter-ndim ! Keep track of function evaluations.
        psum(:)=sum(p(:,:),dim=1)
      end if
    end if
  end do
  return
  !
106 call FFE_ERR(6)
    return
END SUBROUTINE amebsa_private



FUNCTION amotsa(fac)
  IMPLICIT NONE
  REAL*4, INTENT(IN) :: fac
  REAL*4 :: amotsa
! Extrapolates by a factor fac through the face of the simplex across from the high point,
! tries it, and replaces the high point if the new point is better.
  REAL*4 :: fac1,fac2,yflu,ytry,harv
  REAL*4, DIMENSION(FFE_PRM_MAX) :: ptry
  fac1=(1.0-fac)/ndim
  fac2=fac1-fac
  ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
  ytry=FFE_FCHI(ptry)
  if (ytry <= yb) then ! Save the best-ever.
    pb(:)=ptry(:)
    yb=ytry
  end if
  harv = UniRand()
  yflu=ytry+temptr*log(harv) ! We added a thermal fluctuation to all the current
                             ! vertices, but we subtract it here, so
                             ! as to give the simplex a thermal Brownian
                             ! motion: It likes to accept any suggested change.
  if (yflu < yhi) then
    y(ihi)=ytry
    yhi=yflu
    psum(:)=psum(:)-p(ihi,:)+ptry(:)
    p(ihi,:)=ptry(:)
  end if
  amotsa=yflu
END FUNCTION amotsa

END SUBROUTINE FFE_AMEBSA


!**********************************************************************!
!
! SUBROUTINE FFE_FINDMIN
!
SUBROUTINE FFE_FINDMIN(p, iter, ctol, method, chisqr)
!
  IMPLICIT NONE
!
  REAL*4, INTENT(INOUT) :: p(FFE_PRM_MAX), ctol, chisqr
  INTEGER*4, INTENT(INOUT) :: iter
  INTEGER*4, INTENT(IN) :: method
!
  INTEGER*4 :: n, nx, ny, ny1, np
  REAL*4, DIMENSION(FFE_PRM_MAX) :: ptmp
!
  iter = 0
  nx = 0
  ny = 0
  ny1 = 0
  np = FFE_PRM_MAX
  ptmp(1:np) = p(1:np)
  if (method<1.or.method>FFE_MET_MAX) goto 105
  if (allocated(FFE_S)) nx = size(FFE_S)
  if (allocated(FFE_F)) ny = size(FFE_F)
  n = min( nx, ny )
  if (n>0) then
    iter = 0
    select case (method)
    case (1) ! BFGS minimization
      call FFE_DFP( ptmp(1:np), ctol, iter, chisqr)
    case (2) ! Simulated Annealing Simplex
      call FFE_SAS( ptmp(1:np), ctol, iter, chisqr)
    end select ! case (method)
  else
    goto 104
  end if
  p(1:np) = ptmp(1:np)
  return
103 call FFE_ERR(3)
    return
104 call FFE_ERR(4)
    return
105 call FFE_ERR(5)
    return
END SUBROUTINE FFE_FINDMIN



!**********************************************************************!
!
! subroutine FFE_FUNCPRM
!
! Calculation of function values and derivatives for a
! function of N Lorentzians and N Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,N [   ai^2 / ( bi^2 + x^2 + ceps)
!                    + ci^2 * Exp( - (di^2 + ceps) * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, b1, c1, d1, a2, b2, c2, d2, ..., aN, bN, cN, dN )
!
SUBROUTINE FFE_FUNCPRM(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given == 4*n
! -------------------------------------------------------------------- !

  IMPLICIT NONE
  
  real*4, parameter :: ceps = 1.0E-10

  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: n, i, i1, i2, i3, i4
  real*4 :: aa, ab, ac, ad, fexp, fearg, flor, flarg

! initialization
  y = 0.0
  dyda = 0.0
  n = int(na/4)

! function value calculation
  do i=1, n
    
    ! - iterators
    i1 = 4*(i-1) + 1
    i2 = 4*(i-1) + 2
    i3 = 4*(i-1) + 3
    i4 = 4*(i-1) + 4
    ! - current parameters
    aa = a(i1)
    ab = a(i2)
    ac = a(i3)
    ad = a(i4)
    ! - current lorentzian values
    flarg = x*x + ab*ab + ceps
    flor = 1.0 / flarg
    ! - current gaussian values
    fearg = x*x * (ceps+ad*ad)
    fexp = exp(-fearg)
    
    ! sum up functions
    y = y + aa*aa * flor + ac*ac * fexp
    
    ! sum up derivatives
    dyda(i1) = 2.0*aa*flor
    dyda(i2) = -2.0*aa*aa*ab*flor*flor
    dyda(i3) = 2.0*ac*fexp
    dyda(i4) = -2.0*ac*ac*ad*x*x*fexp
    
  end do

  return

END SUBROUTINE FFE_FUNCPRM


!**********************************************************************!
!
! subroutine FFE_FEPRM
!
! Calculation of function values and derivatives for a
! function of N Lorentzians and N Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,N [   ai / ( bi + x^2)
!                    + ci * Exp( - di * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, b1, c1, d1, a2, b2, c2, d2, ..., aN, bN, cN, dN )
!
SUBROUTINE FFE_FEPRM(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given == 4*n
! -------------------------------------------------------------------- !

  IMPLICIT NONE
  
  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: n, i, i1, i2, i3, i4
  real*4 :: aa, ab, ac, ad, fexp, fearg, flor, flarg

! initialization
  y = 0.0
  dyda = 0.0
  n = int(na/4)

! function value calculation
  do i=1, n
    
    ! - iterators
    i1 = 4*(i-1) + 1
    i2 = 4*(i-1) + 2
    i3 = 4*(i-1) + 3
    i4 = 4*(i-1) + 4
    ! - current parameters
    aa = a(i1)
    ab = a(i2)
    ac = a(i3)
    ad = a(i4)
    ! - current lorentzian values
    flarg = x*x + ab
    flor = 1.0 / flarg
    ! - current gaussian values
    fearg = x*x * ad
    fexp = exp(-fearg)
    
    ! sum up functions
    y = y + aa * flor + ac * fexp
    
    ! sum up derivatives
    dyda(i1) = flor
    dyda(i2) = -aa*flor*flor
    dyda(i3) = fexp
    dyda(i4) = -ac*x*x*fexp
    
  end do

  return

END SUBROUTINE FFE_FEPRM



!**********************************************************************!
!
! subroutine FFE_FEPRMY
!
! Calculation of function values for a
! function of N Lorentzians and N Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,N [   ai / ( bi + x^2)
!                    + ci * Exp( - di * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, b1, c1, d1, a2, b2, c2, d2, ..., aN, bN, cN, dN )
!
SUBROUTINE FFE_FEPRMY(x, a, y, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            na: integer*4: size of parameter array given == 4*n
! -------------------------------------------------------------------- !

  IMPLICIT NONE
  
  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: y
  integer*4 :: n, i, i1
  !real*4 :: fexp, fearg, flor, flarg
  real*4 :: x2, a1, b1, a2, b2

! initialization
  y = 0.0
  n = int(na/4)
  x2 = x*x

! function value calculation
  do i=1, n
    
    ! - iterator
    i1 = 4*(i-1)
    ! - parameters
    a1 = a(i1+1)
    b1 = a(i1+2)
    a2 = a(i1+3)
    b2 = a(i1+4)
    !! - current lorentzian values
    !flarg = x*x + a(i1+2)
    !flor = 1.0 / flarg
    ! - current gaussian values
    !fearg = x*x * a(i1+4)
    !fexp = exp(-fearg)
    
    ! sum up functions
    !y = y + a(i1+1) * flor + a(i1+3) * fexp
    y = y + a1 / (x2+b1) + a2 * exp( -b2*x2 )
    
  end do

  return

END SUBROUTINE FFE_FEPRMY


END MODULE fitfeprm




