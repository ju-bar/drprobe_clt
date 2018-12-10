!**********************************************************************!
!**********************************************************************!
!                                                                      !
!                        file   "BasicFuncs".f95                       !
!                                                                      !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                  !
!    Version  :  1.0.0, Mai 03, 2003                                   !
!                                                                      !
!                                                                      !
!**********************************************************************!


!**********************************************************************!
!                                                                      !
!   Purpose: implementing basic helping subroutines and functions      !
!            accessable from F90/95 AND F77                            !
!                                                                      !
!   Includes: none                                                     !
!                                                                      !
!**********************************************************************!
!                                                                      !
!   CONTAINS:                                                          !
!      1) system functions                                             !
!           SUBROUTINE getscalestr(rval,rrsval,scstr)                  !
!           SUBROUTINE gettimestring(stime)                            !
!           SUBROUTINE getclockms(ncount,soption,msecs)                !
!           INTEGER*4 FUNCTION IsSysBigEndian()                        !
!           INTEGER*2 FUNCTION SwapInt2(nin)                           !
!           INTEGER*4 FUNCTION SwapInt4(nin)                           !
!           SUBROUTINE CharTransfer(sarr,sstr,n,ndir)                  !
!      2) geometric functions                                          !
!           SUBROUTINE CrossProd(a,b,c)                                !
!           SUBROUTINE ScalProd(a,b,v,n)                               !
!           SUBROUTINE cart2aradial(x,dx,y,dy,a,r,dr,phi,dphi)         !
!           SUBROUTINE radial2acart(r,dr,phi,dphi,a,x,dx,y,dy)         !
!           LOGICAL FUNCTION nROIValidate(pt,roi,n,nwrap)              !
!           LOGICAL FUNCTION nInROI(pt,roi,n)                          !
!           INTEGER*4 FUNCTION get_nroi_area(roi,n)                    !
!           real*4 FUNCTION Val2Unity(rval,rv1,rv2,nMethod)            !
!           SUBROUTINE EllipTransform2D(rin,rout,n,ra1,rb1,rr1,ra2,rb2,rr2)!
!           SUBROUTINE EllipStretch2D(rin,rout,n,ra,rb,rang)           !
!           SUBROUTINE Rotate2D(rin,rout,n,rang)                       !
!           SUBROUTINE Rotate3D(rin,rout,n,rang,axis)                  !
!           SUBROUTINE DrawLine(rimg,nx,ny,px1,py1,px2,py2,rval)       !
!           SUBROUTINE DrawRect(rimg,nx,ny,px1,py1,px2,py2,bval,bpx, & !
!     &                          fval,fpx)                              !
!           SUBROUTINE DrawCircle(rimg,nx,ny,mx,my,rrad,rval)          !
!           SUBROUTINE DrawCircleR(rimg,nx,ny,rmx,rmy,rrad, &          !
!     &                             rval,nbflag,rfval,nfflag)           !
!           SUBROUTINE DrawLegend24b(bmDest,nPX,nPy,nSX,nSY,nCMethod,rLB,rUB,&
!     &                   nNOL,nDrawLabels,bfcLabel,nPMethod)
!           SUBROUTINE InvMat2x2(a11,a12,a21,a22,b11,b12,b21,b22)      !
!      3) file functions                                               !
!           SUBROUTINE raw_open(nunit,sname)                           !
!           SUBROUTINE raw_c_read(nunit,pos,c,n)                       !
!           SUBROUTINE raw_uread16(sfile, datablock, nsize, fswap, errorcode)
!           SUBROUTINE raw_uwrite16(sfile, datablock, nsize, fswap, errorcode)
!           SUBROUTINE raw_uread32(sfile, datablock, ndatasize, fswap, errorcode)
!           CorrectEndian16(ndata,nx,ny)
!           CorrectEndian32(ndata,nx,ny)
!      4) numerical functions
!           real*4 HT2WL(ht)                                        !
!           real*4 WL2HT(ht)                                        !
!           SUBROUTINE getMeanandSigma(dat,n,mean,sigma)  !
!           SUBROUTINE checkprimefactors(n,pf,npf,nresid) !
!           integer*4 FUNCTION nextprime(n,inc)                        !
!           real*4 FUNCTION gammp(a,x)                                 !
!           real*4 FUNCTION gammq(a,x)                                 !
!           SUBROUTINE gser(gamser,a,x,gln)                            !
!           SUBROUTINE gcf(gammcf,a,x,gln)                             !
!           real*4 FUNCTION gammln(xx)                                 !
!           real*4 FUNCTION rgauss2d(x,y,A0,A,x0,y0,w)                 !
!           integer*4 FUNCTION factorial(n)                            !
!           integer*8 FUNCTION dfactorial(n)                           !
!           integer*4 FUNCTION binomial(n,k)                           !
!           real*4 FUNCTION sigmoid(x,x0,dx)                           !
!      5) string functions                                             !
!           SUBROUTINE SetVarString(carray,string)                     !
!           SUBROUTINE GetVarString(string,carray)                     !
!                                                                      !
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
!************************* BASIC FUNCTIONS ****************************!
!**********************************************************************!
!**********************************************************************!






































!**********************************************************************!
!************************* SYSTEM FUNCTIONS ***************************!
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE getscalestr(rval,rrsval,scstr)
! function: returns scaled value and scaling dimension string
! -------------------------------------------------------------------- !
! parameter: rval : real*4 : value to scale [nm]
!            rrsval : real*4 : rescaling factor
!            scstr : character*1 : scaling string (f,p,n,u,m, ,k,M,G,T)
! -------------------------------------------------------------------- !


  implicit none

! ------------
! declaration
  real*4, intent(in) :: rval
  real*4, intent(out) :: rrsval
  character*1 :: scstr
  
  real*4 :: rabs 
! ------------


! ------------
! init
!  write(*,*) "getscalestr: INIT."
  scstr = ' '
  rrsval = 1.
  rabs = abs(rval)
!  write(*,*) "getscalestr: rabs:",rabs
  if (rabs==0.0) then
    rrsval = 1.e6
    scstr = 'n'
    return
  end if
! ------------


! ------------
  if (rabs<1.e-3) then
    rrsval = 1.e6
    scstr = 'f'
    return
  end if
  if (rabs<1.) then
    rrsval = 1.e3
    scstr = 'p'
    return
  end if
  if (rabs<1.e3) then
    rrsval = 1.
    scstr = 'n'
    return
  end if
  if (rabs<1.e6) then
    rrsval = 1.e-3
    scstr = 'u'
    return
  end if
  if (rabs<1.e9) then
    rrsval = 1.e-6
    scstr = 'm'
    return
  end if
  if (rabs<1.e12) then
    rrsval = 1.e-9
    scstr = ' '
    return
  end if
  if (rabs<1.e15) then
    rrsval = 1.e-12
    scstr = 'k'
    return
  end if
  if (rabs<1.e18) then
    rrsval = 1.e-15
    scstr = 'M'
    return
  end if
  if (rabs<1.e21) then
    rrsval = 1.e-18
    scstr = 'G'
    return
  end if
  if (rabs<1.e24) then
    rrsval = 1.e-21
    scstr = 'T'
    return
  end if

! ------------


! ------------
!  write(*,*) "getscalestr: EXIT."
  return

END SUBROUTINE getscalestr
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE gettimestring(stime)
! function: returns current time, format "YY/MM/DD hh:mm:ss"
! -------------------------------------------------------------------- !
! parameter: stime :: character(len=*) :: return reference
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character (len=*) :: stime

  character*8 :: sysdate
  character*10 :: systime
  character*5 :: syszone
  integer, dimension(8) :: values
! ------------



! ------------
! init
!  write(*,*) "gettimestring: INIT."
! ------------


! ------------
  call date_and_time(sysdate,systime,syszone,values)
  stime=sysdate(3:4)//"/"//sysdate(5:6)//"/"//sysdate(7:8)//" "&
     & //systime(1:2)//":"//systime(3:4)//":"//systime(5:6)
! ------------


! ------------
!  write(*,*) "gettimestring: EXIT.: ",trim(stime)
  return

END SUBROUTINE gettimestring
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE getclockms(ncount,soption,msecs)
! function: calculates milliseconds between calls
!           call with "START" option to start clock
!           call with "TAKE" option only requests current time differ.
!           call with "STOP" option to stop clock
! -------------------------------------------------------------------- !
! parameter:
!            ncount : integer*4 : is a counter that will recieve
!                                 system clock count on start-run
!                                 and transfer old count on stop-run
!            soption : character(*) : option selector
!                                     "START" -> start-run
!                                     "STOP" -> stop-run
!            msecs : real*4 : milliseconds (only valid after stop-run)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: ncount, nrate, nmax, ocount, start_clock, ndelta
  real*4, intent(out) :: msecs
  character(len=*), intent(in) :: soption
! ------------



! ------------
! init
!  write(*,*) "getclockms: INIT."
  msecs = 0.0
  ocount = 0
  start_clock = -1
  if ( trim(soption)=="START".or.&
     & trim(soption)=="Start".or.&
     & trim(soption)=="start") then
    start_clock = 1
  end if
  if ( trim(soption)=="TAKE".or.&
     & trim(soption)=="Take".or.&
     & trim(soption)=="take") then
    start_clock = 2
    ocount = ncount
  end if
  if ( trim(soption)=="STOP".or.&
     & trim(soption)=="Stop".or.&
     & trim(soption)=="stop") then
    start_clock = 0
    ocount = ncount
  end if
  if (start_clock==-1) then
    write(*,*) "getclockms: ERROR in parameter list (soption)."
    return
  end if
! ------------


! ------------
! request clock
  call system_clock(ncount, nrate, nmax)
  if (start_clock /= 1) then
    ndelta = ncount-ocount
    if (ndelta<0) then
      ndelta = ndelta + nmax
    end if
    msecs = 1000.0*real(ndelta)/real(nrate)
  end if
  if ( start_clock == 2 ) ncount = ocount
! ------------


! ------------
!  write(*,*) "getclockms: EXIT."
  return

END SUBROUTINE getclockms
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
INTEGER*4 FUNCTION IsSysBigEndian()
! function: determines data structure of current machine
!           returns '1' if Bigendian and '0' if LittleEndian
! -------------------------------------------------------------------- !

  implicit none

  INTEGER*2 :: nTest
  INTEGER*1 :: nB
  INTEGER*1 :: nByte(2)
!  character(len=2) :: stest

!  write(*,*) "IsSysBigEndian: INIT."

  nTest = 17


!  write(*,*) "IsSysBigEndian: ",nTest

  nByte=TRANSFER(nTest,nB,2)
!  stest=transfer(ntest,stest)
  
!  write(*,*) "IsSysBigEndian: ",nByte,ichar(stest(1:1)),ichar(stest(2:2))

  IsSysBigEndian = 1

  IF (nByte(1)==0) THEN
   IsSysBigEndian = 0
  END IF

!  write(*,*) "IsSysBigEndian: EXIT.", IsSysBigEndian

END FUNCTION IsSysBigEndian
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
INTEGER*2 FUNCTION SwapInt2(nin)
! function: returns byteswapped input of type integer*2
! -------------------------------------------------------------------- !
! parameter: nin: integer*2: nin
!            
! -------------------------------------------------------------------- !

  implicit none

  INTEGER*2,intent(in) :: nin
  INTEGER*2 :: dummy


  dummy = 0

!  do i=15,0,-1
!    if (mod(i+1,8)==0) write(*,'(A1,$)') "|"
!    if (BTEST(nin,i)) then
!      write(*,'(A1,$)') "1"
!    else
!      write(*,'(A1,$)') "0"
!    end if
!  end do
!  write(*,'(A1,I10)') "|",nin

  call MVBITS(nin,0,8,dummy,8)
  call MVBITS(nin,8,8,dummy,0)

!  do i=15,0,-1
!    if (mod(i+1,8)==0) write(*,'(A1,$)') "|"
!    if (BTEST(dummy,i)) then
!      write(*,'(A1,$)') "1"
!    else
!      write(*,'(A1,$)') "0"
!    end if
!  end do
!  write(*,'(A1,I10)') "|",dummy
  
  SwapInt2 = dummy

END FUNCTION SwapInt2
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
INTEGER*4 FUNCTION SwapInt4(nin)
! function: returns byteswapped input of type integer*4
! -------------------------------------------------------------------- !
! parameter: nin: integer*4: nin
!            
! -------------------------------------------------------------------- !

  implicit none

  INTEGER*4,intent(in) :: nin
  INTEGER*4 :: dummy

  dummy = 0
  
!  do i=31,0,-1
!    if (mod(i+1,8)==0) write(*,'(A1,$)') "|"
!    if (BTEST(nin,i)) then
!      write(*,'(A1,$)') "1"
!    else
!      write(*,'(A1,$)') "0"
!    end if
!  end do
!  write(*,'(A1,I10)') "|",nin

  
  call MVBITS(nin,0,8,dummy,24)
  call MVBITS(nin,8,8,dummy,16)
  call MVBITS(nin,16,8,dummy,8)
  call MVBITS(nin,24,8,dummy,0)

!  do i=31,0,-1
!    if (mod(i+1,8)==0) write(*,'(A1,$)') "|"
!    if (BTEST(dummy,i)) then
!      write(*,'(A1,$)') "1"
!    else
!      write(*,'(A1,$)') "0"
!    end if
!  end do
!  write(*,'(A1,I10)') "|",dummy
!  write(*,*)
  
  SwapInt4 = dummy

END FUNCTION SwapInt4
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE CharTransfer(sarr,sstr,n,ndir)
! function: Transfers chars from array structure to string structure
!           the number of the transfered chars n
!           ndir == 0: sarr -> sstr
!           ndir == 1: sstr -> sarr
! -------------------------------------------------------------------- !
! parameter: sarr(n) : character : char array
!            sstr : character(len=n) : char string
!            n : integer*4 : char string and array size
!            ndir : integer*4 : conversion direction flag
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n, ndir
  character(len=n), intent(inout) :: sstr
  character, intent(inout) :: sarr(1:n)
  
  integer*4 :: i
! ------------



! ------------
! init
!  write(*,*) " > CharTransfer: INIT."
  if ((n>0).and.(ndir==0)) then
    do i=1,n
      sstr(i:i) = sarr(i)
    end do
  end if
  if ((n>0).and.(ndir==1)) then
    do i=1,n
      sarr(i) = sstr(i:i)
    end do
  end if
! ------------



! ------------
!  write(*,*) " > CharTransfer: EXIT."
  return

END SUBROUTINE CharTransfer
!**********************************************************************!































































!**********************************************************************!
!************************* GEOMETRIC FUNCTIONS ************************!
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE InvMat2x2(a11,a12,a21,a22,b11,b12,b21,b22)
  implicit none
  real*4, intent(in) :: a11,a12,a21,a22
  real*4, intent(inout) :: b11,b12,b21,b22
  real*4 :: da
  da = a11*a22-a12*a21
  b11 = a22/da
  b12 = -a21/da
  b21 = -a12/da
  b22 = a11/da
  return
END SUBROUTINE InvMat2x2
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE InvMat2x2R8(a11,a12,a21,a22,b11,b12,b21,b22)
  implicit none
  real*8, intent(in) :: a11,a12,a21,a22
  real*8, intent(inout) :: b11,b12,b21,b22
  real*8 :: da
  da = a11*a22-a12*a21
  b11 = a22/da
  b12 = -a21/da
  b21 = -a12/da
  b22 = a11/da
  return
END SUBROUTINE InvMat2x2R8
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE CrossProd(a,b,c)
! function: calculates the cross product of 2 vectors in 3d-space
! -------------------------------------------------------------------- !
! parameter: a(3),b(3) : input vectors
!            c(3) : output vector = a x b
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: a(3),b(3)
  real*4, intent(out) :: c(3)
! ------------


! ------------
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
! ------------


! ------------
!  write(*,*) "CrossProd: EXIT."
  return

END SUBROUTINE CrossProd
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE ScalProd(a,b,v,n)
! function: calculates the scalar product of 2 vectors
! -------------------------------------------------------------------- !
! parameter: a(n),b(n) : input vectors
!            v : output result a*b
!            n : dimension of vectors
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: a(n),b(n)
  integer*4, intent(in) :: n
  real*4, intent(out) :: v
  integer*4 :: i
! ------------


! ------------
  v = 0.0
  do i=1,n
    v = v + a(i)*b(i)
  end do
! ------------


! ------------
!  write(*,*) "ScalProd: EXIT."
  return

END SUBROUTINE ScalProd
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE cart2aradial(x,dx,y,dy,a,r,dr,phi,dphi)
! function: converts catresian measuremnts + errorbars to radial system
! -------------------------------------------------------------------- !
! parameter: all real*4, names say everything
!            a : real*4 : for x=r*cos(a*phi), y=r*sin(a*phi)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: x,dx,y,dy,a
  real*4, intent(out) :: r,dr,phi,dphi
  real*4 :: rpii
! ------------



! ------------
! init
  rpii = atan(1.0)*4.0
  dphi = 0.0
  dr = 0.0
  r = 0.0
  phi = 0.0
!  write(*,*) "cart2aradial: INIT.",x,dx,y,dy,a
  r = sqrt(x*x+y*y)
  dr = sqrt(0.5*(dx*dx+dy*dy))
  if (a==0.0) return
  phi = atan2(y,x)/a
  if (r==0.0) then
    dphi = 2*rpii
    return
  end if
  dphi = dr/r/a
  dphi = MIN(dphi,rpii*2.0) 
! ------------


! ------------
!  write(*,*) "cart2aradial: EXIT.",r,dr,phi,dphi
  return

END SUBROUTINE cart2aradial
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE radial2acart(r,dr,phi,dphi,a,x,dx,y,dy)
! function: converts radial measueremnts + errorbars to cartesian system
! -------------------------------------------------------------------- !
! parameter: all real*4, names say everything
!            a : real*4 : for x=r*cos(a*phi), y=r*sin(a*phi)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(out) :: x,dx,y,dy
  real*4, intent(in) :: r,dr,phi,dphi,a
  real*4 :: si, co
! ------------



! ------------
! init
  x = 0.0
  y = 0.0
  dx = 0.0
  dy = 0.0
!  write(*,*) "radial2acart: INIT.",r,dr,phi,dphi,a
  si = SIN(a*phi)
  co = COS(a*phi)
  x = r*co
  y = r*si
  dx = sqrt(co*co*dr*dr+a*a*y*y*dphi*dphi)
  dy = sqrt(si*si*dr*dr+a*a*x*x*dphi*dphi)
! ------------



! ------------
!  write(*,*) "radial2acart: EXIT.",x,dx,y,dy
  return

END SUBROUTINE radial2acart
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
LOGICAL FUNCTION nROIValidate(pt,roi,n,nwrap)
! function: returns .true. if pt is a valid point in roi or could be
!           wrapped
! -------------------------------------------------------------------- !
! parameter: pt(n) : integer*4 : point of interest
!            roi(2,n) : integer*4 : roi
!            n : integer*4 : dimension
!            nwrap : integer*4 : wrap around flag (0=no, 1=yes)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n, nwrap, roi(2,n)
  integer*4, intent(inout) :: pt(n)
  logical :: lresult

  integer*4 :: i, pd, nl(n)
! ------------



! ------------
! init
!  write(*,*) " > nROIValidate: INIT."
  lresult = .false.
  if (n<1) then
    nROIValidate = lresult
    return
  end if
  if (nwrap==1) then
    lresult = .true.
    nl(:)=abs(roi(2,:)-roi(1,:))
  end if
! ------------


! ------------
  do i=1, n
    if (roi(1,i)>pt(i).or.pt(i)>roi(2,i)) then
      if (nwrap==0) then
        lresult = .false.
        pd = pt(i)
        pt(i) = MAX(roi(1,i),MIN(roi(2,i),pd))
      else
        pt(i) = roi(1,i) + MOD(pt(i)-roi(1,i),nl(i))
      end if
    end if
  end do
! ------------

! ------------
!  write(*,*) " > nROIValidate: EXIT."
  nROIValidate = lresult
  return

END FUNCTION nROIValidate
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
LOGICAL FUNCTION nInROI(pt,roi,n)
! function: returns .true. if pt is within the n dimensional ROI roi
!           integer parameters
! -------------------------------------------------------------------- !
! parameter: pt(n) : integer*4 : point in n dimensional integer array
!            roi(2,n) : integer*4 : region of interest
!            n : integer*4 : number of dimensions
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n, pt(n), roi(2,n)

  integer*4 :: i
  logical :: lResult
! ------------



! ------------
! init
!  write(*,*) "nInROI: INIT."
  lResult = .true.
  if (n<1) then
    nInROI = .false.
    return
  end if
! ------------


! ------------
  do i=1, n
    if (roi(1,i)>pt(i).or.pt(i)>roi(2,i)) lResult = .false.
  end do
! ------------


! ------------
!  write(*,*) "nInROI: EXIT."
  nInROI = lResult
  return

END FUNCTION nInROI
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
INTEGER*4 FUNCTION get_nroi_area(roi,n)
! function: returns area of the n dimensional ROI roi
!           integer parameters
! -------------------------------------------------------------------- !
! parameter: roi(2,n) : integer*4 : region of interest
!            n : integer*4 : number of dimensions
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n, roi(2,n)

  integer*4 :: i, nResult
! ------------



! ------------
! init
!  write(*,*) "get_nroi_area: INIT."
  nResult = 1
! ------------


! ------------
  do i=1, n
    nResult = nResult*(roi(2,i)-roi(1,i))
  end do
! ------------


! ------------
!  write(*,*) "get_nroi_area: EXIT."
  get_nroi_area = nResult
  return

END FUNCTION get_nroi_area
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
real*4 FUNCTION Val2Unity(rval,rv1,rv2,nMethod)
! function: transforms the value rval within the interval (rv1,rv2) to
!           to the interval (0,1) using given functional method nMethod
! -------------------------------------------------------------------- !
! parameter: rval : real*4 : input value
!            rv1 : real*4 : lower value limit
!            rv2 : real*4 : upper value limit
!            nMethod : integer*4 : conversion method selector
!                  0 : linear (min/max)
!                  1 : log (min/max)
!                  2 : sqrt (min/max)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: rval, rv1, rv2
  integer*4, intent(in) :: nMethod
  
  real*4 :: ruval
! ------------



! ------------
! init
!  write(*,*) "Val2Unity: INIT."
  if (rv1>=rv2) then
    write(*,*) "Val2Unity: ERROR invalid parameters (rv1,rv2)."
    stop
  end if
  if (nMethod<0.or.nMethod>2) then
    write(*,*) "Val2Unity: ERROR invalid parameter (nMethod)."
    stop
  end if

  if (rval<=rv1) then
    Val2Unity = 0.0
!    write(*,*) "Val2Unity: rv1 kick.",rval
    return
  end if
  if (rval>=rv2) then
    Val2Unity = 1.0
!    write(*,*) "Val2Unity: rv2 kick.",rval
    return
  end if
! ------------


! ------------
  select case (nMethod)

    case (0)
!     linear
      ruval = (rval-rv1)/(rv2-rv1)

    case (1)
!     log
      ruval = LOG(1.0+EXP(1.0)*(rval-rv1)/(rv2-rv1))

    case (2)
!     sqrt
      ruval = SQRT((rval-rv1)/(rv2-rv1))

  end select
! ------------


! ------------
!  write(*,*) "Val2Unity: EXIT.", rval,ruval
  Val2Unity = ruval
  return

END FUNCTION Val2Unity
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE EllipTransform2D(rin,rout,n,ra1,rb1,rr1,ra2,rb2,rr2)
! function: multiplies the modulus of n 2d vectors by an elliptical
!           stretch factor, ellipse is defined by ra (long half axis),
!           rb (short half axis), rang (angle between x-axis and long
!           half axis) [deg]
!           This Function transforms from an elliptical system 1 to
!           another elliptical system 2
! -------------------------------------------------------------------- !
! parameter: rin(2,n) : real*4 : input vectors
!            rout(2,n) : real*4 : output vectors
!            n : integer*4 : number of vectors
!            ra1 , rb1, rr1 : real*4 : ellipse of current system rin
!            ra2 , rb2, rr2 : real*4 : ellipse of current system rin
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n
  real*4, intent(in) :: rin(2,n), ra1, rb1, rr1, ra2, rb2, rr2
  real*4, intent(out) :: rout(2,n)

  integer*4 :: i
  real*8 :: deg2rad, dpx, dpy
  real*8 :: dcos1, dsin1, dra1, drb1, d2cos1, d2sin1
  real*8 :: dcos2, dsin2, dra2, drb2, d2cos2, d2sin2
  real*8 :: d2c12, d2s12
  real*8 :: Mxx, Myx, Mxy, Myy
! ------------



! ------------
! init
!  write(*,*) " > EllipTransform2D: INIT."
  rout(:,:) = 0.0
  dra1 = DBLE(ra1)
  drb1 = DBLE(rb1)
  dra2 = DBLE(ra2)
  drb2 = DBLE(rb2)
  deg2rad = DATAN(1.0d+00)/45.0
  dcos1 = DCOS(DBLE(rr1)*deg2rad)
  d2cos1 = DCOS(DBLE(2.0*rr1)*deg2rad)
  dsin1 = DSIN(DBLE(rr1)*deg2rad)
  d2sin1 = DSIN(DBLE(2.0*rr1)*deg2rad)
  
  dcos2 = DCOS(DBLE(rr2)*deg2rad)
  d2cos2 = DCOS(DBLE(2.0*rr2)*deg2rad)
  dsin2 = DSIN(DBLE(rr2)*deg2rad)
  d2sin2 = DSIN(DBLE(2.0*rr2)*deg2rad)

  d2c12 = DCOS(DBLE(2.0*(rr1-rr2))*deg2rad)
  d2s12 = DSIN(DBLE(2.0*(rr1-rr2))*deg2rad)
  
  Mxx = -( (dra1 - drb1)*(dra2 + drb2)*d2cos1 &
     & + (dra1 - drb1)*(dra2 - drb2)*d2c12 &
     & - (dra1 + drb1)*(dra2 + drb2 + (dra2 - drb2)*d2cos2)) &
     & / (4.d+00*dra1*drb1)
  Myx = ( -( (dra1 - drb1)*(dra2 + drb2)*d2sin1 ) &
     & + (dra2 - drb2)*((-dra1 + drb1)*d2s12 + (dra1 + drb1)*d2sin2)) &
     & / (4.d+00*dra1*drb1)
  Mxy = ( -( (dra1 - drb1)*(dra2 + drb2)*d2sin1 ) &
     & + (dra2 - drb2)*(( dra1 - drb1)*d2s12 + (dra1 + drb1)*d2sin2)) &
     & / (4.d+00*dra1*drb1)
  Myy = -( (drb1 - dra1)*(dra2 + drb2)*d2cos1 &
     & + (dra1 - drb1)*(dra2 - drb2)*d2c12 &
     & + (dra1 + drb1)*(-dra2 - drb2 + (dra2 - drb2)*d2cos2)) &
     & / (4.d+00*dra1*drb1)
! ------------


! ------------
! transform
  do i=1, n

    dpx = DBLE(rin(1,i))
    dpy = DBLE(rin(2,i))

    if (dpx==0.0.and.dpy==0.0) cycle

    rout(1,i) = real(dpx*Mxx+dpy*Myx)
    rout(2,i) = real(dpx*Mxy+dpy*Myy)

  end do
! ------------


! ------------
!  write(*,*) " > EllipTransform2D: EXIT."
  return

END SUBROUTINE EllipTransform2D
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE EllipStretch2D(rin,rout,n,ra,rb,rang)
! function: multiplies the modulus of n 2d vectors by an elliptical
!           stretch factor, ellipse is defined by ra (long half axis),
!           rb (short half axis), rang (angle between x-axis and long
!           half axis) [deg]
! -------------------------------------------------------------------- !
! parameter: rin(2,n) : real*4 : input vectors
!            rout(2,n) : real*4 : output vectors
!            n : integer*4 : number of vectors
!            ra , rb : real*4 : ellipse axis of multipliers
!            rang : real*4 : ellipse rotation angle [deg]
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n
  real*4, intent(in) :: rin(2,n), ra, rb, rang
  real*4, intent(out) :: rout(2,n)
! ------------



! ------------
! init
!  write(*,*) " > EllipStretch2D: INIT."
  call EllipTransform2D(rin,rout,n,ra,rb,rang,1.0,1.0,0.0)
! ------------


! ------------
!  write(*,*) " > EllipStretch2D: EXIT."
  return

END SUBROUTINE EllipStretch2D
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE Rotate2D(rin,rout,n,rang)
! function: rotates a set of vectors in 2d space
! -------------------------------------------------------------------- !
! parameter: rin(2,n) : real*4 : input vectors
!            rout(2,n) : real*4 : output vectors
!            n : integer*4 : number of vectors
!            rang : real*4 : rotation angle [deg]
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, intent(in) :: n
  real*4, intent(in) :: rin(2,n), rang
  real*4, intent(out) :: rout(2,n)

  integer*4 :: i
  real*8 :: deg2rad, dcosa, dsina, dpx, dpy
! ------------


! ------------
! INIT
!  write(*,*) " > Rotate2D: INIT."
  rout = 0.0
  deg2rad = DATAN(1.d+00)/45.0
  dcosa = DCOS(dble(rang)*deg2rad)
  dsina = DSIN(dble(rang)*deg2rad)
! ------------


! ------------
  do i=1,n
  
    dpx = dble(rin(1,i))
    dpy = dble(rin(2,i))
    rout(1,i) = real(dcosa*dpx-dsina*dpy)
    rout(2,i) = real(dsina*dpx+dcosa*dpy)

!    write(*,*) i
!    write(*,*) sqrt(real(dpx*dpx+dpy*dpy)),datan2(dpy,dpx)/deg2rad
!    write(*,*) sqrt(real(rout(1,i)*rout(1,i)+rout(2,i)*rout(2,i))), &
!     &         atan2(rout(2,i),rout(1,i))/deg2rad
  end do
! ------------


! ------------
! write(*,*) " > Rotate2D: EXIT."
  return

END SUBROUTINE Rotate2D
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE Rotate2DR8(rin,rout,n,rang)
! function: rotates a set of vectors in 2d space
!           real*8 version
! -------------------------------------------------------------------- !
! parameter: rin(2,n) : real*8 : input vectors
!            rout(2,n) : real*8 : output vectors
!            n : integer*4 : number of vectors
!            rang : real*8 : rotation angle [deg]
! -------------------------------------------------------------------- !

  implicit none

  integer*4, intent(in) :: n
  real*8, intent(in) :: rin(2,n), rang
  real*8, intent(out) :: rout(2,n)
  integer*4 :: i
  real*8 :: deg2rad, dcosa, dsina, dpx, dpy

  rout = 0.0d+0
  deg2rad = DATAN(1.d+00)/45.0d+0
  dcosa = DCOS(rang*deg2rad)
  dsina = DSIN(rang*deg2rad)
  do i=1,n
    dpx = rin(1,i)
    dpy = rin(2,i)
    rout(1,i) = dcosa*dpx-dsina*dpy
    rout(2,i) = dsina*dpx+dcosa*dpy
  end do
  return

END SUBROUTINE Rotate2DR8
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE Rotate3D(rin,rout,n,rang,axis)
! function: rotates a set of vectors in 3d space
! -------------------------------------------------------------------- !
! parameter: rin(3,n) : real*4 : input vectors
!            rout(3,n) : real*4 : output vectors
!            n : integer*4 : number of vectors
!            rang : real*4 : rotation angle around axis [deg]
!            axis : integer*4 : rotation axis (1,2,3)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, intent(in) :: n, axis
  real*4, intent(in) :: rin(3,n), rang
  real*4, intent(out) :: rout(3,n)

  integer*4 :: i, j, k, naxis
  real*8 :: deg2rad, dcosa, dsina, dp(3)
  real*8 :: rmatrix(3,3)
! ------------


! ------------
! INIT
!  write(*,*) " > Rotate3D: INIT."
  rout = 0.0
  deg2rad = DATAN(1.d+00)/45.0
  dcosa = DCOS(dble(rang)*deg2rad)
  dsina = DSIN(dble(rang)*deg2rad)
  naxis = axis
  if (axis<1) naxis = 1
  if (axis>3) naxis = 3
! ------------


! ------------
! select axis to define rotation matrix
  rmatrix(:,:) = 0.0D+00
  select case (naxis)
    case (1)
      rmatrix(1,1) = 1.0D+00
      rmatrix(2,2) = dcosa
      rmatrix(3,2) = dsina
      rmatrix(2,3) = -dsina
      rmatrix(3,3) = dcosa
    case (2)
      rmatrix(1,1) = dcosa
      rmatrix(3,1) = -dsina
      rmatrix(2,2) = 1.0D+00
      rmatrix(1,3) = dsina
      rmatrix(3,3) = dcosa
    case (3)
      rmatrix(1,1) = dcosa
      rmatrix(2,1) = dsina
      rmatrix(1,2) = -dsina
      rmatrix(2,2) = dcosa
      rmatrix(3,3) = 1.0D+00
  end select
! ------------


! ------------
  do i=1,n
      
    dp(1) = dble(rin(1,i))
    dp(2) = dble(rin(2,i))
    dp(3) = dble(rin(3,i))
    
    do k=1,3
      rout(k,i) = 0.0
      do j=1,3
        rout(k,i) = rout(k,i) + real(rmatrix(j,k)*dp(j))
      end do
    end do
    
  end do
! ------------


! ------------
! write(*,*) " > Rotate3D: EXIT."
  return

END SUBROUTINE Rotate3D
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE DrawLine(rimg,nx,ny,px1,py1,px2,py2,rval)
! function: sets the pixels in rimg on coordinates of a line with
!           to value rval
! -------------------------------------------------------------------- !
! parameter: rimg(nx,ny) : real*4 : data array (image memory)
!            nx, ny : integer*4 : data array dimensions
!            px1, py1 : real*4 : start point of line
!            px2, py2 : real*4 : stop point of line
!            rval : real*4 : value to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: nx, ny
  real*4, intent(in) :: rval, px1, py1, px2, py2
  real*4, intent(inout) :: rimg(nx,ny)
  
  integer*4 :: l, npx, npy
  real*4 :: dx, dy, dlx, dly, dldx, dldy
! ------------



! ------------
! init
!  write(*,*) " > DrawLine: INIT."
  dx = px2 - px1
  dy = py2 - py1
! ------------



! ------------
! validity checks
  if (nx<=0 .or. ny<=0) return
! ------------



! ------------
! check point issue
  if (dx==0.0.and.dy==0.0) then
    npx = nint(px1)
    npy = nint(px2)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
    return ! all is done in this case
  end if
! ------------


! ------------
! select line path
  if (abs(dx)>=abs(dy)) then
    dlx = dx/abs(dx)
    dly = dy/abs(dx)
  else
    dly = dy/abs(dy)
    dlx = dx/abs(dy)
  end if
! ------------



! ------------
  l = 0
  dldx = 0.0
  dldy = 0.0
! loop along the line
  do

!   check the way already gone and exit if stop point is reached
    if (abs(dldx)>abs(dx).or.abs(dldy)>abs(dy)) exit
    
!   translate to pixel coordinates
    npx = nint(px1+dldx)
    npy = nint(py1+dldy)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval

!   advance one step
    dldx = dldx + dlx
    dldy = dldy + dly
    
  end do
! ------------



! ------------
!  write(*,*) " > DrawLine: EXIT."
  return

END SUBROUTINE DrawLine
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE DrawRect(rimg,nx,ny,px1,py1,px2,py2,bval,bpx,fval,fpx)
! function: sets the pixels in rimg on coordinates of a rect
!           to value rval
! -------------------------------------------------------------------- !
! parameter: rimg(nx,ny) : real*4 : data array (image memory)
!            nx, ny : integer*4 : data array dimensions
!            px1, py1 : real*4 : start point of rect
!            px2, py2 : real*4 : stop point of rect
!            bval, bpx : real*4 : border color value and thickness to set
!            fval, fpx : real*4 : fill color value and flag to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: nx, ny
  real*4, intent(in) :: px1, py1, px2, py2, bval, bpx, fval, fpx
  real*4, intent(inout) :: rimg(nx,ny)
  
  integer*4 :: i, j, i1, i2, j1, j2
! ------------



! ------------
! init
!  write(*,*) " > DrawRect: INIT."
! ------------



! ------------
! validity checks
  if (nx<=0 .or. ny<=0) return
! ------------


! ------------
! filling
  if (fpx>0.0) then
    j1 = NINT(py1)
    j2 = NINT(py2)
    i1 = NINT(px1)
    i2 = NINT(px2)
    do j=j1, j2
      if (j<1.or.j>ny) cycle ! border skip
      do i=i1, i2
        if (i<1.or.j>nx) cycle ! border skip
        rimg(i,j) = fval ! set color
      end do
    end do
  end if
! ------------


! ------------
  if (bpx>0.0) then
!   bottom border
    do j=NINT(py1-bpx*0.5),NINT(py1+bpx*0.4999)
      if (j<1.or.j>ny) cycle ! border skip
      do i=NINT(px1-bpx*0.5),NINT(px2+bpx*0.4999)
        if (i<1.or.j>nx) cycle ! border skip
          rimg(i,j) = bval ! set color
      end do
    end do
!   top border
    do j=NINT(py2-bpx*0.5),NINT(py2+bpx*0.4999)
      if (j<1.or.j>ny) cycle ! border skip
      do i=NINT(px1-bpx*0.5),NINT(px2+bpx*0.4999)
        if (i<1.or.j>nx) cycle ! border skip
          rimg(i,j) = bval ! set color
      end do
    end do
!   left border
    do j=NINT(py1-bpx*0.5),NINT(py2+bpx*0.4999)
      if (j<1.or.j>ny) cycle ! border skip
      do i=NINT(px1-bpx*0.5),NINT(px1+bpx*0.4999)
        if (i<1.or.j>nx) cycle ! border skip
          rimg(i,j) = bval ! set color
      end do
    end do
!   right border
    do j=NINT(py1-bpx*0.5),NINT(py2+bpx*0.4999)
      if (j<1.or.j>ny) cycle ! border skip
      do i=NINT(px2-bpx*0.5),NINT(px2+bpx*0.4999)
        if (i<1.or.j>nx) cycle ! border skip
          rimg(i,j) = bval ! set color
      end do
    end do
  end if
! ------------



! ------------
!  write(*,*) " > DrawRect: EXIT."
  return

END SUBROUTINE DrawRect
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE DrawCircle(rimg,nx,ny,mx,my,rrad,rval)
! function: sets the pixels in rimg on coordinates of a circle with
!           radius rrad and centre (mx,my) to value rval
! -------------------------------------------------------------------- !
! parameter: rimg(nx,ny) : real*4 : data array (image memory)
!            nx, ny : integer*4 : data array dimensions
!            mx, my : integer*4 : circle centre coordinates
!            rrad : real*4 : circle radius [pixel]
!            rval : real*4 : value to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: nx, ny, mx, my
  real*4, intent(in) :: rval, rrad
  real*4, intent(inout) :: rimg(nx,ny)
  
  integer*4 :: l, npx, npy
  real*4 :: cx, cy, rx, ry, rr2, ry2
! ------------



! ------------
! init
!  write(*,*) " > DrawCircle: INIT."
  cx = real(mx)
  cy = real(my)
! ------------


! ------------
! validity checks
  if (nx<=0 .or. ny<=0) return
! ------------



! ------------
  l = 0
  rr2 = rrad*rrad
  do
    ry = real(l)
    ry2 = ry*ry
    rx = sqrt(rr2-ry2)
    if (rx<ry) exit

!   q1: first quater
    npx = nint(cx + rx)
    npy = nint(cy + ry)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q2: x mirror q1
    npx = nint(cx + rx)
    npy = nint(cy - ry)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q3: y mirror q1
    npx = nint(cx - rx)
    npy = nint(cy + ry)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q4: y mirror q2
    npx = nint(cx - rx)
    npy = nint(cy - ry)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q5: reciprocal q1
    npx = nint(cx + ry)
    npy = nint(cy + rx)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q6: x mirror q5
    npx = nint(cx - ry)
    npy = nint(cy + rx)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q7: y mirror q5
    npx = nint(cx + ry)
    npy = nint(cy - rx)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
!   q8: y mirror q6
    npx = nint(cx - ry)
    npy = nint(cy - rx)
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval

    l = l + 1
  end do
! ------------



! ------------
!  write(*,*) " > DrawCircle: EXIT."
  return

END SUBROUTINE DrawCircle
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE DrawCircleR(rimg,nx,ny,rmx,rmy,rrad,rval,nbflag,rfval,nfflag)
! function: sets the pixels in rimg on coordinates of a circle with
!           radius rrad and centre (rmx,rmy) to value rval and fills
!           the circle with value rfval
! -------------------------------------------------------------------- !
! parameter: rimg(nx,ny) : real*4 : data array (image memory)
!            nx, ny : integer*4 : data array dimensions
!            rmx, rmy : real*4 : circle centre coordinates
!            rrad : real*4 : circle radius [pixel]
!            rval : real*4 : value to set at circle radius
!            nbflag : integer*4 : activates drawing at radius with rval
!                                 if /=0
!            rfval : real*4 : in-circle value
!            nfflag : integer*4 : activates drawing inside radius
!                                 with rfval if /=0
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: nx, ny, nbflag, nfflag
  real*4, intent(in) :: rval, rrad, rmx, rmy, rfval
  real*4, intent(inout) :: rimg(nx,ny)
  
  integer*4 :: l, npx, npy, npx2, npy2, nl1, nl2, nl
  real*4 :: cx, cy, rx, ry, rr2, ry2
! ------------



! ------------
! init
!  write(*,*) " > DrawCircleR: INIT."
  if (nbflag==0.and.nfflag==0) return
  cx = rmx
  cy = rmy
! ------------


! ------------
! validity checks
  if (nx<=0 .or. ny<=0) return
! ------------



! ------------
  l = 0
  rr2 = rrad*rrad
  do
    ry = real(l)
    ry2 = ry*ry
    rx = sqrt(rr2-ry2)
    if (rx<ry) exit

!   q1: first quater
    npx = nint(cx + rx)
    npy = nint(cy + ry)
!   q2: y mirror q1
    npx2 = nint(cx - rx)
    npy2 = nint(cy + ry)
!   horizontal line q2->q1)
    nl1 = MIN(nx,MAX(1,npx2))
    nl2 = MIN(nx,MAX(1,npx))
    nl = MIN(ny,MAX(1,npy))
!   fill
    if (nfflag/=0) then
    rimg(nl1:nl2,nl) = rfval
    end if
!   border draw
    if (nbflag/=0) then
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
    if (npx2>0.and.npx2<=nx.and.npy2>0.and.npy2<=ny) rimg(npx2,npy2) = rval
    end if
    
!   q3: x mirror q1
    npx = nint(cx + rx)
    npy = nint(cy - ry)
!   q4: y mirror q3
    npx2 = nint(cx - rx)
    npy2 = nint(cy - ry)
!   horizontal line q4->q3)
    nl1 = MIN(nx,MAX(1,npx2))
    nl2 = MIN(nx,MAX(1,npx))
    nl = MIN(ny,MAX(1,npy))
!   fill
    if (nfflag/=0) then
    rimg(nl1:nl2,nl) = rfval
    end if
!   border draw
    if (nbflag/=0) then
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
    if (npx2>0.and.npx2<=nx.and.npy2>0.and.npy2<=ny) rimg(npx2,npy2) = rval
    end if

!   q5: reciprocal q1
    npx = nint(cx + ry)
    npy = nint(cy + rx)
!   q6: y mirror q5
    npx2 = nint(cx + ry)
    npy2 = nint(cy - rx)
!   vertical line q6->q5)
    nl1 = MIN(ny,MAX(1,npy2))
    nl2 = MIN(ny,MAX(1,npy))
    nl = MIN(nx,MAX(1,npx))
!   fill
    if (nfflag/=0) then
    rimg(nl,nl1:nl2) = rfval
    end if
!   border draw
    if (nbflag/=0) then
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
    if (npx2>0.and.npx2<=nx.and.npy2>0.and.npy2<=ny) rimg(npx2,npy2) = rval
    end if

!   q7: x mirror q5
    npx = nint(cx - ry)
    npy = nint(cy + rx)
!   q8: y mirror q7
    npx2 = nint(cx - ry)
    npy2 = nint(cy - rx)
!   vertical line q8->q7)
    nl1 = MIN(ny,MAX(1,npy2))
    nl2 = MIN(ny,MAX(1,npy))
    nl = MIN(nx,MAX(1,npx))
!   fill
    if (nfflag/=0) then
    rimg(nl,nl1:nl2) = rfval
    end if
!   border draw
    if (nbflag/=0) then
    if (npx>0.and.npx<=nx.and.npy>0.and.npy<=ny) rimg(npx,npy) = rval
    if (npx2>0.and.npx2<=nx.and.npy2>0.and.npy2<=ny) rimg(npx2,npy2) = rval
    end if


    l = l + 1
  end do
! ------------



! ------------
!  write(*,*) " > DrawCircleR: EXIT."
  return

END SUBROUTINE DrawCircleR
!**********************************************************************!


!!!!**********************************************************************!
!!!!**********************************************************************!
!!!SUBROUTINE DrawLegend24b(bmDest,nPX,nPy,nSX,nSY,nCMethod,rLB,rUB,&
!!!     &                   nNOL,nDrawLabels,bfcLabel,nPMethod)
!!!
!!!
!!!
!!!! ------------------
!!!  USE BasicTypes
!!!  USE BitFont
!!!  USE RGBPalette
!!!! ------------------
!!!
!!!  implicit none
!!!
!!!
!!!! ------------------
!!!  INTEGER*1, DIMENSION(:,:,:) :: bmDest
!!!  INTEGER*4 :: nPX,nPY,nSX,nSY,nCMethod,nNOL,nPMethod,nDrawLabels
!!!  REAL*4 :: rLB,rUB,rVal,rLabVal
!!!  TYPE(bf_color) :: bfcLabel
!!!
!!!  INTEGER*1, DIMENSION(:,:,:), ALLOCATABLE :: bmLegend
!!!  INTEGER*4 :: nXdim, nYdim, alloc, i, npos, nStep
!!!  TYPE(bf_color) :: bfcOld
!!!  TYPE(t_rect) :: rcDest,rcLegend,rcEffective,rcE2
!!!  CHARACTER*12 :: sText
!!!  TYPE(rgb_col) :: color
!!!  real*4, external :: Val2Unity
!!!! ------------------
!!!
!!!
!!!
!!!
!!!
!!!! ------------------
!!!!    image dimension
!!!  nXdim = SIZE(bmDest,1)
!!!  nYdim = SIZE(bmDest,2)
!!!    
!!!!    parameterised rectangles
!!!  rcDest = t_rect(t_point(1,1),t_point(nXdim,nYdim))
!!!  rcLegend = t_rect(t_point(nPX,nPY),t_point(nPX+nSX-1,nPY+nSY-1))
!!!
!!!!    clip legend rect to image rect
!!!  rcEffective = InfltRect(rcDest,rcLegend)
!!!!   clip transform to legend rect
!!!  rcE2=t_rect(t_point(rcEffective%lb%x-nPX+1,rcEffective%lb%y-nPY+1),&
!!!     &        t_point(rcEffective%rt%x-nPX+1,rcEffective%rt%y-nPY+1))
!!!
!!!!      write(*,*)"DrawLegend24b: Dest: ([",rcDest%lb%x,",",rcDest%lb%y,"],[",rcDest%rt%x,",",rcDest%rt%y,"])"
!!!!      write(*,*)"DrawLegend24b: Legd: ([",rcLegend%lb%x,",",rcLegend%lb%y,"],[",rcLegend%rt%x,",",rcLegend%rt%y,"])"
!!!!      write(*,*)"DrawLegend24b: Effc: ([",rcEffective%lb%x,",",rcEffective%lb%y,"],[",rcEffective%rt%x,",",rcEffective%rt%y,"])"
!!!
!!!!    valid sizes required, exit if not          
!!!  IF ( (nXdim*nYdim<=0).OR.(nSX<7).OR.(nSY<5) &
!!!   & .OR.(.NOT.(IsRectPositive(rcEffective))) ) THEN
!!!    RETURN
!!!  END IF
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!  ALLOCATE(bmLegend(nSX,nSY,3),STAT=alloc)
!!!  IF (alloc /= 0) THEN
!!!    write(*,*) "DrawLegend24d: Allocation error."
!!!    STOP
!!!  END IF
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!! copy current image content to legend bitmap
!!!  bmLegend(rcE2%lb%x:rcE2%rt%x,rcE2%lb%y:rcE2%rt%y,:) &
!!!     & = bmDest(rcEffective%lb%x:rcEffective%rt%x,&
!!!     &          rcEffective%lb%y:rcEffective%rt%y,:)
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!! save current line/font color and set new
!!!  bfcOld = bfcCurrent
!!!  bfcCurrent = bfcLabel
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!! draw borders of legend
!!!!    left
!!!  bmLegend(1,:,1) = bfcCurrent%r
!!!  bmLegend(1,:,2) = bfcCurrent%g
!!!  bmLegend(1,:,3) = bfcCurrent%b
!!!!    right
!!!  bmLegend(nSX-4,:,1) = bfcCurrent%r
!!!  bmLegend(nSX-4,:,2) = bfcCurrent%g
!!!  bmLegend(nSX-4,:,3) = bfcCurrent%b
!!!!    bottom
!!!  bmLegend(:,1,1) = bfcCurrent%r
!!!  bmLegend(:,1,2) = bfcCurrent%g
!!!  bmLegend(:,1,3) = bfcCurrent%b
!!!!    top
!!!  bmLegend(:,nSY,1) = bfcCurrent%r
!!!  bmLegend(:,nSY,2) = bfcCurrent%g
!!!  bmLegend(:,nSY,3) = bfcCurrent%b
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!! draw+write labels
!!!  IF (nNOL>0) THEN
!!!
!!!!      min. number of labels is 2
!!!    IF (nNOL==1) THEN
!!!
!!!      nNOL=2
!!!
!!!    END IF
!!!
!!!    npos = 1
!!!
!!!! ------------------
!!!! create label text
!!!    DO i=1,nNOL
!!!
!!!      npos = INT((REAL((i-1)*(nSY-1))/REAL(nNOL-1))+0.5)+1
!!!
!!!      bmLegend(nSX-3:nSX,npos,1) = bfcCurrent%r
!!!      bmLegend(nSX-3:nSX,npos,2) = bfcCurrent%g
!!!      bmLegend(nSX-3:nSX,npos,3) = bfcCurrent%b
!!!
!!!      IF (nDrawLabels/=0) THEN
!!!
!!!        rLabVal = rLB + (rUB-rLB)*REAL(i-1)/REAL(nNOL-1)
!!!
!!!        IF ( (rLabVal<(-1.0E1)).OR. &
!!!          &  (rLabVal>(1.0E1)).OR. &
!!!          &  ((rLabVal>(-1.0E-1)).AND.(rLabVal<(1.0E-1))) ) THEN
!!!          WRITE (sText,FMT='(E9.3)') ( rLabVal )
!!!        ELSE
!!!          WRITE (sText,FMT='(F6.3)') ( rLabVal )
!!!        END IF
!!!
!!!        CALL WriteXY(nPX+nSX+2,nPY+npos-1,bmDest,sText)
!!!
!!!      END IF          
!!!
!!!    END DO
!!!! ------------------
!!!
!!!  END IF
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!! prepare legend color drawing
!!!  CALL SetCurrentPal(nPMethod)
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------    
!!!! draw legend color
!!!  DO i=1,nSY-2
!!!
!!!    rVal = rLB + (rUB-rLB)*REAL(i-1)/REAL(nSY-3)
!!!
!!!    CALL GetRGBfromNormVal(Val2Unity(rVal,rLB,rUB,nCMethod),color)
!!!
!!!    bmLegend(2:nSX-5,(1+i),1) = color%r
!!!    bmLegend(2:nSX-5,(1+i),2) = color%g
!!!    bmLegend(2:nSX-5,(1+i),3) = color%b
!!!
!!!  END DO
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!! copy from legend bitmap to image bitmap
!!!  bmDest(rcEffective%lb%x:rcEffective%rt%x, &
!!!     &   rcEffective%lb%y:rcEffective%rt%y,:) &
!!!     &  = bmLegend(rcE2%lb%x:rcE2%rt%x,rcE2%lb%y:rcE2%rt%y,:)
!!!!      write(*,*)"DrawLegend24b: BltDest: ([",rcEffective%lb%x,",",rcEffective%lb%y,"],[",rcEffective%rt%x,",",rcEffective%rt%y,"])"
!!!!      write(*,*)"DrawLegend24b: BltLeg: ([  1 ,  1 ],[",(rcEffective%rt%x-rcEffective%lb%x),",",(rcEffective%rt%y-rcEffective%lb%y),"])"
!!!! ------------------
!!!
!!!
!!!
!!!
!!!! ------------------
!!!  DEALLOCATE(bmLegend,STAT=alloc)
!!!  IF (alloc /= 0) THEN
!!!    write(*,*) "DrawLegend24b: DeAllocation error."
!!!    RETURN
!!!  END IF
!!!! ------------------
!!!
!!!
!!!
!!!! ------------------
!!!! reset color
!!!  bfcCurrent = bfcOld
!!!! ------------------
!!!
!!!
!!!    
!!!  END SUBROUTINE DrawLegend24b
!!!!**********************************************************************!


































































!**********************************************************************!
!************************* FILE FUNCTIONS *****************************!
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE CorrectEndian16(ndata,nx,ny)
! function: corrects wrong data endian of an image
!           works on the assumtion, that for the right endian an image
!           has lower deviation values
!           check both endians of ndata
!    !!!    ndata will be changed if endian is wrong
! -------------------------------------------------------------------- !
! parameter: ndata(nx,ny) : integer*2 : data array
!            nx, ny : integer*4 : array size
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*2, intent(inout) :: ndata(nx,ny)
  integer*4, intent(in) :: nx, ny

  integer*2 :: nval, nswapval, nswapdata(nx,ny)
  integer*4 :: i, j, n
  real*4 :: mean, swapmean, dev, swapdev, rval
  
  integer*2, external :: SwapInt2
! ------------



! ------------
! init
!  write(*,*) " > CorrectEndian16: INIT."
  n = nx*ny
  if (n<=0 .or. nx<=0) return ! invalid size input catched
! ------------


! ------------
! set swap array, and get means
  mean = 0.0
  swapmean = 0.0

  do j=1, ny
    do i=1, nx
      nval = ndata(i,j)
      nswapval = SwapInt2(nval)
      nswapdata(i,j) = nswapval
      mean = mean + real(nval)
      swapmean = swapmean + real(nswapval)
    end do
  end do
  mean = mean / real(n)
  swapmean = swapmean / real(n)
! ------------


! ------------
! get deviation values
  dev = 0.0
  swapdev = 0.0

  do j=1, ny
    do i=1, nx
      rval = real(ndata(i,j)) - mean
      dev = dev + rval * rval
      rval = real(nswapdata(i,j)) - swapmean
      swapdev = swapdev + rval * rval
    end do
  end do
  dev = sqrt(dev / real(n-1))
  swapdev = sqrt(swapdev / real(n-1))
! ------------


! ------------
! set right endian data
  if (swapdev<dev) ndata = nswapdata
! ------------

! ------------
!  write(*,*) " > CorrectEndian16: EXIT."
  return

END SUBROUTINE CorrectEndian16
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE CorrectEndian32(ndata,nx,ny)
! function: corrects wrong data endian of an image
!           works on the assumtion, that for the right endian an image
!           has lower deviation values
!           check both endians of ndata
!    !!!    ndata will be changed if endian is wrong
! -------------------------------------------------------------------- !
! parameter: ndata(nx,ny) : integer*4 : data array
!            nx, ny : integer*4 : array size
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(inout) :: ndata(nx,ny)
  integer*4, intent(in) :: nx, ny

  integer*4 :: nval, nswapval, nswapdata(nx,ny)
  integer*4 :: i, j, n
  real*4 :: mean, swapmean, dev, swapdev, rval
  
  integer*4, external :: SwapInt4
! ------------



! ------------
! init
!  write(*,*) " > CorrectEndian32: INIT."
  n = nx*ny
  if (n<=0 .or. nx<=0) return ! invalid size input catched
! ------------


! ------------
! set swap array, and get means
  mean = 0.0
  swapmean = 0.0

  do j=1, ny
    do i=1, nx
      nval = ndata(i,j)
      nswapval = SwapInt4(nval)
      nswapdata(i,j) = nswapval
      mean = mean + real(nval)
      swapmean = swapmean + real(nswapval)
    end do
  end do
  mean = mean / real(n)
  swapmean = swapmean / real(n)
! ------------


! ------------
! get deviation values
  dev = 0.0
  swapdev = 0.0

  do j=1, ny
    do i=1, nx
      rval = real(ndata(i,j)) - mean
      dev = dev + rval * rval
      rval = real(nswapdata(i,j)) - swapmean
      swapdev = swapdev + rval * rval
    end do
  end do
  dev = sqrt(dev / real(n-1))
  swapdev = sqrt(swapdev / real(n-1))
! ------------


! ------------
! set right endian data
  if (swapdev<dev) ndata = nswapdata
! ------------

! ------------
!  write(*,*) " > CorrectEndian32: EXIT."
  return

END SUBROUTINE CorrectEndian32
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE raw_open(nunit,sname)
! function: connects nunit to the file sname using formatted direct
!           access with record length 1
! -------------------------------------------------------------------- !
! parameter:
!           nunit : integer*4, input : associated file unit
!           sanme : character(len=*), input : file name to open
! -------------------------------------------------------------------- !


  implicit none


! ------------
! declaration
  integer*4, intent(in) :: nunit
  character(len=*), intent(in) :: sname

  integer*4 :: ios
  logical :: isthere
! ------------



! ------------
! init
!  write(*,*) "raw_open: INIT."
! ------------


! ------------
! EXISTENCE CHECK
  inquire(file=TRIM(sname),exist=isthere)
  
  if (.not.isthere) then
    write(*,*) "raw_open: ",sname," does not exist."
    stop
  end if
! ------------



! ------------
! connecting
  open ( unit = nunit, file = sname, form = 'formatted', &
     &   access = 'direct', recl = 1, status = 'old', iostat = ios )
  if (ios/=0) then
    write(*,*) "raw_open: ERROR - opening ",sname," failed."
    stop
  end if
! ------------


! ------------
!  write(*,*) "raw_open: EXIT."
  return

END SUBROUTINE raw_open
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE raw_c_read(nunit,pos,c,n)
! function: reads n characters to c starting from file position pos
!           in nunit
! -------------------------------------------------------------------- !
! parameter:
!            nunit : integer*4 : file unit
!            pos : integer*4 : position in file
!            c(n) : character(n) : character buffer
!            n : integer*4 : number of chars to read
! -------------------------------------------------------------------- !

  implicit none


! ------------
! declaration
  integer*4 :: nunit, pos, n
  character :: c(n)
  
  integer*4 :: i, ios
  logical :: unit_open
! ------------



! ------------
! init
!  write(*,*) "raw_c_read: INIT."
  inquire(unit=nunit,opened=unit_open)
  if (.not.(unit_open)) then
    write(*,*) "raw_c_read: file not opened to read.!"
    stop
  end if
  if (pos<=0) then
    write(*,*) "raw_c_read: illegal record number:", pos
    stop
  end if
! ------------


! ------------
  do i = 1, n

    read ( nunit, rec = pos, fmt = '(a)', iostat = ios ) c(i)

    if ( ios /= 0 ) then
      c(i:n) = CHAR(0)
      pos = -1
      return
    end if

    pos = pos + 1

  end do
! ------------


! ------------
!  write(*,*) "raw_c_read: EXIT."
  return

END SUBROUTINE raw_c_read
!**********************************************************************!






!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE raw_uread16(sfile, datablock, nsize, fswap, errorcode)
!! function: reads 16bit data from file to array
!! -------------------------------------------------------------------- !
!! parameter: sfile : character(*) : name of the file to read
!!            datablock(ndatasize) : integer*2 : data reference
!!            ndatasize=nsize*nsize : integer*4 : number of 16bit values to read
!!            fswap : integer*4 : byte swap flag (0=noswap, 1=swap)
!!            errorcode : integer*4 : returncode
!! -------------------------------------------------------------------- !
!! link : ImgTranfs.f95
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! declaration
!  integer*4, parameter :: uf = 99
!  integer*4, parameter :: ndb = 2
!
!  character(len=*) :: sfile
!  integer*1 :: nidata(ndb*nsize*nsize), nbuffer(ndb*nsize)
!  integer*1 :: checker(ndb*nsize), nt(2)
!  integer*2 :: datablock(nsize,nsize), mold
!  integer*4 :: nsize, fswap, errorcode, nfilesize, ndatasize
!  integer*4 :: nshift, nbsize
!  
!  integer*4 :: i,j,i1,l
!  logical :: isok
!  
!  external :: CorrectEndian16
!! ------------
!
!
!
!! ------------
!! init, open file
!!  write(*,*) " > raw_uread16: INIT."
!! check if requested file exists  
!  inquire(file=trim(sfile), exist=isok)
!  
!! set block size of records to read (in bytes)
!  nbsize = ndb*nsize
!
!! set return array to zero
!  datablock(:,:)=0
!  
!! if the file exists try to open
!  if (isok) then
!
!!    open (unit=uf, file=trim(sfile), action="read", iostat=errorcode, &
!!     & form='unformatted', status='old', position='rewind')
!    
!    ! open in unformatted mode direct access for reading with nbsize record length
!    open(uf,FILE=trim(sfile),STATUS='old', ACCESS='direct', &
!     & form='unformatted',action = 'read', &
!     & iostat=errorcode, RECL=nbsize)
!
!    ! check opening for success
!    if (errorcode/=0) then
!      write(*,*) " > raw_uread16: some error occured during opening ", trim(sfile)
!      write(*,*) " >              Errorcode:",errorcode
!      return
!    end if
!  else
!    ! the file does not exist, throw error message to console and return with 
!    ! -3 as error code
!    write(*,*) " > raw_uread16: ERROR: File ",trim(sfile)," does not exist."
!    errorcode = -3
!    return
!  end if
!
!! ------------
!
!
!! ------------
!! read the data 
!  i = 1
!  j = 1 ! dataline position counter
!  l = 1 ! file length counter
!  nidata(:) = 0
!  nfilesize = 0 ! file length to remember
!  ndatasize = nsize*nbsize ! size of data that has to be read
!
!! buffer-wise until eof
!! !! AS no error code will be generated if the data remaining in the file is
!!    read into a larger buffer (so eof is reached by array-read internally) 
!!    we first have to check-read a record. If the check read is ok, the previous
!!    buffer-read is valid and can be transfered fully into the datablock/nidata
!!    if the checkread failes with error code, the requested record is fully outside
!!    the file size, BUT eventually the previously read record is not completely
!!    full with data. So we skip and reopen later for single byte reading.
!  checker(:) = 0
!  nbuffer(:) = 0
!! primary check-read
!  read(unit=uf,iostat=errorcode,REC=l) checker
!! if check-read is ok, transfer data to nbuffer
!  if (errorcode==0) nbuffer(1:nbsize) = checker(1:nbsize)
!! if check read is ok, try next check read and proceed until error occures  
!  if (errorcode==0) then
!    l = l + 1 ! next record
!    do
!      ! check-read next record
!      read(unit=uf, iostat=errorcode, REC=l) checker
!      ! if check-read failed, exit the do-loop ...
!      if (errorcode/=0) exit
!    
!      ! check-read was ok, transfer previous buffer data to full data line
!      nidata(j:j+nbsize-1) = nbuffer(1:nbsize)
!      ! transfer current check data to buffer
!      nbuffer(1:nbsize) = checker(1:nbsize)
!    
!      nfilesize = nfilesize + nbsize ! nbsize bytes were successfully transfered
!    
!      l = l + 1 ! next record
!      j = j + nbsize ! next position in full-data line
!      ! The following is not very clever, keep an eye on it if altering the code!!
!      ! It is checked, whether the internal data size counter is raised above the
!      ! requested size of the data array. If so, the data should be wrapped around
!      ! to the beginning of the datablok/line. In the current case the file size is
!      ! a multiple of the block size nbsize, so wrapping is consistent with
!      ! setting the next position to 1 if the blocksize was reached previously. A
!      ! more detailed handling of the data should be implemented in cases with variable
!      ! data sizes and variable read-buffer sizes!!
!      if (j>ndatasize) j = 1
!    
!    end do
!  end if
!  !write(*,*) errorcode
!  !write(*,*) nfilesize, l, nbsize
!  errorcode = 0
!
!! eof reached, reopen to read last bytes bytewise
!  close(unit=uf) ! first close
!  ! then open again
!  open(uf,FILE=trim(sfile),STATUS='old', ACCESS='direct', &
!     & form='unformatted',action = 'read', &
!     & iostat=errorcode, RECL=1)
!
!  l = nfilesize + 1
!  i = 1
!  do
!    read(unit=uf, iostat=errorcode, REC=l) nbuffer(1)
!    
!    if (errorcode/=0) exit
!    
!    ! transfer from buffer to data line
!    nidata(j:j) = nbuffer(1:1)
!    
!    l = l + 1
!    j = j + 1
!    
!    ! wrap to file begin
!    if (j>ndatasize) j = 1
!    
!  end do
!  !write(*,*) errorcode
!  
!  errorcode = 0
!  l = l - 1
!  nfilesize = l
!  !write(*,*) nfilesize,j
!  !write(*,*) nfilesize, ndatasize, nfilesize-ndatasize
!
!! ------------
!
!! ------------
!! close file
!  inquire(unit=uf, opened=isok)
!  if (isok) then
!    close(uf)
!  else
!    write(*,*) " > raw_uread16: warning: File ",trim(sfile)," could not be closed."
!  end if
!! ------------
!
!
!
!! ------------
!! reorder
!
!  if (nfilesize<ndatasize) then
!    write(*,*) " > raw_uread16: ERROR: file corrupt or given size is wrong."
!    return
!  end if
!
!  nshift = mod(nfilesize-ndatasize,ndatasize)
!  !write(*,*) nshift
!
!  i1 = 1+nshift
!
!  do j=1, nsize
!    do i=1, nsize
!    
!      nt(1) = nidata(i1)
!      i1 = i1 + 1
!      if (i1>ndatasize) i1 = 1
!      nt(2) = nidata(i1)      
!      i1 = i1 + 1
!      if (i1>ndatasize) i1 = 1
!      
!      datablock(i,j) = TRANSFER(nt(1:2),mold)
!
!    end do
!  end do
!! ------------
!
!
!
!! ------------
!! swap if fswap internal mem
!  call CorrectEndian16(datablock,nsize,nsize)
!! ------------
!
!
!! ------------
!!  write(*,*) " > raw_uread16: EXIT."
!  return
!
!END SUBROUTINE raw_uread16
!!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE raw_uwrite16(sfile, datablock, nsize, fswap, errorcode)
! function: writess 16bit data from array to file
! -------------------------------------------------------------------- !
! parameter: sfile : character(*) : name of the file to read
!            datablock(nsize,nsize) : integer*2 : data reference
!            nsize : integer*4 : quadr. array size
!            fswap : integer*4 : byte swap flag (0=noswap, 1=swap)
!            errorcode : integer*4 : returncode
! -------------------------------------------------------------------- !

  implicit none


! ------------
! declaration
  integer*4, parameter :: uf = 99
  integer*4, parameter :: nb = 2

  character(len=*) :: sfile
  integer*2 :: datablock(nsize,nsize)
  integer*2 :: idatablock(nsize*nsize), iwdata(nsize*nsize)
  integer*4 :: nsize, fswap, errorcode
  
  integer*4 :: i,j,k,ndatasize
  logical :: isok
  
  integer*2, external :: SwapInt2
! ------------



! ------------
! init, open file
  ndatasize = nb*nsize*nsize
  k = 1
  do j=1, nsize
    do i=1, nsize
      idatablock(k) = datablock(i,j)
      k=k+1
    end do
  end do
!  write(*,*) "raw_uwrite16: INIT."
  inquire(file=trim(sfile), exist=isok)

  if (isok) then
    write(*,*) "raw_uwrite16: file ",trim(sfile)," will be overwritten."
  end if

!  open (unit=uf, file=trim(sfile), action="read", iostat=errorcode, &
!   & form='unformatted' )
  open(uf,FILE=trim(sfile), ACCESS='direct', form='unformatted', &
     & ACTION='write', iostat=errorcode, RECL=ndatasize)

  if (errorcode/=0) then
    write(*,*) "raw_uwrite16: some error occured during opening ", trim(sfile)
    write(*,*) "              program stopped. Errorcode:",errorcode
    stop
  end if
! ------------


! ------------
! swap if fswap
  if (fswap==1) then
    do i=1, nsize*nsize
      iwdata(i) = SwapInt2(idatablock(i))
    end do
  else
    do i=1, nsize*nsize
      iwdata(i) = idatablock(i)
    end do
  end if
! ------------



! ------------
! writing
!  write(*,*) "raw_uwrite16: try writing", nsize*nsize," 16-bit values."
   write(uf,REC=1,iostat=errorcode) iwdata
!  write(uf,iostat=errorcode) datablock

!  write(*,*) wdata(1,1),datablock(1,1)
  if (errorcode/=0) then
    write(*,*) "raw_uwrite16: some error occured during writing from ", trim(sfile)
    write(*,*) "             Errorcode:",errorcode
!    stop
  end if

! ------------



! ------------
! close file
  inquire(unit=uf, opened=isok)
  if (isok) then
    close(uf)
  else
    write(*,*) "raw_uwrite16: minor warning: File ",trim(sfile)," could not be closed."
    stop
  end if
! ------------


! ------------
!  write(*,*) "raw_uwrite16: EXIT."
  return

END SUBROUTINE raw_uwrite16
!**********************************************************************!





!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE raw_uread32(sfile, datablock, nsize, fswap, errorcode)
!! function: reads 32bit data from file to array
!! -------------------------------------------------------------------- !
!! parameter: sfile : character(*) : name of the file to read
!!            datablock(nsize,nsize) : integer*4 : data reference
!!            nsize : integer*4 : quadratic array size
!!            fswap : integer*4 : byte swap flag (0=noswap, 1=swap)
!!            errorcode : integer*4 : returncode
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! declaration
!  integer*4, parameter :: uf = 99
!  integer*4, parameter :: ndb = 4
!
!  character(len=*) :: sfile
!  integer*1 :: nidata(ndb*nsize*nsize), nbuffer(ndb*nsize)
!  integer*1 :: checker(ndb*nsize), nt(4)
!  integer*4 :: datablock(nsize,nsize), mold
!  integer*4 :: fswap, errorcode, nfilesize, nsize, ndatasize
!  integer*4 :: nshift, nbsize
!  
!  integer*4 :: i,j,i1,l
!  logical :: isok
!  
!  external :: CorrectEndian16
!! ------------
!
!
!
!! ------------
!! init, open file
!!  write(*,*) " > raw_uread32: INIT."
!! check if requested file exists  
!  inquire(file=trim(sfile), exist=isok)
!  
!! set block size of records to read (in bytes)
!  nbsize = ndb*nsize
!
!! set return array to zero
!  datablock(:,:)=0
!  
!! if the file exists try to open
!  if (isok) then
!
!!    open (unit=uf, file=trim(sfile), action="read", iostat=errorcode, &
!!     & form='unformatted', status='old', position='rewind')
!    
!    ! open in unformatted mode direct access for reading with nbsize record length
!    open(uf,FILE=trim(sfile),STATUS='old', ACCESS='direct', &
!     & form='unformatted',action = 'read', &
!     & iostat=errorcode, RECL=nbsize)
!
!    ! check opening for success
!    if (errorcode/=0) then
!      write(*,*) " > raw_uread32: some error occured during opening ", trim(sfile)
!      write(*,*) " >              Errorcode:",errorcode
!      return
!    end if
!  else
!    ! the file does not exist, throw error message to console and return with 
!    ! -3 as error code
!    write(*,*) " > raw_uread32: ERROR: File ",trim(sfile)," does not exist."
!    errorcode = -3
!    return
!  end if
!
!! ------------
!
!
!! ------------
!! read the data 
!  i = 1
!  j = 1 ! dataline position counter
!  l = 1 ! file length counter
!  nidata(:) = 0
!  nfilesize = 0 ! file length to remember
!  ndatasize = nsize*nbsize ! size of data that has to be read
!
!! buffer-wise until eof
!! !! AS no error code will be generated if the data remaining in the file is
!!    read into a larger buffer (so eof is reached by array-read internally) 
!!    we first have to check-read a record. If the check read is ok, the previous
!!    buffer-read is valid and can be transfered fully into the datablock/nidata
!!    if the checkread failes with error code, the requested record is fully outside
!!    the file size, BUT eventually the previously read record is not completely
!!    full with data. So we skip and reopen later for single byte reading.
!  checker(:) = 0
!  nbuffer(:) = 0
!! primary check-read
!  read(unit=uf,iostat=errorcode,REC=l) checker
!! if check-read is ok, transfer data to nbuffer
!  if (errorcode==0) nbuffer(1:nbsize) = checker(1:nbsize)
!! if check read is ok, try next check read and proceed until error occures  
!  if (errorcode==0) then
!    l = l + 1 ! next record
!    do
!      ! check-read next record
!      read(unit=uf, iostat=errorcode, REC=l) checker
!      ! if check-read failed, exit the do-loop ...
!      if (errorcode/=0) exit
!    
!      ! check-read was ok, transfer previous buffer data to full data line
!      nidata(j:j+nbsize-1) = nbuffer(1:nbsize)
!      ! transfer current check data to buffer
!      nbuffer(1:nbsize) = checker(1:nbsize)
!    
!      nfilesize = nfilesize + nbsize ! nbsize bytes were successfully transfered
!    
!      l = l + 1 ! next record
!      j = j + nbsize ! next position in full-data line
!      ! The following is not very clever, keep an eye on it if altering the code!!
!      ! It is checked, whether the internal data size counter is raised above the
!      ! requested size of the data array. If so, the data should be wrapped around
!      ! to the beginning of the datablok/line. In the current case the file size is
!      ! a multiple of the block size nbsize, so wrapping is consistent with
!      ! setting the next position to 1 if the blocksize was reached previously. A
!      ! more detailed handling of the data should be implemented in cases with variable
!      ! data sizes and variable read-buffer sizes!!
!      if (j>ndatasize) j = 1
!    
!    end do
!  end if
!  !write(*,*) errorcode
!  !write(*,*) nfilesize, l, nbsize
!  errorcode = 0
!
!! eof reached, reopen to read last bytes bytewise
!  close(unit=uf) ! first close
!  ! then open again
!  open(uf,FILE=trim(sfile),STATUS='old', ACCESS='direct', &
!     & form='unformatted',action = 'read', &
!     & iostat=errorcode, RECL=1)
!
!  l = nfilesize + 1
!  i = 1
!  do
!    read(unit=uf, iostat=errorcode, REC=l) nbuffer(i)
!    
!    if (errorcode/=0) exit
!    
!    i = i + 1
!    l = l + 1
!    
!  end do
!  !write(*,*) errorcode
!  
!  errorcode = 0
!  i = i - 1 
!  l = l - 1
!  nfilesize = l
!  !write(*,*) nfilesize,i
!
!  if (i>0) then
!    nidata(j:j+i-1) = nbuffer(1:i)
!  end if
!  
!  
!  !write(*,*) nfilesize, ndatasize, nfilesize-ndatasize
!
!! ------------
!
!
!
!! ------------
!! close file
!  inquire(unit=uf, opened=isok)
!  if (isok) then
!    close(uf)
!  else
!    write(*,*) " > raw_uread32: warning: File ",trim(sfile)," could not be closed."
!  end if
!! ------------
!
!
!
!! ------------
!! reorder
!
!  if (nfilesize<ndatasize) then
!    write(*,*) " > raw_uread32: ERROR: file corrupt or given size is wrong."
!    return
!  end if
!
!  nshift = mod(nfilesize-ndatasize,ndatasize)
!  !write(*,*) nshift
!
!  i1 = 1+nshift
!
!  do j=1, nsize
!    do i=1, nsize
!    
!      nt(1) = nidata(i1)
!      i1 = i1 + 1
!      if (i1>ndatasize) i1 = 1
!      nt(2) = nidata(i1)      
!      i1 = i1 + 1
!      if (i1>ndatasize) i1 = 1
!      
!      datablock(i,j) = TRANSFER(nt(1:4),mold)
!
!    end do
!  end do
!! ------------
!
!
!
!! ------------
!! swap if fswap internal mem
!  call CorrectEndian32(datablock,nsize,nsize)
!! ------------
!
!
!! ------------
!!  write(*,*) " > raw_uread32: EXIT."
!  return
!
!END SUBROUTINE raw_uread32
!!**********************************************************************!




















































































!**********************************************************************!
!************************* NUMERCIAL FUNCTIONS ************************!
!**********************************************************************!

!**********************************************************************!
! transforms high tension [kV] to electron wavelength [nm]
real*4 function HT2WL(ht)
implicit none
real*8, parameter :: te0    = 1.021997856D+003 ! 2 * m_e * c^2 [keV]
real*8, parameter :: hc     = 1.239842041D+000 ! h*c [nm*keV]
real*4,intent(in) :: ht ! [kV]
real*8 :: dht
dht = dble(ht)
HT2WL = real( hc / dsqrt( dht*( te0 + dht ) ), kind=4 )
return
end function HT2WL

!**********************************************************************!
! transforms electron wavelength [nm] to high tension [kV]
real*4 function WL2HT(wl)
implicit none
real*8, parameter :: e0     = 0.510998928D+003 ! m_e * c^2 [keV]
real*8, parameter :: hcoe0  = 2.426310454D-002 ! h*c / (m_e * c^2) [nm]
real*4, intent(in):: wl ! [nm]
real*8 :: hcole0
hcole0 = hcoe0/dble(wl) ! = h*c/e0/lambda [1]
WL2HT = real( e0 * ( dsqrt( 1.D0 + hcole0*hcole0 ) - 1.D0 ) , kind=4)
return
end function WL2HT

!**********************************************************************!
!**********************************************************************!
SUBROUTINE getMeanandSigma(dat,n,mean,sigma)
! function: calculates mean and sigma of array dat
! -------------------------------------------------------------------- !
! parameter: dat(n) : real*4 : input array
!            n : integer*4 : size of dat
!            mean : real*4 : return ref. for mean
!            sigma : real*4 : return ref. for sigma
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n
  real*4, intent(in) :: dat(n)
  real*4, intent(out) :: mean, sigma
  integer*4 :: i
  real*4 :: sqs, s, tmp
! ------------

! ------------
! init
!  write(*,*) " > getMeanandSigma: INIT."
  sqs = 0.0
  s = 0.0
  mean = 0.0
  sigma = 0.0
! ------------

! ------------
! loop over array and calculate sum and square sum
  do i=1, n
    tmp = dat(i)
    s = s + tmp
    sqs = sqs + tmp*tmp
  end do
! ------------

! ------------
! calculate mean and sigma from the sums
  mean = s / real(n)
  sigma = sqrt( sqs / real(n) - mean*mean )
! ------------

! ------------
!  write(*,*) " > getMeanandSigma: EXIT."
  return

END SUBROUTINE getMeanandSigma
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE getMeanandSigma2(dat,sig,n,mean,sigma)
! function: calculates mean and sigma of array dat
!           including sample sigmas sig to calculate total sigma
! -------------------------------------------------------------------- !
! parameter: dat(n) : real*4 : input data array
!            sig(n) : real*4 : input sigma array
!            n : integer*4 : size of dat
!            mean : real*4 : return ref. for mean
!            sigma : real*4 : return ref. for sigma
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n
  real*4, intent(in) :: dat(n), sig(n)
  real*4, intent(out) :: mean, sigma
  integer*4 :: i
  real*4 :: s, tmp1, tmp2
! ------------

! ------------
! init
!  write(*,*) " > getMeanandSigma2: INIT."
  mean = 0.0
  sigma = 0.0
! ------------

! ------------
! loop over array and calculate mean from sum
  s = 0.0
  do i=1, n
    tmp1 = dat(i)
    s = s + tmp1
  end do
  mean = s / real(n)
! loop over array and calculate sigma from mean and sqr sum
  s = 0.0
  do i=1, n
    tmp1 = dat(i) - mean
    tmp2 = sig(i)
    tmp1 = tmp1 + tmp2 / ( tmp1 + (2.*tmp2)**(1./3.) )**2.
    s = s + tmp1*tmp1
  end do
  sigma = sqrt( s / real(n) )
! ------------

! ------------
!  write(*,*) " > getMeanandSigma2: EXIT."
  return

END SUBROUTINE getMeanandSigma2
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE checkprimefactors(n,pf,npf,nresid)
! function: removes all primefactors from n which are contained in
!           pf(npf) and returns the residuum number composed of other
!           prime factors.
!           pf must be a list of prime numbers
! -------------------------------------------------------------------- !
! parameter: n : integer*4 : input number start
!            pf(npf) : integer*4 : input prime factors to check
!            npf : integer*4 : input number of prime factors to check
!            nresid : integer*4 : residual of the factorization with pf
! -------------------------------------------------------------------- !
!
!!! RECURSIVE !!!
!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n, npf, pf(npf)
  integer*4, intent(out) :: nresid
  integer*4 :: i, j

! ------------
! init
!  write(*,*) " > checkprimefactors: INIT."
  nresid = n
  if (n<2) return
  if (npf<1) return

! ------------
  do i=1, npf
    if ( 0==modulo(n,pf(i)) ) then
      j = FLOOR( real(n)/real(pf(i)) )
      call checkprimefactors(j,pf,npf,nresid)
    end if
  end do
! ------------

! ------------
!  write(*,*) " > checkprimefactors: EXIT."
  return

END SUBROUTINE checkprimefactors
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
integer*4 FUNCTION nextprime(n,inc)
! function: calculate next prime
! -------------------------------------------------------------------- !
! parameter: n : integer*4 : input number start
!            inc : integer*4 : search increment
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, intent(in) :: n, inc
  integer*4 :: nr, i, j
! ------------



! ------------
! init
!  write(*,*) " > nextprime: INIT."
  if (inc==0) return
! ------------


! ------------
  nr = n
  do
    j = int(sqrt(real(nr)))+1
    do i=2,j
      if (modulo(nr,i)==0) exit
    end do
    if (i>=j) exit
    nr = nr - 1
  end do
  nextprime = nr
! ------------

! ------------
!  write(*,*) " > nextprime: EXIT."
  return

END FUNCTION nextprime
!**********************************************************************!


!!**********************************************************************!
!!**********************************************************************!
!real*4 FUNCTION gamma(x)
!! function: Returns the complete gamma function G(x)
!! -------------------------------------------------------------------- !
!! parameters: x : independent variable
!! -------------------------------------------------------------------- !
!
!  implicit none
!  
!  i
!!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
real*4 FUNCTION gammp(a,x)
! function: Returns the incomplete gamma function P(a, x).
! -------------------------------------------------------------------- !
! parameter: a : real*4 : dimension parameter
!            x : minimized merit
! -------------------------------------------------------------------- !
!            USES gcf,gser
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: a,x

  real*4 gammcf,gamser,gln
! ------------

  
! ------------
  gammp = 1.0
  if(x<0.0 .or. a<=0.) return ! bad arguments in gammp
  if(x<a+1.) then ! Use the series representation.
    call gser(gamser,a,x,gln)
    gammp=gamser
  else ! Use the continued fraction representation
    call gcf(gammcf,a,x,gln)
    gammp=1.-gammcf ! and take its complement.
  end if
! ------------


! ------------
  return

END FUNCTION gammp
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
real*4 FUNCTION gammq(a,x)
! function: Returns the incomplete gamma function Q(a, x) = 1 - P(a, x)
! -------------------------------------------------------------------- !
! parameter: a : real*4 : dimension parameter
!            x : minimized merit
! -------------------------------------------------------------------- !
!            USES gcf,gser
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: a,x

  real*4 gammcf,gamser,gln
! ------------


! ------------
  gammq = 0.0
  if (x<0.0 .or. a<=0.0) return ! bad arguments in gammq
  if (x<a+1.) then ! Use the series representation
    call gser(gamser,a,x,gln)
    gammq=1.-gamser ! and take its complement.
  else ! Use the continued fraction representation.
    call gcf(gammcf,a,x,gln)
    gammq=gammcf
  end if
! ------------


! ------------
  return

END FUNCTION gammq
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE gser(gamser,a,x,gln)
! function: Returns the incomplete gamma function P(a, x) evaluated
!           by its series representation as gamser. Also returns 
!           lnGAMMA(a) as gln.
! -------------------------------------------------------------------- !
! parameter: gammser : real*4 : series representation of P(a,x)
!            a : real*4 : dimension parameter
!            x : real*4 : minimized merit
!            gln : real*4 : log of GAMMA(a)
! -------------------------------------------------------------------- !
!            USES gammln
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, parameter :: ITMAX = 100
  real*4, parameter :: EPS = 3.e-7

  real*4, intent(in) :: a,x
  real*4, intent(inout) :: gamser,gln
  
  integer*4 :: n
  real*4 :: ap,del,dsum,gammln
! ------------


! ------------
! error catch
  gln=gammln(a)
  if(x<=0.0) then
    if(x<0.0) return ! x < 0 in gser
    gamser=0.0
    return
  end if
! ------------


! ------------
  ap=a
  dsum=1./a
  del=dsum
  do n=1,ITMAX
    ap=ap+1.
    del=del*x/ap
    dsum=dsum+del
    if(abs(del) < abs(dsum)*EPS) goto 1
  end do
  write(*,*) " > gser: <a> too large, ITMAX too small."
  
1 gamser=dsum*exp(-x+a*log(x)-gln)
! ------------

! ------------
  return

END SUBROUTINE gser
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE gcf(gammcf,a,x,gln)
! function: Returns the incomplete gamma function Q(a, x) evaluated
!           by its continued fraction representation as gammcf.
!           Also returns lnGAMMA(a) as gln.
! -------------------------------------------------------------------- !
! parameter: gammser : real*4 : series representation of P(a,x)
!            a : real*4 : dimension parameter
!            x : real*4 : minimized merit
!            gln : real*4 : log of GAMMA(a)
! -------------------------------------------------------------------- !
!            USES gammln
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4, parameter :: ITMAX = 100
  real*4, parameter :: EPS = 3.e-7
  real*4, parameter :: FPMIN = 1.e-30

  real*4, intent(in) :: a, x
  real*4, intent(inout) :: gammcf, gln

  integer*4 :: i
  real*4 :: an,b,c,d,del,h
  real*4, external :: gammln
! ------------


! ------------
  gln=gammln(a)
  b=x+1.-a ! Set up for evaluating continued fraction by modified
           ! Lentz’s method (§5.2) with b0 = 0.
  c=1./FPMIN
  d=1./b
  h=d
! ------------

! ------------
  do i=1,ITMAX ! Iterate to convergence.
    an=-i*(i-a)
    b=b+2.
    d=an*d+b
    if(abs(d) < FPMIN) d = FPMIN
    c=b+an/c
    if(abs(c) < FPMIN) c = FPMIN
    d=1./d
    del=d*c
    h=h*del
    if(abs(del-1.) < EPS) goto 1
  end do
  write(*,*) " > gcf: <a> too large, ITMAX to small."
! ------------


! ------------
1 gammcf=exp(-x+a*log(x)-gln)*h ! Put factors in front.
! ------------


! ------------
  return

END SUBROUTINE gcf
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
real*4 FUNCTION gammln(xx)
! function: Returns the value ln[GAMMA(xx)] for xx > 0.
! -------------------------------------------------------------------- !
! parameter: xx : real*4 : >0
!            
! -------------------------------------------------------------------- !

  implicit none
  
! ------------
! declaration
  real*4, intent(in) :: xx

  integer*4 :: j
  real*8 :: ser,stp,tmp,x,y,cof(6)

! ------------


! ------------
  cof = (/ 76.18009172947146d0, -86.50532032941677d0, &
     &  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
     &  -.5395239384953d-5 /)
  stp = 2.5066282746310005d0
! ------------


! ------------
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! ------------


! ------------
  do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  end do
! ------------


! ------------
  gammln=tmp+log(stp*ser/x)
! ------------

! ------------
return

END FUNCTION gammln
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
real*4 FUNCTION rgauss2d(x,y,A0,A,x0,y0,w)
! function: calculates values of the gaussian function
!           f(x,y) = A0+A*exp( - ( (x-x0)**2 + (y-y0)**2 )/ w**2)
!           integrated over area F if radius of F is much larger than w
!           gives F*A0 + A*Pi*w^2
! -------------------------------------------------------------------- !
! parameter:
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4, intent(in) :: x,y,A0,A,x0,y0,w
  real*4 :: dx, dy
! ------------



! ------------
! init
!  write(*,*) " > rgauss2d: INIT."
! ------------


! ------------
  dx = x-x0
  dy = y-y0
  if (abs(w)>0.0) then
    rgauss2d = A0 + A*exp(-(dx*dx+dy*dy)/(w*w))
  else
    rgauss2d = A0
  end if
! ------------


! ------------
!  write(*,*) " > rgauss2d: EXIT."
!  return

END FUNCTION rgauss2d
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION factorial(n)
! function: calculates the factorial of n -> n!
! -------------------------------------------------------------------- !
! parameter: integer*4 :: n
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: factorial
  integer*4, intent(in) :: n
  integer*4 :: i
! ------------

! ------------
! init
!  write(*,*) " > factorial: INIT."
  factorial = 0 ! precheck default -> this means ERROR!
  if (n<0) return
  factorial = 1 ! postcheck default -> this means NO ERROR!
  i=2
! ------------

! ------------
  do while (n>=i)
    factorial = factorial * i
    i = i + 1
  end do
! ------------

! ------------
!  write(*,*) " > factorial: EXIT."
  return

END FUNCTION factorial
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION dfactorial(n)
! function: calculates the factorial of n -> n!
! -------------------------------------------------------------------- !
! parameter: integer*8 :: n
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*8 :: dfactorial
  integer*8, intent(in) :: n
  integer*8 :: i
! ------------

! ------------
! init
!  write(*,*) " > dfactorial: INIT."
  dfactorial = 0 ! precheck default -> this means ERROR!
  if (n<0) return
  dfactorial = 1 ! postcheck default -> this means NO ERROR!
  i=2
! ------------

! ------------
  do while (n>=i)
    dfactorial = dfactorial * i
    i = i + 1
  end do
! ------------

! ------------
!  write(*,*) " > dfactorial: EXIT."
  return

END FUNCTION dfactorial
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION binomial(n,k)
! function: calculates the binomial coefficient of (n over k), which is
!           equal to (n!)/( (n-k)! * k! )
! -------------------------------------------------------------------- !
! parameter: integer*4 :: n,k
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: binomial
  integer*4, intent(in) :: n, k
  integer*4, external :: factorial
! ------------

! ------------
! init
!  write(*,*) " > binomial: INIT."
  binomial = 0 ! precheck default -> this means ERROR!
  if (n<0.or.k<0.or.n<k) return
  binomial = factorial(n)/( factorial(n-k) * factorial(k) )
! ------------

! ------------
!  write(*,*) " > binomial: EXIT."
  return

END FUNCTION binomial
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION sigmoid(x,x0,dx)
! function: 0.5*(tanh((x-x0)/dx)+1)
! -------------------------------------------------------------------- !
! parameter: all real*4
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4 :: sigmoid
  real*4, intent(in) :: x, x0, dx
! ------------

! ------------
! init
!  write(*,*) " > sigmoid: INIT."
  sigmoid = 0.5*(tanh((x-x0)/dx)+1.0)
! ------------



! ------------
!  write(*,*) " > sigmoid: EXIT."
  return

END FUNCTION sigmoid
!**********************************************************************!



!**********************************************************************!
FUNCTION ERFDP(x)
  implicit none
  real*8 :: ERFDP
  real*8, intent(in) :: x
  real*8, parameter :: prec = 1.0D-6 ! precision level
  real*8 :: pi ! 3.14159265359
  real*8 :: erf ! erf(x)
  integer*4, parameter :: maxit = 10000
  integer*4 :: k ! iterator
  real*8 :: sum, fac, exp, eps ! temp vars
  !
  pi = 4.D0*datan(1.0D0)
  fac = 2.D0*x/dsqrt(pi)
  eps = dabs(prec / fac)
  !
  if (dabs(x)<5.0D0) then ! series calculation
    sum = 1.0D0
    exp = 1.0D0
    k = 0
    iterate: do
      k = k + 1
      exp = -real(2*k-1,kind=8)/real(2*k+1,kind=8)/real(k,kind=8)*x*x*exp
      ! convergence check
      if (dabs(exp)<eps) then
        ! converged
        erf = fac * sum
        exit iterate
      else if (k>=maxit) then
        ! not converged
        erf = -1.0D0
        exit iterate
      else
        ! keep on integrating
        sum = sum + exp
      end if
    end do iterate
  else
    erf = 1.0D0
    if (x<0.0D0) erf = -1.0D0
  end if
  !
  ERFDP = erf
  return
END FUNCTION ERFDP
!**********************************************************************!


!**********************************************************************!
FUNCTION ERFSP(x)
  implicit none
  real*4 :: ERFSP
  real*8, external :: ERFDP
  real*4, intent(in) :: x
  real*8 :: z
  z = real( x, kind=8 )
  ERFSP = real( ERFDP(z) , kind=4 )
  return
END FUNCTION ERFSP
!**********************************************************************!

































!**********************************************************************!
!************************* STRING FUNCTIONS ***************************!
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE SetVarString(carray,n,string)
! function: copy data from string to integer*1 array
! -------------------------------------------------------------------- !
! parameter: carray - integer*1(n) array
!            n      - array size
!            string - character(len=*)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character(len=*) :: string
  integer*1, dimension(1:n) :: carray
  integer*4 :: i, n, slen, alen, mlen
! ------------

! ------------
! init
!  write(*,*) " > SetVarString: INIT."
  alen = size(carray,dim=1)
  slen = len(string)
  mlen = min(alen,slen)
  if (mlen<=0) return
! ------------

! ------------
! copy
  do i=1, mlen
    carray(i)=mod(ichar(string(i:i)),256)
  end do
  if (alen>mlen) then
    do i=1+mlen, alen
      carray(i) = 0
    end do
  end if
! ------------

! ------------
!  write(*,*) " > SetVarString: EXIT."
  return

END SUBROUTINE SetVarString
!**********************************************************************!






!**********************************************************************!
!**********************************************************************!
SUBROUTINE GetVarString(string,carray,n)
! function: copy n data from character array to string
! -------------------------------------------------------------------- !
! parameter: carray - integer*1(n) array
!            n      - size of carray
!            string - character(len=*)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character(len=n) :: string
  integer*1, dimension(1:n) :: carray
  integer*4 :: n, i, slen, alen, mlen
! ------------

! ------------
! init
!  write(*,*) " > GetVarString: INIT."
  alen = size(carray,dim=1)
  slen = len(string)
  mlen = min(alen,slen)
  if (mlen<=0) return
! ------------

! ------------
! copy
  do i=1, mlen
    string(i:i) = achar(carray(i))
  end do
! ------------

! ------------
!  write(*,*) " > GetVarString: EXIT."
  return

END SUBROUTINE GetVarString
!**********************************************************************!







































!**********************************************************************!
!******************** routine comment structure ***********************!
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
!SUBROUTINE <NAME>(<params>)
! function: 
! -------------------------------------------------------------------- !
! parameter:
!            
! -------------------------------------------------------------------- !

!  implicit none

! ------------
! declaration
! ------------



! ------------
! init
!  write(*,*) " > <NAME>: INIT."
! ------------



! ------------
!  write(*,*) " > <NAME>: EXIT."
!  return

!END SUBROUTINE <NAME>
!**********************************************************************!

