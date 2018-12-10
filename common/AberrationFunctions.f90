!*********************************************************************!
!                                                                     !
!    MODULE AberrationFunctions                                       !
!    --------------------------                                       !
!                                                                     !
!    Purpose  : Implementation of functions and memory for            !
!               handling of coherent aberrations                      !
!               MAX. order of aberrations is defined at compilation   !
!               manipulations of simulation and analysis              !
!    File     : AberrationFunctions.f90                               !
!                                                                     !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                 !
!    Version  :  1.0.0, November 08, 2006                             !
!                                                                     !
!    To Link  : BasicFuncs.f90                                        !
!                                                                     !
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
!                                                                     !
!    IMORTANT REMARKS                                                 !
!    ----------------                                                 !
!                                                                     !
!    (1) CALL AF_INIT() before usage                                  !
!    (2) CALL AF_UNINIT() before re-usage or when halting             !
!                                                                     !
!*********************************************************************!



!*********************************************************************!
!*********************************************************************!
!    MODULE DECALRATIONS                                              !
!*********************************************************************!
!*********************************************************************!

MODULE AberrationFunctions

! Use AF as acronym for public parameters funcs and subs!
    
! Global module dependencies
!   USE ...

    
  implicit none
  
  save
    
  ! declare internal data types

  ! accessibility of subroutines or functions
!  private :: AF_***
  private :: AF_ERROR                               ! error message handling
  
  public :: AF_AberrationFunction                   ! calculates aberration function from module variables
  public :: AF_EvenAberrationFunction               ! calculates even part of the aberration function (AF)
  public :: AF_OddAberrationFunction                ! calculates odd part of the AF
  public :: AF_AberrationFunctionGradient           ! calculates gradient of the AF
  public :: AF_EvenAberrationFunctionGradient       ! calculates even part of the AF
  public :: AF_OddAberrationFunctionGradient        ! calculates odd part of the AF gradient
  public :: AF_AberrationFunctionRMSError           ! calculates root mean square error of AF
  public :: AF_AberrationFunctionRMSRingError       ! calculates azimuthally averaged r.m.s. error of AF
  public :: AF_AberrationFunctionSecDeriv           ! calculates second derivative of AF
  public :: AF_IntegratedSecDeriv                   ! calculates second derivative of AF integrated over area
  public :: AF_GetEffectiveAberration               ! calculates effective aberration due to beam tilt
  
  public :: AF_GetBasis                             ! calculates basis functions of AF
  public :: AF_GetDBasis                            ! calculates basis functions of AF gradient
  public :: AF_GetDDBasis                           ! calculates basis functions of AF second derivative
  public :: AF_GetIDDBasis                          ! calculates basis functions of AF second derivative integrals
  public :: AF_GetTiltInductionBasis                ! calculates basis functions of tilt induction
  
  public :: AF_RotateField                          ! rotates aberrations according to rotation of the image
  
  
!  public :: AF_***
  public :: AF_INIT
  public :: AF_INITBINOC
  public :: AF_UNINIT
  
  public :: AF_GetAberrationIndexByName
  public :: AF_SetAberrationByName
  public :: AF_GetAberrationByName
  public :: AF_SetAberration
  public :: AF_GetAberration
  public :: AF_SetAberrationLName
  public :: AF_GetAberrationLName
  public :: AF_SetAberrationSName
  public :: AF_GetAberrationSName
  public :: AF_SetAberrationEnable
  public :: AF_GetAberrationEnable
  public :: AF_SetAberrationError
  public :: AF_GetAberrationError
  public :: AF_WriteAberrationList
  public :: AF_BackupAllAberrations
  public :: AF_SetAberrationBackup
  public :: AF_GetAberrationBackup
  

!  declare module global (but private) variables, params and arrays
!  file units interval
  integer*4, private, parameter :: AF_minunit = 81
  integer*4, private, parameter :: AF_maxunit = 90

!   length of names and internal strings
  integer*4, public, parameter :: AF_ll = 1024

!   enpty line placeholder (emptied at startup)
  character(len=AF_ll), private :: AF_el  
  save
! pi
  real*4, private :: AF_pi
  DATA AF_pi /3.1415927/

! scale degree to radian
  real*4, private :: AF_rd2r
  DATA AF_rd2r /57.2957795/

! error counter  
  integer*4, public :: AF_err_num
  DATA AF_err_num /0/
  
! aberration power threshold
  real*4, private, parameter :: AF_POWERTHRESH_WA = 1.0E-30
  
! wavelength [nm]
  real*4, public :: AF_lamb
  DATA AF_lamb /0.001969/

! ----------->
! ----------->
! IMPORTANT PARAMETER: AF_maxaberration_order
!   controls size of dynamic allocation, length of aberration lists
! ------------------O
! -----------------O|
!                  ||
!                  vv
!              CHANGE THIS PARAMETER
!              TO INCLUDE MORE HIGHER ORDER ABERRATIONS !!!
!              + CAREFULLY CHECK AF_INIT() and AF_UNINIT()
!                FOR CORRECT MEMORY ALLOCATION
! max number of aberration coefficients
  integer*4, public, parameter :: AF_maxaberration_order = 8
! length of name strings  
  integer*4, public, parameter :: AF_asnl = 6
  integer*4, public, parameter :: AF_alnl = 40
! number of short name tables/versions
  integer*4, public, parameter :: AF_snv = 3 ! wxy, Wxy, An
! to CEOS vector translation conde
  integer*4, public :: AF_ceostcode
  DATA AF_ceostcode /0/
     
! aberration tables, dynamically allocated
! size holder, dynamically determined in AF_INIT()
  integer*4, public :: AF_maxaberration
  DATA AF_maxaberration /0/
! allocatable tables
  real*4, dimension(:,:), public, allocatable :: AF_wa ! coeffictients
  real*4, dimension(:,:), public, allocatable :: AF_wabk ! coeffictients backup
  real*4, dimension(:,:), public, allocatable :: AF_waerr ! coeffictient errors
  real*4, dimension(:,:), public, allocatable :: AF_wacov ! coeffictient covariance matrix
  integer*4, dimension(:,:), public, allocatable :: AF_waidx ! coeffictients index hash
  integer*4, dimension(:), public, allocatable :: AF_waact ! aberration activation flags
  character(len=AF_snv*AF_asnl), dimension(:), public, allocatable :: AF_asn ! short names
  character(len=AF_alnl), dimension(:), public, allocatable :: AF_aln ! long names
     
! binomial coefficient tables
  integer*4, dimension(0:2*AF_maxaberration_order,0:2*AF_maxaberration_order) :: AF_BINOMIAL_TAB
  
! last error message
  character(len=AF_ll), public :: AF_LastErrorMessage
  
  

!*********************************************************************!
!*********************************************************************!
!    MODULE BODY                                                      !
!*********************************************************************!
!*********************************************************************!

  CONTAINS





























!*********************************************************************!
!*********************************************************************!
!    INFRASTRUCTURE                                                   !
!*********************************************************************!
!*********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_INIT()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 100
  integer*4, parameter :: asl = AF_snv*AF_asnl
  integer*4, parameter :: all = AF_alnl
  
  integer*4 :: m, n, l, err
  
  real*4, external :: sigmoid ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_INIT: INIT."
! ------------

! ------------
  AF_pi = atan(1.0)*4.0
  AF_rd2r = AF_pi/180.0
  AF_err_num = 0
  AF_el = REPEAT(" ",AF_ll)
! ------------

! ------------
! preset binomial coefficient table
  call AF_INITBINOC()
! ------------

! ------------
! preset aberration tables

! (1) determine length of aberration list -> AF_maxaberration
  ! count up to max. order
  l = 0
  AF_maxaberration = 0
  do m=1, AF_maxaberration_order
    do n=0,m
      if (mod(m+n,2)==1) cycle
      l = l + 1
    end do
  end do
  if (l>0) then
    AF_maxaberration = l
  else
    call AF_ERROR("AF_INIT failed: Invalid number of aberrations.",subnum+1)
    return
  end if
  
  ! allocate arrays for tables
  allocate(AF_wa(2,AF_maxaberration),AF_wabk(2,AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+2)
    return
  end if
  allocate(AF_waerr(2,AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+3)
    return
  end if
  allocate(AF_waidx(2,AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+4)
    return
  end if
  allocate(AF_asn(AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+5)
    return
  end if
  allocate(AF_aln(AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+6)
    return
  end if
  allocate(AF_waact(AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+7)
    return
  end if
   allocate(AF_wacov(2*AF_maxaberration,2*AF_maxaberration),stat=err)
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Allocation error.",subnum+8)
    return
  end if
  
  
  ! index hash up to max. order
  l = 1
  do m=1, AF_maxaberration_order
    do n=0,m
      if (mod(m+n,2)==1) cycle
      AF_waidx(1, l) = m
      AF_waidx(2, l) = n
      l = l + 1
    end do
  end do
  
  ! coefficients
  AF_wa(:,:) = 0.0
  AF_wabk(:,:) = 0.0
  AF_waerr(:,:) = 0.0
  AF_wacov(:,:) = 0.0
  AF_waact(:) = 1 ! preactivate all
  AF_asn = ""
  AF_aln = ""
  
  ! ----------->
  ! ----------->
  ! IMPORTANT -> ADJUST THE FOLLOWING LISTS IF IF AF_maxaberration_order /= 8
  ! long names
  call AF_SetAberrationLName( 1, "Image shift ")
  call AF_SetAberrationLName( 2, "Defocus ")
  call AF_SetAberrationLName( 3, "Twofold astigmatism ")
  call AF_SetAberrationLName( 4, "3rd order axial coma ")
  call AF_SetAberrationLName( 5, "Threefold astigmatism ")
  call AF_SetAberrationLName( 6, "Spherical aberration (Cs) ")
  call AF_SetAberrationLName( 7, "Star aberration ")
  call AF_SetAberrationLName( 8, "Fourfold astigmatism ")
  call AF_SetAberrationLName( 9, "5th order axial coma ")
  call AF_SetAberrationLName(10, "Three lobe aberration ")
  call AF_SetAberrationLName(11, "Fivefold astigmatism ")
  call AF_SetAberrationLName(12, "Spherical aberration (C5) ")
  call AF_SetAberrationLName(13, "6th order Star aberration ")
  call AF_SetAberrationLName(14, "Rosette aberration ")
  call AF_SetAberrationLName(15, "Sixfold astigmatism ")
  call AF_SetAberrationLName(16, "7th order axial coma " )
  call AF_SetAberrationLName(17, "Threefoldness of the 7th order ")
  call AF_SetAberrationLName(18, "Fivefoldness of the 7th order ")
  call AF_SetAberrationLName(19, "Sevenfold astigmatism ")
  call AF_SetAberrationLName(20, "Spherical aberration (C7) ")
  call AF_SetAberrationLName(21, "Twofoldness of the 8th order ")
  call AF_SetAberrationLName(22, "Fourfoldness of the 8th order " )
  call AF_SetAberrationLName(23, "Sixfoldness of the 8th order ")
  call AF_SetAberrationLName(24, "Eightfold astigmatism ")
  ! short names / identifier
  call AF_SetAberrationSName( 1, 1, "w11 ")
  call AF_SetAberrationSName( 2, 1, "w20 ")
  call AF_SetAberrationSName( 3, 1, "w22 ")
  call AF_SetAberrationSName( 4, 1, "w31 ")
  call AF_SetAberrationSName( 5, 1, "w33 ")
  call AF_SetAberrationSName( 6, 1, "w40 ")
  call AF_SetAberrationSName( 7, 1, "w42 ")
  call AF_SetAberrationSName( 8, 1, "w44 ")
  call AF_SetAberrationSName( 9, 1, "w51 ")
  call AF_SetAberrationSName(10, 1, "w53 ")
  call AF_SetAberrationSName(11, 1, "w55 ")
  call AF_SetAberrationSName(12, 1, "w60 ")
  call AF_SetAberrationSName(13, 1, "w62 ")
  call AF_SetAberrationSName(14, 1, "w64 ")
  call AF_SetAberrationSName(15, 1, "w66 ")
  call AF_SetAberrationSName(16, 1, "w71 ")
  call AF_SetAberrationSName(17, 1, "w73 ")
  call AF_SetAberrationSName(18, 1, "w75 ")
  call AF_SetAberrationSName(19, 1, "w77 ")
  call AF_SetAberrationSName(20, 1, "w80 ")
  call AF_SetAberrationSName(21, 1, "w82 ")
  call AF_SetAberrationSName(22, 1, "w84 ")
  call AF_SetAberrationSName(23, 1, "w86 ")
  call AF_SetAberrationSName(24, 1, "w88 ")
  
  call AF_SetAberrationSName( 1, 2, "a11 ")
  call AF_SetAberrationSName( 2, 2, "a20 ")
  call AF_SetAberrationSName( 3, 2, "a22 ")
  call AF_SetAberrationSName( 4, 2, "a31 ")
  call AF_SetAberrationSName( 5, 2, "a33 ")
  call AF_SetAberrationSName( 6, 2, "a40 ")
  call AF_SetAberrationSName( 7, 2, "a42 ")
  call AF_SetAberrationSName( 8, 2, "a44 ")
  call AF_SetAberrationSName( 9, 2, "a51 ")
  call AF_SetAberrationSName(10, 2, "a53 ")
  call AF_SetAberrationSName(11, 2, "a55 ")
  call AF_SetAberrationSName(12, 2, "a60 ")
  call AF_SetAberrationSName(13, 2, "a62 ")
  call AF_SetAberrationSName(14, 2, "a64 ")
  call AF_SetAberrationSName(15, 2, "a66 ")
  call AF_SetAberrationSName(16, 2, "a71 ")
  call AF_SetAberrationSName(17, 2, "a73 ")
  call AF_SetAberrationSName(18, 2, "a75 ")
  call AF_SetAberrationSName(19, 2, "a77 ")
  call AF_SetAberrationSName(20, 2, "a80 ")
  call AF_SetAberrationSName(21, 2, "a82 ")
  call AF_SetAberrationSName(22, 2, "a84 ")
  call AF_SetAberrationSName(23, 2, "a86 ")
  call AF_SetAberrationSName(24, 2, "a88 ")
  
  call AF_SetAberrationSName( 1, 3, "A0  ")
  call AF_SetAberrationSName( 2, 3, "C1  ")
  call AF_SetAberrationSName( 3, 3, "A1  ")
  call AF_SetAberrationSName( 4, 3, "B2  ")
  call AF_SetAberrationSName( 5, 3, "A2  ")
  call AF_SetAberrationSName( 6, 3, "C3  ")
  call AF_SetAberrationSName( 7, 3, "S3  ")
  call AF_SetAberrationSName( 8, 3, "A3  ")
  call AF_SetAberrationSName( 9, 3, "B4  ")
  call AF_SetAberrationSName(10, 3, "D4  ")
  call AF_SetAberrationSName(11, 3, "A4  ")
  call AF_SetAberrationSName(12, 3, "C5  ")
  call AF_SetAberrationSName(13, 3, "S5  ")
  call AF_SetAberrationSName(14, 3, "R5  ")
  call AF_SetAberrationSName(15, 3, "A5  ")
  call AF_SetAberrationSName(16, 3, "B6  ")
  call AF_SetAberrationSName(17, 3, "D6  ")
  call AF_SetAberrationSName(18, 3, "F6  ")
  call AF_SetAberrationSName(19, 3, "A6  ")
  call AF_SetAberrationSName(20, 3, "C7  ")
  call AF_SetAberrationSName(21, 3, "S7  ")
  call AF_SetAberrationSName(22, 3, "R7  ")
  call AF_SetAberrationSName(23, 3, "H7  ")
  call AF_SetAberrationSName(24, 3, "A7  ")
! ------------

! ------------
  AF_LastErrorMessage = "";
! ------------
  
! ------------
!  write(*,*) " > AF_INIT: EXIT."
  return

END SUBROUTINE AF_INIT
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_INITBINOC()
! function: initializes binomial coefficient table
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3900
  integer*4 :: n, m ! iterators
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_INITBINOC: INIT."
  if (AF_maxaberration_order<=0) return
! ------------

! ------------
! preset binomial coefficient table
  do m=0, 2*AF_maxaberration_order
    do n=0, 2*AF_maxaberration_order
      AF_BINOMIAL_TAB(n,m) = binomial(n,m)
    end do
  end do
! ------------

! ------------
!  write(*,*) " > AF_INITBINOC: EXIT."
  return

END SUBROUTINE AF_INITBINOC
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_UNINIT()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 200
  integer*4 :: err
! ------------

! ------------
! INIT
!  write(*,*) " > AF_INIT: AF_UNINIT."
  err = 0
! ------------

! ------------
! deallocate tables
  deallocate(AF_wa,AF_wabk,AF_waerr,AF_waidx,AF_asn,AF_aln, &
     & AF_waact,AF_wacov,stat=err)
  AF_maxaberration = 0
  if (err/=0) then
    call AF_ERROR("AF_INIT failed: Deallocation error.",subnum+1)
    return
  end if
! reset aberration list length
! ------------

! ------------
!  write(*,*) " > AF_UNINIT: EXIT."
  return

END SUBROUTINE AF_UNINIT
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_ERROR(sTxt,nErr)
! function: print error message
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
!            INTEGER*4 :: nErr : error code, displayed when /= 0
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  integer*4, parameter :: subnum = 300
  
  character(len=*), intent(in) :: sTxt
  integer*4, intent(in) :: nErr
  character(len=AF_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = AF_el
! ----------

! ----------
! messaging (to screen)
  write(unit=sinfo,fmt=*) "ERROR: ",trim(sTxt)," Error code:",nErr
  AF_LastErrorMessage = trim(sinfo)
!  call SE_event(trim(sinfo), SE_err)
  write(unit=6,fmt=*) "SFT_ERROR: "//trim(sinfo)
! ----------

  AF_err_num = AF_err_num + 1
  return

END SUBROUTINE AF_ERROR
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
FUNCTION AF_SetAberrationByName(sname,wax,way)
! function: set aberration coefficients for aberration identified by
!           its short name string
! -------------------------------------------------------------------- !
! parameter: character*4 :: sname
!            real*4 :: wax, way
! return value: integer*4 = 1 success
!               integer*4 = 0 error
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 400
  integer*4, parameter :: asl = AF_asnl
  
  integer*4 :: AF_SetAberrationByName
  character(len=*), intent(in) :: sname
  real*4, intent(in) :: wax, way
  
  character(len=asl) :: searchstring
  character(len=AF_snv*asl) :: sourcestring
  integer*4 :: i, j, isl
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberrationByName: INIT."
  AF_SetAberrationByName = 0 ! start with unsuccessful state
  searchstring(1:asl) = AF_el(1:asl)
  isl = min(len(sname),asl)
  searchstring(1:isl) = sname(1:isl)
! ------------

! ------------
! try to locate aberration
  do i=1, AF_maxaberration
    sourcestring = AF_asn(i)
    j = INDEX(sourcestring,trim(searchstring))
    if (j>0) then ! aberration identified
      AF_wa(1,i) = wax
      AF_wa(2,i) = way
      AF_SetAberrationByName = 1
      exit
    end if
  end do
! ------------

! ------------
!  write(*,*) " > AF_SetAberrationByName: EXIT."
  return

END FUNCTION AF_SetAberrationByName
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION AF_GetAberrationByName(sname,wax,way)
! function: Get aberration coefficients for aberration identified by
!           its short name string
! -------------------------------------------------------------------- !
! parameter: character*4 :: sname
!            real*4 :: wax, way
! return value: integer*4 = 1 success
!               integer*4 = 0 error
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 500
  integer*4, parameter :: asl = AF_asnl
  
  integer*4 :: AF_GetAberrationByName
  character(len=*), intent(in) :: sname
  real*4, intent(out) :: wax, way
  
  character(len=AF_snv*asl) :: sourcestring
  character(len=asl) :: searchstring
  integer*4 :: i, j, isl
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationByName: INIT."
  AF_GetAberrationByName = 0 ! start with unsuccessful state
  searchstring(1:asl) = AF_el(1:asl)
  isl = min(len(sname),asl)
  searchstring(1:isl) = sname(1:isl)
! ------------

! ------------
! try to locate aberration
  do i=1, AF_maxaberration
    sourcestring = AF_asn(i)
    j = INDEX(sourcestring,trim(searchstring))
    if (j>0) then ! aberration identified
      wax = AF_wa(1,i)
      way = AF_wa(2,i)
      AF_GetAberrationByName = 1
      exit
    end if
  end do
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationByName: EXIT."
  return

END FUNCTION AF_GetAberrationByName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_SetAberration(nidx,wax,way)
! function: set aberration coefficients for aberration identified by
!           its table index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx
!            real*4 :: wax, way
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 600
  
  integer*4, intent(in) :: nidx
  real*4, intent(in) :: wax, way
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberration: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  AF_wa(1,nidx) = wax
  AF_wa(2,nidx) = way
! ------------

! ------------
!  write(*,*) " > AF_SetAberration: EXIT."
  return

END SUBROUTINE AF_SetAberration
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetAberration(nidx,wax,way)
! function: Get aberration coefficients for aberration identified by
!           its table index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx
!            real*4 :: wax, way
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 700
  
  integer*4, intent(in) :: nidx
  real*4, intent(out) :: wax, way
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberration: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  wax = AF_wa(1,nidx)
  way = AF_wa(2,nidx)
! ------------

! ------------
!  write(*,*) " > AF_GetAberration: EXIT."
  return

END SUBROUTINE AF_GetAberration
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_BackupAllAberrations
! function: Save current aberrations to backup list
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4600
! ------------

! ------------
! INIT
!  write(*,*) " > AF_BackupAllAberrations: INIT."
  if (AF_maxaberration<=1) return ! bounce
! ------------

! ------------
  AF_wabk = AF_wa
! ------------

! ------------
!  write(*,*) " > AF_BackupAllAberrations: EXIT."
  return

END SUBROUTINE AF_BackupAllAberrations
!**********************************************************************!

  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_SetAberrationBackup(nidx,wax,way)
! function: set aberration coefficients for aberration identified by
!           its table index to backup list
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx
!            real*4 :: wax, way
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4500
  
  integer*4, intent(in) :: nidx
  real*4, intent(in) :: wax, way
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberrationBackup: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  AF_wabk(1,nidx) = wax
  AF_wabk(2,nidx) = way
! ------------

! ------------
!  write(*,*) " > AF_SetAberrationBackup: EXIT."
  return

END SUBROUTINE AF_SetAberrationBackup
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetAberrationBackup(nidx,wax,way)
! function: Get aberration coefficients for aberration identified by
!           its table index from backup list
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx
!            real*4 :: wax, way
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4400
  
  integer*4, intent(in) :: nidx
  real*4, intent(out) :: wax, way
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationBackup: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  wax = AF_wabk(1,nidx)
  way = AF_wabk(2,nidx)
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationBackup: EXIT."
  return

END SUBROUTINE AF_GetAberrationBackup
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_SetAberrationLName(nIdx,sLongName)
! function: sets data to long name array
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            character(len=*) :: sLongName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1800
  integer*4, parameter :: all = AF_alnl

  integer*4, intent(in) :: nIdx
  character(len=*), intent(in) :: sLongName
  
  integer*4 :: nLen
!  external :: SetVarString
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberrationLName: INIT."
  nLen = min(all,len(sLongName))
! ------------

! ------------
  AF_aln(nIdx) = sLongName(1:nLen)
! ------------

! ------------
!  write(*,*) " > AF_SetAberrationLName: EXIT."
  return

END SUBROUTINE AF_SetAberrationLName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetAberrationLName(nIdx,sLongName)
! function: retrieves name data
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            character(len=*) :: sLongName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1900
  integer*4, parameter :: all = AF_alnl

  integer*4, intent(in) :: nIdx
  character(len=*), intent(inout) :: sLongName
  
  integer*4 :: nLen
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationLName: INIT."
  if (.not.allocated(AF_aln)) return
  nLen = min(all,len(sLongName))
! ------------

! ------------
  write(unit=sLongName,fmt='(A)') trim(AF_aln(nIdx))
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationLName: EXIT."
  return

END SUBROUTINE AF_GetAberrationLName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_SetAberrationSName(nIdx,nSub,sShortName)
! function: sets data to long name array
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            integer*4 :: nSub : sub table index
!            character(len=*) :: sShortName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2000
  integer*4, parameter :: asl = AF_asnl

  integer*4, intent(in) :: nIdx, nSub
  character(len=*), intent(in) :: sShortName
  
  integer*4 :: nLen, nOff, nEnd
  character(len=AF_snv*asl) :: currentnames
!  external :: SetVarString
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberrationSName: INIT."
  if (nIdx<1.or.nIdx>AF_maxaberration.or.nSub<1.or.nSub>AF_snv) return
  nLen = min(asl,len(sShortName))
  nOff = 1+(nSub-1)*asl
  nEnd = nOff+asl-1
! ------------

! ------------
  currentnames = AF_asn(nIdx)
  currentnames(nOff:nEnd) = trim(sShortName(1:nLen))
  AF_asn(nIdx) = currentnames
! ------------

! ------------
!  write(*,*) " > AF_SetAberrationSName: EXIT."
  return

END SUBROUTINE AF_SetAberrationSName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetAberrationSName(nIdx,nSub,sShortName)
! function: retrieves short name data
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            integer*4 :: nSub : sub table index
!            character(len=*) :: sShortName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2100
  integer*4, parameter :: asl = AF_asnl

  integer*4, intent(in) :: nIdx, nSub
  character(len=*), intent(inout) :: sShortName
  
  integer*4 :: nLen, nOff, nEnd
  character(len=AF_snv*asl) :: currentnames
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationSName: INIT."
  if (nIdx<1.or.nIdx>AF_maxaberration.or.nSub<1.or.nSub>AF_snv) return
  nLen = min(asl,len(sShortName))
  nOff = 1+(nSub-1)*asl
  nEnd = nOff+nLen-1
! ------------

! ------------
  currentnames = AF_asn(nIdx)
  write(unit=sShortName,fmt='(A)') trim(currentnames(nOff:nEnd))
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationSName: EXIT."
  return

END SUBROUTINE AF_GetAberrationSName
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_SetAberrationEnable(nidx,nEnable)
! function: set activation flag for aberration identified by
!           its table index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx, nEnable
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2200
  
  integer*4, intent(in) :: nidx, nEnable
  integer*4 :: tmp
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberrationEnable: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
  tmp = 0
  if (nEnable /= 0) tmp = 1
! ------------

! ------------
! try to locate aberration
  AF_waact(nidx) = tmp
! ------------

! ------------
!  write(*,*) " > AF_SetAberrationEnable: EXIT."
  return

END SUBROUTINE AF_SetAberrationEnable
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetAberrationEnable(nidx,nEnable)
! function: Get activation flag for aberration identified by
!           its table index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx, nEnable
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2300
  
  integer*4, intent(in) :: nidx
  real*4, intent(inout) :: nEnable
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationEnable: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  nEnable = AF_waact(nidx)
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationEnable: EXIT."
  return

END SUBROUTINE AF_GetAberrationEnable
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION AF_GetAberrationIndexByName(sname)
! function: Get aberration list index for aberration identified by
!           its short name string
! -------------------------------------------------------------------- !
! parameter: character*4 :: sname
! return value: integer*4 > 0 index, success
!               integer*4 <= 0 error
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2400
  integer*4, parameter :: asl = AF_asnl
  
  integer*4 :: AF_GetAberrationIndexByName
  character(len=*), intent(in) :: sname
  
  character(len=AF_snv*asl) :: sourcestring
  character(len=asl) :: searchstring
  integer*4 :: i, j, isl
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationIndexByName: INIT."
  AF_GetAberrationIndexByName = 0 ! start with unsuccessful state
  searchstring(1:asl) = AF_el(1:asl)
  isl = min(len(sname),asl)
  searchstring(1:isl) = sname(1:isl)
! ------------

! ------------
! try to locate aberration
  do i=1, AF_maxaberration
    sourcestring = AF_asn(i)
    j = INDEX(sourcestring,trim(searchstring))
    if (j>0) then ! aberration identified
      AF_GetAberrationIndexByName = i
      exit
    end if
  end do
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationIndexByName: EXIT."
  return

END FUNCTION AF_GetAberrationIndexByName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_SetAberrationError(nidx,waxe,waye)
! function: set aberration coefficient errors for aberration identified by
!           its table index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx
!            real*4 :: waxe, waye
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2500
  
  integer*4, intent(in) :: nidx
  real*4, intent(in) :: waxe, waye
! ------------

! ------------
! INIT
!  write(*,*) " > AF_SetAberrationError: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  AF_waerr(1,nidx) = waxe
  AF_waerr(2,nidx) = waye
! ------------

! ------------
!  write(*,*) " > AF_SetAberrationError: EXIT."
  return

END SUBROUTINE AF_SetAberrationError
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetAberrationError(nidx,waxe,waye)
! function: Get aberration coefficient errors for aberration identified by
!           its table index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nidx
!            real*4 :: waxe, waye
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2600
  
  integer*4, intent(in) :: nidx
  real*4, intent(inout) :: waxe, waye
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetAberrationError: INIT."
  if ((nidx<1).or.(nidx>AF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  waxe = AF_waerr(1,nidx)
  waye = AF_waerr(2,nidx)
! ------------

! ------------
!  write(*,*) " > AF_GetAberrationError: EXIT."
  return

END SUBROUTINE AF_GetAberrationError
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_TranslateAberrations(noti,tval,terr,tcode,nerr)
! function: translates aberration values and errors to different
!           notation
! -------------------------------------------------------------------- !
! parameter:
! in:
!   integer*4 :: noti                   = aberration notation
!                                         0: internal
!                                         1: CEOS
!   integer*4 :: tcode                  = translation codes
! in/out:
!   real*4 :: tval(2,AF_maxaberration)  = translated values
!   real*4 :: terr(2,AF_maxaberration)  = translated errors
!   integer*4 :: nerr                   = error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4200

  integer*4, intent(in) :: noti, tcode
  real*4, intent(inout) :: tval(2,AF_maxaberration)
  real*4, intent(inout) :: terr(2,AF_maxaberration)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: i                    ! iterators
  real*4 :: wax, way, waex, waey    ! placeholders
  real*4 :: scfx, scfy              ! scaling factors
! ------------

! ------------
! INIT
!  write(*,*) " > AF_TranslateAberrations: INIT."
  nerr = 0
! ------------

! ------------
  select case (noti)
  case (1) ! translate to CEOS notation
    do i=1, AF_maxaberration
      if (AF_waact(i)==0) cycle
      scfx = 1.0
      scfy = 1.0
      ! set scf to 1/m for all non-isotropic and non-full-fold aberrations
      if (AF_waidx(1,i)/=AF_waidx(2,i) .and. AF_waidx(2,i)/=0) then
        scfx = scfx / real(AF_waidx(1,i))
        scfy = scfy / real(AF_waidx(1,i))
      end if
      ! set sign flips from tcode variable
      ! bit 1 = 1 : flip Ax coefficients (A1x, A2x, A3x, ...)
      ! bit 2 = 1 : flip Ay coefficients (A1y, A2y, A3y, ...)
      ! bit 3 = 1 : flip (B,S,D,R)x coefficients (B2x, S3x, B4x, R4x, ..)
      ! bit 4 = 1 : flip (B,S,D,R)y coefficients (B2y, S3y, B4y, R4y, ..)
      if (AF_waidx(1,i)==AF_waidx(2,i)) then ! sign change in case if A-aberrations
        if (AND(tcode,1)==1) scfx = -scfx
        if (AND(tcode,2)==2) scfy = -scfy
      else
        if (AF_waidx(2,i)/=0) then ! sign change in case of non-A- and non-C-aberrations
          if (AND(tcode,4)==4) scfx = -scfx
          if (AND(tcode,8)==8) scfy = -scfy
        end if
      end if
      wax = AF_wa(1,i)
      way = AF_wa(2,i)
      waex = AF_waerr(1,i)
      waey = AF_waerr(2,i)
      tval(1,i) = scfx*wax
      tval(2,i) = scfy*way
      terr(1,i) = scfx*waex
      terr(2,i) = scfy*waey
    end do
  case default ! unknown or default cases -> no translation, just copy
    do i=1, AF_maxaberration
      if (AF_waact(i)==0) cycle
      tval(1,i) = AF_wa(1,i)
      tval(2,i) = AF_wa(2,i)
      terr(1,i) = AF_waerr(1,i)
      terr(2,i) = AF_waerr(2,i)
    end do
  end select ! case (noti)
! ------------

! ------------
!  write(*,*) " > AF_TranslateAberrations: EXIT."
  return

END SUBROUTINE AF_TranslateAberrations
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_WriteAberrationList(nu,snv,tlv,wae,prefix,suffix,nerr)
! function: writes strings to specified output unit with information
!           about current aberrations
! -------------------------------------------------------------------- !
! parameter:
! in:
!   integer*4 :: nu                 = logical file unit (0 -> stdout)
!   integer*4 :: snv                = shortname version index (1 .. AF_snv)
!   integer*4 :: tlv                = value translation flag
!                                     0: as is
!                                     1: CEOS aberrations
!   integer*4 :: wae                = aberration error flag (0: off, 1: on)
!   character(len=*) :: prefix      = line prefix string
!   character(len=*) :: suffix      = line suffix string
! in/out:
!   integer*4 :: nerr               = error code (0 -> success)
! -------------------------------------------------------------------- !

  implicit none

! DECLARATION
  integer*4, parameter :: subnum = 4100
  
  integer*4, intent(in) :: nu, snv, tlv, wae
  character(len=*), intent(in) :: prefix, suffix
  integer*4, intent(inout) :: nerr
  
  logical :: isopen                 ! open flag
  integer*4 :: i                    ! iterators
  integer*4 :: lfu                  ! internal variable for logical file unit
  integer*4 :: snvi, tlvi, waei, isfx ! internal flags
  character(len=AF_asnl) :: asn     ! aberration short name
  character(len=4) :: sdim          ! dimension string
  character(len=12) :: sv1, sv2, se1, se2, svtmp ! value strings
  real*4 :: dimscl                  ! dimension scaling factor
  real*4 :: rval                    ! temp. value
  real*4, dimension(:,:), allocatable :: waval ! internal values
  real*4, dimension(:,:), allocatable :: waerr ! internal errors

! INIT
!  write(*,*) " > AF_WriteAberrationList: INIT."
  nerr = 0
  snvi = 1
  tlvi = 0
  waei = 0
  ! checks
  if (AF_maxaberration<=0) then
    nerr = 1
    call AF_ERROR("Module not initialized.",subnum+nerr)
    return
  end if
  lfu = 0
  if (nu/=0) then ! check open state of logical file unit
    INQUIRE(unit=nu,opened=isopen)
    if (isopen) then
      lfu = nu
    end if
  end if
  if (snv>=1 .and. snv<=AF_snv) snvi = snv
  if (tlv>=0 .and. tlv<=1) tlvi = tlv
  if (wae>=0 .and. wae<=1) waei = wae
  isfx = len(suffix)
  
  ! allocate internal aberration values and errors
  allocate(waval(2,AF_maxaberration),waerr(2,AF_maxaberration),stat=nerr)
  if (nerr/=0) then
    nerr = 2
    call AF_ERROR("Heap allocation failed.",subnum+nerr)
    return
  end if
  ! preset with zeroes
  waval = 0.0
  waerr = 0.0
  
  ! translate values and errors
  call AF_TranslateAberrations(tlvi,waval,waerr,AF_ceostcode,nerr)
  if (nerr/=0) then
    nerr = 3
    call AF_ERROR("Aberration translation failed.",subnum+nerr)
    return
  end if
  
  ! loop over all aberrations and do the export
  do i=1, AF_maxaberration
    if (AF_waact(i)==0) cycle ! skip inactive aberration
    call AF_GetAberrationSName(i,snvi,asn) ! get aberration short name of wanted version
    rval = sqrt(waval(1,i)**2.0+waval(2,i)**2.0) ! get aberration modulus
    sdim = "nm  " ! standard dimension
    dimscl = 1.0
    if (rval<0.1) then ! below 0.1 nm use pm
      sdim = "pm  "
      dimscl = 1000.0
    else if (rval>1000.0 .and. rval<100000.0) then ! µm
      sdim = "um  "
      dimscl = 0.001
    else if (rval>=100000.0) then ! mm
      sdim = "mm  "
      dimscl = 1.0E-6
    end if
    if (lfu>0) then ! output to logical unit
      if (wae==0) then ! no error info
        if (AF_waidx(2,i)==0) then ! round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=lfu,fmt=501) prefix//trim(asn),sv1,trim(sdim),suffix
        else ! non-round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waval(2,i)*dimscl
          sv2 = adjustr(svtmp)
          write(unit=lfu,fmt=500) prefix//trim(asn),sv1,trim(sdim),suffix, &
     &                            trim(asn),sv2,trim(sdim),suffix
        end if
      else ! with error info
        if (AF_waidx(2,i)==0) then ! round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waerr(1,i)*dimscl
          se1 = adjustr(svtmp)
          write(unit=lfu,fmt=502) prefix//trim(asn),sv1,se1,trim(sdim),suffix
        else ! non-round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waerr(1,i)*dimscl
          se1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waval(2,i)*dimscl
          sv2 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waerr(2,i)*dimscl
          se2 = adjustr(svtmp)
          write(unit=lfu,fmt=503) prefix//trim(asn),sv1,se1,trim(sdim),suffix, &
     &                            trim(asn),sv2,se2,trim(sdim),suffix
        end if
      end if
    else ! output to std. device
      if (wae==0) then ! no error info
        if (AF_waidx(2,i)==0) then ! round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=6,fmt=501) prefix//trim(asn),sv1,trim(sdim),suffix
        else ! non-round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waval(2,i)*dimscl
          sv2 = adjustr(svtmp)
          write(unit=6,fmt=500) prefix//trim(asn),sv1,trim(sdim),suffix, &
     &                          trim(asn),sv2,trim(sdim),suffix
        end if
      else ! with error info
        if (AF_waidx(2,i)==0) then ! round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waerr(1,i)*dimscl
          se1 = adjustr(svtmp)
          write(unit=6,fmt=502) prefix//trim(asn),sv1,se1,trim(sdim),suffix
        else ! non-round aberration
          write(unit=svtmp,fmt=504) waval(1,i)*dimscl
          sv1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waerr(1,i)*dimscl
          se1 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waval(2,i)*dimscl
          sv2 = adjustr(svtmp)
          write(unit=svtmp,fmt=504) waerr(2,i)*dimscl
          se2 = adjustr(svtmp)
          write(unit=6,fmt=503) prefix//trim(asn),sv1,se1,trim(sdim),suffix, &
     &                          trim(asn),sv2,se2,trim(sdim),suffix
        end if
      end if
      
    end if
  end do ! loop over all aberrations
  
  ! memory deallocation
  deallocate(waval,waerr,stat=nerr)
  if (nerr/=0) then
    nerr = 2
    return
  end if

!  write(*,*) " > AF_WriteAberrationList: EXIT."
  return

500 format(" ",A,"x ",A12,A,A,"  ",A,"y ",A12,A,A)
501 format(" ",A,"  ",A12,A,A)
502 format(" ",A,"x ",A12," (",A12,")",A,A, &
     &    "  ",A,"y ",A12," (",A12,")",A,A)
503 format(" ",A,"  ",A12," (",A12,")",A,A)
504 format(G12.4)

END SUBROUTINE AF_WriteAberrationList
!**********************************************************************!




































































!*********************************************************************!
!*********************************************************************!
!    CALCULATIONS                                                     !
!*********************************************************************!
!*********************************************************************!





! -------> COHERENT THINGS
! -------> COHERENT THINGS
! -------> COHERENT THINGS

!**********************************************************************!
!**********************************************************************!
FUNCTION AF_AberrationFunction(wx, wy)
! function: Calculates the aberration function from current aberration
!           vector at given diffraction angle (wx, wy) in the diffraction
!           plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4 :: AF_AberrationFunction

  real*4, intent(in) :: wx, wy
  
  integer*4 :: j,k,l,m,n
  real*4 :: wfield(2,0:AF_maxaberration_order)
  real*4 :: wabs(0:AF_maxaberration_order)
  real*4 :: pwx, pwy, prefac, w, pw
  real*4 :: rtmp, wax, way, ttmp, ttmp1, tsgn
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_AberrationFunction: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  rtmp = 0.0
  AF_AberrationFunction = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  pw = 1.0
  w = sqrt(wx*wx+wy*wy)
  do j = 0,AF_maxaberration_order
    if (j>0) then
      pwy = pwy * wy
      pwx = pwx * wx
      pw = pw * w
    end if
    wabs(j) = pw      ! -> w^a
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
  
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)>AF_POWERTHRESH_WA) then ! use this aberration
      ttmp = 0.0  
      ! get aberration order m and rotational symmetry n
      m = AF_waidx(1,k)
      n = AF_waidx(2,k)
      ! sum up binomial terms
      do l = 0, n
        
        ttmp1 = 0.
        ! get first term sign
        tsgn = sgnprm(mod(l,4))
        ttmp1 = ttmp1 + tsgn*wax
        ! get second term sign
        tsgn = sgnprm(mod(l+3,4))
        ttmp1 = ttmp1 + tsgn*way
     
        ! get exponents
        j = n-l
        !ttmp = ttmp + ttmp1 * binomial(n,l)*wfield(1,j)*wfield(2,l)
        ttmp = ttmp + ttmp1 * AF_BINOMIAL_TAB(n,l)*wfield(1,j)*wfield(2,l)
        
      end do
      
      j = m - n
      ttmp = ttmp * wabs(j) / real(m) ! final scaling for current term
      rtmp = rtmp + ttmp
    end if
    
  end do
  AF_AberrationFunction = rtmp * prefac ! final scaling for aberration function
! ------------


! ------------
!  write(*,*) " > AF_AberrationFunction: EXIT."
  return

END FUNCTION AF_AberrationFunction
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION AF_EvenAberrationFunction(wx, wy)
! function: Calculates the even part of the aberration function
!           from current aberration vector at given diffraction
!           angle (wx, wy) in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 900
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4 :: AF_EvenAberrationFunction

  real*4, intent(in) :: wx, wy
  
  integer*4 :: j,k,l,m,n
  real*4 :: wfield(2,0:AF_maxaberration_order)
  real*4 :: wabs(0:AF_maxaberration_order)
  real*4 :: pwx, pwy, prefac, w, pw
  real*4 :: rtmp, wax, way, ttmp, ttmp1, tsgn
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_EvenAberrationFunction: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  rtmp = 0.0
  AF_EvenAberrationFunction = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  pw = 1.0
  w = sqrt(wx*wx+wy*wy)
  do j = 0,AF_maxaberration_order
    if (j>0) then
      pwy = pwy * wy
      pwx = pwx * wx
      pw = pw * w
    end if
    wabs(j) = pw      ! -> w^a
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
    
    m = AF_waidx(1,k)
    if (modulo(m,2)==1) cycle ! skip odd aberrations
    
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)>AF_POWERTHRESH_WA) then ! use this aberration
      ttmp = 0.0  
      ! get aberration order m and rotational symmetry n
      n = AF_waidx(2,k)
      ! sum up binomial terms
      do l = 0, n
        
        ttmp1 = 0.
        ! get first term sign
        tsgn = sgnprm(mod(l,4))
        ttmp1 = ttmp1 + tsgn*wax
        ! get second term sign
        tsgn = sgnprm(mod(l+3,4))
        ttmp1 = ttmp1 + tsgn*way
     
        ! get exponents
        j = n-l
        !ttmp = ttmp + ttmp1 * binomial(n,l)*wfield(1,j)*wfield(2,l)
        ttmp = ttmp + ttmp1 * AF_BINOMIAL_TAB(n,l)*wfield(1,j)*wfield(2,l)
        
      end do
      
      j = m - n
      ttmp = ttmp * wabs(j) / real(m) ! final scaling for current term
      rtmp = rtmp + ttmp
    end if
    
  end do
  AF_EvenAberrationFunction = rtmp * prefac ! final scaling for aberration function
! ------------


! ------------
!  write(*,*) " > AF_EvenAberrationFunction: EXIT."
  return

END FUNCTION AF_EvenAberrationFunction
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
FUNCTION AF_OddAberrationFunction(wx, wy)
! function: Calculates the odd part of the aberration function
!           from current aberration vector at given diffraction
!           angle (wx, wy) in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1000
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4 :: AF_OddAberrationFunction

  real*4, intent(in) :: wx, wy
  
  integer*4 :: j,k,l,m,n
  real*4 :: wfield(2,0:AF_maxaberration_order)
  real*4 :: wabs(0:AF_maxaberration_order)
  real*4 :: pwx, pwy, prefac, w, pw
  real*4 :: rtmp, wax, way, ttmp, ttmp1, tsgn
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_OddAberrationFunction: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  rtmp = 0.0
  AF_OddAberrationFunction = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  pw = 1.0
  w = sqrt(wx*wx+wy*wy)
  do j = 0,AF_maxaberration_order
    if (j>0) then
      pwy = pwy * wy
      pwx = pwx * wx
      pw = pw * w
    end if
    wabs(j) = pw      ! -> w^a
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  do k = 1, AF_maxaberration ! sum over list!

    if (AF_waact(k)==0) cycle ! skip inactive aberration
        
    m = AF_waidx(1,k)
    if (modulo(m,2)==0) cycle ! skip even aberrations
    
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)>AF_POWERTHRESH_WA) then ! use this aberration
      ttmp = 0.0  
      ! get aberration order m and rotational symmetry n
      n = AF_waidx(2,k)
      ! sum up binomial terms
      do l = 0, n
        
        ttmp1 = 0.
        ! get first term sign
        tsgn = sgnprm(mod(l,4))
        ttmp1 = ttmp1 + tsgn*wax
        ! get second term sign
        tsgn = sgnprm(mod(l+3,4))
        ttmp1 = ttmp1 + tsgn*way
     
        ! get exponents
        j = n-l
        !ttmp = ttmp + ttmp1 * binomial(n,l)*wfield(1,j)*wfield(2,l)
        ttmp = ttmp + ttmp1 * AF_BINOMIAL_TAB(n,l)*wfield(1,j)*wfield(2,l)
        
      end do
      
      j = m - n
      ttmp = ttmp * wabs(j) / real(m) ! final scaling for current term
      rtmp = rtmp + ttmp
    end if
    
  end do
  AF_OddAberrationFunction = rtmp * prefac ! final scaling for aberration function
! ------------


! ------------
!  write(*,*) " > AF_OddAberrationFunction: EXIT."
  return

END FUNCTION AF_OddAberrationFunction
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_AberrationFunctionGradient(wx, wy, afgradx, afgrady)
! function: Calculates the aberration function gradient from current
!           aberration vector at given diffraction angle (wx, wy)
!           in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
!            real*4 :: afgradx, afgrady : gradient of the aber. func.
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3000
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: afgradx, afgrady
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-1:AF_maxaberration_order)
  real*4 :: rtmpx, rtmpy
  real*4 :: pwx, pwy, prefac, wax, way, rm
  real*4 :: betax, betay, xi
  real*4 :: ttmpxx, ttmpxy, ttmpyx, ttmpyy, ttmp1x, ttmp1y
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_AberrationFunctionGradient: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  
  afgradx = 0.0
  afgrady = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------

! ------------
! sum up relevant gradients in x- and y-direction
  rtmpx = 0.0
  rtmpy = 0.0
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
  
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)<AF_POWERTHRESH_WA) cycle ! skip, aberration to weak
        
    ! get aberration order m and rotational symmetry n
    m = AF_waidx(1,k)
    if (m<1) cycle ! skip aberrations of order less than 1, they have no 1st derivative
    rm = 1.0 / real(m)
    n = AF_waidx(2,k)
    im = m-n
    
    ttmpxx = 0.0  
    ttmpxy = 0.0
    ttmpyx = 0.0  
    ttmpyy = 0.0
  
    ! ------------  
    ! do bx terms
    ! sum up binomial terms
    do l = 0, n
    
      ttmp1x = 0.0
      ttmp1y = 0.0
    
      ! get x-term sign and factor
      betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
      ! get y-term sign and factor
      betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
      ! calculate inner sum
      do i=0, im, 2
        xi = AF_BINOMIAL_TAB(im/2,i/2)
        ttmp1x = ttmp1x + xi*(m-l-i)*wfield(1,m-l-i-1)*wfield(2,l+i)
        ttmp1y = ttmp1y + xi*(l+i)  *wfield(1,m-l-i)  *wfield(2,l+i-1)
      end do
    
      ! add to outer sum
      ttmpxx = ttmpxx + ttmp1x * betax
      ttmpxy = ttmpxy + ttmp1x * betay
      ttmpyx = ttmpyx + ttmp1y * betax
      ttmpyy = ttmpyy + ttmp1y * betay
    
    end do
    ! ------------
  
    rtmpx = rtmpx + rm * (wax * ttmpxx + way * ttmpxy)
    rtmpy = rtmpy + rm * (wax * ttmpyx + way * ttmpyy)

  end do
  
  afgradx = rtmpx * prefac ! final scaling for aberration function
  afgrady = rtmpy * prefac ! final scaling for aberration function
! ------------

! ------------
!  write(*,*) " > AF_AberrationFunctionGradient: EXIT."
  return

END SUBROUTINE AF_AberrationFunctionGradient
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_EvenAberrationFunctionGradient(wx, wy, afgradx, afgrady)
! function: Calculates the gradient of the even part of the aberration
!           function from current aberration vector at given diffraction
!           angle (wx, wy) in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
!            real*4 :: afgradx, afgrady : gradient of the aber. func.
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3100
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: afgradx, afgrady
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-1:AF_maxaberration_order)
  real*4 :: rtmpx, rtmpy
  real*4 :: pwx, pwy, prefac, wax, way
  real*4 :: betax, betay, xi
  real*4 :: rm
  real*4 :: ttmpxx, ttmpxy, ttmpyx, ttmpyy, ttmp1x, ttmp1y
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_EvenAberrationFunctionGradient: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  
  afgradx = 0.0
  afgrady = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------

! ------------
! sum up relevant gradients in x- and y-direction
  rtmpx = 0.0
  rtmpy = 0.0
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
  
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)<AF_POWERTHRESH_WA) cycle ! skip, aberration to weak
        
    ! get aberration order m and rotational symmetry n
    m = AF_waidx(1,k)
    if (modulo(m,2)==1) cycle ! skip odd aberrations
    if (m<1) cycle ! skip aberrations of order less than 1, they have no 1st derivative
    rm = 1.0 / real(m)
    n = AF_waidx(2,k)
    im = m-n
    
    ttmpxx = 0.0  
    ttmpxy = 0.0
    ttmpyx = 0.0  
    ttmpyy = 0.0
  
    ! ------------  
    ! do bx terms
    ! sum up binomial terms
    do l = 0, n
    
      ttmp1x = 0.0
      ttmp1y = 0.0
    
      ! get x-term sign and factor
      betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
      ! get y-term sign and factor
      betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
      ! calculate inner sum
      do i=0, im, 2
        xi = AF_BINOMIAL_TAB(im/2,i/2)
        ttmp1x = ttmp1x + xi*(m-l-i)*wfield(1,m-l-i-1)*wfield(2,l+i)
        ttmp1y = ttmp1y + xi*(l+i)  *wfield(1,m-l-i)  *wfield(2,l+i-1)
      end do
    
      ! add to outer sum
      ttmpxx = ttmpxx + ttmp1x * betax
      ttmpxy = ttmpxy + ttmp1x * betay
      ttmpyx = ttmpyx + ttmp1y * betax
      ttmpyy = ttmpyy + ttmp1y * betay
    
    end do
    ! ------------
  
    rtmpx = rtmpx + rm * (wax * ttmpxx + way * ttmpxy)
    rtmpy = rtmpy + rm * (wax * ttmpyx + way * ttmpyy)

  end do
  
  afgradx = rtmpx * prefac ! final scaling for aberration function
  afgrady = rtmpy * prefac ! final scaling for aberration function
! ------------

! ------------
!  write(*,*) " > AF_EvenAberrationFunctionGradient: EXIT."
  return

END SUBROUTINE AF_EvenAberrationFunctionGradient
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_OddAberrationFunctionGradient(wx, wy, afgradx, afgrady)
! function: Calculates the gradient of the odd part of the aberration
!           function from current aberration vector at given diffraction
!           angle (wx, wy) in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
!            real*4 :: afgradx, afgrady : gradient of the aber. func.
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3200
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: afgradx, afgrady
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-1:AF_maxaberration_order)
  real*4 :: rtmpx, rtmpy
  real*4 :: pwx, pwy, prefac, wax, way, rm
  real*4 :: betax, betay, xi
  real*4 :: ttmpxx, ttmpxy, ttmpyx, ttmpyy, ttmp1x, ttmp1y
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_OddAberrationFunctionGradient: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  
  afgradx = 0.0
  afgrady = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------

! ------------
! sum up relevant gradients in x- and y-direction
  rtmpx = 0.0
  rtmpy = 0.0
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
  
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)<AF_POWERTHRESH_WA) cycle ! skip, aberration to weak
        
    ! get aberration order m and rotational symmetry n
    m = AF_waidx(1,k)
    if (modulo(m,2)==0) cycle ! skip even aberrations
    if (m<1) cycle ! skip aberrations of order less than 1, they have no 1st derivative
    rm = 1.0 / real(m)
    n = AF_waidx(2,k)
    im = m-n
    
    ttmpxx = 0.0  
    ttmpxy = 0.0
    ttmpyx = 0.0  
    ttmpyy = 0.0
  
    ! ------------  
    ! do bx terms
    ! sum up binomial terms
    do l = 0, n
    
      ttmp1x = 0.0
      ttmp1y = 0.0
    
      ! get x-term sign and factor
      betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
      ! get y-term sign and factor
      betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
      ! calculate inner sum
      do i=0, im, 2
        xi = AF_BINOMIAL_TAB(im/2,i/2)
        ttmp1x = ttmp1x + xi*(m-l-i)*wfield(1,m-l-i-1)*wfield(2,l+i)
        ttmp1y = ttmp1y + xi*(l+i)  *wfield(1,m-l-i)  *wfield(2,l+i-1)
      end do
    
      ! add to outer sum
      ttmpxx = ttmpxx + ttmp1x * betax
      ttmpxy = ttmpxy + ttmp1x * betay
      ttmpyx = ttmpyx + ttmp1y * betax
      ttmpyy = ttmpyy + ttmp1y * betay
    
    end do
    ! ------------
  
    rtmpx = rtmpx + rm * (wax * ttmpxx + way * ttmpxy)
    rtmpy = rtmpy + rm * (wax * ttmpyx + way * ttmpyy)

  end do
  
  afgradx = rtmpx * prefac ! final scaling for aberration function
  afgrady = rtmpy * prefac ! final scaling for aberration function
! ------------

! ------------
!  write(*,*) " > AF_OddAberrationFunctionGradient: EXIT."
  return

END SUBROUTINE AF_OddAberrationFunctionGradient
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_AberrationFunctionSecDeriv(wx, wy, afdxx, afdxy, afdyy)
! function: Calculates the aberration function second derivatives
!           from current aberration vector at given diffraction
!           angle (wx, wy) in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
!            real*4 :: afdxx, afdxy, afdyy : second derivatives of AF
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3700
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: afdxx, afdxy, afdyy
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-2:AF_maxaberration_order)
  real*4 :: rtmpxx, rtmpxy, rtmpyy
  real*4 :: pwx, pwy, prefac, wax, way, rm
  real*4 :: betax, betay, xi
  real*4 :: ttmpxxx, ttmpxxy, ttmpxyx, ttmpxyy, ttmpyyx, ttmpyyy
  real*4 :: ttmp1xx, ttmp1xy, ttmp1yy
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_AberrationFunctionSecDeriv: INIT."
  afdxx = 0.0
  afdxy = 0.0
  afdyy = 0.0
  prefac = 2.0*AF_pi/AF_lamb
! ------------


! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-2) = 0.0
  wfield(2,-2) = 0.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------

! ------------
! sum up relevant gradients in x- and y-direction
  rtmpxx = 0.0
  rtmpxy = 0.0
  rtmpyy = 0.0
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
  
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)<AF_POWERTHRESH_WA) cycle ! skip, aberration to weak
        
    ! get aberration order m and rotational symmetry n
    m = AF_waidx(1,k)
    if (m<2) cycle ! skip aberrations of order less than 2, they have no 2nd derivatives
    rm = 1.0 / real(m)
    n = AF_waidx(2,k)
    im = m-n
    
    ttmpxxx = 0.0  
    ttmpxxy = 0.0
    ttmpxyx = 0.0  
    ttmpxyy = 0.0
    ttmpyyx = 0.0  
    ttmpyyy = 0.0
  
    ! ------------  
    ! do bx terms
    ! sum up binomial terms
    do l = 0, n
    
      ttmp1xx = 0.0
      ttmp1xy = 0.0
      ttmp1yy = 0.0
      
      ! get x-term sign and factor
      betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
      ! get y-term sign and factor
      betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
      ! calculate inner sum
      do i=0, im, 2
        xi = AF_BINOMIAL_TAB(im/2,i/2)
        ttmp1xx = ttmp1xx + xi*(m-l-i)*(m-l-i-1)*wfield(1,m-l-i-2)*wfield(2,l+i)
        ttmp1xy = ttmp1xy + xi*(m-l-i)*(l+i)    *wfield(1,m-l-i-1)*wfield(2,l+i-1)
        ttmp1yy = ttmp1yy + xi*(l+i)  *(l+i-1)  *wfield(1,m-l-i)  *wfield(2,l+i-2)
      end do
    
      ! add to outer sum
      ttmpxxx = ttmpxxx + ttmp1xx * betax
      ttmpxxy = ttmpxxy + ttmp1xx * betay
      ttmpxyx = ttmpxyx + ttmp1xy * betax
      ttmpxyy = ttmpxyy + ttmp1xy * betay
      ttmpyyx = ttmpyyx + ttmp1yy * betax
      ttmpyyy = ttmpyyy + ttmp1yy * betay
    
    end do
    ! ------------
  
    rtmpxx = rtmpxx + rm * (wax * ttmpxxx + way * ttmpxxy)
    rtmpxy = rtmpxy + rm * (wax * ttmpxyx + way * ttmpxyy)
    rtmpyy = rtmpyy + rm * (wax * ttmpyyx + way * ttmpyyy)

  end do
  
  afdxx = rtmpxx * prefac ! final scaling for aberration function
  afdxy = rtmpxy * prefac ! final scaling for aberration function
  afdyy = rtmpyy * prefac ! final scaling for aberration function
! ------------

! ------------
!  write(*,*) " > AF_AberrationFunctionSecDeriv: EXIT."
  return

END SUBROUTINE AF_AberrationFunctionSecDeriv
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_IntegratedSecDeriv(wx, wy, afdxx, afdxy, afdyy)
! function: Calculates the area integral of the aberration function
!           second derivatives from current aberration vector at
!           diffraction angle (wx, wy).
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
!            real*4 :: afdxx, afdxy, afdyy : integrated 2nd derivatives
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4500
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: afdxx, afdxy, afdyy
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-2:AF_maxaberration_order)
  real*4 :: rtmpxx, rtmpxy, rtmpyy
  real*4 :: pwx, pwy, prefac, wax, way, rm
  real*4 :: betax, betay, xi
  real*4 :: ttmpxxx, ttmpxxy, ttmpxyx, ttmpxyy, ttmpyyx, ttmpyyy
  real*4 :: ttmp1xx, ttmp1xy, ttmp1yy
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_IntegratedSecDeriv: INIT."
  afdxx = 0.0
  afdxy = 0.0
  afdyy = 0.0
  prefac = 2.0*AF_pi/AF_lamb
! ------------


! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-2) = 0.0
  wfield(2,-2) = 0.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------

! ------------
! sum up relevant gradients in x- and y-direction
  rtmpxx = 0.0
  rtmpxy = 0.0
  rtmpyy = 0.0
  do k = 1, AF_maxaberration ! sum over list!
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
  
    ! get power of aberration
    wax = AF_wa(1,k)
    way = AF_wa(2,k)
    
    if ((wax*wax+way*way)<AF_POWERTHRESH_WA) cycle ! skip, aberration to weak
        
    ! get aberration order m and rotational symmetry n
    m = AF_waidx(1,k)
    if (m<2) cycle ! skip aberrations of order less than 2, they have no 2nd derivatives
    rm = 1.0 / real(m)
    n = AF_waidx(2,k)
    im = m-n
    
    ttmpxxx = 0.0  
    ttmpxxy = 0.0
    ttmpxyx = 0.0  
    ttmpxyy = 0.0
    ttmpyyx = 0.0  
    ttmpyyy = 0.0
  
    ! ------------  
    ! do bx terms
    ! sum up binomial terms
    do l = 0, n
    
      ttmp1xx = 0.0
      ttmp1xy = 0.0
      ttmp1yy = 0.0
      
      ! get x-term sign and factor
      betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
      ! get y-term sign and factor
      betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
      ! calculate inner sum
      do i=0, im, 2
        xi = AF_BINOMIAL_TAB(im/2,i/2)
        ttmp1xx = ttmp1xx + xi*(m-l-i)/(l+i+1)  *wfield(1,m-l-i-1)*wfield(2,l+i+1)
        ttmp1xy = ttmp1xy + xi                  *wfield(1,m-l-i)  *wfield(2,l+i)
        ttmp1yy = ttmp1yy + xi*(l+i)  /(m-l-i+1)*wfield(1,m-l-i+1)*wfield(2,l+i-1)
      end do
    
      ! add to outer sum
      ttmpxxx = ttmpxxx + ttmp1xx * betax
      ttmpxxy = ttmpxxy + ttmp1xx * betay
      ttmpxyx = ttmpxyx + ttmp1xy * betax
      ttmpxyy = ttmpxyy + ttmp1xy * betay
      ttmpyyx = ttmpyyx + ttmp1yy * betax
      ttmpyyy = ttmpyyy + ttmp1yy * betay
    
    end do
    ! ------------
  
    rtmpxx = rtmpxx + rm * (wax * ttmpxxx + way * ttmpxxy)
    rtmpxy = rtmpxy + rm * (wax * ttmpxyx + way * ttmpxyy)
    rtmpyy = rtmpyy + rm * (wax * ttmpyyx + way * ttmpyyy)

  end do
  
  afdxx = rtmpxx * prefac ! final scaling for aberration function
  afdxy = rtmpxy * prefac ! final scaling for aberration function
  afdyy = rtmpyy * prefac ! final scaling for aberration function
! ------------

! ------------
!  write(*,*) " > AF_IntegratedSecDeriv: EXIT."
  return

END SUBROUTINE AF_IntegratedSecDeriv
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetBasis(nIdx,wx,wy,bx,by)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1200
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)
  
  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: bx, by
  integer*4, intent(in) :: nIdx
  
  integer*4 :: j,k,l,m,n
  real*4 :: wfield(2,0:AF_maxaberration_order)
  real*4 :: wabs(0:AF_maxaberration_order)
  real*4 :: pwx, pwy, prefac, w, pw
  real*4 :: ttmpx, ttmpy, ttmp1x, ttmp1y
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetBasis: INIT."
  bx = 0.0
  by = 0.0
  prefac = 2.0*AF_pi/AF_lamb
! ------------


! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  pw = 1.0
  w = sqrt(wx*wx+wy*wy)
  do j = 0,AF_maxaberration_order
    if (j>0) then
      pwy = pwy * wy
      pwx = pwx * wx
      pw = pw * w
    end if
    wabs(j) = pw      ! -> w^a
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  k = nIdx
  ttmpx = 0.0  
  ttmpy = 0.0
  ! get aberration order m and rotational symmetry n
  m = AF_waidx(1,k)
  n = AF_waidx(2,k)
  ! sum up binomial terms
  do l = 0, n
    ! get x-term sign
    ttmp1x = sgnprm(mod(l,4))
    ! get y-term term sign
    ttmp1y = sgnprm(mod(l+3,4))
    ! get exponents
    j = n-l
    ttmpx = ttmpx + ttmp1x * AF_BINOMIAL_TAB(n,l)*wfield(1,j)*wfield(2,l)
    ttmpy = ttmpy + ttmp1y * AF_BINOMIAL_TAB(n,l)*wfield(1,j)*wfield(2,l)
  end do
      
  j = m - n

  ttmpx = ttmpx * wabs(j) / real(m) ! final scaling for current term
  ttmpy = ttmpy * wabs(j) / real(m) 
  
  bx = ttmpx * prefac ! final scaling for aberration function
  by = ttmpy * prefac
! ------------

! ------------
!  write(*,*) " > AF_GetBasis: EXIT."
  return

END SUBROUTINE AF_GetBasis
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetDBasis(nIdx,wx,wy,bxx,bxy,byx,byy)
! function: calculates basis values of the first derivative
!           of the aberration function for a given aberration by index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx    = aberration index
!            real*4 :: wx, wy     = diffraction angle [rad]
!            real*4 :: bxx, bxy   = basis values for x derivative
!                                   (bxx for Wx and bxy for Wy)
!            real*4 :: byx, byy   = basis values for y derivative
!                                   (byx for Wx and byy for Wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3600
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  integer*4, intent(in) :: nIdx
  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: bxx, bxy, byx, byy
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-1:AF_maxaberration_order)
  real*4 :: betax, betay, xi
  real*4 :: pwx, pwy, prefac
  real*4 :: ttmpxx, ttmpxy, ttmpyx, ttmpyy, ttmp1x, ttmp1y
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetDBasis: INIT."
  bxx = 0.0
  bxy = 0.0
  byx = 0.0
  byy = 0.0
  prefac = 2.0*AF_pi/AF_lamb
! ------------


! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  k = nIdx
  ttmpxx = 0.0  
  ttmpxy = 0.0
  ttmpyx = 0.0  
  ttmpyy = 0.0
  ! get aberration order m and rotational symmetry n
  m = AF_waidx(1,k)
  if (m<1) return ! skip aberrations of order less than 1, they have no 1st derivative
  prefac = prefac / real(m)
  n = AF_waidx(2,k)
  im = m-n
  
  ! ------------  
  ! sum up binomial terms
  do l = 0, n
    
    ttmp1x = 0.0
    ttmp1y = 0.0
    
    ! get x-term sign and factor
    betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
    ! get y-term sign and factor
    betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
    ! calculate inner sum
    do i=0, im, 2
      xi = AF_BINOMIAL_TAB(im/2,i)
      ttmp1x = ttmp1x + xi*(m-l-i)*wfield(1,m-l-i-1)*wfield(2,l+i)
      ttmp1y = ttmp1y + xi*(l+i)  *wfield(1,m-l-i)  *wfield(2,l+i-1)
    end do
    
    ! add to outer sum
    ttmpxx = ttmpxx + ttmp1x * betax
    ttmpxy = ttmpxy + ttmp1x * betay
    ttmpyx = ttmpyx + ttmp1y * betax
    ttmpyy = ttmpyy + ttmp1y * betay
    
  end do
  ! ------------
  
  bxx = ttmpxx * prefac ! final scaling for aberration function
  bxy = ttmpxy * prefac
  byx = ttmpyx * prefac
  byy = ttmpyy * prefac
! ------------

! ------------
!  write(*,*) " > AF_GetDBasis: EXIT."
  return

END SUBROUTINE AF_GetDBasis
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetDDBasis(nIdx,wx,wy,bxxx,bxxy,bxyx,bxyy,byyx,byyy)
! function: calculates basis values of the second derivative
!           of the aberration function for a given aberration by index
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx    = aberration index
!            real*4 :: wx, wy     = diffraction angle [rad]
!            real*4 :: bxxx, bxxy = basis values for double x derivative
!                                   (bxxx for Wx and bxxy for Wy)
!            real*4 :: bxyx, bxyy = basis values for x & y derivative
!                                   (bxyx for Wx and bxyy for Wy)
!            real*4 :: bxyy, byyy = basis values for double y derivative
!                                   (byyx for Wx and byyy for Wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3500
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  integer*4, intent(in) :: nIdx
  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: bxxx, bxxy, bxyx, bxyy, byyx, byyy
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-2:AF_maxaberration_order)
  real*4 :: betax, betay, xi
  real*4 :: pwx, pwy, prefac
  real*4 :: ttmpxxx, ttmpxxy, ttmpxyx, ttmpxyy, ttmpyyx, ttmpyyy
  real*4 :: ttmp1xx, ttmp1xy, ttmp1yy
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetDDBasis: INIT."
  bxxx = 0.0
  bxxy = 0.0
  bxyx = 0.0
  bxyy = 0.0
  byyx = 0.0
  byyy = 0.0
  prefac = 2.0*AF_pi/AF_lamb
! ------------


! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-2) = 0.0
  wfield(2,-2) = 0.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  k = nIdx
  ttmpxxx = 0.0  
  ttmpxxy = 0.0
  ttmpxyx = 0.0  
  ttmpxyy = 0.0
  ttmpyyx = 0.0  
  ttmpyyy = 0.0
  ! get aberration order m and rotational symmetry n
  m = AF_waidx(1,k)
  if (m<2) return ! skip aberrations of order less than 2, they have no 2nd derivative
  prefac = prefac / real(m)
  n = AF_waidx(2,k)
  im = m-n
  
  ! ------------  
  ! sum up binomial terms
  do l = 0, n
    
    ttmp1xx = 0.0
    ttmp1xy = 0.0
    ttmp1yy = 0.0
    
    ! get x-term sign and factor
    betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
    ! get y-term sign and factor
    betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
    ! calculate inner sum
    do i=0, im, 2
      xi = AF_BINOMIAL_TAB(im/2,i/2)
      ttmp1xx = ttmp1xx + xi*(m-l-i)*(m-l-i-1)*wfield(1,m-l-i-2)*wfield(2,l+i)
      ttmp1xy = ttmp1xy + xi*(m-l-i)*(l+i)    *wfield(1,m-l-i-1)*wfield(2,l+i-1)
      ttmp1yy = ttmp1yy + xi*(l+i)  *(l+i-1)  *wfield(1,m-l-i)  *wfield(2,l+i-2)
    end do
    
    ! add to outer sum
    ttmpxxx = ttmpxxx + ttmp1xx * betax
    ttmpxxy = ttmpxxy + ttmp1xx * betay
    ttmpxyx = ttmpxyx + ttmp1xy * betax
    ttmpxyy = ttmpxyy + ttmp1xy * betay
    ttmpyyx = ttmpyyx + ttmp1yy * betax
    ttmpyyy = ttmpyyy + ttmp1yy * betay
    
  end do
  ! ------------
  
  bxxx = ttmpxxx * prefac ! final scaling for aberration function
  bxxy = ttmpxxy * prefac
  bxyx = ttmpxyx * prefac
  bxyy = ttmpxyy * prefac
  byyx = ttmpyyx * prefac
  byyy = ttmpyyy * prefac
! ------------

! ------------
!  write(*,*) " > AF_GetDDBasis: EXIT."
  return

END SUBROUTINE AF_GetDDBasis
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetIDDBasis(nIdx,wx,wy,bxxx,bxxy,bxyx,bxyy,byyx,byyy)
! function: calculates basis values of the antiderivatives of the
!           second derivatives of the aberration function for a given
!           aberration by index
!           the antiderivatives are integrals of the 2nd derivatives
!           over x and y
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx    = aberration index
!            real*4 :: wx, wy     = diffraction angle [rad]
!            real*4 :: bxxx, bxxy = basis values for double x derivative
!                                   (bxxx for Wx and bxxy for Wy)
!            real*4 :: bxyx, bxyy = basis values for x & y derivative
!                                   (bxyx for Wx and bxyy for Wy)
!            real*4 :: bxyy, byyy = basis values for double y derivative
!                                   (byyx for Wx and byyy for Wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4400
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  integer*4, intent(in) :: nIdx
  real*4, intent(in) :: wx, wy
  real*4, intent(inout) :: bxxx, bxxy, bxyx, bxyy, byyx, byyy
  
  integer*4 :: i,j,k,l,m,n,im
  real*4 :: wfield(2,-2:AF_maxaberration_order)
  real*4 :: betax, betay, xi
  real*4 :: pwx, pwy, prefac
  real*4 :: ttmpxxx, ttmpxxy, ttmpxyx, ttmpxyy, ttmpyyx, ttmpyyy
  real*4 :: ttmp1xx, ttmp1xy, ttmp1yy
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetIDDBasis: INIT."
  bxxx = 0.0
  bxxy = 0.0
  bxyx = 0.0
  bxyy = 0.0
  byyx = 0.0
  byyy = 0.0
  prefac = 2.0*AF_pi/AF_lamb
! ------------


! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  wfield(1,-2) = 0.0
  wfield(2,-2) = 0.0
  wfield(1,-1) = 0.0
  wfield(2,-1) = 0.0
  wfield(1,0) = 1.0
  wfield(2,0) = 1.0
  do j = 1,AF_maxaberration_order
    pwy = pwy * wy
    pwx = pwx * wx
    wfield(1,j) = pwx ! -> wx^a
    wfield(2,j) = pwy ! -> wy^a
  end do
! ------------


! ------------
! sum up relevant aberrations
  k = nIdx
  ttmpxxx = 0.0  
  ttmpxxy = 0.0
  ttmpxyx = 0.0  
  ttmpxyy = 0.0
  ttmpyyx = 0.0  
  ttmpyyy = 0.0
  ! get aberration order m and rotational symmetry n
  m = AF_waidx(1,k)
  if (m<2) return ! skip aberrations of order less than 2, they have no 2nd derivative
  prefac = prefac / real(m)
  n = AF_waidx(2,k)
  im = m-n
  
  ! ------------  
  ! sum up binomial terms
  do l = 0, n
    
    ttmp1xx = 0.0
    ttmp1xy = 0.0
    ttmp1yy = 0.0
    
    ! get x-term sign and factor
    betax = sgnprm(mod(l,4)) * AF_BINOMIAL_TAB(n,l)
    ! get y-term sign and factor
    betay = sgnprm(mod(l+3,4)) * AF_BINOMIAL_TAB(n,l)
    
    ! calculate inner sum
    do i=0, im, 2
      xi = AF_BINOMIAL_TAB(im/2,i/2)
      ttmp1xx = ttmp1xx + xi*(m-l-i)/(i+l+1)  *wfield(1,m-l-i-1)*wfield(2,l+i+1)
      ttmp1xy = ttmp1xy + xi                  *wfield(1,m-l-i)  *wfield(2,l+i)
      ttmp1yy = ttmp1yy + xi*(l+i)  /(m-l-i+1)*wfield(1,m-l-i+1)*wfield(2,l+i-1)
    end do
    
    ! add to outer sum
    ttmpxxx = ttmpxxx + ttmp1xx * betax
    ttmpxxy = ttmpxxy + ttmp1xx * betay
    ttmpxyx = ttmpxyx + ttmp1xy * betax
    ttmpxyy = ttmpxyy + ttmp1xy * betay
    ttmpyyx = ttmpyyx + ttmp1yy * betax
    ttmpyyy = ttmpyyy + ttmp1yy * betay
    
  end do
  ! ------------
  
  bxxx = ttmpxxx * prefac ! final scaling for aberration function
  bxxy = ttmpxxy * prefac
  bxyx = ttmpxyx * prefac
  bxyy = ttmpxyy * prefac
  byyx = ttmpyyx * prefac
  byyy = ttmpyyy * prefac
! ------------

! ------------
!  write(*,*) " > AF_GetIDDBasis: EXIT."
  return

END SUBROUTINE AF_GetIDDBasis
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetTiltInductionBasis(nIdx1,nIdx2,tx,ty,bxx,byx,bxy,byy)
! function: Caluclates the tilt function for induction of lower order
!           aberration nIdx1 by higher-order aberration nIdx2 for
!           a given beam tilt (tx,ty).
!           The lower order inductions can then be calculated by
!           WLx = WHx * bxx + WHy * byx, and
!           WLy = WHx * bxy + WHy * byy,
!           where WL is the lower-order aberration induction
!           and WH is the higher-order aberration
!           Aberrations are identified by list indices
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx1          = list index of lower-order aberration
!            integer*4 :: nIdx2          = list index of higher-order aberration
!            real*4 :: tx, ty            = beam tilt in [rad]
!            real*4 :: bxx               = induction function from WHx to WLx
!            real*4 :: byx               = induction function from WHy to WLx
!            real*4 :: bxy               = induction function from WHx to WLy
!            real*4 :: byy               = induction function from WHy to WLy
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3800

  integer*4, intent(in) :: nIdx1, nIdx2
  real*4, intent(in) :: tx, ty
  real*4, intent(inout) :: bxx, byx, bxy, byy
  
  integer*4 :: i                                ! iterators
  integer*4 :: m1, m2, n1, n2                   ! aberration order and foldness
  integer*4 :: a, b, c, d                       ! aberration alternate indices
  integer*4 :: bf1, bf2                         ! binomial factors
  
  real*4 :: fg                                  ! global factor
  
  complex*8 :: tn, ts                           ! beam tilt and  c.c. as complex number
  complex*8 :: tn1, ts1, tn2, ts2               ! beam-tilt powers
  complex*8 :: cbp1, cbp2                       ! complex basis parts
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetTiltInductionBasis: INIT."
! get aberration indices
  m1 = AF_waidx(1,nIdx1)
  n1 = AF_waidx(2,nIdx1)
  m2 = AF_waidx(1,nIdx2)
  n2 = AF_waidx(2,nIdx2)
  c = (m1 + n1)/2
  d = (m1 - n1)/2
  a = (m2 + n2)/2
  b = (m2 - n2)/2
  fg = 1.0
  if (n1==0) fg = 0.5
  fg = fg * real(m1) / real(m2)
  bf1 = 0
  if ( a-c >= 0 .and. b-d >= 0 ) bf1 = AF_BINOMIAL_TAB(a,a-c) * AF_BINOMIAL_TAB(b,b-d)
  bf2 = 0
  if ( a-d >= 0 .and. b-c >= 0 ) bf1 = AF_BINOMIAL_TAB(a,a-d) * AF_BINOMIAL_TAB(b,b-c)
  tn = cmplx(tx, ty)
  ts = cmplx(tx,-ty)
! ------------

! ------------
! get tilt powers
  tn1 = cmplx(1.0,0.0)
  tn2 = cmplx(1.0,0.0)
  ts1 = cmplx(1.0,0.0)
  ts2 = cmplx(1.0,0.0)
  if ( a-c > 0 ) then
    do i=1, a-c
      ts1 = ts1 * ts
    end do
  end if
  if ( a-d > 0 ) then
    do i=1, a-d
      tn2 = tn2 * tn
    end do
  end if
  if ( b-c > 0 ) then
    do i=1, b-c
      ts2 = ts2 * ts
    end do
  end if
  if ( b-d > 0 ) then
    do i=1, b-d
      tn1 = tn1 * tn
    end do
  end if
! ------------


! ------------
! calculate complex basis parts
  cbp1 = fg*bf1*ts1*tn1
  cbp2 = fg*bf2*ts2*tn2
! ------------

! ------------
! calculate basis terms
  bxx = real(cbp1)+real(cbp2)
  byx = imag(cbp2)-imag(cbp1)
  bxy = imag(cbp1)+imag(cbp2)
  byy = real(cbp1)-real(cbp2)
! ------------


! ------------
!  write(*,*) " > AF_GetTiltInductionBasis: EXIT."
  return

END SUBROUTINE AF_GetTiltInductionBasis
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_GetEffectiveAberration(nIdx,tx,ty,weax,weay)
! function: Calculates the effective aberration (index nIdx) for given
!           beam tilt (tx, ty) on the basis of current aberration set.
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx         = aberration index
!            real*4 :: tx, ty          = beam tilt (rad)
!            real*4 :: weax, weay      = return ref. for eff. aberration
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4000
  integer*4, intent(in) :: nIdx
  real*4, intent(in) :: tx, ty
  real*4, intent(inout) :: weax, weay
  
  integer*4 :: i                        ! iterators
  integer*4 :: m, n                     ! aberration power and symmetry
  integer*4 :: mi, ni                   ! induction aberration power and symmetry
  real*4 :: tibxx, tibyx, tibxy, tibyy  ! tilt induction basis values
  real*4 :: waix, waiy                  ! temp. aberration coefficients
! ------------

! ------------
! INIT
!  write(*,*) " > AF_GetEffectiveAberration: INIT."
  weax = 0.0
  weay = 0.0
! ------------

! ------------
! pre-checks
  if (nIdx<1.or.nIdx>AF_maxaberration) return
  m = AF_waidx(1,nIdx)
  n = AF_waidx(2,nIdx)
! ------------

! ------------
! loop over all aberrations and calculate induction from this aberration
! to aberration(nIdx)
! default preset with axial aberration
  weax = AF_wa(1,nIdx)
  weay = AF_wa(2,nIdx)
  do i=1, AF_maxaberration
    mi = AF_waidx(1,i)
    if (mi<=m) cycle ! skip all aberrations of lower or same order than nIdx (no induction from there)
    ni = AF_waidx(2,i)
    waix = AF_wa(1,i)
    waiy = AF_wa(2,i)
    if ((waix*waix+waiy*waiy)<AF_POWERTHRESH_WA) cycle ! skip, aberration to weak
    ! get basis function values for aberration i on nIdx and given beam tilt
    call AF_GetTiltInductionBasis(nIdx,i,tx,ty,tibxx,tibyx,tibxy,tibyy)
    ! add induction
    weax = weax + waix*tibxx + waiy*tibyx
    weay = weay + waix*tibxy + waiy*tibyy
  end do
! ------------

! ------------
!  write(*,*) " > AF_GetEffectiveAberration: EXIT."
  return

END SUBROUTINE AF_GetEffectiveAberration
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION AF_AberrationFunctionRMSError(wx, wy)
! function: Calculates the RMS Error of the aberration function from
!           current aberration covariance matrix at given diffraction
!           angle (wx, wy) in the diffraction plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3300
  
  real*4 :: AF_AberrationFunctionRMSError

  real*4, intent(in) :: wx, wy
  
  integer*4 :: k, l
  real*4 :: prefac, covfac
  real*4 :: rtmp, bkx, bky, blx, bly
  real*4 :: wackxlx, wackxly, wackylx, wackyly
! ------------

! ------------
! INIT
!  write(*,*) " > AF_AberrationFunctionRMSError: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  rtmp = 0.0
  AF_AberrationFunctionRMSError = 0.0
! ------------


! ------------
! sum up relevant aberrations
  do k = 1, AF_maxaberration ! sum over covariance matrix
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
    
    ! Get basis functions
    call AF_GetBasis(k,wx,wy,bkx,bky)
    
    do l =k, AF_maxaberration
  
      if (AF_waact(l)==0) cycle ! skip inactive aberration
      
      if (l==k) then 
        covfac = 1.0
      else
        covfac = 2.0
      end if
  
      ! get aberration variance or covariance
      wackxlx = AF_wacov(2*(k-1)+1,2*(l-1)+1)
      wackxly = AF_wacov(2*(k-1)+1,2*(l-1)+2)
      wackylx = AF_wacov(2*(k-1)+2,2*(l-1)+1)
      wackyly = AF_wacov(2*(k-1)+2,2*(l-1)+2)
    
      ! get basis functions
      call AF_GetBasis(l,wx,wy,blx,bly)
    
      rtmp = rtmp + covfac*wackxlx*bkx*blx &
     &            +    2.0*wackxly*bkx*bly &
     &            +    2.0*wackylx*bky*blx &
     &            + covfac*wackyly*bky*bly
     
     end do
      
  end do
  
  AF_AberrationFunctionRMSError = sqrt(rtmp) * prefac ! final scaling for aberration function
! ------------


! ------------
!  write(*,*) " > AF_AberrationFunctionRMSError: EXIT."
  return

END FUNCTION AF_AberrationFunctionRMSError
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION AF_AberrationFunctionRMSRingError(wmax)
! function: Calculates the RMS Error of the aberration function from
!           current aberration covariance matrix at given diffraction
!           angle (wmax) in the diffraction plane, by averaging along the
!           wmax circle
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely AF_maxaberration_order and AF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wmax : diffraction angle of circle
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3400
  
  real*4 :: AF_AberrationFunctionRMSRingError

  real*4, intent(in) :: wmax
  
  integer*4 :: k, l, mk, ml, nk, nl
  real*4 :: prefac, covfac, roundfac
  real*4 :: rtmp
  real*4 :: wakp, walp, wfield(AF_maxaberration_order)
  real*4 :: wackxlx, wackyly
! ------------

! ------------
! INIT
!  write(*,*) " > AF_AberrationFunctionRMSRingError: INIT.
  prefac = 2.0*AF_pi/AF_lamb
  rtmp = 0.0
  AF_AberrationFunctionRMSRingError = 0.0
! prepare aberration order power list
  wfield(1) = wmax
  do k=2, AF_maxaberration_order
    wfield(k) = wfield(k-1)*wmax
  end do
! ------------


! ------------
! sum up relevant aberrations
  do k = 1, AF_maxaberration ! sum over covariance matrix
  
    if (AF_waact(k)==0) cycle ! skip inactive aberration
    
    ! get aberration rotational symmetry
    mk = AF_waidx(1,k)
    nk = AF_waidx(2,k)
    
    ! Get basis functions
    wakp = wfield(mk)/real(mk)
    
    do l =k, AF_maxaberration
  
      if (AF_waact(l)==0) cycle ! skip inactive aberration
      
      ! get aberration rotational symmetry
      ml = AF_waidx(1,l)
      nl = AF_waidx(2,l)
      
      
      if (nl/=nk) cycle ! rotational symmetry must be equal to produce RMS error average on
                        ! a circle
      
      wackxlx = 0.0
      wackyly = 0.0
      
      if (nl==0) then   ! round aberrations only have one one parameter, it counts double
        roundfac = 2.0
        ! get aberration variance or covariance
        wackxlx = AF_wacov(2*(k-1)+1,2*(l-1)+1)
      else
        roundfac = 1.0
        ! get aberration variance or covariance
        wackxlx = AF_wacov(2*(k-1)+1,2*(l-1)+1)
        wackyly = AF_wacov(2*(k-1)+2,2*(l-1)+2)
      end if
                        
      ! Get basis functions
      walp = wfield(ml)/real(ml)
      
      if (l==k) then    ! covariances appear double, variances only once
        covfac = 0.5
      else
        covfac = 1.0
      end if
      
      ! sum up
      rtmp = rtmp + roundfac*covfac*(wackxlx+wackyly)*wakp*walp
     
     end do
      
  end do
  
  rtmp = rtmp
  
  AF_AberrationFunctionRMSRingError = sqrt(rtmp) * prefac ! final scaling for aberration function
! ------------


! ------------
!  write(*,*) " > AF_AberrationFunctionRMSRingError: EXIT."
  return

END FUNCTION AF_AberrationFunctionRMSRingError
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE AF_RotateField(acin,n,rot,acout)
! function: changes aberration coefficients due to a rotation of the
!           image frame
!           the order of the aberrations is assumed to be according to
!           the current module setup
!           aberration coefficients are in cartesian notation
! -------------------------------------------------------------------- !
! parameter:
!   IN:     real*4 :: acin(2,n)     = input aberration coefficients
!           integer*4 :: n          = number of input aberrations
!           real*4 :: rot           = rotation angle [radian]
!   IN/OUT: real*4 :: acout(2,n)    = output aberration coefficients
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 4300

  integer*4, intent(in) :: n
  real*4, intent(in) :: acin(2,n), rot
  real*4, intent(inout) :: acout(2,n)
  
  integer*4 :: i, j             ! iterators
  integer*4 :: na               ! used number of aberrations
  real*4    :: arot             ! aberration rotation angle
  real*4    :: ax, ay           ! current aberration values
  real*4    :: ca, sa           ! cosine and sine values
! ------------

! ------------
! INIT
!  write(*,*) " > AF_RotateField: INIT."
  na = min(AF_maxaberration,n)
  if (na<1) then
    call AF_ERROR("Invalid input parameters or module not set up.",subnum + 1)
    return ! module not set up or no data
  end if
! precopy
  acout = acin
! ------------

! ------------
! loop over aberrations
  do i=1, na
    j = AF_waidx(2,i) ! save rotational symmetry number of aberration
    if (j==0) cycle ! do not change round aberrations
    arot = real(j)*rot ! aberration must be rotated j-times stronger than field of view
    ax = acin(1,i)
    ay = acin(2,i)
    ca = cos(arot)
    sa = sin(arot)
    ! apply rotation
    acout(1,i) = ax*ca-ay*sa
    acout(2,i) = ax*sa+ay*ca
  end do
! ------------

! ------------
!  write(*,*) " > AF_RotateField: EXIT."
  return

END SUBROUTINE AF_RotateField
!**********************************************************************!




END MODULE AberrationFunctions

!*********************************************************************!
!*********************************************************************!
!    MODULE END                                                       !
!*********************************************************************!
!*********************************************************************!



















!*********************************************************************!
!*********************************************************************!
!    COMMENTED TEMPLATE                                               !
!*********************************************************************!
!*********************************************************************!




!**********************************************************************!
!**********************************************************************!
!SUBROUTINE <NAME>(<PARAMS>)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

!  implicit none

! ------------
! DECLARATION
!  integer*4, parameter :: subnum = 4600
!
! ------------

! ------------
! INIT
!  write(*,*) " > <NAME>: INIT."
! ------------

! ------------
! 
! ------------

! ------------
!  write(*,*) " > <NAME>: EXIT."
!  return

!END SUBROUTINE <NAME>
!**********************************************************************!