!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  STEMfunction                                          !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2009-2018   !
!                                                                      !
!**********************************************************************!
!                                                                      !
!    MODULE STEMfunctions                                              !
!    --------------------                                              !
!                                                                      !
!    Purpose  : Implementation of STEM imaging functions and data      !
!               manipulations of simulation and analysis               !
!    Version  :  1.0.1, Dec 05, 2018                                   !
!    To Link  : SFFTs.f                                                !
!               BasicFuncs.f90                                         !
!                                                                      !
!**********************************************************************!
!                                                                       
!  Author:  Juri Barthel                                
!           Ernst Ruska-Centre                                          
!           Forschungszentrum Jülich GmbH, 52425 Jülich, Germany        
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
! This program is free software: you can redistribute it and/or modify  
! it under the terms of the GNU General Public License as published by  
! the Free Software Foundation, either version 3 of the License, or     
! (at your option) any later version. See gpl-3.txt                                   
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




!*********************************************************************!
!                                                                     !
!    IMORTANT REMARKS                                                 !
!    ----------------                                                 !
!                                                                     !
!    (1) CALL STF_INIT() before usage                                 !
!    (2) CALL STF_UNINIT() before re-usage or when halting            !
!                                                                     !
!*********************************************************************!



!*********************************************************************!
!*********************************************************************!
!    MODULE DECALRATIONS                                              !
!*********************************************************************!
!*********************************************************************!

MODULE STEMfunctions

! Use STF as acronym for public parameters funcs and subs!
    
! Global module dependencies
!   USE ...
!  USE iso_varying_string
    
  implicit none
  
  ! declare internal data types

  ! accessibility of subroutines or functions
!  private :: STF_***
  private :: STF_ERROR, STF_WARN
  private :: STF_TABBED_EXP
  private :: STF_TABBED_SIGMOID
  private :: STF_SETTAB_USC
  private :: STF_SETTAB_SCR
  private :: STF_SETTAB_USC2
  private :: STF_SETTAB_SCR2
  
  public :: STF_HT2WL
  public :: STF_WL2HT
  public :: STF_FFT
  public :: STF_PrepareProbeWaveFourier
  public :: STF_PrepareVortexProbeWaveFourier
  public :: STF_PreparePlaneWaveFourier
  public :: STF_AberrateWaveFourier
  public :: STF_ApplyObjectiveAperture
  public :: STF_AberrationFunction
  public :: STF_GetCohProbeWave
  public :: STF_GetCohProbeWaveRe
  public :: STF_GetCohProbeWaveIm
  public :: STF_GetPhasePlate
  
!  public :: STF_***
  public :: STF_INIT ! call first!
                     ! set STF_FFT_BOUND
  public :: STF_INIT_ALLOC ! call then
  public :: STF_UNINIT
  
  public :: STF_SetAberrationByName
  public :: STF_GetAberrationByName
  public :: STF_SetAberration
  public :: STF_GetAberration
  public :: STF_SetAberrationLName
  public :: STF_GetAberrationLName
  public :: STF_SetAberrationSName
  public :: STF_GetAberrationSName
  

!  declare module global (but private) variables, params and arrays
!  file units interval
  integer*4, private, parameter :: STF_minunit = 81
  integer*4, private, parameter :: STF_maxunit = 90

!   length of names and internal strings
  integer*4, public, parameter :: STF_ll = 1024
  
!   standard output unit
  integer*4, public, parameter :: STF_stdout = 6

!   enpty line placeholder (emptied at startup)
  character(len=STF_ll), private :: STF_el  
! pi
  real*4, private :: STF_pi
  DATA STF_pi /3.1415927/

! scale degree to radian
  real*4, private :: STF_rd2r
  DATA STF_rd2r /57.2957795/

! error counter  
  integer*4, public :: STF_err_num
  DATA STF_err_num /0/
  
! Max FFT size
  !integer*4, private, parameter :: FFT_BOUND = 8192
  !integer*4, private, parameter :: FFT_NYQ = 4096
  !integer*4, private, parameter :: FFT_BOUND = 4096
  !integer*4, private, parameter :: FFT_NYQ = 2048
  !integer*4, private, parameter :: FFT_BOUND = 2048
  !integer*4, private, parameter :: FFT_NYQ = 1024
  integer*4, public :: STF_FFT_BOUND
  DATA STF_FFT_BOUND /2048/
  integer*4, public :: STF_FFT_NYQ
  DATA STF_FFT_NYQ /1024/

! aberration power threshold
  real*4, private, parameter :: STF_POWERTHRESH_WA = 1.0E-30
  
! aperture power threshold
  real*4, private, parameter :: STF_APERTURETHRESH = 1.0E-03
  
! aperture angle thresholds
  real*4, private, parameter :: STF_RELANGTHRESH = 1.0E-3
  real*4, private, parameter :: STF_RELANGWIDTHTHRESH = 1.0E-2
  real*4, private, parameter :: STF_APSMOOTHPIX = 1.0
  

! explicit defocus average parameter defaults
  integer*4, private, parameter :: STF_DEFOCUS_KERNEL_STEPS_DEFAULT = 13
  real*4, private, parameter :: STF_DEFOCUS_KERNEL_SPREAD_DEFAULT = 3.0
! explicit defocus average parameters
  integer*4, public :: STF_DEFOCUS_KERNEL_STEPS
  real*4, public :: STF_DEFOCUS_KERNEL_SPREAD
  DATA STF_DEFOCUS_KERNEL_STEPS /13/
  DATA STF_DEFOCUS_KERNEL_SPREAD /3.0/

! data size in one dimension  
  integer*4, public :: STF_dim
  DATA STF_dim /256/

! wavelength [nm]
  real*4, public :: STF_lamb
  DATA STF_lamb /0.001969/
! high tension [kV]
  real*4, public :: STF_ht
  DATA STF_ht /300.0/

! current real-space sampling [nm/pix]
  real*4, public :: STF_sampling
  DATA STF_sampling /0.05/

! current Fourier space sampling [(pix*nm)^-1]
  real*4, public :: STF_itog, STF_itow

! ----------->
! ----------->
! IMPORTANT PARAMETER: STF_maxaberration_order
!   controls size of dynamic allocation, length of aberration lists
! ------------------O
! -----------------O|
!                  ||
!                  vv
!              CHANGE THIS PARAMETER
!              TO INCLUDE MORE HIGHER ORDER ABERRATIONS !!!
!              + CAREFULLY CHECK STF_INIT() and STF_UNINIT()
!                FOR CORRECT MEMORY ALLOCATION
! max number of aberration coefficients
  integer*4, public, parameter :: STF_maxaberration_order = 8
  
  integer*4, public, parameter :: STF_aberration_shortname_length = 4
  integer*4, public, parameter :: STF_aberration_longname_length = 40
   
! aberration tables, dynamically allocated
! size holder, dynamically determined in STF_INIT()
  integer*4, public :: STF_maxaberration
  DATA STF_maxaberration /0/
! allocatable tables
  !real*4, public :: STF_wa(2,STF_maxaberration) ! coeffictients
  !integer*4, public :: STF_waidx(2,STF_maxaberration) ! coeffictients index hash
  !character(len=2*STF_aberration_shortname_length), public &
  !   & :: STF_asn(STF_maxaberration) ! short names
  !character*STF_aberration_longname_length, public &
  !   & :: STF_aln(STF_maxaberration) ! long names
  real*4, dimension(:,:), public, allocatable :: STF_wa ! coeffictients
  integer*4, dimension(:,:), public, allocatable :: STF_waidx ! coeffictients index hash
  integer*1, dimension(:,:), public, allocatable :: STF_asn ! short names
  integer*1, dimension(:,:), public, allocatable :: STF_aln ! long names
     
! defocus spread
  real*4, public :: STF_defocusspread ! [nm]
  DATA STF_defocusspread /3.0/
  
! source radius
  real*4, public :: STF_srcradius ! [nm]
  DATA STF_srcradius /0.025/
  
! condenser aperture
  real*4, public :: STF_caperture ! [mrad]
  DATA STF_caperture /30.0/
  real*4, public :: STF_caperture_movex ! [mrad]
  DATA STF_caperture_movex /0.0/
  real*4, public :: STF_caperture_movey ! [mrad]
  DATA STF_caperture_movey /0.0/
  integer*4, public :: STF_cap_type ! 0 = sigmoid(top-hat), 1 = gaussian
  Data STF_cap_type /0/
  
! beam tilt
  real*4, public :: STF_beam_tiltx ! [mrad]
  DATA STF_beam_tiltx /0.0/
  real*4, public :: STF_beam_tilty ! [mrad]
  DATA STF_beam_tilty /0.0/
  
! data array used for calculation
  complex*8, allocatable, dimension(:,:), public :: STF_sc !, STF_sc2
  real*4, public :: STF_PreparedWavePower
  
! anglular function tables
  integer*4, parameter, private :: STF_ANGTAB_SIZE = 1023
  complex*8, dimension(STF_ANGTAB_SIZE), private :: STF_ANGTAB_EXP
! sigmoid function table
  real*4, parameter, private :: STF_SIGMOID_EXT = 3.0
  real*4, dimension(STF_ANGTAB_SIZE), private :: STF_ANGTAB_SIGMOID
! binomial coefficient tables
  integer*4, dimension(0:2*STF_maxaberration_order,0:2*STF_maxaberration_order) :: STF_BINOMIAL_TAB
! scramble and unscramble lists
  integer*4, allocatable, dimension(:), private :: STF_TABBED_SCR, STF_TABBED_USC  
  integer*4, allocatable, dimension(:), private :: STF_TABBED_SCR2, STF_TABBED_USC2  
  
! function names for general interface calling index
! ----------->
! ----------->
! IMPORTANT PARAMETER: STF_FUNC_NUM
!   controls number of callable functions
! ------------------O
! -----------------O|
!                  ||
!                  vv
!              CHANGE THIS PARAMETER
!              TO INCLUDE MORE FUNCTIONS !!!
!              + CAREFULLY CHECK STF_INIT()
!                FOR CORRECT NAME-PRESETS in STF_FuncNames
  integer*4, public, parameter :: STF_FUNC_NUM = 10 ! number of implemented functions
  integer*4, public, parameter :: STF_FUNC_NAME_LENGTH = 80 ! max. length of function names
  integer*1, dimension(STF_FUNC_NAME_LENGTH,STF_FUNC_NUM), public :: STF_FuncNames


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
SUBROUTINE STF_INIT()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 100
  integer*4, parameter :: asl = 2*STF_aberration_shortname_length
  integer*4, parameter :: all = STF_aberration_longname_length
  integer*4, parameter :: fnl = STF_FUNC_NAME_LENGTH
  
!
  integer*4 :: m, n, l, err
  real*4 :: angle, ascale
  integer*4, external :: binomial ! BasicFuncs.f90
  real*4, external :: sigmoid ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_INIT: INIT."
! ------------

! ------------
  STF_pi = atan(1.0)*4.0
  STF_rd2r = STF_pi/180.0
  STF_err_num = 0
  STF_el = REPEAT(" ",STF_ll)
! ------------

! ------------
  STF_DEFOCUS_KERNEL_STEPS = STF_DEFOCUS_KERNEL_STEPS_DEFAULT
  STF_DEFOCUS_KERNEL_SPREAD = STF_DEFOCUS_KERNEL_SPREAD_DEFAULT
! ------------

! ------------
! preset binomial coefficient table
  do m=0, STF_maxaberration_order
    do n=0,STF_maxaberration_order
      STF_BINOMIAL_TAB(n,m) = binomial(n,m)
    end do
  end do
! ------------

! ------------
! preset aberration tables

! (1) determine length of aberration list -> STF_maxaberration
  ! count up to max. order
  l = 0
  STF_maxaberration = 0
  do m=1, STF_maxaberration_order
    do n=0,m
      if (mod(m+n,2)==1) cycle
      l = l + 1
    end do
  end do
  if (l>0) then
    STF_maxaberration = l
    
    ! prepare angular tables now
    ascale = 2.0*STF_pi/real(STF_ANGTAB_SIZE)
    do m=1,STF_ANGTAB_SIZE
      angle = real(m-1)*ascale
      STF_ANGTAB_EXP(m) = cmplx(cos(angle),sin(angle))
    end do
    ascale = 2.0*STF_SIGMOID_EXT/real(STF_ANGTAB_SIZE-1)
    do m=1,STF_ANGTAB_SIZE
      angle = real(m-1)*ascale-STF_SIGMOID_EXT
      STF_ANGTAB_SIGMOID(m) = sigmoid(angle,0.0,1.0)
    end do
    
  else
    call STF_ERROR("STF_INIT failed: Invalid number of aberrations.",subnum+1)
    return
  end if
  
  ! allocate arrays for tables
  allocate(STF_wa(2,STF_maxaberration),stat=err)
  if (err/=0) then
    call STF_ERROR("STF_INIT failed: Allocation error.",subnum+2)
    return
  end if
  allocate(STF_waidx(2,STF_maxaberration),stat=err)
  if (err/=0) then
    call STF_ERROR("STF_INIT failed: Allocation error.",subnum+3)
    return
  end if
  allocate(STF_asn(asl,STF_maxaberration),stat=err)
  if (err/=0) then
    call STF_ERROR("STF_INIT failed: Allocation error.",subnum+4)
    return
  end if
  allocate(STF_aln(all,STF_maxaberration),stat=err)
  if (err/=0) then
    call STF_ERROR("STF_INIT failed: Allocation error.",subnum+5)
    return
  end if
  
  ! index hash up to max. order
  l = 1
  do m=1, STF_maxaberration_order
    do n=0,m
      if (mod(m+n,2)==1) cycle
      STF_waidx(1, l) = m
      STF_waidx(2, l) = n
      l = l + 1
    end do
  end do
  
  ! coefficients
  STF_wa(:,:) = 0.0
  STF_asn(:,:) = ' '
  STF_aln(:,:) = ' '
  
  ! ----------->
  ! ----------->
  ! IMPORTANT -> ADJUST THE FOLLOWING LISTS IF IF STF_maxaberration_order /= 8
  ! long names
  call STF_SetAberrationLName( 1, "Image shift ")
  call STF_SetAberrationLName( 2, "Defocus ")
  call STF_SetAberrationLName( 3, "Twofold astigmatism ")
  call STF_SetAberrationLName( 4, "3rd order axial coma ")
  call STF_SetAberrationLName( 5, "Threefold astigmatism ")
  call STF_SetAberrationLName( 6, "Spherical aberration (Cs) ")
  call STF_SetAberrationLName( 7, "Star aberration ")
  call STF_SetAberrationLName( 8, "Fourfold astigmatism ")
  call STF_SetAberrationLName( 9, "5th order axial coma ")
  call STF_SetAberrationLName(10, "Three lobe aberration ")
  call STF_SetAberrationLName(11, "Fivefold astigmatism ")
  call STF_SetAberrationLName(12, "Spherical aberration (C5) ")
  call STF_SetAberrationLName(13, "6th order Star aberration ")
  call STF_SetAberrationLName(14, "Rosette aberration ")
  call STF_SetAberrationLName(15, "Sixfold aberration ")
  call STF_SetAberrationLName(16, "7th order axial coma " )
  call STF_SetAberrationLName(17, "Threefoldness of the 7th order ")
  call STF_SetAberrationLName(18, "Fivefoldness of the 7th order ")
  call STF_SetAberrationLName(19, "Sevenfold astigmatism ")
  call STF_SetAberrationLName(20, "Spherical aberration (C7) ")
  call STF_SetAberrationLName(21, "Twofoldness of the 8th order ")
  call STF_SetAberrationLName(22, "Fourfoldness of the 8th order " )
  call STF_SetAberrationLName(23, "Sixfoldness of the 8th order ")
  call STF_SetAberrationLName(24, "Eightfold astigmatism ")
  ! short names / identifier
  call STF_SetAberrationSName( 1, "w11 W11 ")
  call STF_SetAberrationSName( 2, "w20 W20 ")
  call STF_SetAberrationSName( 3, "w22 W22 ")
  call STF_SetAberrationSName( 4, "w31 W31 ")
  call STF_SetAberrationSName( 5, "w33 W33 ")
  call STF_SetAberrationSName( 6, "w40 W40 ")
  call STF_SetAberrationSName( 7, "w42 W42 ")
  call STF_SetAberrationSName( 8, "w44 W44 ")
  call STF_SetAberrationSName( 9, "w51 W51 ")
  call STF_SetAberrationSName(10, "w53 W53 ")
  call STF_SetAberrationSName(11, "w55 W55 ")
  call STF_SetAberrationSName(12, "w60 W60 ")
  call STF_SetAberrationSName(13, "w62 W62 ")
  call STF_SetAberrationSName(14, "w64 W64 ")
  call STF_SetAberrationSName(15, "w66 W66 ")
  call STF_SetAberrationSName(16, "w71 W71 ")
  call STF_SetAberrationSName(17, "w73 W73 ")
  call STF_SetAberrationSName(18, "w75 W75 ")
  call STF_SetAberrationSName(19, "w77 W77 ")
  call STF_SetAberrationSName(20, "w80 W80 ")
  call STF_SetAberrationSName(21, "w82 W82 ")
  call STF_SetAberrationSName(22, "w84 W84 ")
  call STF_SetAberrationSName(23, "w86 W86 ")
  call STF_SetAberrationSName(24, "w88 W88 ")
! ------------

! ------------
  call SetVarString(STF_FuncNames(1:fnl,1),fnl,"Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,2),fnl,"Partially Spatial Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,3),fnl,"Partially Temporal Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,4),fnl,"Partially Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,5),fnl,"Partially Spatial Quasi-Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,6),fnl,"Partially Temporal Quasi-Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,7),fnl,"Partially Quasi-Coherent Probe Intensity")
  call SetVarString(STF_FuncNames(1:fnl,8),fnl,"Real Part of Probe Wave Function")
  call SetVarString(STF_FuncNames(1:fnl,9),fnl,"Imaginary Part of Probe Wave Function")
  call SetVarString(STF_FuncNames(1:fnl,10),fnl,"Phase Plate")
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_INIT: EXIT."
  return

END SUBROUTINE STF_INIT
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_INIT_ALLOC()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 105
!
  integer*4 :: nsdim, nalloc
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_INIT_ALLOC: INIT."
  nsdim = STF_FFT_BOUND ! set before calling to initialize FFT arrays
! ------------


  ! ------------
! allocate module arrays
  if (allocated(STF_sc)) deallocate(STF_sc,stat=nalloc)
  allocate(STF_sc(nsdim,nsdim), stat=nalloc)
  if (allocated(STF_TABBED_SCR)) deallocate(STF_TABBED_SCR,stat=nalloc)
  allocate(STF_TABBED_SCR(1:nsdim),stat=nalloc)
  if (allocated(STF_TABBED_USC)) deallocate(STF_TABBED_USC,stat=nalloc)
  allocate(STF_TABBED_USC(1:nsdim),stat=nalloc)
  if (allocated(STF_TABBED_SCR2)) deallocate(STF_TABBED_SCR2,stat=nalloc)
  allocate(STF_TABBED_SCR2(1:nsdim),stat=nalloc)
  if (allocated(STF_TABBED_USC2)) deallocate(STF_TABBED_USC2,stat=nalloc)
  allocate(STF_TABBED_USC2(1:nsdim),stat=nalloc)
! ------------


! ------------
!  write(unit=*,fmt=*) " > STF_INIT_ALLOC: EXIT."
  return

END SUBROUTINE STF_INIT_ALLOC
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_UNINIT()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 200
  integer*4 :: err, nalloc
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_INIT: STF_UNINIT."
  err = 0
! ------------

! ------------
! deallocate tables
  deallocate(STF_wa,STF_waidx,STF_asn,STF_aln,stat=err)
  STF_maxaberration = 0
  if (err/=0) then
    call STF_ERROR("STF_INIT failed: Deallocation error.",subnum+1)
    return
  end if
  if (allocated(STF_sc)) deallocate(STF_sc,stat=nalloc)
  if (allocated(STF_TABBED_SCR)) deallocate(STF_TABBED_SCR,stat=nalloc)
  if (allocated(STF_TABBED_USC)) deallocate(STF_TABBED_USC,stat=nalloc)
  if (allocated(STF_TABBED_SCR2)) deallocate(STF_TABBED_SCR2,stat=nalloc)
  if (allocated(STF_TABBED_USC2)) deallocate(STF_TABBED_USC2,stat=nalloc)
! reset aberration list length
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_UNINIT: EXIT."
  return

END SUBROUTINE STF_UNINIT
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_ERROR(sTxt,nErr)
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
  character(len=STF_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = STF_el
! ----------

! ----------
! messaging (to screen)
  write(unit=sinfo,fmt=*) "ERROR: ",trim(sTxt)," Error code:",nErr
!  call SE_event(trim(sinfo), SE_err)
  write(unit=STF_stdout,fmt='(A)') "STF_ERROR: "//trim(sinfo)
! ----------

  STF_err_num = STF_err_num + 1
  return

END SUBROUTINE STF_ERROR
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_WARN(sTxt)
! function: print warning message
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  integer*4, parameter :: subnum = 300
  
  character(len=*), intent(in) :: sTxt
  character(len=STF_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = STF_el
! ----------

! ----------
! messaging (to screen)
  write(unit=sinfo,fmt=*) "Warning: ",trim(sTxt)
!  call SE_event(trim(sinfo), SE_err)
  write(unit=STF_stdout,fmt='(A)') trim(sinfo)
! ----------

  return

END SUBROUTINE STF_WARN
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION STF_SetAberrationByName(sname,wax,way)
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
  integer*4, parameter :: asl = 2*STF_aberration_shortname_length
  
  integer*4 :: STF_SetAberrationByName
  character(len=*), intent(in) :: sname
  real*4, intent(in) :: wax, way
  
  character(len=asl) :: sourcestring, searchstring
  integer*4 :: i, j, isl
  
  external :: GetVarString
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SetAberrationByName: INIT."
  STF_SetAberrationByName = 0 ! start with unsuccessful state
  searchstring(1:asl) = STF_el(1:asl)
  isl = len(sname)
  searchstring(1:isl) = sname(1:isl)
! ------------

! ------------
! try to locate aberration
  do i=1, STF_maxaberration
    call GetVarString(sourcestring,STF_asn(1:asl,i),asl)
    j = INDEX(sourcestring,trim(searchstring))
    if (j>0) then ! aberration identified
      STF_wa(1,i) = wax
      STF_wa(2,i) = way
      STF_SetAberrationByName = 1
      exit
    end if
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SetAberrationByName: EXIT."
  return

END FUNCTION STF_SetAberrationByName
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION STF_GetAberrationByName(sname,wax,way)
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
  integer*4, parameter :: asl = 2*STF_aberration_shortname_length
  
  integer*4 :: STF_GetAberrationByName
  character(len=*), intent(in) :: sname
  real*4, intent(out) :: wax, way
  
  character(len=asl) :: sourcestring, searchstring
  integer*4 :: i, j, isl
  
  external :: GetVarString
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetAberrationByName: INIT."
  STF_GetAberrationByName = 0 ! start with unsuccessful state
  searchstring(1:asl) = STF_el(1:asl)
  isl = len(sname)
  searchstring(1:isl) = sname(1:isl)
! ------------

! ------------
! try to locate aberration
  do i=1, STF_maxaberration
    call GetVarString(sourcestring,STF_asn(1:asl,i),asl)
    j = INDEX(sourcestring,trim(searchstring))
    if (j>0) then ! aberration identified
      wax = STF_wa(1,i)
      way = STF_wa(2,i)
      STF_GetAberrationByName = 1
      exit
    end if
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetAberrationByName: EXIT."
  return

END FUNCTION STF_GetAberrationByName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SetAberration(nidx,wax,way)
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
!  write(unit=*,fmt=*) " > STF_SetAberration: INIT."
  if ((nidx<1).or.(nidx>STF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  STF_wa(1,nidx) = wax
  STF_wa(2,nidx) = way
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SetAberration: EXIT."
  return

END SUBROUTINE STF_SetAberration
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetAberration(nidx,wax,way)
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
!  write(unit=*,fmt=*) " > STF_GetAberration: INIT."
  if ((nidx<1).or.(nidx>STF_maxaberration)) return ! bounce invalid index
! ------------

! ------------
! try to locate aberration
  wax = STF_wa(1,nidx)
  way = STF_wa(2,nidx)
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetAberration: EXIT."
  return

END SUBROUTINE STF_GetAberration
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SetAberrationLName(nIdx,sLongName)
! function: sets data to long name array
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            character(len=*) :: sLongName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1800
  integer*4, parameter :: all = STF_aberration_longname_length

  integer*4, intent(in) :: nIdx
  character(len=*), intent(in) :: sLongName
  
  integer*4 :: nLen
  external :: SetVarString
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SetAberrationLName: INIT."
  nLen = min(all,len(sLongName))
! ------------

! ------------
  call SetVarString(STF_aln(1:nLen, nIdx),nLen,sLongName(1:nLen))
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SetAberrationLName: EXIT."
  return

END SUBROUTINE STF_SetAberrationLName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetAberrationLName(nIdx,sLongName)
! function: retrieves name data
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            character(len=*) :: sLongName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1900
  integer*4, parameter :: all = STF_aberration_longname_length

  integer*4, intent(in) :: nIdx
  character(len=*), intent(inout) :: sLongName
  
  integer*4 :: nLen
  external :: GetVarString
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetAberrationLName: INIT."
  if (.not.allocated(STF_aln)) return
  nLen = min(all,len(sLongName))
! ------------

! ------------
  call GetVarString(sLongName(1:nLen),STF_aln(1:nLen, nIdx),nLen)
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetAberrationLName: EXIT."
  return

END SUBROUTINE STF_GetAberrationLName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SetAberrationSName(nIdx,sShortName)
! function: sets data to long name array
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            character(len=*) :: sShortName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2000
  integer*4, parameter :: asl = 2*STF_aberration_shortname_length

  integer*4, intent(in) :: nIdx
  character(len=*), intent(in) :: sShortName
  
  integer*4 :: nLen
  external :: SetVarString
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SetAberrationSName: INIT."
  nLen = min(asl,len(sShortName))
! ------------

! ------------
  call SetVarString(STF_asn(1:nLen, nIdx),nLen, sShortName(1:nLen))
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SetAberrationSName: EXIT."
  return

END SUBROUTINE STF_SetAberrationSName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetAberrationSName(nIdx,sShortName)
! function: retrieves short name data
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nIdx : index
!            character(len=*) :: sShortName : name to set
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2100
  integer*4, parameter :: asl = 2*STF_aberration_shortname_length

  integer*4, intent(in) :: nIdx
  character(len=*), intent(inout) :: sShortName
  
  integer*4 :: nLen
  external :: GetVarString
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetAberrationSName: INIT."
  nLen = min(asl,len(sShortName))
! ------------

! ------------
  call GetVarString(sShortName(1:nLen),STF_asn(1:nLen, nIdx),nLen)
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetAberrationSName: EXIT."
  return

END SUBROUTINE STF_GetAberrationSName
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_TABBED_EXP(rangle,cval)
! function: returns tabbed complex exp(I*rangle)
! -------------------------------------------------------------------- !
! parameter: real*4 :: rangle [radian]
!            complex*8 :: cval
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL STF_INIT() BEFORE CALLING THIS ROUTINE !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2400
  real*4 :: rangle
  complex*8 :: cval
  real*4 :: iangle, modangle
  integer*4 :: idx
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_TABBED_EXP: INIT."
  modangle = 2.0 * STF_pi
! ------------

! ------------
! parameter preset
  if (rangle>=0.0) then
    iangle = modulo(rangle,modangle)
    iangle = iangle/modangle*real(STF_ANGTAB_SIZE)
  else
    iangle = modulo(rangle,modangle)
    iangle = iangle/modangle*real(STF_ANGTAB_SIZE)
  end if
  idx = 1+int(iangle) ! round to next integer
  if (idx==0) then
    idx = STF_ANGTAB_SIZE
  end if
  if (idx==STF_ANGTAB_SIZE+1) then
    idx = 1
  end if
! clip angle to interval 0.. STF_ANGTAB_SIZE
! ------------

! ------------
! get data
  cval = STF_ANGTAB_EXP(idx)
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_TABBED_EXP: EXIT."
  return

END SUBROUTINE STF_TABBED_EXP
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_TABBED_SIGMOID(x,x0,dx,rval)
! function: returns tabbed real*4 sigmoid(x,x0,dx)
! -------------------------------------------------------------------- !
! parameter: real*4 :: x,x0,dx,rval
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL STF_INIT() BEFORE CALLING THIS ROUTINE !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2500
  real*4 :: x,x0,dx,rval
  real*4 :: tx, ix, maxx
  integer*4 :: idx
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_TABBED_EXP: INIT."
  maxx = STF_SIGMOID_EXT
  tx = (x-x0)/dx
! ------------

! ------------
! parameter preset
  if (tx>maxx) then
    rval = 1.0
    return
  end if
  if (tx<-maxx) then
    rval = 0.0
    return
  end if
  ix = 0.5*(tx/STF_SIGMOID_EXT+1.0)*real(STF_ANGTAB_SIZE-1)
  idx = 1+int(ix) ! round to next integer
  if (idx==0) then
    idx = STF_ANGTAB_SIZE
  end if
  if (idx==STF_ANGTAB_SIZE+1) then
    idx = 1
  end if
! clip angle to interval 0.. STF_ANGTAB_SIZE
! ------------

! ------------
! get data
  rval = STF_ANGTAB_SIGMOID(idx)
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_TABBED_EXP: EXIT."
  return

END SUBROUTINE STF_TABBED_SIGMOID
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SETTAB_SCR(ndim)
! function: presets scramble tab
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndim
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL BEFORE USING STF_TABBED_SCR !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2900
  integer*4 :: ndim, ndim2, ndim2m1
  integer*4 :: i
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SETTAB_SCR: INIT."
  ndim2 = ndim/2
  ndim2m1 = ndim2-1
! ------------

! ------------
! parameter preset
  STF_TABBED_SCR = 1
  do i=1,ndim
    STF_TABBED_SCR(i)=mod((i+ndim2m1),ndim)-ndim2
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SETTAB_SCR: EXIT."
  return

END SUBROUTINE STF_SETTAB_SCR
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SETTAB_USC(ndim)
! function: presets unscramble tab
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndim
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL BEFORE USING STF_TABBED_USC !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3000
  integer*4 :: ndim, ndim2
  integer*4 :: i
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SETTAB_USC: INIT."
  ndim2 = ndim/2
! ------------

! ------------
! parameter preset
  STF_TABBED_USC = 1
  do i=1,ndim
    STF_TABBED_USC(i)=mod(i-1+ndim2,ndim)+1
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SETTAB_USC: EXIT."
  return

END SUBROUTINE STF_SETTAB_USC
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SETTAB_SCR2(ndim)
! function: presets scramble tab
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndim
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL BEFORE USING STF_TABBED_SCR !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 29000
  integer*4 :: ndim, ndim2, ndim2m1
  integer*4 :: i
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SETTAB_SCR2: INIT."
  ndim2 = ndim/2
  ndim2m1 = ndim2-1
! ------------

! ------------
! parameter preset
  STF_TABBED_SCR2 = 1
  do i=1,ndim
    STF_TABBED_SCR2(i)=mod((i+ndim2m1),ndim)-ndim2
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SETTAB_SCR2: EXIT."
  return

END SUBROUTINE STF_SETTAB_SCR2
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SETTAB_USC2(ndim)
! function: presets unscramble tab
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndim
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL BEFORE USING STF_TABBED_USC !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 30000
  integer*4 :: ndim, ndim2
  integer*4 :: i
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_SETTAB_USC2: INIT."
  ndim2 = ndim/2
! ------------

! ------------
! parameter preset
  STF_TABBED_USC2 = 1
  do i=1,ndim
    STF_TABBED_USC2(i)=mod(i-1+ndim2,ndim)+1
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_SETTAB_USC2: EXIT."
  return

END SUBROUTINE STF_SETTAB_USC2
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_FFT(cdata,nx,ny,dir)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3400
  integer*4, intent(in) :: nx,ny
  character*(*), intent(in) :: dir
  complex*8, intent(inout) :: cdata(STF_FFT_BOUND,STF_FFT_BOUND)
  
  integer*4 :: transformed
  external :: ODDCC128S, ODDCC256S, ODDCC512S, ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_FFT: INIT."
  transformed = 0
! ------------

! ------------
  select case (STF_FFT_BOUND)
  
  case (128)
    call ODDCC128S(cdata,nx,ny,dir)
    transformed = 1
    
  case (256)
    call ODDCC256S(cdata,nx,ny,dir)
    transformed = 1
    
  case (512)
    call ODDCC512S(cdata,nx,ny,dir)
    transformed = 1
    
  case (1024)
    call ODDCC1024S(cdata,nx,ny,dir)
    transformed = 1
  
  case (2048)
    call ODDCC2048S(cdata,nx,ny,dir)
    transformed = 1
    
  case (4096)
    call ODDCC4096S(cdata,nx,ny,dir)
    transformed = 1
    
  case (8192)
    call ODDCC8192S(cdata,nx,ny,dir)
    transformed = 1
  
  end select ! (STF_FFT_BOUND)
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_FFT: EXIT."
  return

END SUBROUTINE STF_FFT
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
FUNCTION STF_AberrationFunction(wx, wy)
! function: Calculates the aberration function from current aberration
!           vector at given diffraction angle (wx, wy) in the diffraction
!           plane
!           limits of the aberration function are determined by
!           the parameters defined in the declaration section of this
!           module, namely STF_maxaberration_order and STF_maxaberration
! -------------------------------------------------------------------- !
! parameter: real*4 :: wx, wy : diffraction angle (wx, wy)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  real*4, dimension(0:3), parameter :: sgnprm = (/ 1., 0., -1., 0. /)

  real*4 :: STF_AberrationFunction

  real*4, intent(in) :: wx, wy
  
  integer*4 :: j,k,l,m,n
  real*4 :: wfield(2,0:STF_maxaberration_order)
  real*4 :: wabs(0:STF_maxaberration_order)
  real*4 :: pwx, pwy, prefac, w, pw
  real*4 :: rtmp, wax, way, ttmp, ttmp1, tsgn
  
  integer*4, external :: binomial ! BasicFuncs.f90
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_AberrationFunction: INIT.
  prefac = 2.0*STF_pi/STF_lamb
  rtmp = 0.0
  STF_AberrationFunction = 0.0
! ------------

! ------------
! preset wfield and wabs with powers of diffraction angle
  pwx = 1.0
  pwy = 1.0
  pw = 1.0
  w = sqrt(wx*wx+wy*wy)
  do j = 0,STF_maxaberration_order
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
  do k = 1, STF_maxaberration ! sum over list!
    ! get power of aberration
    wax = stf_wa(1,k)
    way = stf_wa(2,k)
    
    if ((wax*wax+way*way)>STF_POWERTHRESH_WA) then ! use this aberration
      ttmp = 0.0  
      ! get aberration order m and rotational symmetry n
      m = STF_waidx(1,k)
      n = STF_waidx(2,k)
      ! sum up binimial terms
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
        ttmp = ttmp + ttmp1 * STF_BINOMIAL_TAB(n,l)*wfield(1,j)*wfield(2,l)
        
      end do
      
      j = m - n
      ttmp = ttmp * wabs(j) / real(m) ! final scaling for current term
      rtmp = rtmp + ttmp
    end if
    
  end do
  STF_AberrationFunction = rtmp * prefac ! final scaling for aberration function
! ------------


! ------------
!  write(unit=*,fmt=*) " > STF_AberrationFunction: EXIT."
  return

END FUNCTION STF_AberrationFunction
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetCohProbeWave(re,im,ndim,sampling)
! function: calculates the coherent probe wavefunction from the current
!           coherent aberrations in real-space, with re=real part
!           and im=imaginary part
! -------------------------------------------------------------------- !
! parameter: real*4 :: re(ndim,ndim) : reference to array recieving the
!                                     probe wavefunction real part
!            real*4 :: im(ndim,ndim) : reference to array recieving the
!                                     probe wavefunction imaginary part
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !

  implicit none


! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1700
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4, intent(out) :: re(ndim,ndim),im(ndim,ndim)
  
  integer*4 :: i, j, ndim2, i1, j1
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetCohProbeWave: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_GetCohProbeWave: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  STF_dim = ndim
  call STF_SETTAB_USC(ndim)
  ndim2 = int(ndim/2)
  STF_sampling = sampling
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
! ------------

! ------------
! preset Fourier-space array with exp[-i*chi]
  call STF_PrepareProbeWaveFourier(ndim,ndim,sampling,sampling)
! ------------

! ------------
! transform to real space (scramble first!)
  call STF_FFT(STF_sc,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate wave components by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = STF_sc(i,j)
      re(i1,j1) = real(cval)
      im(i1,j1) = imag(cval)
!      re(i,j) = real(cval)
!      im(i,j) = imag(cval)
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetCohProbeWave: EXIT."
  return

END SUBROUTINE STF_GetCohProbeWave
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetCohProbeWaveRe(re,ndim,sampling)
! function: calculates the coherent probe wavefunction from the current
!           coherent aberrations in real-space, real-part!!
! -------------------------------------------------------------------- !
! parameter: real*4 :: re(ndim,ndim) : reference to array recieving the
!                                     probe wavefunction real part
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !

  implicit none


! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2700
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4 :: re(ndim,ndim)
  
  integer*4 :: i, j, ndim2, i1, j1, ndim2m1
  real*4 :: wx, wy, wy2, chi, power, anglethresh, threshwidth
  real*4 :: tpower, tpowscal, rval
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetCohProbeWave: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_GetCohProbeWaveRe: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  STF_dim = ndim
  ndim2 = int(ndim/2)
  ndim2m1 = ndim2-1
  STF_sampling = sampling
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  anglethresh = STF_caperture*0.001
  threshwidth = anglethresh*0.1
! ------------

! ------------
! preset Fourier-space array with exp[-i*chi]
!!!  call STF_GetProbeWaveFourier(sc,ndim,sampling)
! ------------
! preset Fourier-space array with exp[-i*chi]
  STF_sc(:,:) = cmplx(0.0,0.0)
  tpower = 0.0
  do j=1,ndim
    wy = STF_TABBED_SCR(j)*STF_itow
    wy2 = wy*wy
    do i=1,ndim
      wx = STF_TABBED_SCR(i)*STF_itow
!      power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      chi = STF_AberrationFunction(wx,wy)
!      wave(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      call STF_TABBED_SIGMOID(sqrt(wx*wx+wy2),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = STF_AberrationFunction(wy,wx)  ! TRANSPOSED
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      STF_sc(i,j) = cval*power
!      STF_sc(i,j) = cmplx(cos(chi),-sin(chi))*power
      tpower = tpower + power*power
      
    end do
  end do
!  write(unit=*,fmt=*) "FS-power:",tpower
  tpowscal = 1.0/tpower/real(ndim*ndim)
! ------------
! ------------

! ------------
! transform to real space (scramble first!)
  call STF_FFT(STF_sc,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate wave components by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = STF_sc(i,j)
      re(i1,j1) = real(cval)*tpowscal
!      im(j1,i1) = imag(cval)
!      re(i,j) = real(cval)
!      im(i,j) = imag(cval)
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetCohProbeWaveRe: EXIT."
  return

END SUBROUTINE STF_GetCohProbeWaveRe
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetCohProbeWaveIm(im,ndim,sampling)
! function: calculates the coherent probe wavefunction from the current
!           coherent aberrations in real-space, imaginary-part!!
! -------------------------------------------------------------------- !
! parameter: real*4 :: im(ndim,ndim) : reference to array recieving the
!                                     probe wavefunction real part
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !

  implicit none


! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2800
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4 :: im(ndim,ndim)
  
  integer*4 :: i, j, ndim2, i1, j1, ndim2m1
  real*4 :: wx, wy, wy2, chi, power, anglethresh, threshwidth
  real*4 :: tpower, tpowscal, rval
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetCohProbeWaveIm: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_GetCohProbeWaveIm: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  STF_dim = ndim
  ndim2 = int(ndim/2)
  ndim2m1 = ndim2-1
  STF_sampling = sampling
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  anglethresh = STF_caperture*0.001
  threshwidth = anglethresh*0.1
! ------------

! ------------
! preset Fourier-space array with exp[-i*chi]
!!!  call STF_GetProbeWaveFourier(sc,ndim,sampling)
! ------------
! preset Fourier-space array with exp[-i*chi]
  tpower = 0.0
  STF_sc(:,:) = cmplx(0.0,0.0)
  do j=1,ndim
    wy = STF_TABBED_SCR(j)*STF_itow !(mod((j+ndim2m1),ndim)-ndim2)*STF_itow
    wy2 = wy*wy
    do i=1,ndim
      wx = STF_TABBED_SCR(i)*STF_itow !(mod((i+ndim2m1),ndim)-ndim2)*STF_itow
!      power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      chi = STF_AberrationFunction(wx,wy)
!      wave(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      call STF_TABBED_SIGMOID(sqrt(wx*wx+wy2),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = STF_AberrationFunction(wy,wx) ! TRANSPOSED
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      STF_sc(i,j) = cval*power
!      STF_sc(i,j) = cmplx(cos(chi),-sin(chi))*power
      tpower = tpower + power*power
      
    end do
  end do
!  write(unit=*,fmt=*) "FS-power:",tpower
  tpowscal = 1.0/tpower/real(ndim*ndim)
! ------------
! ------------

! ------------
! transform to real space (scramble first!)
  call STF_FFT(STF_sc,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate wave components by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = STF_sc(i,j)
!      re(j1,i1) = real(cval)*tpowscal
      im(i1,j1) = imag(cval)*tpowscal
!      re(i,j) = real(cval)
!      im(i,j) = imag(cval)
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetCohProbeWaveIm: EXIT."
  return

END SUBROUTINE STF_GetCohProbeWaveIm
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_PrepareProbeWaveFourier(ndimx,ndimy,samplingx,samplingy)
! function: calculates the prove wavefunction in fourier space
!           from aberration and saves it in global public array STF_sc
!           complex*8 array of size (FFT_BOUND,FFT_BOUND)
!           Data is scrambled and transposed
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndimx,ndimy : wave result dimension
!            real*4 :: samplingx,samplingy : wave result sampling x & y
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1700
  
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy
  
  integer*4 :: i, j, ndim2x, ndim2y
  real*4 :: wx, wy, wx2, chi, power, anglethresh, threshwidth
  real*4 :: tpower, rval, itogx, itogy, itowx, itowy
  real*4 :: camx, camy, lbtx, lbty, wap
  complex*8 :: cval
  
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_PrepareProbeWaveFourier: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_PrepareProbeWaveFourier: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndimx)
  call STF_SETTAB_SCR2(ndimy)
  !STF_dim = ndim
  ndim2x = int(ndimx/2)
  ndim2y = int(ndimy/2)
  !STF_sampling = sampling
! set Fourier-space sampling
  !STF_itog = 0.5/sampling/real(ndim2)
  !STF_itow = STF_itog*STF_lamb
  itogx = 0.5/samplingx/real(ndim2x)
  itowx = itogx*STF_lamb
  itogy = 0.5/samplingy/real(ndim2y)
  itowy = itogy*STF_lamb
  !wave = cmplx(0.0,0.0)
  anglethresh = max(STF_caperture*STF_RELANGTHRESH,0.25*(itowx+itowy))
  threshwidth = 0.5*(itowx+itowy)*STF_APSMOOTHPIX
  STF_PreparedWavePower = 0.0
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------


! ------------
! preset Fourier-space array with exp[-i*chi]
  STF_sc(:,:) = 0.0
  tpower = 0.0
  do j=1,ndimx
    wx = STF_TABBED_SCR(j)*itowx+lbtx
    wx2 = (wx-camx)**2
    do i=1,ndimy
      wy = STF_TABBED_SCR2(i)*itowy+lbty
      wap = wx2 + (wy-camy)**2
      if (STF_cap_type==0) then ! set sigmoid aperture function
        call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)
        power = 1.0-rval
      else if (STF_cap_type==1) then ! set gaussian aperture function
        power = exp(-wap/(anglethresh*anglethresh))
      else ! set zero beam only
        power = 0.0
        if (sqrt(wap)<min(itowx,itowy)) then
          power = 1.0
        end if
      end if
      
      
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = STF_AberrationFunction(wx,wy)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
      STF_sc(i,j) = cval*power
      STF_PreparedWavePower = STF_PreparedWavePower + power*power
      
    end do
  end do
!  write(unit=*,fmt=*) "FS-power:",tpower
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_PrepareProbeWaveFourier: EXIT."
  return

END SUBROUTINE STF_PrepareProbeWaveFourier
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_PrepareVortexProbeWaveFourier(oam,ndimx,ndimy,samplingx,samplingy)
! function: calculates the prove wavefunction in fourier space
!           from aberration and saves it in global public array STF_sc
!           complex*8 array of size (FFT_BOUND,FFT_BOUND)
!           Data is scrambled and transposed
!           This creates a STEM vortex probe.
! -------------------------------------------------------------------- !
! parameter: integer*4 :: oam : vortex orbital angular momentum
!            integer*4 :: ndimx,ndimy : wave result dimension
!            real*4 :: samplingx,samplingy : wave result sampling x & y
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1710
  
  integer*4 :: oam
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy
  
  integer*4 :: i, j, ndim2x, ndim2y
  real*4 :: wx, wy, wx2, chi, power, anglethresh, threshwidth
  real*4 :: tpower, rval, itogx, itogy, itowx, itowy
  real*4 :: camx, camy, lbtx, lbty, wap, vtx
  complex*8 :: cval
  
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_PrepareVortexProbeWaveFourier: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_PrepareVortexProbeWaveFourier: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndimx)
  call STF_SETTAB_SCR2(ndimy)
  !STF_dim = ndim
  ndim2x = int(ndimx/2)
  ndim2y = int(ndimy/2)
  !STF_sampling = sampling
! set Fourier-space sampling
  !STF_itog = 0.5/sampling/real(ndim2)
  !STF_itow = STF_itog*STF_lamb
  itogx = 0.5/samplingx/real(ndim2x)
  itowx = itogx*STF_lamb
  itogy = 0.5/samplingy/real(ndim2y)
  itowy = itogy*STF_lamb
  !wave = cmplx(0.0,0.0)
  anglethresh = max(STF_caperture*STF_RELANGTHRESH,0.25*(itowx+itowy))
  threshwidth = 0.5*(itowx+itowy)*STF_APSMOOTHPIX
  STF_PreparedWavePower = 0.0
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------


! ------------
! preset Fourier-space array with exp[-i*chi]
  STF_sc(:,:) = 0.0
  tpower = 0.0
  do j=1,ndimx
    wx = STF_TABBED_SCR(j)*itowx+lbtx
    wx2 = (wx-camx)**2
    do i=1,ndimy
      wy = STF_TABBED_SCR2(i)*itowy+lbty
      wap = wx2 + (wy-camy)**2
      if (STF_cap_type==0) then ! set sigmoid aperture function
        call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)
        power = 1.0-rval
      else if (STF_cap_type==1) then ! set gaussian aperture function
        power = exp(-wap/(anglethresh*anglethresh))
      else ! set zero beam only
        power = 0.0
        if (sqrt(wap)<min(itowx,itowy)) then
          power = 1.0
        end if
      end if
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      vtx = 0. ! preset vortex phase with as zero
      if (wx==0.) then ! calculate the vortex component of the aberration function
        if (wy>0.) vtx = -real(oam)*0.5*STF_pi + 0.5*STF_pi
        if (wy<0.) vtx =  real(oam)*0.5*STF_pi + 0.5*STF_pi
      else ! wx values other than 0
        vtx = -real(oam)*atan2(wy,wx) + 0.5*STF_pi
      end if
      chi = STF_AberrationFunction(wx,wy) + vtx ! combine lens aberrations and a pure vortex aberration
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
      STF_sc(i,j) = cval*power
      STF_PreparedWavePower = STF_PreparedWavePower + power*power
      
    end do
  end do
!  write(unit=*,fmt=*) "FS-power:",tpower
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_PrepareVortexProbeWaveFourier: EXIT."
  return

END SUBROUTINE STF_PrepareVortexProbeWaveFourier
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_PreparePlaneWaveFourier(ndimx,ndimy,samplingx,samplingy)
! function: prepares a plane wavefunction in fourier space
!           and saves it in global public array STF_sc
!           complex*8 array of size (FFT_BOUND,FFT_BOUND)
!           Data is scrambled and transposed
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndimx,ndimy : wave result dimension
!            real*4 :: samplingx,samplingy : wave result sampling x & y
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3500
  
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy
  
  integer*4 :: i, j, ndim2x, ndim2y
  real*4 :: wx, wy, chi, anglethresh, threshwidth
  real*4 :: tpower, itogx, itogy, itowx, itowy
  real*4 :: lbtx, lbty
  
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_PreparePlaneWaveFourier: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_PreparePlaneWaveFourier: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndimx)
  call STF_SETTAB_SCR2(ndimy)
  !STF_dim = ndim
  ndim2x = int(ndimx/2)
  ndim2y = int(ndimy/2)
  !STF_sampling = sampling
! set Fourier-space sampling
  !STF_itog = 0.5/sampling/real(ndim2)
  !STF_itow = STF_itog*STF_lamb
  itogx = 0.5/samplingx/real(ndim2x)
  itowx = itogx*STF_lamb
  itogy = 0.5/samplingy/real(ndim2y)
  itowy = itogy*STF_lamb
  !wave = cmplx(0.0,0.0)
  anglethresh = 0.25*(itowx+itowy)
  threshwidth = anglethresh*0.25
  STF_PreparedWavePower = 0.0
  lbtx = STF_beam_tiltx*0.001/STF_lamb
  lbty = STF_beam_tilty*0.001/STF_lamb
! ------------


! ------------
! preset Fourier-space with 0-beam of 1.0
  STF_sc(:,:) = 0.0
  STF_sc(1,1) = cmplx(1.0,0.0)
  tpower = 1.0
  STF_PreparedWavePower = tpower
  
! transfer to real space and apply a physe wedge to tilt the wave
  call STF_FFT(STF_sc,ndimx,ndimy,"bac") ! external call from SFFTs.f
  
  do j=1,ndimy
    wy = real(j-ndim2y-1)*samplingy
    do i=1,ndimx
      wx = real(i-ndim2x-1)*samplingx
      chi = 2.0*STF_pi*(wx*lbtx+wy*lbty)
      STF_sc(i,j) = cmplx(cos(chi),sin(chi))
    end do
  end do
!  write(unit=*,fmt=*) "FS-power:",tpower
! ------------

! ------------
! transform to Fourier space
  call STF_FFT(STF_sc,ndimx,ndimy,"for") ! external call from SFFTs.f
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_PreparePlaneWaveFourier: EXIT."
  return

END SUBROUTINE STF_PreparePlaneWaveFourier
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetPhasePlate(p,ndim,sampling)
! function: calculates the aberration function in fourier space
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(FFT_BOUND,FFT_BOUND) : pphase plate ref
!            integer*4 :: ndim : wave result dimension
!            real*4 :: sampling
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3200
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4 :: p(ndim,ndim)
  
  integer*4 :: i, j, ndim2
  real*4 :: wx, wy
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_GetPhasePlate: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_GetPhasePlate: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  STF_dim = ndim
  ndim2 = int(ndim/2)
  STF_sampling = sampling
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
! ------------


! ------------
! preset Fourier-space array with exp[-i*chi]
  do j=1,ndim
    wy = (j-ndim2-1)*STF_itow !(mod((j+ndim2m1),ndim)-ndim2)*STF_itow
    do i=1,ndim
      wx = (i-ndim2-1)*STF_itow !(mod((j+ndim2m1),ndim)-ndim2)*STF_itow
      
      p(i,j) = STF_AberrationFunction(wx,wy)
      
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_GetPhasePlate: EXIT."
  return

END SUBROUTINE STF_GetPhasePlate
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_AberrateWaveFourier(wave,ndimx,ndimy,samplingx,samplingy)
! function: applies aberrations to the given wave function
! -------------------------------------------------------------------- !
! parameter: complex*8 :: wave(ndimy,ndimx)
!            integer*4 :: ndimx,ndimy : wave dimension
!            real*4 :: samplingx,samplingy : wave result sampling x & y
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3600
  
  complex*8 :: wave(ndimy,ndimx)
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy
  
  integer*4 :: i, j, ndim2x, ndim2y
  real*4 :: wx, wy, wx2, chi
  real*4 :: itogx, itogy, itowx, itowy
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_PrepareProbeWaveFourier: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_PrepareProbeWaveFourier: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndimx)
  call STF_SETTAB_SCR2(ndimy)
  !STF_dim = ndim
  ndim2x = int(ndimx/2)
  ndim2y = int(ndimy/2)
  !STF_sampling = sampling
! set Fourier-space sampling
  !STF_itog = 0.5/sampling/real(ndim2)
  !STF_itow = STF_itog*STF_lamb
  itogx = 0.5/samplingx/real(ndim2x)
  itowx = itogx*STF_lamb
  itogy = 0.5/samplingy/real(ndim2y)
  itowy = itogy*STF_lamb
  !wave = cmplx(0.0,0.0)
! ------------


! ------------
  do j=1,ndimx
    wx = STF_TABBED_SCR(j)*itowx
    wx2 = wx*wx
    do i=1,ndimy
      wy = STF_TABBED_SCR2(i)*itowy
      
      chi = STF_AberrationFunction(wx,wy)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsic calculation method
      wave(i,j) = cval*wave(i,j)
      
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_AberrateWaveFourier: EXIT."
  return

END SUBROUTINE STF_AberrateWaveFourier
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_ApplyObjectiveAperture(wave,ndimx,ndimy,samplingx,samplingy,acx, acy)
! function: applies aberrations to the given wave function
! -------------------------------------------------------------------- !
! parameter: complex*8 :: wave(ndimy,ndimx)
!            integer*4 :: ndimx,ndimy : wave dimension
!            real*4 :: samplingx,samplingy : wave result sampling x & y
!            real*4 :: acx, acy : aperture center x & y
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3600
  
  complex*8 :: wave(ndimy,ndimx)
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy, acx, acy
  
  integer*4 :: i, j, ndim2x, ndim2y
  real*4 :: wx, wy, wx2, wap
  real*4 :: itogx, itogy, itowx, itowy
  real*4 :: anglethresh, threshwidth, rval, power
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_ApplyObjectiveAperture: INIT."
! check init-status
  if (STF_maxaberration<=0) then
    call STF_ERROR("STF_ApplyObjectiveAperture: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  call STF_SETTAB_SCR(ndimx)
  call STF_SETTAB_SCR2(ndimy)
  !STF_dim = ndim
  ndim2x = int(ndimx/2)
  ndim2y = int(ndimy/2)
  !STF_sampling = sampling
! set Fourier-space sampling
  !STF_itog = 0.5/sampling/real(ndim2)
  !STF_itow = STF_itog*STF_lamb
  itogx = 0.5/samplingx/real(ndim2x)
  itowx = itogx*STF_lamb
  itogy = 0.5/samplingy/real(ndim2y)
  itowy = itogy*STF_lamb
  !wave = cmplx(0.0,0.0)
  anglethresh = max(STF_caperture*STF_RELANGTHRESH,0.25*(itowx+itowy))
  threshwidth = 0.5*(itowx+itowy)*STF_APSMOOTHPIX
! ------------


! ------------
  do j=1,ndimx
    wx = STF_TABBED_SCR(j)*itowx-acx
    wx2 = wx*wx
    do i=1,ndimy
      wy = STF_TABBED_SCR2(i)*itowy-acy
      wap = wx2 + wy**2
      
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)
      power = 1.0-rval
      
      if (power<STF_APERTURETHRESH) power = 0.0
      
      wave(i,j) = power*wave(i,j)
      
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_ApplyObjectiveAperture: EXIT."
  return

END SUBROUTINE STF_ApplyObjectiveAperture
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION STF_HT2WL(ht)
! function: calculates wavelength [nm] from high-tension [kV]
! -------------------------------------------------------------------- !
! parameter: real*4 :: ht
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3700
  real*4, intent(in) :: ht
  real*4 :: STF_HT2WL
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_HT2WL: INIT."
! ------------

! ------------
  STF_HT2WL = 1.239842447 / sqrt( ht * ( 1022.0 + ht ) )
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_HT2WL: EXIT."
  return

END FUNCTION STF_HT2WL
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
FUNCTION STF_WL2HT(wl)
! function: calculates heigh tension in kV from wavelength
! -------------------------------------------------------------------- !
! parameter: real*4 :: wl
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3800
  real*4, intent(in) :: wl
  real*4 :: STF_WL2HT
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_WL2HT: INIT."
! ------------

! ------------
  STF_WL2HT = ( sqrt( 1.0 + (0.0024263/wl)**2 ) - 1.0 ) * 511.0
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_WL2HT: EXIT."
  return

END FUNCTION STF_WL2HT
!**********************************************************************!



END MODULE STEMfunctions

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
!  integer*4, parameter :: subnum = 3900
!
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > <NAME>: INIT."
! ------------

! ------------
! 
! ------------

! ------------
!  write(unit=*,fmt=*) " > <NAME>: EXIT."
!  return

!END SUBROUTINE <NAME>
!**********************************************************************!