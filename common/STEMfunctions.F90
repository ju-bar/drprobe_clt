!*********************************************************************!
!                                                                     !
!    MODULE STEMfunctions                                             !
!    --------------------                                             !
!                                                                     !
!    Purpose  : Implementation of STEM imaging functions and data     !
!               manipulations of simulation and analysis              !
!    File     : STEMfunctions.f90                                     !
!                                                                     !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                 !
!    Version  :  1.0.0, March 28, 2007                                !
!                                                                     !
!    To Link  : SFFTs.f                                               !
!               BasicFuncs.f90                                        !
!                                                                     !
!*********************************************************************!



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
  USE AberrationFunctions
    
  implicit none
  
  save
    
  ! declare internal data types

  ! accessibility of subroutines or functions
!  private :: STF_***
  private :: STF_ERROR
  private :: STF_TABBED_EXP
  private :: STF_TABBED_SIGMOID
  private :: STF_SETTAB_USC
  private :: STF_SETTAB_SCR
  private :: STF_SETTAB_USC2
  private :: STF_SETTAB_SCR2
  private :: STF_SETAMORPH
  private :: STF_FFT
  !private :: STF_SETAMORPHPHASE
  !private :: STF_SETAMORPHSCATAMP
  !private :: STF_VARYAMORPHPHASE
  integer*4, private :: STF_GETAMOPRHVARIANT
  
  
  !public :: STF_RESETAMORPH
  public :: STF_PrepareProbeWaveFourier
  public :: STF_GetCoherentProbe
  public :: STF_GetCoherentProbePhase
!  public :: STF_GetPSQCohProbe
!  public :: STF_GetPTQCohProbe
  public :: STF_SetupPSCKernel
  public :: STF_GetPSCohProbe
  public :: STF_GetPTCohProbe
!  public :: STF_GetQCohProbe
  public :: STF_GetPCohProbe
  public :: STF_GetCohProbeWave
  public :: STF_GetCohProbeWaveRe
  public :: STF_GetCohProbeWaveIm
  public :: STF_GetPhasePlate
  public :: STF_GetRonchigram
  
!  public :: STF_***
  public :: STF_INIT
  public :: STF_UNINIT
  public :: STF_INIT_FFT
  
  public :: STF_GenerateProbe

!  declare module global (but private) variables, params and arrays
!  file units interval
  integer*4, private, parameter :: STF_minunit = 81
  integer*4, private, parameter :: STF_maxunit = 90

!   length of names and internal strings
  integer*4, public, parameter :: STF_ll = 1024

!   enpty line placeholder (emptied at startup)
  character(len=STF_ll), private :: STF_el  
  save
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
  !integer*4, public, parameter :: FFT_BOUND = 8192
  !integer*4, public, parameter :: FFT_BOUND = 4096
  integer*4, public, parameter :: FFT_BOUND_MIN = 128
  integer*4, public, parameter :: FFT_BOUND_MAX = 8192
  integer*4, public :: FFT_BOUND
  DATA FFT_BOUND /2048/
  integer*4, private :: FFT_NYQ
  DATA FFT_NYQ /1024/

! aperture power threshold
  real*4, private, parameter :: STF_APERTURETHRESH = 1.0E-02
  
! aperture angle thresholds
  real*4, private, parameter :: STF_RELANGTHRESH = 1.0E-3
  real*4, private, parameter :: STF_RELANGWIDTHTHRESH = 2.0E-2
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

! current real-space sampling [nm/pix]
  real*4, public :: STF_sampling
  DATA STF_sampling /0.05/

! current Fourier space sampling [(pix*nm)^-1]
  real*4, public :: STF_itog, STF_itow

! defocus spread
  real*4, public :: STF_defocusspread ! [nm]
  DATA STF_defocusspread /3.0/
  
! source radius
  real*4, public :: STF_srcradius ! [nm]
  DATA STF_srcradius /0.025/
  
! source distribution function
  integer*4, public :: STF_srctype ! 0: ideal point, 1: gaussian, 2: lorentz, 3: disk
  DATA STF_srctype /1/
  
! condenser aperture
  real*4, public :: STF_caperture ! [mrad]
  DATA STF_caperture /30.0/
  real*4, public :: STF_caperture_movex ! [mrad]
  DATA STF_caperture_movex /0.0/
  real*4, public :: STF_caperture_movey ! [mrad]
  DATA STF_caperture_movey /0.0/
  
! beam tilt
  real*4, public :: STF_beam_tiltx ! [mrad]
  DATA STF_beam_tiltx /0.0/
  real*4, public :: STF_beam_tilty ! [mrad]
  DATA STF_beam_tilty /0.0/
  
! data array used for calculation
  complex*8, allocatable, dimension(:,:), public :: sc, sc2
  real*4, public :: STF_PreparedWavePower
  
! anglular function tables
  integer*4, parameter, private :: STF_ANGTAB_SIZE = 2047
  complex*8, dimension(STF_ANGTAB_SIZE), private :: STF_ANGTAB_EXP
! sigmoid function table
  real*4, parameter, private :: STF_SIGMOID_EXT = 3.0
  real*4, dimension(STF_ANGTAB_SIZE), private :: STF_ANGTAB_SIGMOID
! scramble and unscramble lists
  integer*4, allocatable, dimension(:), private :: STF_TABBED_SCR, STF_TABBED_USC  
  integer*4, allocatable, dimension(:), private :: STF_TABBED_SCR2, STF_TABBED_USC2

! amorphous object functions
  complex*8, allocatable, dimension(:,:), private :: STF_amorph_prop
  complex*8, allocatable, dimension(:,:,:,:), private :: STF_amorph_phgr
! amorphous material phase maximum
  real*4, parameter, private :: STF_AMORPH_PHASEMAX = 0.4 ! approx. pi/8
! amorphous material scattering function parameters
  real*4, parameter, private :: STF_AMORPH_SCATANG1 = 3.0 ! [1/nm]
  real*4, parameter, private :: STF_AMORPH_SCATANG2 = 5.0 ! [1/nm]
  real*4, parameter, private :: STF_AMORPH_SCATTAIL = 0.2 ! relative power of large angle tails
  real*4, parameter, private :: STF_AMORPH_DWPRM = 0.08 ! Vibration amplitude parameter for Debye-Waller factors
! relative scattering strength parameter for amorph object
  real*4, public :: STF_amorph_relstr
  DATA STF_amorph_relstr /1.0/
! relative variation parameter for amorph object
  real*4, public :: STF_amorph_relvary, STF_amorph_relvarybk
  DATA STF_amorph_relvary /0.05/
  DATA STF_amorph_relvarybk /0.0/
! number of amorphous variants
  integer*4, parameter, private :: STF_AMORPH_VARIANTNUM = 21
! number of used amorphous variants
  real*4, public :: STF_amorph_variants
  DATA STF_amorph_variants /0/
! total thickness of amorphous object
  real*4, public :: STF_AMORPH_THICKTOT
  DATA STF_AMORPH_THICKTOT /0.0/ ! [nm]
! number of object slices
  integer*4, public :: STF_AMORPH_NSLC
  DATA STF_AMORPH_NSLC /1/  ! int
! backup parameters for amorphous functions
  real*4, private :: STF_amorph_sampling
  DATA STF_amorph_sampling /0.0/
  integer*4, private :: STF_amorph_ndim
  DATA STF_amorph_ndim /0/
  integer*4, private :: STF_amorph_nslcbk
  DATA STF_amorph_nslcbk /0/
  real*4, private :: STF_amorph_dz
  DATA STF_amorph_dz /0.0/
  real*4, private :: STF_amorph_sc1, STF_amorph_sc2
  DATA STF_amorph_sc1 /0.0/
  DATA STF_amorph_sc2 /0.0/
!  integer*4, private :: STF_amorph_newphase
!  DATA STF_amorph_newphase /1/
!  integer*4, private :: STF_amorph_newscamp
!  DATA STF_amorph_newscamp /1/
  
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
  character(len=STF_FUNC_NAME_LENGTH), dimension(STF_FUNC_NUM), public :: STF_FuncNames
  
  character(len=STF_ll), public :: STF_LastErrorMessage
  
! status flag
  integer*4, public :: STF_status
  DATA STF_status /0/



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
  integer*4, parameter :: fnl = STF_FUNC_NAME_LENGTH
  
!
  integer*4 :: m
  real*4 :: angle, ascale
  real*4, external :: sigmoid ! BasicFuncs.f90
  external :: InitRand ! random.f90
! ------------

! ------------
! INIT
!  write(*,*) " > STF_INIT: INIT."
! ------------


! ------------
  STF_status = 0
  STF_pi = atan(1.0)*4.0
  STF_rd2r = STF_pi/180.0
  STF_err_num = 0
  STF_el = REPEAT(" ",STF_ll)
  AF_LastErrorMessage = ""
  call InitRand()
! ------------


! ------------
  STF_DEFOCUS_KERNEL_STEPS = STF_DEFOCUS_KERNEL_STEPS_DEFAULT
  STF_DEFOCUS_KERNEL_SPREAD = STF_DEFOCUS_KERNEL_SPREAD_DEFAULT
! ------------


! ------------
  STF_defocusspread = 0.0
  STF_srcradius = 0.0
  STF_caperture = 30.0
! ------------


! ------------
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
! ------------
 

! ------------
  STF_FuncNames( 1) = "Coherent Probe Intensity"
  STF_FuncNames( 2) = "Coherent Probe Phase"
  STF_FuncNames( 3) = "Partially Spatial Coherent Probe Intensity"
  STF_FuncNames( 4) = "Partially Temporal Coherent Probe Intensity"
  STF_FuncNames( 5) = "Partially Coherent Probe Intensity"
  STF_FuncNames( 6) = "Real Part of Probe Wave Function"
  STF_FuncNames( 7) = "Imaginary Part of Probe Wave Function"
  STF_FuncNames( 8) = "Phase Plate"
  STF_FuncNames( 9) = "Ronchigram, weak scattering"
  STF_FuncNames(10) = "Ronchigram, strong scattering"
! ------------

! ------------
!  call STF_SETAMORPHPHASE()
  STF_amorph_sampling = 0.0
  STF_amorph_ndim = 0
  STF_amorph_variants = 0
  STF_amorph_nslcbk = 0
! ------------

! ------------
  if (STF_err_num==0) STF_status = 1
! ------------

! ------------
!  write(*,*) " > STF_INIT: EXIT."
  return

END SUBROUTINE STF_INIT
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_INIT_FFT(ndim, nerr)
! function: initializes module square arrays for FFT handling
! -------------------------------------------------------------------- !
! parameter: integer*4 :: ndim = the largest size of the input data
!                                input the max of ndimx and ndimy
!            integer*4 :: nerr = error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 150
  integer*4, intent(in) :: ndim
  integer*4, intent(out) :: nerr
  
!
  integer*4 :: nec    ! error count
  integer*4 :: nft    ! target FFT size
  integer*4 :: nalloc ! allocation status
! ------------

! ------------
! INIT
!  write(*,*) " > STF_INIT_FFT: INIT."
  nerr = 0
  nalloc = 0
  if (STF_status==0) then
    nerr = 1
    call STF_ERROR("Module not initialized, STF_INIT_FFT failed.", subnum+nerr)
    return
  end if
  if (ndim<2) then
    nerr = 2
    call STF_ERROR("Invalid input (ndim<2), STF_INIT_FFT failed.", subnum+nerr)
    return
  end if
! ------------

! ------------
! determine the next best FFT size
  nft = 2**CEILING( LOG( real( ndim ) )/LOG(2.0) ) ! next 2^N above max(nx,ny)
  ! apply lower limit
  nft = MAX( FFT_BOUND_MIN, nft )
  ! check upper limit
  if (nft > FFT_BOUND_MAX) then
    nerr = 3
    call STF_ERROR("Invalid FFT size (max. 8192), STF_INIT_FFT failed.", subnum+nerr)
    return
  end if
  ! check previous allocations
  if (nft/=FFT_BOUND) then
    nec = STF_err_num
    if (allocated(sc)) deallocate(sc, stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
    if (allocated(sc2)) deallocate(sc2, stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
    if (allocated(STF_TABBED_SCR)) deallocate(STF_TABBED_SCR, stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
    if (allocated(STF_TABBED_SCR2)) deallocate(STF_TABBED_SCR2, stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
    if (allocated(STF_TABBED_USC)) deallocate(STF_TABBED_USC, stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
    if (allocated(STF_TABBED_USC2)) deallocate(STF_TABBED_USC2, stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
    if (STF_err_num>nec) then
      nerr = 4
      call STF_ERROR("Memory deallocation in STF_INIT_FFT failed.", subnum+nerr)
      return
    end if
  end if
  ! set dependent sizes
  FFT_BOUND = nft
  FFT_NYQ = FFT_BOUND/2
! ------------

! ------------
! allocate arrays
  nec = STF_err_num
  if (.not.allocated(sc)) then
    allocate(sc(FFT_BOUND,FFT_BOUND), stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
  end if
  if (.not.allocated(sc2)) then
    allocate(sc2(FFT_BOUND,FFT_BOUND), stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
  end if
  if (.not.allocated(STF_TABBED_SCR)) then
    allocate(STF_TABBED_SCR(FFT_BOUND), stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
  end if
  if (.not.allocated(STF_TABBED_SCR2)) then
    allocate(STF_TABBED_SCR2(FFT_BOUND), stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
  end if
  if (.not.allocated(STF_TABBED_USC)) then
    allocate(STF_TABBED_USC(FFT_BOUND), stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
  end if
  if (.not.allocated(STF_TABBED_USC2)) then
    allocate(STF_TABBED_USC2(FFT_BOUND), stat=nalloc)
    if (nalloc/=0) STF_err_num = STF_err_num + 1
  end if
  if (STF_err_num>nec) then
    nerr = 5
    call STF_ERROR("Memory allocation in STF_INIT_FFT failed.", subnum+nerr)
    return
  end if
! ------------

! ------------
!  write(*,*) " > STF_INIT_FFT: EXIT."
  return

END SUBROUTINE STF_INIT_FFT
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
  integer*4 :: err, m
! ------------

! ------------
! INIT
!  write(*,*) " > STF_INIT: STF_UNINIT."
! ------------
  err = STF_err_num
!  if (allocated(tmpdata)) deallocate(tmpdata,stat=m)
  !if (m/=0) STF_err_num = STF_err_num + 1
  if (allocated(sc)) deallocate(sc,stat=m)
  if (m/=0) STF_err_num = STF_err_num + 1
  if (allocated(sc2)) deallocate(sc2,stat=m)
  if (m/=0) STF_err_num = STF_err_num + 1
  if (allocated(STF_amorph_prop)) deallocate(STF_amorph_prop,stat=m)
  if (m/=0) STF_err_num = STF_err_num + 1
  if (allocated(STF_amorph_phgr)) deallocate(STF_amorph_phgr,stat=m)
  if (m/=0) STF_err_num = STF_err_num + 1
! ------------
  err = 0
  STF_status = 0
! ------------

! ------------
!  write(*,*) " > STF_UNINIT: EXIT."
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
  write(unit=sinfo,fmt=*) trim(sTxt)," Error code:",nErr
  AF_LastErrorMessage = trim(sinfo)
!  call SE_event(trim(sinfo), SE_err)
  write(*,*) " > SFT_ERROR: "//trim(sinfo)
! ----------

  STF_err_num = STF_err_num + 1
  return

END SUBROUTINE STF_ERROR
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
!  write(*,*) " > STF_TABBED_EXP: INIT."
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
!  write(*,*) " > STF_TABBED_EXP: EXIT."
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
!  write(*,*) " > STF_TABBED_EXP: INIT."
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
! clip angle to interval 1 .. STF_ANGTAB_SIZE
  if (idx<1) then
    idx = 1
  end if
  if (idx>STF_ANGTAB_SIZE) then
    idx = STF_ANGTAB_SIZE
  end if
! ------------

! ------------
! get data
  rval = STF_ANGTAB_SIGMOID(idx)
! ------------

! ------------
!  write(*,*) " > STF_TABBED_EXP: EXIT."
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
!  write(*,*) " > STF_SETTAB_SCR: INIT."
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
!  write(*,*) " > STF_SETTAB_SCR: EXIT."
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
!  write(*,*) " > STF_SETTAB_USC: INIT."
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
!  write(*,*) " > STF_SETTAB_USC: EXIT."
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
!  write(*,*) " > STF_SETTAB_SCR2: INIT."
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
!  write(*,*) " > STF_SETTAB_SCR2: EXIT."
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
!  write(*,*) " > STF_SETTAB_USC2: INIT."
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
!  write(*,*) " > STF_SETTAB_USC2: EXIT."
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
  integer*4, parameter :: subnum = 4000
  integer*4 :: nx,ny,transformed
  character*(*) :: dir
  complex*8 :: cdata(FFT_BOUND,FFT_BOUND)   
  character(len=400) :: sinfo
  external :: ODDCC128S, ODDCC256S, ODDCC512S, ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STF_FFT: INIT."
  transformed = 0
! ------------

! ------------
  select case (FFT_BOUND)
  
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
  
  end select ! (FFT_BOUND)
! ------------

! ------------
  if (transformed/=1) then
    write(unit=sinfo,fmt='(A,I5,A)') &
     & "Failed to call FFT routines, data size (",FFT_BOUND,") not supported."
    call STF_ERROR(trim(sinfo), subnum+1)
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > STF_FFT: EXIT."
  return

END SUBROUTINE STF_FFT
!**********************************************************************!







!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SETAMORPH(nerr)
! function: calculates an amorphous phase grating from random object
!           and scattering amplitudes using the current setup
!           STF_dim
! -------------------------------------------------------------------- !
! parameter:
!    IN/OUT: integer*4 :: nerr      ! error code
! -------------------------------------------------------------------- !
! link: SFFTs.f
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3700
  integer*4, intent(inout) :: nerr
  integer*4 :: i, j, k, l, ndim, nslc, nv, nalloc, nyq
  integer*4 :: nrealloc, ndealloc, newphgr, newprop
  integer*4 :: nvarypix
  real*4 :: gx, gy, gx2, g2, itog
  real*4 :: phase, dz, rstr, rvary
  real*4 :: dwf, sca, fd
  real*4 :: rndpha_range
  complex*8 :: czero, cval
  real*4, allocatable :: rndpha(:,:) ! random phase data
  real*4, allocatable :: rndphavar(:,:) ! random phase variant data
  real*4, allocatable :: sctamp(:,:) ! scattering amplitude data
  complex*8, allocatable :: objdat(:,:) ! object data
  
  !external :: ODDCC2048S
  real*4, external :: UniRand 
! ------------

! ------------
! INIT
!  write(*,*) " > STF_SETAMORPH: INIT."
  nerr = 0
  czero = cmplx(0.0,0.0)
  ndim = STF_dim
  nyq = ndim/2
  nslc = STF_AMORPH_NSLC
  nv = 1
  if (STF_amorph_relvary>0.0) nv = STF_AMORPH_VARIANTNUM
  nrealloc = 1
  ndealloc = 0
  newphgr = 0
  newprop = 0
! ------------

! ------------
! check need for allocation
  if (ALLOCATED(STF_amorph_prop).and.ALLOCATED(STF_amorph_phgr)) then ! array are already allocated
    if (     (ndim==STF_amorph_ndim)  .and. (nslc==STF_amorph_nslcbk) &
     &  .and.(nv==STF_amorph_variants) ) then ! arrays have already the correct size
      ndealloc = 0 ! no need to deallocate
      nrealloc = 0 ! no need to reallocate
    else
      ndealloc = 1 ! we need to deallocate
      nrealloc = 1 ! .. and reallocate
    end if
  else
    ndealloc = 0 ! ... no need for deallocation
    nrealloc = 1 ! we need to allocate
  end if
! REMARK: ndealloc can only be =1 if nrealloc is ==1
!
! do the allocation if required
  if (ndealloc==1) then
    deallocate(STF_amorph_prop, STF_amorph_phgr, stat=nalloc)
    if (nalloc/=0) then
      nerr = 1
      call STF_ERROR("Deallocation of amorphous object data arrays failed.",subnum+nerr)
      return
    end if
  end if
  if (nrealloc==1) then
    allocate(STF_amorph_prop(ndim,ndim), STF_amorph_phgr(ndim,ndim,nv,nslc), stat=nalloc)
    if (nalloc/=0) then
      nerr = 2
      call STF_ERROR("Allocation of amorphous object data arrays failed.",subnum+nerr)
      return
    end if
    newphgr = 1 ! recalculate the phasegratings
    if (nslc>1) newprop = 1 ! recalculate the propagator
  end if
! ------------


! ------------
! update local object data
  dz = 0.0
  rstr = STF_amorph_relstr
  rvary = STF_amorph_relvary
  sca = STF_AMORPH_SCATANG1*rstr
  dwf = STF_AMORPH_DWPRM/sqrt(rstr)
  fd = 1.0/sca**2
  if (nslc>1) dz = STF_AMORPH_THICKTOT / real(nslc-1)
  call STF_SETTAB_SCR(ndim)
  itog = 0.5/STF_sampling/real(nyq) ! Fourier-space sampling
  STF_itog = itog
  !sa = STF_AMORPH_SCATTAIL
  !ta = 1.0/(1+sa)
  
  
  
  
  rndpha_range = 2.0*STF_pi
!
! check for other reasons to update the phasegrating data
! (different scattering amplitudes or sampling)
  if (     (STF_amorph_sampling/=STF_sampling) &
     & .or.(STF_amorph_sc1/=sca) &
     & .or.(STF_amorph_relvarybk/=rvary) &
     & .or.(STF_amorph_nslcbk/=nslc) ) then
    newphgr = newphgr + 1 ! recalculate the phasegratings
  end if
!
! check for other reasons to update the propagator data
! (different slice size or sampling)
  if (     (STF_amorph_ndim/=ndim) &
     & .or.(STF_amorph_dz/=dz) &
     & .or.(STF_amorph_sampling/=STF_sampling) ) then
    newprop = newprop + 1 ! recalculate the phasegratings
  end if
! ------------


! ------------
! CALCULATION OF PHASE GRATINGS AND PROPAGATORS
! This is a multi-step routine:
! - preliminary steps
! - preparation of scattering amplitudes and propgator in fourier space, (sctamp) and STF_amorph_prop
! - preparation of nslc random phase data arrays (real*4) in real space (rndpha)
! - - creation of nv variants of each of the nslc data arrays
! - - - transformation of the phase data to fourier space
! - - - application of the scattering amplitudes
! - - - transformation back to real-space
! - - - calculation of the phase object - phase grating
! - - - storing in STF_amorph_phgr
! - -
! - post calculation steps
! 
!
! - preliminary steps
  if ( (newprop/=0) .or. (newphgr/=0) ) then
    allocate(sctamp(ndim,ndim), rndpha(ndim,ndim), &
     &       rndphavar(ndim,ndim), objdat(FFT_BOUND,FFT_BOUND), stat=nalloc)
    if (nalloc/=0) then
      nerr = 3
      call STF_ERROR("Allocation of temporary data arrays failed.",subnum+nerr)
      return
    end if
  end if
!
!
! - preparation of scattering amplitudes and propgator in fourier space, (sctamp) and STF_amorph_prop
  if ( (newprop/=0) .or. (newphgr/=0) ) then
    sctamp(:,:) = 0.0
    do j=1,ndim
      gx = STF_TABBED_SCR(j)*itog
      gx2 = gx*gx
      do i=1,ndim
        gy = STF_TABBED_SCR(i)*itog
        g2 = gx2+gy*gy
        phase = STF_AMORPH_PHASEMAX/(g2*fd+1.0)*exp(-dwf*g2) ! (max. amplitude)*(lorentz of width sc)*(dwf)
        sctamp(i,j) = phase
        if (nslc>1) then
          phase = STF_pi*dz*g2*STF_lamb
          cval = cmplx(0.0,-phase)
          STF_amorph_prop(i,j) = exp( cval )
        end if ! (nslc>1) 
      end do
    end do  
  end if
!
!
! *** begin of pure phase grating actions
  if (newphgr/=0) then
!
!   individual phases for every slice
    do k=1, nslc
!   
! - preparation of nslc random phase data arrays (real*4) in real space (rndpha)
!   the random phase range is 0 ... 2*Pi
      do j=1, ndim
        do i=1, ndim
          rndpha(i,j) = rndpha_range*UniRand()
        end do 
      end do
!
!
! - - creation of nv variants of each of the nslc data arrays
      do l=1, nv
        rndphavar = rndpha ! copy current offset
        nvarypix = int(real(ndim*ndim)*rvary)
        do
          if (nvarypix<=0) exit
          i = modulo( int( real(ndim-1)*UniRand() ), ndim ) + 1
          j = modulo( int( real(ndim-1)*UniRand() ), ndim ) + 1
          rndphavar(i,j) = rndpha_range*UniRand()
          nvarypix = nvarypix - 1
        end do
!
!
! - - - transformation of the phase data to fourier space
        do j=1, ndim
          do i=1, ndim
            objdat(i,j) = rndphavar(i,j)
          end do 
        end do
        call STF_FFT(objdat,ndim,ndim,'forward')
!
!
! - - - application of the scattering amplitudes
        do j=1, ndim
          do i=1, ndim
            objdat(i,j) = objdat(i,j)*sctamp(i,j)
          end do 
        end do
!
!
! - - - transformation back to real-space
        call STF_FFT(objdat,ndim,ndim,'backward')
!
!
! - - - calculation of the phase object - phase grating
! - - - and storing in STF_amorph_phgr
        do j=1, ndim
          do i=1, ndim
            cval = cexp( cmplx( 0.0, real(objdat(i,j)) ) )
            STF_amorph_phgr(i,j,l,k) = cval
          end do 
        end do
!
!
      end do ! l=1, nv
! - -
    end do ! k=1, nslc
!
! *** end of pure phase grating actions
  end if !(newphgr/=0)
!
!
! - post calculation steps
  if ( allocated(sctamp) ) then
    deallocate(sctamp, stat=nalloc)
    if (nalloc/=0) then
      nerr = 4
      call STF_ERROR("Deallocation of temporary data arrays failed.",subnum+nerr)
      return
    end if
  end if
  if ( allocated(rndpha) ) then
    deallocate(rndpha, stat=nalloc)
    if (nalloc/=0) then
      nerr = 5
      call STF_ERROR("Deallocation of temporary data arrays failed.",subnum+nerr)
      return
    end if
  end if
  if ( allocated(rndphavar) ) then
    deallocate(rndphavar, stat=nalloc)
    if (nalloc/=0) then
      nerr = 6
      call STF_ERROR("Deallocation of temporary data arrays failed.",subnum+nerr)
      return
    end if
  end if
  if ( allocated(objdat) ) then
    deallocate(objdat, stat=nalloc)
    if (nalloc/=0) then
      nerr = 7
      call STF_ERROR("Deallocation of temporary data arrays failed.",subnum+nerr)
      return
    end if
  end if
! ------------


! ------------
! store calculation parameters
  STF_amorph_ndim = STF_dim
  STF_amorph_nslcbk = nslc
  STF_amorph_variants = nv
  STF_amorph_relvarybk = rvary
  STF_amorph_sampling = STF_sampling
  STF_amorph_sc1 = sca
  STF_amorph_sc2 = 0.0
  STF_amorph_dz = dz
! ------------  


! ------------
!  write(*,*) " > STF_SETAMORPH: EXIT."
  return

END SUBROUTINE STF_SETAMORPH
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
integer*4 FUNCTION STF_GETAMORPHVARIANT()
! function: returns the number of a random amorphous variant
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3900
  integer*4 :: nv, nr
  real*4, external :: UniRand
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GETAMORPHVARIANT: INIT."
  nv = STF_amorph_variants
  nr = 1
! ------------

! ------------
  if (nv>1) then
    nr = modulo( int(real(nv)*UniRand()), nv)+1
  end if
  STF_GETAMORPHVARIANT = nr
! ------------

! ------------
!  write(*,*) " > STF_GETAMORPHVARIANT: EXIT."
  return

END FUNCTION STF_GETAMORPHVARIANT
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
SUBROUTINE STF_GetCoherentProbe(p,ndim,sampling)
! function: calculates the coherent probefunction from the current
!           coherent aberrations in real-space
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!                                     probefunction result
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !


  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 900
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4, intent(out) :: p(ndim,ndim)
  
  integer*4 :: i, j, ndim2, i1, j1
  real*4 :: wx, wy, wy2, wap, chi, power, anglethresh, threshwidth
  real*4 :: tpower, tpowscal, rval, camx, camy, lbtx, lbty
  complex*8 :: cval
  
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetCoherentProbe: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetCoherentProbe: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  STF_dim = ndim
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  ndim2 = int(ndim/2)
  STF_sampling = sampling
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  sc = cmplx(0.0,0.0)
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  tpower = 0.0
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------


! ------------
! preset Fourier-space array with exp[-i*chi]
  tpower = 0.0
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow
    wy2 = (wx-camx)*(wx-camx)
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      wap = wy2+(wy-camy)*(wy-camy)
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx+lbtx,wy+lbty) ! get aberration function with tilt offset
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      sc(i,j) = cval*power
!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
    
    end do
  end do
!  write(*,*) "FS-power:",tpower
!  call STF_PrepareProbeWaveFourier(ndim,sampling)
!  tpower = STF_PreparedWavePower
!  tpowscal = 1.0/real(ndim*ndim)/tpower
! ------------



! ------------
! transform to real space
  call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
! ------------

! ------------
! calculate probe by taking absolute square in real space
! + unscramble (shifts origin)
  tpower = 0.0
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
      rval = real( cval * conjg(cval) )
      tpower = tpower + rval
      p(i1,j1) = rval
    end do
  end do
  if (tpower>0.0) p = p / tpower
! ------------

! ------------
!  write(*,*) " > STF_GetCoherentProbe: EXIT."
  return

END SUBROUTINE STF_GetCoherentProbe
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetCoherentProbePhase(p,ndim,sampling)
! function: calculates the coherent probefunction from the current
!           coherent aberrations in real-space
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!                                     probefunction result
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !


  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 930
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4, intent(out) :: p(ndim,ndim)
  
  integer*4 :: i, j, ndim2, i1, j1
  real*4 :: wx, wy, wy2, wap, chi, power, anglethresh, threshwidth
  real*4 :: tpower, tpowscal, rval, camx, camy, lbtx, lbty
  complex*8 :: cval
  
  !external :: ODDCC2048S
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetCoherentProbePhase: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetCoherentProbePhase: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  STF_dim = ndim
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  ndim2 = int(ndim/2)
  STF_sampling = sampling
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  sc = cmplx(0.0,0.0)
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  tpower = 0.0
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------


! ------------
! preset Fourier-space array with exp[-i*chi]
  tpower = 0.0
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow
    wy2 = (wx-camx)*(wx-camx)
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      wap = wy2+(wy-camy)*(wy-camy)
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx+lbtx,wy+lbty) ! get aberration function with tilt offset
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      sc(i,j) = cval*power
!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
    
    end do
  end do
!  write(*,*) "FS-power:",tpower
!  call STF_PrepareProbeWaveFourier(ndim,sampling)
!  tpower = STF_PreparedWavePower
!  tpowscal = 1.0/real(ndim*ndim)/tpower
! ------------



! ------------
! transform to real space
  call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
! ------------

! ------------
! calculate probe phase in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
      p(i1,j1) = atan2(imag(cval),real(cval))
    end do
  end do
! ------------

! ------------
!  write(*,*) " > STF_GetCoherentProbePhase: EXIT."
  return

END SUBROUTINE STF_GetCoherentProbePhase
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetCoherentProbeFPS(p,ndim,sampling)
! function: calculates the coherent probefunction from the current
!           coherent aberrations in real-space and transforms then
!           to fourier space - amplitude
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!                                     probefunction result
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !


  implicit none
  
! ------------
! DECLARATION
  integer*4, parameter :: subnum = 970
  
  integer*4 :: ndim
  real*4 :: sampling
  real*4, intent(out) :: p(ndim,ndim)
  
  integer*4 :: i, j, ndim2, i1, j1
  real*4 :: powscl
!  real*4 :: o(ndim,ndim)
  complex*8 :: cval
  
  !external :: ODDCC2048S
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetCoherentProbeFPS: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetCoherentProbeFPS: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  STF_dim = ndim
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  ndim2 = int(ndim/2)
! set Fourier-space sampling
  sc = cmplx(0.0,0.0)
  powscl = real(ndim*ndim)**2
! ------------


! ------------
! get real space probe intensity
  call STF_GetCoherentProbe(p,ndim,sampling)
! ------------


! ------------
! transfer to Fourier space
! set complex array
  do j=1,ndim
    do i=1,ndim
      cval = cmplx(p(i,j),0.0)
      sc(i,j) = cval
    end do
  end do
! transform
  call STF_FFT(sc,ndim,ndim,"for") ! external call from SFFTs.f
! collect data
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
      p(j1,i1) = real( cval * conjg(cval) ) !*powscl
    end do
  end do
  powscl = sum(p)
  if (powscl>0.0) p = p / powscl
! ------------


! ------------
!  write(*,*) " > STF_GetCoherentProbeFPS: EXIT."
  return

END SUBROUTINE STF_GetCoherentProbeFPS
!**********************************************************************!












!! -------> QUASI COHERENT THINGS
!! -------> QUASI COHERENT THINGS
!! -------> QUASI COHERENT THINGS
!
!
!
!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE STF_GetPSQCohProbe(p,ndim,sampling,srcradius)
!! function: calculates the probefunction including partial spatial
!!           coherence from the current coherent aberrations
!!           and from the demagnified source radius in quasi-coherent
!!           formulation
!! -------------------------------------------------------------------- !
!! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!!                                     probefunction result
!!            integer*4 :: ndim : size of array
!!            real*4 :: sampling : real-space sampling for probefunction
!!            real*4 :: srcradius : demagnified source radius [nm]
!! -------------------------------------------------------------------- !
!! required link: SFFTs.f
!! -------------------------------------------------------------------- !
!
!
!  implicit none
!
!! ------------
!! DECLARATION
!  integer*4, parameter :: subnum = 1000
!  
!  integer*4 :: ndim
!  real*4 :: sampling, srcradius
!  real*4, intent(out) :: p(ndim,ndim)
!  
!  integer*4 :: i, j, ndim2, i1, j1, ndim2m1
!  real*4 :: w2, wx, wy, wap, chi, power, anglethresh, threshwidth, wy2, wyap
!  real*4 :: tpower, tpowscal, PSfac, rval, camx, camy, lbtx, lbty
!  complex*8 :: cval
!  
!  !external :: ODDCC2048S
!  real*4, external :: sigmoid
!! ------------
!
!! ------------
!! INIT
!!  write(*,*) " > STF_GetPSQCohProbe: INIT."
!! check init-status
!  if (STF_status<=0.or.AF_maxaberration<=0) then
!    call STF_ERROR("STF_GetPSQCohProbe: Module not initialized. Aborting.",subnum+1)
!    return
!  end if
!! save data size
!  STF_dim = ndim
!  call STF_SETTAB_SCR(ndim)
!  call STF_SETTAB_USC(ndim)
!  ndim2 = int(ndim/2)
!  ndim2m1 = ndim2-1
!  STF_sampling = sampling
!! set Fourier-space sampling
!  STF_itog = 0.5/sampling/real(ndim2)
!  STF_itow = STF_itog*STF_lamb
!  STF_srcradius = srcradius
!  sc = cmplx(0.0,0.0)
!  anglethresh = STF_caperture*STF_RELANGTHRESH
!  threshwidth = STF_itow*STF_APSMOOTHPIX
!  PSfac = STF_pi*STF_srcradius/STF_lamb
!  PSfac = -0.5*PSfac*PSfac
!  tpower = 0.0
!  camx = STF_caperture_movex*0.001
!  camy = STF_caperture_movey*0.001
!  lbtx = STF_beam_tiltx*0.001
!  lbty = STF_beam_tilty*0.001
!! ------------
!
!! ------------
!! preset Fourier-space array with exp[-i*chi]
!  tpower = 0.0
!  do j=1,ndim
!    wx = STF_TABBED_SCR(j)*STF_itow + lbtx
!    wy2 = wx*wx
!    wyap = (wx-camx)*(wx-camx)
!    do i=1,ndim
!      wy = STF_TABBED_SCR(i)*STF_itow + lbty
!      wap = (wy-camy)*(wy-camy)+wyap
!      w2 = wy2+wy*wy
!      
!      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
!      power = 1.0-rval
!      
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      
!      chi = AF_AberrationFunction(wx,wy)
!      power = power * Exp(PSfac*w2)
!      
!      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      sc(i,j) = cval*power
!      tpower = tpower + power*power
!      
!    end do
!  end do
!!  write(*,*) "FS-power:",tpower
!  tpowscal = 1.0/real(ndim*ndim)/tpower
!! ------------
!
!! ------------
!! transform to real space (scramble first!)
!  call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
!! ------------
!
!! ------------
!! calculate probe by taking absolute square in real space
!! + unscramble (shifts origin)
!  do j=1,ndim
!    j1 = STF_TABBED_USC(j)
!    do i=1,ndim
!      i1 = STF_TABBED_USC(i)
!      cval = sc(i,j)
!      p(i1,j1) = real( cval * conjg(cval) )*tpowscal
!    end do
!  end do
!! ------------
!
!! ------------
!!  write(*,*) " > STF_GetPSQCohProbe: EXIT."
!  return
!
!END SUBROUTINE STF_GetPSQCohProbe
!!**********************************************************************!
!
!
!
!
!
!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE STF_GetPTQCohProbe(p,ndim,sampling,defocusspread)
!! function: calculates the probefunction including partial temporal
!!           coherence from the current coherent aberrations
!!           and from the defocus-spread parameter
!! -------------------------------------------------------------------- !
!! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!!                                     probefunction result
!!            integer*4 :: ndim : size of array
!!            real*4 :: sampling : real-space sampling for probefunction
!!            real*4 :: defocusspread : defocus spread [nm]
!! -------------------------------------------------------------------- !
!! required link: SFFTs.f
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! DECLARATION
!  integer*4, parameter :: subnum = 1200
!  
!  integer*4 :: ndim
!  real*4 :: sampling, defocusspread
!  real*4, intent(out) :: p(ndim,ndim)
!  
!  integer*4 :: i, j, ndim2, i1, j1, ndim2m1
!  real*4 :: w2, wx, wy, chi, power, anglethresh, threshwidth, wy2, wyap, wap
!  real*4 :: tpower, tpowscal, PTfac, rval, camx, camy, lbtx, lbty
!  complex*8 :: cval
!  
!  !external :: ODDCC2048S
!  real*4, external :: sigmoid
!! ------------
!
!! ------------
!! INIT
!!  write(*,*) " > STF_GetPTQCohProbe: INIT."
!! check init-status
!  if (STF_status<=0.or.AF_maxaberration<=0) then
!    call STF_ERROR("STF_GetPTQCohProbe: Module not initialized. Aborting.",subnum+1)
!    return
!  end if
!! save data size
!  STF_dim = ndim
!  call STF_SETTAB_SCR(ndim)
!  call STF_SETTAB_USC(ndim)
!  ndim2 = int(ndim/2)
!  ndim2m1 = ndim2-1
!  STF_sampling = sampling
!! set Fourier-space sampling
!  STF_itog = 0.5/sampling/real(ndim2)
!  STF_itow = STF_itog*STF_lamb
!  STF_defocusspread = defocusspread
!  sc = cmplx(0.0,0.0)
!  anglethresh = STF_caperture*STF_RELANGTHRESH
!  threshwidth = STF_itow*STF_APSMOOTHPIX
!  PTfac = STF_pi*STF_defocusspread/STF_lamb
!  PTfac = -0.25*PTfac*PTfac
!  tpower = 0.0
!  camx = STF_caperture_movex*0.001
!  camy = STF_caperture_movey*0.001
!  lbtx = STF_beam_tiltx*0.001
!  lbty = STF_beam_tilty*0.001
!! ------------
!
!! ------------
!! preset Fourier-space array with exp[-i*chi]
!  tpower = 0.0
!  do j=1,ndim
!    wx = STF_TABBED_SCR(j)*STF_itow + lbtx
!    wy2 = wx*wx
!    wyap = (wx-camx)*(wx-camx)
!    do i=1,ndim
!      wy = STF_TABBED_SCR(i)*STF_itow + lbty
!      w2 = wy*wy+wy2
!      wap = wyap+(wy-camy)*(wy-camy)
!
!      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
!      power = 1.0-rval
!      
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      
!      chi = AF_AberrationFunction(wx,wy)
!      power = power * Exp(PTfac*w2*w2)
!      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      sc(i,j) = cval*power
!
!      tpower = tpower + power*power
!
!    end do
!  end do
!!  write(*,*) "FS-power:",tpower
!  tpowscal = 1.0/real(ndim*ndim)/tpower
!! ------------
!
!! ------------
!! transform to real space (scramble first!)
!  call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
!! ------------
!
!! ------------
!! calculate probe by taking absolute square in real space
!! + unscramble (shifts origin)
!  do j=1,ndim
!    j1 = STF_TABBED_USC(j)
!    do i=1,ndim
!      i1 = STF_TABBED_USC(i)
!      cval = sc(i,j)
!      p(i1,j1) = real( cval * conjg(cval) )*tpowscal
!    end do
!  end do
!! ------------
!
!! ------------
!!  write(*,*) " > STF_GetPTQCohProbe: EXIT."
!  return
!
!END SUBROUTINE STF_GetPTQCohProbe
!!**********************************************************************!
!
!
!
!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE STF_GetQCohProbe(p,ndim,sampling,srcradius,defocusspread)
!! function: calculates the probefunction including partial spatial
!!           and partial temporal coherence from the current coherent
!!           aberrations and from the demagnified source radius
!!           in quasi-coherent formulation
!! -------------------------------------------------------------------- !
!! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!!                                     probefunction result
!!            integer*4 :: ndim : size of array
!!            real*4 :: sampling : real-space sampling for probefunction
!!            real*4 :: srcradius : demagnified source radius [nm]
!!            real*4 :: defocusspread : defocus spread [nm]
!! -------------------------------------------------------------------- !
!! required link: SFFTs.f
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! DECLARATION
!  integer*4, parameter :: subnum = 1400
!  
!  integer*4 :: ndim
!  real*4 :: sampling, srcradius, defocusspread
!  real*4, intent(out) :: p(ndim,ndim)
!  
!  integer*4 :: i, j, ndim2, i1, j1, ndim2m1
!  real*4 :: w2, wx, wy, chi, power, anglethresh, threshwidth, wy2, wyap, wap
!  real*4 :: tpower, tpowscal, PSfac, PTfac, rval, camx, camy, lbtx, lbty
!  complex*8 :: cval
!  
!  !external :: ODDCC2048S
!  real*4, external :: sigmoid
!! ------------
!
!! ------------
!! INIT
!!  write(*,*) " > STF_GetQCohProbe: INIT."
!! check init-status
!  if (STF_status<=0.or.AF_maxaberration<=0) then
!    call STF_ERROR("STF_GetQCohProbe: Module not initialized. Aborting.",subnum+1)
!    return
!  end if
!! save data size
!  STF_dim = ndim
!  call STF_SETTAB_SCR(ndim)
!  call STF_SETTAB_USC(ndim)
!  ndim2 = int(ndim/2)
!  ndim2m1 = ndim2-1
!  STF_sampling = sampling
!! set Fourier-space sampling
!  STF_itog = 0.5/sampling/real(ndim2)
!  STF_itow = STF_itog*STF_lamb
!  STF_srcradius = srcradius
!  STF_defocusspread = defocusspread
!  sc = cmplx(0.0,0.0)
!  anglethresh = STF_caperture*STF_RELANGTHRESH
!  threshwidth = STF_itow*STF_APSMOOTHPIX
!  PSfac = STF_pi*STF_srcradius/STF_lamb
!  PSfac = -0.5*PSfac*PSfac
!  PTfac = STF_pi*STF_defocusspread/STF_lamb
!  PTfac = -0.25*PTfac*PTfac
!  tpower = 0.0
!  camx = STF_caperture_movex*0.001
!  camy = STF_caperture_movey*0.001
!  lbtx = STF_beam_tiltx*0.001
!  lbty = STF_beam_tilty*0.001
!! ------------
!
!! ------------
!! preset Fourier-space array with exp[-i*chi]
!  tpower = 0.0
!  do j=1,ndim
!    wx = STF_TABBED_SCR(j)*STF_itow+lbtx
!    wy2 = wx*wx
!    wyap = (wx-camx)*(wx-camx)
!    do i=1,ndim
!      wy = STF_TABBED_SCR(i)*STF_itow+lbty
!      w2 = wy*wy+wy2
!      wap = wyap + (wy-camy)*(wy-camy)
!!      power = 1.0-sigmoid(sqrt(w2),anglethresh,threshwidth)
!!      if (power<0.001) cycle ! apply condenser aperture
!!      chi = STF_AberrationFunction(wx,wy)
!!      power = power * Exp(PSfac*w2) * Exp(PTfac*w2*w2)
!!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
!!      tpower = tpower + power*power
!      
!!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!!     use tabbed version here, (speed increased by approx 15%)
!      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
!      power = 1.0-rval
!      !power = 1.0-sigmoid(sqrt(w2),anglethresh,threshwidth)
!      
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      
!      chi = AF_AberrationFunction(wx,wy)  ! TRANSPOSED
!      power = power * Exp(PSfac*w2) * Exp(PTfac*w2*w2)
!!     use intrinsic version here, ( no speed gain with tabbed version)
!      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
!      sc(i,j) = cval*power
!!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
!      
!    end do
!  end do
!!  write(*,*) "FS-power:",tpower
!  tpowscal = 1.0/tpower/real(ndim*ndim)
!! ------------
!
!! ------------
!! transform to real space (scramble first!)
!  call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
!! ------------
!
!! ------------
!! calculate probe by taking absolute square in real space
!! + unscramble (shifts origin)
!  do j=1,ndim
!    j1 = STF_TABBED_USC(j)
!    do i=1,ndim
!      i1 = STF_TABBED_USC(i)
!      cval = sc(i,j)
!      p(i1,j1) = real( cval * conjg(cval) )*tpowscal
!    end do
!  end do
!! ------------
!
!! ------------
!!  write(*,*) " > STF_GetQCohProbe: EXIT."
!  return
!
!END SUBROUTINE STF_GetQCohProbe
!!**********************************************************************!
!
!










! -------> PARTIAL COHERENCE EXPLICITLY
! -------> PARTIAL COHERENCE EXPLICITLY
! -------> PARTIAL COHERENCE EXPLICITLY



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_SetupPSCKernel(rkernel,nkx2,nky2,sx,sy,rad,klim,ntype)
! function: prepares the convolution kernel in real space
! -------------------------------------------------------------------- !
! parameter: 
!       real*4 :: rkernel(-nkx2:nkx2,-nky2:nky2) ! kernel data
!       integer*4 :: nkx2, nky2 ! kernel half sizes
!       real*4 :: sx, sy        ! kernel sampling rate (nm)
!       real*4 :: rad           ! kernel radius parameter (nm)
!       real*4 :: klim          ! kernel radial limit (nm)
!       integer*4 :: ntype      ! kernel function selector
!                               ! 0 = Delta function = no convolution
!                               ! 1 = Gaussian
!                               ! 2 = Cauchy (Lorentzian)
!                               ! 3 = Disk (Sigmoidal)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1150
  real*4, intent(inout)             :: rkernel(-nkx2:nkx2,-nky2:nky2)
  integer*4, intent(in)             :: nkx2, nky2
  real*4, intent(in)                :: sx, sy, rad, klim
  integer*4, intent(in)             :: ntype
  
  integer*4 :: i, j, ki, kj, i2, j2, ntyuse
  real*4 :: kernprm, kernsum, kernval, rprm, aprm, kap, r2pi
  real*4 :: pkx, pky, pkys, pks

! ---------------
! INITIALIZATION
  rkernel = 0.0
  kernsum = 0.0
! limit the convolution type
  ntyuse = max(0,min(3,ntype))
  aprm = 1.3047660265
  rprm = aprm*rad
  kap = 1.5
  r2pi = 6.2831853072

! -------------------
! KERNEL PREPARATION
  select case (ntyuse)
  
    case (0) ! DELTA FUNCTION
      rkernel(0,0) = 1.0
      kernsum = 1.0
  
    case (1) ! GAUSSIAN
      !write(*,*) "- Gaussian kernel setup:"
      ! Exp[ -(x^2 + y^2) /HWHM^2 *Log[2]]
      kernprm = -LOG(2.0)/(rad*rad)
      do j = -nky2, nky2
        pky = real(j)*sy
        pkys = pky*pky
        do i = -nkx2, nkx2
          pkx = real(i)*sx
          pks = pkx*pkx + pkys
          if (pks>klim*klim) cycle ! cut-off
          kernval = EXP( pks*kernprm )
          kernsum = kernsum + kernval
          rkernel(i,j) = kernval
!          write(*,*) "  - kernel at (",pkx,",",pky,")nm = ",kernval
        end do
      end do
      ! normalization
      rkernel = rkernel / kernsum
      
    case (2) ! MODIFIED CAUCHY (LORENTZIAN)
      !  A * HWHM / ( 2*Pi *( ( A * HWHM )^2 + ( x - x0 )^2 + ( y - y0 )^2 )^[Kappa] )
      !  A = Sqrt[1/3 (1 + 2*2^(1/3) + 2^(2/3))] = 1.3047660265041066
      !  Kappa = 1.5
      !
      kernprm = rprm * rprm
      do j = -nky2, nky2
        pky = real(j)*sy
        pkys = pky*pky
        do i = -nkx2, nkx2
          pkx = real(i)*sx
          pks = pkx*pkx + pkys
          if (pks>klim*klim) cycle ! cut-off
          kernval = rprm / r2pi / ( kernprm + pks )**kap
          kernsum = kernsum + kernval
          rkernel(i,j) = kernval
        end do
      end do
      ! normalization
      rkernel = rkernel / kernsum
    
    case (3) ! DISC (SIGMOID)
      kernprm = rad
      do j = -nky2, nky2
        pky = real(j)*sy
        pkys = pky*pky
        do i = -nkx2, nkx2
          pkx = real(i)*sx
          pks = sqrt(pkx*pkx + pkys)
          if (pks>klim) cycle ! cut-off
          kernval = 0.5-0.5*tanh((pks/kernprm - 1.0)*100.0)
          kernsum = kernsum + kernval
          rkernel(i,j) = kernval
        end do
      end do
      ! normalization
      rkernel = rkernel / kernsum
      
  end select ! case (ntype)
! ------------
  
  return
  
END SUBROUTINE STF_SetupPSCKernel
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetPSCohProbe(p,ndim,sampling,srcradius,srctype)
! function: calculates the probefunction including partial spatial
!           coherence from the current coherent aberrations
!           and from the demagnified source radius
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!                                     probefunction result
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
!            real*4 :: srcradius : demagnified source radius [nm]
!            integer*4 :: srctype : source distribution type
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1100
  
  integer*4 :: ndim, srctype
  real*4 :: sampling, srcradius
  real*4, intent(out) :: p(ndim,ndim)
  
  integer*4 :: i, j, ndim2, ki, kj, i1, j1, ntyuse, nky2, nkx2, nalloc
  real*4 :: w2, wx, wy, chi, power, anglethresh, threshwidth, wy2, wap
  real*4 :: tpower, tpowscal, PSfac, rval, camx, camy, lbtx, lbty, klim
  real*4, allocatable :: rkernel(:,:)
  !real*4 :: pcoh(ndim,ndim)
  complex*8 :: cval
  
  !external :: ODDCC2048S
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetPSCohProbe: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetPSCohProbe: Module not initialized. Aborting.",subnum+1)
    return
  end if
! save data size
  STF_dim = ndim
  ndim2 = int(ndim/2)
  STF_sampling = sampling
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
! set Fourier-space sampling
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  STF_srcradius = srcradius
  sc = cmplx(0.0,0.0)
  sc2 = cmplx(0.0,0.0)
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  PSfac = STF_pi*STF_srcradius/STF_lamb
  PSfac = PSfac*PSfac
  tpower = 0.0
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
  ! limit the convolution type
  ntyuse = max(0,min(3,srctype))
  ! determine the size of the kernel to use
  select case (ntyuse)
    case (1)
      ! in case of a gaussian use 
      klim = 3.0*srcradius
    case (2)
      klim = 10.0*srcradius
    case (3)
      klim = 1.5*srcradius
  end select ! case (ntype)
  nkx2 = min( ndim2, CEILING(klim/sampling) )
  nky2 = nkx2
  ! klim = sampling*nkx2
  ! allocate memory for the working array and the kernel, zero both arrays
  allocate(rkernel(-nkx2:nkx2,-nky2:nky2), stat=nalloc)
  if (nalloc/=0) then
    call STF_ERROR("Failed to allocate memory for image convolution.", subnum+2)
    return
  end if
  rkernel = 0.0
! ------------

! ------------
! preset Fourier-space array with exp[-i*chi]
  tpower = 0.0
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow
    wy2 = (wx-camx)*(wx-camx)
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow
      wap = wy2 + (wy-camy)*(wy-camy)
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx+lbtx,wy+lbty)
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      sc(i,j) = cval*power
!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
    
    end do
  end do
!  write(*,*) "FS-power:",tpower
!  tpowscal = 1.0/real(ndim*ndim)/tpower
! ------------


! ------------
! transform coherent wavefunction to real space
  call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
! ------------

! ------------
! calculate probe intensity by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
      sc2(i1,j1) = cval * conjg(cval) ! * tpowscal
    end do
  end do
! ------------


! ------------
! Transform coherent probe to Fourier-space
  call STF_FFT(sc2,ndim,ndim,"for") ! external call from SFFTs.f
! ------------


! ------------
! Prepare the convolution kernel in real space
  call STF_SetupPSCKernel(rkernel, nkx2, nky2, sampling, sampling, &
     &                    srcradius, klim, ntyuse)
! Transfer real space kernel to working array for DFT
  sc = cmplx(0.0,0.0)
  do j=-nky2, nky2
    j1 = j + 1
    if (j1<1) j1 = j1 + ndim
    do i=-nkx2, nkx2
      i1 = i + 1
      if (i1<1) i1 = i1 + ndim
      sc(i1,j1) = rkernel(i,j)
    end do
  end do
! Transform kernel to fourier space
  call STF_FFT(sc,ndim,ndim,"for") ! external call from SFFTs.f
! ------------


! ------------
! multiply probe and kernel in Fourier-space for convolution
! scrambled and transposed
  do i=-ndim2, ndim2-1
    !wx = i*STF_itow+lbtx
    !wy2 = wx*wx
    ki = mod( i+ndim, ndim )+1
    do j=-ndim2, ndim2-1
      !wy = j*STF_itow+lbty
      kj = mod( j+ndim, ndim )+1
      !w2 = wy*wy+wy2
      sc2(kj,ki) = sc2(kj,ki) * conjg(sc(kj,ki))
    end do
  end do
! ------------

! ------------
! transform back to real space
  call STF_FFT(sc2,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate probe by taking absolute square in real space
  do j=1, ndim
    do i=1, ndim
      cval = sc2(i,j)
      rval = abs( real( cval ) )
      tpower = tpower + rval
      p(i,j) = rval
    end do
  end do
! rscale incident probe intensity to 1
  if (tpower>0.0) p = p / tpower
! ------------

! ------------
!  write(*,*) " > STF_GetPSCohProbe: EXIT."
  deallocate(rkernel,stat=nalloc)
  return

END SUBROUTINE STF_GetPSCohProbe
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetPTCohProbe(p,ndim,sampling,defocusspread)
! function: calculates the probefunction including partial temporal
!           coherence from the current coherent aberrations
!           and from the defocus spread
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!                                     probefunction result
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
!            real*4 :: srcradius : demagnified source radius [nm]
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !


  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1400
  
  integer*4 :: ndim
  real*4 :: sampling, defocusspread
  real*4, intent(out) :: p(ndim,ndim)
  
  integer*4 :: i, j, l, ndim2, i1, j1, NEF
  real*4 :: tpower, tpowscal, PTfac, PTscal, meandefocus, PTweight, wap
  real*4 :: defocus, dz, z0, ddz, rval, wx, wy, wy2, w2, chi, zfac
  real*4 :: anglethresh, threshwidth, power, wyap, camx, camy, lbtx, lbty
  real*4, dimension(:,:), allocatable :: pcoh, pdef, pcap
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetPTCohProbe: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetPTCohProbe: Module not initialized. Aborting.",subnum+1)
    return
  end if
! allocate pcoh data array
  allocate(pcoh(ndim,ndim),pdef(ndim,ndim),pcap(ndim,ndim),stat=i)
  if (i/=0) then
    call STF_ERROR("STF_GetPTCohProbe: Allocation failed. Aborting.",subnum+2)
    return
  end if
! save data size
  STF_dim = ndim
  ndim2 = int(ndim/2)
  STF_sampling = sampling
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  sc=cmplx(0.0,0.0)
!  sc2=cmplx(0.0,0.0)
! set Fourier-space sampling
  STF_defocusspread = defocusspread
  if (defocusspread==0.0) then
    PTfac = -1.0
    PTscal = 1.0
    NEF = 1
  else
    PTfac = -1.0/defocusspread/defocusspread
    PTscal = 1.0/(sqrt(STF_pi)*defocusspread)
    NEF = STF_DEFOCUS_KERNEL_STEPS
  end if
  
  
  zfac = STF_pi/STF_lamb
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  tpower = 0.0
  p = 0.0
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
  
  
! ------------




! ------------
! Get current mean defocus
  meandefocus = AF_wa(1,2)
! Get start defocus
  z0 = meandefocus - STF_DEFOCUS_KERNEL_SPREAD*defocusspread
! Get defocus kernel step
  dz = 2.0*STF_DEFOCUS_KERNEL_SPREAD*defocusspread/real(STF_DEFOCUS_KERNEL_STEPS-1)
! Set internal defocus register to zero
  AF_wa(1,2) = 0.0
! ------------

! ------------
! preset mean aberration function without defocus and defocus aberration field
! in temporary complex*8 datafield sc2
  pcoh = 0.0
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow
    wy2 = (wx+lbtx)**2
    wyap = (wx-camx)*(wx-camx)
    
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow
      w2 = (wy+lbty)**2 + wy2
      wap = (wy-camy)*(wy-camy) + wyap
      
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      
      pcap(i,j) = rval
      power = 1.0-rval
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      pdef(i,j) = w2
      pcoh(i,j) = AF_AberrationFunction(wx+lbtx,wy+lbty)  ! TRANSPOSED
    
    end do
  end do
!  write(*,*) "FS-power:",tpower
! ------------


! ------------
! sum up defocus panes
  ! loop over all defoci of kernel
  do l=1,NEF
  
    ! get current defocus
    defocus = z0 + real(l-1)*dz
    ddz = defocus - meandefocus
    
    ! calculate coherent probe function
    ! ------------
    ! preset Fourier-space array with exp[-i*chi]
    tpower = 0.0
    sc = 0.0
    do j=1,ndim
      do i=1,ndim
        rval = pcap(i,j)
        power = 1.0-rval
        if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
        chi = pcoh(i,j)+defocus*zfac*pdef(i,j)
        cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
        sc(i,j) = cval*power
!        tpower = tpower + power*power
      end do
    end do
    !  write(*,*) "FS-power:",tpower
    tpowscal = 1.0/real(ndim*ndim)/tpower
    
    ! ------------
    ! transform to real space
    call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
    
    ! ------------
    ! calculate probe by taking absolute square in real space
    ! + unscramble (shifts origin)
    ! + weight with distribution function and power
    PTweight = PTscal*EXP(PTfac*ddz*ddz) !*tpowscal  
    do j=1,ndim
      j1 = STF_TABBED_USC(j)
      do i=1,ndim
        i1 = STF_TABBED_USC(i)
        cval = sc(i,j)
        p(i1,j1) = p(i1,j1) + real( cval * conjg(cval) )*PTweight
      end do
    end do
    
  end do
! ------------

! ------------
! rescale incident probe intensity to 1
  tpower = sum(p)
  if (tpower>0.0) p = p / tpower
! ------------

! ------------
! reset mean defocus
  AF_wa(1,2) = meandefocus
  deallocate(pcoh,pcap,pdef,stat=i)
! ------------

! ------------
!  write(*,*) " > STF_GetPTCohProbe: EXIT."
  return

END SUBROUTINE STF_GetPTCohProbe
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetPCohProbe(p,ndim,sampling,srcradius,srctype,defocusspread)
! function: calculates the probefunction including partial temporal
!           and partial spatial coherence from the current coherent
!           aberrations, from the source size, and from the defocus spread
! -------------------------------------------------------------------- !
! parameter: real*4 :: p(ndim,ndim) : reference to array recieving the
!                                     probefunction result
!            integer*4 :: ndim : size of array
!            real*4 :: sampling : real-space sampling for probefunction
!            real*4 :: srcradius : demagnified source radius [nm]
!            integer*4 :: srctype : source distribution type
!            real*4 :: defocusspread : defocus spread [nm]
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !


  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1500
  
  integer*4 :: ndim
  real*4 :: sampling, defocusspread, srcradius
  real*4, intent(out) :: p(ndim,ndim)
  
  integer*4 :: i, j, l, ndim2, i1, j1, ki, kj, NEF, srctype
  integer*4 :: nkx2, nky2, ntyuse, nalloc
  real*4 :: tpower, tpowscal, PTfac, PTscal, meandefocus, PTweight
  real*4 :: defocus, dz, z0, ddz, rval, wx, wy, w2, chi, zfac, wy2
  real*4 :: anglethresh, threshwidth, power, PSFac, wyap, wap
  real*4 :: camx, camy, lbtx, lbty, klim
  real*4, dimension(:,:), allocatable :: ptcoh, pcap, rkernel
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetPTCohProbe: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetPCohProbe: Module not initialized. Aborting.",subnum+1)
    return
  end if
! allocate pcoh data array
  allocate(ptcoh(ndim,ndim),pcap(ndim,ndim),stat=i)
  if (i/=0) then
    call STF_ERROR("STF_GetPCohProbe: Allocation failed. Aborting.",subnum+2)
    return
  end if
! save data size
  STF_dim = ndim
  ndim2 = int(ndim/2)
  STF_sampling = sampling
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
  sc=cmplx(0.0,0.0)
  sc2=cmplx(0.0,0.0)
! set Fourier-space sampling
  STF_defocusspread = defocusspread
  if (defocusspread==0.0) then
    PTfac = -1.0
    PTscal = 1.0
    NEF = 1
  else
    PTfac = -1.0/defocusspread/defocusspread
    PTscal = 1.0/(sqrt(STF_pi)*defocusspread)
    NEF = STF_DEFOCUS_KERNEL_STEPS
  end if
  zfac = STF_pi/STF_lamb
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  STF_srcradius = srcradius
  STF_srctype = srctype
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  PSfac = STF_pi*STF_srcradius/STF_lamb
  PSfac = PSfac*PSfac
  tpower = 0.0
! preset outputarray  
  p = 0.0
  ptcoh = 0.0
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
  
  ! limit the convolution type
  ntyuse = max(0,min(3,srctype))
  ! determine the size of the kernel to use
  select case (ntyuse)
    case (1)
      ! in case of a gaussian use 
      klim = 3.0*srcradius
    case (2)
      klim = 10.0*srcradius
    case (3)
      klim = 1.5*srcradius
  end select ! case (ntype)
  nkx2 = min( ndim2, CEILING(klim/sampling) )
  nky2 = nkx2
  ! klim = sampling * nkx2
  ! allocate memory for the working array and the kernel, zero both arrays
  allocate(rkernel(-nkx2:nkx2,-nky2:nky2), stat=nalloc)
  if (nalloc/=0) then
    call STF_ERROR("Failed to allocate memory for image convolution.", subnum+3)
    return
  end if
  rkernel = 0.0
! ------------




! ------------
! Get current mean defocus
  meandefocus = AF_wa(1,2)
! Get start defocus
  z0 = meandefocus - STF_DEFOCUS_KERNEL_SPREAD*defocusspread
! Get defocus kernel step
  dz = 2.0*STF_DEFOCUS_KERNEL_SPREAD*defocusspread/real(STF_DEFOCUS_KERNEL_STEPS-1)
! Set internal defocus register to zero
  AF_wa(1,2) = 0.0
! ------------

! ------------
! preset mean aberration function without defocus and defocus aberration field
! in temporary complex*8 datafield sc2
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow+lbtx
    wy2 = (STF_TABBED_SCR(j)*STF_itow)**2
    wyap = (wx-camx)**2
    
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow+lbty
      w2 = (STF_TABBED_SCR(i)*STF_itow)**2 + wy2
      wap = wyap + (wy-camy)**2
      
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      pcap(i,j) = rval
      power = 1.0-rval
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx,wy)  ! TRANSPOSED
      
      sc2(i,j) = cmplx(chi,zfac*w2);
    
    end do
  end do
!  write(*,*) "FS-power:",tpower
! ------------


! ------------
! sum up defocus panes
  ! loop over all defoci of kernel
  do l=1,NEF
  
    ! get current defocus
    defocus = z0 + real(l-1)*dz
    ddz = defocus - meandefocus
    
    ! calculate coherent probe function
    ! ------------
    ! preset Fourier-space array with exp[-i*chi]
    sc=0.0
    tpower = 0.0
    do j=1,ndim
        do i=1,ndim
            rval = pcap(i,j)
            power = 1.0-rval
            if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
            chi = real(sc2(i,j))+defocus*imag(sc2(i,j))
            cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
            sc(i,j) = cval*power
!            tpower = tpower + power*power
        end do
    end do
    !  write(*,*) "FS-power:",tpower
!    tpowscal = 1.0/real(ndim*ndim)/tpower
    
    ! ------------
    ! transform to real space
    call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
    
    ! ------------
    ! calculate probe by taking absolute square in real space
    ! + unscramble (shifts origin)
    ! + weight with distribution function and power
    !! PTfac = -1.0/defocusspread/defocusspread
    !! PTscal = 1.0/(sqrt(STF_pi)*defocusspread)
    PTweight = PTscal*EXP(PTfac*ddz*ddz) !*tpowscal  
    do j=1,ndim
        j1 = STF_TABBED_USC(j)
        do i=1,ndim
            i1 = STF_TABBED_USC(i)
            cval = sc(i,j)
            ptcoh(i1,j1) = ptcoh(i1,j1) + real( cval * conjg(cval) )*PTweight
        end do
    end do
    
  end do
! ------------

! ------------
! reset mean defocus
  AF_wa(1,2) = meandefocus
! ------------

! NOW APPLY SPATIAL INCOHERENCE ON PARTIALLY TEMPORAL INCOHERENT PROBE
! ------------
! sc2 can be used again
!  sc2=cmplx(0.0,0.0)
  do j=1,ndim
    do i=1,ndim
      sc2(i,j) = ptcoh(i,j)
    end do
  end do
! Transform current probe to Fourier-space
  call STF_FFT(sc2,ndim,ndim,"for") ! external call from SFFTs.f
! ------------

! ------------
! Prepare the convolution kernel in real space
  call STF_SetupPSCKernel(rkernel, nkx2, nky2, sampling, sampling, &
     &                    srcradius, klim, ntyuse)
! Transfer real space kernel to working array for DFT
  sc = cmplx(0.0,0.0)
  do j=-nky2, nky2
    j1 = j + 1
    if (j1<1) j1 = j1 + ndim
    do i=-nkx2, nkx2
      i1 = i + 1
      if (i1<1) i1 = i1 + ndim
      sc(i1,j1) = rkernel(i,j)
    end do
  end do
! Transform kernel to fourier space
  call STF_FFT(sc,ndim,ndim,"for") ! external call from SFFTs.f
! ------------

! ------------
! multiply probe and kernel in Fourier-space for convolution
! scrambled and transposed
  do i=-ndim2, ndim2-1
    !wx = i*STF_itow+lbtx
    !wy2 = wx*wx
    ki = mod( i+ndim, ndim )+1
    do j=-ndim2, ndim2-1
      !wy = j*STF_itow+lbty
      kj = mod( j+ndim, ndim )+1
      !w2 = wy*wy+wy2
      sc2(kj,ki) = sc2(kj,ki) * conjg(sc(kj,ki))
    end do
  end do
! ------------

! ------------
! transform back to real space
  call STF_FFT(sc2,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate probe in real space
  tpower = 0.0
  do j=1, ndim
    do i=1, ndim
      cval = sc2(i,j)
      rval = abs(real(cval))
      p(i,j) = rval
      tpower = tpower + rval
    end do
  end do
  if (tpower>0.0) p = p / tpower ! rescale to probe intensity 1
! ------------


! ------------
!  write(*,*) " > STF_GetPTCohProbe: EXIT."
  deallocate(ptcoh,pcap,rkernel,stat=i)
  return

END SUBROUTINE STF_GetPCohProbe
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GetRonchigram(p, ndim, sampling, defocusspread)
! function: calculates a ronchigram in p, using current aberrations
!           a random object (amorph), slightly randomly changed,
!           and partial temporal coherence
! -------------------------------------------------------------------- !
! parameter: p(ndim,ndim) : real*4 : ronchigram
!            ndim : integer*4 : data size
!            sampling : real*4 : real-space sampling
!            defocusspread : real*4 : defocus spread
! -------------------------------------------------------------------- !
! required link: SFFTs.f
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3300
  
  integer*4, intent(in) :: ndim
  real*4, intent(in) :: sampling, defocusspread
  real*4, dimension(ndim,ndim), intent(inout) :: p
  
  integer*4 :: i, j, l, ndim2, i1, j1, NEF, k, nk, nerr, navar
  real*4 :: tpower, tpowscal, PTfac, PTscal, meandefocus, PTweight
  real*4 :: defocus, dz, z0, ddz, rval, wx, wy, w2, chi, zfac, wy2, tsc, itsc
  real*4 :: anglethresh, threshwidth, power, wyap, wap, camx, camy, lbtx, lbty
  real*4, allocatable :: pcap(:,:)
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetRonchigram: INIT."
  nerr = 0
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
    call STF_ERROR("STF_GetRonchigram: Module not initialized. Aborting.",subnum+1)
    return
  end if
! allocate pcap data array
  allocate(pcap(ndim,ndim),stat=i)
  if (i/=0) then
    call STF_ERROR("STF_GetRonchigram: Allocation failed. Aborting.",subnum+2)
    return
  end if
! save data size
  STF_dim = ndim
  ndim2 = int(ndim/2)
  tsc = 1/real(ndim)
  itsc = 1.0/tsc
  nk = STF_AMORPH_NSLC
  STF_sampling = sampling
  call STF_SETTAB_SCR(ndim)
  call STF_SETTAB_USC(ndim)
!  sc=cmplx(0.0,0.0)
  sc2=cmplx(0.0,0.0)
! set Fourier-space sampling
  STF_defocusspread = defocusspread
  if (defocusspread==0.0) then
    PTfac = -1.0
    PTscal = 1.0
    NEF = 1
  else
    PTfac = -1.0/defocusspread/defocusspread
    PTscal = 1.0/(sqrt(STF_pi)*defocusspread)
    NEF = STF_DEFOCUS_KERNEL_STEPS
  end if
  zfac = STF_pi/STF_lamb
  STF_itog = 0.5/sampling/real(ndim2)
  STF_itow = STF_itog*STF_lamb
  tpower = 0.0
! preset outputarray  
  p = 0.0
!  ptcoh = 0.0
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
! prepare real-space object function (amorph)
  call STF_SETAMORPH(nerr)
  if (nerr/=0) return
  navar = STF_GETAMORPHVARIANT()
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------


! ------------
! Get current mean defocus
  meandefocus = AF_wa(1,2)
! Get start defocus
  z0 = meandefocus - STF_DEFOCUS_KERNEL_SPREAD*defocusspread
! Get defocus kernel step
  dz = 2.0*STF_DEFOCUS_KERNEL_SPREAD*defocusspread/real(STF_DEFOCUS_KERNEL_STEPS-1)
! Set internal defocus register to zero
  AF_wa(1,2) = 0.0
! ------------

! ------------
! preset mean aberration function without defocus and defocus aberration field
! in temporary complex*8 datafield sc2
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow+lbtx
    wy2 = wx**2
    wyap = (wx-camx)*(wx-camx)
    
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow+lbty
      w2 = wy**2 + wy2
      wap = (wy-camy)*(wy-camy) + wyap
      
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      pcap(i,j) = rval
      power = 1.0-rval
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx,wy)  ! TRANSPOSED
      sc2(i,j) = cmplx(chi,zfac*w2) ! save mean phase in real part and focus phase basis in imaginary part
!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
    end do
  end do
!  write(*,*) "FS-power:",tpower
! ------------


! ------------
! sum up defocus panes
  ! loop over all defoci of kernel
  do l=1,NEF
  
    ! get current defocus
    defocus = z0 + real(l-1)*dz
    ddz = defocus - meandefocus
    
    ! calculate coherent probe function
    ! ------------
    ! preset Fourier-space array with exp[-i*chi]
    sc=0.0
    tpower = 0.0
    do j=1,ndim
      do i=1,ndim
        rval = pcap(i,j)
        power = 1.0-rval
        if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
            
        chi = real(sc2(i,j))+defocus*imag(sc2(i,j))
            
        cval = exp(cmplx(0.0,-chi)) !cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
        sc(i,j) = cval*power
        tpower = tpower + power*power
        
      end do
    end do
    !  write(*,*) "FS-power:",tpower
    tpowscal = 1.0/tpower
    
    ! ------------
    ! transform to real space
    call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
    
    ! ------------
    ! do the elastic scattering (multislice)
    if (STF_AMORPH_NSLC>1)  then
      do k=1,nk
        do j=1,ndim
          do i=1,ndim
            sc(i,j) = sc(i,j)*STF_amorph_phgr(i,j,navar,k)
          end do
        end do
        call STF_FFT(sc,ndim,ndim,"for") ! external call from SFFTs.f
        do j=1,ndim
          do i=1,ndim
            sc(i,j) = sc(i,j)*STF_amorph_prop(i,j)
          end do
        end do
        if (k==STF_AMORPH_NSLC) exit
        call STF_FFT(sc,ndim,ndim,"back") ! external call from SFFTs.f
        !sc = sc * tsc
        navar = STF_GETAMORPHVARIANT() ! update variant number for the next slice
      end do
    else ! (STF_AMORPH_NSLC<=1)
      do j=1,ndim
        do i=1,ndim
          sc(i,j) = sc(i,j)*STF_amorph_phgr(i,j,navar,1)
        end do
      end do
      ! transform to fourier space
      call STF_FFT(sc,ndim,ndim,"for") ! external call from SFFTs.f
    end if ! (STF_AMORPH_NSLC>1) 
    
    
    
    ! ------------
    ! calculate ronchigram by taking absolute square in fourier space
    ! + unscramble (shifts origin)
    ! + weight with distribution function and power
    !! PTfac = -1.0/defocusspread/defocusspread
    !! PTscal = 1.0/(sqrt(STF_pi)*defocusspread)
    PTweight = PTscal*EXP(PTfac*ddz*ddz)*tpowscal
    do j=1,ndim
      j1 = STF_TABBED_USC(j)
      do i=1,ndim
        i1 = STF_TABBED_USC(i)
        cval = sc(i,j)
        p(j1,i1) = p(j1,i1) + real( cval * conjg(cval) )*PTweight
      end do
    end do
    
  end do
  
!  tpower = 0.0
!  do j=1,ndim
!    do i=1,ndim
!      tpower = tpower + p(j,i)
!    end do
!  end do
! ------------

! ------------
! reset mean defocus
  AF_wa(1,2) = meandefocus
! ------------

! ------------
!  write(*,*) " > STF_GetRonchigram: EXIT."
  return

END SUBROUTINE STF_GetRonchigram
!**********************************************************************!





































! -----------> FULL COHERENT THINGS
! -----------> FULL COHERENT THINGS
! -----------> FULL COHERENT THINGS



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
!  write(*,*) " > STF_GetCohProbeWave: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
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
  call STF_FFT(sc,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate wave components by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
      re(i1,j1) = real(cval)
      im(i1,j1) = imag(cval)
!      re(i,j) = real(cval)
!      im(i,j) = imag(cval)
    end do
  end do
! ------------

! ------------
!  write(*,*) " > STF_GetCohProbeWave: EXIT."
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
  real*4 :: wx, wy, chi, power, anglethresh, threshwidth
  real*4 :: tpower, tpowscal, rval, wyap, wap, camx, camy, lbtx, lbty
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetCohProbeWave: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
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
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------

! ------------
! preset Fourier-space array with exp[-i*chi]
!!!  call STF_GetProbeWaveFourier(sc,ndim,sampling)
! ------------
! preset Fourier-space array with exp[-i*chi]
  sc(:,:) = cmplx(0.0,0.0)
  tpower = 0.0
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow+lbtx
    wyap = (wx-camx)*(wx-camx)
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow+lbty
      wap = (wy-camy)*(wy-camy)+wyap
!      power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      chi = STF_AberrationFunction(wx,wy)
!      wave(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx,wy)
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      sc(i,j) = cval*power
!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
      tpower = tpower + power*power
      
    end do
  end do
!  write(*,*) "FS-power:",tpower
  tpowscal = 1.0/tpower !/real(ndim*ndim)
! ------------
! ------------

! ------------
! transform to real space (scramble first!)
  call STF_FFT(sc,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate wave components by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
      re(i1,j1) = real(cval)*tpowscal
!      im(j1,i1) = imag(cval)
!      re(i,j) = real(cval)
!      im(i,j) = imag(cval)
    end do
  end do
! ------------

! ------------
!  write(*,*) " > STF_GetCohProbeWaveRe: EXIT."
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
  real*4 :: wx, wy, chi, power, anglethresh, threshwidth
  real*4 :: tpower, tpowscal, rval, wyap, wap, camx, camy, lbtx, lbty
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetCohProbeWaveIm: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
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
  anglethresh = STF_caperture*STF_RELANGTHRESH
  threshwidth = STF_itow*STF_APSMOOTHPIX
  camx = STF_caperture_movex*0.001
  camy = STF_caperture_movey*0.001
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------

! ------------
! preset Fourier-space array with exp[-i*chi]
!!!  call STF_GetProbeWaveFourier(sc,ndim,sampling)
! ------------
! preset Fourier-space array with exp[-i*chi]
  tpower = 0.0
  sc(:,:) = cmplx(0.0,0.0)
  do j=1,ndim
    wx = STF_TABBED_SCR(j)*STF_itow+lbtx
    wyap = (wx-camx)*(wx-camx)
    do i=1,ndim
      wy = STF_TABBED_SCR(i)*STF_itow+lbty
      wap = (wy-camy)*(wy-camy) + wyap
!      power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
!      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
!      chi = STF_AberrationFunction(wx,wy)
!      wave(i,j) = cmplx(cos(chi),-sin(chi))*power
!      tpower = tpower + power*power
      
!      rval = sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)  ! intrinsict calculation method
!     use tabbed version here, (speed increased by approx 15%)
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)  ! tabbed calculation method
      power = 1.0-rval
      !power = 1.0-sigmoid(sqrt(wx*wx+wy*wy),anglethresh,threshwidth)
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx,wy)
!     use intrinsic version here, ( no speed gain with tabbed version)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
!      call STF_TABBED_EXP(-chi,cval) ! tabbed calculation method
      sc(i,j) = cval*power
!      sc(i,j) = cmplx(cos(chi),-sin(chi))*power
      tpower = tpower + power*power
      
    end do
  end do
!  write(*,*) "FS-power:",tpower
  tpowscal = 1.0/tpower !/real(ndim*ndim)
! ------------
! ------------

! ------------
! transform to real space (scramble first!)
  call STF_FFT(sc,ndim,ndim,"bac") ! external call from SFFTs.f
! ------------

! ------------
! calculate wave components by taking absolute square in real space
! + unscramble (shifts origin)
  do j=1,ndim
    j1 = STF_TABBED_USC(j)
    do i=1,ndim
      i1 = STF_TABBED_USC(i)
      cval = sc(i,j)
!      re(j1,i1) = real(cval)*tpowscal
      im(i1,j1) = imag(cval)*tpowscal
!      re(i,j) = real(cval)
!      im(i,j) = imag(cval)
    end do
  end do
! ------------

! ------------
!  write(*,*) " > STF_GetCohProbeWaveIm: EXIT."
  return

END SUBROUTINE STF_GetCohProbeWaveIm
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_PrepareProbeWaveFourier(ndimx,ndimy,samplingx,samplingy)
! function: calculates the prove wavefunction in fourier space
!           from aberration and saves it in global public array sc
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
  real*4 :: wx, wy, wx2, chi, power, anglethresh, threshwidth, wap
  real*4 :: tpower, rval, itogx, itogy, itowx, itowy
  real*4 :: camx, camy, lbtx, lbty
  complex*8 :: cval
  
  real*4, external :: sigmoid
! ------------

! ------------
! INIT
!  write(*,*) " > STF_PrepareProbeWaveFourier: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
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
  sc(:,:) = 0.0
  tpower = 0.0
  do j=1,ndimx
    wx = STF_TABBED_SCR(j)*itowx+lbtx
    wx2 = (wx-camx)*(wx-camx)
    do i=1,ndimy
      wy = STF_TABBED_SCR2(i)*itowy+lbty
      wap = wx2 + (wy-camy)*(wy-camy)
      call STF_TABBED_SIGMOID(sqrt(wap),anglethresh,threshwidth,rval)
      power = 1.0-rval
      
      if (power<STF_APERTURETHRESH) cycle ! apply condenser aperture
      
      chi = AF_AberrationFunction(wx,wy)
      cval = cmplx(cos(chi),-sin(chi)) ! intrinsict calculation method
      sc(i,j) = cval*power
      STF_PreparedWavePower = STF_PreparedWavePower + power*power
      
    end do
  end do
!  write(*,*) "FS-power:",tpower
! ------------

! ------------
!  write(*,*) " > STF_PrepareProbeWaveFourier: EXIT."
  return

END SUBROUTINE STF_PrepareProbeWaveFourier
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
  real*4 :: wx, wy, lbtx, lbty
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GetPhasePlate: INIT."
! check init-status
  if (STF_status<=0.or.AF_maxaberration<=0) then
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
  lbtx = STF_beam_tiltx*0.001
  lbty = STF_beam_tilty*0.001
! ------------


! ------------
! preset Fourier-space array with exp[-i*chi]
  do j=1,ndim
    wy = (j-ndim2-1)*STF_itow + lbty
    do i=1,ndim
      wx = (i-ndim2-1)*STF_itow + lbtx
      
      p(i,j) = AF_AberrationFunction(wx,wy)
      
    end do
  end do
! ------------

! ------------
!  write(*,*) " > STF_GetPhasePlate: EXIT."
  return

END SUBROUTINE STF_GetPhasePlate
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE STF_GenerateProbe(nProbeType, rData, ndim, sampling)
! function: General Interface Function for Probe Generation
!           Be aware to set all other parameters befor calling
!           STF_lamb, STF_caperture,  STF_srcradius, STF_defocusspread
!           STF_DEFOCUS_KERNEL_STEPS, STF_ STF_DEFOCUS_KERNEL__SPREAD
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nProbeType : switch for different function
!            real*4 :: rData(ndim,ndim) : reference to data array
!            integer*4 :: ndim : size of data array
!            real*4 :: sampling : real space sampling [nm/pix]
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2200
  
  integer*4, intent(in) :: nProbeType, ndim
  real*4, dimension(ndim,ndim), intent(inout) :: rData
  real*4, intent(inout) :: sampling
  
  integer*4 :: ierr
! ------------

! ------------
! INIT
!  write(*,*) " > STF_GenerateProbe: INIT."
  ierr = 0
  if (STF_status<1.or.AF_maxaberration<1) then
    call STF_ERROR("AF or STF module not ready.",subnum+1)
    return
  end if
  if (nProbeType<1.or.nProbeType>STF_FUNC_NUM) then
    call STF_ERROR("STF_GenerateProbe: Invalid parameter (nProbeType).",subnum+2)
    return
  end if
  if (ndim<16.or.ndim>FFT_BOUND_MAX) then
    call STF_ERROR("STF_GenerateProbe: Invalid parameter (ndim).",subnum+3)
    return
  end if
  if (sampling<=0.0) then
    call STF_ERROR("STF_GenerateProbe: Invalid parameter (sampling).",subnum+4)
    return
  end if
  STF_dim = ndim
  STF_sampling = sampling
! ------------

! ------------
! select the type
  select case ( nProbeType )
    case (1)
        ! "Coherent Probe Intensity"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetCoherentProbe(rData, ndim, STF_sampling)
    case (2)
        ! "Coherent Probe Phase"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetCoherentProbePhase(rData, ndim, STF_sampling)
    case (3)
        ! "Partially Spatial Coherent Probe Intensity"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetPSCohProbe(rData, ndim, STF_sampling, STF_srcradius, STF_srctype)
    case (4)
        ! "Partially Temporal Coherent Probe Intensity"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetPTCohProbe(rData, ndim, STF_sampling, STF_defocusspread)
    case (5)
        ! "Partially Coherent Probe Intensity"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetPCohProbe(rData, ndim, STF_sampling, STF_srcradius, STF_srctype, STF_defocusspread)
    case (6)
        ! "Real Part of Probe Wave Function"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetCohProbeWaveRe(rData, ndim, STF_sampling)
    case (7)
        ! "Imaginary Part of Probe Wave Function"
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetCohProbeWaveIm(rData, ndim, STF_sampling)
    case (8)
        ! "Phase plate"
        call STF_GetPhasePlate(rData, ndim, STF_sampling)
    case (9)
        ! "Ronchigram, weak scattering"
        STF_amorph_relstr = 1.0
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetRonchigram(rData, ndim, STF_sampling, STF_defocusspread)
    case (10)
        ! "Ronchigram, strong scattering"
        STF_amorph_relstr = 1.67
        call STF_INIT_FFT(ndim, ierr )
        if (ierr==0) call STF_GetRonchigram(rData, ndim, STF_sampling, STF_defocusspread)
  end select
! ------------

! ------------
!  write(*,*) " > STF_GenerateProbe: EXIT."
  return

END SUBROUTINE STF_GenerateProbe
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
!  integer*4, parameter :: subnum = 4100
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