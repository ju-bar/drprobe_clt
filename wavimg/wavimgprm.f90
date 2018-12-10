!**********************************************************************!
!**********************************************************************!
!
! FILE: "wavimgprm.f90"
!
! AUTHOR: Dr. J. Barthel
!         Ernst Ruska-Centre
!         RWTH Aachen University, Aachen, Germany
!         Forschungszentrum, Jülich, Germany
!              
!         ju.barthel@fz-juelich.de
!
! PURPOSE: Global parameter declaration for program WAVIMG
!          and parameter file handling
!
! USAGE:  use wavimgprm
!         call wavimgprm_init(nerr) ! to allocate module arrays
!         ...
!         call wavimgprm_uninit(nerr) ! to deallocate module arrays
!
! VERSION: 1.1, J.B., 13.12.2017
!
!**********************************************************************!
!**********************************************************************!


MODULE wavimgprm

  implicit none
  
! --------- [ module routines ] ----------------------------------------------------------
  
  public :: wavimgprm_init
  public :: wavimgprm_uninit
  public :: wavimgprm_init2
  public :: wavimgprm_initcl
  public :: wavimgprm_getcls

! --------- [ global parameters ] --------------------------------------------------------
  integer*4, parameter :: stdout = 6        ! standard output unit
  integer*4, parameter :: namax = 50        ! max. number of aberration definitions       
!  integer*4, parameter :: nft = 2048        ! max. size for ft
  integer*4, parameter :: nloopmax = 5      ! max. number of parameter loops
  integer*4, parameter :: nalmax = 2048     ! max. number of list-directed aberration variations
  integer*4, parameter :: ntfmax = 10240    ! max. number of tf values
  real*4, parameter :: d2r = 0.0174533      ! degree to radian factor
  real*4, parameter :: twopi = 6.28318531   ! 2*Pi
  real*4, parameter :: pi = 3.14159265      ! Pi
  integer*4, parameter :: nft_min = 128     ! min. fft plan size
  integer*4, parameter :: nft_max = 8192    ! max. fft plan size
  
! --------- [ global variables ] ---------------------------------------------------------
  integer*4 :: nft                          ! max. size for ft
  integer*4 :: nerr                         ! error code
  integer*4 :: ndbg                         ! debug export flag
  integer*4 :: nsil                         ! message flag
  integer*4 :: numw                         ! number of warnings
  integer*4 :: nume                         ! number of errors
  integer*4 :: nsre                         ! flag for supressing fokal averaging resonances
  integer*4 :: nnli                         ! flag for suppressing loop image output in map output mode
  real*4 :: wl                              ! wavelength [nm]
  real*4 :: gres1                           ! resonance frequency of focal kernel [1/nm]
  real*4 :: ginfo                           ! information limit [1/nm]
  complex*8, allocatable :: cft(:,:)        ! transformation array
  complex*8, allocatable :: cft2(:,:)       ! transformation array 2
  character(len=600) :: sdbgmsg             ! debug message
  
! --------- [ input variables ] ----------------------------------------------------------
  character(len=600) :: sprmfile            ! input parameter file name
  character(len=600) :: swavfile            ! wave function file name
  character(len=600) :: simgfile            ! image file name
  character(len=600) :: simgfile_ex         ! image file name, external parameter
  character(len=600) :: smtffile            ! detector mtf file name
  real*4 :: ht                              ! TEM high tension
  integer*4 :: nwx, nwy                     ! wave discretization
  real*4 :: swx, swy                        ! wave sampling
  integer*4 :: notype                       ! output type option
  integer*4 :: nx, ny                       ! image discretization
  real*4 :: ox, oy, ax, simg                ! image frame offset, rotation and sampling
  integer*4 :: nfoc_ex                      ! external focus usage flag
  integer*4 :: ncohm                        ! index of coherence model
  integer*4 :: doptc                        ! flags usage of partial temporal coherence
  integer*4 :: dopsc                        ! flags usage of partial spatial coherence
  integer*4 :: domtf                        ! flags usage of detector mtf
  integer*4 :: dofrm                        ! flags usage of image frame extraction
  integer*4 :: dovib                        ! flags usage of vibration envelope (0 = OFF, 1 = ISO, 2 = ANISO)
  integer*4 :: doint                        ! flags creation of integer image
  integer*4 :: dornsb                       ! flags the renormalization of output waves created with side-band shifted MTF 
  real*4 :: foc_ex                          ! focus value [nm] specified via command-line argument
  real*4 :: btx, bty                        ! beam tilt for the transfer terms [mrad]
  real*4 :: oapr, oapx, oapy                ! objective aperture parameters radius, decnter x and decenter y [mrad]
  real*4 :: sbshx, sbshy                    ! side-band shift with respect to the zero beam [pix]
  real*4 :: oapr_ex                         ! objective aperture radius external (command-line)
  real*4 :: fs                              ! focus spread parameter [nm]
  integer*4 :: NKFS                         ! number of focal kernel samples offside central focus
  real*4 :: fkw                             ! relative width of the focal kernel wrt. fs
  real*4 :: sc                              ! semi-angle of convergence [mrad]
  integer*4 :: NKCB                         ! number of convergent beam kernel samples offside central beam
  real*4 :: cbkw                            ! relative width of the convergent beam kernel wrt. sc
  integer*4 :: nal                          ! number of aberrations loaded
  integer*4 :: nmtf                         ! number of mtf file pixels
  real*4, allocatable :: mtfdata(:)         ! mtf data
  real*4, allocatable :: ntfdata(:)         ! noise-tf data
  real*4 :: mtfscal                         ! mtf scaling
  real*4 :: vibamp, vibamp2, vibdir         ! vibration amplitudes and orientation
  real*4 :: iimean, dark_noise, el_conv     ! integer image: vacuum intensity mean, readout noise, counts/electron converision

! --------- [ variable loop definitions ] ------------------------------------------------
  integer*4 :: nloop                            ! number of loops (input)
  DATA nloop /0/
  integer*4, allocatable :: lpidx(:)            ! loop indices (input)
  integer*4, allocatable :: lpsz(:)             ! loop sizes (input)
  integer*4, allocatable :: lpcl(:)             ! loop class (input)
  integer*4, allocatable :: lpvr(:)             ! loop variable (input)
  integer*4, allocatable :: lpvf(:)             ! loop variation form (input)
  real*4, allocatable    :: lpv0(:)             ! loop variable start (input)
  real*4, allocatable    :: lpv1(:)             ! loop variable stop (input)
  real*4, allocatable    :: lpvd(:)             ! loop variable delta (input)
  character(len=2048), allocatable :: lpid(:)   ! loop ID string (input) (0 terminated)
  integer*4, allocatable :: lpali(:,:)          ! loop aberration parameter list index (header)
  real*4, allocatable    :: lpalv(:,:,:)        ! loop aberration parameter list values
  integer*4 :: nloopsteps                       ! total number of loop steps to be performed
  DATA nloopsteps /1/
  integer*4 :: iloopstep                        ! current loop step
  DATA iloopstep /0/
  
      
! ---------- [ temporary and other variables ] -------------------------------------------
  character(len=2048) :: stmp, stmpx , stmpy ! temp. string
  character(len=2048) :: siout, sioutbk      ! modified image output file name
  character(len=2048) :: siwav               ! modified input wave file name
  

! ---------- [ image caclulation variables ] ---------------------------------------------
  integer*4, allocatable :: ic_iwx(:), ic_iwy(:)        ! wave frequency hash
  real*4, allocatable :: ic_twx(:), ic_twy(:)           ! wave transfer frequencies (+beam tilt)
  real*4, allocatable :: ic_objaper(:,:)                ! objective aperture (on wave beams)
  
! ---------- [ run timing flag and variables ] -------------------------------------------
  integer*4, public :: wicl_runtimes  ! flags run time info  0: OFF (default) 1: ON (/rti)
  DATA wicl_runtimes /0/
  integer*8, public :: wicl_init      ! stores initial clock count 
  DATA wicl_init /0/
  integer*8, public :: wicl_rate      ! stores clock count rate rate (cts/s)
  DATA wicl_rate /1/
  integer*8, public :: wicl_max       ! stores clock max. count
  DATA wicl_max /1/
  integer*8, public :: wicl_measure   ! stores recent measured count
  DATA wicl_measure /0/
  
  



CONTAINS



subroutine wavimgprm_init(ierr)

  implicit none
  
  integer*4, intent(inout) :: ierr
  integer*4 :: nalloc, m, n
  
  nalloc = 0
  ierr = 0
  m = 1
  n = 1
  
! MOVED TO wavimgprm_init2
! JB 2017-12-13 for version 0.66
!  allocate(cft(nft,nft), cft2(nft,nft), stat=nalloc)
!  if (nalloc/=0) then
!    ierr = 1
!    return
!  end if
!  cft = cmplx(0.0,0.0)
!  cft2 = cmplx(0.0,0.0)
!  
!  allocate(mtfdata(ntfmax), ntfdata(ntfmax), stat=nalloc)
!  if (nalloc/=0) then
!    ierr = 2
!    return
!  end if
!  mtfdata = 0.0
!  ntfdata = 0.0
  
  n = nloopmax
  allocate(lpidx(n), lpsz(n), lpcl(n), lpvr(n), lpvf(n) , stat=nalloc)
  if (nalloc/=0) then
    ierr = 3
    return
  end if
  nloop = 0
  lpidx = 0
  lpidx = 0
  lpsz = 0
  lpcl = 0
  lpvr = 0
  lpvf = 0
  
  n = nloopmax
  allocate(lpv0(n), lpv1(n), lpvd(n), stat=nalloc)
  if (nalloc/=0) then
    ierr = 4
    return
  end if
  lpv0 = 0.0
  lpv1 = 0.0
  lpvd = 0.0
  
  n = nloopmax
  allocate(lpid(n), stat=nalloc)
  if (nalloc/=0) then
    ierr = 5
    return
  end if
  lpid = ""
  
  m = namax
  n = nloopmax
  allocate(lpali(0:m,0:n), stat=nalloc)
  if (nalloc/=0) then
    ierr = 6
    return
  end if
  lpali = 0
  
  m = namax
  n = nloopmax
  allocate(lpalv(1:m,1:nalmax,1:n), stat=nalloc)
  if (nalloc/=0) then
    ierr = 7
    return
  end if
  lpalv = 0.0
  
  return

end subroutine wavimgprm_init




subroutine wavimgprm_uninit(ierr)

  implicit none
  
  integer*4, intent(inout) :: ierr
  integer*4 :: nalloc
  
  nalloc = 0
  ierr = 0
  
  if (allocated(cft)) deallocate(cft, stat=nalloc)
  if (allocated(cft2)) deallocate(cft2, stat=nalloc)
  if (allocated(mtfdata)) deallocate(mtfdata, stat=nalloc)
  if (allocated(ntfdata)) deallocate(ntfdata, stat=nalloc)
  if (allocated(lpidx)) deallocate(lpidx, stat=nalloc)
  if (allocated(lpsz)) deallocate(lpsz, stat=nalloc)
  if (allocated(lpcl)) deallocate(lpcl, stat=nalloc)
  if (allocated(lpvr)) deallocate(lpvr, stat=nalloc)
  if (allocated(lpvf)) deallocate(lpvf, stat=nalloc)
  if (allocated(lpv0)) deallocate(lpv0, stat=nalloc)
  if (allocated(lpv1)) deallocate(lpv1, stat=nalloc)
  if (allocated(lpvd)) deallocate(lpvd, stat=nalloc)
  if (allocated(lpid)) deallocate(lpid, stat=nalloc)
  if (allocated(lpali)) deallocate(lpali, stat=nalloc)
  if (allocated(lpalv)) deallocate(lpalv, stat=nalloc)
  
  return

end subroutine wavimgprm_uninit



subroutine wavimgprm_init2(ierr)

  implicit none
  
  integer*4, intent(inout) :: ierr
  integer*4 :: nalloc
  
  ierr = 0
  nalloc = 0
  
  allocate(cft(nft,nft), cft2(nft,nft), stat=nalloc)
  if (nalloc/=0) then
    ierr = 1
    return
  end if
  cft = cmplx(0.0,0.0)
  cft2 = cmplx(0.0,0.0)
  
  allocate(mtfdata(ntfmax), ntfdata(ntfmax), stat=nalloc)
  if (nalloc/=0) then
    ierr = 2
    return
  end if
  mtfdata = 0.0
  ntfdata = 0.0
 
  return

end subroutine wavimgprm_init2



subroutine wavimgprm_initcl()
! Initializes the clock timer counts.

  implicit none

  if (wicl_runtimes==1) then ! init the run-time clock state
    call system_clock(wicl_init, wicl_rate, wicl_max)
  end if

  return

end subroutine wavimgprm_initcl




subroutine wavimgprm_getcls(cls)
! Determines the time span to the previous clock init in seconds.

  implicit none

  integer*8 :: clnow, cldelta
  real*4, intent(out) :: cls

  cls = 0.0
  if (wicl_runtimes==1) then 
    call system_clock(clnow)
    if (clnow<wicl_init) then
      cldelta = wicl_max - wicl_init + clnow
      wicl_init = clnow
    else
      cldelta = clnow - wicl_init
    end if
    cls = real(cldelta,kind=4)/wicl_rate
  end if

  return

end subroutine wavimgprm_getcls
!**********************************************************************!




END MODULE wavimgprm