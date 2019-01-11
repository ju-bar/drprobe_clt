!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  MultiSlice.f90                                        !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2009-2019   !
!                                                                      !
!**********************************************************************!
!                                                                      !
!    MODULE MultiSlice                                                 !
!    -----------------                                                 !
!                                                                      !
!    Purpose  : Implementation of a multi-slice algorithm              !
!    Version  : 1.0.0, May 07, 2007                                    !
!               1.1.0, Sept. 01, 2010                                  !
!               1.2.0, March 08, 2017                                  !
!               1.3.0, June 21, 2017                                   !
!               1.3.1, June 26, 2017                                   !
!               1.3.2, Nov 09, 2017                                    !
!               1.4.0, Jan 08, 2019                                    !
!                                                                      !
!   Linked Libs: libfftwf-3.3.lib                                      !
!   Includes   : fftw3.f03.in                                          !
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


!*********************************************************************!
!                                                                     !
!    IMORTANT REMARKS                                                 !
!    ----------------                                                 !
!                                                                     !
!    (1) CALL MS_INIT() before usage                                  !
!    (2) CALL MS_UNINIT() before re-usage or when halting             !
!                                                                     !
!*********************************************************************!



!*********************************************************************!
!*********************************************************************!
!    MODULE DECALRATIONS                                              !
!*********************************************************************!
!*********************************************************************!

MODULE MultiSlice

! Use MS as acronym for public parameters funcs and subs!
    
  ! Global module dependencies
  use, intrinsic :: iso_c_binding   
   
  implicit none

  ! Declare internal data types
  integer, parameter :: C_FFTW_R2R_KIND = C_INT ! missing declare in fftw3 include
  
  INCLUDE 'fftw3.f03.in'
  
  

  ! accessibility of subroutines or functions
!  private :: MS_***
  private :: MS_ERROR
  private :: MS_TABBED_SIGMOID
  private :: MS_ALLOCSTATICS
  private :: MS_DEALLOCSTATICS
  private :: MS_SETTAB_SCR
  private :: MS_SETTAB_USC
!  private :: MS_PropWaveElastic
  
!  public :: MS_***
  public :: MS_INIT
  public :: MS_UNINIT
  public :: MS_FFT
  public :: MS_FFT_WORK
!  public :: MS_SetAnnularDetector
  public :: MS_SetIncomingWave
  public :: MS_ShiftWave
  public :: MS_OffsetIncomingWave
  public :: MS_PrepareSlice
  public :: MS_SliceApplyBuni
  public :: MS_SlicePot2Pgr
  public :: MS_PreparePropagators
  public :: MS_SetStack
  public :: MS_Start
  public :: MS_Stop
  public :: MS_CalculateNextSlice
  public :: MS_ApplyPSpatialCoh
!  public :: MS_Calculate
  

!  declare module global (but private) variables, params and arrays
!  file units interval
  integer*4, private, parameter :: MS_minunit = 81
  integer*4, private, parameter :: MS_maxunit = 90

!   length of names and internal strings
  integer*4, public, parameter :: MS_ll = 1024
  
!   standard output unit
  integer*4, public, parameter :: MS_stdout = 6

!   enpty line placeholder (emptied at startup)
  character(len=MS_ll), private :: MS_el  
  save
! pi
  real*4, public :: MS_pi
  DATA MS_pi /3.1415927/

! scale degree to radian
  real*4, public :: MS_rd2r
  DATA MS_rd2r /0.01745329/

! error counter  
  integer*4, public :: MS_err_num
  DATA MS_err_num /0/

! error message string
  character(len=MS_ll), public :: MS_err_msg
! message string
  character(len=MS_ll), public :: MS_msg
    
! status
  integer*4, public :: MS_status
  DATA MS_status /0/
  
! debug flag
  integer*4, public :: MS_DEBUG_EXPORT
  DATA MS_DEBUG_EXPORT /0/
  integer*4, public :: MS_VERBO_EXPORT
  DATA MS_VERBO_EXPORT /0/
  
! aberration power threshold
  real*4, private, parameter :: MS_POWERTHRESH_WA = 1.0E-30
  
! aperture power threshold
  real*4, private, parameter :: MS_APERTURETHRESH = 1.0E-02
  
! hard Fourier-space aperture for multi slice ! used in the propagators
! relative to w_max
  real*4, private, parameter :: MS_RELAPERTURE = 2.0 / 3.0
  

! wave data size in one dimension  
  integer*4, public :: MS_dimx, MS_dimy
  DATA MS_dimx /256/
  DATA MS_dimy /256/
  
! slice data size in one dimension

! wavelength [nm]
  real*4, public :: MS_lamb ! this value is used by all routines
  DATA MS_lamb /0.001969/
! electron energy [keV]
  real*4, public :: MS_ht ! only kept for fun
  DATA MS_ht /300.0/


! current real-space sampling [nm/pix]
  real*4, public :: MS_samplingx, MS_samplingy
  DATA MS_samplingx /0.05/
  DATA MS_samplingy /0.05/
  
! object tilt angle to be considered by propagators [deg]
  real*4, public :: MS_objtiltx, MS_objtilty
  DATA MS_objtiltx /0.0/
  DATA MS_objtilty /0.0/
  
! detector area [mrad]
  real*4, public :: MS_detminang, MS_detmaxang  
  DATA MS_detminang /0.0/
  DATA MS_detmaxang /100.0/

! current Fourier space sampling [(pix*nm)^-1]
!  real*4, public :: MS_itog, MS_itow

! condenser aperture
  real*4, public :: MS_caperture ! [mrad]
  DATA MS_caperture /30.0/
  
  
! angular function tables
  integer*4, parameter, private :: MS_ANGTAB_SIZE = 1023
! sigmoid function table
  real*4, parameter, private :: MS_SIGMOID_EXT = 3.0
  !real*4, dimension(MS_ANGTAB_SIZE), private :: MS_ANGTAB_SIGMOID
  real*4, dimension(:), allocatable, private :: MS_ANGTAB_SIGMOID
! scramble and unscramble lists
  integer*4, dimension(:), allocatable, public :: MS_TABBED_SCR, MS_TABBED_USC
  integer*4, dimension(:), allocatable, public :: MS_TABBED_SCR2, MS_TABBED_USC2

! incoming wavefunction, and its backup
  complex*8, dimension(:,:), allocatable, public :: MS_wave_in
  complex*8, dimension(:,:), allocatable, public :: MS_wave_in_bk
! current wavefunctions
!  complex*8, dimension(:,:), allocatable, public :: MS_wave !, MS_wavei
  complex(C_FLOAT_COMPLEX), public, pointer :: MS_wave(:,:)
! average wavefunctions at different exit planes
  complex*8, dimension(:,:,:), allocatable, public :: MS_wave_avg
  integer*4, public :: MS_wave_avg_num ! numer of exit-plane waves
  DATA MS_wave_avg_num /0/
  integer*4, public :: MS_wave_avg_idx ! index of plane for average wave export (volatile and set on use)
  DATA MS_wave_avg_idx /0/
  integer*4, dimension(:), allocatable, public :: MS_wave_avg_nac ! accumulations in the average
  
! fftw data handlers
  integer*4, public, dimension(2) :: MS_fft_dims
  DATA MS_fft_dims /0, 0/
  integer(C_INT), public :: MS_fft_flags ! fft flags (0: FFTW_MEASURE, 64: FFTW_ESTIMATE)
  DATA MS_fft_flags /FFTW_MEASURE/
  type(C_PTR), public :: MS_fft_planf ! module internal forward plan
  type(C_PTR), public :: MS_fft_planb ! module internal backward plan
  type(C_PTR), public :: MS_fft_pdata ! module internal pointers to aligned data
  type(C_PTR), public :: MS_fft_work_planf ! module internal forward plan (working array)
  type(C_PTR), public :: MS_fft_work_planb ! module internal backward plan (working array)
  type(C_PTR), public :: MS_fft_work_pdata ! module internal pointers to aligned data (working array)
  complex(C_FLOAT_COMPLEX), public, pointer :: MS_work(:,:) ! pointer to working array

  
! Multi-slice stack arrays
! number of different slices
  integer*4, public :: MS_slicenum
  DATA MS_slicenum /0/
! number of different propagators
  integer*4, public :: MS_propnum
  DATA MS_propnum /0/
! size of multi-slice stack
  integer*4, public :: MS_stacksize
  DATA MS_stacksize /0/
! number of digits for slice indexing in output
  integer*4, public :: MS_nslid
  DATA MS_nslid /3/
! detection flags
  integer*4, public :: MS_useinelastic
  DATA MS_useinelastic /0/
  integer*4, public :: MS_excludeelastic
  DATA MS_excludeelastic /0/
  integer*4, public :: MS_propagateall
  DATA MS_propagateall /0/
! stack index list connection to slices
  integer*4, dimension(:), allocatable, public :: MS_slicestack
  integer*4, dimension(:), allocatable, public :: MS_propstack
! stack thickness per slice
  real*4, dimension(:), allocatable, public :: MS_slicethick
! fresnel propagators (fourier-space, scrambled and transposed)
  complex*8, dimension(:,:,:), allocatable, public :: MS_propagator

! current slice counter for external state checking
  integer*4, public :: MS_slicecur, MS_lastmaxslice
  real*4, public :: MS_calcthick
  
! flag: export of wave function data
  integer*4, public :: MS_wave_export
  DATA MS_wave_export /0/   ! OFF by default, set to >=1 to activate.
                            ! Warning: Switching ON this option will cause
                            !          an increase of calulation time.

! flag: export form of the wave function data (RS / FS)
  integer*4, public :: MS_wave_export_form
  DATA MS_wave_export_form /0/ ! 0: default -> real space wave functions
                            ! 1: Fourier space wave functions

! flag: export of average wave function data
  integer*4, public :: MS_wave_avg_export
  DATA MS_wave_avg_export /0/ ! OFF by default, set to >=1 to activate.
                            ! Warning: Switching ON this option will cause
                            !          an increase of calulation time.
                            
! flag: export of integrated probe intensities
  integer*4, public:: MS_pint_export
  DATA MS_pint_export /0/ ! OFF by default, set to >=1 to activate.
                            ! Warning: Switching ON this option will cause
                            !          an increase of calulation time.

! index: current index in the intensity array (volatile, set on use)
  integer*4, public :: MS_pint_idx
  DATA MS_pint_idx /0/

! flag: export of the incident wave function data
  integer*4, public :: MS_incwave_export
  DATA MS_incwave_export /0/  ! OFF by default, set to >=1 to activate.

! wave function file names 
  character(len=MS_ll), public :: MS_wave_filenm, MS_wave_filenm_bk
  character(len=MS_ll), public :: MS_wave_filenm_avg
  DATA MS_wave_filenm /"epw"/
  DATA MS_wave_filenm_avg /"epw"/
! period of wave export in number of slices  
  integer*4, public :: MS_wave_export_pzp
  DATA MS_wave_export_pzp /1/
  

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
SUBROUTINE MS_INIT()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 100
  
!
  integer*4 :: m , err
  real*4 :: angle, ascale
  complex*8 :: c0
  character(len=400) :: smsg
  real*4, external :: sigmoid ! BasicFuncs.f90
  
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_INIT: INIT."
! ------------

! ------------
  MS_pi = atan(1.0)*4.0
  MS_rd2r = MS_pi/180.0
  MS_err_num = 0
  MS_el = REPEAT(" ",MS_ll)
  MS_err_msg = MS_el
  MS_status = 0
  MS_slicecur = -1
  c0 = cmplx(0.0,0.0)
  call MS_ALLOCSTATICS(err)
  if (err/=0) then
    call MS_ERROR("MultiSlice init failed.",subnum+1)
    return
  end if
  MS_propnum = 0
  if (allocated(MS_propagator)) deallocate(MS_propagator,stat=err)
! ------------


! ------------
! prepare angular tables now
  ascale = 2.0*MS_SIGMOID_EXT/real(MS_ANGTAB_SIZE-1)
  do m=1,MS_ANGTAB_SIZE
    angle = real(m-1)*ascale-MS_SIGMOID_EXT
    MS_ANGTAB_SIGMOID(m) = sigmoid(angle,0.0,1.0)
  end do
! ------------
  
  

! ------------
  if (MS_slicenum>0.and.MS_stacksize>0.and.MS_dimx>0.and.MS_dimy>0) then
    ! stack index list connection to slices
    err = 0
    allocate(MS_slicestack(1:MS_stacksize),STAT=err)
    if (err/=0) then
      call MS_ERROR("MultiSlice init, failed to allocate slice stack.",subnum+1)
      return
    end if
    MS_slicestack = 0
    ! stack hash of propagators
    allocate(MS_propstack(1:MS_slicenum),STAT=err)
    if (err/=0) then
      call MS_ERROR("MultiSlice init, failed to allocate propagator stack.",subnum+3)
      return
    end if
    MS_propstack = 0
    ! thickness array
    allocate(MS_slicethick(1:MS_slicenum),STAT=err)
    if (err/=0) then
      call MS_ERROR("MultiSlice init failed.",subnum+3)
      return
    end if
    MS_slicethick = 0.0;

    ! call here, be aware that changing MS_dim needs REINIT!
    ! Thus, no need to call MS_SETTAB_*** elsewhere    
    call MS_SETTAB_SCR(MS_dimx, MS_dimy)
    call MS_SETTAB_USC(MS_dimx, MS_dimy)
    
    ! prepare plans for DFT
    if (MS_fft_dims(1)/=MS_dimy .or. MS_fft_dims(2)/=MS_dimx) then
      if ( MS_fft_dims(1)>0 .or. MS_fft_dims(2)>0 ) then
        call fftwf_free(MS_fft_pdata)
        call fftwf_destroy_plan(MS_fft_planf)
        call fftwf_destroy_plan(MS_fft_planb)
        call fftwf_free(MS_fft_work_pdata)
        call fftwf_destroy_plan(MS_fft_work_planf)
        call fftwf_destroy_plan(MS_fft_work_planb)
        MS_fft_dims = (/0, 0/)
      end if
      MS_fft_pdata = fftwf_alloc_complex( int(MS_dimy*MS_dimx, C_SIZE_T) )
      MS_fft_work_pdata = fftwf_alloc_complex( int(MS_dimy*MS_dimx, C_SIZE_T) )
      call c_f_pointer(MS_fft_pdata, MS_wave, [MS_dimx,MS_dimy])
      call c_f_pointer(MS_fft_work_pdata, MS_work, [MS_dimx,MS_dimy])
      MS_fft_planf = fftwf_plan_dft_2d(MS_dimy, MS_dimx, MS_wave, MS_wave, FFTW_FORWARD, MS_fft_flags)
      MS_fft_planb = fftwf_plan_dft_2d(MS_dimy, MS_dimx, MS_wave, MS_wave, FFTW_BACKWARD, MS_fft_flags)
      MS_fft_work_planf = fftwf_plan_dft_2d(MS_dimy, MS_dimx, MS_work, MS_work, FFTW_FORWARD, MS_fft_flags)
      MS_fft_work_planb = fftwf_plan_dft_2d(MS_dimy, MS_dimx, MS_work, MS_work, FFTW_BACKWARD, MS_fft_flags)
      MS_fft_dims = (/ MS_dimy, MS_dimx /)
    end if
    ! set ready state
    MS_status = 1
  else
    smsg = ""
    if (MS_slicenum<=0) then
      smsg = trim(smsg)//" number of slice files"
    end if
    if (MS_stacksize<=0) then
      if (len_trim(smsg)>0) smsg = trim(smsg)//","
      smsg = trim(smsg)//" number of object slices"
    end if
    if (MS_dimx<=0) then
      if (len_trim(smsg)>0) smsg = trim(smsg)//","
      smsg = trim(smsg)//" slice x-dimension"
    end if
    if (MS_dimy<=0) then
      if (len_trim(smsg)>0) smsg = trim(smsg)//","
      smsg = trim(smsg)//" slice y-dimension"
    end if
    call MS_ERROR("MultiSlice init failed, invalid parameter:"//trim(smsg),subnum+8)
    return
  end if
! ------------

  return

END SUBROUTINE MS_INIT
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_UNINIT()
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
!  write(unit=*,fmt=*) " > MS_INIT: MS_UNINIT."
  err = 0
  MS_err_num = 0
  MS_status = 0
! ------------

! ------------
  call MS_DEALLOCSTATICS(err)
  if (allocated(MS_slicestack)) then
    deallocate(MS_slicestack)
  end if
  if (allocated(MS_propstack)) then
    deallocate(MS_propstack,stat=err)
  end if
  if (allocated(MS_slicethick)) then
    deallocate(MS_slicethick)
  end if
  if (allocated(MS_propagator)) then
    deallocate(MS_propagator)
  end if
  if ( MS_fft_dims(1)>0 .or. MS_fft_dims(2)>0 ) then
    call fftwf_free(MS_fft_pdata)
    call fftwf_destroy_plan(MS_fft_planf)
    call fftwf_destroy_plan(MS_fft_planb)
    call fftwf_free(MS_fft_work_pdata)
    call fftwf_destroy_plan(MS_fft_work_planf)
    call fftwf_destroy_plan(MS_fft_work_planb)
    MS_fft_dims = (/0, 0/)
  end if
! ------------
!  write(unit=*,fmt=*) " > MS_UNINIT: EXIT."
  return

END SUBROUTINE MS_UNINIT
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_ALLOCSTATICS(nerr)
! function: allocates static size arrays in heap
! -------------------------------------------------------------------- !
! parameter:
!       integer*4 :: nerr       = error code (0=success)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1700
  integer*4, intent(inout) :: nerr
  integer*4 :: ierr                         ! internal error counter
  integer*4 :: asc                          ! allocation status check
  complex*8 :: c0
! ------------

! ------------
! INIT
!  write(*,*) " > MS_ALLOCSTATICS: INIT."
  nerr = 0
  ierr = 0
  asc = 0
  c0 = cmplx(0.0,0.0)
! ------------

! ------------
! since these arrays are always of same size, no deallocation is required
! if they are already allocated, use them as is
! --------------
  if (.not.allocated(MS_ANGTAB_SIGMOID)) then
    allocate(MS_ANGTAB_SIGMOID(MS_ANGTAB_SIZE),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (.not.allocated(MS_TABBED_SCR)) then
    allocate(MS_TABBED_SCR(1:MS_dimx),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (.not.allocated(MS_TABBED_USC)) then
    allocate(MS_TABBED_USC(1:MS_dimx),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (.not.allocated(MS_TABBED_SCR2)) then
    allocate(MS_TABBED_SCR2(1:MS_dimy),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (.not.allocated(MS_TABBED_USC2)) then
    allocate(MS_TABBED_USC2(1:MS_dimy),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (.not.allocated(MS_wave_in)) then
    allocate(MS_wave_in(MS_dimx,MS_dimy),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (.not.allocated(MS_wave_in_bk)) then
    allocate(MS_wave_in_bk(MS_dimx,MS_dimy),stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! ------------

! ------------
  if (ierr/=0) then
    nerr = 1
    return
  end if
! ------------

! ------------
  MS_ANGTAB_SIGMOID = 0.0
  MS_TABBED_SCR = 0
  MS_TABBED_USC = 0
  MS_TABBED_SCR2 = 0
  MS_TABBED_USC2 = 0
  MS_wave_in = c0
  MS_wave_in_bk = c0
  MS_wave = c0
  !MS_wavei = c0
! ------------

! ------------
!  write(*,*) " > MS_ALLOCSTATICS: EXIT."
  return

END SUBROUTINE MS_ALLOCSTATICS
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_DEALLOCSTATICS(nerr)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1800
  integer*4, intent(inout) :: nerr
  integer*4 :: ierr, asc
! ------------

! ------------
! INIT
!  write(*,*) " > MS_DEALLOCSTATICS: INIT."
  nerr = 0
  ierr = 0
  asc = 0
! ------------

! ------------
  if (allocated(MS_ANGTAB_SIGMOID)) then
    deallocate(MS_ANGTAB_SIGMOID,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (allocated(MS_TABBED_SCR)) then
    deallocate(MS_TABBED_SCR,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (allocated(MS_TABBED_USC)) then
    deallocate(MS_TABBED_USC,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (allocated(MS_TABBED_SCR2)) then
    deallocate(MS_TABBED_SCR2,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (allocated(MS_TABBED_USC2)) then
    deallocate(MS_TABBED_USC2,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (allocated(MS_wave_in)) then
    deallocate(MS_wave_in,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! --------------
  if (allocated(MS_wave_in_bk)) then
    deallocate(MS_wave_in_bk,stat=asc)
    if (asc/=0) then
      ierr = ierr + 1
      asc = 0
    end if
  end if
! ------------
  

! ------------
  if (ierr/=0) then
    nerr = 1
    return
  end if
! ------------

! ------------
!  write(*,*) " > MS_DEALLOCSTATICS: EXIT."
  return

END SUBROUTINE MS_DEALLOCSTATICS
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_ERROR(sTxt,nErr)
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
  character(len=MS_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = MS_el
! ----------

! ----------
! messaging (to screen)
  write(unit=sinfo,fmt=*) "ERROR: ",trim(sTxt)," Error code:",nErr
  MS_err_msg = sinfo ! save last error message
!  call SE_event(trim(sinfo), SE_err)
  write(unit=MS_stdout,fmt='(A)') "MS_ERROR: "//trim(sinfo)
! ----------

  MS_err_num = MS_err_num + 1
  return

END SUBROUTINE MS_ERROR
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_TABBED_SIGMOID(x,x0,dx,rval)
! function: returns tabbed real*4 sigmoid(x,x0,dx)
! -------------------------------------------------------------------- !
! parameter: real*4 :: x,x0,dx,rval
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL MS_INIT() BEFORE CALLING THIS ROUTINE !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 400
  real*4 :: x,x0,dx,rval
  real*4 :: tx, ix, maxx
  integer*4 :: idx
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_TABBED_EXP: INIT."
  maxx = MS_SIGMOID_EXT
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
  ix = 0.5*(tx/MS_SIGMOID_EXT+1.0)*real(MS_ANGTAB_SIZE-1)
  idx = 1+int(ix) ! round to next integer
  if (idx==0) then
    idx = MS_ANGTAB_SIZE
  end if
  if (idx==MS_ANGTAB_SIZE+1) then
    idx = 1
  end if
! clip angle to interval 0.. MS_ANGTAB_SIZE
! ------------

! ------------
! get data
  rval = MS_ANGTAB_SIGMOID(idx)
! ------------

! ------------
!  write(unit=*,fmt=*) " > MS_TABBED_EXP: EXIT."
  return

END SUBROUTINE MS_TABBED_SIGMOID
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_SETTAB_SCR(nx, ny)
! function: presets scramble tabs
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL BEFORE USING MS_TABBED_SCR !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 500
  integer*4, intent(in) :: nx, ny
  integer*4 :: n2, n2m1
  integer*4 :: i
! ------------

! ------------
! HORIZONTALS (X)
  if (nx > 0 .and. nx <= MS_dimx) then
    n2 = nx/2
    n2m1 = n2-1
    MS_TABBED_SCR = 1
    do i=1, nx
      MS_TABBED_SCR(i)=mod((i+n2m1),nx)-n2
    end do
  end if
! ------------
  
! ------------
! VERTICALS (Y)
  if (ny > 0 .and. ny <= MS_dimy) then
    n2 = ny/2
    n2m1 = n2-1
    MS_TABBED_SCR2 = 1
    do i=1, ny
      MS_TABBED_SCR2(i)=mod((i+n2m1),ny)-n2
    end do
  end if
! ------------

! ------------
  return

END SUBROUTINE MS_SETTAB_SCR
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_SETTAB_USC(nx, ny)
! function: presets unscramble tabs
! -------------------------------------------------------------------- !
!!!!!!! BE SHURE TO CALL BEFORE USING MS_TABBED_SCR !!!!!!!!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 500
  integer*4, intent(in) :: nx, ny
  integer*4 :: n2
  integer*4 :: i
! ------------

! ------------
! HORIZONTALS (X)
  if (nx > 0 .and. nx <= MS_dimx) then
    n2 = nx/2
    MS_TABBED_USC = 1
    do i=1, nx
      MS_TABBED_USC(i)=mod(i-1+n2,nx)+1
    end do
  end if
! ------------
  
! ------------
! VERTICALS (Y)
  if (ny > 0 .and. ny <= MS_dimy) then
    n2 = ny/2
    MS_TABBED_USC2 = 1
    do i=1, ny
      MS_TABBED_USC2(i)=mod(i-1+n2,ny)+1
    end do
  end if
! ------------

! ------------
  return

END SUBROUTINE MS_SETTAB_USC
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_SetIncomingWave(wave)
! function: copies incoming wave data to internal backup memory
! -------------------------------------------------------------------- !
! parameter: complex*8, dimension (MS_dimx, MS_dimy) :: wave
! -------------------------------------------------------------------- !

  implicit none

! DECLARATION
  integer*4, parameter :: subnum = 700
  complex*8, intent(in), dimension(MS_dimx,MS_dimy) :: wave

  MS_wave_in_bk = wave

  return

END SUBROUTINE MS_SetIncomingWave
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_OffsetIncomingWave(dx,dy,dz)
! function: generates new MS_wave_in from MS_wave_in_bk with given
!           offset in nanometres
! -------------------------------------------------------------------- !
! parameter: real*4 :: dx, dy , offset wave shift
!            real*4 :: dz , defocus
!            integer*4 :: thread number
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  real*4 :: dx, dy, dz
  real*4 :: pfac, wx, wy, w2y, w2
  real*4 :: chi, chiy, itowx, itowy, ffac
  integer*4 :: i, j, nx, ny
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_OffsetIncomingWave: INIT."
! catch zero offset
  if (dx==0.0 .and. dy==0.0 .and. dz==0.0) then
    MS_wave_in = MS_wave_in_bk
    return
  end if
! clear old data
  MS_wave_in = cmplx(0.0,0.0)
! prepare parameters
  nx = MS_dimx
  ny = MS_dimy
! set Fourier-space sampling
  itowx = MS_lamb / (MS_samplingx*real(nx))
  itowy = MS_lamb / (MS_samplingy*real(ny))
  pfac = 2.0*MS_pi/MS_lamb
  ffac = 0.5*dz
! ------------

! ------------
! multiply shifting and defocussing phase plate
  do j=1, ny
    wy = MS_TABBED_SCR2(j)*itowy
    chiy = wy*dy
    w2y = wy*wy
    do i=1, nx
      wx = MS_TABBED_SCR(i)*itowx
      w2 = w2y + wx*wx
      chi = pfac*(wx*dx+chiy+ffac*w2)
      cval = MS_wave_in_bk(i,j)*cmplx(cos(chi),sin(-chi))
      MS_wave_in(i,j) = cval
    end do
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > MS_OffsetIncomingWave: EXIT."
  return

END SUBROUTINE MS_OffsetIncomingWave
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_ShiftWave(wave,dx,dy)
! function: shifts the wave function by the given ammount in real space
!           offset in nanometres
! -------------------------------------------------------------------- !
! parameter: real*4 :: dx, dy , offset wave shift
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  real*4, intent(in) :: dx, dy
  complex*8, intent(inout) :: wave(MS_dimx,MS_dimy)
  
  real*4 :: rsca, pfac, wx, wy
  real*4 :: chi, chiy, itowx, itowy
  integer*4 :: i, j, nx, ny
  complex*8 :: cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_ShiftWave: INIT."
! catch zero offset
  if (dx==0.0 .and. dy==0.0) return
  if (MS_status==0) return
! prepare parameters
  nx = MS_dimx
  ny = MS_dimy
! set Fourier-space sampling
  itowx = MS_lamb / (MS_samplingx*real(nx))
  itowy = MS_lamb / (MS_samplingy*real(ny))
  pfac = 2.0*MS_pi/MS_lamb
  rsca = 1.0 / real(nx*ny)
! ------------

! copy data
  MS_work = wave
! transform to fourier space
  call fftwf_execute_dft(MS_fft_work_planf,MS_work,MS_work)

! ------------
! apply the shifting phase plate
  do j=1, MS_dimy
    wy = MS_TABBED_SCR2(j)*itowy
    chiy = wy*dy
    do i=1, MS_dimx
      wx = MS_TABBED_SCR(i)*itowx
      chi = pfac*(wx*dx+chiy)
      cval = cmplx(cos(chi),sin(-chi))*rsca
      MS_work(i,j) = MS_work(i,j)*cval
    end do
  end do
! ------------

! transform to real space
  call fftwf_execute_dft(MS_fft_work_planf,MS_work,MS_work)
  
! copy data
  wave(i,j) = MS_work(i,j)

! ------------
!  write(unit=*,fmt=*) " > MS_ShiftWave: EXIT."
  return

END SUBROUTINE MS_ShiftWave
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
!SUBROUTINE MS_PrepareSlice(nidx, sdata, sthick)
SUBROUTINE MS_PrepareSlice(nidx, sthick)
! function: Prepares array data for a slice
! -------------------------------------------------------------------- !
! parameter: ndix : integer*4 : slice index in memory
!!            sdata(MS_dimy,MS_dimx) : complex*8 : slice data EXP[i*V]
!            sthick : real*4 : slice thickness in [nm]
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 900
  integer*4 :: nidx
  real*4 :: sthick
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_PrepareSlice: INIT."
  if (MS_status<1) then
    call MS_ERROR("Module not initialised.",subnum+1)
    return
  end if
  if (nidx<1.or.nidx>MS_slicenum) then
    call MS_ERROR("Wrong parameter (nidx) calling MS_PrepareSlice.",subnum+2)
    return
  end if
! set slice thickness
  MS_slicethick(nidx) = sthick
! ------------


! ------------
!  write(unit=*,fmt=*) " > MS_PrepareSlice: EXIT."
  return

END SUBROUTINE MS_PrepareSlice
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_SliceApplyBuni( vbuni, ndimx, ndimy, cpot, nerr)
! funtion: applies dampening by DWF to a projected potential
! -------------------------------------------------------------------- !
! parameters:
! real*4 :: vbuni = Debye parameter for Debye-Waller factor in nm^2
! integer*4 :: ndimx, ndimy = potential data x,y dimension
! complex*8 :: cpot(ndimx,ndimy) = projected potential
!              (input) = projected potential
!              (output) = dampened projected potential
! integer*4 :: nerr = routine error code
! -------------------------------------------------------------------- !
  
  implicit none
  
! constant parameters
  integer*4, parameter :: subnum = 2200
  real*4, parameter :: twopi = 6.2831853
! function parameters
  real*4, intent(in) :: vbuni
  integer*4, intent(in) :: ndimx, ndimy
  complex*8, intent(inout) :: cpot(1:ndimx,1:ndimy)
  integer*4, intent(inout) :: nerr
! working variables
  integer*4 :: i, j, i1, j1 ! iterators
  integer*4 :: nx2, ny2 ! nyquist numbers
  real*4 :: itogx, itogy ! Fourier-space sampling rate of potentials
  real*4 :: gx, gy, gy2, g2 ! spatial frequencies
  real*4 :: dwa, dwf, rsca ! Debye-Waller factor and temp values

! ---------------
! init
  nerr = 0
  if (MS_status<1) goto 101
  if (ndimx<=0 .or. ndimy<=0) goto 102
  nx2 = (ndimx-modulo(ndimx,2))/2
  ny2 = (ndimy-modulo(ndimy,2))/2
  itogx = 1./(MS_samplingx*real(ndimx))
  itogy = 1./(MS_samplingy*real(ndimy))
  dwa = -0.25 * vbuni
  rsca = 1.0 / real(ndimx*ndimy)
! ---------------

! ---------------
! transform the potential to Fourier space
! (using the slow variant since the case is rare and cpot must not be on wave dimensions)
  call MS_FFT(cpot, ndimx, ndimy, 1)
! ---------------

! ---------------
! apply the DWF
! FT-potential is scrambled and transposed in MS_sc
  do j=1, ndimy ! get y-frequency
    j1 = modulo(j - 1 + ny2, ndimy) - ny2
    gy = itogy*real(j1)
    gy2 = gy*gy
    do i=1, ndimx ! get x-frquency
      i1 = modulo(i - 1 + nx2, ndimx) - nx2
      gx = itogx*real(i1)
      g2 = gy2+gx*gx
      dwf = exp( dwa * g2 ) * rsca
      cpot(i,j) = cpot(i,j) * dwf
    end do
  end do
! ---------------

! ---------------
! transform the potential back to real space
  call MS_FFT(cpot, ndimx, ndimy, -1)
! ---------------

  goto 1000

! ---------------
! error handling
101 nerr = 1
  call MS_ERROR("MS_SliceApplyBuni: Module not initialized.", subnum+nerr)
  goto 1000
102 nerr = 2
  call MS_ERROR("MS_SliceApplyBuni: Invalid size of potential grid.", subnum+nerr)
  goto 1000

! ---------------
! routine exit
1000 continue
  return
  
END SUBROUTINE MS_SliceApplyBuni
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_SlicePot2Pgr( ht, dz, nabf, vabf, nbuni, vbuni, &
                          & ndimx, ndimy, ndimz, slcdat, nerr )
! funtion: transforms projected potentials into phase gratings
!          applies absorption depending on nabf
!          applies debye-waller factors on nbuni
! -------------------------------------------------------------------- !
! parameters:
! real*4 :: ht = electron energy keV
! real*4 :: sz = slice thickness [nm]
! integer*4 :: nabf = flag to apply absorption factor vabf
! integer*4 :: nbuni = flag to apply universal Debye-Waller factor (buni)
! real*4 :: vabf = absorption factor
! real*4 :: vbuni = Debye parameter for Debye-Waller factor in nm^2
! integer*4 :: ndimx, ndimy, ndimz = slice data x,y,z dimension
! complex*8 :: slcdat(ndimx,ndimy,ndimz) = slice data
!              (input) = projected potential
!              (output) = phase grating
! integer*4 :: nerr = routine error code
! -------------------------------------------------------------------- !
  
  implicit none
  
! constants
  integer*4, parameter :: subnum = 2100
  real*4, parameter :: elwl  = 1.23984197396 ! c*h/e [nm/kV]
  real*4, parameter :: erest = 510.998946269 ! m0*c^2 [keV] electron rest energy in keV
  real*4, parameter :: psig = 2.0886573708  ! m0 * e / (2*Pi * hbar^2) * (10^-9)^2  [ V^-1 nm^-2 ]
! function parameters
  real*4, intent(in) :: ht, dz, vabf, vbuni
  integer*4, intent(in) :: nabf, nbuni, ndimx, ndimy, ndimz
  complex*8, intent(inout) :: slcdat(1:ndimx,1:ndimy,1:ndimz)
  integer*4, intent(inout) :: nerr
! working parameters
  integer*4 :: i ! iterators
  integer*4 :: nalloc ! allocation status
  real*4 :: wl ! electon wavelength [nm]
  real*4 :: relati ! reltivistic correction factor
  real*4 :: sigmae ! interaction constant
  complex*8 :: cimag ! imaginary number
  complex*8 :: cabsorb ! absorption factor
  complex*8 :: prefac ! cross-section and other factors.
  complex*8, allocatable :: ctmp(:,:) ! working array

! ---------------  
! init
  nerr = 0
  nalloc = 0
  cabsorb = cmplx(1.0, 0.0) ! default, no absorption
  if (nabf==1) cabsorb = cmplx(1.0, vabf) ! absorption potential transfer to imaginary part
  cimag = cmplx(0.0, 1.0) ! factor to transform *I (imaginary constant)
  wl = elwl / sqrt( ht * ( 2.0*erest + ht ) ) ! lambda, electron wavelength [nm]
  sigmae = psig * wl ! interaction constant [nm-1]
  relati = 1.0 + ht / erest
  prefac = sigmae * dz * relati * cabsorb * cimag ! = gamma * sigma * t * (1+I*abf) * (I)
! ---------------

! ---------------
! allocate
  allocate( ctmp(1:ndimx,1:ndimy), stat=nalloc )
  if (nalloc/=0) goto 101
  ctmp = cmplx(0.0,0.0)
! ---------------

! ---------------
! calculations
  do i=1, ndimz ! loop over all slices in the stack
    ctmp(1:ndimx,1:ndimy) = slcdat(1:ndimx,1:ndimy,i)
  ! -----------------
    if (nbuni/=0) then ! apply Debye-Waller factor
      write(unit=MS_msg,fmt=*) vbuni
      call PostDebugMessage( "Applying global Debye-Waller factors with Biso = "// &
         & trim(adjustl(MS_msg))//" nm^2." )
      call MS_SliceApplyBuni( vbuni, ndimx, ndimy, ctmp, nerr)
      if (nerr/=0) goto 102
    end if 
  ! -----------------

  ! -----------------
  ! transform & transfer back
  ! here is the actual potential-to-phase-grating transformation
    slcdat(1:ndimx,1:ndimy,i) = EXP( prefac * ctmp(1:ndimx,1:ndimy) )
  ! -----------------

  end do ! loop (i) over all slices
! ---------------

  goto 1000

! ---------------
! error handling
101 nerr = 1
  call MS_ERROR("MS_SlicePot2Pgr: Memory allocation failed.", subnum+nerr)
  goto 1000
102 nerr = 2
  call MS_ERROR("MS_SlicePot2Pgr: Failed to apply universal Debye-Waller factor.", subnum+nerr)
  goto 1000

! ---------------
! routine exit
1000 continue
  if (allocated(ctmp)) deallocate(ctmp,stat=nalloc)
  return
  
END SUBROUTINE MS_SlicePot2Pgr
!**********************************************************************!


!------------------------------------------------------------------------
!
! This is the old version of Fresnel propagation.
!
!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_PreparePropagators()
! function: Prepares all propagators for calculation
!           do not forget to call this function pefore starting the
!           multisclice algorithm
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2000
  real*8, parameter :: dd2r = 0.1745329252D-001
  real*8, parameter :: dpi  = 3.141592653589D+000
  
  integer*4 :: nx, ny, i, j, k, err, addid, matchid
  real*8 :: pfac, wx, wy, wy2, w2, ux, uy, uy2, u2
  real*8 :: chi, itowx, itowy, itogx, itogy !, power
  real*8 :: wthresh, wmax, otx, oty, pt0, rsca
  complex*8 :: cval, cval0
  character(len=400) :: stmp
  real*4 :: isthick(MS_slicenum)
! ------------

! ------------
! INIT
!  write(*,*) " > MS_PreparePropagators: INIT."
  if (MS_status<1) goto 101
  if (MS_slicenum<1) goto 102
! prepare parameters
  nx = MS_dimx
  ny = MS_dimy
! set Fourier-space sampling
  itogx = 1.0D+0/dble(MS_samplingx*nx)
  itowx = itogx*dble(MS_lamb)
  itogy = 1.0D+0/dble(MS_samplingy*ny)
  itowy = itogy*dble(MS_lamb)
  wmax = min(itowx*nx/2,itowy*ny/2)
  wthresh = dble(MS_RELAPERTURE)*wmax
  wthresh = wthresh*wthresh
  cval0 = cmplx(0.0,0.0)
  otx = dble(MS_objtiltx)* dd2r ! object tilt x in radian
  oty = dble(MS_objtilty)* dd2r ! object tilt y in radian
  !MS_slicethick(nidx) = sthick
  pfac = dpi/dble(MS_lamb)*dcos(dsqrt(otx*otx+oty*oty))
  rsca = 1.0D+0/dble(nx*ny)
! ------------


! ------------
! determine number of distinguished propagators
  isthick = 0.0
  MS_propstack = 0
  MS_propnum = 0
  addid = 0
  if (allocated(MS_propagator)) deallocate(MS_propagator,stat=err)
  do i=1, MS_slicenum  ! loop through distinguished slices
    addid = 0 ! reset id to add
    
    if (MS_propnum==0) then ! add current slice since stack is empty
      addid = i
    else ! stack is not empty try to find match
      k = 0
      matchid = 0 ! reset matching id
      do j=1, MS_propnum
        if ( MS_slicethick(i)==isthick(j) ) then
          matchid = j
          exit ! match found exit now, j holds matching index in propagator list
        else
          k = k + 1 ! not a match, raise count
        end if
      end do
      if (matchid==0) then ! no match found
        addid = i ! save current slice id for adding a new propagator
      end if
    end if
    if (addid>0) then ! add a new propagator
      MS_propnum = MS_propnum + 1
      matchid = MS_propnum
      isthick(matchid) = MS_slicethick(i)
    end if
    MS_propstack(i) = matchid ! set index in propagator stack
    
  end do  ! loop through distinguished slices
  if (MS_propnum<=0) goto 103
! ------------

! ------------
! allocate propagator array
! fresnel propagators (fourier-space, scrambled and transposed)
  err = 0
  allocate(MS_propagator(1:MS_dimx,1:MS_dimy,1:MS_propnum),STAT=err)
  if (err/=0) goto 104
  MS_propagator = cval0
! ------------

! ------------
! setup fresnel propagators
  pt0 = pfac*(otx*otx+oty*oty) ! Pi/lambda * Cos[tilt] * t^2
  do k=1, MS_propnum ! L+ all propagators
  do j=1, ny ! L+ all lines ! y frequencies
    wy = dble(MS_TABBED_SCR2(j))*itowy ! wy = ky*lambda
    wy2 = wy*wy
    uy = wy - oty ! uy -> tilted
    uy2 = uy*uy
    do i=1, nx ! L+ all columns ! x frequencies
      wx = dble(MS_TABBED_SCR(i))*itowx ! wx = kx*lambda
      w2 = wy2+wx*wx
      ux = wx - otx ! ux -> tilted
      u2 = uy2+ux*ux
      if (w2>wthresh) then ! hard aperture in wave frame
        MS_propagator(i,j,k) = cval0
      else
        chi = dble(isthick(k))*(pfac*u2-pt0) ! chi(k) = dZ * Pi/lambda * Cos[tilt] * (w^2 - 2 w.t)
                                             !                                        ^Prop   ^Tilt
                                             ! The tilted propagator effectively shifts the slices
                                             ! against each other along the tilt direction.
        cval = cmplx(dcos(chi),-dsin(chi),4)*rsca
        MS_propagator(i,j,k) = cval !cmplx(cos(chi),sin(chi),4)
      end if
    end do ! L- all columns
  end do ! L- all lines
  if (MS_DEBUG_EXPORT/=0) then
    ! debug output of fresnel propagator functions
    write(unit=stmp,fmt='("fre_",I3.3,".dat")') k
    call SaveDataC8(trim(stmp), MS_propagator(1:nx,1:ny,k),nx*ny,err)
  end if
  end do ! L- all propagators
! ------------

  goto 1000

! ------------
! error handling
101 call MS_ERROR("Module not initialised.",subnum+1)
  goto 1000
102 call MS_ERROR("No slices specified.", subnum+2)
  goto 1000
103 call MS_ERROR("No propagators, zero object thickness.",subnum+3)
  goto 1000
104 call MS_ERROR("MultiSlice propagator allocation failed.",subnum+4)
  goto 1000

! ------------
1000 continue
!  write(*,*) " > MS_PreparePropagators: EXIT."
  return

END SUBROUTINE MS_PreparePropagators
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_PreparePropagators2()
! function: Prepares all propagators for calculation
!           do not forget to call this function pefore starting the
!           multisclice algorithm
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2000
  real*8, parameter :: dd2r = 0.1745329252D-001
  real*8, parameter :: dpi  = 3.141592653589D+000
  
  integer*4 :: nx,ny, i, j, k, err, addid, matchid
  real*8 :: rsca, pfac, gx, gy, gy2, g2
  real*8 :: chi, itogx, itogy ! , power !, itowx, itowy, 
  real*8 :: gthresh, gmax, otx, oty !, pt0
  real*8 :: ot, od, dt, dd ! object and diffraction tilt angles and directions [rad]
  real*8 :: cot, sot ! cosine and sine of the object tilt magnitude
  complex*8 :: cval, cval0
  character(len=400) :: stmp
  real*4 :: isthick(MS_slicenum)
! ------------

! ------------
! INIT
!  write(*,*) " > MS_PreparePropagators2: INIT."
  if (MS_status<1) goto 101
  if (MS_slicenum<1) goto 102
! prepare parameters
  nx = MS_dimx
  ny = MS_dimy
! set Fourier-space sampling
  itogx = 1.0D+0 / dble(MS_samplingx*nx)
!  itowx = itogx*dble(MS_lamb)
  itogy = 1.0D+0 / dble(MS_samplingy*ny)
!  itowy = itogy*dble(MS_lamb)
  gmax = min(itogx*nx/2,itogy*ny/2)
  gthresh = dble(MS_RELAPERTURE)*gmax
  gthresh = gthresh*gthresh
  cval0 = cmplx(0.0,0.0)
  otx = dble(MS_objtiltx)* dd2r ! object tilt x [rad]
  oty = dble(MS_objtilty)* dd2r ! object tilt y [rad]
  !MS_slicethick(nidx) = sthick
  !pfac = dpi/dble(MS_lamb)*dcos(dsqrt(otx*otx+oty*oty))
  pfac = 2.0D+0 * dpi / dble(MS_lamb)
  rsca = 1.0D+0 / dble(nx*ny) ! rescaling to keep the norm
  ot = dsqrt(otx*otx+oty*oty) ! object tilt magnitude [rad]
  cot = dcos(ot)
  sot = dsin(ot)
  od = 0.0
  if (ot>0.D+0) od = datan2(oty,otx) ! object tilt direction [rad]
! ------------


! ------------
! determine number of distinguished propagators
  isthick = 0.0
  MS_propstack = 0
  MS_propnum = 0
  addid = 0
  if (allocated(MS_propagator)) deallocate(MS_propagator,stat=err)
  do i=1, MS_slicenum  ! loop through distinguished slices
    addid = 0 ! reset id to add
    if (MS_propnum==0) then ! add current slice since stack is empty
      addid = i
    else ! stack is not empty try to find match
      k = 0
      matchid = 0 ! reset matching id
      do j=1, MS_propnum
        if ( MS_slicethick(i)==isthick(j) ) then
          matchid = j
          exit ! match found exit now, j holds matching index in propagator list
        else
          k = k + 1 ! not a match, raise count
        end if
      end do
      if (matchid==0) then ! no match found
        addid = i ! save current slice id for adding a new propagator
      end if
    end if
    if (addid>0) then ! add a new propagator
      MS_propnum = MS_propnum + 1
      matchid = MS_propnum
      isthick(matchid) = MS_slicethick(i)
    end if
    MS_propstack(i) = matchid ! set index in propagator stack
    
  end do  ! loop through distinguished slices
  if (MS_propnum<=0) goto 103
! ------------

! ------------
! allocate propagator array
! fresnel propagators (fourier-space, scrambled and transposed)
  err = 0
  allocate(MS_propagator(1:MS_dimx,1:MS_dimy,1:MS_propnum),STAT=err)
  if (err/=0) goto 104
  MS_propagator = cval0
! ------------

! ------------
! setup fresnel propagators
!  pt0 = pfac*(otx*otx+oty*oty) ! Pi/lambda * Cos[tilt] * t^2
  do k=1, MS_propnum ! L+ all propagators
  do j=1, ny ! L+ all lines ! y frequencies
    gy = dble(MS_TABBED_SCR2(j))*itogy ! gy
    gy2 = gy*gy
    do i=1, nx ! L+ all columns ! x frequencies
      gx = dble(MS_TABBED_SCR(i))*itogx ! gx
      g2 = gy2+gx*gx
      if (g2>gthresh) then ! hard aperture in wave frame
        MS_propagator(i,j,k) = cval0 ! no transfer beyond aperture
      else
        !dt = 2.D+0*dasin( 0.5D+0 * dble(MS_lamb) * dsqrt(g2) ) ! |theta|
        dt = dasin( dble(MS_lamb) * dsqrt(g2) ) ! |theta|
        dd = 0.D+0
        if (dt>0.D+0) dd = datan2(gy,gx) ! theta_phi
        ! calculate phase shift to non-diffracted beam
        chi = pfac*dble(isthick(k))*(  &
              &   1.D+0/(cot*dcos(dt) + sot*dcos(od-dd)*dsin(dt) ) &
              & - 1.D+0/cot )
        ! calculate complex phase factor to wavefunction, exp(-I*chi)
        cval = cmplx( dcos(chi), -dsin(chi), 4) * rsca
        ! store the phase factor
        MS_propagator(i,j,k) = cval
        !
      end if
    end do ! L- all columns
  end do ! L- all lines
  if (MS_DEBUG_EXPORT/=0) then
    ! debug output of propagator functions
    write(unit=stmp,fmt='("pro_",I3.3,".dat")') k
    call SaveDataC8(trim(stmp), MS_propagator(1:nx,1:ny,k),nx*nx,err)
  end if
  end do ! L- all propagators
! ------------

  goto 1000

! ------------
! error handling
101 call MS_ERROR("Module not initialised.",subnum+1)
  goto 1000
102 call MS_ERROR("No slices specified.", subnum+2)
  goto 1000
103 call MS_ERROR("No propagators, zero object thickness.",subnum+3)
  goto 1000
104 call MS_ERROR("MultiSlice propagator allocation failed.",subnum+4)
  goto 1000

! ------------
1000 continue
!  write(*,*) " > MS_PreparePropagators2: EXIT."
  return

END SUBROUTINE MS_PreparePropagators2
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_SetStack(nstack, nstacklen)
! function: set content of the multislice stack list
! -------------------------------------------------------------------- !
! parameter: nstack : integer*4, dimension(nstacklen) : stack list
!            nstacklen : integer*4 : length of stack
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1000
  integer*4, dimension(nstacklen) :: nstack
  integer*4 :: nstacklen
  integer*4 :: i, n
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_SetStack: INIT."
  if (MS_status<1) then
    call MS_ERROR("Module not initialised.",subnum+1)
    return
  end if
  MS_slicestack(:) = 0
! ------------

! ------------
  n = min(MS_stacksize,nstacklen)
  do i=1,n
    MS_slicestack(i) = nstack(i)
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > MS_SetStack: EXIT."
  return

END SUBROUTINE MS_SetStack
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_FFT(cdata,nx,ny,dir)
! function: calculates a complex-to-complex 2D Fourier transform of
!           cdata with given dimensions nx, ny in directions dir
!           dir >=0: forwards, <0: backwards
! remarks : The transformation never applyies a re-normalization of
!           the data. With each call, there is a factor of
!           sqrt(nx*ny) applied to each coefficient.
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 3400
  integer*4, intent(in) :: nx, ny, dir
  complex*8, intent(inout) :: cdata(nx,ny)   
! fft variables
  type(C_PTR) :: fft_plan ! local plan
  type(C_PTR) :: fft_pdata ! local pointers to aligned data
  complex(C_FLOAT_COMPLEX), pointer :: work(:,:) ! pointer to working array
! ------------

! ------------
  if (nx <= 0 .or. ny <= 0) return ! skip invalid dimensions
  fft_pdata = fftwf_alloc_complex( int(ny*nx, C_SIZE_T) )
  call c_f_pointer(fft_pdata, work, [nx,ny])
  if (dir >= 0) then ! forward plan
    fft_plan = fftwf_plan_dft_2d(ny, nx, work, work, FFTW_FORWARD, FFTW_ESTIMATE)
  else
    fft_plan = fftwf_plan_dft_2d(ny, nx, work, work, FFTW_BACKWARD, FFTW_ESTIMATE)
  end if
  work(1:nx,1:ny) = cdata(1:nx,1:ny)
! ------------

! ------------
  call fftwf_execute_dft(fft_plan, work, work)
! ------------

! ------------
  cdata(1:nx,1:ny) = work(1:nx,1:ny)
! clean up fft
  call fftwf_free(fft_pdata)
  call fftwf_destroy_plan(fft_plan)
  
! ------------
  return

END SUBROUTINE MS_FFT
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_FFT_WORK(dir)
! function: Calculates a DFT of MS_work with direction dir
! -------------------------------------------------------------------- !
! parameter: integer*4 :: dir :: transform direction
!                                >=0 : forwards
!                                < 0 : backwards
! -------------------------------------------------------------------- !

  implicit none

  integer*4, parameter :: subnum = 3450
  integer*4, intent(in) :: dir

  if (MS_fft_dims(1) > 0 .and. MS_fft_dims(2) > 0) then
    if (dir>=0) then ! forwards fft
      call fftwf_execute_dft(MS_fft_work_planf, MS_work, MS_work)
    else ! backwards fft
      call fftwf_execute_dft(MS_fft_work_planb, MS_work, MS_work)
    end if
  end if
  
  
  return

END SUBROUTINE MS_FFT_WORK
!**********************************************************************!














!*********************************************************************!
!*********************************************************************!
!    CALCULATIONS                                                     !
!*********************************************************************!
!*********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_Start(istart)
! function: Starts Multislice algorithm
! -------------------------------------------------------------------- !
! parameter: integer*4 :: calculation start index
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1100
  complex*8 :: c0
  integer*4, optional, intent(in) :: istart
  integer*4 :: nslice, ncurslice, nstart, nstop
  character(len=1024) :: stmp
  external :: ExportWave
  external :: sinsertslcidx ! (idx,idxlen,sfnin,sfnadd,sfnext,sfnout) ! msasub.f90
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_Start: INIT."
  if (MS_status<1) then
    call MS_ERROR("Module not initialised.",subnum+1)
    return
  end if
  if (MS_slicecur>=0) then
    call MS_ERROR("Algorithm still running.",subnum+3)
    return
  end if
  c0 = cmplx(0.0,0.0)
  nstart = 0
  nstop = 0
  if (present(istart)) then
    nstop = min(istart,MS_stacksize)
  end if
! ------------

! ------------
! setup starting wave
  MS_wave(:,:) = MS_wave_in(:,:)
!  MS_absorbed(:,:,:) = c0
!
! setup the multi-slice index and counters
! (default)
  MS_slicecur = 0
  MS_lastmaxslice = 0
  MS_calcthick = 0.0
! (inserted wave function) - version 0.60b - 14.06.2012 (JB)
  if (nstart<nstop) then
    ! simulate an empty calculation down to the insert position MSP_extinwslc
    do MS_slicecur=nstart, nstop
      ncurslice = MS_slicecur
      nslice = MS_slicestack(ncurslice+1)+1
      MS_calcthick = MS_calcthick + MS_slicethick(nslice)
    end do
    MS_slicecur = nstop
  end if
! ------------

! ------------
! optional wave export after each slice, also saves incoming wave
  if (MS_incwave_export==1) then
    ! prepare file name
    call sinsertslcidx(MS_slicecur,MS_nslid,trim(MS_wave_filenm),"",".wav",stmp)
    MS_wave_avg_idx = 0 ! store incoming wave in channel 0 regardless of its actual plane.
    MS_pint_idx = 0
    !
    call ExportWave(trim(stmp))
    !
  end if

! ------------
!  write(unit=*,fmt=*) " > MS_Start: EXIT."
  return

END SUBROUTINE MS_Start
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_Stop()
! function: Stops Multislice algorithm, leaves current MS_wave data as is
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1200
  integer*4 :: nslmax
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_Stop: INIT."
  if (MS_status<1) then
    call MS_ERROR("Module not initialised.",subnum+1)
    return
  end if
  nslmax = MS_lastmaxslice
! ------------

! ------------
! reset slice id
  MS_slicecur = -1
! ------------

! ------------
!  write(unit=*,fmt=*) " > MS_Stop: EXIT."
  return

END SUBROUTINE MS_Stop
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_CalculateNextSlice(slc, nx,ny)
! function: Calculates wave with next slice in stack
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1300
  integer*4, intent(in) :: nx, ny
  complex*8, intent(in) :: slc(nx,ny)
  
  integer*4 :: i,j, nslice, ncurslice, nprop, i1, j1
  character(len=MS_ll) :: stmp
  
  external :: ExportWave
  external :: sinsertslcidx ! (idx,idxlen,sfnin,sfnadd,sfnext,sfnout) ! msasub.f90
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_CalculateNextSlice: INIT."
  if (MS_status<1) then
    call MS_ERROR("Module not initialised.",subnum+1)
    return
  end if
  if (MS_slicecur<0) then
    call MS_ERROR("No multislice calculation is running.",subnum+3)
    return
  end if
! ------------

! ------------
! count up slice
  MS_slicecur = MS_slicecur + 1
  ncurslice = MS_slicecur
  if (ncurslice>MS_stacksize) then ! final step reached, stop algorithm automatically
    MS_slicecur = -1
    return
  end if
  MS_lastmaxslice = ncurslice
  nslice = MS_slicestack(ncurslice)+1
  nprop = MS_propstack(nslice)
! ------------

! ------------
! in the beginning, the wave is given in Fourier-space
!   scrambled and transposed
! ------------

! ------------
! transform to real space
  !call MS_FFT(MS_wave(:,:),MS_dimx,MS_dimy,"back") ! external call from SFFTs.f
  call fftwf_execute_dft(MS_fft_planb,MS_wave,MS_wave)
! ------------


! ------------
! apply the phase grating
! distinguish two cases:
! case 1: slice data size is equal or larger than wave size, no wrap around / repetition, this is faster
! case 2: slice data size is smaller thn wave size, repetition of data is required
!         the repetition is done by taking the iterator (modulo) slice size
  if (nx>=MS_dimx .and. ny>=MS_dimy) then ! case 1
  
    do j=1, MS_dimy
      do i=1, MS_dimx
        MS_wave(i,j) = MS_wave(i,j)*slc(i,j)
      end do
    end do
  
  else ! case 2 (one of the slice sizes is smaller than the wave frame
    
    do j=1, MS_dimy
      j1 = modulo(j-1,ny)+1
      do i=1, MS_dimx
        i1 = modulo(i-1,nx)+1
        MS_wave(i,j) = MS_wave(i,j)*slc(i1,j1)
      end do
    end do
    
  end if
! ------------
  

! ------------
! transform to fourier space
  !call MS_FFT(MS_wave(:,:),MS_dimx,MS_dimy,"for") ! external call from SFFTs.f
  call fftwf_execute_dft(MS_fft_planf,MS_wave,MS_wave)
! ------------



! ------------
! apply the fresnel propagator
  do j=1, MS_dimy
    do i=1, MS_dimx
      MS_wave(i,j) = MS_wave(i,j)*MS_propagator(i,j,nprop)
    end do
  end do
  MS_calcthick = MS_calcthick + MS_slicethick(nslice)
! ------------


! ------------
! in the end, the wave is given in Fourier-space
!   scrambled
! ------------
!
! optional wave export at selected slices, version 0.52b
! - moved here with special export routine ExportWave
! - the export routine does another inverse FT before saving
! - by this way the wave function is multiplied by the propagator at the end
! - this version is used for CTEM mode only, where we want the full thickness effect
  if (      (MS_wave_export==1 .or. MS_wave_avg_export==1 .or. &
             MS_pint_export==1 ) &
     & .and. 0==modulo(ncurslice,MS_wave_export_pzp) ) then
    ! prepare file name
    call sinsertslcidx(ncurslice,MS_nslid,trim(MS_wave_filenm),"",".wav",stmp)
    MS_wave_avg_idx = ncurslice/MS_wave_export_pzp ! this should always give an integer number
    MS_pint_idx = MS_wave_avg_idx
    !
    call ExportWave(trim(stmp))
    !
  end if
! ------------
  
! ------------
!  write(unit=*,fmt=*) " > MS_CalculateNextSlice: EXIT."
  return

END SUBROUTINE MS_CalculateNextSlice
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE MS_ApplyPSpatialCoh(rdata,nx,ny,sx,sy,srcrad,ntype)
! function: convolution with source function on rdata
!           by multiplication with source geometry spectrum in fourier-space
!           and apodization with a cos**2 half-period in the outer 1/4
! -------------------------------------------------------------------- !
! parameter: 
!       real*4 :: rdata(nx,ny) ! input and output image array
!       integer*4 :: nx, ny    ! input and output image dimensions
!       real*4 :: sx, sy       ! image sampling rate
!       real*4 :: srcrad       ! source radius parameter [nm]
!       integer*4 :: ntype     ! distribution type selector (optional)
!                              ! 0 = delta function = no convolution
!                              ! 1 = Gaussian
!                              ! 2 = Cauchy (Lorentzian)
!                              ! 3 = Disk (Sigmoidal)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1400
  real*4, intent(inout)             :: rdata(nx,ny)
  integer*4, intent(in)             :: nx, ny
  real*4, intent(in)                :: sx, sy, srcrad
  integer*4, intent(in), optional   :: ntype
  
  integer*4 :: i, j, ki, kj, i2, j2, ntyuse
  integer*4 :: nkdimx2, nkdimy2, nalloc
  
  real*4 :: kernprm, kernrad, kernsum, kernval
  real*4 :: pkx, pky, pkys, pks
  !real*4 :: itowx, itowy, pfac, wx, wy, w2x, w2, wmax, wl
  !real*4 :: apod, apodprm1, apodprm2
  real*4, allocatable :: rwork(:,:), rkernel(:,:)
  !complex*8 :: cval0, cval
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MS_ApplyPSpatialCoh: INIT."
! prepare parameters
!  ndimx = nx
!  ndim2x = ndimx/2
!  ndimy = ny
!  ndim2y = ndimy/2
  ntyuse = 1
! get optional type selector
  if (present(ntype)) then
    ntyuse = max(0,min(3,ntype))
  end if
!  write(*,*) "- using kernel type ", ntyuse
! return in case of delta distribution
  if (ntyuse==0) return  
! determine the size of the kernel to use
  select case (ntyuse)
    case (1)
      kernrad = 3.0*srcrad
    case (2)
      kernrad = 10.0*srcrad
    case (3)
      kernrad = 1.5*srcrad
  end select ! case (ntype)
  nkdimx2 = CEILING(kernrad/sx)
  nkdimy2 = CEILING(kernrad/sy)
! ------------
  

! ------------
! allocate memory for the working array and the kernel, zero both arrays
  allocate(rwork(nx,ny), rkernel(-nkdimx2:nkdimx2,-nkdimy2:nkdimy2), stat=nalloc)
  if (nalloc/=0) then
    call MS_ERROR("Failed to allocate memory for finite source convolution.", nalloc)
    return
  end if
  rwork = 0.0
  rkernel = 0.0
  kernsum = 0.0
! ------------


! ------------
! KERNEL PREPARATION
  select case (ntyuse)
  
    case (1) ! GAUSSIAN
      !write(*,*) "- Gaussian kernel setup:"
      kernprm = -1.0/(srcrad*srcrad)
!      write(*,*) "  - 1/e HW: ",srcrad," nm"
!      write(*,*) "  - cut-off: ",kernrad," nm"
!      write(*,*) "  - nkdimx2: ",nkdimx2
!      write(*,*) "  - nkdimy2: ",nkdimy2
      do j = -nkdimy2, nkdimy2
        pky = real(j)*sy
        pkys = pky*pky
        do i = -nkdimx2, nkdimx2
          pkx = real(i)*sx
          pks = pkx*pkx + pkys
          if (pks>kernrad*kernrad) cycle ! cut-off
          kernval = EXP( pks*kernprm )
          kernsum = kernsum + kernval
          rkernel(i,j) = kernval
!          write(*,*) "  - kernel at (",pkx,",",pky,")nm = ",kernval
        end do
      end do
      ! normalization
      rkernel = rkernel / kernsum
      
    case (2) ! CAUCHY (LORENTZIAN)
      kernprm = srcrad*srcrad
      do j = -nkdimy2, nkdimy2
        pky = real(j)*sy
        pkys = pky*pky
        do i = -nkdimx2, nkdimx2
          pkx = real(i)*sx
          pks = pkx*pkx + pkys
          if (pks>kernrad*kernrad) cycle ! cut-off
          kernval = 1.0 / ( pks + kernprm )
          kernsum = kernsum + kernval
          rkernel(i,j) = kernval
        end do
      end do
      ! normalization
      rkernel = rkernel / kernsum
    
    case (3) ! DISC (SIGMOID)
      kernprm = srcrad
      do j = -nkdimy2, nkdimy2
        pky = real(j)*sy
        pkys = pky*pky
        do i = -nkdimx2, nkdimx2
          pkx = real(i)*sx
          pks = sqrt(pkx*pkx + pkys)
          if (pks>kernrad) cycle ! cut-off
          kernval = 0.5-0.5*tanh((pks/kernprm - 1.0)*100.0)
          kernsum = kernsum + kernval
          rkernel(i,j) = kernval
        end do
      end do
      ! normalization
      rkernel = rkernel / kernsum
      
  end select ! case (ntype)
! ------------


! ------------
! CONVOLUTION with wrap around
  do j=1, ny ! loop over output rows
    do i=1, nx ! loop over output columns
      do kj = -nkdimy2, nkdimy2 ! loop over kernel rows
        ! get input row index
        j2 = modulo(j + kj - 1, ny) + 1 ! wrap to array bounds
        do ki = -nkdimx2, nkdimx2 ! loop over kernel columns
          ! get input coulmn index
          i2 = modulo(i + ki -1, nx) + 1 ! wrap to array bounds
          !
          ! convultion integral
          rwork(i,j) = rwork(i,j) + rkernel(ki,kj)*rdata(i2,j2)
          !
        end do
      end do
    end do
  end do
! ------------

  
! ------------
! transfer from rwork to rdat
  rdata = rwork
! ------------


! ------------
! deallocate working arrays
  deallocate(rwork, rkernel, stat=nalloc)
  if (nalloc/=0) then
    call MS_ERROR("Failed to allocate memory for finite source convolution.", nalloc)
    return
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > MS_ApplyPSpatialCoh: EXIT."
  return

END SUBROUTINE MS_ApplyPSpatialCoh
!**********************************************************************!





END MODULE MultiSlice

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
!  integer*4, parameter :: subnum = 2300
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