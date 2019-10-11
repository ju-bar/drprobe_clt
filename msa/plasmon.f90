!**********************************************************************!
!**********************************************************************!
!                                                                      !
!                    file   "plasmon".f90                              !
!                                                                      !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                  !
!    Version  :  1.1.0, Oct.   11, 2019                                !
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
!   Purpose: implementing subroutines for electron plasmon scattering  !
!            using a Monte-Carlo Approach                              !
!            B. Mendis, Ulramic. 206 (2019) 112816.                    !
!            doi:10.1016/j.ultramic.2019.112816                        !
!                                                                      !
!   Contains corrections and alternative code based on discussions     !
!   between J. Barthel, B. Mendis, S. Findlay, and  L.J. Allen         !
!   (Aug.-Oct. 2019)                                                   !
!                                                                      !
!   External sources to link:                                          !
!            random.f90 : provides random number generators            !
!              UniRand() : float, uniform distribution [0, 1[          !
!              PoissonRand(m): integer, Poisson distr. with mean m     !
!                                                                      !
!**********************************************************************!
  
!**********************************************************************!
!
! How to use this module ?
!
! 1) Set required input parameters
!    PL_ek = electron kinetic energy [eV]
!    PL_ep = energy of single plasmon excitations [eV]
!    PL_lp = mean-free path for single plasmon excitations [nm]
!    PL_tmax = max. sample thickness [nm]
!    PL_wthr = probability threshold for neglecting higher-excitation
!              plasmon excitations, keep between 0.001 and 0.1
!              smaller values open higher excitation levels
!    (then call PL_init)
!  
!    Alternatively, set
!    PL_lp = mean-free path for single plasmon excitations [nm]
!    PL_qe = characteristic angle for single scattering [rad]
!    PL_qc = critical angle for single scattering [rad]
!    PL_npemax = number of excitations allowed per electron
!    (then call PL_init2)
!
! 2) Call PL_init(errorcode)
!      before starting the Monte-Carlo
!      the returned errorcode indicates problems with the parameters
!      set in step (1), errorcode == 0 means success
!
!    Alternative, for the alternative parameter setup, see (1)
!    Call PL_init2(errorcode)
!      before starting the Monte-Carlo
!      the returned errorcode indicates problems with the parameters
!      set in step (1), errorcode == 0 means success

!
! 3) Call PL_reset()
!      before each individual Monte-Carlo pass
!
! 4) Call PL_scatt_slc(dt, sqx, sqy, wave, nx, ny)
!      for each slice of the object
!      input:
!        dt = slice thickness [nm]
!        sqx, sqy = wave function sampling rate [rad/pixel]
!        wave(:,:) = wave function samples
!        nx, ny = number of wave function samples
!      output:
!        wave(:,:) = possibly modified wave function samples
!
! 5) Call PL_populate(num_exc)
!      after each individual Monte-Carlo pass to register excitation
!      channel populations
!      input:
!        num_exc = number of plasmon-excitations in previous run
!
! 6) Call PL_deinit()
!      after all is done and the module is no longer needed.
!
!**********************************************************************!


!*********************************************************************!
!*********************************************************************!
!    MODULE DECALRATIONS                                              !
!*********************************************************************!
!*********************************************************************!

MODULE Plasmon

! Use PL_ as prefix for public parameters funcs and subs!
    
  ! Global module dependencies

  implicit none
  
  public :: PL_init ! call to initialize the module from input
                    ! see variable decl. below for required input
                    ! do this at least once before running the whole
                    ! Monte-Carlo calculation
  public :: PL_init2 ! call to initialize the module from alternative input
                     ! see variable decl. below for required input
                     ! do this at least once before running the whole
                     ! Monte-Carlo calculation
                     !
  public :: PL_deinit ! call to deallocate module arrays, once the whole
                      ! calculation is finished
                      !
  public :: PL_mc_slc ! call to decide for plasmon scattering in a thin slice.
                      ! /!\ Use this to decide scattering on-the-fly.
                      ! /!\ This is ment to be called by module functions only.
                      !
  public :: PL_reset ! call to reset the variables householding scattering per run
                     ! This should be done before each single Monte-Carlo pass
                     !
  public :: PL_populate ! call to update the plasmon-loss channel population statistics
                        ! This should be done after each single Monte-Carlo pass
                        !
  public :: PL_scatt_slc ! call to apply plasmon scattering in a slice
                         ! /!\ this changes the wave function depending on what
                         !     is returned by PL_mc_slc.
  
  ! 
  ! physical parameters - to be defined before calling any routine
  !
  real*4, public, parameter :: PL_e0 = 510999. ! electron rest energy [eV]
  !
  ! /!\ REQUIRED INPUT PARAMETERS for PL_INIT
  ! /!\ ALTERNATIVE INPUT PARAMETERS FOR PL_INIT2 ARE MARKED BY -> (*AI*)
  !
  real*4, public :: PL_ek ! kinetic energy of the probing electron [eV]
  DATA PL_ek /200000./
  real*4, public :: PL_ep ! plasmon energy [eV]
  DATA PL_ep /20./
  real*4, public :: PL_lp ! inelastic mean-free path [nm] (*AI*)
  DATA PL_lp /100./
  real*4, public :: PL_tmax ! max. sample thickness [nm]
  DATA PL_tmax /0./
  real*4, public :: PL_wthr ! probability threshold for neglecting higher
                            ! plasmon excitations
  DATA PL_wthr /0.01/       ! default: 1%
  !
  ! /!\ END OF INPUT PARAMETERS
  !

  ! derived parameters - calculated by PL_init from the above values
  ! /!\ Modify these variables with care, as they determine how the
  !     plasmon excitation Monte-Carlo performs.
  real*4, public :: PL_tol ! t / Lambda
  real*4, public :: PL_qe, PL_qe2 ! characteristic angle [rad] (*AI*)
  real*4, public :: PL_qc, PL_qc2 ! critical angle [rad] (*AI*)
  integer*4, public :: PL_npemax ! max number of plasmon excitations (*AI*)
  DATA PL_npemax /0/ ! turn this off by default
  
  ! results of the Monte-Carlo for one run
  ! /!\ Allocations sizes on last dimension: (0:PL_npemax) /!\
  !     corresponds to the number of plasmon excitations.
  integer*4, public :: PL_exc_num ! number of excitations for the run <= PL_npemax
  DATA PL_exc_num /0/
  real*4, public, allocatable :: PL_exc_dq(:,:) ! excitation scattering angle [rad]
  real*4, public, allocatable :: PL_exc_dqtot(:,:) ! accumulated scattering angle [rad]
  
  ! statistics of the whole Monte-Carlo
  integer*4, public :: PL_mc_num ! number of MC runs
  DATA PL_mc_num /0/
  integer*4, public, allocatable :: PL_mc_exc_pop(:) ! population numbers
  
  ! infrastructure
  integer*4, public :: PL_num_err ! number of errors
  DATA PL_num_err /0/ 
  character(len=2048), public :: PL_msg_err ! last error message
  DATA PL_msg_err /""/
  
  
  
!*********************************************************************!
!*********************************************************************!
!    MODULE BODY                                                      !
!*********************************************************************!
!*********************************************************************!

  CONTAINS
  
  !*******************************************************************!
  !
  ! PL_init
  !
  ! input: none
  !
  ! output: 
  !          integer*4 :: nerr = error code
  !
  ! purpose: (re-)initialize module variables based on current input
  !          parameters PL_ek, PL_ep, PL_lp, PL_tmax, PL_nslc, PL_wthr
  !
  !          This is the initialization assuming real bulk plasmon
  !          scattering as determined by plasmon energy and mean free
  !          path.
  !
  subroutine PL_init(nerr)
    
    implicit none
    
    integer*4, intent(inout) :: nerr
    integer*4 :: nalloc
    real*4 :: wt, facn
    
    nerr = 0
    nalloc = 0
    call PL_deinit() ! reset
    
    if (PL_ep<=1.) goto 102
    if (PL_ek<PL_ep) goto 103
    if (PL_wthr<0. .or. PL_wthr>1.) goto 104
    if (PL_lp<=0.) goto 105
    if (PL_tmax<0.) goto 106
    
    ! calculate derived parameters
    ! - characteristic angle [rad] (Edgerton (2011) Electron Energy-Loss Spectroscopy in the Electron Microscope)
    PL_qe = PL_ep / PL_ek * (PL_ek + PL_e0) / (PL_ek + 2.*PL_e0)
    PL_qe2 = PL_qe * PL_qe
    ! - critical angle [rad] ! This parameter is critical.
    !   Egerton Ultramic. 107 (2007) 575-586.: approximation
    PL_qc = sqrt(PL_ep / PL_ek) 
    PL_qc2 = PL_qc * PL_qc
    ! - t over lambda
    PL_tol = PL_tmax/PL_lp
    ! max number number of registered plasmon excitations
    PL_npemax = 0
    wt = 1. - exp(-PL_tol)
    facn = 1.
    do while (wt > PL_wthr) ! include more excitation levels until remaining
	                        ! total probability of all higher levels is below
							! the user defined threshold PL_wthr
      PL_npemax = PL_npemax + 1
      facn = facn * real(PL_npemax)
      wt = wt - (PL_tol**PL_npemax * exp(-PL_tol) / facn)
    end do
	!
	! allocations
    allocate(PL_exc_dq(1:2,0:PL_npemax), stat=nalloc)
    if (nalloc/=0) goto 107
    allocate(PL_exc_dqtot(1:2,0:PL_npemax), stat=nalloc)
    if (nalloc/=0) goto 107
    allocate(PL_mc_exc_pop(0:PL_npemax), stat=nalloc)
    if (nalloc/=0) goto 107
    PL_exc_dq = 0.0
    PL_exc_dqtot = 0.0
    PL_mc_exc_pop = 0
    
100 return
    
102 nerr = 2
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_ep<1) too small plasmon energy"
    goto 100
103 nerr = 3
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_ek<PL_ep) too small electron energy"
    goto 100
104 nerr = 4
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): PL_wthr out of bounds (>0 & <1)"
    goto 100
105 nerr = 5
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_lp<=0) mean-free path not positive"
    goto 100
106 nerr = 6
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_tmax<0) negative sample thickness"
    goto 100
107 nerr = 7
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): memory allocation failed"
    goto 100
  
  end subroutine PL_init
  
  
  !*******************************************************************!
  !
  ! PL_init2
  !
  ! input: none
  !
  ! output: 
  !          integer*4 :: nerr = error code
  !
  ! purpose: (re-)initialize module variables based on current input
  !          parameters PL_lp, PL_qe, PL_qc, PL_npemax
  !
  !          This is the initialization for faking inelastic scattering
  !          for low-loss intraband transitions for a given mean-free
  !          path, characteristic angle and maximum number of
  !          excitations per incident electron.
  !
  subroutine PL_init2(nerr)
    
    implicit none
    
    integer*4, intent(inout) :: nerr
    integer*4 :: nalloc
    
    nerr = 0
    nalloc = 0
    ! reset (alternative)
    PL_num_err = 0
    PL_msg_err = ""
    PL_exc_num = 0
    PL_mc_num = 0
    if (allocated(PL_exc_dq)) deallocate(PL_exc_dq, stat=nalloc)
    if (allocated(PL_exc_dqtot)) deallocate(PL_exc_dqtot, stat=nalloc)
    if (allocated(PL_mc_exc_pop)) deallocate(PL_mc_exc_pop, stat=nalloc)
    !
    PL_ep = 10. ! set some fake value (this will not be used)
    if (PL_qe<=0.) goto 102
    if (PL_qc<=PL_qe) goto 103
    if (PL_npemax<=0) goto 104
    if (PL_lp<=0.) goto 105
    !
    ! calculate derived parameters
    PL_qe2 = PL_qe * PL_qe
    PL_qc2 = PL_qc * PL_qc
    PL_tol = PL_tmax / PL_lp ! just in case someone defined tmax, not required
	  !
	  ! allocations
    allocate(PL_exc_dq(1:2,0:PL_npemax), stat=nalloc)
    if (nalloc/=0) goto 107
    allocate(PL_exc_dqtot(1:2,0:PL_npemax), stat=nalloc)
    if (nalloc/=0) goto 107
    allocate(PL_mc_exc_pop(0:PL_npemax), stat=nalloc)
    if (nalloc/=0) goto 107
    PL_exc_dq = 0.0
    PL_exc_dqtot = 0.0
    PL_mc_exc_pop = 0
    
100 return
    
102 nerr = 2
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_qe<=0) too small characteristic angle"
    goto 100
103 nerr = 3
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_qc<PL_qe) too small critical angle"
    goto 100
104 nerr = 4
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init): (PL_npemax<=0) too low number of allowed excitations"
    goto 100
105 nerr = 5
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init):(PL_lp<=0) mean-free path not positive"
    goto 100
107 nerr = 7
    PL_num_err = PL_num_err + 1
    PL_msg_err = trim(PL_msg_err)//" (init):memory allocation failed"
    goto 100
  
  end subroutine PL_init2
  
  
  
  !*******************************************************************!
  !
  ! PL_deinit
  !
  ! input: none
  ! output: none
  ! purpose: de-initialize module variables
  !
  subroutine PL_deinit()
    
    implicit none
    
    integer*4 :: nalloc
    
    nalloc = 0
    PL_num_err = 0
    PL_msg_err = ""
    PL_npemax = 0
    PL_exc_num = 0
    PL_mc_num = 0
    
    !if (allocated(PL_exc_islc)) deallocate(PL_exc_islc, stat=nalloc)
    if (allocated(PL_exc_dq)) deallocate(PL_exc_dq, stat=nalloc)
    if (allocated(PL_exc_dqtot)) deallocate(PL_exc_dqtot, stat=nalloc)
    if (allocated(PL_mc_exc_pop)) deallocate(PL_mc_exc_pop, stat=nalloc)
    
    return
  
  end subroutine PL_deinit
  
  
  
  !*******************************************************************!
  !
  ! PL_mc_slc
  !
  ! input:
  !     real*4 :: dt : slice thickness
  !     integer*4 :: iexc_max : max. allowed excitations
  !
  ! output:
  !     integer*4 :: iexc : excitation flag (0 = no, >0: number of excitations)
  !     real*4 :: output dqx, dqy scattering angle [rad]
  !
  ! purpose: calculate whether inelastic scattering takes place
  !          in a slice, and if so, determines the scattering angle.
  !
  ! results: None of the other module variables are changed. These
  !          changes need to be applied by the calling routine.
  !
  subroutine PL_mc_slc(dt, iexc_max, iexc, dqx, dqy)
  
    implicit none
    
    real*4, intent(in) :: dt ! slice thickness
    integer*4, intent(in) :: iexc_max ! max. allowed excitations
    integer*4, intent(out) :: iexc ! excitation flag
    real*4, intent(out) :: dqx, dqy ! scattering angle [rad]
    
    integer*4 :: i
    real*4 :: pcur, qcur
    
    real*4, external :: UniRand ! from "random.f90" -> RNG: 0 ... 1
    integer*4, external :: PoissonRand ! from "random.f90" -> Poissonian RNG
    
    iexc = 0
    dqx = 0.
    dqy = 0.
    
    ! calculate a poissonian random number (limited to max. allowed excitation level)
    iexc = min(PoissonRand(dt/PL_lp), iexc_max)
    
    if (iexc > 0) then ! excitation happens
      do i=1, iexc
        pcur = 6.283185 * UniRand() ! uniform RNG phi
        qcur = sqrt( PL_qe2 * (PL_qc2 / PL_qe2 + 1.)**UniRand() - PL_qe2) ! Lorentzian Theta
        dqx = dqx + qcur * cos(pcur)
        dqy = dqy + qcur * sin(pcur)
      end do
    end if
    
100 return
  
  end subroutine PL_mc_slc
  
  
  !*******************************************************************!
  !
  ! PL_reset
  !
  ! input:
  !     none
  !
  ! output:
  !     none
  !
  ! purpose: resets the housholder arrays keeping track of tilts and
  !          number of excitations per run through the specimen
  !
  subroutine PL_reset()
  
    implicit none
    
    PL_exc_num = 0
    PL_exc_dq = 0.
    PL_exc_dqtot = 0.
    
  end subroutine PL_reset
  
  
  !*******************************************************************!
  !
  ! PL_populate
  !
  ! input:
  !     integer*4 :: num_exc ! number of plasmon excitations to register
  !                          ! for the finished run in the population
  !                          ! statistics
  !
  ! output:
  !     none
  !
  ! purpose: updates the population statistics in plasmon-loss channels
  !          call this when a multislice run is finished
  !
  subroutine PL_populate(num_exc)
  
    implicit none
    
    integer*4, intent(in) :: num_exc
    integer*4 :: n
    
    n = max(0, min( PL_npemax, num_exc) ) ! plasmon-loss channel 0 .. PL_npemax
    PL_mc_exc_pop(n) = PL_mc_exc_pop(n) + 1 ! increment excitation level count
	  PL_mc_num = PL_mc_num + 1 ! increment number of MC passes
	
    
  end subroutine PL_populate
  
  
  !*******************************************************************!
  !
  ! PL_scatt_slc
  !
  ! input:
  !     real*4 :: dt ! current slice thickness [nm]
  !     real*4 :: sqx, sqy ! sampling rates of the wave function [rad/pixel]
  !     complex*8 :: wave(:,:) ! wave function before plasmon excitat.
  !     integer*4 :: nx, ny ! size of the wave function grid
  !
  ! output:
  !     complex*8 :: wave(:,:) ! wave function after plasmon excitat.
  !
  ! purpose: resets the housholder arrays keeping track of tilts and
  !          number of excitations per run through the specimen
  !
  subroutine PL_scatt_slc(dt, sqx, sqy, wave, nx, ny)
  
    implicit none
    
    real*4, intent(in) :: dt ! slice thickness [nm]
    real*4, intent(in) :: sqx, sqy ! sampling rates of the wave function [rad/pixel]
    complex*8, intent(inout) :: wave(nx,ny) ! wave function
    integer*4, intent(in) :: nx, ny ! grid size
    
    integer*4 :: iexc, iexc_max, nexc ! number of excitations
    integer*4 :: isx0, isy0, isx1, isy1, isxd, isyd ! pixel shifts
    integer*4 :: i, j, i1, j1
    real*4 :: dqx, dqy
    complex*8, allocatable :: wtmp(:,:)
    
    iexc = 0
    iexc_max = PL_npemax - PL_exc_num ! max. further allowed excitations
    dqx = 0.
    dqy = 0.
    
    call PL_mc_slc(dt, iexc_max, iexc, dqx, dqy)
    
    if (iexc > 0) then ! excitations have happened
      
      ! new excitation level
      nexc = PL_exc_num + iexc
      
      ! set the change of scattering angle with this excitation
      PL_exc_dq(1, nexc ) = dqx
      PL_exc_dq(2, nexc ) = dqy
      
      ! set the total scattering angle with this excitation
      do i=PL_exc_num+1, nexc ! loop in case that more than 1 excitation happened (iexc>=1)
        PL_exc_dqtot(1, i) = PL_exc_dqtot(1, PL_exc_num ) + dqx
        PL_exc_dqtot(2, i) = PL_exc_dqtot(2, PL_exc_num ) + dqy
      end do
      
      ! get current pixel shift ( dqtot old -> pixels )
      isx0 = nint( PL_exc_dqtot(1, PL_exc_num) / sqx )
      isy0 = nint( PL_exc_dqtot(2, PL_exc_num) / sqy )
      ! get new pixel shift ( dqtot new -> pixels )
      isx1 = nint( PL_exc_dqtot(1, nexc) / sqx )
      isy1 = nint( PL_exc_dqtot(2, nexc) / sqy )
      ! get change of pixel shift (pixel new - pixel old)
      isxd = isx1 - isx0
      isyd = isy1 - isy0
      
      if (isxd/=0 .or. isyd/=0) then
		! apply wave function tilt by shifting in Fourier space
		! wave is expected to be in this representation
        ! write(*,*) "PL_scatt_slc: tilt shift: (", isxd, ",", isyd, ")"
		!
		! The implementation below looks somewhat inefficient,
		! however it allows to interface "wave" by C pointers
		! and allocatable arrays without access violations.
		!
		! allocate a dummy buffer for the wave function
        allocate(wtmp(nx,ny), stat=i)
		! copy current w.f. to the dummy
        do j=1, ny
          do i=1, nx
            wtmp(i,j) = wave(i,j)
          end do
        end do
		! copy back from the dummy by re-addressing with shift
        do j=1, ny
          j1 = 1 + modulo(j-isyd, ny)
          do i=1, nx
            i1 = 1 + modulo(i-isxd, nx)
            wave(i,j) = wtmp(i1,j1)
          end do
        end do
		! de-alloc the dummy
        deallocate(wtmp, stat=i)
      end if
      
      ! update the number of excitations for this run
      PL_exc_num = PL_exc_num + iexc
      
    end if
    
    return
    
  end subroutine PL_scatt_slc
  
  
END MODULE Plasmon

!*********************************************************************!
!*********************************************************************!
!    MODULE END                                                       !
!*********************************************************************!
!*********************************************************************!

