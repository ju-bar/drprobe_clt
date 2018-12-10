!**********************************************************************!
!
! FILE "celslcprm.f90"
! FORTRAN MODULE celslcprm
! 
! PURPOSE: global parameters for program CELSLC
! 
! AUTHOR: Juri Barthel, ju.barthel@fz-juelich.de
!         RWTH Aachen University, Aachen, Germany
!         21.09.2015
!
! last modified: J.B. 29.05.2018
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


MODULE celslcprm


! ---[ routine names ]------------------------------------------------------------------------------

  public :: csprm_initcl
  public :: csprm_getcls

! ---[ global parameters ]--------------------------------------------------------------------------

  integer*4, parameter :: stdout = 6        ! standard output unit
  integer*4, parameter :: fft_dmin = 32     ! min. size of fft
  integer*4, parameter :: fft_dmax = 8192   ! max. size of fft

! ---[ global variables ]---------------------------------------------------------------------------

  integer*4 :: nerr                         ! error code
  integer*4 :: nverbose                     ! verbosity level of std out
  real*4 :: wl                              ! wavelength
  real*4 :: wl1                             ! wavelength corrected for refraction in the sample
  real*4 :: u0cell                          ! mean inner potential of the input cell
  real*4 :: sdx, sdy, sdz                   ! and cell sampling
  
  ! defaults
  DATA nerr /0/
  DATA nverbose /1/
  DATA wl   /0.001969/
  DATA wl1  /0.001969/
  DATA u0cell /0.0/
  DATA sdx  /0.005/
  DATA sdy  /0.005/
  DATA sdz  /0.1/
  
  ! run time measurement
  integer*4 :: csprm_runtimes               ! switch for run time measurements
  integer*8 :: csprm_clinit                 ! initialized clock state
  integer*8 :: csprm_clrate                 ! system_clock rate (ticks per second)
  integer*8 :: csprm_clmax                  ! system_clock maximum tick count
  integer*8 :: csprm_clmeasure              ! last measured system_clock tick
  
  DATA csprm_runtimes /0/                   ! runtime measurement is off by default
  DATA csprm_clinit /0/                     ! initial clock state
  DATA csprm_clrate /1/                     ! initial clock rate
  DATA csprm_clmax /0/                      ! initial clock max. count
  DATA csprm_clmeasure /0/                  ! initial clock measure
      
! ---[ input variables ]----------------------------------------------------------------------------

  character(len=1024) :: sinputprm          ! input parameter file name
  character(len=1024) :: scellfile          ! super-cell file name
  character(len=1024) :: sfxfile            ! file name for x-ray scattering factors
  character(len=1024) :: sfefile            ! file name for el. scattering factor parameters
  integer*4 :: nfin                         ! input format: 0 = cel, 1 = cif, 10 = asc
  real*4 :: ht, abf, buniv                  ! el. energy, abs., universal Biso
  integer*4 :: nx, ny, nz, nv               ! super-cell discretization
  integer*4 :: nfl, ndwf, nabs, nabf        ! flags frozen lattice, Debye-Waller factor and abs.
  integer*4 :: nrev                         ! flag for slicing the supercell in reverse order
  integer*4 :: npot                         ! flag for potentials export
  integer*4 :: npps                         ! flag for potential output to slice files
  integer*4 :: n3dp                         ! flag for 3d potential data from input structure
  integer*4 :: nfx                          ! flag for using external x-ray scattering factors
  integer*4 :: nfe                          ! flag for using external el. scattering factors
  integer*4 :: nsca                         ! switch for using specific implemented form factors
                                            ! 0 = default -> 1,
                                            ! 1 = Weickenmeier and Kohl (EMS)
                                            ! 2 = Waasmaier and Kirfel (muSTEM)
  integer*4 :: ssc                          ! index of slice (0 = all, >0 only this slice)
  integer*4 :: buni                         ! flag for using a universal DWF
  integer*4 :: block                        ! flag for cutting a new block from the input structure
  integer*4 :: ntla                         ! flag for translation of atoms
  integer*4 :: nffdec                       ! flag for form factor decay output
  integer*4 :: nf2dec                       ! flag for form factor loss output
  real*4 :: bloh, blok, blol                ! block orientation in hkl of the input
  real*4 :: blyh, blyk, blyl                ! block y-axis direction in hkl of the input
  real*4 :: blsa, blsb, blsc                ! block dimension in nm
  real*4 :: tlax, tlay, tlaz                ! atom translation vector (fractional coordinates)
  real*4 :: vffdec                          ! form-factor decay value (for output only)
  real*4 :: vf2dec                          ! form-factor loss value (for output only)
  
  
  
  ! defaults
  DATA nfin /0/
  DATA ht /300.0/
  DATA abf /0.1/
  DATA buniv /0.005/
  DATA nx /64/
  DATA ny /64/
  DATA nz /0/
  DATA nv /1/
  DATA nfl /0/
  DATA ndwf /0/
  DATA nabs /0/
  DATA nabf /0/
  DATA nrev /0/
  DATA npot /0/
  DATA npps /0/
  DATA n3dp /0/
  DATA nfx /0/
  DATA nfe /0/
  DATA nsca /0/
  DATA ssc /0/
  DATA buni /0/
  DATA block /0/
  DATA bloh /0.0/
  DATA blok /0.0/
  DATA blol /0.0/
  DATA blyh /0.0/
  DATA blyk /0.0/
  DATA blyl /0.0/
  DATA blsa /0.0/
  DATA blsb /0.0/
  DATA blsc /0.0/
  DATA nffdec /0/
  DATA vffdec /0.04321392/ ! default decay of form factors exp(-pi)
  DATA nf2dec /0/
  DATA vf2dec /0.001/ ! default tolerated loss of scattering power

! ---[ output variables ]---------------------------------------------------------------------------

  integer*4 :: ndigsl, ndigvr               ! number of digits for slice files and variant files
  character(len=256) :: sslcfile            ! slice file name
  
  ! defaults
  DATA ndigsl /3/
  DATA ndigvr /3/
  
  
CONTAINS



SUBROUTINE csprm_initcl()
! function: Initializes the clock timer counts.
! parameter: none

  implicit none

  if (csprm_runtimes==1) then ! init the run-time clock state
    call system_clock(csprm_clinit, csprm_clrate, csprm_clmax)
  end if

  return

END SUBROUTINE csprm_initcl
!**********************************************************************!



SUBROUTINE csprm_getcls(cls)
! function: Determines the time span to the previous clock init in
!           seconds.
! parameter: cls : real*4 : time span in seconds

  implicit none

  integer*8 :: clnow, cldelta
  real*4, intent(out) :: cls

  cls = 0.0

  if (csprm_runtimes==1) then 
    call system_clock(clnow)
    if (clnow<csprm_clinit) then
      cldelta = csprm_clmax - csprm_clinit + clnow
      csprm_clinit = clnow
    else
      cldelta = clnow - csprm_clinit
    end if
    cls = real(cldelta,kind=4)/csprm_clrate
  end if

  return

END SUBROUTINE csprm_getcls



END MODULE celslcprm