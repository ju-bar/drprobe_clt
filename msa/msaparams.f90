!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  msaparams.F90                                         !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2009-2018   !
!                                                                      !
!**********************************************************************!
!                                                                      !
!    MODULE MSAparams                                                  !
!    -----------------                                                 !
!                                                                      !
!    Purpose  : parameters, parameter I/O and memory management for    !
!               the program MSA (see msa.f90)                          !
!    Version  : 1.2.3, Dec 05, 2018                                    !
!    To Link  : MultiSlice.F90                                         !
!               STEMfunctions.F90                                      !
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

!**********************************************************************!
!                                                                      !
!   Purpose: MODULE holding global parameters program msa              !
!                                                                      !
!**********************************************************************!
!                                                                      !
!   CONTAINS:                                                          !
!      1) SETUP & INIT ROUTINES                                        !
!      2) DATA MANAGEMENT ROUTINES                                     !
!      3) CALCULATIONS                                                 !
!      4) CALLING & INTERFACES                                         !
!                                                                      !
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
!**********************************************************************!












!**********************************************************************!
!*********************  MODULE MAIN  **********************************!
!**********************************************************************!
MODULE MSAparams

  use STEMfunctions
  use MultiSlice

  IMPLICIT NONE

! declare internal data types

! accessibility of subroutines or functions
!  private :: 
  private :: MSP_ERROR, MSP_WARN
  
!  public :: 
  public :: MSP_HALT
  public :: MSP_INIT, MSP_UNINIT
  public :: MSP_INITCL
  public :: MSP_GETCLS
  public :: MSP_READparams
  public :: MSP_READdetdef
  public :: MSP_READdetprofiles
  public :: MSP_READBLOCK_microscope
  public :: MSP_READBLOCK_multislice
  public :: MSP_ALLOCPGR, MSP_PREALLOCPGR, MSP_ALLOCDET
  public :: MSP_GetRandomProbeDefocus
  public :: MSP_GetRandomProbeShift
  public :: MSP_SetAnnularDetectors
  public :: MSP_SetRectangularDetectors
!  public :: MSP_SetCalibratedDetectors
  public :: MSP_GetPGRIndex
  public :: MSP_GetNumberOfDigits
  public :: MSP_InitTextOutput
  public :: MSP_WriteTextOutput
  public :: MSP_FinishTextOutput
  
  
  external :: PostMessage
  
! declare module global (but PRIVATE) variables, params and arrays

!   length of names and internal strings
  integer*4, public, parameter :: MSP_ll = 1024

!   enpty line placeholder (emptied at startup)
  character(len=MSP_ll), private :: MSP_el

!   number of known systems
  integer*4, private, parameter :: MSP_sys_max = 2
  
!   system identifiers
  integer*4, private, parameter :: MSP_sys_WIN  = 1
  integer*4, private, parameter :: MSP_sys_UNIX = 2
  
 
! declare module global (but PUBLIC) variables, params and arrays

!   type of compile version
  integer*4, public, parameter :: MSP_VERSION_NORMAL = 0
  integer*4, public, parameter :: MSP_VERSION_SERVICE = 99

!   standard output unit
  integer*4, public, parameter :: MSP_stdout = 6

!   current compile version
  integer*4, public :: MSP_compile_version

!   system identifier index
  integer*4, public :: MSP_sys
!   system identifier strings  
  character(len=MSP_ll), public :: MSP_sys_names(MSP_sys_max)
!   system directory delimiters
  character(len=MSP_sys_max), public :: MSP_sys_pathdelim
!   system count byte size
  integer*4, public :: MSP_sys_nszcnt(MSP_sys_max)

  
!   min data division value
  real*4, public, parameter :: MSP_div0 = 1.0E-15
  
!   run timing flag and variables
  integer*4, public :: MSP_runtimes
  DATA MSP_runtimes /0/ ! flags run time information  0: OFF (default) 1: ON (/rti)
  integer*8, public :: MSP_clinit
  DATA MSP_clinit /0/
  integer*8, public :: MSP_clrate
  DATA MSP_clrate /1/
  integer*8, public :: MSP_clmax
  DATA MSP_clmax /1/
  integer*8, public :: MSP_clmeasure
  DATA MSP_clmeasure /0/


 
! error and warning counters 
  integer*4, public :: MSP_err_num, MSP_warn_num
  
!   application info
  character(len=MSP_ll), public :: MSP_callApp
!   application version
  character(len=MSP_ll), public :: MSP_verApp
!   application author
  character(len=MSP_ll), public :: MSP_authApp
  
! flag for using the text output form (non-default)
  integer*4, public :: MSP_txtout
  DATA MSP_txtout /0/ ! 0: OFF (default), 1: ON
  
! flag for using 3d result files (no _sl suffixes added)
  integer*4, public :: MSP_3dout
  DATA MSP_3dout /0/ ! 0: OFF (default) = one file per export slice plane, 1: ON = all export planes in one file

! spatial coherence convolution switch
  integer*4, public :: MSP_ApplyPSC
  integer*4, public :: MSP_ExplicitPSC
  DATA MSP_ApplyPSC /0/
  DATA MSP_ExplicitPSC /0/

! multislice scan parameters
  integer*4, public :: MSP_ScanPixelX
  integer*4, public :: MSP_ScanPixelY
  
  integer*4, public :: MSP_LastScanPixelX
  integer*4, public :: MSP_LastScanPixelY
  
! incoming beam tilt [mrad]
  real*4, public :: MSP_BeamTiltX
  real*4, public :: MSP_BeamTiltY
  
! absorption parameters (only in case of potential input from slice files)
  integer*4, public :: MSP_nabf
  real*4, public :: MSP_Absorption
  DATA MSP_nabf /0/
  DATA MSP_Absorption /0.0/
  
! universal Debye parameter (only applied in case of potential input from slice files)
  integer*4, public :: MSP_nbuni
  real*4, public :: MSP_Buni
  DATA MSP_nbuni /0/
  DATA MSP_Buni /0.01/
  
! global supercell data
  integer*4, public :: MSP_dimcellx
  integer*4, public :: MSP_dimcelly
  real*4, public :: MSP_szcellx
  real*4, public :: MSP_szcelly


! CALCULATION RESULT
  real*4, public :: MSP_TheResult

! debug flag
  integer*4, public :: DEBUG_EXPORT
  DATA DEBUG_EXPORT /0/
  integer*4, public :: VERBO_EXPORT
  DATA VERBO_EXPORT /0/
  
! calculation modes
  !
  ! CTEM mode, switch type of input wave, (0: STEM Probe, 1: TEM Plane Wave)
  integer*4, public :: MSP_ctemmode
  DATA MSP_ctemmode /0/ ! off by default
  !
  ! IMAGE EXTRACTION, extract real space intensity distribution, (0: OFF, 1: ON)
  integer*4, public :: MSP_pimgmode
  DATA MSP_pimgmode /0/ ! off by default
  !
  ! DIFFRACTION EXTRACTION, extract Fourier-space intensity distribution, (0: OFF, 1: ON)
  integer*4, public :: MSP_pdifmode
  DATA MSP_pdifmode /0/ ! off by default
  !
  ! REMARK: If either of the two above modes is active in STEM simulations, the calculation will be done with
  !         an explicit random probe offset, simulating the chosen type of partial spatial coherence, and with
  !         an explicit random probe defocus, simulating partial temporal coherence.
  !
  ! VORTEX PROBE FLAG AND ORBITAL ANGULAR MOMENTUM
  integer*4, public :: MSP_Vortex
  DATA MSP_Vortex /0/ ! off by default
  !
  ! SWITCH between large angle propagators (0: default) or Fresnel propagators (1)
  integer*4, public :: MSP_use_fre
  DATA MSP_use_fre /1/ ! on by default
  !
  ! external override defocus parameters
  integer*4, public :: MSP_use_extdefocus
  real*4, public :: MSP_extdefocus
  DATA MSP_use_extdefocus /0/
  DATA MSP_extdefocus /0.0/
  
  ! external override object tilt parameters
  integer*4, public :: MSP_use_extot
  real*4, public :: MSP_OTExX, MSP_OTExY
  DATA MSP_use_extot /0/
  DATA MSP_OTExX /0.0/
  DATA MSP_OTExY /0.0/
  
  ! external override source profile parameters
  integer*4, public :: MSP_use_extsrcprm
  real*4, public :: MSP_extsrcrad
  DATA MSP_use_extsrcprm /0/
  DATA MSP_extsrcrad /0.01/
  
  ! external input of a wave function inserted in the multislice
  ! - flag: import activated (loads file named by MSP_inwfile)
  integer*4, public :: MSP_use_extinwave
  ! - index: thickness of insertion (this (0 based) slice index of the object stack will be applied)
  integer*4, public :: MSP_extinwslc
  ! - flag: import form of the wave function data (0: RS [default], 1: FS)
  integer*4, public :: MSP_extinwform
  DATA MSP_use_extinwave /0/
  DATA MSP_extinwslc /0/
  DATA MSP_extinwform /0/
  
!   global file names

!   measurement input file
  character(len=MSP_ll), public :: MSP_infile

!   measurement output file
  character(len=MSP_ll), public :: MSP_outfile

!   parameter file
  character(len=MSP_ll), public :: MSP_prmfile
  
!   input wave function file
  character(len=MSP_ll), public :: MSP_inwfile
  
!   PI
  real*4, public :: MSP_pi
!   conversion factor degree to radian
  real*4, public :: MSP_d2r
!   conversion factor radian to degree
  real*4, public :: MSP_r2d


! parameters defined outsides the calculation modules
!
!

! slice file title length
  integer*4, public, parameter :: MSP_SF_TITLE_LENGTH = 40

! scane frame parameters
  real*4, public :: MSP_SF_offsetx, MSP_SF_offsety, MSP_SF_sizex, MSP_SF_sizey
  real*4, public :: MSP_SF_rot, MSP_SF_rotcos, MSP_SF_rotsin
  integer*4, public :: MSP_SF_ndimx, MSP_SF_ndimy

! partial coherence parameters
!   MSP_PC_temporal | Action
!   ------------------------------------------------------------------
!   0               | Off
!   1               | Gaussian focus distribution (1/e-HW definition)
!   ------------------------------------------------------------------
!   MSP_PC_spatial  | Action
!   ------------------------------------------------------------------
!   0               | Off
!   1               | Gaussian source distribution (HWHM definition)
!   2               | Cauchy source distribution (HWHM definition)
!   3               ! Cylindrical source distribution (cut-off radius)
!   ------------------------------------------------------------------
  integer*4, public :: MSP_PC_temporal, MSP_PC_spatial

! supercell params
  integer*4, public :: MSP_SC_repeatx, MSP_SC_repeaty, MSP_SC_repeatz
  DATA MSP_SC_repeatx /1/
  DATA MSP_SC_repeaty /1/
  DATA MSP_SC_repeatz /1/
  
! frozen lattice variant data
  integer*4, public :: MSP_SLI_filenamestruct
  DATA MSP_SLI_filenamestruct /0/ ! file name structure for slice files, determines the loading type
                                  ! 0 = all slice variants in one file
                                  ! 1 = all slice variants in separate files
  integer*4, public :: MSP_FL_varnum ! number of frozen lattice variants, minimum 1
  DATA MSP_FL_varnum /1/
  integer*4, public :: MSP_FL_varcalc ! min. number of frozen lattice variants used for calculating one scan pixel or exit plane waves
  DATA MSP_FL_varcalc /1/
  
! 
  
  
! slice data
  character(len=MSP_SF_TITLE_LENGTH), dimension(:), allocatable, public :: MSP_SLC_title
  character(len=MSP_ll), public :: MSP_SLC_filenames
  integer*4, dimension(:), allocatable, public :: MSP_SLC_object
  integer*4, dimension(:,:), allocatable, public :: MSP_SLC_setup
  complex*8, dimension(:,:,:), allocatable, public :: MSP_phasegrt
  integer*4, public :: MSP_SLC_num ! number of allocated phasegratings in MSP_phasegrt
  DATA MSP_SLC_num /0/
  
! detector parameters
  integer*4, public :: MSP_detslc                           ! detector read out period [slices]
  integer*4, public :: MSP_usedetdef                        ! detector definition file, usage flag
  character(len=MSP_ll), public :: MSP_detfile              ! detecor definition file name
  integer*4, public :: MSP_detnum                           ! number of defined detectors
  integer*4, public :: MSP_detpln                           ! number of readout planes (for output handling)
  integer*4, dimension(:), allocatable, public :: MSP_detplnflg ! flag list for detector readout for each object slice (0 = no readout, 1 = readout)
  real*4, dimension(:,:), allocatable, public :: MSP_detdef ! detector definitions
  character(len=MSP_ll), dimension(:), allocatable, public :: MSP_detname ! detector names
  ! detector area data
  integer*4, public, parameter :: MSP_detrspnhdr = 3        ! size of the radial sensitivity profile header (use flag, number of data points, theta 1 pixel)
  character(len=MSP_ll), dimension(:), allocatable, public :: MSP_detrspfile ! file name defining a detector relative radial sensitivity profile
  real*4, dimension(:,:), allocatable, public :: MSP_detrspdat ! detector relative radial sensitivity profile data
  real*4, dimension(:,:), allocatable, public :: MSP_detrsphdr ! detector relative radial sensitivity profile header data
  real*4, dimension(:,:,:), allocatable, public :: MSP_detarea ! detector area setup
  integer*4, dimension(:,:,:), allocatable, public :: MSP_detcols ! detector column limit data
  real*4, dimension(:,:), allocatable, public :: MSP_detresult ! detector readout results
  DATA MSP_usedetdef /0/
  DATA MSP_detnum /0/
  DATA MSP_detpln /0/
  DATA MSP_detslc /0/
  ! detector image output flag
  integer*4, public :: MSP_detimg_output                    ! flag detector image output
  DATA MSP_detimg_output /0/
  
  ! probe image and diffraction data
  real*4, dimension(:,:,:), allocatable, public :: MSP_pimg ! probe image data for each plane
  real*4, dimension(:,:,:), allocatable, public :: MSP_pdif ! probe diffraction data for each plane
  integer*4, dimension(:), allocatable, public :: MSP_pint_nac ! # accumulations of probe intensities
  integer*4, public :: MSP_pint_num ! number of probe image intensities per scan position
  DATA MSP_pint_num /0/
  
  ! output indexing digits
  integer*4 :: MSP_nvard, MSP_nslid, MSP_nslcd, MSP_nn1d, MSP_nn2d
  DATA MSP_nvard /3/
  DATA MSP_nslid /3/
  DATA MSP_nslcd /3/
  DATA MSP_nn1d /3/
  DATA MSP_nn2d /3/
  
! temp. string
  character(len=MSP_ll), public :: MSP_stmp, MSP_stmp2
  
  
  

  CONTAINS




!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!
!* >>
!* >>
!* >> SETUP & INIT ROUTINES


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_INIT()
! function: inits the module
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 100
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_INIT: INIT."
  MSP_compile_version = MSP_VERSION_NORMAL  !MSP_VERSION_SERVICE !MSP_VERSION_NORMAL
! ------------


! ------------
! empty a line
  MSP_el = REPEAT(" ",MSP_ll)
! empty all strings, so be careful with initializing
  MSP_verApp = MSP_el
  MSP_callApp = MSP_el
  MSP_authApp = MSP_el
  MSP_prmfile = MSP_el
  MSP_outfile = MSP_el
  MSP_infile = MSP_el
  MSP_ctemmode = 0
  MSP_TheResult = 0.0
  MSP_ApplyPSC = 0
  DEBUG_EXPORT = 0
  MSP_err_num = 0
  MSP_warn_num = 0
! ------------


! ------------
!  default system select
!   system directory delimiters
  MSP_sys_pathdelim = "\/"
!   system identifier index
  
  MSP_sys_names(MSP_sys_WIN) = MSP_el
  MSP_sys_names(MSP_sys_WIN) = "WIN "
  MSP_sys_nszcnt(MSP_sys_WIN) = 1

  MSP_sys_names(MSP_sys_UNIX) = MSP_el
  MSP_sys_names(MSP_sys_UNIX) = "UNIX "
  MSP_sys_nszcnt(MSP_sys_UNIX) = 1

! *************************************
! **** SET CURRENT SYSTEM IDENT HERE!!!
  MSP_sys = MSP_sys_WIN
! **** SET CURRENT SYSTEM IDENT HERE!!!
! *************************************


! ------------
  MSP_pi = atan(1.0)*4.0
  MSP_d2r = MSP_pi / 180.0
  MSP_r2d = 1.0 / MSP_d2r
! ------------


! ------------
  MSP_ScanPixelX = 0
  MSP_ScanPixelY = 0
  MSP_BeamTiltX = 0.0
  MSP_BeamTiltY = 0.0
  MSP_OTExX = 0.0
  MSP_OTExY = 0.0
  MSP_SF_offsetx = 0.0
  MSP_SF_offsety = 0.0
  MSP_SF_sizex = 0.0
  MSP_SF_sizey = 0.0
  MSP_SF_rot = 0.0
  MSP_SF_rotcos = 1.0
  MSP_SF_rotsin = 0.0
  MSP_SF_ndimx = 0
  MSP_SF_ndimy = 0
  MSP_PC_temporal = 0
  MSP_PC_spatial = 0
  MSP_SC_repeatx = 1
  MSP_SC_repeaty = 1
  if (allocated(MSP_SLC_object)) then
    deallocate(MSP_SLC_object)
  end if
  if (allocated(MSP_phasegrt)) then
    deallocate(MSP_phasegrt)
    MSP_SLC_num = 0
  end if
  if (allocated(MSP_SLC_setup)) then
    deallocate(MSP_SLC_setup)
  end if
  if (allocated(MSP_SLC_title)) then
    deallocate(MSP_SLC_title)
  end if
  MSP_usedetdef = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_INIT: EXIT."
  return

END SUBROUTINE MSP_INIT
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_UNINIT()
! function: uninits the module
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 200
  integer*4 :: nalloc
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_UNINIT: INIT."
! ------------

! ------------
  if (allocated(MSP_SLC_object)) then
    deallocate(MSP_SLC_object,stat=nalloc)
  end if
  if (allocated(MSP_phasegrt)) then
    deallocate(MSP_phasegrt,stat=nalloc)
    MSP_SLC_num = 0
  end if
  if (allocated(MSP_SLC_title)) then
    deallocate(MSP_SLC_title,STAT=nalloc)
  end if
  if (allocated(MSP_SLC_setup)) then
    deallocate(MSP_SLC_setup,STAT=nalloc)
  end if
  if (allocated(MSP_detdef)) then
    deallocate(MSP_detdef,stat=nalloc)
  end if
  if (allocated(MSP_detname)) then
    deallocate(MSP_detname,stat=nalloc)
  end if
  if (allocated(MSP_detresult)) then
    deallocate(MSP_detresult,stat=nalloc)
  end if
  if (allocated(MSP_detarea)) then
    deallocate(MSP_detarea,stat=nalloc)
  end if
  if (allocated(MSP_detcols)) then
    deallocate(MSP_detcols,stat=nalloc)
  end if
  if (allocated(MSP_detrspdat)) then
    deallocate(MSP_detrspdat,stat=nalloc)
  end if
  if (allocated(MSP_detrspfile)) then
    deallocate(MSP_detrspfile,stat=nalloc)
  end if
  if (allocated(MSP_detrsphdr)) then
    deallocate(MSP_detrsphdr,stat=nalloc)
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_UNINIT: EXIT."
  return

END SUBROUTINE MSP_UNINIT
!**********************************************************************!


!* << SETUP & INIT ROUTINES
!* <<
!* <<
!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!



























!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!
!* >>
!* >>
!* >> DATA MANAGEMENT ROUTINES


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_READparams(nu)
! function: reads parameter data from file
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 300
  integer*4, intent(in) :: nu
  logical :: isopen
  character*256 :: buffer
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_READparams: INIT."
  INQUIRE(unit=nu,opened=isopen)
  if (.not.isopen) then
    call MSP_ERROR("File unit not connected.",subnum+1)
  end if
! ------------


! ------------
! load data from parameter file
! current parameter file format:

! ----------> Input file structure: example

! --------------
! read the data
   do
     read(unit=nu,fmt=*,err=15) buffer
     
     if (trim(buffer)=="[Microscope Parameters]") then
       if (DEBUG_EXPORT==1) then
         call PostMessage( "Start reading [Microscope Parameters]" )
       end if
       call MSP_READBLOCK_microscope(nu)
       if (DEBUG_EXPORT==1) then
         call PostMessage( "Finished reading [Microscope Parameters]" )
       end if
     end if
     
     if (trim(buffer)=="[Multislice Parameters]") then
       if (DEBUG_EXPORT==1) then
         call PostMessage( "Start reading [Multislice Parameters]" )
       end if
       call MSP_READBLOCK_multislice(nu)
       if (DEBUG_EXPORT==1) then
         call PostMessage( "Finished reading [Multislice Parameters]" )
       end if
       exit ! terminate then
     end if
   
   end do 

! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_READparams: EXIT."
  return
  
15 call CriticalError("Failed reading parameter string from file.")
  return

END SUBROUTINE MSP_READparams
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_READBLOCK_microscope(nunit)
! function: reads filename parameters from connected file
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : connected file unit
! -------------------------------------------------------------------- !
! block example and comments:
!    30.000000						(semi angle of incident wave (mrad) in STEM mode, ignored in CTEM mode)
!    70.000000						(lower semi angle of fourier-space detector (mrad), ignored in CTEM mode)
!    120.00000						(upper semi angle of fourier-space detector (mrad), ignored in CTEM mode)
!    1, 'detectors.prm'				(added 100901, switch for using a detector definition file, and the name of the detector definition file, ignored in CTEM mode)
!    3.000000						(defocus spread (nm), irgnored in CTEM mode)
!    2.000000						(defocus spread kernel width, irgnored in CTEM mode)
!    0.001969						(wavelength (nm))
!    0.010000						(source radius (nm), irgnored in CTEM and STEM mode, used only for application of partial spatial coherence to an input file)
!    7								(defocus spread kernel steps/size, ignored in CTEM mode)
!    24								([NOA] = number of aberration definitions following this line, they MUST follow)
!    0 0.000000 0.000000			(aberration definition: (index) (aberration.x) (aberration.y) )
!    1 0.000000 0.000000			( ... same up to num=NOA required)
!    2 0.000000 0.000000			(ATTENTION: unformatted reading expects the dot "." as decimal)
!    3 0.000000 0.000000            (           delimeter. DO NOT USE COMMA ","!!! It bugs the reading.)
!    4 0.000000 0.000000			(The space char " " is used to separate the 3 values read as [integer real real])
!    5 0.000000 0.000000			(aberration numbering starts from 0 up to max. aberration)
!    6 0.000000 0.000000			(current implementation supports max. 24 aberrations (1st to 8th order))
!    7 0.000000 0.000000			(beginning with 0=image shift, 1=defocus, 2=2-fold astigmatism, 3=axial coma, )
!    8 0.000000 0.000000			(4=3-fold astigm., 5=CS, 6=star aberration, 7=4-fold astigm., 8=5th order coma, )
!    9 0.000000 0.000000			(9=3-lobe aberration, )
!    10 0.000000 0.000000    		(10=5-fold astigmatism, )
!    11 0.000000 0.000000    		(11=spherical aberration of 6th order, )
!    12 0.000000 0.000000    		(12=2-fold aberration of of 6th order, )
!    13 0.000000 0.000000    		(13=rosette aberration, )
!    14 0.000000 0.000000   		(14=6-fold astigmatism, )
!    15 0.000000 0.000000    		(15=7th order coma, )
!    16 0.000000 0.000000    		(16=3-fold aberration of 7th order, )
!    17 0.000000 0.000000    		(17=5-fold aberration of 7th order, )
!    18 0.000000 0.000000    		(18=7-fold astigmatism, )
!    19 0.000000 0.000000    		(19=Spherical aberration of 8th order, )
!    20 0.000000 0.000000    		(20=2-fold aberration of of 8th order, )
!    21 0.000000 0.000000    		(21=4-fold aberration of of 8th order, )
!    22 0.000000 0.000000    		(22=6-fold aberration of of 8th order, )
!    23 0.000000 0.000000    		(23=8-fold astigmatism )
! -------------------------------------------------------------------- !
  
  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 400
  integer*4, intent(in) :: nunit
  integer*4 :: anum, idx, anumint, i
  real*4 :: ax, ay, rlamb
  logical :: isopen
  character*STF_aberration_longname_length :: aname
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_READBLOCK_microscope: INIT."
! check open status
  INQUIRE(unit=nunit,opened=isopen)
  if (.not.isopen) then
    call MSP_ERROR("File unit not connected.",subnum+1)
  end if
! ------------


! ------------
! 
  read(unit=nunit,fmt=*,err=17) STF_caperture  ! (semi angle of incident wave (mrad)
  read(unit=nunit,fmt=*,err=17) MS_detminang ! (lower semi angle of fourier-space detector (mrad))
  read(unit=nunit,fmt=*,err=17) MS_detmaxang ! (upper semi angle of fourier-space detector (mrad))
  read(unit=nunit,fmt=*,err=16) MSP_usedetdef, MSP_detfile ! (switch for using a detector definition file, and the name of the detector definition file, ignored in CTEM mode)
  read(unit=nunit,fmt=*,err=17) rlamb ! (wavelength (nm)) or (high tension (kV))
  !read(unit=nunit,fmt=*) ! (not used)
  read(unit=nunit,fmt=*,err=17) STF_srcradius ! (source radius (nm))
  read(unit=nunit,fmt=*,err=17) STF_defocusspread ! (defocus spread (nm))
  read(unit=nunit,fmt=*,err=17) STF_DEFOCUS_KERNEL_SPREAD ! (defocus spread kernel width)
  read(unit=nunit,fmt=*,err=16) STF_DEFOCUS_KERNEL_STEPS ! (defocus spread kernel steps/size)
  read(unit=nunit,fmt=*,err=16) anum ! (noa = number of aberration definitions)
  anumint = min(anum,STF_maxaberration)
  if (anumint>0) then
    do i = 1, anumint
      read(unit=nunit,fmt=*,err=16) idx,ax,ay ! (aberration definition: (index) (aberration.x) (aberration.y) )
      call STF_SetAberration(idx+1,ax,ay)
    end do
  end if
    
! ***********************
!   CHECKS
  if (rlamb>0.1) then ! assume that the the HT is specified
    STF_ht = rlamb
    ! transform to wavelength
    STF_lamb = 1.239842447 / sqrt( rlamb * ( 1022.0 + rlamb ) )
  else ! assume that wavelength is specified
    STF_lamb = rlamb ! just set it
    ! transform to high tension
    STF_ht = sqrt( 511.0**2 + (1.239842447 / rlamb)**2 ) - 511.0
  end if
  MS_lamb = STF_lamb
  MS_ht = STF_ht
  if (STF_lamb < 0.0001 .or. STF_lamb > 0.1 ) goto 18
!   CHECKS
! ***********************
    
! ***********************
!   PSEUDO DEBUG
  if (DEBUG_EXPORT==1) then
    call PostMessage("[Microscope Parameters] input report:")
    write(unit=MSP_stmp,fmt='(A,G12.4)') "semi angle of incident wave (mrad):",STF_caperture
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "lower semi angle of fourier-space detector (mrad):",MS_detminang
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "upper semi angle of fourier-space detector (mrad):",MS_detmaxang
    call PostMessage(trim(MSP_stmp))
    if (MSP_usedetdef/=0) then
      call PostMessage( "! The previous two parameters will be ignored." )
      call PostMessage( "Using detector definition from file ["//trim(MSP_detfile)//"]." )
    end if
    write(unit=MSP_stmp,fmt='(A,G12.4)') "defocus spread (nm):",STF_defocusspread
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "defocus spread kernel width:",STF_DEFOCUS_KERNEL_SPREAD
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "electron wavelength (nm):",STF_lamb
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "electron energy (keV):",STF_ht
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "source radius (nm):",STF_srcradius
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4)') "defocus spread kernel steps/size",STF_DEFOCUS_KERNEL_STEPS
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4)') "number of aberration definitions:",anumint
    call PostMessage(trim(MSP_stmp))
    if (anumint>0) then
      do idx = 1, anumint
        call STF_GetAberration(idx+1,ax,ay)
        call STF_GetAberrationLName(idx+1,aname)
        write(unit=MSP_stmp,fmt='(A,2G12.4)') "- Set "//trim(aname)//" (nm):",ax,ay
        call PostMessage(trim(MSP_stmp))
      end do
    end if
  end if
! ***********************
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_READBLOCK_microscope: EXIT."
  return
  
15 call CriticalError("Failed reading string parameter from file.")
16 call CriticalError("Failed reading integer parameter from file.")
17 call CriticalError("Failed reading float parameter from file.")
18 call CriticalError("Wavelength parameter is out of range (0.0001 ... 0.1).")
  return

END SUBROUTINE MSP_READBLOCK_microscope
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_READBLOCK_multislice(nunit)
! function: reads filename parameters from connected file
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : connected file unit
! -------------------------------------------------------------------- !
! block example and comments:
!    0.000000						(object tilt x [deg])
!    0.000000						(object tilt y [deg])
!    0.000000						(horizontal scan frame offset [nm])
!    0.000000						(vertical scan frame offset [nm])
!    1.104640						(horzontal scan frame size [nm])
!    1.170920						(vertical scan frame size [nm])
!    0.0               				(frame rotation w.r.t. slice x-axis [deg]) (added 091119)
!    100							(number of scan columns = number of pixels on horizontal image axis)
!    106							(number of scan rows = number of pixels on vertical image axis)
!    0								(de/activates partial temporal coherence calculation -> set to 1 to activate)
!    0								(de/activates partial spatial coherence calculation -> set to 1 to activate, ignored by msa)
!    1								(integer factor to internally repeat supercell data in horizontal direction, x)
!    1								(integer factor to internally repeat supercell data in vertical direction, y)
!    1								(integer factor to internally repeat supercell data in depth, z)
!    'D:\Data\Juri\STO\sto110c5x5\sli\sto110c5x5' 	(slice file name string <SFN>, slice files will be searched with names <SFN>+"_###.sli" where ### is a three digit number identifying the indiviual slices in numbered order from entrance slice 001 to exit slice [NOSD]. )
!    4								([NOSD] = number of slice definitions, SAME NUMBER OF SLICE DEFINITIONS MUST FOLLOW BELOW)
!    1								(number of frozen lattice alternatives)
!    1                              (number of frozen lattice variations for one scan pixel in STEM mode, number of frozen lattice variations generating exit plane waves in CTEM mode. )
!    0								(periodic detector readout positions, number of slices after which detectors are read out, a sequence of output files will be generated, adding index "_t###" where ### is a 3 digit number specifying the slice number of the detection.
!    8								([NOS] = number of slices in object, same number of slice IDs MUST FOLLOW!)
!    0								(slice id, must be a valid slice index from the slice list 0 = slice 001, <NOSD>-1 = slice <NOSD> as set above)
!    1								(...)
!    2
!    3
!    4
!    5
!    6
!    7
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 500

  integer*4, intent(in) :: nunit

  logical :: isopen, fex
  integer*4 :: i, status, ntmp, nerr
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_READBLOCK_multislice: INIT."
! check open status
  INQUIRE(unit=nunit,opened=isopen)
  if (.not.isopen) then
    call MSP_ERROR("File unit not connected.", subnum+1)
  end if
! ------------


! ------------
  read(unit=nunit,fmt=*,err=17) MS_objtiltx
  read(unit=nunit,fmt=*,err=17) MS_objtilty
  read(unit=nunit,fmt=*,err=17) MSP_SF_offsetx
  read(unit=nunit,fmt=*,err=17) MSP_SF_offsety
  read(unit=nunit,fmt=*,err=17) MSP_SF_sizex
  read(unit=nunit,fmt=*,err=17) MSP_SF_sizey
  read(unit=nunit,fmt=*,err=17) MSP_SF_rot
  read(unit=nunit,fmt=*,err=16) MSP_SF_ndimx
  read(unit=nunit,fmt=*,err=16) MSP_SF_ndimy
  read(unit=nunit,fmt=*,err=16) MSP_PC_temporal
  read(unit=nunit,fmt=*,err=16) MSP_PC_spatial
  read(unit=nunit,fmt=*,err=16) MSP_SC_repeatx
  read(unit=nunit,fmt=*,err=16) MSP_SC_repeaty
  read(unit=nunit,fmt=*,err=16) MSP_SC_repeatz
  read(unit=nunit,fmt=*,err=15) MSP_SLC_filenames  ! slice file name string
  read(unit=nunit,fmt=*,err=16) MS_slicenum        ! number of slices
  read(unit=nunit,fmt=*,err=16) MSP_FL_varnum      ! number of slice variants
  read(unit=nunit,fmt=*,err=16) MSP_FL_varcalc     ! number of variant loops per pixel
  read(unit=nunit,fmt=*,err=16) MSP_detslc         ! detector positions in slice stack
  read(unit=nunit,fmt=*,err=16) MS_stacksize       ! read the slice stacking order and size
  if (MS_stacksize>0) then
    if (allocated(MSP_SLC_object)) then
      deallocate(MSP_SLC_object)
    end if
    allocate(MSP_SLC_object(MS_stacksize),stat=status)
    if (status/=0) then
      call MSP_ERROR("Memory allocation failed.", subnum+2)
    end if
    do i = 1, MS_stacksize
      read(unit=nunit,fmt=*,err=16) ntmp
      MSP_SLC_object(i) = ntmp
    end do
    call MSP_GetNumberOfDigits(MS_stacksize,MSP_nslid)
    MSP_nslid = max(3,MSP_nslid)
  else
    goto 18
  end if
  call MSP_GetNumberOfDigits(MS_slicenum,MSP_nslcd)
  MSP_nslcd = max(3,MSP_nslcd) ! min. number of expected digits defining the slice number in the file name
  call MSP_GetNumberOfDigits(max(MSP_SF_ndimx,MSP_SF_ndimy),MSP_nn1d)
  MSP_nn1d = max(3,MSP_nn1d)
  
! ***********************
!   CHECKS

  
  if (sqrt(MS_objtiltx**2+MS_objtilty**2)>10.0) then
    call MSP_WARN("Very large object tilt applied (>10 degree).")
  end if
  if (MSP_SC_repeatx<1) then
    MSP_SC_repeatx = 1
    call MSP_WARN("Invalid parameter: super cell repeat along x. Set to 1.")
  end if
  if (MSP_SC_repeaty<1) then
    MSP_SC_repeaty = 1
    call MSP_WARN("Invalid parameter: super cell repeat along y. Set to 1.")
  end if
  if (MSP_SC_repeatz<1) then
    MSP_SC_repeatz = 1
    call MSP_WARN("Invalid parameter: super cell repeat along z. Set to 1.")
  end if
  if (MS_stacksize>0) then
    do i = 1, MS_stacksize
      if (MSP_SLC_object(i)>=MS_slicenum) then
        call MSP_ERROR("Detected invalid slice definition index in object slice list.",subnum+4)
      end if
    end do
  end if
  if (MSP_FL_varcalc<1) MSP_FL_varcalc = 1
  if (MSP_FL_varnum<1) MSP_FL_varnum = 1
  call MSP_GetNumberOfDigits(MSP_FL_varcalc,MSP_nvard)
  MSP_nvard = max(3,MSP_nvard) ! min. number of expected digits defining the variation number in the file name
  ! check if the slice file for slice 1 and variant 1 can be found
  MSP_SLI_filenamestruct = 0 ! test for multi-variant files structure
  call GetSliceFileName(1,1,MSP_stmp,nerr)
  inquire(file=trim(MSP_stmp),exist=fex)
  if (.not.fex) then ! file not found
    ! test for single variant file structure
    MSP_SLI_filenamestruct = 1
    call GetSliceFileName(1,1,MSP_stmp2,nerr)
    inquire(file=trim(MSP_stmp2),exist=fex)
    if (.not.fex) then
      call MSP_ERROR("First slice file ["//trim(MSP_stmp)// &
           & "] or ["//trim(MSP_stmp2)//"] not found.",subnum+5)
      goto 19
    end if
  end if

!   CHECKS
! ***********************



! ***********************
!   CALCULATIONS

  MSP_SF_rotcos = cos(MSP_SF_rot*MSP_d2r)
  MSP_SF_rotsin = sin(MSP_SF_rot*MSP_d2r)

!   CALCULATIONS
! ***********************



! ***********************
!   PSEUDO DEBUG
  if (DEBUG_EXPORT==1) then
    call PostMessage( "[Multislice Parameters] input report:" )
    write(unit=MSP_stmp,fmt='(A,G12.4)') "object tilt x [deg]:",MS_objtiltx
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "object tilt y [deg]:",MS_objtilty
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "horizontal scan frame offset [nm]:",MSP_SF_offsetx
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "vertical scan frame offset [nm]:",MSP_SF_offsety
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "horzontal scan frame size [nm]:",MSP_SF_sizex
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "vertical scan frame size [nm]:",MSP_SF_sizey
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4)') "scan frame rotation [nm]:",MSP_SF_rot
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I5)') "number of scan columns:",MSP_SF_ndimx
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I5)') "number of scan rows:",MSP_SF_ndimy
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I3)') "consider partial temporal coherence:",MSP_PC_temporal
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I3)') "consider partial spatial coherence:",MSP_PC_spatial
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4)') "repetition of supercell along x:",MSP_SC_repeatx
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4)') "repetition of supercell along y:",MSP_SC_repeaty
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4)') "repetition of supercell along z:",MSP_SC_repeatz
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I5)') "min. number of frozen lattice variations per cal.:",MSP_FL_varcalc
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I5)') "number of slice definitions:",MS_slicenum
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I5)') "number of frozen lattice variants:",MSP_FL_varnum
    call PostMessage(trim(MSP_stmp))
    call GetSliceFileName(1,1,MSP_stmp,nerr)
    call PostMessage("slice files:    ["//trim(MSP_stmp)//"]")
    call GetSliceFileName(MS_slicenum,MSP_FL_varnum,MSP_stmp,nerr)
    call PostMessage("            ... ["//trim(MSP_stmp)//"]")
    write(unit=MSP_stmp,fmt='(A,I5)') "number of slices in object:",MS_stacksize
    call PostMessage(trim(MSP_stmp))
!    if (MS_stacksize>0) then
!      write(unit=MSP_stmp,fmt='(A,<MS_stacksize>I5)') "slice order:",MSP_SLC_object(1:MS_stacksize)
!      call PostMessage(trim(MSP_stmp))
!    end if
    if (MSP_detslc<=0) then
      write(unit=MSP_stmp,fmt='(A,I4,A)') "detector readout at last slice (",MS_stacksize,")."
      call PostMessage(trim(MSP_stmp))
    else
      write(unit=MSP_stmp,fmt='(A,I4,A)') "detector readout every ",MSP_detslc," slices."
      call PostMessage(trim(MSP_stmp))
    end if
    
  end if
! ***********************
    
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_READBLOCK_multislice: EXIT."
  return
  
15 call CriticalError("Failed reading string parameter from file.")
16 call CriticalError("Failed reading integer parameter from file.")
17 call CriticalError("Failed reading float parameter from file.")
18 call CriticalError("Invalid number of object slices.")
19 call CriticalError("Failed to locate slice files.")

END SUBROUTINE MSP_READBLOCK_multislice
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_READdetdef(nu)
! function: reads parameter data from file
! -------------------------------------------------------------------- !
! parameter: integer*4 :: nu     ! logical file unit opened outsides
! -------------------------------------------------------------------- !
! Example detetctor definition file:
!    '[Detector parameters]'
!    2013041001
!    3
!    80.0, 200.0, 0.0, 0.0, 0.0, 0.0, 'HAADF', 'prm/HAADFprofile.dat'
!    30.0, 90.0, 0.0, 0.0, 0.0, 0.0, 'ADF', ''
!    12.0, 23.0, 0.0, 0.0, 0.0, 0.0, 'ABF', ''
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  integer*4, intent(in) :: nu
  logical :: isopen, fexists
  character(len=800) :: buffer
  integer*4 :: i, ddt, nalloc, nversion, nversion_found
  real*4 :: dummy, d2r
! ------------


! ------------
! INIT
  nalloc = 0
  nversion_found = 0
  d2r = 0.01745329
!  write(unit=*,fmt=*) " > MSP_READdetdef: INIT."
  INQUIRE(unit=nu,opened=isopen)
  if (.not.isopen) then
    call MSP_ERROR("File unit not connected.",subnum+1)
  end if
! ------------


! ------------
! load data from parameter file
! current parameter file format:

! ----------> Input file structure: example

! --------------
! read the data
   do
     read(unit=nu,fmt=*,err=666) buffer
     
     if (INDEX(buffer,"[Detector Parameters]")>0) then
       if (DEBUG_EXPORT==1) then
         call PostMessage( "Start reading [Detector parameters]" )
       end if
       
       read(unit=nu,fmt=*,err=555) nversion
       
       ! FOR ALL VERSIONS
       read(unit=nu,fmt=*,err=555) MSP_detnum
       write(unit=MSP_stmp,fmt='(A,I3,A)') "- reading definitions for ",MSP_detnum," detectors."
       call PostDebugMessage(trim(MSP_stmp))
       
       if (MSP_detnum>0) then ! check number of definitions
       
         call PostDebugMessage("  - allocating memory for detector parameter arrays.")
       
         ! allocate detector parameter arrays
         if (allocated(MSP_detdef)) deallocate(MSP_detdef,stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to deallocate memory of previous detector definitions.",subnum+2)
            goto 777
         end if
         allocate(MSP_detdef(0:10,MSP_detnum), stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to allocate memory for new detector definitions.",subnum+3)
            goto 777
         end if
         MSP_detdef = 0.0
         ! ***
         if (allocated(MSP_detname)) deallocate(MSP_detname,stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to deallocate memory of previous detector names.",subnum+2)
            goto 777
         end if
         allocate(MSP_detname(MSP_detnum), stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to allocate memory for new detector names.",subnum+3)
            goto 777
         end if
         MSP_detname = ""
         ! ***
         if (allocated(MSP_detrspfile)) deallocate(MSP_detrspfile,stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to deallocate memory of previous detector sensitivity profile file names.",subnum+2)
            goto 777
         end if
         allocate(MSP_detrspfile(MSP_detnum), stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to allocate memory for new detector detector sensitivity profile file names.",subnum+3)
            goto 777
         end if
         MSP_detrspfile = ""
         ! ***
         if (allocated(MSP_detrsphdr)) deallocate(MSP_detrsphdr,stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to deallocate memory of previous detector sensitivity profile headers.",subnum+2)
            goto 777
         end if
         allocate(MSP_detrsphdr(MSP_detrspnhdr,MSP_detnum), stat=nalloc)
         if (nalloc/=0) then
            call MSP_ERROR("Failed to allocate memory for new detector detector sensitivity profile headers.",subnum+3)
            goto 777
         end if
         MSP_detrsphdr = 0.0
         
         
         ! read the definitions       
         DETDEF_PRMREAD: do i=1, MSP_detnum
           ! read the complete detector definition line
           read(unit=nu,fmt='(A)',err=555) buffer
! REMOVED 0.62 !
!           ! get the detector type switch
!           read(unit=buffer,fmt=*,err=555) MSP_detdef(0,i)
!           !
!           ddt = nint(MSP_detdef(0,i))
!           !
           MSP_detdef(0,i) = 1
           ddt = 1 ! support only annular detectors
           ! process the line depending on the type switch
           TYPE_SWITCH: select case (ddt)
           case(1)

             call PostDebugMessage("  - reading an annular detector &
               &segment definition.")
             ! read the line using the correct sequence of parameters
! REMOVED 0.62 !             
!             read(unit=buffer,fmt=*,err=555) &
!               &  MSP_detdef(0,i), MSP_detdef(1,i), MSP_detdef(2,i), &
!               &  MSP_detdef(3,i), MSP_detdef(4,i), MSP_detdef(5,i), &
!               &  MSP_detdef(6,i), MSP_detname(i)

             VERSION_SWITCH: select case (nversion)
             case(2013041001)
               nversion_found = nversion
               read(unit=buffer,fmt=*,err=555) &
                 &  MSP_detdef(1,i), MSP_detdef(2,i), &
                 &  MSP_detdef(3,i), MSP_detdef(4,i), MSP_detdef(5,i), &
                 &  MSP_detdef(6,i), MSP_detname(i)
             case(2016021801)
               nversion_found = nversion
               read(unit=buffer,fmt=*,err=555) &
                 &  MSP_detdef(1,i), MSP_detdef(2,i), &
                 &  MSP_detdef(3,i), MSP_detdef(4,i), MSP_detdef(5,i), &
                 &  MSP_detdef(6,i), MSP_detname(i), MSP_detrspfile(i)
               if (len_trim( MSP_detrspfile(i) ) > 0) then
                 inquire(file=trim(MSP_detrspfile(i)),exist=fexists)
                 if (fexists) MSP_detrsphdr(1,i) = 0.1 ! indicate that a profile should be loaded
               end if
             end select VERSION_SWITCH ! case (nversion)
             
             
             ! version check
             if (0==nversion_found) goto 556
             ! preliminary parameter checks
             if (MSP_detdef(1,i)==MSP_detdef(2,i)) then
               write(unit=MSP_stmp,fmt='(A,I4,A,G12.4,A)') &
                 & "Detector #", i," has two identical radius parameters&
                 &:  R1 = R2 =", MSP_detdef(1,i)," (mrad)."
               call MSP_WARN( trim(MSP_stmp) )
             end if
             if (MSP_detdef(1,i)>MSP_detdef(2,i)) then
               write(unit=MSP_stmp,fmt='(A,I4,A,G12.4,A)') &
                 & "Found R1>R2 for annular detector #", i, &
                 & ", the values will be switched now.."
               call MSP_WARN( trim(MSP_stmp) )
               dummy = MSP_detdef(1,i)
               MSP_detdef(1,i) = MSP_detdef(2,i)
               MSP_detdef(2,i) = dummy
             end if
             ! save max. detection angle
             MSP_detdef(10,i) = MSP_detdef(2,i)
             ! parameter reading report
             write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
               & "- annular segment detector #", i,": "//trim(MSP_detname(i))
             call PostDebugMessage( trim(MSP_stmp) )
             write(unit=MSP_stmp,fmt='(A,F8.2,A,F8.2,A)') &
               & "  + R1 =", MSP_detdef(1,i),", R2 =", MSP_detdef(2,i)," (mrad)"
             call PostDebugMessage( trim(MSP_stmp) )
             write(unit=MSP_stmp,fmt='(A,F8.2,A,F8.2,A)') &
               & "  + A1 =", MSP_detdef(3,i),", A2 =", MSP_detdef(4,i)," (deg)"
             call PostDebugMessage( trim(MSP_stmp) )
             write(unit=MSP_stmp,fmt='(A,F8.2,A,F8.2,A)') &
               & "  + C1 =", MSP_detdef(5,i),", C2 =", MSP_detdef(6,i)," (mrad)"
             call PostDebugMessage( trim(MSP_stmp) )
             if (MSP_detrsphdr(1,i)>0.0) then
               call PostDebugMessage( "  + Sensitivity profile: "//trim( MSP_detrspfile(i) ) )
             end if
             cycle DETDEF_PRMREAD
           
!           case(2) ! DEACTIVATED 0.62 !
!             call PostDebugMessage("  - reading a rectangular detector segment definition.")
!             ! read the line using the correct sequence of parameters
!             read(unit=buffer,fmt=*,err=555) &
!               &  MSP_detdef(0,i), MSP_detdef(1,i), MSP_detdef(2,i), &
!               &  MSP_detdef(3,i), MSP_detdef(4,i), MSP_detdef(5,i), &
!               &  MSP_detname(i)
!             ! preliminary parameter checks
!             if (MSP_detdef(3,i)<=0.0 .or. MSP_detdef(4,i)<=0.0) then
!               write(unit=MSP_stmp,fmt='(A,I4,A,G12.4,A)') &
!                 & "Detector #", i," has invalid or zero area. &
!                 & the detector will not be used."
!               call MSP_WARN( trim(MSP_stmp) )
!             end if
!             ! save max. detection angle, using the max. abs of the rectangle corners
!             dummy = 0.0
!             ! get rotated size vectors (vx = (vxx, vxy), |vx| = DX ... )
!             vxx = MSP_detdef(3,i)*cos(d2r*MSP_detdef(5,i))
!             vxy = MSP_detdef(3,i)*sin(d2r*MSP_detdef(5,i))
!             vyx = -MSP_detdef(4,i)*sin(d2r*MSP_detdef(5,i))
!             vyy =  MSP_detdef(4,i)*cos(d2r*MSP_detdef(5,i))
!             ! offset point
!             px = MSP_detdef(1,i)
!             py = MSP_detdef(2,i)
!             dummy = max(dummy, sqrt(px*px + py*py) )
!             ! offset point moved right along the rectangle x direction
!             px = px + vxx
!             py = py + vxy
!             dummy = max(dummy, sqrt(px*px + py*py) )
!             ! ... moved up along the rect y direction
!             px = px + vyx
!             py = py + vyy
!             dummy = max(dummy, sqrt(px*px + py*py) )
!             ! .. moved left along the rect x direction
!             px = px - vxx
!             py = py - vxy
!             dummy = max(dummy, sqrt(px*px + py*py) )
!             MSP_detdef(10,i) = dummy
!             ! parameter reading report
!             write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
!               & "- rectangular segment detector #", i,": "//trim(MSP_detname(i))
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,F8.2,A,F8.2,A)') &
!               & "  + Offset position (X0,Y0) = (", MSP_detdef(1,i), &
!               & ",", MSP_detdef(2,i),") (mrad)"
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,F8.2,A,F8.2,A)') &
!               & "  + Size (DX,DY) = (", MSP_detdef(3,i),", ", &
!               & MSP_detdef(4,i),") (mrad)"
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,F8.2,A)') &
!               & "  + Orientation = ", MSP_detdef(5,i)," deg"
!             call PostDebugMessage( trim(MSP_stmp) )
!             cycle DETDEF_PRMREAD
!           
!           case(3) ! DEACTIVATED 0.62 !
!             call PostDebugMessage("  - reading a detector calibration parameter set.")
!             ! read the line using the correct sequence of parameters
!             ! <S1> = 3: <DEF>=(FN,PN,R1,CI,SL,X0,Y0,A0,DN) parameter sequence defining a calibrated detector.
!             ! Structure=(String,Int,Float,Float,Int,Float,Float,Float,String)
!             read(unit=buffer,fmt=*,err=555) &
!               &  MSP_detdef(0,i), &
!               &  MSP_detstr1(i),  MSP_detdef(1,i), MSP_detdef(2,i), &
!               &  MSP_detdef(3,i), MSP_detdef(4,i), MSP_detdef(5,i), &
!               &  MSP_detdef(6,i), MSP_detdef(7,i), MSP_detname(i)
!             ! preliminary parameter checks ? TODO elsewhere or here
!             ! saving max. detection angle will be done elsewhere, we need the wave sampling first
!             ! parameter reading report
!             write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
!               & "- calibrated detector #", i,": "//trim(MSP_detname(i))
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A)') &
!               & "  + image file = ["//trim(MSP_detstr1(i))//"]"
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,I6,A)') &
!               & "  + image size = ", nint(MSP_detdef(1,i)), &
!               & " pixels"
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
!               & "  + image scale = ",abs(MSP_detdef(2,i)), &
!               & " (mrad/pix)"
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,G12.4)') &
!               & "  + calibration # electrons = ",abs(MSP_detdef(3,i))
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,I6)') &
!               & "  + sensitivity levels = ", MSP_detdef(4,i)
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,F8.2,A,F8.2,A)') &
!               & "  + center X0 =", MSP_detdef(5,i),", Y0 =", &
!               & MSP_detdef(6,i)," (pix)"
!             call PostDebugMessage( trim(MSP_stmp) )
!             write(unit=MSP_stmp,fmt='(A,F8.2,A)') &
!               & "  + Orientation = ", MSP_detdef(7,i)," deg"
!             call PostDebugMessage( trim(MSP_stmp) )
!             cycle DETDEF_PRMREAD
!           
!           case default
!             call MSP_WARN("Unsupported detector type switch. The expected type switch range is (1, 2, 3).")
!             MSP_detdef(0,i) = 0.0 ! switch the detector off
!             cycle DETDEF_PRMREAD
           
           end select TYPE_SWITCH ! case (ddt)
           
           
           
         end do DETDEF_PRMREAD
         
         ! read detector profiles
         call MSP_READdetprofiles()
       
       else ! invalid number of detector definitions, switch to standard mode
         call MSP_WARN( "Detector definition file contains no valid definitions." )
         goto 777
       end if
       
       call PostDebugMessage( "Finished reading [Detector Parameters]" )
       
       exit ! terminate the line reading do-loop here
     
     end if
     
   end do 
   


! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_READdetdef: EXIT."
  return
  
555 call MSP_ERROR("Failed to read from detector parameter file.",subnum+3)
    goto 666
556 call MSP_ERROR("Unsupported format of detector parameter file.",subnum+4)
    goto 666
666 call MSP_ERROR("Failed to identify detector parameters.",subnum+2)
    goto 777
777 call MSP_WARN( "Switching OFF detector definition by parameter file." )
    call PostMessage("Using the standard detector definition instead." )
    call PostMessage("Use the example parameters from 'msa howto.txt' as template.")
    MSP_detnum = 0
    MSP_usedetdef = 0
  return

END SUBROUTINE MSP_READdetdef
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_READdetprofiles()
! function: Reads information from all detector profiles.
!           Allocates the detector profile array on the fly.
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1700
  integer*4 :: ioerr
  integer*4 :: nalloc
  integer*4 :: i, j
  integer*4 :: maxnrs ! max. number of radial samples in all profiles
  integer*4 :: lfu
  real*4 :: rtmp1, rtmp2
  character(len=MSP_ll) :: smsg
  external :: GetFreeLFU ! (lfu,lfu0,lfumax)
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_READdetprofiles: INIT."
  nalloc = 0
  maxnrs = 0
  lfu = 0
  ioerr = 0
! CHECKS
  if ( MSP_detnum<1 .or. .not.allocated(MSP_detrsphdr) &
     & .or. .not.allocated(MSP_detrspfile) ) return ! required data is missing, do nothing
  call GetFreeLFU(lfu,20,99)
  if ( lfu<20 ) then
    call MSP_ERROR("Failed to obtain logical file unit for loading detector profiles.", subnum+1)
    return ! can't open profiles
  end if
! ------------

! ------------
! DETERMINE NUMBER OF SAMPLES IN EACH FILE AND MAX. NUMBER OF SAMPLES
  do i=1, MSP_detnum
    if (MSP_detrsphdr(1,i)<0.05) cycle ! file doesn't exist or is not specified, try next
    open(unit=lfu,file=trim(MSP_detrspfile(i)),iostat=ioerr,action="read",share='DENYNONE')
    if (ioerr/=0) then
      call MSP_ERROR("Failed to open detector sensitivity profile ["// &
                    & trim(MSP_detrspfile(i))//"].", subnum+2)
      cycle ! try next
    end if
    !
    read(unit=lfu, fmt=*, iostat=ioerr) MSP_detrsphdr(2,i), MSP_detrsphdr(3,i) ! read the header line
    if (ioerr==0) then ! header read complete
      maxnrs = max(maxnrs, nint(MSP_detrsphdr(2,i))) ! store max number of items
      if ( nint(MSP_detrsphdr(2,i)) > 1 .and. &
         & nint(MSP_detrsphdr(2,i)) > ceiling(MSP_detrsphdr(3,i)) ) then ! reasonable data
        MSP_detrsphdr(1,i) = 1.0 ! indicate valid header
      end if
    end if
    !
    close(lfu)
  end do
! ------------


! ------------
! ALLOCATE THE PROFILE DATA ARRAY
  if (maxnrs<=0) return ! do nothing, nothing to load
  if (allocated(MSP_detrspdat)) deallocate(MSP_detrspdat, stat=nalloc)
  allocate(MSP_detrspdat(maxnrs+10,MSP_detnum), stat=nalloc)
  MSP_detrspdat = 0.0 ! initialize
! ------------


! ------------
! LOAD THE PROFILE DATA
  do i=1, MSP_detnum
    if (MSP_detrsphdr(1,i)<0.9) cycle ! file is invalid, try next
    open(unit=lfu,file=trim(MSP_detrspfile(i)),iostat=ioerr,action="read",share='DENYNONE')
    if (ioerr/=0) then
      call MSP_ERROR("Failed to open detector sensitivity profile ["// &
                    & trim(MSP_detrspfile(i))//"].", subnum+2)
      cycle ! try next
    end if
    !
    read(unit=lfu, fmt=*, iostat=ioerr) rtmp1, rtmp2 ! read the header line ! this should always work OK
    !
    ioerr = 0
    do j=1, nint(MSP_detrsphdr(2,i)) ! read all data items
      read(unit=lfu, fmt=*, iostat=ioerr) rtmp1 ! read the next data
      if (ioerr/=0) exit ! stop reading data due to error
      MSP_detrspdat(j,i) = rtmp1 ! allow also negative sensitivity ?
    end do
    if (ioerr/=0) then ! deactivate profile due to incomplete data
      MSP_detrsphdr(1,i) = 0.1
      write(unit=smsg,fmt=*) j
      call MSP_ERROR( "Failed to read data item #"//trim(adjustl(smsg))// &
                    & " of profile ["//trim(MSP_detrspfile(i))//"].", subnum+3)
    end if
    !
    close(lfu)
  end do
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_READdetprofiles: EXIT."
  return

END SUBROUTINE MSP_READdetprofiles
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_ALLOCPGR(nx,ny,nerr)
! function: allocates memory for phase gratings
! -------------------------------------------------------------------- !
! parameter: 
! INPUT:     integer*4 :: nx, ny ! size of phase gratings
! IN/OUT:    integer*4 :: nerr
! REMARKS:
! - requires a previous call of MSP_PREALLOCPGR
! - requires MSP_SLC_setup to be set up completely
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 600
  integer*4, intent(in) :: nx, ny
  integer*4, intent(inout) :: nerr
  integer*4 :: err, i, j, npgrnum
  real*4 :: mre
  
  character(len=600) :: smsg
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_ALLOCPGR: INIT."
  nerr = 0
! ------------

! ------------
  if (.not.allocated(MSP_SLC_setup)) then
    nerr = subnum+1
    call MSP_ERROR("Phase grating allocation failed, pre-allocation was not finished properly.",nerr)
    return
  end if
  if (nx<=0 .or. ny<=0 ) then
    nerr = subnum+2
    call MSP_ERROR("Phase grating allocation failed, invalid size input parameters.",nerr)
    return
  end if
! ------------

! ------------
  if (allocated(MSP_phasegrt)) then
    deallocate(MSP_phasegrt,STAT=err)
    MSP_SLC_num = 0
  end if
! ------------

! ------------
! determine the number of phasegratings to be allocated from MSP_SLC_setup
  npgrnum = 0
  do j=1, MS_slicenum
    i = MSP_SLC_setup(0,j)
    if (i<=0) then
      nerr = subnum+3
      write(unit=smsg,fmt='(A,I4)') "No phase grating set up for slice #",j
      call MSP_ERROR(trim(smsg),nerr)
      return
    end if
    npgrnum = npgrnum + i
  end do
  if (npgrnum<=0) then
    nerr = subnum+4
    call MSP_ERROR("Phase grating allocation failed, invalid number of phase gratings",nerr)
    return
  end if
! ------------

! ------------
  allocate(MSP_phasegrt(1:nx,1:ny,1:npgrnum),STAT=err)
  mre = real(nx)*real(ny)*real(npgrnum)*8.0/1024.0/1024.0
  if (err/=0) then
    nerr = subnum+5
    write(unit=smsg,fmt='(A,F8.2,A)') "Phase grating allocation failed ( ",mre," MB)."
    call MSP_ERROR(trim(smsg),nerr)
    return
  end if
  MSP_phasegrt = cmplx(0.0,0.0)
  MSP_dimcellx = nx
  MSP_dimcelly = ny
  MSP_SLC_num = npgrnum
  write(unit=smsg,fmt='(A,I4)') "Total number of allocated phase gratings: ", MSP_SLC_num
  call PostDebugMessage(trim(smsg))
  write(unit=smsg,fmt='(A,F8.2,A)') "Total allocated phase grating memory: ",mre," MB)."
  call PostDebugMessage(trim(smsg))
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_ALLOCPGR: EXIT."
  return

END SUBROUTINE MSP_ALLOCPGR
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_PREALLOCPGR(nz,nv,nerr)
! function: allocates structure arrays for phase gratings
! -------------------------------------------------------------------- !
! parameter: 
! INPUT:     integer*4 :: nz     ! number of phase gratings
!            integer*4 :: nv     ! number of phase grating variants
! IN/OUT:    integer*4 :: nerr
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  integer*4, intent(in) :: nz, nv
  integer*4, intent(inout) :: nerr
  integer*4 :: err
  
  character(len=600) :: smsg
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_PREALLOCPGR: INIT."
  nerr = 0
! ------------

! ------------
  if (nz<=0 .or. nv<=0) then
    write(unit=smsg,fmt='(A,F8.2,A)') "Phase grating pre-allocation failed, invalid size input parameters."
    call MSP_ERROR(trim(smsg),subnum+1)
    return
  end if
! ------------

! ------------
  if (allocated(MSP_SLC_setup)) then
    deallocate(MSP_SLC_setup,STAT=err)
  end if
  if (allocated(MSP_SLC_title)) then
    deallocate(MSP_SLC_title,STAT=err)
  end if
! ------------

! ------------
  allocate(MSP_SLC_setup(0:nv,1:nz),STAT=err)
  if (err/=0) then
    write(unit=smsg,fmt='(A,F8.2,A)') "Phase grating pre-allocation failed."
    call MSP_ERROR(trim(smsg),subnum+2)
    return
  end if
  MSP_SLC_setup = 0
  MS_slicenum = nz
  MSP_FL_varnum = nv
  allocate(MSP_SLC_title(1:nz),STAT=err)
  if (err/=0) then
    call MSP_ERROR("Name string allocation failed.",subnum+3)
    return
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_PREALLOCPGR: EXIT."
  return

END SUBROUTINE MSP_PREALLOCPGR
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
FUNCTION MSP_GetPGRIndex(nv, nz, nerr)
! function: returns the index of phasegrating variant nv of slice nz
!           in the central phase grating data MSP_phasegrt
! -------------------------------------------------------------------- !
! parameter:
! INPUT:    
!   integer*4, intent(in) :: nv ! variant index
!   integer*4, intent(in) :: nz ! slice index
! IN/OUTPUT:
!   integer*4, intent(inout) :: nerr ! error code, 0 = success 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 900
  integer*4 :: MSP_GetPGRIndex
  integer*4, intent(in) :: nv ! variant index
  integer*4, intent(in) :: nz ! slice index
  integer*4, intent(inout) :: nerr ! error code, 0 = success 
  integer*4 :: i, j, inv, nnv, vcur
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_GetPGRIndex: INIT."
  MSP_GetPGRIndex = 0
  write(unit=MSP_stmp,fmt='(A,I<MSP_nslcd>.<MSP_nslcd>,A,I<MSP_nvard>.<MSP_nvard>)') &
     & "   Request for phase grating index: slice #",nz,", variant #", nv
  call PostDebugMessage( trim(MSP_stmp) )
! ------------

! ------------
  if (nz<0 .or. nz>MS_slicenum .or. nv<0 .or. nv>MSP_FL_varnum) then
    nerr = subnum+1
    call MSP_ERROR("Phase grating index retrieval failed, invalid input parameters.",nerr)
    return
  end if
! ------------

! ------------
  inv = nv
  ! 1. get number of present variants for the selected slice
  j = MSP_SLC_setup(0, nz)
  if (j<=0) then
    nerr = subnum+2
    call MSP_ERROR("The requested slice contains no variant.",nerr)
    return
  end if
  ! 2. walk periodically through the variant setup channels until variant number nv is found
  nnv = 0
  i = 1
  vcur = 0
  do
    vcur = MSP_SLC_setup(i, nz)
    if (vcur>0) nnv = nnv + 1
    if (nnv==nv) exit
    i = i + 1
    if (i>MSP_FL_varnum) i = 1
  end do
! ------------

! ------------
  MSP_GetPGRIndex = vcur
  write(unit=MSP_stmp,fmt='(A,I5)') "   Returned phase grating index: ", vcur
  call PostDebugMessage( trim(MSP_stmp) )
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_GetPGRIndex: EXIT."
  return

END FUNCTION MSP_GetPGRIndex
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_ALLOCDET(nerr)
! function: allocates the detector arrays MSP_detarea and MSP_detcols
!           using the current number of detectors MSP_detnum
! -------------------------------------------------------------------- !
! parameter:
!   integer*4, intent(inout) :: nerr = error code
!                                      0 = success
!                                      or allocation error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1200
  
  integer*4, intent(inout) :: nerr
  integer*4 :: ncur1, ncur2, nnew, nalloc
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_ALLOCDET: INIT."
  nerr = 0
  ncur1 = 0
  ncur2 = 0
  nnew = MSP_detnum
! ------------


! ------------
  if (allocated(MSP_detarea)) then
    ncur1 = size(MSP_detarea,3)
  end if
  if (allocated(MSP_detcols)) then
    ncur2 = size(MSP_detcols,3)
  end if
! ------------

! ------------
  write(unit=MSP_stmp,fmt='(A,I3)') "   New number of detectors: ", nnew
  call PostDebugMessage( trim(MSP_stmp) )
  write(unit=MSP_stmp,fmt='(A,I3)') "   Old number of detectors 1: ", ncur1
  call PostDebugMessage( trim(MSP_stmp) )
  write(unit=MSP_stmp,fmt='(A,I3)') "   Old number of detectors 2: ", ncur2
  call PostDebugMessage( trim(MSP_stmp) )
  if ( ncur1/=nnew ) then ! the sizes do not match, we need new arrays
    if (ncur1/=0) then
      ! 1) deallocate existing arrays
      deallocate(MSP_detarea, stat=nalloc)
      if (nalloc/=0) then
        nerr = nalloc
        call MSP_ERROR( &
          & "Failed to deallocate memory of old detector area arrays.", &
          & subnum+nerr)
        return
      end if
    end if
    ncur1 = 0
    ! 2) allocate for new arrays
    if (nnew/=0) then
      allocate(MSP_detarea(FFT_BOUND,FFT_BOUND,nnew), stat=nalloc)
      if (nalloc/=0) then
        nerr = nalloc
        call MSP_ERROR( &
          & "Failed to allocate memory for detector area arrays.", &
          & subnum+nerr)
        return
      end if
      ncur1 = nnew
      MSP_detarea = 0
    end if
  end if
! ------------


! ------------
  if ( ncur2/=nnew ) then ! the sizes do not match, we need new arrays
    if (ncur2/=0) then
      ! 1) deallocate existing arrays
      deallocate(MSP_detcols, stat=nalloc)
      if (nalloc/=0) then
        nerr = nalloc
        call MSP_ERROR( &
          & "Failed to deallocate memory of old detector col arrays.", &
          & subnum+nerr)
        return
      end if
    end if
    ncur2 = 0
    ! 2) allocate for new arrays
    if (nnew/=0) then
      allocate(MSP_detcols(3,FFT_BOUND,nnew), stat=nalloc)
      if (nalloc/=0) then
        nerr = nalloc
        call MSP_ERROR( &
          & "Failed to allocate memory for detector col arrays.", &
          & subnum+nerr)
        return
      end if
      ncur2 = nnew
      MSP_detcols = 0
    end if
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_ALLOCDET: EXIT."
  return

END SUBROUTINE MSP_ALLOCDET
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_GetRandomProbeDefocus(fs, dz)
! function: calculates a random defocus from focus-spread parameters
! -------------------------------------------------------------------- !
! parameter:
!  real*4 :: fs (in) = focus spread parameter (1/e half width) in nm
!  real*4 :: dz (out) = random probe defocus in nm
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2000
  real*4, intent(in) :: fs
  real*4, intent(out) :: dz
  
  real*4, external :: GaussRand ! link "random.f90"
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_GetRandomProbeDefocus: INIT."
  dz = 0.0
  if (MSP_PC_temporal == 0) return ! not active, return
  if (abs(fs)==0.0) return ! not relevant, return
! ------------

! ------------
  dz = abs(fs) * GaussRand() * 0.7071068
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_GetRandomProbeDefocus: EXIT."
  return

END SUBROUTINE MSP_GetRandomProbeDefocus
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_GetRandomProbeShift(sc, dx, dy)
! function: returns a random probe shift in nm depending on source
!           size parameters
! -------------------------------------------------------------------- !
! parameter:
!   real*4 :: sc (in) = source size (HWHM) in nm
!   real*4 :: dx, dy (out) = random probe shift
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 2100
  real*4, intent(in) :: sc
  real*4, intent(out) :: dx, dy
  real*4 :: tmp1, tmp2, tmp3
  real*4, external :: UniRand, GaussRand ! link "random.f90"
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_GetRandomProbeShift: INIT."
  dx = 0.0
  dy = 0.0
  if (MSP_PC_spatial==0) return ! not active, return
  if (abs(sc)==0.0) return ! not relevant, return
! ------------

! ------------
! switch depending on source distribution type
  select case (MSP_PC_spatial)
  case (1) ! Gaussian
    dx = sc * GaussRand() * 0.849322
    dy = sc * GaussRand() * 0.849322
  case (2) ! Cauchy
    tmp1 = UniRand()*0.9 ! limit the distribution to approx 10*HWHM
    tmp2 = UniRand()*6.283 ! azimuth
    tmp3 = 1.3*sc*sqrt( ((2.0 - tmp1)*tmp1) / ((tmp1 - 1.0)*(tmp1 - 1.0)) ) ! radial
    dx = tmp3 * cos( tmp2 )
    dy = tmp3 * sin( tmp2 )
  case (3) ! Disk
    do
      dx = sc*(2.0*UniRand()-1.0) ! uniform x
      dy = sc*(2.0*UniRand()-1.0) ! uniform y
      if (sqrt(dx*dx+dy*dy)<=sc) exit ! point is in disk, good
    end do
  end select ! case (MSP_PC_spatial)
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_GetRandomProbeShift: EXIT."
  return

END SUBROUTINE MSP_GetRandomProbeShift
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_GetNumberOfDigits(maxidx,nd)
! function: resturns number of digits required to print maxidx
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1300
  integer*4, intent(in) :: maxidx
  integer*4, intent(inout) :: nd
  character(len=40) :: snum
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_GetNumberOfDigits: INIT."
! ------------


! ------------
  write(unit=snum,fmt='(i)') maxidx
  nd = LEN_TRIM( adjustl ( snum ) )
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_GetNumberOfDigits: EXIT."
  return

END SUBROUTINE MSP_GetNumberOfDigits
!**********************************************************************!


  
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_InitTextOutput(nerr)
! function: Initializes text output by creating a new file with the
!           specified output file name and by writing the global header.
!           Initializes program parameters used for saving records to
!           this file.
! -------------------------------------------------------------------- !
! parameter:
!  integer*4 :: nerr        = error code (0 = success)
!
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1400
  integer*4, intent(inout) :: nerr

  integer*4 :: i, k, lfu, ioerr, ndet, nslcmax, nslc, ndetect
  
  character(len=MSP_ll) :: sline
  
  external :: GetFreeLFU ! (lfu,lfu0,lfumax)
  external :: createfilefolder
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_InitTextOutput: INIT."
  nerr = 0
  lfu = 0
  if (MSP_txtout==0) return ! Do not initialize this output option.
                            ! Just leave.
! Determine the number of detectors
  ndet = max(1,MSP_detnum)
! Determine the number of output slices -> nslc
  nslcmax = MS_stacksize
  ndetect = nslcmax
  if (MSP_detslc>0) ndetect = min(MSP_detslc,nslcmax)
  nslc = 0
  do i=1, nslcmax
    if (0==modulo(i,ndetect)) nslc = nslc + 1
  end do
! ------------

! ------------
! Create a new file
  call GetFreeLFU(lfu,20,99)
  call createfilefolder(trim(MSP_outfile),ioerr)
  open(unit=lfu, file=trim(MSP_outfile), iostat=ioerr, status="replace", action="write" )
  if (ioerr/=0) then
    call CriticalError("InitTextOutput: Failed to create file ["//trim(MSP_outfile)//"].")
  end if
! ------------

! ------------
! Write the detector list to file
  write(unit=sline, fmt='(I)') ndet                            ! line 1: number of detectors
  write(unit=lfu, fmt='(A)') trim(adjustl(sline))
  do k=1, ndet                                                 ! lines 2 - 1+ndet: detector names
    write(unit=lfu, fmt='(A)') "'"//trim(MSP_detname(k))//"'"
  end do
! ------------

! ------------
! Write the slice list to file
  write(unit=sline, fmt='(I)') nslc                            ! line 2+ndet: number of detection planes ( = length of detector arrays)
  write(unit=lfu, fmt='(A)') trim(adjustl(sline))
  do i=1, nslcmax                                              ! line 3+ndet - 2+ndet+nslc: detection plane indices
    if (0/=modulo(i,ndetect)) cycle ! skip this slice
    write(unit=sline, fmt='(I)') i
    write(unit=lfu, fmt='(A)') trim(adjustl(sline))
  end do
! ------------

! ------------
! Close the file
  close(unit=lfu, iostat=ioerr)
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_InitTextOutput: EXIT."
  return

END SUBROUTINE MSP_InitTextOutput
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_WriteTextOutput(nerr)
! function: Writes (adds) a new record to a text output file.
! -------------------------------------------------------------------- !
! parameter: 
!  integer*4 :: nerr        = error code (0 = success)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1500
  integer*4, intent(inout) :: nerr

  integer*4 :: i, k, lfu, ioerr, ndet, nslcmax, ndetect
  
  character(len=MSP_ll) :: sline, stmp1, stmp2
  
  external :: GetFreeLFU ! (lfu,lfu0,lfumax)
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_WriteTextOutput: INIT."
  nerr = 0
  lfu = 0
  if (MSP_txtout==0) return ! Do not output. Just leave.
! Determine the number of detectors
  ndet = max(1,MSP_detnum)
! Determine the number of output slices
  nslcmax = MS_stacksize
  ndetect = nslcmax
  if (MSP_detslc>0) ndetect = min(MSP_detslc,nslcmax)
! ------------


! ------------
! Open the existing file
  call GetFreeLFU(lfu,20,99)
  ! no folder creation here, the file MUST exist, thus the folder
  open(unit=lfu, file=trim(MSP_outfile), iostat=ioerr, status="old", &
     & position="append", action="write" )
  if (ioerr/=0) then
    call CriticalError("WriteTextOutput: Failed to open file ["//trim(MSP_outfile)//"].")
  end if
! ------------


! ------------
! Write the current record to the file
  write(unit=lfu, fmt='(A)') "1" ! new record indicator
  write(unit=stmp1, fmt='(I)') MSP_ScanPixelX
  write(unit=stmp2, fmt='(I)') MSP_ScanPixelY
  write(unit=lfu, fmt='(A)') trim(adjustl(stmp1))//", "//trim(adjustl(stmp2)) ! scan pixel numbers
  do k=1, ndet ! loop over detectors
    do i=1, nslcmax ! loop over thickness
      if (0/=modulo(i,ndetect)) cycle
      write(unit=sline, fmt='(E14.6)') MSP_detresult(k,i)
      write(unit=lfu, fmt='(A)') trim(adjustl(sline)) ! the intensity data
    end do
  end do
! ------------


! ------------
! Close the file
  close(unit=lfu, iostat=ioerr)
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_WriteTextOutput: EXIT."
  return

END SUBROUTINE MSP_WriteTextOutput
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_FinishTextOutput(nerr)
! function: Adss the finish record to the output text file.
! -------------------------------------------------------------------- !
! parameter: 
!  integer*4 :: nerr        = error code (0 = success)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1600
  integer*4, intent(inout) :: nerr

  integer*4 :: lfu, ioerr
  
  external :: GetFreeLFU ! (lfu,lfu0,lfumax)
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_FinishTextOutput: INIT."
  nerr = 0
  lfu = 0
  if (MSP_txtout==0) return ! Do not output. Just leave.
! ------------


! ------------
! Open the existing file
  call GetFreeLFU(lfu,20,99)
  ! no folder creation here, the file MUST exist, so the folder
  open(unit=lfu, file=trim(MSP_outfile), iostat=ioerr, status="old", &
     & position="append", action="write" )
  if (ioerr/=0) then
    call CriticalError("FinishTextOutput: Failed to open file ["//trim(MSP_outfile)//"].")
  end if
! ------------

! ------------
! Write the final record to the file
  write(unit=lfu, fmt='(A)') "0" ! final record indicator
! ------------

! ------------
! Close the file
  close(unit=lfu, iostat=ioerr)
! ------------

! ------------
!  write(unit=*,fmt=*) " > MSP_FinishTextOutput: EXIT."
  return

END SUBROUTINE MSP_FinishTextOutput
!**********************************************************************!



!* << DATA MANAGEMENT ROUTINES
!* <<
!* <<
!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!














!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!
!* >>
!* >>
!* >> CALCULATIONS & PROCESSING

!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_SetAnnularDetectors(nerr)
! function: Setup of an annular detectors for the current slice dimension
! -------------------------------------------------------------------- !
! parameter:
!   INPUT:
!   IN/OUTPUT:
!     integer*4 :: nerr         ! error code
!                               ! 0 = success
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 700

  integer*4, intent(inout) :: nerr
  
  integer*4 :: nalloc
  integer*4 :: i, j, k, ndimx, ndim2x, ndimy, ndim2y, m, l, i1, j1
  real*4 :: r0, r1, a0, a1, c1, c2, itogx, itogy, itowx, itowy, itowr
  real*4 :: wt0, wt02, wt1, wt12, wx, wy, wx2, w2, ac, wc1, wc2
  real*4 :: rcur, rk, scur
  character(len=1024) :: stmp
  real*4, allocatable :: detfunc(:,:)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_SetAnnularDetectors: INIT."
  nerr = 0
  !if (MSP_usedetdef/=1) return
  !if (MSP_detnum<2) return
  ! check allocation status
  if (MS_status<1) then
    nerr = 1
    call MSP_ERROR("Multislice Module not initialized.",subnum+nerr)
    return
  end if
  if (MS_dimx<=0.or.MS_dimy<=0) then
    nerr = 2
    call MSP_ERROR("Invalid slice array size.",subnum+nerr)
    return
  end if
  if (.not.(allocated(MSP_detarea).and.allocated(MSP_detcols))) then
    nerr = 3
    call MSP_ERROR("Memory for detector data is not allocated.",subnum+nerr)
    return
  end if
  if (MSP_detnum<1) then
    nerr = 4
    call MSP_ERROR("No detector is defined.",subnum+nerr)
    return
  end if
! ------------

! ------------
! loop through detectors
  do k=1, MSP_detnum
  
  if (1/=nint(MSP_detdef(0,k))) cycle ! this is not an annular detector
  
  r0 = MSP_detdef(1,k)
  r1 = MSP_detdef(2,k)
  a0 = MSP_detdef(3,k)
  a1 = MSP_detdef(4,k)
  c1 = MSP_detdef(5,k)
  c2 = MSP_detdef(6,k)
  
  ! add 360° to the azimuths until both are definitely positive
  do while (a0<0.0 .or. a1<0.0)
    a0 = a0 + 360.0
    a1 = a1 + 360.0
  end do
  
  ! a1 should alsway be larger than a0, if not, than add 360° to a1
  if (a0>=a1) then
    a1 = a1 + 360.0
  end if
    
  ! now bring both angles back to the first rotation cycle, at least for a0
  do while (a0>=360.0)
    a0 = a0 - 360.0
    a1 = a1 - 360.0
  end do
  
  if (DEBUG_EXPORT==1) then
    write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
      & "Setting annular segment detector data for detector #",k, &
      & " "//trim(MSP_detname(k))
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- inner detector radius: ", r0," (mrad)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- outer detector radius: ", r1," (mrad)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- start azimuth angle: ", a0," (deg)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- stop azimuth angle: ", a1," (deg)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- de-center x: ", c1," (mrad)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- de-center y: ", c2," (mrad)"
    call PostDebugMessage(trim(MSP_stmp))
    if (MSP_detrsphdr(1,k)>0.9) then
      call PostDebugMessage("- detector radial sensitivity profile: "//trim(MSP_detrspfile(k)))
    else
      call PostDebugMessage("- 100% sensitivity")
    end if
  end if

! ------------
! prepare parameters
  wt0 = r0*0.001
  wt02 = wt0*wt0
  wt1 = r1*0.001
  wt12 = wt1*wt1
  wc1 = c1*0.001
  wc2 = c2*0.001
  ndimx = MS_dimx
  ndim2x = ndimx/2
  ndimy = MS_dimy
  ndim2y = ndimy/2
  
  ! set Fourier-space sampling
  itogx = 0.5/MS_samplingx/real(ndim2x)
  itowx = itogx*MS_lamb
  itogy = 0.5/MS_samplingy/real(ndim2y)
  itowy = itogy*MS_lamb
  
  ! pre-set Fourier-space sampling for the radial sensitivity profile
  itowr = 0.0
  if (MSP_detrsphdr(1,k)>0.9) then ! use this detector with sensitivity profile
    itowr = 0.001*r0 / (MSP_detrsphdr(3,k) - 1.0 ) ! theta-1 / (associated pixel number - 1) ... and from mrad to rad / pixel
  end if
  
  ! preset sensitivity
  scur = 1.0
  
  ! loop through calculation array in fourier space
  do j=1,ndimx
    MSP_detcols(2,j,k) = ndimy ! set left detector column to maximum
    wx = itowx*MS_TABBED_SCR(j) - wc1 ! get wave frequency-x
    wx2 = wx*wx
    do i=1, ndimy
      wy = itowy*MS_TABBED_SCR2(i) - wc2 ! get wave frequency-y
      w2 = wx2+wy*wy
      ac = atan2(wy, wx)*MSP_r2d
      ! try to find azimuthal cycle beyond a0
      do while (ac<a0)
        ac = ac + 360.0
      end do
      ! range checks on the following general interval : min <= x < max
      if ((w2>=wt02).and.(w2<wt12).and.(ac>=a0).and.(ac<a1)) then ! check if current fourier pixel is inside the annular segement
        if (itowr>0.0) then ! use radial sensitivity profile
          ! get the sensitivity for this theta angle
          rcur = max( 0.0, min( MSP_detrsphdr(2,k)-1.0, sqrt(w2)/itowr ) ) ! zero based index in the sensitivity table
          i1 = floor(rcur) ! lower integer (index)
          rk = rcur-real(i1) ! fraction between actual radius and integer radius
          ! access the sensitivity curve at (i1 + 1)
          scur = (1.0-rk)*MSP_detrspdat(i1+1,k) + rk*MSP_detrspdat(i1+2,k) ! linear interpolation of the sensitivity curve
        end if
        MSP_detarea(i,j,k) = scur ! set detector pixel weight to relative sensitivity
        MSP_detcols(1,j,k) = MSP_detcols(1,j,k) + 1 ! increase pixel sum
        MSP_detcols(2,j,k) = min(MSP_detcols(2,j,k),i) ! update minimum column index
        MSP_detcols(3,j,k) = max(MSP_detcols(3,j,k),i) ! update maximum column index
      end if
    end do ! loop i
  end do ! loop j

  if (MSP_detimg_output==1) then ! export the detector function to binary file
    if (.not.allocated(detfunc)) then
      allocate(detfunc(MS_dimx, MS_dimy), stat=nalloc)
      if (nalloc/=0) then
        nerr = 5
        call MSP_ERROR("Allocation of detector function array failed.",subnum+nerr)
        goto 667
      end if
    end if
    ! create the detector function from prepared arrays, but in non-transposed and non-scrambled order
    do j=1, ndimy
      j1 = MS_TABBED_SCR2(j) + ndim2y+1
      do i=1, ndimx
        i1 = MS_TABBED_SCR(i) + ndim2x+1
        detfunc(i,j) = MSP_detarea(j1,i1,k)
      end do 
    end do
    
    m = LEN_TRIM(MSP_outfile) ! get length of the standard output file
    l = INDEX(trim(MSP_outfile),".",back=.TRUE.) ! search for extension point in given output file
    if (MSP_usedetdef/=0) then ! update output file name with detetctor name
      if (l<1) then ! no extension wanted
        write(unit=stmp,fmt='(A)') trim(MSP_outfile)//"_detector_"//trim(MSP_detname(k))
      else ! last extension starts at position l, insert the index string before
        write(unit=stmp,fmt='(A)') MSP_outfile(1:l-1)//"_detector_"//trim(MSP_detname(k))//MSP_outfile(l:m)
      end if
    else ! keep current output file name
      if (l<1) then ! no extension wanted
        write(unit=stmp,fmt='(A)') trim(MSP_outfile)//"_detector"
      else ! last extension starts at position l, insert the index string before
        write(unit=stmp,fmt='(A)') MSP_outfile(1:l-1)//"_detector"//MSP_outfile(l:m)
      end if
    end if
    
    write(unit=MSP_stmp,fmt='(A,I3.3,A)') "Saving detector function #",k," to file ["//trim(stmp)//"]."
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4,A,I4)') "- detector function: 32-bit float, size: ",MS_dimx," x ",MS_dimy
    call PostDebugMessage(trim(MSP_stmp))
    call SaveDataR4(trim(stmp),detfunc,MS_dimx*MS_dimy,nerr)
    if (nerr/=0) then
      call MSP_ERROR("Failed to save detector function.",subnum+6)
    end if
    
    
    deallocate(detfunc, stat=nalloc)
    
667 end if ! (MSP_detimg_output==1)
  
  
  end do ! loop k for detectors
  if (allocated(detfunc)) deallocate(detfunc, stat=nalloc)
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_SetAnnularDetectors: EXIT."
  return

END SUBROUTINE MSP_SetAnnularDetectors
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_SetRectangularDetectors(nerr)
! function: Sets data for rectangular segment detectors
!           in MSP_detarea and MSP_detcols
!           Requires a prior allocation of the two arrays
! -------------------------------------------------------------------- !
! parameter: 
!  integer*4, intent(in) :: nerr = error code
!                                  0 = success
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1000
!
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nalloc
  integer*4 :: i, j, k, ndimx, ndim2x, ndimy, ndim2y, m, l, i1, j1
  real*4 :: x0, y0, dx, dy, a0, ca, sa, r0
  real*4 :: itogx, itogy, itowx, itowy
  real*4 :: wx, wy, tx, ty
  character(len=1024) :: stmp
  real*4, allocatable :: detfunc(:,:)
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_SetRectangularDetectors: INIT."
  nerr = 0
  !if (MSP_usedetdef/=1) return
  !if (MSP_detnum<2) return
  ! check allocation status
  if (MS_status<1) then
    nerr = 1
    call MSP_ERROR("Multislice Module not initialized.",subnum+nerr)
    return
  end if
  if (MS_dimx<=0.or.MS_dimy<=0) then
    nerr = 2
    call MSP_ERROR("Invalid slice array size.",subnum+nerr)
    return
  end if
  if (.not.(allocated(MSP_detarea).and.allocated(MSP_detcols))) then
    nerr = 3
    call MSP_ERROR("Memory for detector data is not allocated.",subnum+nerr)
    return
  end if
  if (MSP_detnum<1) then
    nerr = 4
    call MSP_ERROR("No detector is defined.",subnum+nerr)
    return
  end if
! ------------


! ------------
! loop through detectors
  do k=1, MSP_detnum
  
  if (2/=nint(MSP_detdef(0,k))) cycle ! this is not a rectangular detector
  
  x0 = MSP_detdef(1,k)
  y0 = MSP_detdef(2,k)
  dx = MSP_detdef(3,k)
  dy = MSP_detdef(4,k)
  a0 = MSP_detdef(5,k)
  
  if (DEBUG_EXPORT==1) then
    write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
      & "Setting rectangular segment detector data for detector #",k, &
      & " "//trim(MSP_detname(k))
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A,G12.4,A)') &
      & "- offset point: ( ",x0,", ",y0,") (mrad)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A,G12.4,A)') &
      & "- size: (", dx,", ", dy, ") (mrad)"
    call PostDebugMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
      & "- orientation: ", a0," (deg)"
    call PostDebugMessage(trim(MSP_stmp))
  end if

! ------------
! prepare parameters
  ! rectangle unit vectors (matrix) = ( (  exx,  exy ), (  eyx,  eyy ) )
  !                        inverse = ( (  eyy, -exy ), ( -eyx,  exx ) )
  ! exx = eyy = ca = cos(a0)
  ! exy = -eyx = sa = sin(a0)
  r0 = MSP_d2r*a0
  ca = cos(r0)
  sa = sin(r0)
  
  ndimx = MS_dimx
  ndim2x = ndimx/2
  ndimy = MS_dimy
  ndim2y = ndimy/2
  
  ! set Fourier-space sampling
  itogx = 0.5/MS_samplingx/real(ndim2x)
  itowx = itogx*MS_lamb
  itogy = 0.5/MS_samplingy/real(ndim2y)
  itowy = itogy*MS_lamb
  
  ! loop through calculation array in fourier space
  do j=1, ndimx
    MSP_detcols(2,j,k) = ndimy ! set left detector column to maximum
    wx = itowx*MS_TABBED_SCR(j) - x0 ! get diffraction angle wx and substract rect offset
    do i=1, ndimy
      wy = itowy*MS_TABBED_SCR2(i) - y0 ! get diffraction angle wy and substract rect offset
      !
      ! decompose along the orientation vectors by
      ! wx = tx*exx + ty*eyx
      ! wy = tx*exy + ty*eyy
      ! =>
      ! tx = wx*eyy - wy*eyx = wx*ca + wy*sa
      ! ty = wx*exy - wy*exx = wx*sa - wy*ca
      tx = wx*ca + wy*sa
      ty = wx*sa - wy*ca
      !
      ! check for position inside the rectangle
      ! range checks on the following general interval : min <= x < max
      if ((tx>=0.0).and.(tx<dx).and.(ty>=0.0).and.(ty<dy)) then ! yes = point is in rect
        MSP_detarea(i,j,k) = 1 ! set detector pixel weight to 1
        MSP_detcols(1,j,k) = MSP_detcols(1,j,k) + 1 ! increase pixel sum
        MSP_detcols(2,j,k) = min(MSP_detcols(2,j,k),i) ! update minimum column index
        MSP_detcols(3,j,k) = max(MSP_detcols(3,j,k),i) ! update maximum column index
      end if
    end do ! loop i
  end do ! loop j

  if (MSP_detimg_output==1) then ! export the detector function to binary file
    if (.not.allocated(detfunc)) then
      allocate(detfunc(MS_dimx, MS_dimy), stat=nalloc)
      if (nalloc/=0) then
        nerr = 5
        call MSP_ERROR("Allocation of detector function array failed.",subnum+nerr)
        goto 667
      end if
    end if
    ! create the detector function from prepared arrays, but in non-transposed and non-scrambled order
    do j=1, ndimy
      j1 = MS_TABBED_SCR2(j) + ndim2y+1
      do i=1, ndimx
        i1 = MS_TABBED_SCR(i) + ndim2x+1
        detfunc(i,j) = MSP_detarea(j1,i1,k)
      end do 
    end do
    
    m = LEN_TRIM(MSP_outfile) ! get length of the standard output file
    l = INDEX(trim(MSP_outfile),".",back=.TRUE.) ! search for extension point in given output file
    if (MSP_usedetdef/=0) then ! update output file name with detetctor name
      if (l<1) then ! no extension wanted
        write(unit=stmp,fmt='(A)') trim(MSP_outfile)//"_detector_"//trim(MSP_detname(k))
      else ! last extension starts at position l, insert the index string before
        write(unit=stmp,fmt='(A)') MSP_outfile(1:l-1)//"_detector_"//trim(MSP_detname(k))//MSP_outfile(l:m)
      end if
    else ! keep current output file name
      if (l<1) then ! no extension wanted
        write(unit=stmp,fmt='(A)') trim(MSP_outfile)//"_detector"
      else ! last extension starts at position l, insert the index string before
        write(unit=stmp,fmt='(A)') MSP_outfile(1:l-1)//"_detector"//MSP_outfile(l:m)
      end if
    end if
    
    write(unit=MSP_stmp,fmt='(A,I3.3,A)') "Saving detector function #",k," to file ["//trim(stmp)//"]."
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt='(A,I4,A,I4)') "- detector function: 32-bit float, size: ",MS_dimx," x ",MS_dimy
    call PostDebugMessage(trim(MSP_stmp))
    call SaveDataR4(trim(stmp),detfunc,MS_dimx*MS_dimy,nerr)
    if (nerr/=0) then
      call MSP_ERROR("Failed to save detector function.",subnum+6)
    end if
    
    
    deallocate(detfunc, stat=nalloc)
    
667 end if ! (MSP_detimg_output==1)
  
  
  end do ! loop k for detectors
  if (allocated(detfunc)) deallocate(detfunc, stat=nalloc)
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_SetRectangularDetectors: EXIT."
  return

END SUBROUTINE MSP_SetRectangularDetectors
!**********************************************************************!




!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE MSP_SetCalibratedDetectors(nerr)
!! function: Sets data for calibrated detectors
!!           in MSP_detarea and MSP_detcols
!!           Requires a prior allocation of the two arrays
!! -------------------------------------------------------------------- !
!! parameter: 
!!  integer*4, intent(in) :: nerr = error code
!!                                  0 = success
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! DECLARATION
!  integer*4, parameter :: subnum = 1300
!!
!  integer*4, intent(inout) :: nerr
!  
!  logical :: file_exists
!  integer*4 :: nalloc
!  integer*4 :: i, j, k, ndimx, ndim2x, ndimy, ndim2y, m, l, i1, j1
!  integer*4 :: ndimd, ndimd2, nlev
!  real*4 :: itodw, dscl, a0, x0, y0, ci, r0, ca, sa
!  real*4 :: itogx, itogy, itowx, itowy
!  real*4 :: wx, wy, tx, ty, ac, di, dj, dv
!  character(len=1024) :: stmp, dfile, dname
!  real*4, allocatable :: detdata(:,:), detused(:,:)
!  real*4, allocatable :: detlev(:,:), detfunc(:,:)
!! ------------
!
!
!! ------------
!! INIT
!!  write(unit=*,fmt=*) " > MSP_SetCalibratedDetectors: INIT."
!  nerr = 0
!  !if (MSP_usedetdef/=1) return
!  !if (MSP_detnum<2) return
!  ! check allocation status
!  if (MS_status<1) then
!    nerr = 1
!    call MSP_ERROR("Multislice Module not initialized.",subnum+nerr)
!    return
!  end if
!  if (MS_dimx<=0.or.MS_dimy<=0) then
!    nerr = 2
!    call MSP_ERROR("Invalid slice array size.",subnum+nerr)
!    return
!  end if
!  if (.not.(allocated(MSP_detarea).and.allocated(MSP_detcols))) then
!    nerr = 3
!    call MSP_ERROR("Memory for detector data is not allocated.",subnum+nerr)
!    return
!  end if
!  if (MSP_detnum<1) then
!    nerr = 4
!    call MSP_ERROR("No detector is defined.",subnum+nerr)
!    return
!  end if
!! ------------
!
!
!! ------------
!! loop through detectors
!  do k=1, MSP_detnum
!  
!  if (3/=nint(MSP_detdef(0,k))) cycle ! this is not the correct type of detector
!  
!  dname = trim(MSP_detname(k))
!  dfile = trim(MSP_detstr1(k))
!  ndimd = nint(abs(MSP_detdef(1,k)))
!  dscl = abs(MSP_detdef(2,k))
!  ci = abs(MSP_detdef(3,k))
!  nlev = nint(abs(MSP_detdef(4,k)))
!  x0 = MSP_detdef(5,k)
!  y0 = MSP_detdef(6,k)
!  a0 = MSP_detdef(7,k)
!  
!  ! parameter check
!  ! - file exists
!  inquire(file=trim(dfile),exist=file_exists)
!  if (.not.file_exists) then
!    write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
!      & "Failed to setup detector #",k,", ("//trim(dname)// &
!      & "), file ["//trim(dfile)//"] not found."
!    nerr = 5
!    call MSP_ERROR(trim(MSP_stmp),subnum+nerr)
!    cycle
!  end if
!  ! - reasonable file size (128 - 16384)
!  if (ndimd<128 .or. ndimd>16384) then
!    write(unit=stmp,fmt='(I)') ndimd
!    write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
!      & "Failed to setup detector #",k,", ("//trim(dname)// &
!      & "), invalid array size ("//trim(adjustl(stmp))// &
!      & "). Valid sizes 128 ... 16384."
!    nerr = 6
!    call MSP_ERROR(trim(MSP_stmp),subnum+nerr)
!    cycle
!  end if
!  ! - non-zero image scale
!  if (dscl<=0.0) then
!    nerr = 7
!    call MSP_ERROR("Invalid detector image scale. "// &
!      & "Expecting a non-negative number in mrad per pixel.", &
!      & subnum+nerr)
!    cycle
!  end if
!  ! - non-zero reference charge
!  if (ci<=0.0) then
!    nerr = 8
!    call MSP_ERROR("Invalid reference charge parameter. "// &
!      & "Expecting a non-negative number of electrons per pixel.", &
!      & subnum+nerr)
!    cycle
!  end if
!  ! - non-negative number of levels
!  if (nlev<=0) nlev = 1
!  
!  
!  
!  if (DEBUG_EXPORT==1) then
!    write(unit=MSP_stmp,fmt='(A,I3.3,A)') &
!      & "Setting calibrated detector data for detector #",k, &
!      & ": "//trim(dname)
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A)') &
!      & "- detector image file name: ["//trim(dfile)//"]"
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,I5,A,I5,A,I5,A)') &
!      & "- detector image size: ",ndimd," x ",ndimd, &
!      & ", expected file size: ", &
!      & nint(real(4*ndimd*ndimd)/1024.0)," kB."
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
!      & "- image scale: ", dscl," (mrad/pix)"
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
!      & "- image ref. charge: ", ci," electrons"
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,I5)') &
!      & "- used # response levels: ", nlev
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,G12.4)') &
!      & "- center pixel x: ", x0
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,G12.4)') &
!      & "- center pixel y: ", y0
!    call PostDebugMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,G12.4,A)') &
!      & "- azimuth angle wrt wave x-axis: ", a0," (deg)"
!    call PostDebugMessage(trim(MSP_stmp))
!  end if
!
!! ------------
!! prepare parameters
!  itodw = dscl
!  ndimd2 = ndimd/2
!  ! rectange unit vectors (matrix) = ( (  exx,  exy ), (  eyx,  eyy ) )
!  !                        inverse = ( (  eyy, -exy ), ( -eyx,  exx ) )
!  ! exx = eyy = ca = cos(a0)
!  ! exy = -eyx = sa = sin(a0)
!  r0 = MSP_d2r*a0
!  ca = cos(r0)
!  sa = sin(r0)
!  !
!  ndimx = MS_dimx
!  ndim2x = ndimx/2
!  ndimy = MS_dimy
!  ndim2y = ndimy/2
!  
!  ! set Fourier-space sampling of the wave function
!  itogx = 0.5/MS_samplingx/real(ndim2x)
!  itowx = itogx*MS_lamb
!  itogy = 0.5/MS_samplingy/real(ndim2y)
!  itowy = itogy*MS_lamb
!  
!  ! allocate arrays for the detector image analysis
!  if (allocated(detdata)) deallocate(detdata, stat=nalloc)
!  if (allocated(detused)) deallocate(detused, stat=nalloc)
!  allocate(detdata(ndimd,ndimd), detused(ndimy, ndimx), stat=nalloc)
!  if (nalloc/=0) then
!    nerr = 9
!    call MSP_ERROR("Failed to allocate memory for det. image data.", &
!      & subnum+nerr)
!    cycle
!  end if
!  detdata = 0.0
!  detused = 0.0
!  
!  ! load the detector data
!  call LoadDataR4(trim(dfile),detdata,ndimd*ndimd,nerr)
!  if (nerr/=0) then
!    nerr = 10
!    call MSP_ERROR("Failed to load detector image data.", &
!      & subnum+nerr)
!    cycle
!  end if
!  
!  ! allocate array for the detector sensitivity levels
!  if (allocated(detlev)) deallocate(detlev, stat=nalloc)
!  allocate(detlev(1:2,0:nlev), stat=nalloc)
!  if (nalloc/=0) then
!    nerr = 11
!    call MSP_ERROR("Failed to allocate memory for det. data.", &
!      & subnum+nerr)
!    cycle
!  end if
!  detlev = 0.0
!  
!  ! determine the detector levels
!  ! - begin with the zero level and its variation
!  !   we expect the corners outside the biggest circle to
!  !   contain only zero-electron data, i.e. background counts
!  !   => Sum up the corner intensities and get mean value
!  !      and standard deviation (RMS)
!  l = 0
!  do j=1, ndimd
!    dj = real(j-ndimd2-1)
!    do i=1, ndimd
!      di = real(i-ndimd2-1)
!      if (di*di+dj*dj>real(ndimd2*ndimd2)) then
!        l = l + 1
!        dv = detdata(i,j)
!        detlev(1,0) = detlev(1,0) + dv
!        detlev(2,0) = detlev(2,0) + dv*dv
!      end if
!    end do
!  end do
!  if (l>0) then
!    detlev(1,0) = detlev(1,0) / real(l)
!    detlev(2,0) = sqrt( abs( detlev(2,0) / real(l) &
!     & - detlev(1,0)*detlev(1,0) ) ) 
!  else
!    nerr = 12
!    call MSP_ERROR("Failed to detect zero level at corners.", &
!      & subnum+nerr)
!    cycle
!  end if
!  
!  
!  ! loop through wave calculation array in fourier space
!  do j=1, ndimx
!    MSP_detcols(2,j,k) = ndimy ! set left detector column to maximum
!    wx = itowx*MS_TABBED_SCR(j)
!    ! translate to 
!    do i=1, ndimy
!      wy = itowy*MS_TABBED_SCR2(i)
!      !
!    end do ! loop i
!  end do ! loop j
!
!  if (MSP_detimg_output==1) then ! export the detector function to binary file
!    if (.not.allocated(detfunc)) then
!      allocate(detfunc(MS_dimx, MS_dimy), stat=nalloc)
!      if (nalloc/=0) then
!        nerr = 5
!        call MSP_ERROR("Allocation of detector function array failed.",subnum+nerr)
!        goto 667
!      end if
!    end if
!    ! create the detector function from prepared arrays, but in non-transposed and non-scrambled order
!    do j=1, ndimy
!      j1 = MS_TABBED_SCR2(j) + ndim2y+1
!      do i=1, ndimx
!        i1 = MS_TABBED_SCR(i) + ndim2x+1
!        detfunc(i,j) = MSP_detarea(j1,i1,k)
!      end do 
!    end do
!    
!    m = LEN_TRIM(MSP_outfile) ! get length of the standard output file
!    l = INDEX(trim(MSP_outfile),".",back=.TRUE.) ! search for extension point in given output file
!    if (MSP_usedetdef/=0) then ! update output file name with detetctor name
!      if (l<1) then ! no extension wanted
!        write(unit=stmp,fmt='(A)') trim(MSP_outfile)//"_detector_"//trim(MSP_detname(k))
!      else ! last extension starts at position l, insert the index string before
!        write(unit=stmp,fmt='(A)') MSP_outfile(1:l-1)//"_detector_"//trim(MSP_detname(k))//MSP_outfile(l:m)
!      end if
!    else ! keep current output file name
!      if (l<1) then ! no extension wanted
!        write(unit=stmp,fmt='(A)') trim(MSP_outfile)//"_detector"
!      else ! last extension starts at position l, insert the index string before
!        write(unit=stmp,fmt='(A)') MSP_outfile(1:l-1)//"_detector"//MSP_outfile(l:m)
!      end if
!    end if
!    
!    write(unit=MSP_stmp,fmt='(A,I3.3,A)') "Saving detector function #",k," to file ["//trim(stmp)//"]."
!    call PostMessage(trim(MSP_stmp))
!    write(unit=MSP_stmp,fmt='(A,I4,A,I4)') "- detector function: 32-bit float, size: ",MS_dimx," x ",MS_dimy
!    call PostDebugMessage(trim(MSP_stmp))
!    call SaveDataR4(trim(stmp),detfunc,MS_dimx*MS_dimy,nerr)
!    if (nerr/=0) then
!      call MSP_ERROR("Failed to save detector function.",subnum+6)
!    end if
!    
!    
!    deallocate(detfunc, stat=nalloc)
!    
!667 end if ! (MSP_detimg_output==1)
!  
!  
!  end do ! loop k for detectors
!! ------------
!
!
!! ------------
!  if (allocated(detdata)) deallocate(detdata, stat=nalloc)
!  if (allocated(detused)) deallocate(detused, stat=nalloc)
!  if (allocated(detlev)) deallocate(detlev, stat=nalloc)
!! ------------
!
!
!! ------------
!!  write(unit=*,fmt=*) " > MSP_SetCalibratedDetectors: EXIT."
!  return
!
!END SUBROUTINE MSP_SetCalibratedDetectors
!!**********************************************************************!




!* << CALCULATIONS & PROCESSING
!* <<
!* <<
!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!








!**********************************************************************!
!********************* SPECIAL STD ROUTINES ***************************!
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_ERROR(sTxt,nErr)
! function: print error message
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
!            INTEGER*4 :: nErr : error code, displayed when /= 0
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  character*(*), intent(in) :: sTxt
  integer*4, intent(in) :: nErr
  character(len=MSP_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = MSP_el
! ----------

! ----------
! messaging (to screen)
  write(unit=MSP_stdout,fmt='(A)')   ""
  write(unit=MSP_stdout,fmt='(A,I)') "ERROR : MS-Parameter : "//trim(sTxt)//" Error code:",nErr
  write(unit=MSP_stdout,fmt='(A)')   ""
  MSP_err_num = MSP_err_num + 1
  call MSP_HALT()
! ----------
  
  return

END SUBROUTINE MSP_ERROR
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_WARN(sTxt)
! function: print warning message
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
!            INTEGER*4 :: nErr : error code, displayed when /= 0
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  character*(*), intent(in) :: sTxt
  character(len=MSP_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = MSP_el
! ----------

! ----------
! messaging (to screen)
  write(unit=MSP_stdout,fmt='(A)') "Warning : MS-Parameter : "//trim(sTxt)
  MSP_warn_num = MSP_warn_num + 1
! ----------
  
  return

END SUBROUTINE MSP_WARN
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_HALT()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_HALT: INIT."
! ------------


! ------------
  write(unit=MSP_stdout,fmt='(A,I)') "Errors  :",MSP_err_num
  write(unit=MSP_stdout,fmt='(A,I)') "Warnings:",MSP_warn_num
  write(unit=MSP_stdout,fmt='(A)')   "Halting program."
  write(unit=MSP_stdout,fmt='(A)')   ""
  stop
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_HALT: EXIT."
  return

END SUBROUTINE MSP_HALT
!**********************************************************************!


 
!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_INITCL()
! function: Initializes the clock timer counts.
! -------------------------------------------------------------------- !
! parameter: none
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1800
!
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_INITCL: INIT."
! ------------


! ------------
  if (MSP_runtimes==1) then ! init the run-time clock state
    call system_clock(MSP_clinit, MSP_clrate, MSP_clmax)
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_INITCL: EXIT."
  return

END SUBROUTINE MSP_INITCL
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSP_GETCLS(cls)
! function: Determines the time span to the previous clock init in
!           seconds.
! -------------------------------------------------------------------- !
! parameter: cls : real*4 : time span in seconds
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1900
  integer*8 :: clnow, cldelta
  real*4, intent(out) :: cls
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > MSP_GETCLS: INIT."
  cls = 0.0
! ------------


! ------------
  if (MSP_runtimes==1) then 
    call system_clock(clnow)
    if (clnow<MSP_clinit) then
      cldelta = MSP_clmax - MSP_clinit + clnow
      MSP_clinit = clnow
    else
      cldelta = clnow - MSP_clinit
    end if
    cls = real(cldelta,kind=4)/MSP_clrate
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > MSP_GETCLS: EXIT."
  return

END SUBROUTINE MSP_GETCLS
!**********************************************************************!




END MODULE MSAparams
!**********************************************************************!
!**********************************************************************!
!**********************************************************************!









!**********************************************************************!
!************************ EXTERNAL ROUTINES ***************************!
!**********************************************************************!








!**********************************************************************!
!*********************SUBROUTINE COMMENTS TEMPLATE*********************!
!**********************************************************************!

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
!  integer*4, parameter :: subnum = 2200
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

