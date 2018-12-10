!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  emsdata.f90                                           !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2009-2018   !
!    Version  :  1.0.0, July 10, 2009                                  !
!    Version  :  2.0.0, July 11, 2012                                  !
!    Version  :  2.1.0, December 16, 2014                              !
!                                                                      !
!                                                                      !
!**********************************************************************!
!
!  Purpose: Implementation of routines to handle EMS type of data files
!           used by the program MSA for I/O of phase gratings and
!           structure data.
!
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
!**********************************************************************!
INTEGER*4 FUNCTION SwapInt4(nin)
! function: returns byteswapped input of type integer*4
! -------------------------------------------------------------------- !
! parameter: nin: integer*4: nin
!            
! -------------------------------------------------------------------- !

  implicit none

  integer*4,intent(in) :: nin
  integer*4 :: dummy, i, j
  character*4 :: ch, chs

  dummy = 0
  ch = transfer(nin,ch)
  do i=1,4
    j=5-i
    chs(i:i) = ch(j:j)
  end do
  dummy = transfer(chs,dummy)
  SwapInt4 = dummy

END FUNCTION SwapInt4
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
REAL*4 FUNCTION SwapReal4(rin)
! function: returns byteswapped input of type real*4
! -------------------------------------------------------------------- !
! parameter: rin: real*4: nin
!            
! -------------------------------------------------------------------- !

  implicit none

  real*4,intent(in) :: rin
  real*4 :: dummy
  integer*4 :: i, j
  character*4 :: ch, chs

  dummy = 0.0
  ch = transfer(rin,ch)
  do i=1,4
    j=5-i
    chs(i:i) = ch(j:j)
  end do
  dummy = transfer(chs,dummy)
  SwapReal4 = dummy

END FUNCTION SwapReal4
!**********************************************************************!








!**********************************************************************!
!                                                                      !
!   Purpose: MODULE helping to handle ems data                         !
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
MODULE EMSdata

!  use binio
  USE IFPORT

  IMPLICIT NONE

! declare internal data types

! accessibility of subroutines or functions
!  private :: 
  private :: EMS_ERROR, EMS_WARN
  private :: EMS_SLI_PRMREAD_INT32
  private :: EMS_SLI_PRMWRITE_INT32
  private :: EMS_SLI_PRMREAD_REAL32
  private :: EMS_SLI_PRMWRITE_REAL32
  private :: EMS_SLI_PRMREAD_STRING
  private :: EMS_SLI_PRMWRITE_STRING
  private :: EMS_SLI_ALLOC_ITAB
  private :: EMS_SLI_DEALLOC_ITAB
  private :: EMS_SLI_LOAD_ITAB
  private :: EMS_SLI_SAVE_ITAB
  
!  public :: 
  public :: EMS_INIT, EMS_UNINIT, EMS_GetFreeLFU
  public :: EMS_SLI_loaddata
  public :: EMS_SLI_loadparams
  public :: EMS_SLI_save
  public :: EMS_SLI_settitle
  public :: EMS_SLI_settab0
  public :: EMS_SLI_settab1
  
! declare module global (but PRIVATE) variables, params and arrays

!   length of names and internal strings
  integer*4, public, parameter :: EMS_ll = 1024

!   empty line placeholder (emptied at startup)
  character(len=EMS_ll), private :: EMS_el
  
!   logical unit limits
  integer*4, public, parameter :: EMS_lfumin = 50
  integer*4, public, parameter :: EMS_lfumax = 100
  
! declare module global (but PUBLIC) variables, params and arrays
  integer*4, private, parameter :: EMS_SLI_VERSION          = 2012071101
  integer*4, private, parameter :: EMS_SLI_IFORM            = 1 ! expected data format = complex 8-byte (complex*8)

  integer*4, private, parameter :: EMS_SLI_POS_DATASIZE_X   = Z'00000000'
  integer*4, private, parameter :: EMS_SLI_POS_DATASIZE_Y   = Z'00000004'
  integer*4, private, parameter :: EMS_SLI_POS_IFORM        = Z'00000010'
  integer*4, private, parameter :: EMS_SLI_POS_TITLE        = Z'0000004C'
  integer*4, private, parameter :: EMS_SLI_POS_EXTHEADER_VER= Z'00000074'
  integer*4, private, parameter :: EMS_SLI_POS_VARIANT_NUM  = Z'00000078'
  integer*4, private, parameter :: EMS_SLI_POS_CONTENT_TYPE = Z'00000080'
  integer*4, private, parameter :: EMS_SLI_POS_SUBSLC_NUM   = Z'00000084'
  integer*4, private, parameter :: EMS_SLI_POS_ALT_OFFSET   = Z'00000088' ! alternative data offset,
                                                                          ! the data table starts with the elastic phase gratings
                                                                          ! after the elastic phase gratings, inelastic transition
                                                                          ! probabilities are saved in logical sequence
                                                                          ! according to the inelastic data table definitions
  integer*4, private, parameter :: EMS_SLI_POS_ITAB_OFFSET  = Z'0000008C' ! offset of the inelastic data table (=1100 by default)
  integer*4, private, parameter :: EMS_SLI_POS_THICKNESS    = Z'00001034'
  integer*4, private, parameter :: EMS_SLI_POS_HIGHTENSION  = Z'00001040'
  integer*4, private, parameter :: EMS_SLI_POS_PHYSSIZE_X   = Z'00001044'
  integer*4, private, parameter :: EMS_SLI_POS_PHYSSIZE_Y   = Z'00001048'
  integer*4, private, parameter :: EMS_SLI_POS_SUBSLC_THICK = Z'00001060'
  integer*4, private, parameter :: EMS_SLI_POS_DATA         = Z'00002000' ! default data offset, may vary in case of inelastic potentials
  
  integer*4, private, parameter :: EMS_SLI_TITLE_SIZE       = 40
  integer*4, private, parameter :: EMS_SLI_ITAB_OFFSET_DEFAULT = Z'00001100' ! = 4352 0x1100, default inelastic data table offset
  
!

!   PI
  real*4, public :: EMS_pi
!   conversion factor degree to radian
  real*4, public :: EMS_d2r
!   conversion factor radian to degree
  real*4, public :: EMS_r2d
  
  integer*4, private :: EMS_SLI_data_swap ! swap flag
  DATA EMS_SLI_data_swap /0/
  
  integer*4, public :: EMS_SLI_data_ver ! file extended header version
  integer*4, public :: EMS_SLI_data_nvar ! file number of variants
  integer*4, public :: EMS_SLI_data_dimx, EMS_SLI_data_dimy ! file data dimension
  integer*4, public :: EMS_SLI_data_iform ! file data format
  integer*4, public :: EMS_SLI_data_ctype ! file data content type (0=phase grating, 1=potential)
  DATA EMS_SLI_data_ctype /0/ ! default is phase grating
  integer*4, public :: EMS_SLI_data_sslnum ! file data sub-slice number
  integer*4, public :: EMS_SLI_data_alt_offset ! alternative data offset (overrides the default)
  DATA EMS_SLI_data_alt_offset /EMS_SLI_POS_DATA/
  integer*4, public :: EMS_SLI_data_itab_offset ! inelastic table data offset
  DATA EMS_SLI_data_itab_offset /EMS_SLI_ITAB_OFFSET_DEFAULT/
  character(len=EMS_SLI_TITLE_SIZE), public :: EMS_SLI_data_title
  real*4, public :: EMS_SLI_data_ht ! high-tension value
  real*4, public :: EMS_SLI_data_szz ! thickness of slice
  real*4, public :: EMS_SLI_data_szx ! slice size horizontal
  real*4, public :: EMS_SLI_data_szy ! slice size vertical
  real*4, public :: EMS_SLI_data_sslz ! sub-slice thickness
  
!
! Variables holding inelastic data tables
!
! the table is allocatable, check the allocation before reading and writing data from/to the tables
! * allocation is done automatically while reading data from a slice file, previous data will be lost
! * allocation must be done from outsides when preparing data for writing to a slice file
!
! USE THE FOLLOWING THREE VARIABLES TO DEFINE OR DETERMINE THE ARRAY SIZES
! - define from outsides when writing slice files (before calling EMS_SLI_save)
! - determine from outsides when reading slice files (after calling EMS_SLI_loadparams and before calling EMS_SLI_loaddata_inel)
!
  integer*4, public :: EMS_SLI_data_niaty ! number of atom types in the inelastic data table
  DATA EMS_SLI_data_niaty /0/
  integer*4, public :: EMS_SLI_data_niatr_max ! max number of atomic transitions
  DATA EMS_SLI_data_niatr_max /0/
  integer*4, public :: EMS_SLI_data_niatp_max ! max number of atomic coordinates
  DATA EMS_SLI_data_niatr_max /0/
!
! DETERMINATION OF THE DATA STRUCTURE SIZE OR DO WE NEED TO MOVE THE HEADER SOMEWHERE ELSE?
! - expected size = 4*(3 + EMS_SLI_data_niaty*( 4 + EMS_SLI_data_niatr_max + 2*EMS_SLI_data_niatp_max )) Bytes
! - available space at the default table position: 0x2000 - 0x1100 = 0x0F00 = 3840 Bytes
!
  integer*4, public, allocatable :: EMS_SLI_data_attyz(:) ! atomic numbers for each atom type (not necessarily different)
                                                          ! allocation size: (EMS_SLI_data_niaty)
  integer*4, public, parameter :: EMS_SLI_data_natdat = 3 ! number of atom data per type: charge, Biso, Occ
  real*4, public, allocatable    :: EMS_SLI_data_atdat(:,:) ! data for each atom type (not necessarily different, ionic charge, Biso, Occ)
                                                          ! allocation size: (EMS_SLI_data_niaty)
  integer*4, public, allocatable :: EMS_SLI_data_niatr(:) ! number of atom transitions/potentials
                                                          ! allocation size: (EMS_SLI_data_niaty)
  integer*4, public, allocatable :: EMS_SLI_data_iatrcd(:,:) ! transistion codes
                                                          ! allocation size: (EMS_SLI_data_niatr_max, EMS_SLI_data_niaty)
  real*4, public, allocatable :: EMS_SLI_data_iatren(:,:) ! transistion energies
                                                          ! allocation size: (EMS_SLI_data_niatr_max, EMS_SLI_data_niaty)
  integer*4, public, allocatable :: EMS_SLI_data_niatp(:) ! number of atoms of this type = number of proj. atom positions
                                                          ! allocation size: (EMS_SLI_data_niaty) 
  real*4, public, allocatable :: EMS_SLI_data_iaatpo(:,:,:) ! fract. slice atom coordinates (equilibrium positions)
                                                          ! allocation size: (3,EMS_SLI_data_niatp_max, EMS_SLI_data_niaty)


  CONTAINS
  
  
  
!**********************************************************************!
!
! Information on the EMS image and slice files
! ----------------------------------------
!
! Structure information was interpreted from the fortran source code
! of EMS by P. Stadelmann, files 'openfi.f' and 'openfi.com'.
!
! The file consists of a header block and a data block.
! The header block is 8192 bytes long.
! The data block consists of nsam*nrow data items, where the
! item size depends on the data type.
!
! Header structure:
!    OFFSET(DEC,HEX) TYPE            SIZE        INFO
!    0      0000     integer*4       4           nsam, image dims, number of samples
!    4      0004     integer*4       4           nrow, image dims, number of rows
!    8      0008     integer*4       4           nsaeff, effective x dim of the image
!    12     000C     integer*4       4           nroeff, effective y dim of the image
!    16     0010     integer*4       4           iform, data format, -1=='F' FT=complex_64, 0=='R' real_32, 1=='C' complex_64, 2=='I' int_8
!    20     0014     integer*4       4           iprot, ? flag for title interpretation ?
!    24     0018     integer*4       4           imami, flags that min. and max. values have been calculated already (0=no, 1=yes)
!    28     001C     integer*4       4           igeom, flags that geometry correction was done (0=no, 1=yes)
!    32     0020     integer*4       4           nslice, number of multislice iteration
!    36     0024     integer*4       4           iu, first zone index
!    40     0028     integer*4       4           iv, second zone index
!    44     002C     integer*4       4           iw, third zone index
!    48     0030     integer*4       4           ispcel, super-cell flag (0: unit cell, 1: super-cell)
!    52     0034     integer*4       24          itr(3,2), 2-d to 3-d indices transform matrix
!    56     0038     integer*4       ||          itr(3,2) (2,1)
!    60     003C     integer*4       ||          itr(3,2) (3,1)
!    64     0040     integer*4       ||          itr(3,2) (1,2)
!    68     0044     integer*4       ||          itr(3,2) (2,2)
!    72     0048     integer*4       ||          itr(3,2) (3,2)
!    76     004C     character*40    40          imatit, title string of image / slice
!    ...    ...      ...             ||
!*** 116    0074     integer*4       4           extended header version (added by J.B. 101202)
!*** 120	0078     integer*4       4           number of data variants
!*** 128    0080     integer*4       4           content type (0=phase grating, 1=potential) (added by J.B. 110407)
!*** 132    0084     integer*4       4           number of sub-slices (>1 if sub-slicing is done) (added by J.B. 110407)
!*** 136    0088     integer*4       4           alternative data offset, overrides default data offset (added by J.B. 120711)
!*** 140    008C     integer*4       4           offset of the inelastic data table [see below for table structure definition] (added by J.B. 120711)
!    ...
!    4096   1000     real*4          4           beamh, 1st beam trace index (3-d)
!    4100   1004     real*4          4           beamk, 2nd beam trace index (3-d)
!    4104   1008     real*4          4           beaml, 3rd beam trace index (3-d)
!    4108   100C     real*4          4           picmin, min. of the image (real part)
!    4112   1010     real*4          4           picmax, max. of the image (real part)
!    4116   1014     real*4          4           commin, min. of the image (imag part)
!    4120   1018     real*4          4           commax, max. of the image (imag part)
!    4124   101C     real*4          4           a, 1st lattice side [nm]
!    4128   1020     real*4          4           b, 2nd lattice side [nm]
!    4132   1024     real*4          4           c, 3rd lattice side [nm]
!    4136   1028     real*4          4           alfa, (b,c) lattice angle [deg]
!    4140   102C     real*4          4           beta, (a,c) lattice angle [deg]
!    4144   1030     real*4          4           gama, (a,b) lattice angle [deg]
!    4148   1034     real*4          4           thick, slice thickness [nm]
!    4152   1038     real*4          4           absor, absorption coefficiant [a.u.]
!    4156   103C     real*4          4           sigma, interaction constant
!    4160   1040     real*4          4           voltag, accelerating voltage  [kv]
!    4164   1044     real*4          4           ca, image geometry in real*4 space [a.u.]
!    4168   1048     real*4          4           cb, image geometry in real*4 space [a.u.]
!    4172   104C     real*4          4           angle, image geometry in real*4 space [a.u.]
!    4176   1050     real*4          4           b2dh, 1st beam trace index (2-d)
!    4180   1054     real*4          4           b2dk, 2nd beam trace index (2-d)
!    4184   1056
!    4188   105C
!    4192   1060     real*4          4           sub-slice thickness [nm] (added by J.B. 110407)
!
!    !!!! Atom & Inelastic potential table (max size to default data offset:
!*** 4352   1100     integer*4       4           naty, number of atom types (added by J.B.) 120711)
!*** 4360   1108     integer*4       4           natpm, max. number of atom positions (for pre-allocation of arrays)
!*** 4356   1104     integer*4       4           natrm, max. number of atom transitions (for pre-allocation of arrays)
!*** 4364   110C     integer*4       4    i-+    attyz(i), atomic number of the current atom type
!*** 4368   1110     3*real*4        12     |    atdat(1:3,i), data of the current atom type
!*** 4380   1118     integer*4       4      |    niatp(i), max. number of atom positions for each atom type in the slice
!*** 4384   111C     3*real*4        12     +-k  iaatpo(1:3,k,i), fractional atom coordinates x, y, z [ca, cb, cz]
!*** 4384+12*niatp   integer*4       4      |    niatr(i), number of atom transitions for each atom type [= number of potentials for the atom type]
!*** 4388+12*niatp   integer*4       4      +-j  iatrcd(j,i), inelastic transistion codes [= 100*shell index + l-index, K -> 100, L2,3-> 223, M4,5 -> 345]
!*** 4388+4*niatr+12*niatp  real*4   4      +-j  iatren(j,i), energy losses [eV]
!
!    TABLE STRUCTURE IN CLEAR WRITING:
!
!    1) number of atom types
!    2) max. number of positions
!    3) max. number of transitions
!    4) LOOP atom types
!    4.1) atomic number Z
!    4.2) atom data (ionic charge, Biso, Occupancy)
!    4.3) number of atom positions
!    4.4) LOOP atom positions
!    4.4.1) atom position (ca, cb, cdz)
!    4.5) number of atom state transitions
!    4.5.1) atom transition code
!    4.5.2) atom transition energy
!
!*** 4392+8*niatr+12*niatp ... next attyz(i+1)
!
! Data offset (default, may change):
!    8192   2000
!
!**********************************************************************!




!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!
!* >>
!* >>
!* >> SETUP & INIT ROUTINES


!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_INIT()
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
!  write(unit=*,fmt=*) " > EMS_INIT: INIT."
! ------------


! ------------
! empty a line
  EMS_el = REPEAT(" ",EMS_ll)
! empty all strings, so be careful with initializing
! ------------


! ------------
  EMS_pi = atan(1.0)*4.0
  EMS_d2r = EMS_pi / 180.0
  EMS_r2d = 1.0 / EMS_d2r
! ------------

! ------------
  EMS_SLI_data_dimx = 0
  EMS_SLI_data_dimy = 0
! ------------

  EMS_SLI_data_alt_offset = EMS_SLI_POS_DATA ! default elastic data offset
  EMS_SLI_data_itab_offset = EMS_SLI_ITAB_OFFSET_DEFAULT ! default header table offset

! ------------
!  write(unit=*,fmt=*) " > EMS_INIT: EXIT."
  return

END SUBROUTINE EMS_INIT
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_UNINIT()
! function: uninits the module
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 200
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_UNINIT: INIT."
  integer*4 :: lerr
! ------------

! ------------
  call EMS_SLI_DEALLOC_ITAB(lerr)
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_UNINIT: EXIT."
  return

END SUBROUTINE EMS_UNINIT
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
SUBROUTINE EMS_GetFreeLFU(lfu)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 500
  integer*4, intent(inout) :: lfu
  logical :: isopen
  integer*4 :: u0, u1
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > GetFreeLFU: INIT."
  lfu = EMS_lfumin
  u0 = EMS_lfumin
  u1 = max(EMS_lfumax,EMS_lfumin)
! ------------

! ------------
  
  do
    INQUIRE(unit=lfu,opened=isopen)
    if (.not.isopen) exit
    lfu = lfu + 1
!   catch to many open trials
    if (lfu>u1) then
      call EMS_ERROR("Failed to acquire logical file unit.",subnum+1)
    end if
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > GetFreeLFU: EXIT."
  return

END SUBROUTINE EMS_GetFreeLFU
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_settitle(stitle)
! function: sets the internal title string
! -------------------------------------------------------------------- !
! parameter: character(len=*) :: stitle     ! new title
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 600
  character(len=*), intent(in) :: stitle
  integer*4 :: outlen
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_settitle: INIT."
! preset the characters to zero
  EMS_SLI_data_title = repeat(char(0),EMS_SLI_TITLE_SIZE)
! determine length
  outlen = min(40,len(stitle))
! ------------

! ------------
! EMS_SLI_data_title
  EMS_SLI_data_title(1:outlen) = stitle(1:outlen)
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_settitle: EXIT."
  return

END SUBROUTINE EMS_SLI_settitle
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_PRMREAD_INT32(lfu, npos, prm, nerr)
! function: reads a 32-bit integer parameter from specified position
!           of an open file unit
! -------------------------------------------------------------------- !
! parameter:
! INPUT:
!  integer*4 :: lfu         ! logical file unit number
!  integer*4 :: npos        ! io-offset from begin of the file
! OUTPUT:
!  integer*4 :: prm         ! the parameter variable reference
!  integer*4 :: nerr        ! error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 700
  integer*4, intent(in) :: lfu, npos
  integer*4, intent(inout) :: prm, nerr
  integer*4, external :: SwapInt4
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_PRMREAD_INT32: INIT."
  nerr = 0
! ------------

! ------------
  nerr = fseek(lfu, npos, 0) ! set file pointer to requested position
  if (nerr/=0) then
    call EMS_ERROR("Sequential positioning failed.", subnum+1)
    return
  end if
  read(unit=lfu,iostat=nerr) prm ! read parameter data
  if (nerr/=0) then
    call EMS_ERROR("Failed to read data.", subnum+2)
    return
  end if
  if (1==EMS_SLI_data_swap) prm = SwapInt4(prm) ! handle byte-swap
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_PRMREAD_INT32: EXIT."
  return

END SUBROUTINE EMS_SLI_PRMREAD_INT32
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_PRMWRITE_INT32(lfu, npos, prm, nerr)
! function: writes a 32-bit integer parameter at specified position
!           to an open file unit
! -------------------------------------------------------------------- !
! parameter:
! INPUT:
!  integer*4 :: lfu         ! logical file unit number
!  integer*4 :: npos        ! io-offset from begin of the file
!  integer*4 :: prm         ! the parameter variable reference
! OUTPUT:
!  integer*4 :: nerr        ! error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  integer*4, intent(in) :: lfu, npos, prm
  integer*4, intent(inout) :: nerr
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_PRMWRITE_INT32: INIT."
  nerr = 0
! ------------

! ------------
  nerr = fseek(lfu, npos, 0) ! position
  if (nerr/=0) then
    call EMS_ERROR("Sequential positioning failed.", subnum+1)
    return
  end if
  write(unit=lfu,iostat=nerr) prm
  if  (nerr/=0) then
    call EMS_ERROR("Failed to write data to file.", subnum+2)
    return
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_PRMWRITE_INT32: EXIT."
  return

END SUBROUTINE EMS_SLI_PRMWRITE_INT32
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_PRMREAD_REAL32(lfu, npos, prm, nerr)
! function: reads a 32-bit float parameter from specified position
!           of an open file unit
! -------------------------------------------------------------------- !
! parameter: 
! INPUT:
!  integer*4 :: lfu         ! logical file unit number
!  integer*4 :: npos        ! io-offset from begin of the file
! OUTPUT:
!  real*4 :: prm            ! the parameter variable reference
!  integer*4 :: nerr        ! error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 900
  integer*4, intent(in) :: lfu, npos
  integer*4, intent(inout) :: nerr
  real*4, intent(inout) :: prm
  real*4, external :: SwapReal4
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_PRMREAD_REAL32: INIT."
  nerr = 0
! ------------

! ------------
  nerr = fseek(lfu, npos, 0) ! set file pointer to requested position
  if (nerr/=0) then
    call EMS_ERROR("Sequential positioning failed.", subnum+1)
    return
  end if
  read(unit=lfu,iostat=nerr) prm ! read parameter data
  if (nerr/=0) then
    call EMS_ERROR("Failed to read data.", subnum+2)
    return
  end if
  if (1==EMS_SLI_data_swap) prm = SwapReal4(prm) ! handle byte-swap
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_PRMREAD_REAL32: EXIT."
  return

END SUBROUTINE EMS_SLI_PRMREAD_REAL32
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_PRMWRITE_REAL32(lfu, npos, prm, nerr)
! function: writes a 32-bit float parameter to specified position
!           of an open file unit
! -------------------------------------------------------------------- !
! parameter:
! INPUT:
!  integer*4 :: lfu         ! logical file unit number
!  integer*4 :: npos        ! io-offset from begin of the file
!  real*4 :: prm            ! the parameter variable reference
! OUTPUT:
!  integer*4 :: nerr        ! error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1000
  integer*4, intent(in) :: lfu, npos
  integer*4, intent(inout) :: nerr
  real*4, intent(in) :: prm
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_PRMWRITE_REAL32: INIT."
  nerr = 0 
! ------------


! ------------
  nerr = fseek(lfu, npos, 0) ! position
  if (nerr/=0) then
    call EMS_ERROR("Sequential positioning failed.", subnum+1)
    return
  end if
  write(unit=lfu,iostat=nerr) prm
  if  (nerr/=0) then
    call EMS_ERROR("Failed to write data to file.", subnum+2)
    return
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_PRMWRITE_REAL32: EXIT."
  return

END SUBROUTINE EMS_SLI_PRMWRITE_REAL32
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_PRMREAD_STRING(lfu, npos, prm, nerr)
! function: reads a character string from specified position
!           of an open file unit
!           the full length of prm is read in
! -------------------------------------------------------------------- !
! parameter: 
! INPUT:
!  integer*4 :: lfu         ! logical file unit number
!  integer*4 :: npos        ! io-offset from begin of the file
! OUTPUT:
!  character(len=*) :: prm  ! the parameter variable reference
!  integer*4 :: nerr        ! error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1100
  integer*4, intent(in) :: lfu, npos
  integer*4, intent(inout) :: nerr
  character(len=*), intent(inout) :: prm
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_PRMREAD_STRING: INIT."
  nerr = 0
! ------------

! ------------
  nerr = fseek(lfu, npos, 0) ! set file pointer to requested position
  if (nerr/=0) then
    call EMS_ERROR("Sequential positioning failed.", subnum+1)
    return
  end if
  read(unit=lfu,iostat=nerr) prm ! read parameter data
  if (nerr/=0) then
    call EMS_ERROR("Failed to read data.", subnum+2)
    return
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_PRMREAD_STRING: EXIT."
  return

END SUBROUTINE EMS_SLI_PRMREAD_STRING
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_PRMWRITE_STRING(lfu, npos, prm, nerr)
! function: writes a character string at specified position
!           to an open file unit
!           the full length of prm is written
! -------------------------------------------------------------------- !
! parameter:
! INPUT:
!  integer*4 :: lfu         ! logical file unit number
!  integer*4 :: npos        ! io-offset from begin of the file
!  character(len=*) :: prm  ! the parameter variable reference
! OUTPUT:
!  integer*4 :: nerr        ! error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1200
  integer*4, intent(in) :: lfu, npos
  integer*4, intent(inout) :: nerr
  character(len=*), intent(in) :: prm
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_PRMWRITE_STRING: INIT."
  nerr = 0
! ------------

! ------------
  nerr = fseek(lfu, npos, 0) ! position
  if (nerr/=0) then
    call EMS_ERROR("Sequential positioning failed.", subnum+1)
    return
  end if
  write(unit=lfu,iostat=nerr) prm
  if  (nerr/=0) then
    call EMS_ERROR("Failed to write data to file.", subnum+2)
    return
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_PRMWRITE_STRING: EXIT."
  return

END SUBROUTINE EMS_SLI_PRMWRITE_STRING
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_ALLOC_ITAB(naty, natr, natp, nerr)
! function: Allocates table memory for the position and inelastic
!           transition data.
! -------------------------------------------------------------------- !
! parameter: 
!  integer*4, intent(in) :: naty    = # individual atom types
!  integer*4, intent(in) :: natr    = # considered transitions
!  integer*4, intent(in) :: natp    = # backup atom positions in slice
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1400
  integer*4, intent(in) :: naty
  integer*4, intent(in) :: natr
  integer*4, intent(in) :: natp
  integer*4, intent(inout) :: nerr
  integer*4 :: lerr, nfail, nadat
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_ALLOC_ITAB: INIT."
  lerr = 0
  nerr = 0
! deallocate the arrays
  call EMS_SLI_DEALLOC_ITAB(lerr)
  if (lerr/=0) then
    nfail = 1
    goto 13
  end if
! ------------

! ------------
! Reallocate the arrays
  if (naty>0) then ! there are atoms in the slice
  
    ! allocate the primary lists depending on the atom type number only
    allocate(EMS_SLI_data_attyz(naty),stat=nerr) ! ATOM TYPES (Z numbers)
    if (nerr/=0) then
      nfail = 2
      goto 14
    end if
    EMS_SLI_data_attyz = 0
    nadat = EMS_SLI_data_natdat
    allocate(EMS_SLI_data_atdat(nadat,naty),stat=nerr) ! ATOM DATA
    if (nerr/=0) then
      nfail = 3
      goto 14
    end if
    EMS_SLI_data_atdat = 0.0
    allocate(EMS_SLI_data_niatr(naty),stat=nerr) ! NUMBERS OF TRANSITIONS
    if (nerr/=0) then
      nfail = 4
      goto 14
    end if
    EMS_SLI_data_niatr = 0
    allocate(EMS_SLI_data_niatp(naty),stat=nerr) ! NUMBERS OF POSITIONS
    if (nerr/=0) then
      nfail = 7
      goto 14
    end if
    EMS_SLI_data_niatp = 0
    
    if (natr>0) then
      ! allocate arrays depending on the number of transitions
      allocate(EMS_SLI_data_iatrcd(natr,naty),stat=nerr) ! TRANSITION CODES
      if (nerr/=0) then
        nfail = 5
        goto 14
      end if
      EMS_SLI_data_iatrcd = 0
      allocate(EMS_SLI_data_iatren(natr,naty),stat=nerr) ! TRANSITION ENERGIES
      if (nerr/=0) then
        nfail = 6
        goto 14
      end if
      EMS_SLI_data_iatren = 0.0
    end if
    
    if (natp>0) then
      ! allocate arrays depending on the number of positions
      allocate(EMS_SLI_data_iaatpo(3,natp,naty),stat=nerr)
      if (nerr/=0) then
        nfail = 8
        goto 14
      end if
      EMS_SLI_data_iaatpo = 0.0
    end if
    
  else ! no atoms in slice, set all to zero and return
  
    EMS_SLI_data_niaty = 0
    EMS_SLI_data_niatr_max = 0
    EMS_SLI_data_niatp_max = 0
    return
    
  end if
! ------------


! ------------
! update the size members
  EMS_SLI_data_niaty = max(0,naty)
  EMS_SLI_data_niatr_max = max(0,natr)
  EMS_SLI_data_niatp_max = max(0,natp)
! ------------


! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_ALLOC_ITAB: EXIT."
  return
  
13 nerr = 1
  call EMS_ERROR("Failed to deallocate array memory.", subnum+nfail)
  return
14 nerr = 2
  call EMS_ERROR("Failed to allocate array memory.", subnum+nfail)
  return
  

END SUBROUTINE EMS_SLI_ALLOC_ITAB
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_DEALLOC_ITAB(nerr)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1500
  integer*4, intent(inout) :: nerr
  integer*4 :: nfail
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_DEALLOC_ITAB: INIT."
  nerr = 0
! ------------


! ------------
! De-Allocate the tables
  if (allocated(EMS_SLI_data_attyz)) deallocate(EMS_SLI_data_attyz,stat=nerr)
  if (nerr/=0) then
    nfail = 1
    goto 13
  end if
  if (allocated(EMS_SLI_data_atdat)) deallocate(EMS_SLI_data_atdat,stat=nerr)
  if (nerr/=0) then
    nfail = 2
    goto 13
  end if
  if (allocated(EMS_SLI_data_niatr)) deallocate(EMS_SLI_data_niatr,stat=nerr)
  if (nerr/=0) then
    nfail = 3
    goto 13
  end if
  if (allocated(EMS_SLI_data_iatrcd)) deallocate(EMS_SLI_data_iatrcd,stat=nerr)
  if (nerr/=0) then
    nfail = 4
    goto 13
  end if
  if (allocated(EMS_SLI_data_iatren)) deallocate(EMS_SLI_data_iatren,stat=nerr)
  if (nerr/=0) then
    nfail = 5
    goto 13
  end if
  if (allocated(EMS_SLI_data_niatp)) deallocate(EMS_SLI_data_niatp,stat=nerr)
  if (nerr/=0) then
    nfail = 6
    goto 13
  end if
  if (allocated(EMS_SLI_data_iaatpo)) deallocate(EMS_SLI_data_iaatpo,stat=nerr)
  if (nerr/=0) then
    nfail = 7
    goto 13
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_DEALLOC_ITAB: EXIT."
  return
  
13 nerr = 1
  call EMS_ERROR("Failed to deallocate array memory.", subnum+nfail)
  return


END SUBROUTINE EMS_SLI_DEALLOC_ITAB
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_LOAD_ITAB(llfu, nerr)
! function: Loads the inelastic data table to the local variables and
!           arrays. Requires slice version 2012071101 or higher and a
!           file opened for binary reading and connected to a logical
!           file unit llfu.
! -------------------------------------------------------------------- !
! parameter: 
!  integer*4, intent(in) :: llfu    = logical file unit for reading
!  integer*4, intent(inout) :: nerr = error codes (0 = no error)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1300
  integer*4, intent(in) :: llfu
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nfail        ! local error number
  logical :: isopen
  integer*4 :: ncurpos      ! current read offset
  integer*4 :: itab_size    ! total table size in bytes
  integer*4 :: i, j, k      ! iterators
  integer*4 :: naty, natr, natrmax, natp, natpmax, nadat ! number of data
  integer*4 :: itmp1, itmp2 ! temp. variables
  real*4 :: rtmp1, rtmp2    ! temp. variables
  integer*4, external :: SwapInt4 ! 4-byte integer byte swap function
  real*4, external :: SwapReal4 ! 4-byte real byte swap function
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_LOAD_ITAB: INIT."
  nfail = 0
  nadat = EMS_SLI_data_natdat
! ------------


! ------------
! check if the file is open
  inquire(unit=llfu, OPENED=isopen)
  if (.not.isopen) then
    nfail = 1
    goto 12
  end if
! ------------


! ------------
! read the table offset from the primary header
  EMS_SLI_data_itab_offset = EMS_SLI_ITAB_OFFSET_DEFAULT
  call EMS_SLI_PRMREAD_INT32(llfu, EMS_SLI_POS_ITAB_OFFSET, EMS_SLI_data_itab_offset, nerr)
  if (nerr/=0) then
    nfail = 2
    goto 13
  end if
  ncurpos = EMS_SLI_data_itab_offset
  if (ncurpos < EMS_SLI_ITAB_OFFSET_DEFAULT) then
    ! invalid position, do not read the table
    return
  end if
! ------------


! ------------
! read the table size parameters
!  integer*4, public :: EMS_SLI_data_niaty ! number of atom types in the inelastic data table
!  integer*4, public :: EMS_SLI_data_niatr_max ! max number of atomic transitions
!  integer*4, public :: EMS_SLI_data_niatp_max ! max number of atomic coordinates
  call EMS_SLI_PRMREAD_INT32(llfu, ncurpos, EMS_SLI_data_niaty, nerr)
  if (nerr/=0) then
    nfail = 3
    goto 13
  end if
  ncurpos = ncurpos + 4
  call EMS_SLI_PRMREAD_INT32(llfu, ncurpos, EMS_SLI_data_niatp_max, nerr)
  if (nerr/=0) then
    nfail = 4
    goto 13
  end if
  ncurpos = ncurpos + 4
  call EMS_SLI_PRMREAD_INT32(llfu, ncurpos, EMS_SLI_data_niatr_max, nerr)
  if (nerr/=0) then
    nfail = 5
    goto 13
  end if
  !
  ncurpos = ncurpos + 4
  naty = max(0,EMS_SLI_data_niaty)
  natrmax = max(0,EMS_SLI_data_niatr_max)
  natpmax = max(0,EMS_SLI_data_niatp_max)
  
! ... and calculate the maximum table size
  itab_size = 4*(3 + naty*( 3 + nadat + 2*natrmax + 3*natpmax ) )
! ------------


! ------------
! Table reading 
!  integer*4, public, allocatable :: EMS_SLI_data_attyz(:) ! atomic numbers for each atom type (not necessarily different)
!  real*4, public, allocatable    :: EMS_SLI_data_atdat(:,:) ! atom type data (not necessarily different)
!  integer*4, public, allocatable :: EMS_SLI_data_niatr(:) ! number of atom transitions
!  integer*4, public, allocatable :: EMS_SLI_data_iatrcd(:,:) ! transistion codes
!  real*4, public, allocatable :: EMS_SLI_data_iatren(:,:) ! transistion energies
!  integer*4, public, allocatable :: EMS_SLI_data_niatp(:) ! number of atoms of this type = number of proj. atom positions
!  real*4, public, allocatable :: EMS_SLI_data_iaatpo(:,:,:) ! projected fract. atom coordinates
!
! Allocate the tables, the routine call will deallocate existing tables first
  call EMS_SLI_ALLOC_ITAB(naty, natrmax, natpmax, nerr)
  if (nerr/=0) then
    nfail = 14
    goto 14
  end if
! 
! Loading ...
! ------------
  nerr = fseek(llfu, ncurpos, 0) ! jump to the data offset
  if (nerr/=0) then
    nfail = 15
    goto 15
  end if
  !
  if (naty>0) then
    ! Loading loop over all atom types in the table
    do i=1, naty
      !
      natp = 0
      natr = 0
      !
      read(unit=llfu,iostat=nerr) EMS_SLI_data_attyz(i) ! ATOM TYPE (Z)
      if (nerr/=0) then
        nfail = 16
        goto 16
      end if
      read(unit=llfu,iostat=nerr) EMS_SLI_data_atdat(1:nadat,i) ! ATOM DATA (dZ,Biso,Occ)
      if (nerr/=0) then
        nfail = 17
        goto 16
      end if
      !
      read(unit=llfu,iostat=nerr) EMS_SLI_data_niatp(i) ! NUMBER OF POSITIONS
      if (nerr/=0) then
        nfail = 21
        goto 16
      end if
      !
      natp = EMS_SLI_data_niatp(i) ! length of position table
      ! integrity of the table is not checked
      !
      if (natp>0) then
        ! loading loop for atom type position table
        do k=1, natp
          !
          read(unit=llfu,iostat=nerr) EMS_SLI_data_iaatpo(1:3, k, i) ! XYZ POSITION
          if (nerr/=0) then
            nfail = 22
            goto 16
          end if
          !
        end do
        !
      end if
      !
      read(unit=llfu,iostat=nerr) EMS_SLI_data_niatr(i) ! NUMBER OF TRANSITIONS
      if (nerr/=0) then
        nfail = 18
        goto 16
      end if
      !
      natr = EMS_SLI_data_niatr(i) ! length of transition table
      ! integrity of the table is not checked
      !
      if (natr>0) then
        ! loading loop for atom type transition table
        do j=1, natr
          !
          read(unit=llfu,iostat=nerr) EMS_SLI_data_iatrcd(j, i) ! TRANS CODE
          if (nerr/=0) then
            nfail = 19
            goto 16
          end if
          read(unit=llfu,iostat=nerr) EMS_SLI_data_iatren(j, i) ! TRANS ENERGY (eV)
          if (nerr/=0) then
            nfail = 20
            goto 16
          end if
          !
        end do
        !
      end if
      !
    end do
    !
  end if
  !
  ! Should we apply byte swap or not?
  if (1==EMS_SLI_data_swap) then
    ! swap the data bytes
    do i=1, EMS_SLI_data_niaty
      itmp1 = EMS_SLI_data_attyz(i)
      itmp2 = SwapInt4(itmp1)
      EMS_SLI_data_attyz(i) = itmp2
      do j=1, nadat
        rtmp1 = EMS_SLI_data_atdat(j,i)
        rtmp2 = SwapReal4(rtmp1)
        EMS_SLI_data_atdat(j,i) = rtmp2
      end do
      itmp1 = EMS_SLI_data_niatr(i)
      itmp2 = SwapInt4(itmp1)
      EMS_SLI_data_niatr(i) = itmp2
      do j=1, EMS_SLI_data_niatr(i)
        itmp1 = EMS_SLI_data_iatrcd(j,i)
        itmp2 = SwapInt4(itmp1)
        EMS_SLI_data_iatrcd(j,i) = itmp2
        rtmp1 = EMS_SLI_data_iatren(j,i)
        rtmp2 = SwapReal4(rtmp1)
        EMS_SLI_data_iatren(j,i) = rtmp2
      end do
      itmp1 = EMS_SLI_data_niatp(i)
      itmp2 = SwapInt4(itmp1)
      EMS_SLI_data_niatp(i) = itmp2
      do k=1, EMS_SLI_data_niatp(i)
        do j=1, 3
          rtmp1 = EMS_SLI_data_iaatpo(j,k,i)
          rtmp2 = SwapReal4(rtmp1)
          EMS_SLI_data_iaatpo(j,k,i) = rtmp2
        end do
      end do
    end do
  end if
! ------------
! 
! ------------


! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_LOAD_ITAB: EXIT."
  return
  
! error handling
12 nerr = 1
  call EMS_ERROR("Failed to read parameters, file not connected.", subnum+nfail)
  return
13 nerr = 2
  call EMS_ERROR("Failed to read parameter from file.", subnum+nfail)
  return
14 nerr = 3
  call EMS_ERROR("Failed to allocate array memory.", subnum+nfail)
  return
15 nerr = 4
  call EMS_ERROR("File positioning failed.", subnum+nfail)
  return
16 nerr = 5
  call EMS_ERROR("Failed to read parameter table from file.", subnum+nfail)
  return

END SUBROUTINE EMS_SLI_LOAD_ITAB
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_SAVE_ITAB(llfu, nerr)
! function: Saves the inelastic data table from the local variables and
!           arrays to a file.
!           Requires slice version 2012071101 or higher and a
!           file opened for binary writing and connected to a logical
!           file unit llfu.
!           The routine assumes that the global variable
!           EMS_SLI_data_itab_offset is already set to the desired
!           position of the table.
! -------------------------------------------------------------------- !
! parameter: 
!  integer*4, intent(in) :: llfu    = logical file unit for reading
!  integer*4, intent(inout) :: nerr = error codes (0 = no error)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1600
  integer*4, intent(in) :: llfu
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nfail        ! local error number
  logical :: isopen
  integer*4 :: ncurpos      ! current read offset
  integer*4 :: i, j, k      ! iterators
  integer*4 :: naty, natr, natrmax, natp, natpmax, nadat ! number of data
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_SAVE_ITAB: INIT."
  nerr = 0
  nfail = 0
  nadat = EMS_SLI_data_natdat
! ------------


! ------------
! check if the file is open
  inquire(unit=llfu, OPENED=isopen)
  if (.not.isopen) then
    nfail = 1
    goto 12
  end if
! ------------


! ------------
! write the table at the preset position in the logical file unit
  ncurpos = EMS_SLI_data_itab_offset 
! ------------
! write the table size parameters
  naty = 0
  natrmax = 0
  natpmax = 0
  ! get number of data-sets
  if (allocated(EMS_SLI_data_attyz)) then
    naty = SIZE(EMS_SLI_data_attyz)
    if (naty>0) then
      if (allocated(EMS_SLI_data_niatr)) then
        natrmax = MAXVAL(EMS_SLI_data_niatr)
      end if
      if (allocated(EMS_SLI_data_niatp)) then
        natpmax = MAXVAL(EMS_SLI_data_niatp)
      end if
    end if
  end if
  !
  ! write all the numbers ...
  ! NUMBER OF ATOM TYPES
  EMS_SLI_data_niaty = naty
  call EMS_SLI_PRMWRITE_INT32(llfu, ncurpos, naty, nerr)
  if (nerr/=0) then
    nfail = 3
    goto 13
  end if
  ncurpos = ncurpos + 4
  ! MAX. NUMBER OF POSITIONS PER ATOM TYPE
  call EMS_SLI_PRMWRITE_INT32(llfu, ncurpos, EMS_SLI_data_niatp_max, nerr)
  if (nerr/=0) then
    nfail = 4
    goto 13
  end if
  ncurpos = ncurpos + 4
  ! MAX. NUMBER OF TRANSITIONS PER ATOM TYPE
  EMS_SLI_data_niatr_max = natrmax
  call EMS_SLI_PRMWRITE_INT32(llfu, ncurpos, EMS_SLI_data_niatr_max, nerr)
  if (nerr/=0) then
    nfail = 5
    goto 13
  end if
  ncurpos = ncurpos + 4
! ------------


! ------------
! Table writing
! ------------
  nerr = fseek(llfu, ncurpos, 0) ! jump to the data offset
  if (nerr/=0) then
    nfail = 15
    goto 15
  end if
  !
  if (naty>0) then
    ! Saving loop over all atom types in the table
    do i=1, naty
      !
      write(unit=llfu,iostat=nerr) EMS_SLI_data_attyz(i) ! ATOM TYPE (Z)
      if (nerr/=0) then
        nfail = 16
        goto 16
      end if
      write(unit=llfu,iostat=nerr) EMS_SLI_data_atdat(1:nadat,i) ! ATOM DATA (dZ, Biso, Occ)
      if (nerr/=0) then
        nfail = 17
        goto 16
      end if
      !
      write(unit=llfu,iostat=nerr) EMS_SLI_data_niatp(i) ! NUMBER OF POSITIONS
      if (nerr/=0) then
        nfail = 21
        goto 16
      end if
      natp = EMS_SLI_data_niatp(i)
      !
      !
      if (natp>0) then
        ! Saving loop for atom type position table
        do k=1, natp
          !
          write(unit=llfu,iostat=nerr) EMS_SLI_data_iaatpo(1:3, k, i) ! XYZ POSITION
          if (nerr/=0) then
            nfail = 22
            goto 16
          end if
          !
        end do
        !
      end if
      !
      write(unit=llfu,iostat=nerr) EMS_SLI_data_niatr(i) ! NUMBER OF TRANSITIONS
      if (nerr/=0) then
        nfail = 18
        goto 16
      end if
      natr = EMS_SLI_data_niatr(i)
      !
      if (natr>0) then
        ! Saving loop for atom type transition table
        do j=1, natr
          !
          write(unit=llfu,iostat=nerr) EMS_SLI_data_iatrcd(j, i) ! TRANS CODE
          if (nerr/=0) then
            nfail = 19
            goto 16
          end if
          write(unit=llfu,iostat=nerr) EMS_SLI_data_iatren(j, i) ! TRANS ENERGY (eV)
          if (nerr/=0) then
            nfail = 20
            goto 16
          end if
          !
        end do
        !
      end if
      !
    end do
    !
  end if
  !
! ------------


! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_SAVE_ITAB: EXIT."
  return
  
! error handling
12 nerr = 1
  call EMS_ERROR("Failed to write parameters, file not connected.", subnum+nfail)
  return
13 nerr = 2
  call EMS_ERROR("Failed to write parameter to file.", subnum+nfail)
  return
14 nerr = 3
  call EMS_ERROR("Failed to allocate array memory.", subnum+nfail)
  return
15 nerr = 4
  call EMS_ERROR("File positioning failed.", subnum+nfail)
  return
16 nerr = 5
  call EMS_ERROR("Failed to write parameter table to file.", subnum+nfail)
  return

END SUBROUTINE EMS_SLI_SAVE_ITAB
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_loaddata(sfile,nx,ny,nv,cdata,szz,nerr)
! function: loads ems sli data from specified file, only elastic slices
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 300
  character*(*), intent(in) :: sfile
  integer*4, intent(in) :: nx, ny, nv
  integer*4, intent(inout) :: nerr
  real*4, intent(inout) :: szz
  complex*8, intent(inout) :: cdata(1:nx,1:ny,1:nv)
  
  integer*4 :: i, j, k, nlfu, nskip, nfail
  integer*4 :: nnx, nny, nnv, mnv
  real*4 :: rtmp1, rtmp2, ht, szx, szy
  real*4, external :: SwapReal4
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_loaddata: INIT."
  EMS_SLI_data_swap = 0 ! preset byte swap to OFF by default
  EMS_SLI_data_alt_offset = EMS_SLI_POS_DATA ! preset the data offset to the default data offset
  nfail = 0
! ------------


! ------------
! re-read the parameters from this file
! - do this before opening the file for reading the data
! - for sake of security
! - for determining byte swap
! - and for determining the actual data offset (will be stored in EMS_SLI_data_alt_offset)
  call EMS_SLI_loadparams(trim(sfile),nnx,nny,nnv,ht,szx,szy,szz,nerr)
  if (nerr/=0) then
    nfail = 1
    goto 13
  end if
! ------------

! ------------
! check data size consistency
  if (nnx/=nx .or. nny/=ny) then
    nfail = 2
    goto 14
  end if
  ! check data size validity
  if (nnx<=0 .or. nny<=0) then
    nfail = 3
    goto 15
  end if
! store minimum number of slice variants
  mnv = min(nv,nnv)
! ------------

! ------------
  nskip = EMS_SLI_data_alt_offset ! in any case EMS_SLI_data_alt_offset holds the true data offset
! check if there is something to read at all 
  if (mnv<1) then
    ! number of variants for loading is zero or completely invalid
    ! - we simply skip reading the elastic data here, 
    goto 12
  end if
! ------------


! ------------
! START THE READING OF ELASTIC DATA
! ------------
! open binary file
  call EMS_GetFreeLFU(nlfu)
  open(unit=nlfu, file=trim(sfile), form='binary', access='sequential', &
     & iostat=nerr, status="old", action="read", share='DENYNONE' )
  if (nerr/=0) then
    call EMS_ERROR("Failed to connect to file", subnum+1)
    return
  end if
! ------------

! ------------
  nerr = fseek(nlfu, nskip, 0) ! jump to the data offset
  if (nerr/=0) then
    nfail = 4
    call EMS_ERROR("Sequential positioning failed.", subnum+nfail)
    return
  end if
  read(unit=nlfu,iostat=nerr) cdata(1:nx,1:ny,1:mnv) ! read the slices
  if (nerr/=0) then
    nfail = 5
    call EMS_ERROR("Failed to read data.", subnum+nfail)
    return
  end if
  if (1==EMS_SLI_data_swap) then
    ! swap the data bytes
    do k=1, nv
      do j=1, ny
        do i=1, nx
          rtmp1 = real(cdata(i,j,k))
          rtmp2 = imag(cdata(i,j,k))
          rtmp1 = SwapReal4(rtmp1)
          rtmp2 = SwapReal4(rtmp2)
          cdata(i,j,k) = cmplx(rtmp1,rtmp2)
        end do
      end do
    end do
  end if
! ------------

! ------------
  close (unit=nlfu, iostat=nerr)
  if (nerr/=0) then
    nfail = 6
    call EMS_ERROR("Failed to disconnect file.", subnum+nfail)
    return
  end if
! ------------
! FINISH THE READING OF ELASTIC DATA
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_loaddata: EXIT."
12 return
  
13 nerr = 1
  call EMS_ERROR("Failed read parameter from file.", subnum+nfail)
  return
14 nerr = 2
  call EMS_ERROR("Inconsitency between function parameters and file parameters.", subnum+nfail)
  return 
15 nerr = 3
  call EMS_ERROR("Invalid data size found in file parameters.", subnum+nfail)
  return 

END SUBROUTINE EMS_SLI_loaddata
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_loadparams(sfile,nx,ny,nv,ht,szx,szy,szz,nerr)
! function: loads ems sli parameters from specified file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 400
  character*(*), intent(in) :: sfile
  integer*4, intent(inout) :: nx, ny, nv, nerr
  real*4, intent(inout) :: szx, szy, szz,ht
  
  integer*4 :: nlfu, nfail
  integer*4, external :: SwapInt4
  real*4, external :: SwapReal4
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_loadparams: INIT."
  EMS_SLI_data_swap = 0
  nv = 1 ! default
  EMS_SLI_data_nvar = 1 ! default
  EMS_SLI_data_alt_offset = EMS_SLI_POS_DATA ! default
  EMS_SLI_data_itab_offset = EMS_SLI_ITAB_OFFSET_DEFAULT ! default
  nfail = 0
! ------------

! ------------
! open binary file
  call EMS_GetFreeLFU(nlfu)
  open(unit=nlfu, file=trim(sfile), form='binary', access='sequential', &
     & iostat=nerr, status="old", action="read", share='DENYNONE' )
  if (nerr/=0) then
    nfail = 1
    call EMS_ERROR("Failed to connect to file ["//trim(sfile)//"]", subnum+nfail)
    return
  end if
! ------------


! ------------
! check byte swap and load slice sampling along x
  call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_DATASIZE_X, EMS_SLI_data_dimx, nerr)
  if (nerr/=0) then
    nfail = 2
    goto 13
  end if
  if (EMS_SLI_data_dimx>32000.or.EMS_SLI_data_dimx<0) then
    EMS_SLI_data_swap = 1
    EMS_SLI_data_dimx = SwapInt4(EMS_SLI_data_dimx)
  end if
  nx = EMS_SLI_data_dimx
! ------------

! ------------
! load slice sampling along y
  call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_DATASIZE_Y, EMS_SLI_data_dimy, nerr)
  if (nerr/=0) then
    nfail = 3
    goto 13
  end if
  ny = EMS_SLI_data_dimy
! ------------

! ------------
! load image format id
  call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_IFORM, EMS_SLI_data_iform, nerr)
  if (nerr/=0) then
    nfail = 4
    goto 13
  end if
  if (EMS_SLI_data_iform/=EMS_SLI_IFORM) then
    ! we expect this one format, even if inelastic potentials are saved in the file
    ! the flag for existence of inelastic data is the number of atom types
    ! saved at the beginning of the inelastic data table
    call EMS_ERROR("Unexpected data format.", subnum+6)
    return
  end if
! ------------

! ------------
! load title string
  call EMS_SLI_PRMREAD_STRING(nlfu, EMS_SLI_POS_TITLE, EMS_SLI_data_title, nerr)
  if (nerr/=0) then
    nfail = 5
    goto 13
  end if
! ------------

! ------------
! load extended header version
  call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_EXTHEADER_VER, EMS_SLI_data_ver, nerr)
  if (nerr/=0) then
    nfail = 6
    goto 13
  end if
  ! ------------
  ! load the extended header components depending on the header version
  select case (EMS_SLI_data_ver)
  case (2010112301) ! first version
    ! ------------
    ! number of variants
    call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_VARIANT_NUM, EMS_SLI_data_nvar, nerr)
    if (nerr/=0) then
      nfail = 7
      goto 13
    end if
    if (EMS_SLI_data_nvar<0) EMS_SLI_data_nvar = 1 ! fail saveness
    nv = EMS_SLI_data_nvar ! transfer to dummy
    ! ------------
  case (2012071101) ! inelastic slices
    ! ------------
    ! number of variants
    call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_VARIANT_NUM, EMS_SLI_data_nvar, nerr)
    if (nerr/=0) then
      nfail = 12
      goto 13
    end if
    if (EMS_SLI_data_nvar<0) EMS_SLI_data_nvar = 1 ! fail saveness
    nv = EMS_SLI_data_nvar ! transfer to dummy
    ! ------------
    ! content type
    call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_CONTENT_TYPE, EMS_SLI_data_ctype, nerr)
    if (nerr/=0) then
      nfail = 13
      goto 13
    end if
    if (EMS_SLI_data_ctype<0) EMS_SLI_data_ctype = 0 ! fail saveness, expect phase grating
    ! ------------
    ! number of sub-slices
    call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_SUBSLC_NUM, EMS_SLI_data_sslnum, nerr)
    if (nerr/=0) then
      nfail = 14
      goto 13
    end if
    if (EMS_SLI_data_sslnum<0) EMS_SLI_data_sslnum = 1 ! fail saveness, there is at least one slice
    ! ------------
    ! sub-slice thickness
    call EMS_SLI_PRMREAD_REAL32(nlfu, EMS_SLI_POS_SUBSLC_THICK, EMS_SLI_data_sslz, nerr)
    if (nerr/=0) then
      nfail = 15
      goto 13
    end if
    if (EMS_SLI_data_sslz<0.0) EMS_SLI_data_sslz = 0.0 ! fail saveness, 0.0 means: no propagator
    ! ------------
    ! alternative elastic data offset (EMS_SLI_POS_ALT_OFFSET)
    call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_ALT_OFFSET, EMS_SLI_data_alt_offset, nerr)
    if (nerr/=0) then
      nfail = 16
      goto 13
    end if
    ! cross-check if the alternative data offset is reasonable
    ! if not: set to default
    if (EMS_SLI_data_alt_offset<EMS_SLI_POS_DATA) then
      EMS_SLI_data_alt_offset = EMS_SLI_POS_DATA
    end if
    ! ------------
    ! inelastic header table table offset
    call EMS_SLI_PRMREAD_INT32(nlfu, EMS_SLI_POS_ITAB_OFFSET, EMS_SLI_data_itab_offset, nerr)
    if (nerr/=0) then
      nfail = 18
      goto 13
    end if
    ! cross-check if the header tables are present
    ! the tables MUST be present from version #2012071101 on, if not set to default
    if (EMS_SLI_data_itab_offset<EMS_SLI_ITAB_OFFSET_DEFAULT) then
      EMS_SLI_data_itab_offset = EMS_SLI_ITAB_OFFSET_DEFAULT
    end if
    ! ------------
    ! inelastic header tables
    call EMS_SLI_LOAD_ITAB(nlfu, nerr)
    if (nerr/=0) then
      nfail = 20
      goto 13
    end if
  end select ! case (EMS_SLI_data_ver)
  ! ------------
! ------------

! ------------
! total slice thickness
  call EMS_SLI_PRMREAD_REAL32(nlfu, EMS_SLI_POS_THICKNESS, EMS_SLI_data_szz, nerr)
  if (nerr/=0) then
    nfail = 50
    goto 13
  end if
  szz = EMS_SLI_data_szz
! ------------

! ------------
! slice data high tension [kV]
  call EMS_SLI_PRMREAD_REAL32(nlfu, EMS_SLI_POS_HIGHTENSION, EMS_SLI_data_ht, nerr)
  if (nerr/=0) then
    nfail = 51
    goto 13
  end if
  ht = EMS_SLI_data_ht
! ------------

! ------------
! slice size along x in [nm]
  call EMS_SLI_PRMREAD_REAL32(nlfu, EMS_SLI_POS_PHYSSIZE_X, EMS_SLI_data_szx, nerr)
  if (nerr/=0) then
    nfail = 52
    goto 13
  end if
  szx = EMS_SLI_data_szx
! ------------

! ------------
! slice size along y in [nm]  
  call EMS_SLI_PRMREAD_REAL32(nlfu, EMS_SLI_POS_PHYSSIZE_Y, EMS_SLI_data_szy, nerr)
  if (nerr/=0) then
    nfail = 53
    goto 13
  end if
  szy = EMS_SLI_data_szy
! ------------

! ------------
  close (unit=nlfu, iostat=nerr)
  if (nerr/=0) then
    nfail = 54
    call EMS_ERROR("Failed to disconnect file.", subnum+nfail)
    return
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_loadparams: EXIT."
  return
  
13 nerr = 1
  call EMS_ERROR("Failed read parameter from file.", subnum+nfail)
  return
14 nerr = 2
  call EMS_ERROR("Invalid header table offset read from primary header.", subnum+nfail)
  return
15 nerr = 5
  call EMS_ERROR("Invalid alternative data offset read from primary header.", subnum+nfail)
  return

END SUBROUTINE EMS_SLI_loadparams
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_save(sfile,nx,ny,nv,szx,szy,szz,ht,cdata,nerr)
! function: saves ems sli parameters and data to specified file
! -------------------------------------------------------------------- !
! parameter: sfile :: character(len=*)  = slice file name (absolute on HD)
!            nx, ny :: integer*4        = slice discretization (x and y)
!            nv :: integer*4            = number of slice variants
!            szx, szy :: real*4         = slice size (x and y) in nm
!            szz :: real*4              = slice thickness [nm]
!            ht :: real*4               = TEM high tension
!            cdata(nx,ny,nv) :: complex*8  = slice data
!            nerr :: integer*4          = error code
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 400
  character*(*), intent(in) :: sfile
  integer*4, intent(in) :: nx, ny, nv
  integer*4, intent(inout) :: nerr
  real*4, intent(in) :: szx, szy, szz, ht
  complex*8, intent(inout) :: cdata(nx,ny,nv)
  
  integer*4 :: nlfu, nskip, nfail
  integer*4 :: data_offset, itab_offset, itab_size, nadat
  integer*4 :: header(2048)
  external :: createfilefolder
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_save: INIT."
  nerr = 0
  nfail = 0
  nadat = EMS_SLI_data_natdat
  if (nv<1) then
    call EMS_ERROR("Invalid parameter (number of variants).", subnum+1)
    return
  end if
  header = 0
! ------------


! ------------
! open binary file
  call EMS_GetFreeLFU(nlfu)
  call createfilefolder(trim(sfile),nerr)
  open(unit=nlfu, file=trim(sfile), form='binary', access='sequential', &
     & iostat=nerr, status="replace", action="write", share='DENYRW' )
  if (nerr/=0) then
    call EMS_ERROR("Failed to connect to file", subnum+2)
    return
  end if
! ------------

! ------------ PAY ATTENTION !
!
! With sequential file access you MUST WRITE THE DATA IN IT'S
! INTENDED SEQUENCE to the file unit. Using FSEEK to jump back
! will erase data writen before and beyond the jump-point.
!
! ------------

! ------------ SLICE SAMPLING X
  nskip = EMS_SLI_POS_DATASIZE_X !                      = Z'00000000'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 1
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) nx
! ------------ SLICE SAMPLING Y
  nskip = EMS_SLI_POS_DATASIZE_Y !                      = Z'00000004'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 1
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) ny
! ------------ DATA FORM
  nskip = EMS_SLI_POS_IFORM !                           = Z'00000010'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 3
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) EMS_SLI_IFORM
! ------------ SLICE TITLE STRING
  nskip = EMS_SLI_POS_TITLE !                           = Z'0000004C'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 4
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) EMS_SLI_data_title
! ------------ EXTENDED HEADER VERSION
  nskip = EMS_SLI_POS_EXTHEADER_VER !                   = Z'00000074'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 5
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) EMS_SLI_VERSION
! ------------ NUMBER OF FROZEN LATTICE VARIANTS
  nskip = EMS_SLI_POS_VARIANT_NUM !                     = Z'00000078'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 6
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) nv
! ------------ DATA CONTENT TYPE
  nskip = EMS_SLI_POS_CONTENT_TYPE !                    = Z'00000080'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 6
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) EMS_SLI_data_ctype
! ------------ DECIDE ON WHERE TO WRITE TABLES AND LARGE DATA!
  ! - atom table offset = default = Z'00001100'
  itab_offset = EMS_SLI_ITAB_OFFSET_DEFAULT
  EMS_SLI_data_itab_offset = itab_offset
  ! - get size of the inelastic data table
  itab_size = 4*(3 + EMS_SLI_data_niaty*( 3 + nadat &
     &        + 2*EMS_SLI_data_niatr_max + 3*EMS_SLI_data_niatp_max ) )
  ! - determine offset for the elastic data
  if (itab_offset + itab_size + 16 > EMS_SLI_POS_DATA) then
    ! set a new / alternative data offset
    data_offset = itab_offset + itab_size + 16
  else
    ! set to default data offset = Z'00002000'
    data_offset = EMS_SLI_POS_DATA 
  end if
! ------------ ELASTIC DATA OFFSET
  nskip = EMS_SLI_POS_ALT_OFFSET !                      = Z'00000088'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 13
  write(unit=nlfu,iostat=nerr,err=13) data_offset
  if (nerr/=0) goto 14
! ------------ INELASTIC AND ATOM POSITION TABLE OFFSET
  nskip = EMS_SLI_POS_ITAB_OFFSET !                     = Z'0000008C'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 13
  
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) itab_offset
! ------------ SLICE THICKNESS
  nskip = EMS_SLI_POS_THICKNESS !                       = Z'00001034'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 7
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) szz
! ------------ DATA HIGH TENSION / ELECTRON ENERGY
  nskip = EMS_SLI_POS_HIGHTENSION !                     = Z'00001040'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 8
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) ht
! ------------ SLICE PHYSICAL X SIZE
  nskip = EMS_SLI_POS_PHYSSIZE_X !                      = Z'00001044'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 9
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) szx
! ------------ SLICE PHYSICAL Y SIZE
  nskip = EMS_SLI_POS_PHYSSIZE_Y !                      = Z'00001048'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 10
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) szy
! ------------ INELASTIC AND ATOM POSITION TABLE
  nskip = itab_offset !                                 = Z'00001100'                                   
  nerr = fseek(nlfu, nskip, 0)
  nfail = 14
  call EMS_SLI_SAVE_ITAB(nlfu, nerr)
  if (nerr/=0) goto 14
! ------------ ELASTIC PHASE GRATING DATA <-------------------------------------------------------
  nskip = data_offset !                                >= Z'00002000'
  nerr = fseek(nlfu, nskip, 0)
  nfail = 15
  if (nerr/=0) goto 14
  write(unit=nlfu,iostat=nerr,err=13) cdata ! EMS_SLI_POS_DATA
! ------------ DONE. CLOSE THE FILE !
  close (unit=nlfu, iostat=nerr)
  if (nerr/=0) then
    call EMS_ERROR("Failed to disconnect file.", subnum+11)
    return
  end if
! -----------------------------------------------------------------------------------------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_save: EXIT."
  return
  
13 nerr = 1
  call EMS_ERROR("Failed to write data.", subnum+51)
  return
14 nerr = 2
  call EMS_ERROR("Sequential positioning failed.", subnum+nfail)
  return

END SUBROUTINE EMS_SLI_save
!**********************************************************************!






!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_settab0()
! function: Deletes all atom tables to be saved to the slice file.
! -------------------------------------------------------------------- !
! parameter: none.
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1800
  
  integer*4 :: nierr
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_settab0: INIT."
  nierr = 0
  call EMS_SLI_ALLOC_ITAB(0, 0, 0, nierr)
  return
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_settab0: EXIT."
  return

END SUBROUTINE EMS_SLI_settab0
!**********************************************************************!





     
!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_SLI_settab1(fz0,fz1,nau,auhash,na,atn,atc,atb,ato,atp)
! function: Creates an atom table to be saved to the slice file
!           from input lists. The created table considers the elastic
!           scattering information only. No transition table will be
!           created.
! -------------------------------------------------------------------- !
! parameter:
!           fz0 = slice fractional z-offset in the super-cell
!           fz1 = slice fractional z-termination in the super-cell
!           nau = number of atoms in this slice
!           auhash(nau) = list of atom indices used in this slice
!           na = number of all atoms in the supercell
!           atn(na) = atomic number list
!           atc(na) = ionic charge list
!           atb(na) = Biso list
!           ato(na) = Occupancy list
!           atp(3,na) = (fx,fy,fz) list of positions in the supercell
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1700
  
  integer*4, intent(in) :: nau, auhash(nau), na, atn(na)
  real*4, intent(in) :: fz0, fz1, atc(na), atb(na), ato(na), atp(3,na)
  
  integer*4 :: nierr, nalloc
  integer*4 :: i, j, k, i1, i2
  integer*4 :: naty, natpmax, natrmax, nadat
  integer*4 :: natp(nau), atty(nau)
  real*4 :: dfz, fps
  real*4, allocatable :: atdat(:,:), atpos(:,:,:)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > EMS_SLI_settab1: INIT."
  nierr = 0
  naty = 0
  natpmax = 0
  natrmax = 0
  natp = 0
  atty = 0
  atdat = 0.0
  atpos = 0.0
  dfz = abs(fz1 - fz0)
  if (nau<=0 .or. na<=0 .or. dfz<=0.0) then
    call EMS_SLI_ALLOC_ITAB(naty, natrmax, natpmax, nierr)
    return
  end if
  nalloc = 0
  nadat = EMS_SLI_data_natdat
  allocate( atdat(nadat,nau), atpos(3,nau,nau), stat=nalloc )
  if (nalloc/=0) return
! ------------

! ------------
! Determine how many different atom types exist in this
! slice and fill the local arrays accordingly
  do i=1, nau ! loop over used atoms in the slice
    ! Get the atom index in the full atom table
    j = auhash(i)
    ! Consistency check
    if (j<1 .or. j>na) cycle ! invalid atom index, skip this atom
    !
    ! Check if the atom type is already defined
    if (naty==0) then ! no atom type defined yet
      ! Create a new atom type + data
      ! add also one atom position to the position table
      goto 100
      !
    else ! there are atom types defined already
      !
      ! Search existing atom data
      i1 = 0
      do i2=1, naty
        ! compare
        if (        atty(i2) == atn(j) & ! atomic number
            & .and. abs( atdat(1,i2)-atc(j) ) < 0.001 & ! ionic charge
            & .and. abs( atdat(2,i2)-atb(j) ) < 0.00001 & ! Biso
            & .and. abs( atdat(3,i2)-ato(j) ) < 0.001 & ! Occupancy
            & ) then ! atom type exists
          i1 = i2
          exit
        end if
      end do
      !
      if (i1>0) then ! the atom type is known ...
        goto 101 ! ... just add the atom position
      else ! atom type doesn't exist yet
        ! Create new atom type and add first position
        goto 100
      end if
      !
    end if
    
    ! the program should never reach this point
    cycle ! if so, skip this atom
    
    ! JUMP POINT: Create a new atom type + data
100 naty = naty + 1
    i1 = naty ! update atom type index for adding new position
    atty(naty) = atn(j)
    atdat(1,naty) = atc(j)
    atdat(2,naty) = atb(j)
    atdat(3,naty) = ato(j)
    ! when a new atom type is created, this type gets also its first position data
    ! ...
    ! JUMP POINT: Add a new position at atom type index i1 
101 natp(i1) = natp(i1) + 1 
    k = natp(i1)
    atpos(1,k,i1) = atp(1,j)
    atpos(2,k,i1) = atp(2,j)
    ! rescale the fractional z coordinate of the atom in the supercell
    ! to a fractional z coordinate of the slice
    fps = (atp(3,j)-fz0)/dfz
    ! ... limit to range [0,1)
    fps = min(0.999999,max(0.0,fps))
    atpos(3,k,i1) = fps
    
  end do
! ------------

! ------------
! Data is now stored in the local arrays.
! Transfer the data to the gobal module arrays
  ! get the max. number of atom positions
  natpmax = maxval(natp)
  ! allocate the global arrays, this also sets the global size numbers
  nierr = 0
  call EMS_SLI_ALLOC_ITAB(naty, natrmax, natpmax, nierr)
  if (nierr/=0) goto 299
  ! copy the arrays
  EMS_SLI_data_attyz(1:naty) = atty(1:naty)
  EMS_SLI_data_atdat(1:nadat,1:naty) = atdat(1:nadat,1:naty)
  EMS_SLI_data_niatp(1:naty) = natp(1:naty)
  EMS_SLI_data_iaatpo(1:3,1:natpmax,1:naty) = atpos(1:3,1:natpmax,1:naty)
! ------------

! ------------
299 continue
  if (allocated(atdat)) deallocate( atdat, stat=nalloc )
  if (allocated(atpos)) deallocate( atpos, stat=nalloc )
! ------------

! ------------
!  write(unit=*,fmt=*) " > EMS_SLI_settab1: EXIT."
  return

END SUBROUTINE EMS_SLI_settab1
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
SUBROUTINE EMS_ERROR(sTxt,nErr)
! function: print error message
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
!            INTEGER*4 :: nErr : error code, displayed when /= 0
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  character*(*), intent(in) :: sTxt
  integer*4, intent(in) :: nErr
  character(len=EMS_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = EMS_el
! ----------

! ----------
! messaging (to screen)
  write(unit=6,fmt=*) " "
  write(unit=6,fmt=*) "ERROR : EMS-Module : ",trim(sTxt)," Error code:",nErr
  write(unit=6,fmt=*) " Halting program."
  write(unit=6,fmt=*) " "
  stop
! ----------
  
  return

END SUBROUTINE EMS_ERROR
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE EMS_WARN(sTxt)
! function: print warning message
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
!            INTEGER*4 :: nErr : error code, displayed when /= 0
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  character*(*), intent(in) :: sTxt
  character(len=EMS_ll) :: sinfo
! ----------

! ----------
! init
  sinfo = EMS_el
! ----------

! ----------
! messaging (to screen)
  write(unit=6,fmt=*) "Warning : EMS-Module : ",trim(sTxt)
! ----------
  
  return

END SUBROUTINE EMS_WARN
!**********************************************************************!




END MODULE EMSdata
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
!  integer*4, parameter :: subnum = 1900
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

