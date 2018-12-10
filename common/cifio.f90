!**********************************************************************!
!**********************************************************************!
!
! F90 Module - cifio (CIF)
!
! Author: Dr. J. Barthel
!         Forschungszentrum Juelich GmbH
!         Jülich, Germany
!         ju.barthel@fz-juelich.de
!         first version: 05.01.2016
!         last version: 21.03.2018
!
! Purpose: Input and output of structure data from CIF files.
!          The most recently loaded data is stored in module variables.
!          CIF file structure is interpreted relatively strict
!          according to the CIF dictionary version 2.4:
!          ftp://ftp.iucr.org/pub/cif_core.dic
!
! Todo:    Application of symmetry operations on the anisotropic
!          thermal vibration tensors.
!
! Linkage: binio2.f90 (general I/O routines)
!          spacegroups.f90 (module with spacegroup data)
!          symops.f90 (module handling symmetry operations)
!
! Remarks:
!   - The unit cell (super cell) is defines by a, b, c, alpha, beta,
!     and gamma, all related lengths and reciprocal lengths are
!     calculated from those numbers. These numbers MUST BE defined
!     by the CIF input. A quick check that this has happened is done
!     by calling CIF_get_cell_volume after calling CIF_READ
!
!
!**********************************************************************!
!**********************************************************************!
!* This program is free software: you can redistribute it and/or modify
!* it under the terms of the GNU General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or   
!* (at your option) any later version.                                 
!*                                                                     
!* This program is distributed in the hope that it will be useful,     
!* but WITHOUT ANY WARRANTY; without even the implied warranty of      
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
!* GNU General Public License for more details.                        
!*                                                                     
!* You should have received a copy of the GNU General Public License   
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
MODULE cifio

  use symops ! uses spacegroups
  
  implicit none
  
  save
  
  !
  ! variable and parameter definitions
  !
  integer, parameter, private :: dp = SELECTED_REAL_KIND(15,307)
  integer, parameter, private :: ssl = 128
  integer, parameter, private :: lsl = 1024
  integer, parameter, private :: CIF_lstat_len = 10
  integer, public, parameter :: CIF_atom_site_namealts = 2
  integer, parameter, public :: CIF_atom_site_items_num = 14
  real(dp), parameter, public :: CIF_8pi2 = 0.789568352d+2 ! 8 * Pi**2 (from Uiso to Biso)
  integer*4, public :: CIF_stdio
  DATA CIF_stdio /6/
  integer*4, public :: CIF_nerr, CIF_warn
  DATA CIF_nerr /0/ ! number of errors
  DATA CIF_warn /0/ ! number of warnings
  integer*4, public :: CIF_verbosity ! verbosity level
  DATA CIF_verbosity /0/ ! default verbosity level is root level (almost nothing)
  !
  ! <<< The content of the CIF file in line sequence >>>
  !     call CIF_CLEAR_LINES to get rid of this.
  !
  character(len=lsl), private, allocatable :: CIF_lines(:) ! the lines of a CIF file
  !
  ! <<< File structure information >>>
  !     To be allocated when needed.
  integer*4, private, allocatable :: CIF_lstats(:,:) ! information on loop_ constructs (tables)
  !
  ! <<< Structure data >>>
  !     When adding more structure data, make sure to update the
  !     CIF_read* and CIF_clear* routines!
  !
  real(dp), public :: CIF_cell_length_a
  real(dp), public :: CIF_cell_length_b
  real(dp), public :: CIF_cell_length_c
  real(dp), public :: CIF_cell_angle_alpha
  real(dp), public :: CIF_cell_angle_beta
  real(dp), public :: CIF_cell_angle_gamma
  character(len=ssl), public :: CIF_symmetry_space_group_name_HM
  integer, public :: CIF_symmetry_Int_Tables_number
  integer*4, public :: CIF_atom_site_number ! number of set atom site data
  DATA CIF_atom_site_number /0/
  real(dp), public, allocatable :: CIF_atom_site(:,:) ! atomic site data
  !
  ! <<< Auto item ID list >>>
  ! 
  ! structure of the internal _atom_site table
  character(len=ssl), public :: CIF_atom_site_items(CIF_atom_site_namealts, CIF_atom_site_items_num)
  DATA CIF_atom_site_items(1:2, 1) /"_atom_site_fract_x", "_atom_site_fract_x"/
  DATA CIF_atom_site_items(1:2, 2) /"_atom_site_fract_y", "_atom_site_fract_y"/
  DATA CIF_atom_site_items(1:2, 3) /"_atom_site_fract_z", "_atom_site_fract_z"/
  DATA CIF_atom_site_items(1:2, 4) /"_atom_site_type_symbol", "_atom_site_aniso_type_symbol"/
  DATA CIF_atom_site_items(1:2, 5) /"_atom_site_type_charge", "_atom_site_type_charge"/
  DATA CIF_atom_site_items(1:2, 6) /"_atom_site_occupancy", "_atom_site_occupancy"/
  DATA CIF_atom_site_items(1:2, 7) /"_atom_site_adp_type", "_atom_site_thermal_displace_type"/
  DATA CIF_atom_site_items(1:2, 8) /"_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv"/
  DATA CIF_atom_site_items(1:2, 9) /"_atom_site_aniso_U_11", "_atom_site_aniso_B_11"/
  DATA CIF_atom_site_items(1:2,10) /"_atom_site_aniso_U_22", "_atom_site_aniso_B_22"/
  DATA CIF_atom_site_items(1:2,11) /"_atom_site_aniso_U_33", "_atom_site_aniso_B_33"/
  DATA CIF_atom_site_items(1:2,12) /"_atom_site_aniso_U_12", "_atom_site_aniso_B_12"/
  DATA CIF_atom_site_items(1:2,13) /"_atom_site_aniso_U_13", "_atom_site_aniso_B_13"/
  DATA CIF_atom_site_items(1:2,14) /"_atom_site_aniso_U_23", "_atom_site_aniso_B_23"/
  ! data types for reading the _atom_site table from file: 0 = string, 1 = integer, 2 = float/double
  integer, public :: CIF_atom_site_dtype(CIF_atom_site_items_num)
  DATA CIF_atom_site_dtype( 1) /2/
  DATA CIF_atom_site_dtype( 2) /2/
  DATA CIF_atom_site_dtype( 3) /2/
  DATA CIF_atom_site_dtype( 4) /0/
  DATA CIF_atom_site_dtype( 5) /2/
  DATA CIF_atom_site_dtype( 6) /2/
  DATA CIF_atom_site_dtype( 7) /0/
  DATA CIF_atom_site_dtype( 8) /2/
  DATA CIF_atom_site_dtype( 9) /2/
  DATA CIF_atom_site_dtype(10) /2/
  DATA CIF_atom_site_dtype(11) /2/
  DATA CIF_atom_site_dtype(12) /2/
  DATA CIF_atom_site_dtype(13) /2/
  DATA CIF_atom_site_dtype(14) /2/
  !
  !
  !  routine definitions
  PRIVATE :: CIF_MSG, CIF_WRN, CIF_ERR
  PRIVATE :: CIF_find_atom_site_name
  PRIVATE :: CIF_get_value
  PRIVATE :: CIF_str_remsub
  PRIVATE :: CIF_read_item
  PRIVATE :: CIF_read_loop1, CIF_read_loop2
  PRIVATE :: CIF_get_iso_from_ani
  PRIVATE :: CIF_create_atom_type_list
  !  interface
  PUBLIC  :: CIF_READ    ! read and extract the content of a CIF file
  PUBLIC  :: CIF_WRITE   ! write the current structure data in a CIF file
  PUBLIC  :: CIF_CLEAR   ! clear the current data, set to default
  PUBLIC  :: CIF_get_atnum ! get the atomic number from an atom-site type symbol string
  PUBLIC  :: CIF_get_atom_type_symbol ! re-creates the atom-site type symbol from atomic number and oxidation state (charge)
  PUBLIC  :: CIF_get_atom_site_type_symbol ! re-creates the atom-site type symbol from the current data
  PUBLIC  :: CIF_get_atom_site_biso ! returns the biso value of an atom-site from the current data
  PUBLIC  :: CIF_get_atom_site_uiso ! returns the uiso value of an atom-site from the current data
  PUBLIC  :: CIF_get_cell_volume ! returns the cell volume in A^3 from a, b, c, alpha, beta, and gamma
  PUBLIC  :: CIF_get_cell_reciprocal ! returns the reciprocal cell parameters from a, b, c, alpha, beta, and gamma
  
  CONTAINS


!********************************************************
! CIF_MSG
! This subroutine posts a message depending on the
! module verbosity level.
! Currently this is a wrapper function calling a program
! specific messaging routine.
!********************************************************
SUBROUTINE CIF_MSG(smsg,nlev)
  implicit none
  character(len=*), intent(in) :: smsg
  integer, optional :: nlev
  integer :: ulev
  ulev = 0 ! default is to post base level messages (almost nothing)
  if (present(nlev)) ulev=nlev
  if (ulev<=CIF_verbosity) then
    call PostMessage(trim(smsg))
  end if
END SUBROUTINE CIF_MSG


!********************************************************
! CIF_WRN
! This subroutine posts a warning.
! Currently this is a wrapper function calling a program
! specific messaging routine.
!********************************************************
SUBROUTINE CIF_WRN(smsg)
  implicit none
  character(len=*), intent(in) :: smsg
  call PostWarning(trim(smsg))
  CIF_warn = CIF_warn+1
END SUBROUTINE CIF_WRN


!********************************************************
! CIF_ERR
! This subroutine posts an error message.
!********************************************************
SUBROUTINE CIF_ERR(smsg, error_code)
  implicit none
  character(len=*), intent(in) :: smsg
  integer, optional :: error_code
  character(len=20) :: scode
  scode=""
  if( present(error_code) ) THEN
    write (unit=scode,fmt=*) error_code
    write (unit=6,fmt='(A)') "ERROR: "//trim(smsg)// &
        & " ("//trim(adjustl(scode))//")"
  else
    write (unit=6,fmt='(A)') "ERROR: "//trim(smsg)
  end if
  CIF_nerr = CIF_nerr+1
END SUBROUTINE CIF_ERR


!********************************************************
! CIF_CLEAR
! This routine clears the current data, resets to default
! values and deallocates all arrays.
!********************************************************
SUBROUTINE CIF_CLEAR()
  implicit none
  integer :: nalloc
  nalloc = 0
  call CIF_CLEAR_LINES() ! the lines of the previous CIF file if still present
  if (allocated(CIF_atom_site)) deallocate(CIF_atom_site,stat=nalloc) ! atomic site data
  CIF_cell_length_a = 0.d+0
  CIF_cell_length_b = 0.d+0
  CIF_cell_length_c = 0.d+0
  CIF_cell_angle_alpha = 9.d+1
  CIF_cell_angle_beta = 9.d+1
  CIF_cell_angle_gamma = 9.d+1
  CIF_symmetry_space_group_name_HM = ""
  CIF_symmetry_Int_Tables_number = 0
END SUBROUTINE CIF_CLEAR


!********************************************************
! CIF_CLEAR_LINES
! This routine clears the current CIF file lines from
! memory by deallocation.
!********************************************************
SUBROUTINE CIF_CLEAR_LINES()
  implicit none
  integer :: nalloc
  nalloc = 0
  if (allocated(CIF_lines)) deallocate(CIF_lines,stat=nalloc)
END SUBROUTINE CIF_CLEAR_LINES


!********************************************************
! CIF_str_remsub
! This routine removes a substring from a string. The
! removal includes possible adjacent quotation marks.
!********************************************************
SUBROUTINE CIF_str_remsub(str,sub)
  implicit none
  character(len=*), intent(inout) :: str
  character(len=*), intent(in) :: sub
  integer :: i1, i2, l1, l2
  integer :: qmark
  qmark = 0
  l1 = len(str)
  l2 = len(sub)
  if (l1==0) return ! emtpy string ! really nothing to delete
  if (l2==0) then   ! empty substring ! there could be quotations marks in str
    call CIF_str_remsub(str,"''")
    call CIF_str_remsub(str,'""')
    return
  end if
  i1 = index(str,sub) ! start index of first occurance
  if (i1<1) return ! substring not found
  i2 = i1 + l2-1 ! stop index of sub in str
  ! replace content by space characters
  ! include the two adjacent characters, they should be
  ! spaces anyway in CIF, or could be quotation marks
  i1 = max(1, i1-1)
  ! check if there is a quotation mark right before the sub-string
  if ( str(i1:i1)=="'" ) qmark = 1 ! type 1 quo-mark
  if ( str(i1:i1)=='"' ) qmark = 2 ! type 2 quo-mark
  ! in case of quatation marks, scan for the next occurance,
  if (qmark>0) then
    if (qmark==1) i2 = i2+index(str(i2:l1),"'")-1
    if (qmark==2) i2 = i2+index(str(i2:l1),'"')-1
  else ! otherwise just include the next character for deletion
    i2 = min(i2+1,l1)
  end if
  str(i1:i2) = repeat(" ",i2-i1+1)
END SUBROUTINE CIF_str_remsub


!********************************************************
! CIF_get_value
! This routine interpretes the string sval as numerical
! data and tries to transform it to dval. Precision
! brackets as often found in CIF files will be ignored.
!********************************************************
SUBROUTINE CIF_get_value(sval, dval, nchk)
  implicit none
  character(len=ssl), intent(in) :: sval
  real(dp), intent(out) :: dval
  integer, intent(out) :: nchk
  integer :: ilen, ipar, ioerr
  character(len=ssl) :: stmp
  nchk = 0
  ipar = 0
  ioerr = 0
  dval = 0.0d+0
  if (0==len_trim(sval)) return ! zero length strings are interpreted as 0
  if ("."==trim(sval)) return ! "." is interpreted as 0
  ilen = len_trim(sval) ! get length of the data string
  ipar = index(sval,"(") ! get position of a possible precision bracket
  if (ipar>0) ilen = ipar-1 ! shorten the string length for interpretation
  stmp = sval(1:ilen) ! prepare the final string for interpretation
  read(unit=stmp,fmt=*,iostat=ioerr) dval ! read from string to numerical variable
  if (ioerr/=0) nchk = ioerr
END SUBROUTINE CIF_get_value


!********************************************************
! CIF_get_atnum
! This routine translates from an atomic site symbol
! string to an atomic number and oxidation state if
! possible
!********************************************************
SUBROUTINE CIF_get_atnum(assymb,atnum,atoxs)
  implicit none
  character(len=*), intent(in) :: assymb
  integer, intent(out) :: atnum
  real(dp), intent(out) :: atoxs
  integer, parameter :: maxato = 116
  integer :: i, inum, ich, isp, iup, ilo
  character(len=10) :: ssy, ssy2
  character*26 :: lowcha, upchar
  character*10 :: numcha
  character*2  :: sigcha, symato(0:maxato), sym2
  data lowcha /'abcdefghijklmnopqrstuvwxyz'/
  data upchar /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
  data numcha /'0123456789'/
  data sigcha /'+-'/
  data symato /'L ', &
     &         'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
     &         'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
     &         'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
     &         'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
     &         'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
     &         'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
     &         'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
     &         'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
     &         'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
     &         'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
     &         'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
     &         'Rg','Cn','Uu','Fl','Up','Lv' /
  atnum = -1
  atoxs = 0.0d+0
  inum = 0
  ssy = trim(adjustl(assymb))
  ! check for reasonable input string
  if (0==len_trim(ssy)) return ! empty string
  ! check for occurance of a number on the input string
  inum = SCAN(ssy,numcha) ! determine the position of the first numerical character
  ! get the pure element symbol part
  if (inum>1) then
    ssy2 = ssy(1:inum-1)
  else
    ssy2 = trim(ssy)
  end if
  sym2 = '  '
  iup = scan(upchar,ssy2(1:1)) ! determine the upper-case character type of the first symbol character
  ilo = scan(lowcha,ssy2(1:1)) ! determine the lower-case character type of the first symbol character
  ich = max(iup,ilo)
  if (ich<=0) return ! invalid symbol, do not analyse further
  sym2(1:1) = upchar(ich:ich)
  iup = scan(upchar,ssy2(2:2)) ! determine the upper-case character type of the second symbol character
  ilo = scan(lowcha,ssy2(2:2)) ! determine the lower-case character type of the second symbol character
  ich = max(iup,ilo)
  if (ich>0) sym2(2:2) = lowcha(ich:ich)
  ! sym2 is now a "Xx" type of string for comparison to symato
  do i=0, maxato ! find the atomic number (atnum remains 0 if not found)
    if (trim(sym2)==trim(symato(i))) then
      atnum = i
      exit
    end if
  end do
  if (atnum<0) return ! invalid symbol, do not analyse further
  isp = SCAN(ssy,sigcha) ! determine the position of the sign character
  if (isp>0) ich = SCAN(sigcha,ssy(isp:isp)) ! determine which sign character is there
  if (inum>1.and.ich>=1.and.isp>inum) then ! there is an oxydation number in the symbol string
    ssy2 = ssy(inum:isp-1)
    read(unit=ssy2,fmt=*,err=801) iup
    if (ich==1) atoxs =  dble(iup)
    if (ich==2) atoxs = -dble(iup)
  end if
  goto 1000
801 continue
  call CIF_ERR("Failed to read oxidation state number from symbol string ("//&
     & trim(ssy)//").")
  goto 1000
1000 continue
END SUBROUTINE CIF_get_atnum


!********************************************************
! CIF_get_atom_site_type_symbol
! This routine re-creates the _atom_site_type_symbol
! string from the current data. Whenever the charge is
! non-zero, a respective oxidation string is added to
! the element symbol (rounded to the neares integer).
! Returns an empty string in case of invalid data.
!********************************************************
SUBROUTINE CIF_get_atom_site_type_symbol(iat,satsym)
  implicit none
  integer, intent(in) :: iat
  character(len=*), intent(out) :: satsym
  integer :: NP, natnum, natcrg, natcrg2
  integer, parameter :: maxato = 116
  character*2  :: symato(0:maxato), symcrg
  data symato /'L ', &
     &         'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
     &         'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
     &         'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
     &         'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
     &         'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
     &         'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
     &         'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
     &         'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
     &         'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
     &         'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
     &         'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
     &         'Rg','Cn','Uu','Fl','Up','Lv' /
  !   
  satsym = ""
  NP = 0
  if (.not.allocated(CIF_atom_site)) return ! no data
  NP = SIZE(CIF_atom_site,2)
  if (iat<1.or.iat>CIF_atom_site_number.or.iat>NP) return ! invalid index
  natnum = nint(CIF_atom_site(4,iat))
  natcrg = nint(CIF_atom_site(5,iat))
  if (natnum<0.or.natnum>maxato) return
  if (natcrg/=0) then
    symcrg = ""
    natcrg2 = abs(natcrg)
    if (natcrg>0) then
      write(unit=symcrg,fmt='(I1,"+")') natcrg2
    else
      write(unit=symcrg,fmt='(I1,"-")') natcrg2
    end if
    satsym = trim(symato(natnum))//trim(symcrg)
  else
    satsym = trim(symato(natnum))
  end if
END SUBROUTINE CIF_get_atom_site_type_symbol


!********************************************************
! CIF_get_atom_type_symbol
! This routine re-creates the _atom_site_type_symbol
! string from the atomic number and oxidation state.
! The oxidation string is added to
! the element symbol (rounded to the neares integer).
! Returns an empty string in case of invalid data.
!********************************************************
SUBROUTINE CIF_get_atom_type_symbol(atnum,charge,satsym)
  implicit none
  integer, intent(in) :: atnum
  real(dp), intent(in) :: charge
  character(len=*), intent(out) :: satsym
  integer :: natnum, natcrg, natcrg2
  integer, parameter :: maxato = 116
  character*2  :: symato(0:maxato), symcrg
  data symato /'L ', &
     &         'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
     &         'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
     &         'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
     &         'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
     &         'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
     &         'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
     &         'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
     &         'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
     &         'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
     &         'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
     &         'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
     &         'Rg','Cn','Uu','Fl','Up','Lv' /
  !   
  satsym = ""
  natnum = atnum
  natcrg = nint(charge)
  if (natnum<0.or.natnum>maxato) return
  if (natcrg/=0) then
    symcrg = ""
    natcrg2 = abs(natcrg)
    if (natcrg>0) then
      write(unit=symcrg,fmt='(I1,"+")') natcrg2
    else
      write(unit=symcrg,fmt='(I1,"-")') natcrg2
    end if
    satsym = trim(symato(natnum))//trim(symcrg)
  else
    satsym = trim(symato(natnum))
  end if
END SUBROUTINE CIF_get_atom_type_symbol


!********************************************************
! CIF_get_cell_volume
! This routine calculates the unit cell volume from the
! cell a, b, c, alpha, beta, and gamma parameters.
! A zero volume indicates an undefined structure.
!********************************************************
SUBROUTINE CIF_get_cell_volume(vol)
  implicit none
  real(dp), parameter :: dd2r = 1.7453292519943295769237d-2 
  real(dp), intent(out) :: vol
  real(dp) :: ca, cb, cc, da
  vol = 0.d+0
  ca = dcos(CIF_cell_angle_alpha*dd2r)
  cb = dcos(CIF_cell_angle_beta *dd2r)
  cc = dcos(CIF_cell_angle_gamma*dd2r)
  da = dsqrt(dabs(1.d+0-ca*ca-cb*cb-cc*cc+2.d+0*ca*cb*cc))
  vol = CIF_cell_length_a*CIF_cell_length_b * &
      & CIF_cell_length_c* da
END SUBROUTINE CIF_get_cell_volume


!********************************************************
! CIF_get_cell_reciprocal
! This routine calculates the reciprocal unit cell
! parameters from the cell a, b, c, alpha, beta, and
! gamma parameters.
! Calculations accordint to
! Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
!      New York: John Wiley & Sons Inc.
!********************************************************
SUBROUTINE CIF_get_cell_reciprocal(ra, rb, rc, ralpha, rbeta, rgamma)
  implicit none
  real(dp), parameter :: dd2r = 1.7453292519943295769237d-2 
  real(dp), intent(out) :: ra, rb, rc, ralpha, rbeta, rgamma
  real(dp) :: vol
  real(dp) :: ca, cb, cc, sa, sb, sc
  real(dp) :: cralpha, crbeta, crgamma
  ra = 0.d+0
  rb = 0.d+0
  rc = 0.d+0
  ralpha = 9.d+1
  rbeta  = 9.d+1
  rgamma = 9.d+1
  call CIF_get_cell_volume(vol)
  if (vol>0.d+0) then
    sa = dsin(CIF_cell_angle_alpha*dd2r)
    sb = dsin(CIF_cell_angle_beta *dd2r)
    sc = dsin(CIF_cell_angle_gamma*dd2r)
    ca = dcos(CIF_cell_angle_alpha*dd2r)
    cb = dcos(CIF_cell_angle_beta *dd2r)
    cc = dcos(CIF_cell_angle_gamma*dd2r)
    ra = CIF_cell_length_b*CIF_cell_length_c * sa / vol
    rb = CIF_cell_length_a*CIF_cell_length_c * sb / vol
    rc = CIF_cell_length_a*CIF_cell_length_b * sc / vol
    cralpha = (cb*cc - ca) / (sb*sc)
    crbeta  = (cc*ca - cb) / (sc*sa)
    crgamma = (ca*cb - cc) / (sa*sb)
    ralpha = dacos(cralpha) / dd2r
    rbeta  = dacos(crbeta ) / dd2r
    rgamma = dacos(crgamma) / dd2r
  end if
END SUBROUTINE CIF_get_cell_reciprocal


!********************************************************
! CIF_get_iso_from_ani
! This routine determines the isotropic or equivalent
! thermal displacement parameter from the current
! anisotropic thermal displacement parameters of the
! requested atomic site
!********************************************************
SUBROUTINE CIF_get_iso_from_ani(iat,viso)
  implicit none
  integer, intent(in) :: iat
  real(dp), intent(out) :: viso
  integer :: NP, adp_type, i, j
  real(dp) :: tmp, a0(3), a1(3), aa1(3), adp(3,3)
  viso = 0.d+0
  a1 = 0.d+0
  aa1 = 0.d+0
  a0 = (/ CIF_cell_length_a, CIF_cell_length_b, CIF_cell_length_c /)
  call CIF_get_cell_reciprocal(a1(1),a1(2),a1(3),aa1(1),aa1(2),aa1(3))
  adp(1:3,1) = (/ CIF_atom_site( 9,iat), CIF_atom_site(12,iat), &
                & CIF_atom_site(13,iat) /)
  adp(1:3,2) = (/ CIF_atom_site(12,iat), CIF_atom_site(10,iat), &
                & CIF_atom_site(14,iat) /)
  adp(1:3,3) = (/ CIF_atom_site(13,iat), CIF_atom_site(14,iat), &
                & CIF_atom_site(11,iat) /)
  tmp = 0.d+0
  do j=1,3
    do i=1,3
      tmp = tmp + adp(i,j)*a0(i)*a0(j)*a1(i)*a1(j)
    end do
  end do
  viso = tmp / 3.d+0
END SUBROUTINE CIF_get_iso_from_ani


!********************************************************
! CIF_get_atom_site_biso
! This routine returns the Biso value in A^2 units for
! the specified atom site (by index).
!********************************************************
SUBROUTINE CIF_get_atom_site_biso(iat,biso)
  implicit none
  integer, intent(in) :: iat
  real(dp), intent(out) :: biso
  integer :: NP, adp_type, i, j
  real(dp) :: tmp
  biso = 0.0d+0
  tmp = 0.0d+0
  NP = 0
  if (.not.allocated(CIF_atom_site)) return ! no data
  NP = SIZE(CIF_atom_site,2)
  if (iat<1.or.iat>CIF_atom_site_number.or.iat>NP) return ! invalid index
  adp_type = nint(CIF_atom_site(7,iat))
  select case (adp_type)
  case (1) ! Uiso
    biso = CIF_atom_site(8,iat)*CIF_8pi2 ! translate from Uiso to Biso
  case (2) ! Biso
    biso = CIF_atom_site(8,iat) ! transfer as is
  case (3) ! Uani
    ! calculate from Uani to Uiso (tmp)
    call CIF_get_iso_from_ani(iat,tmp)
    ! calculate from Uiso to Biso
    biso = tmp*CIF_8pi2
  case (4) ! Bani
    ! calculate from Bani to Biso (tmp)
    call CIF_get_iso_from_ani(iat,tmp)
    biso = tmp
  case default ! assuming Uiso or Uani
    biso = CIF_atom_site(8,iat)*CIF_8pi2 ! translate from Uiso to Biso
    if (biso <= 0.d+0) then ! try Uani to Biso
      call CIF_get_iso_from_ani(iat,tmp)
      biso = tmp*CIF_8pi2
    end if
    if (biso <= 0.d+0) biso = 0.d+0 ! reset to zero
  end select ! (adp_type)
END SUBROUTINE CIF_get_atom_site_biso


!********************************************************
! CIF_get_atom_site_uiso
! This routine returns the Uiso value in A^2 units for
! the specified atom site (by index).
!********************************************************
SUBROUTINE CIF_get_atom_site_uiso(iat,uiso)
  implicit none
  integer, intent(in) :: iat
  real(dp), intent(out) :: uiso
  integer :: NP, adp_type, i, j
  real(dp) :: tmp
  uiso = 0.0d+0
  tmp = 0.0d+0
  NP = 0
  if (.not.allocated(CIF_atom_site)) return ! no data
  NP = SIZE(CIF_atom_site,2)
  if (iat<1.or.iat>CIF_atom_site_number.or.iat>NP) return ! invalid index
  adp_type = nint(CIF_atom_site(7,iat))
  select case (adp_type)
  case (1) ! Uiso
    uiso = CIF_atom_site(8,iat) ! transfer as is
  case (2) ! Biso
    uiso = CIF_atom_site(8,iat)/CIF_8pi2 ! translate from Biso to Uiso
  case (3) ! Uani
    ! calculate from Uani to Uiso (tmp)
    call CIF_get_iso_from_ani(iat,tmp)
    uiso = tmp
  case (4) ! Bani
    ! calculate from Bani to Biso (tmp)
    call CIF_get_iso_from_ani(iat,tmp)
    ! calculate from Biso to Uiso
    uiso = tmp / CIF_8pi2
  case default ! assuming Uiso or Uani
    uiso = CIF_atom_site(8,iat) ! transfer as is
    if (uiso <= 0.d+0) then ! try Uani to Biso
      call CIF_get_iso_from_ani(iat,tmp)
      uiso = tmp
    end if
    if (uiso <= 0.d+0) uiso = 0.d+0 ! reset to zero
  end select ! (adp_type)
END SUBROUTINE CIF_get_atom_site_uiso


!********************************************************
! CIF_find_atom_site_name
! This subroutine returns the index of the parameter name
! in the _atom_site table of this module. It returns 0
! if the name is not in the list
!********************************************************
SUBROUTINE CIF_find_atom_site_name(sname,idx,jdx)
  implicit none
  character(len=*), intent(in) :: sname
  integer, intent(inout) :: idx
  integer, optional, intent(inout) :: jdx
  integer :: i, j
  idx = 0
  if (len_trim(sname)>0) then
    do i=1, CIF_atom_site_items_num
      do j=1, CIF_atom_site_namealts
        if (1==index(CIF_atom_site_items(j,i),trim(adjustl(sname)))) then
          idx = i
          if (present(jdx)) jdx = j
          return ! found, done
        end if
      end do
    end do
  end if
  ! not found, should return 0
END SUBROUTINE CIF_find_atom_site_name


!********************************************************
! CIF_read_item
! This routine reads a name and a value from the input
! string and writes the value to a module variable if
! the name is recognized. The routine assumes that the
! input string contains at least one pair of name and
! value. It will remove the name and the value from the
! input string, leaving possible following items for the
! next call or an empty string.
! The input/output string variable is of fix length=lsl.
!********************************************************
SUBROUTINE CIF_read_item(sline)
  implicit none
  character(len=lsl), intent(inout) :: sline
  character(len=ssl) :: snam, sval
  real(dp) :: dval
  integer :: i, ioerr, ilen, ival, j, nchk
  ! initialize
  ioerr = 0
  snam = ''
  sval = ''
  dval = 0.0d+0
  ! try to read the next name and value
  read (unit=sline,fmt=*,iostat=ioerr)  snam, sval
  if (ioerr/=0) then ! failed
    call CIF_ERR("Failed to read next data name and value "//&
       & "in '"//trim(sline)//"', skipping.", ioerr)
    sline = ""
    goto 1000
  end if
  ilen = len_trim(sline)
  ival = index(sline,trim(sval))
  ! set the value to the variable by name
  if     ("_cell_length_a"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_cell_length_a = dval
  elseif ("_cell_length_b"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_cell_length_b = dval
  elseif ("_cell_length_c"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_cell_length_c = dval
  elseif ("_cell_angle_alpha"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_cell_angle_alpha = dval
  elseif ("_cell_angle_beta"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_cell_angle_beta = dval
  elseif ("_cell_angle_gamma"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_cell_angle_gamma = dval
  elseif ("_symmetry_space_group_name_H-M"==trim(snam) .or. &
       &  "_symmetry_space_group_name_H-M_alt"==trim(snam) ) then
    CIF_symmetry_space_group_name_HM = ''
    j = 0
    do i=1, len_trim(sval)
      ! if (sval(i:i)==" ") cycle ! remove space character
      if (sval(i:i)=="'") cycle ! remove quotation marks
      if (sval(i:i)=='"') cycle ! remove quotation marks
      j = j + 1
      CIF_symmetry_space_group_name_HM(j:j) = sval(i:i)
    end do
    CIF_symmetry_space_group_name_HM(j+1:j+1) = achar(0) ! terminate the string
  elseif ("_symmetry_Int_Tables_number"==trim(snam)) then
    call CIF_get_value(sval,dval,nchk)
    if (nchk/=0) goto 800
    CIF_symmetry_Int_Tables_number = nint(dval)
  end if
  !
  goto 900
  !
800 continue
  call CIF_ERR("Failed to interprete the value of "// &
     & trim(snam)//" from '"//trim(sval)//"', skipping.", nchk)
  goto 900
  !
900 continue
  ! remove the name and value string from sline
  call CIF_str_remsub(sline,adjustl(trim(snam)))
  sline = trim(adjustl(sline))
  call CIF_str_remsub(sline,adjustl(trim(sval)))
  sline = trim(adjustl(sline))
  !
1000 continue
  !
END SUBROUTINE CIF_read_item


!********************************************************
! CIF_get_loop_stats
! This routine determines the table size of a CIF loop_
! construct. It doesn't interprete anything.
! The interface of the routine requires a start line
! number iline and the remaining string in this line
! the string must begin with "loop_".
! When returning iline is changed to the final line of
! the table reading and sline contains the remaining
! unprocessed information in this line.
! Tries to determine the type of data in the loop (cid):
!  (1) : _atom_site data
!  (2) : secondary _atom_site data
!  (3) : _symmetry data
!  (4) : unknown category
!********************************************************
SUBROUTINE CIF_get_loop_stats(iline,sline,lstats)
  implicit none
  integer*4, intent(inout) :: iline
  character(len=lsl), intent(inout) :: sline
  integer*4, intent(out) :: lstats(CIF_lstat_len)
  integer*4 :: ncol, nrow, cid
  integer*4 :: i, j, ml(1)
  integer*4 :: ilen, irow
  integer*4 :: nlines
  integer*4 :: nitems
  integer*4 :: cidvotes(3)
  character(len=ssl) :: sname
  !
  ! Initialize
  nlines = SIZE(CIF_lines)
  lstats = 0
  ncol = 0
  nrow = 0
  ilen = 0
  cidvotes = 0
  sname = ""
  if (iline>nlines) return ! reached EOF, strange, do nothing.
  lstats(1) = iline
  !
  ! The current sline contains "loop_". Remove this!
  ilen = len_trim(sline)
  if (ilen<=5) then
    sline = ""
  else
    sline = adjustl(sline(6:ilen))
  end if
  !
  ! <<< HEADER ANALYSIS >>>
  !
  do ! loop header parsing
    !
    ! further analyse the content of sline
    ! determine the number of data items per table row = ncol
    !
    if (len_trim(sline)<1) goto 20 ! skip empty lines
    !
    ! content parser
    sname = ""
    read (unit=sline, fmt=*) sname ! get the name of the item
    if (sname(1:1)=="_") then ! this is a data name
      ! count
      ncol = ncol+1
      if (0<index(sname,"_atom_site")) cidvotes(2) = cidvotes(2) + 1
      if (0<index(sname,"_atom_site_fract")) cidvotes(1) = cidvotes(1) + 1
      if (0<index(sname,"_symmetry_equiv_pos") .or. 0<index(sname,"_space_group_symop")) cidvotes(3) = cidvotes(3) + 1
    else ! this is no data name 
      ! exit the counting loop, leave sline and iline as is
      exit
    end if
    !
    ! we are here, becaus a name was found, remove it now from sline
    call CIF_str_remsub(sline,adjustl(trim(sname)))
    sline = trim(adjustl(sline))
    cycle ! try next content of sline
    !
20  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! header parsing // counting ncol
  !
  ! determine the type of table from cidvotes
  if (0==maxval(cidvotes)) then
    cid = 4 ! unkown content
  else
    ml = maxloc(cidvotes)
    cid = ml(1) ! this will be either 1, 2, or 3, probably 2 or 3
    ! in case of 2, look for votes in 1, if at least 3 votes exist, this is the primary atom_site table cid=1
    if (cid==2.and.cidvotes(1)>2) cid = 1 
  end if
  !
  ! <<< DATA COUNTING >>>
  !
  nitems = 0
  do ! data counting loop
    !
    ! Read items, check if they can be data. If not, stop!
    ! Fill up rows and count.
    !
    if (len_trim(sline)<1) goto 30 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*) sname ! get the item
    if (sname(1:1)=="_") then ! this is a data name
      exit ! stop counting items
    elseif (sname(1:5)=="loop_") then ! this is another loop
      exit ! stop counting items
    elseif (sname(1:5)=="data_" .or. sname(1:5)=="save_") then ! this is other data
      exit ! stop counting items
    else ! this is no data name 
      nitems = nitems + 1 ! could be an item
    end if
    ! remove name from sline
    call CIF_str_remsub(sline,adjustl(trim(sname)))
    sline = trim(adjustl(sline))
    cycle ! next item
    !
30  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! data counting loop
  !
  ! nitems is the number of items
  nrow = (nitems-modulo(nitems,ncol))/ncol ! count full data rows
  ! the current line #iline already contains content that doesn't belong to the loop
  ! thus reduce the line index by one
  iline = iline -1
  !
  ! set the loop information to the interface output array
  lstats(2) = iline
  lstats(3) = ncol
  lstats(4) = nrow
  lstats(5) = cid
  !
END SUBROUTINE CIF_get_loop_stats


!********************************************************
! CIF_read_loop1
! This routine reads a table of data following a "loop_"
! item in a CIF file.
! The interface of the routine requires a list of values
! as stored in CIF_lstats, including the first and the
! last line of the table data.
! The data is interpreted as primary CIF atomic site
! table and used to fill the modules CIF_atom_site table.
!********************************************************
SUBROUTINE CIF_read_loop1(lstats,nchk)
  implicit none
  integer, intent(in) :: lstats(CIF_lstat_len)
  integer, intent(out) :: nchk
  integer :: ioerr
  integer :: iline
  integer :: nlines
  integer :: idx, jdx
  integer :: icol, irow
  integer :: ncol, nrow
  integer :: ipos, ilen
  integer :: nalloc
  integer :: ichk
  integer :: ival
  integer :: idtype
  integer, allocatable :: idxhdr(:,:)
  real(dp) :: dval
  character(len=ssl) :: snam, sval
  character(len=lsl) :: sline, smsg
  !
  ! init
  nchk = 0
  ioerr = 0
  ichk = 0
  nalloc = 0
  nlines = SIZE(CIF_lines)
  iline = lstats(1)
  sline = CIF_lines(iline)
  idx = 0
  jdx = 0
  idtype = 0
  ival = 0
  dval = 0.0d+0
  icol = 0
  ncol = lstats(3)
  irow = 0
  nrow = lstats(4)
  if (ncol<1.or.nrow<1) return ! table contains no data
  allocate(idxhdr(2,ncol),stat=nalloc) ! allocate the header hash table
  if (nalloc/=0) goto 801
  idxhdr = 0
  !
  ! The current sline contains (it should) "loop_". Remove this!
  ipos = index(sline,"loop_")
  ilen = len_trim(sline)
  if (ipos>0.and.ilen>5) then
    sline = adjustl(sline(6:ilen))
  else
    sline = ""
  end if
  !
  !
  !
  ! <<< HEADER ANALYSIS >>>
  !
  ! - find out, which type of data is in which column of the table
  !
  do ! loop header parsing
    !
    ! further analyse the content of sline
    ! determine the number of data items per table row = ncol
    !
    if (len_trim(sline)<1) goto 20 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*) snam ! get the name of the item
    if (snam(1:1)=="_") then ! this is a data name
      ! count
      icol = icol+1
      if (icol>ncol) exit ! max. number of columns reached, just exit the counting loop
      call CIF_find_atom_site_name(snam, idx, jdx) ! get the location index in the module table
      idxhdr(1,icol) = idx ! store the index
      idxhdr(2,icol) = jdx ! store the alternative index
    else ! this is no data name 
      ! exit the counting loop
      exit
    end if
    !
    ! we are here, because a name was found, remove it now from sline
    call CIF_str_remsub(sline,adjustl(trim(snam)))
    sline = trim(adjustl(sline))
    cycle ! try next content of sline
    !
20  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! header parsing
  !
  if (icol==0) then
    write(unit=smsg,fmt=*) iline
    call CIF_ERR("Parsing table header failed in line #"// &
       & trim(adjustl(smsg))//".")
    goto 1000 ! corrupt table, do nothing more
  end if
  !
  !
  ! <<< Data Reading >>>
  !
25 continue
  icol = 1 ! current column
  irow = 1 ! current row
  do ! data counting loop
    !
    ! Read items, check if they can be data. If not, stop!
    ! Fill up the table.
    !
    if (len_trim(sline)<1) goto 30 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*) sval ! get the item
    if (sval(1:1)=="_") then ! this is a data name
      exit ! stop reading data
    elseif (sval(1:5)=="loop_") then ! this is another loop
      exit ! stop reading data
    elseif (sval(1:5)=="data_" .or. sval(1:5)=="save_") then ! this is other data
      exit ! stop reading data
    else ! this can be data
      !
      ! interprete the data depending on the column type
      idx = idxhdr(1,icol)
      jdx = idxhdr(2,icol)
      ! and the respective data type
      if (idx>0) then ! this is known data, try to get it
        idtype = CIF_atom_site_dtype(idx)
        select case (idtype)
        case (0) ! string transferring
          if (idx==4) then ! symbol interpretation to atomic number
            call CIF_get_atnum(trim(sval),ival,dval)
            CIF_atom_site(idx,CIF_atom_site_number+1) = dble(ival)
            CIF_atom_site(idx+1,CIF_atom_site_number+1) = dval
          elseif (idx==7) then ! explicit thermal displacement type interpretation
            CIF_atom_site(idx,CIF_atom_site_number+1) = 0.0d+0
            if     (trim(sval)=="Uiso") then
              CIF_atom_site(idx,CIF_atom_site_number+1) = 1.0d+0
            elseif (trim(sval)=="Biso") then
              CIF_atom_site(idx,CIF_atom_site_number+1) = 2.0d+0
            elseif (trim(sval)=="Uani") then
              CIF_atom_site(idx,CIF_atom_site_number+1) = 3.0d+0
            elseif (trim(sval)=="Bani") then
              CIF_atom_site(idx,CIF_atom_site_number+1) = 4.0d+0
            end if
          end if
        case (1) ! integer reading
          call CIF_get_value(sval,dval,ichk)
          if (ichk/=0) goto 800
          CIF_atom_site(idx,CIF_atom_site_number+1) = dval
        case (2) ! float reading
          call CIF_get_value(sval,dval,ichk)
          if (ichk/=0) goto 800
          CIF_atom_site(idx,CIF_atom_site_number+1) = dval
          if (idx==8 .and. nint(CIF_atom_site(7,CIF_atom_site_number+1))<1 ) then
            CIF_atom_site(7,CIF_atom_site_number+1) = dble(jdx) ! in case of U_iso or B_iso, preset the adp type depending on the alternative index
          end if
        end select
      end if
      !
    end if
    ! remove value string from sline
    call CIF_str_remsub(sline,adjustl(trim(sval)))
    sline = trim(adjustl(sline))
    icol = icol + 1 ! advance column index
    if (icol>ncol) then ! last column reached
      icol = 1 ! move back to columns one
      CIF_atom_site_number = CIF_atom_site_number + 1 ! increase atom site table content count
      irow = irow + 1 ! advance row index
    end if
    if (irow>nrow) exit ! stop reading data, table is complete
    cycle ! next item
    !
30  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! data counting loop
  !
35 continue
  !
  goto 1000
  !
  !
  ! <<< Error Handling >>>
  !
800 continue
  call CIF_ERR("Failed to interprete the value from '"//trim(sval) &
     & //"', skipping table.", nchk)
  goto 1000
801 continue
  nchk = 2
  call CIF_ERR("Failed to allocate memory.", nalloc)
  goto 1000
802 continue
  nchk = 3
  call CIF_ERR("Failed to read from string '"//trim(sline)// &
     & "'.", ioerr)
  goto 1000
  !
  ! <<< Exit >>>
  !
1000 continue
  if (allocated(idxhdr)) deallocate(idxhdr, stat=nalloc)
  return
  !
END SUBROUTINE CIF_read_loop1


!********************************************************
! CIF_read_loop2
! This routine reads a table of data following a "loop_"
! item in a CIF file.
! The interface of the routine requires a list of values
! as stored in CIF_lstats, including the first and the
! last line of the table data.
! The loop is interpreted as secondary loop, applying
! the table data to all atomic sites of the current
! module table CIF_atom_site. Table rows will be selected
! by the _atom_site_symbol string.
!********************************************************
SUBROUTINE CIF_read_loop2(lstats, nchk)
  implicit none
  integer, intent(in) :: lstats(CIF_lstat_len)
  integer, intent(out):: nchk
  integer :: ioerr
  integer :: iline
  integer :: nlines
  integer :: idx, jdx, i, j
  integer :: icol, irow
  integer :: ncol, nrow
  integer :: ipos, ilen
  integer :: nalloc
  integer :: ichk
  integer :: ival
  integer :: idtype
  integer, allocatable :: idxhdr(:,:)
  real(dp), allocatable :: l_atom_site(:)
  real(dp) :: dval
  character(len=ssl) :: snam, sval
  character(len=lsl) :: sline, smsg
  !
  ! init
  nchk = 0
  ioerr = 0
  ichk = 0
  nalloc = 0
  nlines = SIZE(CIF_lines)
  iline = lstats(1)
  sline = CIF_lines(iline)
  idx = 0
  idtype = 0
  ival = 0
  dval = 0.0d+0
  icol = 0
  ncol = lstats(3)
  irow = 0
  nrow = lstats(4)
  if (ncol<1.or.nrow<1) return ! table contains no data
  allocate(idxhdr(2,ncol),stat=nalloc) ! allocate the header hash table
  if (nalloc/=0) goto 801
  idxhdr = 0
  allocate(l_atom_site(CIF_atom_site_items_num),stat=nalloc) ! allocate the local table row
  if (nalloc/=0) goto 801
  l_atom_site = 0.0d+0
  !
  ! The current sline contains (it should) "loop_". Remove this!
  ipos = index(sline,"loop_")
  ilen = len_trim(sline)
  if (ipos>0.and.ilen>5) then
    sline = adjustl(sline(6:ilen))
  else
    sline = "" ! though strange, "loop_" should be there
  end if
  !
  !
  !
  ! <<< HEADER ANALYSIS >>>
  !
  ! - find out, which type of data is in which column of the table
  !
  do ! loop header parsing
    !
    ! further analyse the content of sline
    ! determine the number of data items per table row = ncol
    !
    if (len_trim(sline)<1) goto 20 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*, iostat=ioerr, err=802) snam ! get the name of the item
    if (snam(1:1)=="_") then ! this is a data name
      ! count
      icol = icol+1
      if (icol>ncol) exit ! max. number of columns reached, just exit the counting loop
      call CIF_find_atom_site_name(snam, idx, jdx) ! get the location index in the module table
      idxhdr(1,icol) = idx ! store the index
      idxhdr(2,icol) = jdx ! store the alternative index
      if (idx==4) ival = icol ! remember column of the atom symbol item
    else ! this is no data name 
      ! exit the counting loop
      exit
    end if
    !
    ! we are here, because a name was found, remove it now from sline
    call CIF_str_remsub(sline,adjustl(trim(snam)))
    sline = trim(adjustl(sline))
    cycle ! try next content of sline
    !
20  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! header parsing
  !
  if (icol==0) then
    write(unit=smsg,fmt=*) iline
    call CIF_ERR("Parsing table header failed in line #"// &
       & trim(adjustl(smsg))//".")
    goto 1000 ! corrupt table, do nothing more
  end if
  !
  if (ival==0) then
    write(unit=smsg,fmt=*) lstats(1)
    call CIF_WRN("Missing item _atom_site_type_symbol in table of line #"// &
       & trim(adjustl(smsg))//".")
    call CIF_WRN("Failed to link to primary atom site data.")
    goto 1000 ! invalid table, do nothing more
  end if
  !
  !
  ! <<< Data Reading >>>
  !
25 continue
  icol = 1 ! current column
  irow = 1 ! current row
  l_atom_site = 0.0d+0
  do ! data reading loop
    !
    ! Read items, check if they can be data. If not, stop!
    ! Fill up the table.
    !
    if (len_trim(sline)<1) goto 30 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*, iostat=ioerr, err=802) sval ! get the item
    if (sval(1:1)=="_") then ! this is a data name
      exit ! stop reading data
    elseif (sval(1:5)=="loop_") then ! this is another loop
      exit ! stop reading data
    elseif (sval(1:5)=="data_" .or. sval(1:5)=="save_") then ! this is other data
      exit ! stop reading data
    else ! this can be data
      !
      ! interprete the data depending on the column type
      idx = idxhdr(1,icol)
      jdx = idxhdr(2,icol)
      ! and the respective data type
      if (idx>0) then ! this is known data, try to get it
        idtype = CIF_atom_site_dtype(idx)
        select case (idtype)
        case (0) ! string transferring
          if (idx==4) then ! symbol interpretation to atomic number
            call CIF_get_atnum(trim(sval),ival,dval)
            l_atom_site(idx) = dble(ival)
            l_atom_site(idx+1) = dval
          elseif (idx==7) then ! explicit thermal displacement type interpretation
            l_atom_site(idx) = 0.0d+0 ! preset with "not specified"
            if     (trim(sval)=="Uiso") then
              l_atom_site(idx) = 1.0d+0
            elseif (trim(sval)=="Biso") then
              l_atom_site(idx) = 2.0d+0
            elseif (trim(sval)=="Uani") then
              l_atom_site(idx) = 3.0d+0
            elseif (trim(sval)=="Bani") then
              l_atom_site(idx) = 4.0d+0
            end if
          end if
        case (1) ! integer reading
          call CIF_get_value(sval,dval,ichk)
          if (ichk/=0) goto 800
          l_atom_site(idx) = dval
        case (2) ! float reading
          call CIF_get_value(sval,dval,ichk)
          if (ichk/=0) goto 800
          l_atom_site(idx) = dval
          if (idx>=9.and.idx<=14) then ! Uani or Bani
            ! override the current adp_type index with the superior data
            l_atom_site(7) = dble(2+jdx) ! 3 = Uani, 4 = Bani
          end if
        end select
      end if
      !
    end if
    ! remove value string from sline
    call CIF_str_remsub(sline,adjustl(trim(sval)))
    sline = trim(adjustl(sline))
    icol = icol + 1 ! advance column index
    if (icol>ncol) then ! last column reached
      icol = 1 ! move back to loop column one
      irow = irow + 1 ! advance row index
      ! transfer current local table row to the main table
      if (CIF_atom_site_number>0) then ! there is table content
        do i=1, CIF_atom_site_number ! loop primary atomic site table
          if (       0.1d+0>dabs(CIF_atom_site(4,i)-l_atom_site(4)) & ! same atomic number
             & .and. 0.1d+0>dabs(CIF_atom_site(5,i)-l_atom_site(5)) & ! same oxidation state or charge
             & ) then
            ! same atom type / species -> transfer secondary data
            do j=1, ncol ! loop local row
              idx = idxhdr(1,j)
              if (idx>5) then ! only transfer / override auxiliary data
                CIF_atom_site(idx,i) = l_atom_site(idx)
              end if
              if (idx>=9.and.idx<=14) then ! Uani or Bani ! do this again, since idx=7 is usually not in idxhdr
                ! override the current adp_type index with the superior data
                 CIF_atom_site(7,i) = l_atom_site(7) ! 3 = Uani, 4 = Bani
              end if
            end do ! loop local row
          end if !
        end do ! loop primary atomic site table
      end if !
      l_atom_site = 0.0d+0 ! prepare to read next table row
    end if
    if (irow>nrow) exit ! stop reading data, table is complete
    cycle ! next item
    !
30  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! data reading loop
  !
35 continue
  !
  goto 1000
  !
  !
  ! <<< Error Handling >>>
  !
800 continue
  nchk = 1
  call CIF_ERR("Failed to interprete the value from '"//trim(sval) &
     & //"', skipping table.", ichk)
  goto 1000
801 continue
  nchk = 2
  call CIF_ERR("Failed to allocate memory.", nalloc)
  goto 1000
802 continue
  nchk = 3
  call CIF_ERR("Failed to read from string '"//trim(sline)// &
     & "'.", ioerr)
  goto 1000
  !
  ! <<< Exit >>>
  !
1000 continue
  if (allocated(l_atom_site)) deallocate(l_atom_site, stat=nalloc)
  if (allocated(idxhdr)) deallocate(idxhdr, stat=nalloc)
  return
  !
END SUBROUTINE CIF_read_loop2


!********************************************************
! CIF_read_loop3
! This routine reads a table of data following a "loop_"
! item in a CIF file.
! The interface of the routine requires a list of values
! as stored in CIF_lstats, including the first and the
! last line of the table data.
! The loop is interpreted as symmetry loop, applying
! the symmetry operations given as string by the CIF
! item name "_symmetry_equiv_pos_as_xyz" to all atomic
! sites of the current module table CIF_atom_site.
!********************************************************
SUBROUTINE CIF_read_loop3(lstats, nchk)
  implicit none
  integer, intent(in) :: lstats(CIF_lstat_len)
  integer, intent(out):: nchk
  integer :: ioerr
  integer :: iline
  integer :: nlines
  integer :: idx, i, j, idxid, idxop
  integer :: icol, irow
  integer :: ncol, nrow
  integer :: ipos, ilen
  integer :: nalloc
  integer :: ichk
  integer :: ival
  integer :: idtype
  integer :: symop_id
  character(len=ssl) :: snam, sval, symop_xyz
  character(len=lsl) :: sline, smsg
  !
  ! init
  nchk = 0
  ioerr = 0
  ichk = 0
  nalloc = 0
  nlines = SIZE(CIF_lines)
  iline = lstats(1)
  sline = CIF_lines(iline)
  idx = 0
  idxid = 0
  idxop = 0
  idtype = 0
  ival = 0
  icol = 0
  ncol = lstats(3)
  irow = 0
  nrow = lstats(4)
  symop_id = 0
  symop_xyz = ""
  if (ncol<1.or.nrow<1) return ! table contains no data
  !
  ! Prepare the module handling the symmetry operations
  if (allocated(symops_trf)) deallocate(symops_trf, stat=nalloc)
  allocate(symops_trf(symops_nltrf,nrow), stat=nalloc)
  if (nalloc/=0) goto 801
  call SYMOPS_INIT()
  !
  ! The current sline contains (it should) "loop_". Remove this!
  ipos = index(sline,"loop_")
  ilen = len_trim(sline)
  if (ipos>0.and.ilen>5) then
    sline = adjustl(sline(6:ilen))
  else
    sline = "" ! though strange, "loop_" should be there
  end if
  !
  !
  !
  ! <<< HEADER ANALYSIS >>>
  !
  ! - find out, which type of data is in which column of the table
  !
  do ! loop header parsing
    !
    ! further analyse the content of sline
    ! determine the number of data items per table row = ncol
    !
    if (len_trim(sline)<1) goto 20 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*, iostat=ioerr, err=802) snam ! get the name of the item
    if (snam(1:1)=="_") then ! this is a data name
      ! count
      icol = icol+1
      if (icol>ncol) exit ! max. number of columns reached, just exit the counting loop
      if (      trim(adjustl(snam))=="_symmetry_equiv_pos_as_xyz" &
         & .or. trim(adjustl(snam))=="_space_group_symop_operation_xyz" ) then
        idxop = icol ! store the index of the operation in the table rows
      end if
      if (      trim(adjustl(snam))=="_symmetry_equiv_pos_site_id" &
         & .or. trim(adjustl(snam))=="_space_group_symop_id" ) then
        idxid = icol
      end if
    else ! this is no data name 
      ! exit the counting loop
      exit
    end if
    !
    ! we are here, because a name was found, remove it now from sline
    call CIF_str_remsub(sline,adjustl(trim(snam)))
    sline = trim(adjustl(sline))
    cycle ! try next content of sline
    !
20  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! header parsing
  !
  if (icol==0) then
    write(unit=smsg,fmt=*) iline
    call CIF_ERR("Parsing table header failed in line #"// &
       & trim(adjustl(smsg))//".")
    goto 1000 ! corrupt table, do nothing more
  end if
  !
  if (idxop==0) then
    write(unit=smsg,fmt=*) lstats(1)
    call CIF_WRN("Missing item _*_xyz in table of line #"// &
       & trim(adjustl(smsg))//".")
    call CIF_WRN("Failed to load symmetry operation data.")
    goto 1000 ! invalid table, do nothing more
  end if
  !
  if (idxid==0) then ! this is a non-lethal condition
    write(unit=smsg,fmt=*) lstats(1)
    call CIF_WRN("Missing item _*_id in table of line #"// &
       & trim(adjustl(smsg))//".")
  end if
  !
  !
  ! <<< Data Reading >>>
  !
25 continue
  icol = 1 ! current column
  irow = 1 ! current row
  do ! data reading loop
    !
    ! Read items, check if they can be data. If not, stop!
    ! Fill up the table.
    !
    if (len_trim(sline)<1) goto 30 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*, iostat=ioerr, err=802) sval ! get the item
    if (sval(1:1)=="_") then ! this is a data name
      exit ! stop reading data
    elseif (sval(1:5)=="loop_") then ! this is another loop
      exit ! stop reading data
    elseif (sval(1:5)=="data_" .or. sval(1:5)=="save_") then ! this is other data
      exit ! stop reading data
    else ! this can be data
      !
      ! interprete the data depending on the column type
      !
      if (icol==idxid) then ! read the ID
        read (unit=sval, fmt=*, iostat=ioerr, err=800) symop_id
      end if
      !
      if (icol==idxop) then ! read the operation
        symop_xyz = trim(adjustl(sval)) ! does this really work?
      end if
      !
    end if
    !
    ! remove value string from sline
    call CIF_str_remsub(sline,adjustl(trim(sval)))
    sline = trim(adjustl(sline))
    !
    icol = icol + 1 ! advance column index
    if (icol>ncol) then ! last column reached
      !
      ! apply the symmetry operation to the current atom sites
      !
      CALL SYMOPS_SET_STR(TRIM(symop_xyz),irow,ichk)
      if (ichk/=1) goto 800
      !
      ! prepare for next row
      icol = 1 ! move back to columns one
      irow = irow + 1 ! advance row index
      symop_id = 0
      symop_xyz = ""
    end if
    if (irow>nrow) exit ! stop reading data, table is complete // alternative: goto 35
    cycle ! next item
    !
30  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do ! data counting loop
  !
35 continue
  !
  goto 1000
  !
  !
  ! <<< Error Handling >>>
  !
800 continue
  nchk = 1
  call CIF_ERR("Failed to interprete the value from '"//trim(sval) &
     & //"', skipping table.", ichk)
  goto 1000
801 continue
  nchk = 2
  call CIF_ERR("Failed to allocate memory.", nalloc)
  goto 1000
802 continue
  nchk = 3
  call CIF_ERR("Failed to read from string '"//trim(sline)// &
     & "'.", ioerr)
  goto 1000
  !
  ! <<< Exit >>>
  !
1000 continue
  !
END SUBROUTINE CIF_read_loop3



!********************************************************
! CIF_READ
! This subroutine reads data from a file and interpretes
! it as structure information.
!********************************************************
SUBROUTINE CIF_READ(sfile, nchk)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nchk ! error code (1=success, <=0 error)
  integer*4 :: ioerr ! I/O error code
  integer*4 :: ichk ! internal check number
  integer*4 :: ilun ! logical file unit
  integer*4 :: iline ! current file line
  integer*4 :: nlines ! number of lines in the CIF file
  integer*4 :: iloop ! identifies the current loop
  integer*4 :: nloops ! number of loops in the CIF file
  integer*4 :: nalloc ! allocation state
  integer*4 :: icomm ! comment location
  integer*4 :: i, j, k, l, n
  integer*4 :: isymm ! structure symmetry operations found
  integer*4 :: nsymm ! number of symmetry operations to apply
  real(dp) :: H(3,3) ! cell transformation matrix
  character(len=lsl) :: sline, smsg ! current line of the file
  character(len=ssl) :: sname ! item name string
  integer*4, external :: getfreelun
  !
  !
  ! <<< Initialization >>>
  !
  call CIF_MSG("begin CIF_READ", 1)
  i = 0
  j = 0
  nchk = 0 ! error code init (0=unknown error)
  ichk = 0 !
  ilun = getfreelun() ! get a unused I/O unit
  ioerr = 0 !
  iline = 0 !
  nlines = 0 !
  iloop = 0 !
  nloops = 0 !
  nalloc = 0 !
  isymm = 0 !
  nsymm = 0 !
  H = 0.0d+0
  call CIF_MSG("- clearing current data", 2)
  call CIF_CLEAR()
  if (ilun<0) goto 101
  !
  !
  ! <<< Open the file >>>
  !
  open (unit=ilun, file=trim(sfile), status='UNKNOWN', iostat=ioerr)
  if (ioerr/=0) goto 102
  rewind(unit=ilun, iostat=ioerr)
  if (ioerr/=0) goto 103
  !
  !
  ! <<< Load the file content as sequence of strings to memory >>>
  !
  iline = 0
  do ! determine the number of lines in the file
     ! try to read the next line
     read (unit=ilun, fmt='(a1024)', iostat=ioerr) sline
     if (-1==ioerr) exit ! could not read due to end of file
     if (0/=ioerr) exit ! could not read due to other error
     iline = iline+1
  end do
  nlines = iline
  write(unit=smsg,fmt=*) nlines
  call CIF_MSG("- reading "//trim(adjustl(smsg))//" lines.", 2)
  if (nlines<1) goto 105
  !
  ! allocate memory for storing the CIF text
  if (allocated(CIF_lines)) deallocate(CIF_lines, stat=nalloc)
  allocate(CIF_lines(nlines), stat=nalloc)
  if (nalloc/=0) goto 106
  CIF_lines = ""
  !
  ! read the lines again, this time store them all in the module array
  rewind(unit=ilun, iostat=ioerr)
  if (ioerr/=0) goto 103
  do iline=1, nlines ! read only the same number of lines as successfully read before
     ! try to read the next line
     read (unit=ilun, fmt='(a1024)', iostat=ioerr) sline
     if (-1==ioerr) exit ! could not read due to end of file
     if (0/=ioerr) exit ! could not read due to other error
     ! remove comments from the line
     icomm = index(sline,"#")
     if (icomm==1) sline = "" ! full comment line, erase
     if (icomm>1) sline = sline(1:icomm-1) ! cut comments
     CIF_lines(iline) = trim(adjustl(sline)) ! store remaining content
  end do
  ! close the file
  close (unit=ilun, iostat=ioerr)
  call CIF_MSG("- finished reading text lines.", 2)
  !
  ! The content of the CIF file is now stored in CIF_lines.
  ! We will extract the relevant data below.
  !
  !
  ! <<< Analyse the file structure >>>
  !
  ! Count the number of "_loop"
  ! - we assume here, that "_loop" doesn't appear twice in a line
  do iline=1, nlines
    i = index(CIF_lines(iline),"loop_")
    if (i>0) nloops = nloops+1
  end do
  write(unit=smsg,fmt=*) nloops
  call CIF_MSG("- found "//trim(adjustl(smsg))//" tables ...", 2)
  iloop = 0 ! reset loop selector to - "no loop present"
  if (nloops>0) then ! there are loops in the CIF ...
    ! get their stats
    if (allocated(CIF_lstats)) deallocate(CIF_lstats,stat=nalloc)
    allocate(CIF_lstats(CIF_lstat_len,nloops),stat=nalloc)
    CIF_lstats = 0
    !
    ! - Retrieve the positions, number of columns, and number of rows of each loop
    j = 0
    iline = 1
    l = 1
    sline = ""
    do
      sline = trim(CIF_lines(l))
      if (len_trim(sline)>0) then
        i = index(sline,"loop_")
        if (i>0) then
          k = len_trim(sline)
          sline = adjustl(sline(i:k))
          j = j + 1
          ! analyse the loop structure
          call CIF_get_loop_stats(l,sline,CIF_lstats(:,j))
        end if
      end if
      l = l + 1
      if (l>nlines) exit
    end do
  end if
  !
  !
  ! <<< Read the basic data items (cell definition) >>>
  !
  call CIF_MSG("- reading cell data.", 2)
  iline = 1
  if (nloops>0) iloop = 1 ! select the first loop in the file for jumping
  if (iloop>0) then ! check if it starts with a loop
    if (iline==CIF_lstats(1,iloop)) then ! it does ...
      iline = CIF_lstats(2,iloop)+1 ! jump over it
      iloop = iloop+1
      if (iloop>nloops) iloop = 0
    end if
  end if
  sline = trim(adjustl(CIF_lines(iline)))
  do ! main reading/parsing loop
    !
    ! - reads primary structure data
    ! - skips loops
    !
    if (len_trim(sline)<1) goto 20 ! skip empty lines
    !
    ! content parser
    read (unit=sline, fmt=*) sname ! get the name of the item
    sname = trim(adjustl(sname))
    if (sname(1:1)=="_") then ! singular data name
      call CIF_read_item(sline) ! read the item, clear it from the string
      cycle ! proceed parsing the rest of this line
    else ! none of the above applied, remove sname from sline
      i = len_trim(sname)
      j = len_trim(sline)
      if (i>=j) then
        sline = ""
      else
        sline = trim(adjustl(sline(i+1:j)))
      end if
      cycle ! proceed parsing the rest of this line
    end if
    !
20  iline = iline+1 ! proceed to next line
    if (iline>nlines) exit ! done, final line reached
    if (iloop>0) then ! check for loop
LSKL: do ! this code loop jumps over consecutive CIF loops
        if (iline==CIF_lstats(1,iloop)) then ! a loop starts here
          iline = CIF_lstats(2,iloop)+1 ! jump over the loop
          iloop = iloop+1 ! mark next loop
          if (iloop>nloops) then
            iloop = 0 ! no more loops coming
            exit LSKL
          end if
        else
          exit LSKL
        end if
      end do LSKL
    end if
    if (iline>nlines) exit ! done, final line reached (check twice here because of loop jumping)
    sline = trim(adjustl(CIF_lines(iline))) ! copy content of next line to sline
    !
  end do
  !
  call CONVMAT( CIF_cell_length_a, CIF_cell_length_b, &
              & CIF_cell_length_c, CIF_cell_angle_alpha, &
              & CIF_cell_angle_beta, CIF_cell_angle_gamma, H)
  !
  ! <<< Read the loop data (tables) >>>
  !
  if (nloops>0) then
    !
    ! - process all primary tables (they contain atomic site coordinates)
    call CIF_MSG("- processing primary atom site tables ...", 2)
    !
    ! - first determine the number of data rows
    n = 0
    do iloop=1, nloops
      if (1==CIF_lstats(5,iloop)) then ! this is a primary loop
        n = n + CIF_lstats(4,iloop) ! sum up number of table rows
      end if
    end do
    if (n<1) then
      call CIF_WRN("The CIF file contains no atom site definitions.")
      goto 100 ! stop reading, and exit
    end if
    !
    ! - prepare local memory for receiving data
    if (allocated(CIF_atom_site)) deallocate(CIF_atom_site,stat=nalloc)
    allocate(CIF_atom_site(CIF_atom_site_items_num,n), stat=nalloc) ! allocate
    CIF_atom_site = 0.0d+0 ! initialize
    CIF_atom_site_number = 0 ! CIF_atom_site_number keeps track of the rows used
    !
    ! - now read process the table data
    do iloop=1, nloops
      if (1==CIF_lstats(5,iloop)) then ! this is a primary loop
        call CIF_read_loop1(CIF_lstats(:,iloop), ichk)
        if (ichk/=0) then
          write(unit=smsg,fmt=*) CIF_lstats(1,iloop)
          call CIF_ERR("Failure while interpreting primary table (line #"// &
             & trim(adjustl(smsg))//").")
        end if
      end if
    end do
    write(unit=smsg,fmt=*) CIF_atom_site_number
    call CIF_MSG("  - found "//trim(adjustl(smsg))// &
       & " primary atom site definitions.", 3)
    !
    !
    ! - process all secondary tables
    !   they contain additional atomic site data
    !   that we needs to fill into the current CIF_atom_site table
    !
    ! - read process the table data
    do iloop=1, nloops
      if (2==CIF_lstats(5,iloop)) then ! this is a secondary loop
        write(unit=smsg,fmt=*) CIF_lstats(1,iloop)
        call CIF_MSG("- processing secondary atom site table of line #"// &
           & trim(adjustl(smsg))//" ...", 2)
        call CIF_read_loop2(CIF_lstats(:,iloop),ichk)
        if (ichk/=0) then
          call CIF_ERR("Failure while interpreting secondary table (line #"// &
             & trim(adjustl(smsg))//").")
        end if
      end if
    end do
    !
    !
    ! - process the symmetry tables
    !   always apply the explicit symmetry table if given in the file.
    !   Check later if symmetry was applied and apply other symmetry
    !   definitions later only if no explicit table was found.
    !
    ! - read process the table data
    do iloop=1, nloops
      if (3==CIF_lstats(5,iloop)) then ! this is a symmetry loop
        write(unit=smsg,fmt=*) CIF_lstats(1,iloop)
        call CIF_MSG("- processing symmetry table of line #"// &
           & trim(adjustl(smsg))//" ...", 2)
        call CIF_read_loop3(CIF_lstats(:,iloop),ichk)
        if (ichk/=0) then
          call CIF_ERR("Failure while interpreting symmetry table (line #"// &
             & trim(adjustl(smsg))//").")
          cycle ! try next symmetry table
        end if
        isymm = isymm + 1 ! indicate that symmetry operations were loaded from file
      end if
    end do
    !
    ! - process the read symmetry data
    !
    ! - initialize the symops module
    nsymm = 0
    if (isymm>0) then ! there was a symmetry table
      !
      ! get the number of symmetry operations
      if (allocated(symops_trf)) nsymm = SIZE(symops_trf,2)
      !
    else if (CIF_symmetry_Int_Tables_number>1) then ! ... a non-trivial space group number
      !
      ! set the symmetry operations from the space group number
      call SYMOPS_SET_SGNUM(CIF_symmetry_Int_Tables_number, ichk)
      if (ichk/=1) goto 107
      ! get the number of symmetry operations
      if (allocated(symops_trf)) nsymm = SIZE(symops_trf,2)
      !
    else if (len_trim(CIF_symmetry_space_group_name_HM)>0) then ! ... a space group name given
      !
      ! set the symmetry operations from the space group name
      call SYMOPS_SET_SGNAME(trim(CIF_symmetry_space_group_name_HM), ichk)
      if (ichk/=1) goto 108
      ! get the number of symmetry operations
      if (allocated(symops_trf)) nsymm = SIZE(symops_trf,2)
      !
    end if
    !
    if (nsymm>0) then ! there are symmetry operations
      !
      ! apply the symmetry operations
      call CIF_MSG("  - applying symmetry operations ...", 3)
      call SYMOPS_APPLY( H, CIF_atom_site, 0.5d+0, ichk)
      if (ichk/=1) goto 108
      CIF_atom_site_number = SIZE(CIF_atom_site, 2)
      write(unit=smsg,fmt=*) CIF_atom_site_number
      call CIF_MSG("  - structure contains "//trim(adjustl(smsg))// &
         & " atom sites after applying symmetry operations.", 3)
      !
    end if 
    !
    !
  end if
  !
  
  !
  !
  !
  ! <<< Clean up and finish. >>> 
  !
100 call CIF_MSG("- finished reading CIF", 2)
  call CIF_CLEAR_LINES()
  nchk = 1
  goto 1000
  !
  !
  ! error handling
101 continue !
  nchk = -1
  call CIF_ERR("No free logical I/O unit.")
  goto 1000
102 continue !
  nchk = -2
  call CIF_ERR("Failed to open file ["//trim(sfile)// &
     & "] for reading.", ioerr)
  goto 1000
103 continue !
  nchk = -3
  call CIF_ERR("File operation failed.", ioerr)
  goto 1000
104 continue !
  nchk = -4
  write(unit=smsg,fmt=*) iline
  call CIF_ERR("Failed to read line #"//trim(adjustl(smsg))//".", ioerr)
  goto 1000
105 continue
  nchk = -5
  call CIF_ERR("CIF file seems to be empty.")
  goto 1000
106 continue
  nchk = -6
  call CIF_ERR("Failed to allocate heap memory.")
  goto 1000
107 continue
  nchk = -7
  call CIF_ERR("Invalid space group number.")
  goto 1000
108 continue
  nchk = -8
  call CIF_ERR("Failed to apply symmetry operations.")
  goto 1000
  
1000 continue
  !
END SUBROUTINE CIF_READ


!********************************************************
! CIF_create_atom_type_list
! This subroutine analyses the current CIF_atom_site 
! table and generates an atom type list from it.
! Atom types are distinguished by atomic number and
! oxidation state. This routine generates in addition
! a list of type assignments, where each atom site is
! assigned to one atom type.
!********************************************************
SUBROUTINE CIF_create_atom_type_list(T, A, ntypes, nchk)
  implicit none
  real(dp), dimension(:,:), allocatable :: T
  real(dp), dimension(:), allocatable :: A
  integer, intent(out) :: ntypes, nchk
  integer :: i, j, k, l
  integer :: nalloc
  real(dp), dimension(:,:), allocatable :: TTMP
  nalloc = 0
  nchk = 0
  ntypes = 0
  if (allocated(T)) deallocate(T,stat=nalloc)
  if (allocated(A)) deallocate(A,stat=nalloc)
  if (CIF_atom_site_number<=0) goto 100 ! nothing to do
  allocate(A(CIF_atom_site_number),stat=nalloc)
  if (nalloc/=0) goto 106
  A = 0.d+0
  allocate(TTMP(3,CIF_atom_site_number),stat=nalloc)
  if (nalloc/=0) goto 106
  TTMP = 0.d+0
  ! determine atom types
  ! - loop through site list and gather information
  ! - create types on the fly
  ! - count types of same atomic number to create label indices
  do i=1, CIF_atom_site_number
    k = 0
    l = 0
    if (ntypes>0) then ! compare to existing types
      do j=1, ntypes
        if ( nint( CIF_atom_site(4,i) ) == nint( TTMP(1,j) ) ) then ! same atomic number
          l = l + 1 ! count them
          if ( 0.01 > dabs( CIF_atom_site(5,i)-TTMP(2,j) ) ) then ! same ionic charge // oxidation state
            k = j ! assign
          end if
        end if
      end do
    end if
    if (k>0) then ! assign existing type
      A(i) = dble(k)
    else ! create new type and assign
      ntypes = ntypes + 1
      TTMP(1,ntypes) = CIF_atom_site(4,i)
      TTMP(2,ntypes) = CIF_atom_site(5,i)
      TTMP(3,ntypes) = dble(l+1)
      A(i) = dble(ntypes)
    end if
  end do
  if (ntypes>0) then ! create the list returned
    allocate(T(3,ntypes),stat=nalloc)
    if (nalloc/=0) goto 106
    T(1:3,1:ntypes) = TTMP(1:3,1:ntypes)
  end if
  deallocate(TTMP,stat=nalloc)
100 continue ! fine exit
  nchk = 1
  return
101 continue ! error exits
  nchk = 0
  goto 1000
106 continue
  nchk = -6
  call CIF_ERR("Failed to allocate heap memory.")
  goto 1000
1000 continue
END SUBROUTINE CIF_create_atom_type_list


!********************************************************
! CIF_WRITE
! This subroutine writes structure data to a file in
! standard CIF format.
!********************************************************
SUBROUTINE CIF_WRITE(sfile, nchk)
  implicit none
  character(len=*), intent(in) :: sfile ! file name
  integer*4, intent(inout) :: nchk ! error code (1==success, <=0 error)
  integer*4 :: ioerr ! I/O error code
  integer*4 :: ichk ! internal check number
  integer*4 :: ilun ! logical file unit
  integer*4 :: iline ! current file line
  integer*4 :: nlines ! number of lines in the CIF file
  integer*4 :: iloop ! identifies the current loop
  integer*4 :: nloops ! number of loops in the CIF file
  integer*4 :: nalloc ! allocation state
  integer*4 :: icomm ! comment location
  integer*4 :: i, j, k, l, n
  integer*4 :: isymm ! structure symmetry operations found
  integer*4 :: nsymm ! number of symmetry operations to apply
  real(dp), allocatable, dimension(:,:) :: atom_types ! atom type list
  real(dp), allocatable, dimension(:) :: atom_site_types ! atom site type assignment list
  real(dp) :: dtmp ! temp. value
  character(len=lsl) :: sline ! current line of the file
  character(len=lsl) :: stmp, snum, smsg ! temp. strings
  character(len=ssl) :: sname ! item name string
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  !
  ! <<< Initialization >>>
  !
  call CIF_MSG("begin CIF_WRITE", 1)
1001 FORMAT(F18.4)
1002 FORMAT(F18.6)
  i = 0
  j = 0
  nchk = 0 ! success code init (0)
  ichk = 0 !
  ilun = getfreelun() ! get a unused I/O unit
  ioerr = 0 !
  iline = 0 !
  nlines = 0 !
  iloop = 0 !
  nloops = 0 !
  nalloc = 0 !
  isymm = 0 !
  nsymm = 0 !
  if (ilun<0) goto 101
  !
  !
  ! <<< Open the file for writing >>>
  !
  call createfilefolder(trim(sfile),ioerr)
  if (ioerr/=0) goto 108
  open (unit=ilun, file=trim(sfile), status='REPLACE', action='WRITE', iostat=ioerr)
  if (ioerr/=0) goto 102
  !
  !
  ! <<< Write the CIF file header >>>
  !
  sline = "#\#CIF_2.0"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  sline = "data_CELSLC"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  sline = ""
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  !
  ! <<< Write the Crystal data >>>
  !
  sline = "# Crystal data"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  sline = "_pd_phase_name                         'Super Cell'"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt=1002) CIF_cell_length_a
  sline = "_cell_length_a                         "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt=1002) CIF_cell_length_b
  sline = "_cell_length_b                         "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt=1002) CIF_cell_length_c
  sline = "_cell_length_c                         "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt=1001) CIF_cell_angle_alpha
  sline = "_cell_angle_alpha                      "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt=1001) CIF_cell_angle_beta
  sline = "_cell_angle_beta                       "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt=1001) CIF_cell_angle_gamma
  sline = "_cell_angle_gamma                      "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  if (CIF_symmetry_Int_Tables_number<1 .or. CIF_symmetry_Int_Tables_number>sg_nummax) goto 107
  call SG_INIT(ichk) ! <- MODULE spacegroups initialized here
  i = CIF_symmetry_Int_Tables_number
  call SG_NUMGETNAME(i,stmp,ichk)
  sline = "_symmetry_space_group_name_H-M_alt     '"//trim(adjustl(stmp))//"'"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  write (unit=stmp, fmt='(I)') CIF_symmetry_Int_Tables_number
  sline = "_symmetry_Int_Tables_number            "//trim(adjustl(stmp))
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  sline = ""
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  !
  ! <<< Write Symmetry operation table >>>
  !
  !
  sline = "# Symmetry operations"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  sline = "loop_"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  sline = "    _space_group_symop_id"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  sline = "    _space_group_symop_operation_xyz"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  call SG_NUMGETSYMNUM(CIF_symmetry_Int_Tables_number,n,ichk)
  do j=1, n
    write (unit=stmp, fmt='(I)') j
    call SG_NUMGETSYMOP(CIF_symmetry_Int_Tables_number,j,snum(1:sg_soplen),ichk)
    sline = "    "//trim(adjustl(stmp))//"    '"//trim(snum(1:sg_soplen))//"'"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  end do
  call SG_UNINIT() ! <- MODULE spacegroups uninitialized here
  !
  sline = ""
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  !
  if (CIF_atom_site_number>0) then
    !
    ! <<< Write the atom type table >>>
    !
    sline = "# Atomic types"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "loop_"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_type_symbol"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_type_oxidation_number"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    !
    ! create a list of atom types and assign the type indices to the atomic site list
    call CIF_create_atom_type_list(atom_types, atom_site_types, k, ichk)
    if (ichk<=0) goto 105
    if (k>0 .and. allocated(atom_types) .and. allocated(atom_site_types) ) then ! write the atom type table
      do i=1, k
        j = nint(atom_types(1,i))
        call CIF_get_atom_type_symbol(j,atom_types(2,i),stmp)
        write(unit=snum,fmt=1001) atom_types(2,i)
        sline = "    "//trim(adjustl(stmp))//"  "//trim(adjustl(snum))
        write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
      end do
    end if
    !
    sline = ""
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    !
    !
    ! <<< Write the atomic site table >>>
    !
    sline = "# Atomic sites"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "loop_"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_label"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_type_symbol"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_occupancy"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_fract_x"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_fract_y"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_fract_z"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_thermal_displace_type"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    sline = "    _atom_site_U_iso_or_equiv"
    write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    !
    do i=1, CIF_atom_site_number ! write a line for each atom site
      ! init
      sline = ""
      j = CIF_atom_site(4,i) ! atomic number
      l = 0
      ! generate label
      call CIF_get_atom_type_symbol(j,0.0D+0,stmp) ! only the symbol here
      snum = ""
      if (k>0 .and. allocated(atom_types) .and. allocated(atom_site_types) ) then ! add label index
        l = nint(atom_site_types(i))
        write(unit=snum,fmt='(I)') nint(atom_types(3,l))
      end if
      sline = "    "//trim(stmp)//trim(adjustl(snum))
      ! add type symbol
      call CIF_get_atom_type_symbol(j,CIF_atom_site(5,i),stmp) ! atom type with charger / oxidation state
      sline = trim(sline)//"  "//trim(adjustl(stmp))
      ! add the occupancy
      dtmp = min(1.D+0,max(0.D+0,CIF_atom_site(6,i)))
      write(unit=stmp,fmt=1001) dtmp
      sline = trim(sline)//"  "//trim(adjustl(stmp))
      ! add fractional coordinates
      write(unit=stmp,fmt=1002) CIF_atom_site(1,i)
      sline = trim(sline)//"  "//trim(adjustl(stmp))
      write(unit=stmp,fmt=1002) CIF_atom_site(2,i)
      sline = trim(sline)//"  "//trim(adjustl(stmp))
      write(unit=stmp,fmt=1002) CIF_atom_site(3,i)
      sline = trim(sline)//"  "//trim(adjustl(stmp))
      ! add Uiso (type and value)
      call CIF_get_atom_site_uiso(i,dtmp)
      write(unit=stmp,fmt=1002) dtmp
      sline = trim(sline)//"  Uiso  "//trim(adjustl(stmp))
      ! write the line
      write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
    end do
    !
  end if
  !
  sline = ""
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  !
  ! <<< Last line >>>
  !
  sline = "# End of data_CELSLC"
  write (unit=ilun, fmt='(A)', iostat=ioerr, err=104) trim(sline)
  !
  ! <<< Close the file >>>
  !
  close (unit=ilun, iostat=ioerr)
  call CIF_MSG("- finished writing CIF.", 2)
  !
  !
  !
  ! <<< Clean up and finish. >>> 
  !
  nchk = 1
100 goto 1000
  !
  !
  ! error handling
101 continue !
  nchk = -1
  call CIF_ERR("No free logical I/O unit.")
  goto 1000
102 continue !
  nchk = -2
  call CIF_ERR("Failed to open file ["//trim(sfile)// &
     & "] for writing.", ioerr)
  goto 1000
103 continue !
  nchk = -3
  call CIF_ERR("File operation failed.", ioerr)
  close (unit=ilun, iostat=ioerr)
  goto 1000
104 continue !
  nchk = -4
  call CIF_ERR("Failed to write to file.", ioerr)
  close (unit=ilun, iostat=ioerr)
  call SG_ISREADY(j)
  if (j/=0) call SG_UNINIT()
  goto 1000
105 continue !
  nchk = -5
  call CIF_ERR("Failed to generate atom type list.", ichk)
  close (unit=ilun, iostat=ioerr)
  goto 1000
106 continue
  nchk = -6
  call CIF_ERR("Failed to allocate heap memory.")
  goto 1000
107 continue
  nchk = -7
  call CIF_ERR("Invalid space group number.")
  close (unit=ilun, iostat=ioerr)
  goto 1000
108 continue
  nchk = -8
  call CIF_ERR("Failed to create folder to store file ["//trim(sfile)// &
     & "].", ioerr)
  goto 1000
  
1000 continue
  if (allocated(atom_types)) deallocate(atom_types,stat=nalloc)
  if (allocated(atom_site_types)) deallocate(atom_site_types,stat=nalloc)
  !
END SUBROUTINE CIF_WRITE



END MODULE cifio
!**********************************************************************!
!**********************************************************************!