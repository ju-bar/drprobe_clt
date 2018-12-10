!**********************************************************************!
!**********************************************************************!
!
! FILE: "fscatab.f90"
!
! AUTHOR: Dr. J. Barthel
!         Forschungszentrum Jülich
!         Jülich, Germany
!
! PURPOSE: Implementations for epw simulations
!          Input and output and evaluation of scattering factor tables
!
! LINK:   ..\basics\binio2.f90
!         ..\basics\levmrq.F90
!
! VERSION: 0.01b, J.B., 13.09.2013
!          - implementation of initial module structure
!            MODULE FSCATAB, shortcut FST
!          1.0, J.B., 13.11.2014
!          - tested all routines, removed some bugs.
!          1.1, J.B., 29.05.2018
!          - added f^2(theta) integrals
!
!**********************************************************************!
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



!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
!
MODULE FSCATAB

  IMPLICIT NONE

  SAVE
  

!**********************************************************************!
! GLOBAL DECLARATIONS
!
! --------------------------------------------------------------------
! parameters
  integer*4, public, parameter :: FST_ll = 1024
  integer*4, public, parameter :: FST_ntab_max = 1024 ! max. number of tables (elements+ions)
  integer*4, public, parameter :: FST_ntab_len = 256 ! max. number of table items (tabulates scattering angles)
  integer*4, public, parameter :: FST_prm_num = 12 ! number of parameters for fitted function
  
  real*4, public, parameter :: FST_twopi = 6.28318530718
  real*4, public, parameter :: FST_fourpi = 12.56637061436
  real*4, public, parameter :: FST_r4pi = 0.079577471546
  real*4, public, parameter :: FST_r8pi2 = 0.0126651479553
  real*4, public, parameter :: FST_elmc2 = 510.998946269 ! m0*c^2 [keV]

! --------------------------------------------------------------------
! variables
  character(len=FST_ll), public :: FST_smsg ! string for messaging
  character(len=FST_ll), public :: FST_sinfo ! string for messaging
  character(len=FST_ll), public :: FST_serr ! string for error messaging
  character(len=FST_ll), public :: FST_swarn ! string for waring messaging
  
  integer*4, public :: FST_dowarn ! flag for warning output
  DATA FST_dowarn /1/
  integer*4, public :: FST_domsg ! flag for message output
  DATA FST_domsg /1/
  integer*4, public :: FST_lasterr ! last error code
  DATA FST_lasterr /0/

! (STACK ALLOCATIONS)
! - listed elements symbols (Symb)_i
  character(len=20) :: FST_sym(FST_ntab_max)
! - listed elements elemental charge and oxidation state (used flag, Z, dZ, ntab)_i
  real*4, public :: FST_crg(4,FST_ntab_max)
! - listed elements electron scattering factor tables (s,fel)_i,j
  real*4, public :: FST_fel(2,FST_ntab_len,FST_ntab_max)
! - listed elements electron scattering factor parameters (a)_l
  real*4, public :: FST_fep(FST_prm_num,FST_ntab_max)
! - listed elements electron scattering factor parameterization errors (sigma,E,R)
  real*4, public :: FST_per(3,FST_ntab_max)

! - listed number of tables
  integer*4, public :: FST_nat
  DATA FST_nat /0/
  
! - current atom index in the current FFE table
! ! Warning: The value maybe volatile, carefully check how it is used!
! ! Warning: This value maybe corrupted in case of multi-threading.
! !          Need to extend it to an array if parallel computing
! !          is applied in one of the related functions.
  integer*4, public :: FST_PRM_J
  DATA FST_PRM_J /0/
! - current atom Debye-Waller parameter B [A^2] in the current FFE table
! ! Warning: same as for FST_PRM_J
  real*4, public :: FST_PRM_J_BISO
  DATA FST_PRM_J_BISO /0.0/
! - current diffraction vector G in the absorptive form factor integration
! ! Warning: same as for FST_PRM_J
  real*4, public :: FST_PRM_G
  DATA FST_PRM_G /0.0/
! - current incident wave vector K0 in the absorptive form factor integration
! ! Warning: same as for FST_PRM_J
  real*4, public :: FST_PRM_K0
  DATA FST_PRM_K0 /0.0/

! - temporaray fit parameters  
  real*4, public :: FST_feptmp(FST_prm_num)

! --------------------------------------------------------------------
! routines
  private :: FST_ERROR, FST_Warning, FST_Message
  private :: FST_GETFREELUN
  public :: FST_INIT
  public :: FST_LoadFX
  public :: FST_LoadFEP
  public :: FST_SaveFEP
  public :: FST_FX2FE
  public :: FST_FITFE
  public :: FST_GETIDX
  public :: FST_GETFE
  public :: FST_GETMU0
  public :: FST_GETMUG
  public :: FST_GETF2G
  public :: FST_GETSCAR1
  public :: FST_GETSCAC
  public :: FST_GETSCAR
  public :: FST_GETSCAABS
  public :: FST_GETSCAF2
  public :: FST_GETSCADWF
  public :: FST_GETCRG
  
  ! fit functions
  public :: FST_FFG3
  public :: FST_FFL3G3
  public :: FST_FFLL3G3
  public :: FST_FFLLGN
  public :: FST_FFQLGN
  
  
  
!
!**********************************************************************!
! IMPLEMENTATIONS
  CONTAINS




!**********************************************************************!
!
! subroutine FST_ERROR
!
! posts an error message to console
!
subroutine FST_ERROR(smessage, ncode)

  implicit none
  
  character(len=*), intent(in) :: smessage
  integer*4, intent(in) :: ncode
  
  FST_lasterr = ncode
  FST_serr = trim(smessage)
  write (unit=*,fmt='(A,I5.5,") ",A)') "ERROR: (FST-", ncode, &
     &   trim(FST_serr)

  return

end subroutine FST_ERROR
!**********************************************************************!



!**********************************************************************!
!
! subroutine FST_Warning
!
! posts a warning message to console.
!
subroutine FST_Warning(smessage)

  implicit none
  
  character(len=*), intent(in) :: smessage
  
  FST_swarn = trim(smessage)
  write (unit=*,fmt='(A)') "Warning: (FST) "//trim(FST_swarn)

  return

end subroutine FST_Warning


!**********************************************************************!
!
! subroutine FST_Message
!
! posts a message to console.
!
subroutine FST_Message(smessage)

  implicit none
  
  character(len=*), intent(in) :: smessage
  
  if (FST_domsg<=0) return

  FST_smsg = trim(smessage)
  write (unit=*,fmt='(A)') trim(FST_smsg)

  return

end subroutine FST_Message
!**********************************************************************!


!**********************************************************************!
!
! subroutine FST_INIT
!
! initializes the scattering tables, removes all content
!
subroutine FST_INIT

  implicit none
  
  FST_lasterr = 0
  
  FST_smsg = ""
  FST_serr = ""
  FST_swarn = ""
  
  FST_nat = 0
  FST_sym = ""
  FST_crg = 0.0
  FST_fel = 0.0
  FST_fep = 0.0
  FST_per = 0.0
  
  FST_sinfo = "FST: initialized."
  call FST_Message(trim(FST_sinfo))
  
  return
  
end subroutine FST_INIT


!**********************************************************************!
!
! FST_GETFREELUN
!
! function returning next free logical unit number
!
function FST_GETFREELUN()

  implicit none
  
  integer*4 :: FST_GETFREELUN
  
  logical :: isopen             ! isopen request from inquire
  integer*4 :: lun, nret, nerr  ! temporary luns
  
  nret = 0
  
  do lun = 20, 99
    inquire(unit=lun,opened=isopen,iostat=nerr)
    if (nerr/=0) then
      FST_GETFREELUN = 0
      return
    end if
    if (.not.isopen) then
      nret = lun
      exit
    end if
  end do
  
  return

end function FST_GETFREELUN




!**********************************************************************!
!
! subroutine FST_LoadFX
!
! loads new x-ray scattering factors from file
! translates to electron scattering factors
!
subroutine FST_LoadFX(sfile, nerr)

  implicit none
  
  integer*4, parameter :: subid = 200 ! error id
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  
  logical :: fexists
  integer*4 :: lun, ioerr, nline, nitem, i
  real*4 :: fz, fdz, fs, fx
  real*4 :: tabFX(2,FST_ntab_len) ! temporary table of x-ray scattering factors
  character(len=FST_ll) :: sitem1, sitem2, snum, ssym, smsg
  
  ! initialize
  nerr = 0
  nline = 0
  nitem = 0
  
  ! clear old tables
  FST_nat = 0
  FST_sym = ""
  FST_crg = 0.0
  FST_fel = 0.0
  FST_fep = 0.0
  FST_per = 0.0
  
  ! check if the file exists
  inquire(file=trim(sfile),exist=fexists)
  if (.not.fexists) goto 201
  ! check for available file unit
  lun = FST_GETFREELUN()
  if (lun<0) goto 202
  
  FST_sinfo = "FST: loading x-ray scattering factors from file [" &
     & //trim(sfile)//"]."
  call FST_Message(trim(FST_sinfo))
  
  ! open the file for reading
  open( unit=lun, file=trim(sfile), &
      & iostat=ioerr, action='READ', status='OLD' )
  if (ioerr/=0) goto 203
  
  ! READING TABLE HEADER LINE (not used)
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt=*) sitem1
  if (ioerr/=0) goto 204
  
  ! READING TABLES
  !
  ! EXPECTED FORMAT:
  ! TABLE
  ! <Atomic symbol> (Cu)
  ! <Atomic core charge, atomic Z number> (29)
  ! <Atomic charge> (0)
  ! <scattering angle s[A-1], scattering factor fX> (0.00, 29.0000)
  ! <scattering angle s[A-1], scattering factor fX> (0.05, 28.4481)
  ! <scattering angle s[A-1], scattering factor fX> (0.05, 27.0853)
  ! ...
  ! <scattering angle s[A-1], scattering factor fX> (6.00,  0.9280)
  ! TABLE  <- used to separate between lists ->
  ! < next list >
  ! END  <- used to mark end of input ->
  !
  ! -------------------------------------------------------------------
  
100 nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt='(A)') sitem1
  if (ioerr/=0) goto 204
  if (INDEX(trim(sitem1),"TABLE")>0) goto 101
  goto 104 ! stop input, unknown or END situation
  
101 FST_nat = FST_nat + 1
  call FST_message("")
  call FST_message("FST: Reading next x-ray scattering factor table:")
  nitem = 0
 
  ! - symbol used to identify the atom/ion type
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt='(A)') ssym
  call FST_message("- Symbol: "//trim(ssym))
  if (ioerr/=0) goto 204
  if (len_trim(ssym)==0) goto 204
  write(unit=FST_sym(FST_nat),fmt='(A)') trim(ssym)
  !
  
  ! - core charge Z
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt=*) fz
  write(unit=FST_smsg,fmt='("-  Z = ",I3)') nint(fz)
  call FST_message(trim(FST_smsg))
  if (ioerr/=0) goto 204
  if (fz<1.0 .or. fz>250.0) goto 204
  FST_crg(2,FST_nat) = fz
  !
  
  ! - ionic charge dZ
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt=*) fdz
  write(unit=FST_smsg,fmt='("- dZ = ",F7.4)') fdz
  call FST_message(trim(FST_smsg))
  if (ioerr/=0) goto 204
  if (fdz<-10.0 .or. fdz>10.0) goto 204
  FST_crg(3,FST_nat) = fdz
  !
  
  ! - clear the temporary fx table
  tabFX = 0.0
  ! - read the next table item
102 nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt='(A)') sitem1
  if (INDEX(trim(sitem1),"TABLE")>0) goto 103 ! next table starts, interpret current table
  if (INDEX(trim(sitem1),"END")>0) goto 103 ! table set ends, interpret current table
  ! next table item needs to be processed
  read(unit=sitem1, iostat=ioerr, fmt=*) fs, fx
  if (ioerr/=0) goto 103 ! assume table is ended
  nitem = nitem + 1
  tabFX(1,nitem) = fs
  tabFX(2,nitem) = fx
  goto 102 ! next item
  !
  !
103 nerr = 0
  !
  ! - check if the table is long enough for this program
  !   we require a certain amount of data points to be able
  !   interpolate correctly
  if (nitem<int(FST_prm_num*1.5)) then
    ! the table is too short, do not use it
    if (INDEX(trim(sitem1),"TABLE")>0) goto 101 ! proceed with next table
    if (INDEX(trim(sitem1),"END")>0) goto 104 ! finish input
  end if
  
  ! table reading has ended here
  if (FST_domsg>1) then
  ! - report fx table
    call FST_Message("- Loaded table of x-ray scattering factors:")
    call FST_Message("  s [1/A]         fx [A]")
    do i=1, nitem
      write(unit=FST_smsg,fmt='("  ",F7.2,"   ",F12.6)') &
     &      tabFX(1,i), tabFX(2,i)
      call FST_Message(trim(FST_smsg))
    end do
  end if
  !
  call FST_Message("Translating table to electron scattering factors:")
  FST_feptmp = 0.0
  ! - interpret the table, translate to electron scattering factors using the
  !   Mott formula
  call FST_FX2FE( tabFX(1:2,1:nitem), FST_fel(1:2,1:nitem, FST_nat), &
     &            fz, fdz, nitem, nerr)
  if (nerr/=0) goto 206
  !
  ! The FE-table is correctly created at this point.
  ! - save the table length
  FST_crg(4,FST_nat) = real(nitem)
  
  !
  call FST_Message("Fitting parameterization to table:")
  ! - further analyse the electron scattering factors by fitting the
  !   Kirkland parameterization to the screened part saved in FST_fel
  call FST_FITFE( FST_fel(1:2, 1:nitem, FST_nat), nitem, fz, fdz, &
     &            FST_fep(1:FST_prm_num, FST_nat), &
     &            FST_per(1:3, FST_nat), nerr )
  if (nerr/=0) goto 207
  !
  ! The FE-parameterization is finished.
  ! - flag the table to be usable
  FST_crg(1,FST_nat) = 1.0
  ! - post a short report
  call FST_message("- Best fit parameters:")
  ! - - atom type
  write(unit=smsg,fmt='(A,": Z = ",I3,", dZ = ",F5.2)') &
     &  trim(FST_sym(FST_nat)), nint(FST_crg(2,FST_nat)), &
     &  FST_crg(3,FST_nat)
  call FST_message("  Processed atom: "//trim(smsg))
  ! - - best fitting parameters
  ! - - - Ai
  call FST_message("  A1             B1             C1             D1")
  write(unit=smsg, &
     &       fmt='("  ",E13.6,"  ",E13.6,"  ",E13.6,"  ",E13.6)') &
     &       FST_fep(1:4, FST_nat)
  call FST_message(trim(smsg))
  ! - - - Bi
  call FST_message("  A2             B2             C2             D2")
  write(unit=smsg, &
     &       fmt='("  ",E13.6,"  ",E13.6,"  ",E13.6,"  ",E13.6)') &
     &       FST_fep(5:8, FST_nat)
  call FST_message(trim(smsg))
  ! - - - Ci
  call FST_message("  A3             B3             C3             D3")
  write(unit=smsg, &
     &       fmt='("  ",E13.6,"  ",E13.6,"  ",E13.6,"  ",E13.6)') &
     &       FST_fep(9:12, FST_nat)
  call FST_message(trim(smsg))
  ! - - DC value
  write(unit=smsg, &
     &       fmt='("  f_el(0) = ",F8.4)') &
     &       100.0*FST_per(1, FST_nat)/FST_per(2, FST_nat)
  call FST_message(trim(smsg))
  ! - - residuals
  write(unit=smsg, &
     &       fmt='("  S = ",G13.6,", E = ",G13.6,", R = ",G13.6)') &
     &  FST_per(1:3, FST_nat)
  call FST_message(trim(smsg))
  !
  !
  ! CONTINUE DEPENDING ON LAST LINE CONTENT
  if (INDEX(trim(sitem1),"TABLE")>0) goto 101 ! proceed with next table
  if (INDEX(trim(sitem1),"END")>0) goto 104 ! finish input
  
  ! -------------------------------------------------------------------
  !
  ! close the file
104 close(unit=lun, iostat=ioerr)
  
  FST_sinfo = "FST: finished loading x-ray scattering factors."
  call FST_Message(trim(FST_sinfo))
  
  return
  
201 nerr = 1
    call FST_ERROR("Failed to load x-ray scattering factors. File [" &
     & //trim(sfile)//"] was not found.",subid+nerr)
    return
    
202 nerr = 2
    call FST_ERROR("Failed to load x-ray scattering factors. " &
     & //"No free logical file unit is available.",subid+nerr)
    return
    
203 nerr = 3
    call FST_ERROR("Failed to load x-ray scattering factors. File [" &
     & //trim(sfile)//"] could not be opened for reading.",subid+nerr)
    return
    
204 nerr = 4
    write(unit=snum,fmt='(I)') nline
    call FST_ERROR("Failed to read data from file ["//trim(sfile)// &
     & "] in line "//trim(adjustl(snum))//".",subid+nerr)
    return
    
205 nerr = 5
    write(unit=snum,fmt='(I)') nline
    call FST_ERROR("Insufficient FX-table length in file ["// &
     & trim(sfile)//"], line "//trim(adjustl(snum))//".",subid+nerr)
    return
    
206 nerr = 6
    call FST_ERROR("Failure while translating from FX-table to "// &
     & "FE-table for atom type ["//trim(FST_sym(FST_nat))//"].", &
     & subid+nerr)
    return

207 nerr = 7
    call FST_ERROR("Failure while fitting the parameterization to "// &
     & "FE-table for atom type ["//trim(FST_sym(FST_nat))//"].", &
     & subid+nerr)
    return
  
end subroutine FST_LoadFX




!**********************************************************************!
!
! subroutine FST_LoadFEP
!
! loads paramaters for electron scattering factors from file
!
subroutine FST_LoadFEP(sfile, nerr)

  use fitfeprm

  implicit none
  
  integer*4, parameter :: subid = 600 ! error id
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  
  logical :: fexists
  integer*4 :: lun, ioerr, nline, nitem, i, np
  real*4 :: fz, fdz, fep(FFE_PRM_MAX), ffe0, ffed(FFE_PRM_MAX), asym
  character(len=FST_ll) :: sitem1, sitem2, snum, ssym, smsg
  
  ! initialize
  nerr = 0
  nline = 0
  nitem = 0
  np = FFE_PRM_MAX
  
  ! clear old tables
  FST_nat = 0
  FST_sym = ""
  FST_crg = 0.0
  FST_fel = 0.0
  FST_fep = 0.0
  FST_per = 0.0
  
  ! check if the file exists
  inquire(file=trim(sfile),exist=fexists)
  if (.not.fexists) goto 201
  ! check for available file unit
  lun = FST_GETFREELUN()
  if (lun<0) goto 202
  
  FST_sinfo = "FST: loading parameterized f_el data from file [" &
     & //trim(sfile)//"]."
  call FST_Message(trim(FST_sinfo))
  
  ! open the file for reading
  open( unit=lun, file=trim(sfile), &
      & iostat=ioerr, action='READ', status='OLD' )
  if (ioerr/=0) goto 203
  
  ! READING TABLE HEADER LINE (not used)
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt=*) sitem1
  if (ioerr/=0) goto 204
  
  ! READING TABLES
  !
  ! EXPECTED FORMAT:
  ! TABLE
  ! <Atomic symbol>   (Cu)
  ! <atomic Z number> (29)
  ! <Atomic charge>   (0.00)
  ! <a1, b1, a2, b2>  (  3.58774531e-001  1.06153463e-001  1.76181348e+000  1.01640995e+000)
  ! <a3, b3, c1, d1>  (  6.36905053e-001  1.53659093e+001  7.44930667e-003  3.85345989e-002)
  ! <c2, d2, c3, d3>  (  1.89002347e-001  3.98427790e-001  2.29619589e-001  9.01419843e-001)
  ! TABLE  <- used to separate between lists ->
  ! < next list >
  ! END  <- used to mark end of input ->
  !
  ! -------------------------------------------------------------------
  
100 nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt='(A)') sitem1
  if (ioerr/=0) goto 204
  if (INDEX(trim(sitem1),"TABLE")>0) goto 101
  goto 104 ! stop input, unknown or END situation
  
101 FST_nat = FST_nat + 1
  call FST_message("")
  call FST_message("FST: Reading next parameter table:")
  nitem = 0
 
  ! - symbol used to identify the atom/ion type
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt='(A)') ssym
  call FST_message("- Symbol: "//trim(ssym))
  if (ioerr/=0) goto 204
  if (len_trim(ssym)==0) goto 204
  write(unit=FST_sym(FST_nat),fmt='(A)') trim(ssym)
  !
  
  ! - core charge Z
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt=*) fz
  write(unit=FST_smsg,fmt='("-  Z = ",I3)') nint(fz)
  call FST_message(trim(FST_smsg))
  if (ioerr/=0) goto 204
  if (FST_sym(FST_nat)(1:1)/="L" .and. (fz<1.0 .or. fz>250.0)) goto 204
  FST_crg(2,FST_nat) = fz
  !
  
  ! - ionic charge dZ
  nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt=*) fdz
  write(unit=FST_smsg,fmt='("- dZ = ",F7.4)') fdz
  call FST_message(trim(FST_smsg))
  if (ioerr/=0) goto 204
  if (fdz<-10.0 .or. fdz>10.0) goto 204
  FST_crg(3,FST_nat) = fdz
  !
  
  ! - clear the temporary fe parameters
  fep = 0.0
  ! - read the next table line
102 nline = nline + 1
  read(unit=lun, iostat=ioerr, fmt='(A)') sitem1
  if (INDEX(trim(sitem1),"TABLE")>0) goto 103 ! next table starts, interpret current parameters
  if (INDEX(trim(sitem1),"END")>0) goto 103 ! table set ends, interpret current parameters
  ! next table item needs to be processed
  read(unit=sitem1, iostat=ioerr, fmt=*) fep(1+4*nitem:4+4*nitem)
  if (ioerr/=0) goto 103 ! assume table is ended
  nitem = nitem + 1
  goto 102 ! next item
  !
  !
103 nerr = 0
  !
  ! - transfer and report fel parameters
  FST_fep(1:2, FST_nat) = fep(1:2)
  FST_fep(5:6, FST_nat) = fep(3:4)
  FST_fep(9:10, FST_nat) = fep(5:6)
  FST_fep(3:4, FST_nat) = fep(7:8)
  FST_fep(7:8, FST_nat) = fep(9:10)
  FST_fep(11:12, FST_nat) = fep(11:12)
  ! table reading has ended here, do some checks
  call FST_Message("- Loaded f_el parameters:")
  write(unit=smsg,fmt='(A,": Z = ",I3,", dZ = ",F5.2)') &
     &  trim(FST_sym(FST_nat)), nint(FST_crg(2,FST_nat)), &
     &  FST_crg(3,FST_nat)
  call FST_message("  - atom: "//trim(smsg))
  call FST_message("  A1             B1             C1             D1")
  write(unit=smsg, &
     &       fmt='("  ",E13.6,"  ",E13.6,"  ",E13.6,"  ",E13.6)') &
     &       FST_fep(1:4, FST_nat)
  call FST_message(trim(smsg))
  call FST_message("  A2             B2             C2             D2")
  write(unit=smsg, &
     &       fmt='("  ",E13.6,"  ",E13.6,"  ",E13.6,"  ",E13.6)') &
     &       FST_fep(5:8, FST_nat)
  call FST_message(trim(smsg))
  call FST_message("  A3             B3             C3             D3")
  write(unit=smsg, &
     &       fmt='("  ",E13.6,"  ",E13.6,"  ",E13.6,"  ",E13.6)') &
     &       FST_fep(9:12, FST_nat)
  call FST_message(trim(smsg))
  call FFE_FEPRM( 0.0, FST_fep(1:np, FST_nat), ffe0, ffed(1:np), np)
  write(unit=smsg, &
     &       fmt='("  - f_el(0) = ",F8.4)') ffe0
  call FST_message(trim(smsg))
  asym = fz - fdz - ( FST_fep(1, FST_nat) + &
                    & FST_fep(5, FST_nat) + &
                    & FST_fep(9, FST_nat) ) / FFE_PREFAC
  write(unit=smsg, &
     &       fmt='("  - physical consistency: ",E13.6)') asym
  call FST_message(trim(smsg))
  
  if (FST_sym(FST_nat)(1:1)/="L" .and. abs(asym)>0.1) then
    call FST_Warning("Loaded scattering factor parameterization "//&
     &               "shows large deviation from expected physics.")
  end if
  
  ! - flag the table to be usable
  FST_crg(1,FST_nat) = 1.0
  !
  ! CONTINUE DEPENDING ON LAST LINE CONTENT
  if (INDEX(trim(sitem1),"TABLE")>0) goto 101 ! proceed with next table
  if (INDEX(trim(sitem1),"END")>0) goto 104 ! finish input
  
  ! -------------------------------------------------------------------
  !
  ! close the file
104 close(unit=lun, iostat=ioerr)
  
  FST_sinfo = "FST: finished loading parameterized f_el data."
  call FST_Message(trim(FST_sinfo))
  
  return
  
201 nerr = 1
    call FST_ERROR("Failed to load f_el parameters. File [" &
     & //trim(sfile)//"] was not found.",subid+nerr)
    return
    
202 nerr = 2
    call FST_ERROR("Failed to load f_el parameters. " &
     & //"No free logical file unit is available.",subid+nerr)
    return
    
203 nerr = 3
    call FST_ERROR("Failed to load f_el parameters. File [" &
     & //trim(sfile)//"] could not be opened for reading.",subid+nerr)
    return
    
204 nerr = 4
    write(unit=snum,fmt='(I)') nline
    call FST_ERROR("Failed to read data from file ["//trim(sfile)// &
     & "] in line "//trim(adjustl(snum))//".",subid+nerr)
    return
    
end subroutine FST_LoadFEP



!**********************************************************************!
!
! subroutine FST_SaveFEP
!
! stores parameters for electron scattering factors to file
!
subroutine FST_SaveFEP(sfile, nerr, noverride)

  use fitfeprm

  implicit none
  
  integer*4, parameter :: subid = 700 ! error id
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(in), optional :: noverride
  integer*4, intent(inout) :: nerr
  
  logical :: fexists
  integer*4 :: lun, ioerr, nline, nitem, i, nov, np
  real*4 :: fz, fdz, fep(FFE_PRM_MAX), ffe0, ffed(FFE_PRM_MAX), asym
  character(len=FST_ll) :: sitem1, sitem2, snum, ssym, smsg
  
  ! initialize
  nerr = 0
  nitem = 0
  np = FFE_PRM_MAX
  nov = 0
  if (present(noverride)) nov = noverride
  
  ! check if the file already exists
  inquire(file=trim(sfile),exist=fexists)
  if ((nov==0).and.fexists) goto 201
  ! check for available file unit
  lun = FST_GETFREELUN()
  if (lun<0) goto 202
  
  if (fexists) call FST_Warning("Replacing content of existing file.")
  
  FST_sinfo = "FST: saving parameterized f_el data to file [" &
     & //trim(sfile)//"]."
  call FST_Message(trim(FST_sinfo))
  
  ! open the file for reading
  call createfilefolder(trim(sfile),ioerr)
  open( unit=lun, file=trim(sfile), &
      & iostat=ioerr, action='WRITE', status='REPLACE' )
  if (ioerr/=0) goto 203
  
  ! HEADER LINE (information string)
  nline = nline + 1
  write(unit=lun, iostat=ioerr, fmt='(A)') &
     & "# Electron scattering factors, parameterization by "// &
     & "E. Kirkland, calculated by CELSLC, J. Barthel (2013)"
  if (ioerr/=0) goto 204
  
  ! WRITING TABLES
  !
  ! EXPECTED FORMAT:
  ! TABLE
  ! <Atomic symbol>   (Cu)
  ! <atomic Z number> (29)
  ! <Atomic charge>   (0.00)
  ! <a1, b1, a2, b2>  (  3.58774531e-001  1.06153463e-001  1.76181348e+000  1.01640995e+000)
  ! <a3, b3, c1, d1>  (  6.36905053e-001  1.53659093e+001  7.44930667e-003  3.85345989e-002)
  ! <c2, d2, c3, d3>  (  1.89002347e-001  3.98427790e-001  2.29619589e-001  9.01419843e-001)
  ! TABLE  <- used to separate between lists ->
  ! < next list >
  ! END  <- used to mark end of input ->
  !
  ! -------------------------------------------------------------------
  
  if (FST_nat>0) then
    do i=1, FST_nat
      if (FST_crg(1,FST_nat) == 0.0) cycle ! disabled table
      ! - write table init with SIGMA info
      write(unit=lun, fmt='(A,F9.6)') &
     &     "TABLE    # SIGMA = ", FST_per(1, i)
      ! - write atom or ion symbol
      write(unit=lun, fmt='(A)') trim(FST_sym(i))
      ! - write atomic number Z
      write(unit=lun, fmt='(I3)') nint(FST_crg(2,i))
      ! - write ionic charge dZ
      write(unit=lun, fmt='(F7.3)') FST_crg(3,i)
      ! - write parameter line 1: a1  b1  a2  b2
      write(unit=lun, fmt='(4E17.9)') FST_fep(1:2,i), FST_fep(5:6,i)
      ! - write parameter line 2: a3  b3  c1  d1
      write(unit=lun, fmt='(4E17.9)') FST_fep(9:10,i), FST_fep(3:4,i)
      ! - write parameter line 2: c2  d2  c3  d3
      write(unit=lun, fmt='(4E17.9)') FST_fep(7:8,i), FST_fep(11:12,i)
    end do
  end if
  
  ! last line
  write(unit=lun, fmt='(A)') "END"
  !
  ! extra lines for documenting the format
  write(unit=lun, fmt='(A)') ""
  write(unit=lun, fmt='(A)') ""
  write(unit=lun, fmt='(A)') &
     & "# LEAVE TWO LINES EMPTY ABOVE TO CANCEL INPUT FEP READING"
  write(unit=lun, fmt='(A)') &
     & "# The parameter table form is similar to that used by "
  write(unit=lun, fmt='(A)') &
     & "# E. Kirkland in 'Advanced Computing in Electron Microscopy'."
  write(unit=lun, fmt='(A)') &
     & "# Information is added to distinguish neutral atoms from ions."
  write(unit=lun, fmt='(A)') &
     & "# A parameter table for one atom or ion has the following form:"
  write(unit=lun, fmt='(A)') "TABLE"
  write(unit=lun, fmt='(A)') "Atom symbol"
  write(unit=lun, fmt='(A)') "Z  = atomic number"
  write(unit=lun, fmt='(A)') "dZ = ionic charge"
  write(unit=lun, fmt='(A)') "a1  b1  a2  b2"
  write(unit=lun, fmt='(A)') "a3  b3  c1  d1"
  write(unit=lun, fmt='(A)') "c2  d2  c3  d3"
  write(unit=lun, fmt='(A)') ""
  
  ! -------------------------------------------------------------------
  !
  ! close the file
104 close(unit=lun, iostat=ioerr)
  
  FST_sinfo = "FST: finished loading parameterized f_el data."
  call FST_Message(trim(FST_sinfo))
  
  return
  
201 nerr = 1
    call FST_ERROR("Refused to save f_el parameters. File [" &
     & //trim(sfile)//"] already exists.",subid+nerr)
    return
    
202 nerr = 2
    call FST_ERROR("Failed to save f_el parameters. " &
     & //"No free logical file unit is available.",subid+nerr)
    return
    
203 nerr = 3
    call FST_ERROR("Failed to save f_el parameters. File [" &
     & //trim(sfile)//"] could not be opened for writing.",subid+nerr)
    return
    
204 nerr = 4
    write(unit=snum,fmt='(I)') nline
    call FST_ERROR("Failed to write data to file ["//trim(sfile)// &
     & "] in line "//trim(adjustl(snum))//".",subid+nerr)
    return
    
end subroutine FST_SaveFEP



!**********************************************************************!
!
! subroutine FST_FX2FE
!
! translates X-ray scattering factors to electron scattering factors
! in a table
!
subroutine FST_FX2FE(fx, fe, z, dz, n, nerr)

  use fitfeprm
  
  implicit none
  
  integer*4, parameter :: subid = 300     ! error id
  
  real*4, parameter :: scentral = 0.5     ! central regime of the atomic
                                          ! scattering factors
  real*4, parameter :: seps = 0.0001      ! small s value
  
  real*4, intent(in) :: fx(2,n), z, dz
  real*4, intent(out) :: fe(2,n)
  integer*4, intent(in) :: n
  integer*4, intent(inout) :: nerr
  
  integer*4 :: i, j, nfit
  real*4 :: fits(n), fitfe(n)
  real*4 :: fs, fcprm(6)
  
  ! initialization
  nerr = 0
  fe = 0.0
  !
  ! initial translation, take only data-points with s>0
  j = 0
  do i=1, n
    fs = fx(1,i)
    if (fs>0.0) then
      j = j + 1
      fe(1,j) = fs
      fe(2,j) = FFE_PREFAC * (z-dz-fx(2,i)) / (fs*fs)
    end if
  end do
  !
!  !
!  ! preset fe(0) with a very rough but secure method
!  !   f(0) = f(0.05) + [ f(0.05)-f(0.1) ] / 2
!  fe(1,1) = 0.0
!  fe(2,1) = (3.0*fe(2,2) - fe(2,3)) * 0.5
!  !
!  ! interpolate 3 Gaussians ( 6 parameters )
!  !   fit(s) = sum_i=1,3 A_i * Exp[ -B_i s^2 ]
!  ! in the low s regime ( 0 > s > 0.5 /A ) of fe
!  ! and take fe(0) = fit(0) = sum_i=1,3 A_i from
!  ! this interpolation
!  !
!  ! - prepare (extract) the input data of the fit
!  nfit = 0
!  fits = 0.0
!  fitfe = 0.0
!  do i=2, n
!    if (fe(1,i)>seps .and. fe(1,i)<scentral+seps) then
!      nfit = nfit + 1
!      fits(nfit) = fe(1,i)
!      fitfe(nfit) = fe(2,i)
!    end if
!  end do
!  ! - do the fit (special subroutine)
!  call FST_FITCENTRAL(fits(1:nfit),fitfe(1:nfit),nfit,fcprm(1:6),nerr)
!  !
!  ! - backup temporary parameters
!  FST_feptmp(1:6) = fcprm(1:6)
!  !
!  ! - calculate fe(0)
!  fe(2,1) = sum(fcprm(1:3))
  !
  if (FST_domsg>1) then
    ! - report fe table
    call FST_Message("- Translated table of electron scattering factors:")
    call FST_Message("  s [1/A]         fe [A]")
    do i=1, j
      write(unit=FST_smsg,fmt='("  ",F7.2,"   ",F12.6)') fe(1,i), fe(2,i)
      call FST_Message(trim(FST_smsg))
    end do
    !
  end if
  !
  return
  

end subroutine FST_FX2FE





!**********************************************************************!
!
! subroutine FST_FITFE
!
! Fit a parameterization of 3 Gaussians and 3 Lorentzians
!   (according to E. Kirkland)
! to a table of electron scattering factors
!
! model function used externally:
! fe(s) = SUM_i=1,3 [ Ai / ( s^2 + Bi ) + Ci * Exp( -Di * s^2 ) ]
!
! OUTPUT ARRAY ORGANIZATION
! fep = (A1, B1, C1, D1, A2, B2, C2, D2, ..., AN, BN, CN, DN )
!
! actual model used in fitting:
! fe(s) = SUM_i=1,3 [ ai^2 / ( s^2 + bi^2 ) + ci^2 * Exp( - di^2 * s^2 ) ]
!
subroutine FST_FITFE(sfe, n, z, dz, fep, per, nerr)

  USE fitfeprm

  implicit none
  
  integer*4, parameter :: subid = 500     ! error id
  integer*4, parameter :: prm_num = FFE_PRM_MAX ! number of fit parameters
  integer*4, parameter :: imp_max = 5     ! number of fit parameters
  integer*4, parameter :: rep_max = 100   ! max. number of repeats to improvement
  integer*4, parameter :: try_min = 100   ! min. number of trials
  real*4, parameter :: inirng = 2.0       ! initial search range
  
  real*4, intent(in) :: sfe(2,n), z, dz
  integer*4, intent(in) :: n
  real*4, intent(inout) :: fep(prm_num), per(3)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: i, j, nfit, nprm, uprm(prm_num), niter, nrep, nimp, ntry
  real*4 :: fits(n), fitfe(n), fitfes(n), rftol, tfemax
  real*4 :: prm(prm_num), lprm(2,prm_num), covprm(prm_num,prm_num)
  real*4 :: fs, fe, chisqr, ffe, fedp(prm_num), sigma, Eval, Rval
  real*4 :: Rden, csqold, dfe, ffe0, maxB, tfep(prm_num)
  
  external :: InitRand
  real*4, external :: UniRand
  
  !
  ! initialize
  nerr = 0
  fep = 0.0
  per = 0.0
  !
  nprm = prm_num
  nfit = 0
  uprm = 1
  prm = 0.0
  tfep = 0.0
  lprm = 0.0
  covprm = 0.0
  !
  ! - transfer the data for the extra side condition
  FFE_ZRED = z - dz ! strength of the x-ray decay at large s
  FFE_W_EXTRA = 100.0/(FFE_ZRED*FFE_ZRED) ! related fit weight, (strong)
  ! - transfer the input data to extra arrays, using only fe's at s>0
  nfit = 0
  do i=1, n
    fs = sfe(1,i)
    fe = sfe(2,i)
    if (abs(fs)>0.0) then
      nfit = nfit + 1
      fits(nfit) = fs
      fitfe(nfit) = fe
      fitfes(nfit) = 1.0/(fe*fe) ! (normal)
    end if
  end do
  ! - chek number of data points
  if (nfit<2*nprm) then
    if (nfit<nint(1.25*real(nprm))) then
      nerr = 1
      call FST_ERROR("Insufficient data for fitting parameterization.",subid+nerr)
      return
    end if
    call FST_Warning("Weak data supply for fitting parameterization.")
  end if
  ! - tfemax = max. table value
  tfemax = maxval(fitfe(1:nfit))
  !
  ! - Setup the initial parameter values (random)
  ! - - parameters (Ai,Bi,Ci,Di)
  do i=1, nprm
    prm(i) = 0.001 + inirng*UniRand()
  end do
  !
  ! - Initialize the fit.
  !call InitRand()
  !
  FFE_SHOW_ERR = 0
  ! - - setup the table s values
  call FFE_SET_S(nfit,fits(1:nfit))
  ! - - setup the table fe values
  call FFE_SET_F(nfit,fitfe(1:nfit))
  ! - - setup the chi-square weight table
  FFE_USE_WEIGHTS = 1 ! activate weights
  call FFE_SET_W(nfit,fitfes(1:nfit))
  ! - Setup the parameter range limits
  FFE_USE_RANGE = 1 ! activate the use of range limits
  ! - - Limit all parameters to the range +/- 100.
  !     This range is more than sufficient for the
  !     Lorentzian and Gaussian amplitudes, which
  !     are realized by squares of the fit parameters.
  !     Limiting the width parameters to +/- 100
  !     seems ok, since the interesting range is on the
  !     order of 10. Here we may se room for a change,
  !     using a larger search range for the width
  !     parameters.
  do i=1, nprm
    call FFE_SET_PRNG(i,-100.0,100.0)
  end do
  ! - number of iterations
  niter = 0
  ! - tolerance limit
  rftol = 1.0E-12
  ! - fractional temperature reduction parameter for
  !   simulated annealing algorithms
  FFE_SAS_TRED = 0.1
  ! - number of trials for 
  !   simulated annealing algorithms
  FFE_SAS_NTRY = 200
  ! - chi-square test number
  csqold = -1.0
  ! - number of fit repeats
  nrep = 0
  ! - number of improvements
  nimp = 0
  ! - number of trials
  ntry = 0
  ! DO THE FIT !
  do
  ! - try new chi-square minimization
    call FFE_FINDMIN(prm(1:nprm), niter, rftol, 2, chisqr ) ! restart - random
    !write(*,fmt='(A,I5,A,E13.6)') " > ", niter, " iterations: ", sqrt(chisqr)
    call FFE_FINDMIN(prm(1:nprm), niter, rftol, 2, chisqr ) ! restart - local
    !write(*,fmt='(A,I5,A,E13.6,$)') " > ", niter, " iterations: ", sqrt(chisqr)
    ntry = ntry + 1
  ! - check for improvement ...
    if (csqold>0.0) then !  ... compare to previous chi-square
      nrep = nrep + 1 ! increase number of repeat counts
      maxB = prm(2)**2
      maxB = max(maxB,1.0/prm(4)**2)
      maxB = max(maxB,prm(6)**2)
      maxB = max(maxB,1.0/prm(8)**2)
      maxB = max(maxB,prm(10)**2)
      maxB = max(maxB,1.0/prm(12)**2)
      if (chisqr<csqold .and. maxB<200.0) then ! ... found better chi-square -> save the new chi-square and try again
        nimp = nimp + 1
        nrep = 0
        tfep(1:nprm) = prm(1:nprm)
        csqold = chisqr
        write(*,fmt='(A,I2,$)') " - improvement #", nimp
        write(*,fmt='(A,E13.6,A,F8.3)') " CHI: ",sqrt(chisqr), &
          & ", Amp-Consist: ",z-dz-(prm(1)**2+prm(5)**2+prm(9)**2)/FFE_PREFAC
      else                    ! ... not a better chi-square -> set back to old and try again or quit
        prm(1:nprm) = tfep(1:nprm)
        !write(*,fmt='(A)') " REJECTED."
      end if
    else                 !  ... no previous chi-square -> save the new chi-square and try again
      tfep(1:nprm) = prm(1:nprm)
      csqold = chisqr
      write(*,fmt='(A,$)')      " - first run ....."
      write(*,fmt='(A,E13.6,A,F8.3)') " CHI: ",sqrt(chisqr), &
        & ", Amp-Consist: ",z-dz-(prm(1)**2+prm(5)**2+prm(9)**2)/FFE_PREFAC
    end if
    if (nimp>=imp_max.and.ntry>=try_min) goto 101 ! exit, max. number of improvements reached
    if (nrep>rep_max) goto 101 ! exit, max. number of repetitions reached
  ! - want to try again, use a new set of initial parameters (random)
    do i=1, nprm
      prm(i) = 0.001 + inirng*UniRand()
    end do
  end do
  !
  !
  ! - Get the residual values, sigma, E and R.
  ! - - sigma and Rval
101 sigma = 0.0
  Rval = 0.0
  Rden = 0.0
  fep = tfep**2
  ! - Get the fe at s=0
  call FFE_FEPRM(0.0, fep(1:nprm), ffe0, fedp(1:nprm), nprm)
  ! - Get test data for all other s
  do i=1, nfit
    call FFE_FEPRM(fits(i), fep(1:nprm), ffe, fedp(1:nprm), nprm)
    dfe = fitfe(i) - ffe
    sigma = sigma + dfe*dfe
    Rval = Rval + abs(dfe)
    Rden = Rden + abs(fitfe(i))
  end do
  sigma = sqrt( sigma / real(nfit) )
  Rval = Rval / Rden
  ! - - E
  Eval = 100.0 * sigma / ffe0
  ! - - store the residuals
  per(1) = sigma
  per(2) = Eval
  per(3) = Rval
  
  return
  
end subroutine FST_FITFE


!**********************************************************************!
!
! integer*4 function FST_GETIDX
!
! returns the list index of atom or ion data 
! identified by its symbol
!
! return value = 0 indicates a failure in identifying the
!                  input atomic symbol string.
!
function FST_GETIDX(sym)
  use fitfeprm
  implicit none
  character(len=*), intent(in) :: sym
  integer*4 :: FST_GETIDX
  integer*4 :: i, j, nat
  FST_GETIDX = 0 ! preset with fail
  j = 0
  nat = FST_nat ! get active number of atom definitions
  if (nat>0) then
    find_at: do i=1, nat ! find atom that has equal symbol
      if (trim(FST_sym(i))==trim(sym)) then
        j = i
        exit find_at
      end if
    end do find_at
    FST_GETIDX = j ! copy to function return
  end if
end function FST_GETIDX

!**********************************************************************!
!
! function FST_GETFE
!
! returns the electron scattering factor for an atom or ion
! identified by its symbol, for a given scattering angle s [1/A]
!
! a negative value returned denotes a failure in identifying the
! input atomic symbol string.
!
function FST_GETFE(sym,s)
  
  use fitfeprm
  
  implicit none
  
  character(len=*), intent(in) :: sym
  real*4, intent(in) :: s
  
  real*4 :: FST_GETFE
  
  integer*4 :: j, np, nat
  real*4 :: fe, dz, feio
  
  FST_GETFE = -1.0 ! preset with fail
  fe = -1.0
  np = FFE_PRM_MAX ! get number of function parameters from module
  j = FST_GETIDX(sym)
  feio = 0.0
  if (j>0) then ! found an entry in the list
    ! get the scattering factor
    call FFE_FEPRMY(s, FST_fep(1:np,j), fe, np)
    ! get ionic charge
    !dz = FST_crg(3,j)
    !if (abs(s)>0.0) feio = dz * FFE_PREFAC / (s*s)
    !
    FST_GETFE = fe ! + feio ! function return value
  end if
  
end function FST_GETFE



!**********************************************************************!
!
! function FST_GETMU0(s)
!
! returns the integrand for the absorptive form factor at g=0
!
! The atom type is identified by an index FST_PRM_J
! The atom Debye-Waller parameter is FST_PRM_J_BISO
! Set the values of these parameters before using this function!
!
function FST_GETMU0(s)
  
  use fitfeprm
  
  implicit none
  
  real*8, parameter :: twopi = 6.283185307179586476925286766559
  
  real*8, intent(in) :: s
  real*8 :: FST_GETMU0
  real*8 :: fs, dwfs, biso2
  
  fs = FST_GETSCAR1(s)
  biso2 = dble(2.0*FST_PRM_J_BISO)
  dwfs = dexp( -biso2*s*s )
  FST_GETMU0 = twopi * fs*fs * (1.d+0 - dwfs) * s
  
  return
    
end function FST_GETMU0


!**********************************************************************!
!
! function FST_GETMUG(theta,phi)
!
! returns the integrand for the absorptive form factor at g=FST_PRM_G
! along scattering path (theta,phi)
!
! The atom type is identified by an index FST_PRM_J
! The incident wave vector is given by FST_PRM_K0
! The atom Debye-Waller parameter is FST_PRM_J_BISO
! Set the values of these parameters before using this function!
!
function FST_GETMUG(theta,phi)
  
  use fitfeprm
  
  implicit none
  
  real*8, parameter :: twopi = 6.283185307179586476925286766559
  
  real*8, intent(in) :: theta, phi
  real*8 :: FST_GETMUG
  real*8 :: k, q, g, qg, biso, ct, st, cp, twok, sg, sq, sqg
  real*8 :: fq, fqg, dwfg, dwfq, dwfqg
  
  ct = dcos(theta)
  st = dsin(theta)
  cp = dcos(phi)
  ! q = scattering vector length of Q in the ewald sphere
  ! Q = (qx, qy, qz) = K' - K
  ! qx = k * sin(theta) * cos(phi)
  ! qy = k * sin(theta) * sin(phi)
  ! qz = k * ( cos(theta) - 1 )
  ! k = |K| = wekok
  k = dble(FST_PRM_K0)
  twok = k+k
  ! q = Sqrt[ 2 k^2 (1 - Cos[q]) ]
  q = k * dsqrt( 2.D+0 - 2.D+0 * ct )
  ! g = length of some reciprocal space vector G = (gx, gy, gz)
  !     with gx = g, gy = 0, gz = 0 (obda)
  g = dble(FST_PRM_G)
  ! qg = length of the difference vector Q - G
  qg = dsqrt( g*g + twok*k - twok*(k*ct + g*cp*st) )
  biso = dble(FST_PRM_J_BISO)
  sg = g * 0.5D+0
  sq = q * 0.5D+0
  sqg = qg * 0.5D+0
  fq = FST_GETSCAR1(sq) ! f(s) --- f(g/2)
  fqg = FST_GETSCAR1(sqg)
  dwfg = dexp( -biso*sg*sg ) ! exp( -1/4 biso*g^2 ) = exp( -biso*s^2)
  dwfq = dexp( -biso*sq*sq )
  dwfqg = dexp( -biso*sqg*sqg )
  FST_GETMUG = fq*fqg * (dwfg - dwfq*dwfqg) * st
  
  return
    
end function FST_GETMUG



!**********************************************************************!
!
! function FST_GETF2G(theta,phi)
!
! returns the integrand for the squared form factor
! along scattering path (theta)
!
! The atom type is identified by an index FST_PRM_J
! The incident wave vector is given by FST_PRM_K0
! The atom Debye-Waller parameter is FST_PRM_J_BISO
! Set the values of these parameters before using this function!
!
function FST_GETF2G(theta)
  
  use fitfeprm
  
  implicit none
  
  real*8, intent(in) :: theta
  real*8 :: FST_GETF2G
  real*8 :: k, q, g, qg, biso, ct, st, sq
  real*8 :: fq, fqg, dwfg, dwfq, dwfqg
  
  ct = dcos(theta)
  st = dsin(theta)
  ! q = scattering vector length of Q in the ewald sphere
  ! Q = (qx, qy, qz) = K' - K
  ! qx = k * sin(theta) * cos(phi)
  ! qy = k * sin(theta) * sin(phi)
  ! qz = k * ( cos(theta) - 1 )
  ! k = |K| = FST_PRM_K0
  k = dble(FST_PRM_K0)
  ! q = Sqrt[ 2 k^2 (1 - Cos[q]) ]
  q = k * dsqrt( 2.D+0 - 2.D+0 * ct )
  biso = dble(abs(FST_PRM_J_BISO))
  sq = q * 0.5D+0
  fq = FST_GETSCAR1(sq) ! f(s) --- f(g/2)
  if ( biso>0.D+0 ) then
    fq = fq * dexp( -biso*sq*sq ) ! exp( -1/4 biso*q^2 ) = exp( -biso*s^2)
  end if
  FST_GETF2G = st*fq*fq ! sin(theta) * ( f(theta)*DWF(theta) )^2
  
  return
    
end function FST_GETF2G



!**********************************************************************!
!
! function FST_GETSCAR1(s)
!
! returns the atomic form factor in [A] depending on
! - diffraction vector s = g/2 [1/A]
! (no relativistic correction applied)
!
! the atom type is identified by an index FST_PRM_J
! Set the value of this parameter before using this function!
!
! DOUBLE PRECISION IMPLEMENTATION
!
function FST_GETSCAR1(s)
  
  use fitfeprm
  
  implicit none
  
  real*8, parameter :: alph2 = 4.0d-04 ! square of inverse yukawa range in A (only used for small s2 in ions) 
  
  real*8, intent(in) :: s
  real*8 :: FST_GETSCAR1
  integer*4 :: i, j, np, nat
  real*4 :: fe, sf
  real*8 :: c1, dz, feio
  
  FST_GETSCAR1 = 0.d+0 ! preset
  j = FST_PRM_J ! get the list index for the input symbol
  c1 = dble(FFE_PREFAC) ! FFE_PREFAC = 0.0239337 = m0 * e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) [A]
  !
  if (j>0) then ! found an entry in the list
    !
    ! init
    fe = 0.0
    feio = 0.d+0
    np = FFE_PRM_MAX ! get number of function parameters from module
    !
    sf = real(s, kind=4 ) ! translate to interface form
	call FFE_FEPRMY(sf, FST_fep(1:np,j), fe, np)
	!
	if (dabs(dz)>0.d+0) then ! an ionic charge potential needs to be added
      dz = dble(FST_crg(3,j))
      feio = dz * c1 / (s*s+alph2)
    end if
    !
	FST_GETSCAR1 = dble(fe) + feio
	!
  end if
  
end function FST_GETSCAR1


!**********************************************************************!
!
! function FST_GETSCAABS(sym, g, B, ht)
!
! returns the absorptive form factor for an atom or ion
! identified by its symbol depending on
! - diffraction vector g [1/nm]
! - Debye-Waller parameter B [nm^2]
! - high-tension value ht [kV]
!
! result: absorptive atomic form factor [A]
!         contains the screened core potential and the
!         ionic charge coulomb potentials
!
!         no relativistic correction applied
!         need to multiply by gamma^2/(2*pi*k0)
!         no DWF applied
!
function FST_GETSCAABS(sym, g, B, ht)
  
  use fitfeprm
  
  implicit none
    
  character(len=*), intent(in) :: sym
  real*4, intent(in) :: g, B, ht
  
  real*4 :: FST_GETSCAABS
  
  integer*4 :: i, j, np, nat
  real*4 :: ga, k0
  real*4 :: dwasp
  real*8 :: si0, si1, fa, pi0, pi1
  
  !external :: dqsimp, dqsimp2d
  real*8, external :: dsgrid2d
  
  FST_GETSCAABS = 0.0 ! preset
  np = FFE_PRM_MAX ! get number of function parameters from module
  j  = FST_GETIDX(sym) ! get the list index for the input symbol
  FST_PRM_J = j ! set the atomic data list index in the module parameter
  k0 = 0.080655587*sqrt((2.0*FST_elmc2 + ht)*ht) ! = ! k[ht] = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
  !
  if (j>0) then ! found an entry in the list
    !
    ! -------[nm] to [A]
    dwasp = 100.0 * B ! parameter of the DWF [A^2]
    FST_PRM_J_BISO = dwasp
	ga	= 0.1 * g ! scattering angle [1/A]
	FST_PRM_G = ga  ! store current g [1/A] for use in integrator function
	FST_PRM_K0 = k0 ! store current k0 [1/A] for use in integrator function
	! set integration range
	si0 = 0.d+0 ! 0
	si1 = dble(0.5*FST_twopi) ! Pi
	pi0 = 0.d+0 ! 0
	pi1 = dble(0.5*FST_twopi) ! Pi
	fa = dsgrid2d(FST_GETMUG,si0,si1,pi0,pi1,2.0d+0,1.0d+0,128,64) * 2.0D+0
	!
	FST_GETSCAABS = real(fa, kind=4 )
	!
  end if
  
end function FST_GETSCAABS


!**********************************************************************!
!
! function FST_GETSCAF2(sym, tmax, B, ht)
!
! returns the integral of the squared form factor for an atom or ion
! identified by its symbol up to scattering angle tmax [rad]
! depending on
! - Debye-Waller parameter B [nm^2]
! - high-tension value ht [kV]
!
! result: integrated squared atomic form factor [A^2]
!         contains the screened core potential and the
!         ionic charge coulomb potentials
!
!         no relativistic correction applied
!         need to multiply by gamma^2
!
function FST_GETSCAF2(sym, tmax, B, ht)
  
  use fitfeprm
  
  implicit none
  
  integer*4, parameter :: ngrid = 2048
  real*8, parameter :: twopi = 6.283185307179586476925286766559
    
  character(len=*), intent(in) :: sym
  real*4, intent(in) :: tmax, B, ht
  
  real*8 :: FST_GETSCAF2
  
  integer*4 :: i, j, np, nat
  real*4 :: k0
  real*4 :: dwasp
  real*8 :: si0, si1, fa
  
  !external :: dqsimp, dqsimp1d
  real*8, external :: dsgrid1d
  
  FST_GETSCAF2 = 0.0 ! preset
  np = FFE_PRM_MAX ! get number of function parameters from module
  j  = FST_GETIDX(sym) ! get the list index for the input symbol
  FST_PRM_J = j ! set the atomic data list index in the module parameter
  k0 = 0.080655587*sqrt((2.0*FST_elmc2 + ht)*ht) ! = ! k[ht] = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
  !
  if (j>0) then ! found an entry in the list
    !
    ! -------[nm] to [A]
    dwasp = 100.0 * B ! parameter of the DWF [A^2]
    FST_PRM_J_BISO = dwasp
	FST_PRM_G = 0.0 ! G is not used
	FST_PRM_K0 = k0 ! store current k0 [1/A] for use in integrator function
	! set integration range
	si0 = 0.d+0 ! 0
	si1 = min( dble(0.5*FST_twopi), tmax ) ! min( tmax , Pi )
	fa = dsgrid1d(FST_GETF2G,si0,si1,2.0d+0,ngrid) * twopi
	!
	FST_GETSCAF2 = fa
	!
  end if
  
end function FST_GETSCAF2





!**********************************************************************!
!
! function FST_GETSCAC(sym, g, B, ht, dwfflg, absflg)
!
! returns the electron scattering power for an atom or ion
! identified by its symbol depending on
! - diffraction vector g [1/nm]
! - Debye-Waller parameter B [nm^2]
! - high-tension value ht [kV]
! - flag for using the Debye-Waller factor
! - flag for calculating an absorptive factor
!
! result: complex atomic form factor [nm]
!         (real part => elastic potential)
!         (imag part => absorptive potential)
!         contains only the screened core potential
!         ionic charge coulomb potentials should be
!         added when forming the structure factor
!
! a negative value returned in the real part denotes a failure in
! identifying the input atomic symbol string.
!
function FST_GETSCAC(sym, g, B, ht, dwfflg, absflg)
  
  use fitfeprm
  
  implicit none
  
  real*4, parameter :: thrima   = 1.0E-8 ! threshold for imaginary part calculation
  
  character(len=*), intent(in) :: sym
  real*4, intent(in) :: g, B, ht
  logical, intent(in) :: dwfflg, absflg
  
  complex*8 :: FST_GETSCAC
  
  integer*4 :: i, j, np, nat
  real*4 :: sa, fe, fr, fa, fi, feio, k0
  real*4 :: dwas, dwasp, rc !, dz
  
  FST_GETSCAC = -1.0 ! preset with fail
  fe = -1.0
  fa = 0.0 ! absorptive form factor without rel. corrections
  fi = 0.0 ! absorptive form factor for the current ht
  np = FFE_PRM_MAX ! get number of function parameters from module
  j = FST_GETIDX(sym) ! get the list index for the input symbol
  feio = 0.0
  k0 = 0.506774*sqrt((2.0*FST_elmc2 + ht)*ht) ! = ! 2*Pi*k ! andere k-Notation hier: Aufpassen!
  if (j>0) then ! found an entry in the list
    !
    ! -------Debye-Waller Temperature Coeff dwa = 8 pi**2 usa**2
    ! -------[nm] to [ang]
    dwas  = 1.0 ! preset of the Debye-Waller factor (DWF)
    dwasp = 100.0 * B ! parameter of the DWF [A^2]
	sa	= 0.050 * g ! scattering angle (s = g/2) [1/A]
	if (dwfflg) dwas = exp( -dwasp * sa * sa ) ! DWF
	rc	= (FST_elmc2 + ht) / FST_elmc2 ! relativistic correction
	call FFE_FEPRMY(sa, FST_fep(1:np,j), fe, np)
    ! The ionic charge potential is not done here. This will be added later
    ! in the form factor summation.
    ! dz = FST_crg(3,j)
    ! if (abs(sa)>0.0) feio = dz * FFE_PREFAC / (sa*sa) ! FFE_PREFAC = 0.0239337 = m0 * e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) [A]
    ! fe = fe + feio
	fr 	= rc * fe * FST_fourpi * dwas  ! real form factor  dampened by DWF
	                                   ! multiplied by relativistic correction and
	                                   !   4*Pi to be consistent with fscatt.f
    !
    ! absorptive form factor calculation
    if (absflg .and. dwfflg .and. (abs(dwasp)>=thrima) .and. (abs(fr)>=thrima) ) then
      ! calculate the absorptive form factor
      fa = FST_GETSCAABS(sym, g, B, ht)
      fi = rc*rc*fa*k0       ! imaginary form factor
    end if
    ! -------[ang] to [nm]
	FST_GETSCAC	= 0.1 * cmplx( fr, fi )
	!
  end if
  
end function FST_GETSCAC



!**********************************************************************!
!
! function FST_GETSCAR(sym, g, B, ht, dwfflg)
!
! returns the electron scattering power for an atom or ion
! identified by its symbol
! depending on
! - diffraction vector g [1/nm]
! - Debye-Waller parameter B [nm^2]
! - high-tension value ht [kV]
! - flag for using the Debye-Waller factor
!
! a negative value returned denotes a failure in identifying the
! input atomic symbol string.
!
function FST_GETSCAR(sym, g, B, ht, dwfflg)
  
  use fitfeprm
  
  implicit none
  
  character(len=*), intent(in) :: sym
  real*4, intent(in) :: g, B, ht
  logical, intent(in) :: dwfflg
  
  real*4 :: FST_GETSCAR
  
  integer*4 :: i, j, np, nat
  real*4 :: sa, fe, fr, feio
  real*4 :: dwas, dwasp, rc !, dz
  
  FST_GETSCAR = -1.0 ! preset with fail
  fe = -1.0
  np = FFE_PRM_MAX ! get number of function parameters from module
  j = FST_GETIDX(sym) ! get the list index for the input symbol
  feio = 0.0
  if (j>0) then ! found an entry in the list
    !
    ! -------Debye-Waller Temperature Coeff dwa = 8 pi**2 usa**2
    ! -------[nm] to [ang]
    dwas  = 1.0 ! preset of the Debye-Waller factor (DWF)
	dwasp = 100.0 * B ! parameter of the DWF [A^2]
	sa	= 0.050 * g ! scattering angle (s = g/2) [1/A]
	rc	= (FST_elmc2 + ht) / FST_elmc2 ! relativistic correction
	call FFE_FEPRMY(sa, FST_fep(1:np,j), fe, np)
	! The ionic charge potential is not done here. This will be added later
    ! in the form factor summation.
    ! dz = FST_crg(3,j)
    ! if (abs(sa)>0.0) feio = dz * FFE_PREFAC / (sa*sa) ! FFE_PREFAC = 0.0239337 = m0 * e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) [A]
    ! fe = fe + feio
	fr 	= rc * fe * FST_fourpi ! real scattering power
	                           ! multiplied by relativistic correction and
	                           !   4*Pi to be consistent with fscatt.f
    if (dwfflg) then
      dwas = exp( -dwasp * sa * sa ) ! DWF
      fr   = fr  * dwas ! dampened scattering power
    end if
    ! -------[ang] to [nm]
	FST_GETSCAR	= 0.1 * fr
	!
  end if
  
end function FST_GETSCAR




!**********************************************************************!
!
! function FST_GETCRG(sym)
!
! returns the charge assigned to the atom/ion item
! identified by its symbol
!
function FST_GETCRG(sym)
  
  use fitfeprm
  
  implicit none
  
  character(len=*), intent(in) :: sym
  
  real*4 :: FST_GETCRG
  
  integer*4 :: j
  real*4 :: dz
  
  FST_GETCRG = 0.0 ! preset with default
  dz = 0.0
  j = FST_GETIDX(sym) ! get the list index for the input symbol
  if (j>0) then ! found an entry in the list
    !
    dz = FST_crg(3,j)
	FST_GETCRG	= dz
	!
  end if
  
end function FST_GETCRG



!**********************************************************************!
!
! function FST_GETSCADWF(g, B, ht, dwfflg)
!
! returns the debye-waller factor to be applied to a projected potential
! depending on
! - diffraction vector g [1/nm]
! - Debye-Waller parameter B [nm^2]
! - high-tension value ht [kV]
! - flag for using the Debye-Waller factor
!
function FST_GETSCADWF(g, B, ht, dwfflg)
  
  use fitfeprm
  
  implicit none
  
  real*4, intent(in) :: g, B, ht
  logical, intent(in) :: dwfflg
  
  real*4 :: FST_GETSCADWF
  
  real*4 :: dwas
  
  FST_GETSCADWF = 1.0 ! preset with default ( no DWF )
  !
  dwas  = 1.0 ! preset of the Debye-Waller factor (DWF)
  if (dwfflg) then
    dwas = exp( -0.25*B*g*g ) ! DWF
    FST_GETSCADWF = dwas ! dampening factor for scattering power
  end if
  
end function FST_GETSCADWF
















!**********************************************************************!
!
! subroutine FST_FFG3
!
! Calculation of function values and derivatives for a
! function of 3 Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,3 [ ai * Exp( - bi * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, a2, a3, b1, b2, b3)
!
SUBROUTINE FST_FFG3(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: gaussian parameters
!                               a(1) : a1
!                               a(2) : a2
!                               a(3) : a3
!                               a(4) : b1
!                               a(5) : b2
!                               a(6) : b3
!                               a(>6) : useless
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given
! -------------------------------------------------------------------- !

  IMPLICIT NONE

  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: i, j, n
  real*4 :: amp, rwi, fexp, farg

  y = 0.0
  if (na/=6) return
  n = 3

! function value
  do i=1, n
    j = n+i
    amp = a(i)
    rwi = a(j)
    farg = x*x*rwi
    fexp = exp(-farg)
    ! sum up functions
    y = y + amp*fexp
    ! set derivatives
    dyda(i) = fexp
    dyda(j) = -amp*fexp*x*x
  end do

  return

END SUBROUTINE FST_FFG3
!**********************************************************************!



!**********************************************************************!
!
! subroutine FST_FFL3G3
!
! Calculation of function values and derivatives for a
! function of 3 Lorentzians and 3 Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,3 [   ai^2 / ( bi^2 + x^2 )
!                    + ci^2 * Exp( - di^2 * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3)
!
SUBROUTINE FST_FFL3G3(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given == 12
! -------------------------------------------------------------------- !

  IMPLICIT NONE

  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: i, j, k, l, n
  real*4 :: aa, ab, ac, ad, fexp, fearg, flor, flarg

! initialization
  y = 0.0
  if (na/=12) return
  n = 3

! function value calculation
  do i=1, n
    
    ! - iterators
    j =   n+i
    k = 2*n+i
    l = 3*n+i
    ! - current parameters
    aa = a(i)
    ab = a(j)
    ac = a(k)
    ad = a(l)
    ! - current lorentzian values
    flarg = x*x + ab*ab
    flor = 1.0 / flarg
    ! - current gaussian values
    fearg = x*x*ad*ad
    fexp = exp(-fearg)
    
    ! sum up functions
    y = y + aa*aa * flor + ac*ac * fexp
    
    ! sum up derivatives
    dyda(i) = 2.0*aa*flor
    dyda(j) = -2.0*aa*aa*ab*flor*flor
    dyda(k) = 2.0*ac*fexp
    dyda(l) = -2.0*ac*ac*ad*x*x*fexp
    
  end do

  return

END SUBROUTINE FST_FFL3G3
!**********************************************************************!



!**********************************************************************!
!
! subroutine FST_FFLL3G3
!
! Calculation of function values and derivatives for a
! function of 3 Lorentzians and 3 Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,3 [   ai / ( bi^2 + x^2 )
!                    + ci * Exp( - di * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3)
!
SUBROUTINE FST_FFLL3G3(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given == 12
! -------------------------------------------------------------------- !

  IMPLICIT NONE

  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: i, j, k, l, n
  real*4 :: aa, ab, ac, ad, fexp, fearg, flor, flarg

! initialization
  y = 0.0
  if (na/=12) return
  n = 3

! function value calculation
  do i=1, n
    
    ! - iterators
    j =   n+i
    k = 2*n+i
    l = 3*n+i
    ! - current parameters
    aa = a(i)
    ab = a(j)
    ac = a(k)
    ad = a(l)
    ! - current lorentzian values
    flarg = x*x + ab
    flor = 1.0 / flarg
    ! - current gaussian values
    fearg = x*x*ad
    fexp = exp(-fearg)
    
    ! sum up functions
    y = y + aa * flor + ac * fexp
    
    ! sum up derivatives
    dyda(i) = flor
    dyda(j) = aa*flor*flor
    dyda(k) = fexp
    dyda(l) = -ac*x*x*fexp
    
  end do

  return

END SUBROUTINE FST_FFLL3G3
!**********************************************************************!



!**********************************************************************!
!
! subroutine FST_FFLLGN
!
! Calculation of function values and derivatives for a
! function of N Lorentzians and N Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,N [   ai / ( bi + x^2 )
!                    + ci * Exp( - di * x^2 ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, b1, c1, d1, a2, b2, c2, d2, ..., aN, bN, cN, dN )
!
SUBROUTINE FST_FFLLGN(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given == 4*n
! -------------------------------------------------------------------- !

  IMPLICIT NONE

  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: n, i, i1, i2, i3, i4
  real*4 :: aa, ab, ac, ad, fexp, fearg, flor, flarg

! initialization
  y = 0.0
  dyda = 0.0
  n = int(na/4)

! function value calculation
  do i=1, n
    
    ! - iterators
    i1 = 4*(i-1) + 1
    i2 = 4*(i-1) + 2
    i3 = 4*(i-1) + 3
    i4 = 4*(i-1) + 4
    ! - current parameters
    aa = a(i1)
    ab = a(i2)
    ac = a(i3)
    ad = a(i4)
    ! - current lorentzian values
    flarg = x*x + ab
    flor = 1.0 / flarg
    ! - current gaussian values
    fearg = x*x*ad
    fexp = exp(-fearg)
    
    ! sum up functions
    y = y + aa * flor + ac * fexp
    
    ! sum up derivatives
    dyda(i1) = flor
    dyda(i2) = -aa*flor*flor
    dyda(i3) = fexp
    dyda(i4) = -ac*x*x*fexp
    
  end do

  return

END SUBROUTINE FST_FFLLGN
!**********************************************************************!




!**********************************************************************!
!
! subroutine FST_FFQLGN
!
! Calculation of function values and derivatives for a
! function of N Lorentzians and N Gaussians
!
! model function used externally:
! f(x) = SUM_i=1,N [   ai^2 / ( bi^2 + x^2 + ceps)
!                    + ci^2 * Exp( - di^2 * ( x^2 + ceps) ) ]
!
! PARAMETER ARRAY ORGANIZATION
! a = (a1, b1, c1, d1, a2, b2, c2, d2, ..., aN, bN, cN, dN )
!
SUBROUTINE FST_FFQLGN(x, a, y, dyda, na)
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: function parameters
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given == 4*n
! -------------------------------------------------------------------- !

  IMPLICIT NONE
  
  real*4, parameter :: ceps = 1.0E-10

  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: n, i, i1, i2, i3, i4
  real*4 :: aa, ab, ac, ad, fexp, fearg, flor, flarg

! initialization
  y = 0.0
  dyda = 0.0
  n = int(na/4)

! function value calculation
  do i=1, n
    
    ! - iterators
    i1 = 4*(i-1) + 1
    i2 = 4*(i-1) + 2
    i3 = 4*(i-1) + 3
    i4 = 4*(i-1) + 4
    ! - current parameters
    aa = a(i1)
    ab = a(i2)
    ac = a(i3)
    ad = a(i4)
    ! - current lorentzian values
    flarg = x*x + ab*ab + ceps
    flor = 1.0 / flarg
    ! - current gaussian values
    fearg = x*x * (ceps+ad*ad)
    fexp = exp(-fearg)
    
    ! sum up functions
    y = y + aa*aa * flor + ac*ac * fexp
    
    ! sum up derivatives
    dyda(i1) = 2.0*aa*flor
    dyda(i2) = -2.0*aa*aa*ab*flor*flor
    dyda(i3) = 2.0*ac*fexp
    dyda(i4) = -2.0*ac*ac*ad*x*x*fexp
    
  end do

  return

END SUBROUTINE FST_FFQLGN
!**********************************************************************!



!
END MODULE FSCATAB
!
!**********************************************************************!
!**********************************************************************!
!**********************************************************************!