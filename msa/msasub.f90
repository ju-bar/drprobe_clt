!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  msasub.f90                                            !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2009-2023   !
!                                                                      !
!**********************************************************************!
!                                                                      !
!    Purpose  : Implementation of important sub-routines used by the   !
!               program MSA (see msa.f90)                              !
!                                                                      !
!    Link with: MSAparams.F90, MultiSlice.F90, emsdata.F90             !
!                                                                      !
!    Uses     : ifport                                                 !
!                                                                      !
!**********************************************************************!
!                                                                       
!  Author:  Juri Barthel                                
!           Ernst Ruska-Centre                                          
!           Forschungszentrum J�lich GmbH, 52425 J�lich, Germany        
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
!
! List of implemented functions: (2011-09-23)
!
!    SUBROUTINE CriticalError(smessage)
!    SUBROUTINE PostWarning(smessage)
!    SUBROUTINE PostMessage(smessage)
!    SUBROUTINE PostDebugMessage(smessage)
!    SUBROUTINE PostSureMessage(smessage)
!    SUBROUTINE Introduce()
!    subroutine CheckLicense()
!    FUNCTION factorial(n)
!    FUNCTION binomial(n,k)
!    FUNCTION sigmoid(x,x0,dx)
!    SUBROUTINE SetVarString(carray,n,string)
!    SUBROUTINE GetVarString(string,carray,n)
!    SUBROUTINE GetFreeLFU(lfu,lfu0,lfumax)
!    SUBROUTINE SaveDataC8(sfile,dat,n,nerr)
!    SUBROUTINE SaveDataR4(sfile,dat,n,nerr)
!    SUBROUTINE RepeatDataComplex(cin, cout, nix, nrepx, nox, niy, nrepy, noy, nerr)
!    SUBROUTINE ExplainUsage()
!    SUBROUTINE ParseCommandLine()
!    SUBROUTINE LoadParameters(sprmfile)
!    SUBROUTINE GetSliceFileName(nslc,nvar,sfname,nerr)
!    SUBROUTINE SetGlobalCellParams()
!    SUBROUTINE PrepareSupercells()
!    SUBROUTINE PrepareWavefunction()
!    SUBROUTINE DetectorReadout(rdata, ndet, nret)
! OUT    SUBROUTINE DetectorReadoutSpecial(rdata, rlost, nret)
!    SUBROUTINE ApplySpatialCoherence()
!    SUBROUTINE MSACalculate()
!    SUBROUTINE STEMMultiSlice()
!    SUBROUTINE CTEMMultiSlice()
!    SUBROUTINE ExportSTEMData(sfile)
!    SUBROUTINE ExportWave(nidx,sfile)
!    SUBROUTINE ExportWave2(nidx,sfile)
!    SUBROUTINE ExportWaveDirect(nidx,sfile)
!    SUBROUTINE SaveResult(soutfile)
!
!**********************************************************************!
!**********************************************************************!






!**********************************************************************!
!**********************************************************************!
SUBROUTINE CriticalError(smessage)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none

! ------------
! DECLARATION
  character*(*) :: smessage
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > CriticalError: INIT."
! ------------

! ------------
  MSP_err_num = MSP_err_num + 1
  write(unit=MSP_stdout,fmt='(A)') "Critical error: "//trim(smessage)
  write(unit=MSP_stdout,fmt='(A)') ""
  call MSP_HALT()
! ------------

! ------------
!  write(unit=*,fmt=*) " > CriticalError: EXIT."
  return

END SUBROUTINE CriticalError
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE PostWarning(smessage)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none

! ------------
! DECLARATION
  character*(*) :: smessage
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PostWarning: INIT."
! ------------

! ------------
  MSP_warn_num = MSP_warn_num + 1
  write(unit=MSP_stdout,fmt='(A)') ""
  write(unit=MSP_stdout,fmt='(A)') "Warning: "//trim(smessage)
  write(unit=MSP_stdout,fmt='(A)') ""
! ------------

! ------------
!  write(unit=*,fmt=*) " > PostWarning: EXIT."
  return

END SUBROUTINE PostWarning
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE PostMessage(smessage)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none

! ------------
! DECLARATION
  character*(*) :: smessage ! text message to be displayed
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PostMessage: INIT."
! ------------

! ------------
  if (DEBUG_EXPORT>=1 .or. VERBO_EXPORT>=1) then
    write(unit=MSP_stdout,fmt='(A)') trim(smessage)
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > PostMessage: EXIT."
  return

END SUBROUTINE PostMessage
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE PostSureMessage(smessage)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none

! ------------
! DECLARATION
  character*(*) :: smessage
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PostSureMessage: INIT."
! ------------

! ------------
  if (DEBUG_EXPORT>=1 .or. VERBO_EXPORT>=0) then
    write(unit=MSP_stdout,fmt='(A)') trim(smessage)
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > PostSureMessage: EXIT."
  return

END SUBROUTINE PostSureMessage
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE PostDebugMessage(smessage)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none

! ------------
! DECLARATION
  character*(*) :: smessage
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PostMessage: INIT."
! ------------

! ------------
  if (DEBUG_EXPORT>=1) then
    write(unit=MSP_stdout,fmt='(A)') trim(smessage)
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > PostMessage: EXIT."
  return

END SUBROUTINE PostDebugMessage
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE PostRuntime(sinfo, ipl)
! function: Posts a run-time message with "sinfo" as add text and
!           ipl setting the output function:
!           ipl = 0 : PostMessage
!           ipl = 1 : PostSureMessage
! -------------------------------------------------------------------- !
! parameter: sinfo : character(len=*) : add info
!            ipl : integer*4 : message function flag
! -------------------------------------------------------------------- !

  use MSAparams
  
  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sinfo
  integer*4, intent(in) :: ipl
  real*4 :: rt
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PostRuntime: INIT."
! ------------

! ------------
  if (MSP_runtimes==1 .and. VERBO_EXPORT>=0) then
    call MSP_GETCLS(rt)
    write(unit=MSP_stmp,fmt='(G15.3)') rt
    if (ipl==1) call PostSureMessage("- "//trim(sinfo)//" ("//trim(adjustl(MSP_stmp))//" s).")
    if (ipl==0) call PostMessage("- "//trim(sinfo)//" ("//trim(adjustl(MSP_stmp))//" s).")
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > PostRuntime: EXIT."
  return

END SUBROUTINE PostRuntime
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE Introduce()
! function: Posts a "start message"
! -------------------------------------------------------------------- !

  use MSAparams
  
  implicit none

! ------------
! DECLARATION
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > Introduce: INIT."
! ------------

! ------------
  call PostSureMessage("")
  call PostSureMessage(" +---------------------------------------------------+")
  !                                "[msa] MultiSlice Algorithm"
  call PostSureMessage(" | Program: "//trim(MSP_callApp)//   "               |")
  !                                "1.1.4 64-bit  -  2023 November 13  -"
  call PostSureMessage(" | Version: "//trim(MSP_verApp)//              "     |")
  !                                "Dr. J. Barthel, ju.barthel@fz-juelich.de"
  call PostSureMessage(" | Author : "//trim(MSP_authApp)//                 " |")
  call PostSureMessage(" |          Forschungszentrum Juelich GmbH, GERMANY  |")
  call PostSureMessage(" | License: GNU GPL 3 <http://www.gnu.org/licenses/> |")
  call PostSureMessage(" +---------------------------------------------------+")
  call PostSureMessage("")
! ------------

! ------------
!  write(unit=*,fmt=*) " > Introduce: EXIT."
  return

END SUBROUTINE Introduce
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION factorial(n)
! function: calculates the factorial of n -> n!
! -------------------------------------------------------------------- !
! parameter: integer*4 :: n
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: factorial
  integer*4, intent(in) :: n
  integer*4 :: i
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > factorial: INIT."
  factorial = 0 ! precheck default -> this means ERROR!
  if (n<0) return
  factorial = 1 ! postcheck default -> this means NO ERROR!
  i=2
! ------------

! ------------
  do while (n>=i)
    factorial = factorial * i
    i = i + 1
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > factorial: EXIT."
  return

END FUNCTION factorial
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION binomial(n,k)
! function: calculates the binomial coefficient of (n over k), which is
!           equal to (n!)/( (n-k)! * k! )
! -------------------------------------------------------------------- !
! parameter: integer*4 :: n,k
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: binomial
  integer*4, intent(in) :: n, k
  integer*4, external :: factorial
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > binomial: INIT."
  binomial = 0 ! precheck default -> this means ERROR!
  if (n<0.or.k<0.or.n<k) return
  binomial = factorial(n)/( factorial(n-k) * factorial(k) )
! ------------

! ------------
!  write(unit=*,fmt=*) " > binomial: EXIT."
  return

END FUNCTION binomial
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
FUNCTION sigmoid(x,x0,dx)
! function: 0.5*(tanh((x-x0)/dx)+1)
! -------------------------------------------------------------------- !
! parameter: all real*4
!            
! -------------------------------------------------------------------- !

  use precision
  
  implicit none

! ------------
! declaration
  real(fpp) :: sigmoid
  real(fpp), intent(in) :: x, x0, dx
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > sigmoid: INIT."
  sigmoid = 0.5*(tanh((x-x0)/dx)+1.0)
! ------------



! ------------
!  write(unit=*,fmt=*) " > sigmoid: EXIT."
  return

END FUNCTION sigmoid
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE SetVarString(carray,n,string)
! function: copy data from string to integer*1 array
! -------------------------------------------------------------------- !
! parameter: carray - integer*1(n) array
!            n      - array size
!            string - character(len=*)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character(len=*) :: string
  integer*1, dimension(1:n) :: carray
  integer*4 :: i, n, slen, alen, mlen
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > SetVarString: INIT."
  alen = size(carray,dim=1)
  slen = len(string)
  mlen = min(alen,slen)
  if (mlen<=0) return
! ------------

! ------------
! copy
  do i=1, mlen
    carray(i)=mod(ichar(string(i:i)),256)
  end do
  if (alen>mlen) then
    do i=1+mlen, alen
      carray(i) = 0
    end do
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > SetVarString: EXIT."
  return

END SUBROUTINE SetVarString
!**********************************************************************!






!**********************************************************************!
!**********************************************************************!
SUBROUTINE GetVarString(string,carray,n)
! function: copy n data from character array to string
! -------------------------------------------------------------------- !
! parameter: carray - integer*1(n) array
!            n      - size of carray
!            string - character(len=*)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character(len=n) :: string
  integer*1, dimension(1:n) :: carray
  integer*4 :: n, i, slen, alen, mlen
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > GetVarString: INIT."
  alen = size(carray,dim=1)
  slen = len(string)
  mlen = min(alen,slen)
  if (mlen<=0) return
! ------------

! ------------
! copy
  do i=1, mlen
    string(i:i) = achar(carray(i))
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > GetVarString: EXIT."
  return

END SUBROUTINE GetVarString
!**********************************************************************!







!**********************************************************************!
!**********************************************************************!
SUBROUTINE GetFreeLFU(lfu,lfu0,lfumax)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, intent(inout) :: lfu
  integer*4, intent(in) :: lfu0, lfumax
  logical :: isopen
  integer*4 :: u0, u1
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > GetFreeLFU: INIT."
  lfu = lfu0
  u0 = lfu0
  u1 = max(lfumax,lfu0)
! ------------

! ------------
  
  do
    INQUIRE(unit=lfu,opened=isopen)
    if (.not.isopen) exit
    lfu = lfu + 1
!   catch to many open trials
    if (lfu>u1) then
      call CriticalError("GetFreeLFU: Failed to acquire logical file unit.")
    end if
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > GetFreeLFU: EXIT."
  return

  END SUBROUTINE GetFreeLFU
!**********************************************************************!

! here are two functions copied from binio2.f90
! take care to remove when linking with binio2.f90
  


!**********************************************************************!
!
! subroutine sinsertslcidx
!
! inserts and additional suffix string with slice index suffix to
! a file name before the last dot character and sets a given file
! extension
!
! idxlen <= 0 will switch off the addition of the index suffix
!             and the "_sl" prefix.
!
subroutine sinsertslcidx(idx,idxlen,sfnin,sfnadd,sfnext,sfnout)
  implicit none
  integer*4, intent(in) :: idx ! index to append
  integer*4, intent(in) :: idxlen ! number of characters to use
  character(len=*), intent(in) :: sfnin ! input file name
  character(len=*), intent(in) :: sfnadd ! additional substring
  character(len=*), intent(in) :: sfnext ! file name extension
  character(len=*), intent(out) :: sfnout ! output file name
  integer*4 :: ipt ! position of last dot in sfnin
  character(len=12) :: snum
  snum = ""
  if (idxlen>0) then
    write(unit=snum,fmt='(A,I<idxlen>.<idxlen>)') "_sl", idx
  end if
  ipt = index(trim(sfnin),".",BACK=.TRUE.) ! remember position of the last dot in file prefix
  if (ipt>0) then
    sfnout = sfnin(1:ipt-1)//trim(sfnadd)//trim(snum)//trim(sfnext)
  else 
    sfnout = trim(sfnin)//trim(sfnadd)//trim(snum)//trim(sfnext)
  end if
  return
  end subroutine sinsertslcidx
  
!**********************************************************************!
!
! subroutine saddsuffix
!
! adds a suffix to the file name sfin before the last file extension
! (last dot) with a pre-string sadd and a number idx written with
! nlid digits (leading zeros might be added). The result is stored
! in sfout. With idln == 0 the numerical suffix will not be added.
!
subroutine saddsuffix(sfin, sadd, idx, idln, sfout)
  implicit none
  integer*4, parameter :: lnummax = 10 ! max length of numerical suffix
  integer*4, intent(in) :: idx ! index to append
  integer*4, intent(in) :: idln ! number of characters to use
  character(len=*), intent(in) :: sfin ! input file name
  character(len=*), intent(in) :: sadd ! additional substring
  character(len=*), intent(out) :: sfout ! output file name
  integer*4 :: ipt ! position of last dot in sfnin
  integer*4 :: inlen ! input string length
  integer*4 :: lin ! internal number length
  character(len=12) snum ! numerical suffix
  inlen = len_trim(sfin) ! length of input string
  snum = "" ! init empty
  lin = max(0,min(lnummax,idln)) ! limit the numerical number suffix length
  if (lin>0) then
    write(unit=snum,fmt='(I<lin>.<lin>)') idx ! write temp. num. suffix
    snum = trim(adjustl(snum)) ! left adjust and trim
  end if
  if (inlen==0) then ! handle case of empty input string
    sfout = trim(sadd)//trim(snum)
  else
    ipt = index(trim(sfin),".",BACK=.TRUE.) ! remember position of the last dot in file prefix
    if (ipt>0) then
      sfout = sfin(1:ipt-1)//trim(sadd)//trim(snum)//sfin(ipt:inlen)
    else 
      sfout = trim(sfin)//trim(sadd)//trim(snum)
    end if
  end if
  return
end subroutine saddsuffix

!**********************************************************************!
!
! subroutine getpath
!
! returns a string "spath" containing the folder part of the input 
! filepathname string "sfilepath"
!
! the length of the string variables should be set by the calling routine
!
subroutine getpath(sfilepath, spath)
  implicit none
  character(len=*), intent(in) :: sfilepath ! input
  character(len=*), intent(out) :: spath ! output
  character(len=2), parameter :: sdelim = "/\" ! path delimiting characters
  integer*4 :: ild ! last delimiter character
  ild = 0
  spath = ""
  if (len_trim(sfilepath)>1) then ! need some characters in the input
    ild = max(ild,index( trim(sfilepath), sdelim(1:1), back=.true. ))
    ild = max(ild,index( trim(sfilepath), sdelim(2:2), back=.true. ))
  end if
  if (ild>0 .and. len(spath)>=ild ) then
    !     the backwards search found a delimiter character
    ! AND the receiving string is long enough to get it
    spath = sfilepath(1:ild)
  end if
  return
end subroutine getpath


!**********************************************************************!
!
! subroutine createfilefolder
!
! creates a folder on disk to contain files for writing. The folder
! name is extracted from the input filepathname string "sfilepath"
!
subroutine createfilefolder(sfilepath,nerr)
  use ifport
  implicit none
  character(len=*), intent(in) :: sfilepath ! input
  integer*4, intent(inout) :: nerr
  character(len=8192) :: spath ! path string
  logical :: pathexists
  external :: getpath
  nerr = 0
  pathexists = .false.
  call getpath(trim(sfilepath),spath)
  if (len_trim(spath)>0) then
  
    !DEC$ IF DEFINED (__INTEL_COMPILER)
    INQUIRE ( DIRECTORY=trim(spath), EXIST=pathexists )
    !DEC$ ELSE
    INQUIRE ( FILE=trim(spath), EXIST=pathexists )
    !DEC$ ENDIF
    
    if (.not.pathexists) then  
      nerr = systemqq('mkdir "'//trim(spath)//'"')
    end if
    
  end if
  return
end subroutine createfilefolder

  
!**********************************************************************!
! opens a file for exclusive (non-shared) binary read&write access  
subroutine fileopenexclrw(sfile, lun, nexist, nerr)
  use IFPORT
  integer*4, parameter :: ntrymax = 100
  integer*4, parameter :: nwaitms = 50
  character(len=*), intent(in) :: sfile ! file name
  integer*4, intent(out) :: lun, nexist, nerr
  integer*4 :: ilun, ioerr, iopentry
  logical :: fexist, fopen
  fopen = .FALSE.
  fexist = .FALSE.
  nerr = 0
  ioerr = 0
  nexist = 0
  lun = 0
  ilun = 0
  iopentry = 0
  do while ((.not.fopen) .and. iopentry < ntrymax)
    ! check whether file already exists
    inquire(file=trim(sfile),exist=fexist)
    call GetFreeLFU(ilun,20,99) ! get a logical file unit
    if (fexist) then ! file exists, open -> load -> add -> store
      ! open an existing file (exclusive access)
      open ( unit=ilun, file=trim(sfile), iostat=ioerr, &
        &    form='BINARY', action='READWRITE', status='OLD', &
        &    share='DENYRW', position='REWIND', access='STREAM' )
    else ! file doesn't exists -> create -> store
      ! open a new file (exclusive access)
      open ( unit=ilun, file=trim(sfile), iostat=ioerr, &
        &    form='BINARY', action='READWRITE', status='NEW', &
        &    share='DENYRW', access='STREAM' )
    end if
    iopentry = iopentry + 1
    if (nerr/=0) then
      call sleepqq(nwaitms)
      cycle ! -> opening loop
    end if
    INQUIRE(unit=ilun, opened=fopen)
  end do ! while (.not.fopen)
  if ((.not.fopen) .and. iopentry>=ntrymax) then
    nerr = 1
    return
  end if
  if (fexist) nexist = 1
  if (fopen) lun = ilun
  return
end subroutine fileopenexclrw
!**********************************************************************!
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE SaveDataC(sfile,dat,n,nerr)
! function: saves complex data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use precision
  
  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  complex(fpp), intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GETFreeLFU
  external :: createfilefolder
! ------------

  nerr = 0
  call GetFreeLFU(lu,20,99)
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    call CriticalError("SaveDataC: Failed to open file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("SaveDataC: Failed to write data to file ["//trim(sfile)//"].")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  return

  END SUBROUTINE SaveDataC
!**********************************************************************!
  
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE SaveDataC8(sfile,dat,n,nerr)
! function: saves complex data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  complex*8, intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GETFreeLFU
  external :: createfilefolder
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > SaveDataC8: INIT."
  nerr = 0
  call GetFreeLFU(lu,20,99)
! ------------

! ------------
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    call CriticalError("SaveDataC8: Failed to open file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("SaveDataC8: Failed to write data to file ["//trim(sfile)//"].")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  !
  ! done
  !
! ------------

! ------------
!  write(unit=*,fmt=*) " > SaveDataC8: EXIT."
  return

END SUBROUTINE SaveDataC8
!**********************************************************************!

  
  !**********************************************************************!
!**********************************************************************!
SUBROUTINE AppendDataC(sfile,dat,n,nerr)
! function: appands complex data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use precision

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  complex(fpp), intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GETFreeLFU
  external :: createfilefolder
! ------------

  nerr = 0
  call GetFreeLFU(lu,20,99)
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu, file=trim(sfile), iostat=nerr, &
     & form='binary', action='write', status='unknown', position='append')
  if (nerr/=0) then
    call CriticalError("AppendDataC: Failed to open file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("AppendDataC: Failed to write data to file ["//trim(sfile)//"].")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  return

END SUBROUTINE AppendDataC
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE AppendDataC8(sfile,dat,n,nerr)
! function: appands complex data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  complex*8, intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GETFreeLFU
  external :: createfilefolder
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > AppendDataC8: INIT."
  nerr = 0
  call GetFreeLFU(lu,20,99)
! ------------

! ------------
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu, file=trim(sfile), iostat=nerr, &
     & form='binary', action='write', status='unknown', position='append')
  if (nerr/=0) then
    call CriticalError("AppendDataC8: Failed to open file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("AppendDataC8: Failed to write data to file ["//trim(sfile)//"].")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  !
  ! done
  !
! ------------

! ------------
!  write(unit=*,fmt=*) " > AppendDataC8: EXIT."
  return

  END SUBROUTINE AppendDataC8
!**********************************************************************!

  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE CreateDataR(sfile,n,nerr)
! function: creates a file with n float zeros
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !
  use precision
  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu, nalloc
  external :: GetFreeLFU
  external :: createfilefolder
  real(fpp), allocatable :: fa(:)
! ------------

  nerr = 0
  nalloc = 0
  allocate(fa(n),stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("CreateDataR: Failed to allocate memory.")
  end if
  fa = 0.
  call GetFreeLFU(lu,20,99)
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    call CriticalError("CreateDataR: Failed to create file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) fa
  if (nerr/=0) then
    call CriticalError("CreateDataR: Failed to initialize file.")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  !
  ! done
  !
  deallocate(fa,stat=nalloc)
  return

END SUBROUTINE CreateDataR
!**********************************************************************!
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE CreateDataR4(sfile,n,nerr)
! function: creates a file with n float zeros
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu, nalloc
  external :: GetFreeLFU
  external :: createfilefolder
  real*4, allocatable :: fa(:)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > CreateDataR4: INIT."
  nerr = 0
  nalloc = 0
  allocate(fa(n),stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("CreateDataR4: Failed to allocate memory.")
  end if
  fa = 0.
  call GetFreeLFU(lu,20,99)
! ------------

! ------------
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    call CriticalError("CreateDataR4: Failed to create file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) fa
  if (nerr/=0) then
    call CriticalError("CreateDataR4: Failed to initialize file.")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  !
  ! done
  !
  deallocate(fa,stat=nalloc)
! ------------

! ------------
!  write(unit=*,fmt=*) " > CreateDataR4: EXIT."
  return

END SUBROUTINE CreateDataR4
!**********************************************************************!
  
  
!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE AddDataR4(sfile,dat,n,np,nerr)
!! function: adss real*4 data array to file at location np
!! -------------------------------------------------------------------- !
!! parameter: 
!! -------------------------------------------------------------------- !
!
!  use IFPORT
!  
!  implicit none
!
!! ------------
!! DECLARATION
!  integer*4, parameter :: nbuf = 1024 ! use a fix size write buffer
!  
!  character(len=*), intent(in) :: sfile
!  integer*4, intent(in) :: n, np
!  real*4, intent(in) :: dat(n)
!  integer*4, intent(inout) :: nerr
!  
!  integer*4 :: lu, nalloc, npos, nposn, nrem, nnow
!  external :: GETFreeLFU
!  real*4 :: buf(nbuf)
!! ------------
!
!! ------------
!! INIT
!!  write(unit=*,fmt=*) " > AddDataR4: INIT."
!  nerr = 0
!  nalloc = 0
!  buf = 0.
!  call GetFreeLFU(lu,20,99)
!! ------------
!
!! ------------
!  !
!  ! open file, connect to lun, using shared access
!  ! Warning: This can cause trouble on multi-process calls
!  !          To avoid or at least keep chances low for simultaneous access:
!  !          - keep nbuf low 
!  !          - call only once per process and file
!  !
!  open(unit=lu, file=trim(sfile), iostat=nerr, &
!     & form='binary', action='write', status='old', share='DENYNONE')
!  if (nerr/=0) goto 100
!  !
!  !
!  nrem = n
!  npos = np
!  nposn = 1
!  !
!  ! while items remain to be written
!  do while (nrem>0)
!    !
!    ! determine number of bytes for next block
!    nnow = min(nrem,nbuf)
!    !
!    ! go to requested position in file
!    nerr = FSEEK(lu, npos, 0)
!    if (nerr/=0) goto 103
!    !
!    ! read to buffer
!    read(unit=lu,iostat=nerr) buf(1:nnow)
!    if (nerr/=0) goto 101
!    !
!    ! add from input
!    buf(1:nnow) = buf(1:nnow) + dat(nposn:nposn+nnow-1)
!    !
!    ! go back to requested position in file
!    nerr = FSEEK(lu, npos, 0)
!    if (nerr/=0) goto 103
!    !
!    ! write data to file, sequential binary
!    write(unit=lu,iostat=nerr) buf(1:nnow)
!    if (nerr/=0) goto 102
!    !
!    nrem = nrem - nnow
!    npos = npos + nnow
!    nposn = nposn + nnow
!    !
!  end do
!  !
!  ! close and disconnect
!  !
!  close(unit=lu)
!  !
!  ! done
!  !
!! ------------
!
!! ------------
!!  write(unit=*,fmt=*) " > AddDataR4: EXIT."
!  return
!  
!100 continue
!  call CriticalError("AddDataR4: Failed to open file ["//trim(sfile)//"].")
!  return
!101 continue
!  call CriticalError("AddDataR4: Failed to read data from file ["//trim(sfile)//"].")
!  return
!102 continue
!  call CriticalError("AddDataR4: Failed to write data to file ["//trim(sfile)//"].")
!  return
!103 continue
!  call CriticalError("AddDataR4: Failed to position in file ["//trim(sfile)//"].")
!  return
!
!END SUBROUTINE AddDataR4
!!**********************************************************************!
  

!**********************************************************************!
!**********************************************************************!
SUBROUTINE SaveDataR(sfile,dat,n,nerr)
! function: saves real*4 data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !
  use precision
  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real(fpp), intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GetFreeLFU
  external :: createfilefolder
! ------------

  nerr = 0
  call GetFreeLFU(lu,20,99)
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    call CriticalError("SaveDataR: Failed to open file.")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("SaveDataR: Failed to write data to file.")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  return

END SUBROUTINE SaveDataR
!**********************************************************************!
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE SaveDataR4(sfile,dat,n,nerr)
! function: saves real*4 data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real*4, intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GetFreeLFU
  external :: createfilefolder
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > SaveDataR4: INIT."
  nerr = 0
  call GetFreeLFU(lu,20,99)
! ------------

! ------------
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    call CriticalError("SaveDataR4: Failed to open file.")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("SaveDataR4: Failed to write data to file.")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  !
  ! done
  !
! ------------

! ------------
!  write(unit=*,fmt=*) " > SaveDataR4: EXIT."
  return

END SUBROUTINE SaveDataR4
!**********************************************************************!
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE AppendDataR(sfile,dat,n,nerr)
! function: appends real(fpp) data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !
  use precision
  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real(fpp), intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GETFreeLFU
  external :: createfilefolder
! ------------

  nerr = 0
  call GetFreeLFU(lu,20,99)
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu, file=trim(sfile), iostat=nerr, &
     & form='binary', action='write', status='unknown', position='append')
  if (nerr/=0) then
    call CriticalError("AppendDataR: Failed to open file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("AppendDataR: Failed to write data to file ["//trim(sfile)//"].")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  !
  ! done
  !
  return

END SUBROUTINE AppendDataR
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE AppendDataR4(sfile,dat,n,nerr)
! function: appends real*4 data array to file
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real*4, intent(in) :: dat(n)
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lu
  external :: GETFreeLFU
  external :: createfilefolder
! ------------

  nerr = 0
  call GetFreeLFU(lu,20,99)
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  open(unit=lu, file=trim(sfile), iostat=nerr, &
     & form='binary', action='write', status='unknown', position='append')
  if (nerr/=0) then
    call CriticalError("AppendDataR4: Failed to open file ["//trim(sfile)//"].")
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lu,iostat=nerr) dat
  if (nerr/=0) then
    call CriticalError("AppendDataR4: Failed to write data to file ["//trim(sfile)//"].")
  end if
  !
  ! close and disconnect
  !
  close(unit=lu)
  return

END SUBROUTINE AppendDataR4
!**********************************************************************!

!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE LoadDataR4(sfile,dat,n,nerr)
!! function: loads real*4 data array from file
!! -------------------------------------------------------------------- !
!! parameter: 
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! DECLARATION
!  character(len=*), intent(in) :: sfile
!  integer*4, intent(in) :: n
!  real*4, intent(inout) :: dat(n)
!  integer*4, intent(inout) :: nerr
!  
!  integer*4 :: lu, i
!  character(len=1000) :: stmp1, stmp2
!! ------------
!
!! ------------
!! INIT
!!  write(unit=*,fmt=*) " > LoadDataR4: INIT."
!  nerr = 0
!  call GetFreeLFU(lu,20,99)
!! ------------
!
!! ------------
!  !
!  ! open file, connect to lun
!  !
!  open(unit=lu,file=trim(sfile),iostat=nerr,&
!     & form='binary',action='read',status='old')
!  if (nerr/=0) then
!    call PostWarning("LoadDataR4: Failed to open file ["// &
!      & trim(sfile)//"].")
!    nerr = 1
!    return
!  end if
!  !
!  ! write data to file, sequential binary
!  !
!  read(unit=lu,iostat=nerr) dat
!  if (nerr/=0) then
!    i = nint(real(n*4)/1024.0)
!    write(unit=stmp1,fmt='(I)') n
!    write(unit=stmp2,fmt='(I)') i
!    call PostWarning("LoadDataR4: Failed to read "// &
!      & trim(adjustl(stmp1))//" 32-bit floats ("// &
!      & trim(adjustl(stmp2))//" kB) from file ["// &
!      & trim(sfile)//"].")
!    return
!  end if
!  !
!  ! close and disconnect
!  !
!  close(unit=lu)
!  !
!  ! done
!  !
!! ------------
!
!! ------------
!!  write(unit=*,fmt=*) " > LoadDataR4: EXIT."
!  return
!
!END SUBROUTINE LoadDataR4
!!**********************************************************************!





!!**********************************************************************!
!!**********************************************************************!
!SUBROUTINE RepeatDataComplex(cin, cout, nix, nrepx, nox, &
!     &  niy, nrepy, noy, nerr)
!! function: repeats complex data cin nrepx times along dimension 1
!!           and nrepy times along dimension 2 and writes result to
!!           array cout
!! -------------------------------------------------------------------- !
!! parameter: all treated as INOUT
!!   IN/OUTPUT:
!!     complex*8 :: cin(nix,niy)     = input data array
!!     complex*8 :: cout(nox,noy)    = output data array
!!     integer*4 :: nix, niy         = size of input array
!!     integer*4 :: nox, noy         = size of output array
!!     integer*4 :: nrepx, nrepy     = repetition rate of dimension 1 & 2
!!     integer*4 :: nerr             = error code, success = 0
!!   REMARK:
!!     MUST BE: nox = nix*nrepx AND noy = niy*nrepy
!!              cin, cout allocated
!!              none of nix, niy, nox, noy, nrepx, nrepy has value 0
!! -------------------------------------------------------------------- !
!
!  implicit none
!
!! ------------
!! DECLARATION
!  complex*8, intent(inout) :: cin(nix,niy)
!  complex*8, intent(inout) :: cout(nox,noy)
!  integer*4, intent(inout) :: nix, niy, nox, noy, nrepx, nrepy, nerr
!  
!  integer*4 :: j, k, i1, i2
!! ------------
!
!! ------------
!! INIT
!!  write(unit=*,fmt=*) " > RepeatDataComplex: INIT."
!  nerr = 0
!! ------------
!
!! ------------
!! initial checks of input parameters
!  if ((nix*nrepx/=nox).or.(niy*nrepy/=noy)) then
!    nerr = 1
!    return
!  end if
!  if ((nix*niy<=0).or.(nrepx*nrepy<=0).or.(nox*noy<=0)) then
!    nerr = 1
!    return
!  end if
!! ------------
!
!
!! ------------
!  do j = 1, niy ! loop through all rows of input data
!    ! repeat row in x
!    do k = 1, nrepx
!      i1 = 1+nix*(k-1)
!      i2 = nix*k
!      cout(i1:i2,j) = cin(1:nix,j)
!    end do
!    ! repeat row in y
!    if (nrepy>1) then
!      do k = 2, nrepy
!        i2 = j+niy*(k-1)
!        cout(1:nox,i2) = cout(1:nox,j)
!      end do
!    end if
!    ! next row
!  end do
!! ------------
!
!
!! ------------
!!  write(unit=*,fmt=*) " > RepeatDataComplex: EXIT."
!  return
!
!END SUBROUTINE RepeatDataComplex
!!**********************************************************************!

















!**********************************************************************!
!**********************************************************************!
SUBROUTINE ExplainUsage()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none

! ------------
! DECLARATION
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > <NAME>: INIT."
! ------------

! ------------
  call Introduce()
  call PostSureMessage("Usage of msa in command line:")
  call PostSureMessage("msa -prm <Parameter filename>")
  call PostSureMessage("    -out <Output filename>")
  call PostSureMessage("    [-in <Input image filename>]")
  call PostSureMessage("    [-inw <Input wave filename> <insert slice number>]")
  call PostSureMessage("    [-px <horizontal scan-pixel number>]")
  call PostSureMessage("    [-py <vertical scan-pixel number>]")
  call PostSureMessage("    [-lx <last horiz. scan pixel>]")
  call PostSureMessage("    [-ly <last vert. scan pixel]")
  call PostSureMessage("    [-pcv <probe convergence semi angle in mrad>]")
  call PostSureMessage("    [-foc <defocus value in nm>]")
  call PostSureMessage("    [-tx <x beam tilt in mrad>]")
  call PostSureMessage("    [-ty <y beam tilt in mrad>]")
  call PostSureMessage("    [-otx <x object tilt in mrad>]")
  call PostSureMessage("    [-oty <y object tilt in mrad>]")
  call PostSureMessage("    [-sr <effective source radius in nm>]")
  call PostSureMessage("    [-abf <factor> apply absorption (potentials only)]")
  call PostSureMessage("    [-buni <Biso in nm^2> apply DWF (potentials only)]")
  call PostSureMessage("    [-npass <number of QEP passes)]")
  call PostSureMessage("    [-slc <slice filename>]")
  call PostSureMessage("    [/ctem, switch to imaging TEM simulation]")
  call PostSureMessage("    [/txtout, switch for text output, STEM only]")
  call PostSureMessage("    [/3dout, switch for 3d data output, STEM only]")
  call PostSureMessage("    [/wave, save wavefunctions]")
  call PostSureMessage("    [/avwave, save average wavefunctions]")
  call PostSureMessage("    [/detimg, save images of detector functions]")
  call PostSureMessage("    [/lapro, use large-angle propagators]")
  call PostSureMessage("    [/dftest, use estimated fft plans]")
  call PostSureMessage("    [/verbose, show processing output]")
  call PostSureMessage("    [/debug, show more processing output]")
  call PostSureMessage("    [/silent, show no processing output]")
  call PostSureMessage("    [] = optional parameter")
! ------------

! ------------
!  write(unit=*,fmt=*) " > ExplainUsage: EXIT."
  return

END SUBROUTINE ExplainUsage
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE ParseCommandLine()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams
  use Plasmon

  implicit none

! ------------
! DECLARATION
  character*512 :: buffer, cmd
  logical :: fex
  integer*4 :: i, cnt, status, clen, plen, nfound, nsil, nwef, nawef
  integer*4 :: nprm, nout, nposx, nposy, nbtx, nbty
  real*4 :: mrad2deg, rtmp
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ParseCommandLine: INIT."
  i = 0
  cnt = command_argument_count()
  if (cnt==0) then
    call ExplainUsage()
    call CriticalError("Missing required program arguments.")
  end if
  nprm = 0
  nposx = 0
  nposy = 0
  nout = 0
  nbtx = 0
  nbty = 0
  nsil = 0
  mrad2deg = 0.05729578
! ------------


! ------------
! LOOP OVER ALL GIVEN ARGUMENTS
  DEBUG_EXPORT = 0
  VERBO_EXPORT = 0
  MSP_ctemmode = 0
  MSP_pimgmode = 0
  MSP_pdifmode = 0
  MSP_padifmode = 0
  MSP_ScanPixelX = -1
  MSP_ScanPixelY = -1
  MSP_LastScanPixelX = -1
  MSP_LastScanPixelY = -1
  MSP_BeamTiltX = 0.0
  MSP_BeamTiltY = 0.0
  MSP_FL_varcalc_ex = 0
  MSP_use_extalpha = 0
  MSP_extalpha = 0.0
  MSP_use_extdefocus = 0
  MSP_extdefocus = 0.0
  MSP_ext_bsx = 0.0
  MSP_ext_bsy = 0.0
  MSP_use_extbsh = 0
  MSP_use_extot = 0
  MSP_OTExX = 0.0
  MSP_OTExY = 0.0
  MSP_detimg_output = 0
  MSP_outfile = "msa_stdout.dat"
  MSP_txtout = 0
  MSP_3dout = 0
  MSP_use_extinwave = 0
  MSP_nabf = 0
  MSP_Absorption = 0.0
  MSP_nbuni = 0
  MSP_Buni = 0.01
  MSP_use_fre = 1
  MSP_SLC_lod = 0
  !MSP_FFTW_FLAG = 0
  MSP_Kmomout = 0
  MSP_KmomMmax = -1
  MSP_KmomRange = 0.0
  MSP_useldet = 0
  MSP_do_plasm = 0
  MSP_use_SLC_filenames_ex = 0
  MSP_SLC_filenames_ex = ""
  do
    i = i + 1
    if (i>cnt) exit
    
    call get_command_argument (i, buffer, clen, status)
    if (status/=0) goto 100
    
    ! CHECK COMMAND
    nfound = 0
    cmd = buffer(1:clen)
    CHECK_COMMAND: select case (cmd(1:clen))
    
    ! ???
    case ("msa")
      nfound = 1
    
    ! THE PARAMETER FILE
    case ("-prm")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_prmfile, fmt='(A)') buffer(1:plen)
      inquire(file=trim(MSP_prmfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": parameter file ["//trim(MSP_prmfile)//"] not found.")
        return
      end if
      nprm = 1
    
    ! THE OUTPUT FILE
    case ("-out")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_outfile, fmt='(A)') buffer(1:plen)
      nout = 1
    
    ! THE INPUT IMAGE FILE (FOR APPLICATION OF PARTIAL SPATIAL COHERENCE ONLY)
    case ("-in")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_infile, fmt='(A)') buffer(1:plen)
      inquire(file=trim(MSP_infile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": input image file ["//trim(MSP_infile)//"] not found.")
        return
      end if
      MSP_ApplyPSC = 1
      
    ! THE INPUT WAVE FUNCTION IMAGE FILE (e.g. FOR FURTHER CHANNELING AFTER INELASTIC TRANSITION)
    ! [REAL SPACE DATA]
    case ("-inw")
      nfound = 1
      ! read the wave file name from the next parameter
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_inwfile, fmt='(A)') buffer(1:plen)
      inquire(file=trim(MSP_inwfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": input wave function file ["//trim(MSP_inwfile)//"] not found.")
        return
      end if
      ! read the insert slice number from the next parameter
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_extinwslc
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read insert slice number.")
        return
      end if
      ! allow the use
      MSP_use_extinwave = 1
      ! set the format to real space
      MSP_extinwform = 0
    
    ! THE INPUT WAVE FUNCTION IMAGE FILE (e.g. FOR FURTHER CHANNELING AFTER INELASTIC TRANSITION)
    ! [FOURIER SPACE DATA]
    case ("-inwft")
      nfound = 1
      ! read the wave file name from the next parameter
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_inwfile, fmt='(A)') buffer(1:plen)
      inquire(file=trim(MSP_inwfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": input wave function file ["//trim(MSP_inwfile)//"] not found.")
        return
      end if
      ! read the insert slice number from the next parameter
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_extinwslc
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read insert slice number.")
        return
      end if
      ! allow the use
      MSP_use_extinwave = 1
      ! set the format to Fourier space
      MSP_extinwform = 1
      
    ! THE CURRENT HORIZONTAL SCAN POSITION (IMAGE PIXEL POS)  
    case ("-px")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_ScanPixelX
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read position index.")
        return
      end if
      nposx = 1
    
    ! THE CURRENT VERTICAL SCAN POSITION  (IMAGE PIXEL POS)
    case ("-py")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_ScanPixelY
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read position index.")
        return
      end if
      nposy = 1
      
    ! THE LAST HORIZONTAL SCAN POSITION (IMAGE PIXEL POS)  
    case ("-lx")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_LastScanPixelX
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read position index.")
        return
      end if
    
    ! THE LAST VERTICAL SCAN POSITION  (IMAGE PIXEL POS)
    case ("-ly")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_LastScanPixelY
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read position index.")
        return
      end if
      
    ! AN OPTION FOR SETTING PROBE CONVERGENCE EXTERNALLY
    case ("-pcv")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_extalpha
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read convergence angle.")
        return
      end if
      MSP_use_extalpha = 1
    
    ! AN OPTION FOR SETTING NUMBER OF QEP PASSES EXTERNALLY
    case ("-pass")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_FL_varcalc_ex
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read number of passes.")
        MSP_FL_varcalc_ex = 0
        return
      end if
      
    ! AN OPTION TO SET THE SLICE FILE NAME
    case ("-slc")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_SLC_filenames_ex, fmt='(A)') buffer(1:plen)
      ! make no checks here, this will be done when loading the parameter file
      MSP_use_SLC_filenames_ex = 1
      
    ! AN OPTION FOR SETTING A FIX DEFOCUS EXTERNALLY
    case ("-foc")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_extdefocus
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read defocus value.")
        return
      end if
      MSP_use_extdefocus = 1
      
    ! AN OPTION FOR SETTING A FIX BEAM SHIFT X EXTERNALLY
    case ("-bsx")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_ext_bsx
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to probe offset x value.")
        return
      end if
      MSP_use_extbsh = 1
      
    ! AN OPTION FOR SETTING A FIX BEAM SHIFT Y EXTERNALLY
    case ("-bsy")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_ext_bsy
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to probe offset y value.")
        return
      end if
      MSP_use_extbsh = 1
    
    ! AN OPTION FOR SETTING OBJECT TILT EXTERNALLY
    case ("-otx")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_OTExX
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read tilt value.")
        return
      end if
      MSP_OTExX = MSP_OTExX * mrad2deg
      MSP_use_extot = 1
    
    ! AN OPTION FOR SETTING OBJECT TILT EXTERNALLY
    case ("-oty")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_OTExY
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read tilt value.")
        return
      end if
      MSP_OTExY = MSP_OTExY * mrad2deg
      MSP_use_extot = 1
      
    ! AN OPTION FOR SETTING SOURCE RADIUS EXTERNALLY
    case ("-sr")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_extsrcrad
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read source radius.")
        return
      end if
      MSP_use_extsrcprm = 1
    
    ! AN OPTION FOR SETTING A BEAM TILT X EXTERNALLY
    case ("-tx")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_BeamTiltX
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read tilt value.")
        return
      end if
    
    ! AN OPTION FOR SETTING A BEAM TILT X EXTERNALLY
    case ("-ty")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_BeamTiltY
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read tilt value.")
        return
      end if
    
    ! APPLY ABSORPTIVE POTENTIALS (ONLY IN CASE OF POTENTIAL SLICE FILE INPUT)
    case ("-abf")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_Absorption
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read absorption parameter.")
        return
      end if
      MSP_nabf = 1
      
    ! APPLY UNIVERSAL DEBYE-WALLER FACTOR (ONLY IN CASE OF POTENTIAL SLICE FILE INPUT)
    case ("-buni")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_Buni
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read Biso value.")
        return
      end if
      MSP_nbuni = 1
      
    ! APPLY UNIVERSAL DEBYE-WALLER FACTOR (ONLY IN CASE OF POTENTIAL SLICE FILE INPUT)
    case ("-uuni")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_Buni
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read Uiso value.")
        return
      end if
      MSP_Buni = MSP_Buni * 78.9568352 * 0.01 ! from U [A^2] to B [nm^2]
      MSP_nbuni = 1
      
    ! CREATE A VORTEX PROBE (STEM ONLY) (added 2017-10-17 JB)
    case ("-vtx")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_Vortex
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read OAM value.")
        return
      end if
      
    ! Calculate k-moments (STEM ONLY) (added 2019-01-11 JB)
    case ("-kmom")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_KmomMmax
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read max. moment index.")
        return
      end if
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) MSP_KmomRange
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read theta range.")
        return
      end if
      MSP_Kmomout = 1
      
    ! THE DETECTION PLANE LIST FILE
    case ("-detslc")
      nfound = 1
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      write(unit = MSP_ldetfile, fmt='(A)') buffer(1:plen)
      inquire(file=trim(MSP_ldetfile),exist=fex)
      if (.not.fex) then
        if (status/=0) then
          call CriticalError("Invalid data for "//cmd(1:clen)// &
            & ": file ["//trim(MSP_ldetfile)//"] not found.")
          return
        end if
      end if
      MSP_useldet = 1
      
    ! THE PLASMON EXCITATION MODE (real plasmons)
    case ("-ple")
      nfound = 1
      ! read the plasmon energy (eV)
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) PL_ep
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read plasmon energy (eV).")
        return
      end if
      ! real the mean free path (nm)
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) PL_lp
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read mean-free path (nm).")
        return
      end if
      MSP_do_plasm = 1
      
    ! THE PLASMON EXCITATION MODE (fake low loss intraband transitions)
    case ("-plp")
      nfound = 1
      ! read the mean free path (nm)
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) PL_lp
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read mean-free path (nm).")
        return
      end if
      ! read the characteristic angle (mrad)
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) rtmp
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read characteristic angle (mrad).")
        return
      end if
      PL_qe = rtmp * 0.001 ! characteristic angle from mrad to rad
      ! read the critical angle (mrad)
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) rtmp
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read critical angle (mrad).")
        return
      end if
      PL_qc = rtmp * 0.001 ! critical angle from mrad to rad
      ! read the maximum number of excitations per electron
      i = i + 1
      if (i>cnt) goto 101
      call get_command_argument (i, buffer, plen, status)
      if (status/=0) goto 102
      read(unit=buffer,fmt=*,iostat=status) PL_npemax
      if (status/=0) then
        call CriticalError("Invalid data for "//cmd(1:clen)// &
          & ": failed to read allowed number of excitations per electron.")
        return
      end if
      MSP_do_plasm = 2

    ! ACTIVATE OUTPUT TO A TEXT LIST FILE
    case ("/txtout")
      nfound = 1
      MSP_txtout = 1
      
    ! ACTIVATE OUTPUT TO ONE 3D DATA FILE
    case ("/3dout")
      nfound = 1
      MSP_3dout = 1
      
    ! ACTIVATE DEBUG/VERBOSE OUTPUT TO CONSOLE  
    case ("/debug")
      DEBUG_EXPORT = 1
      nfound = 1
    
    ! ACTIVATE DEBUG/VERBOSE OUTPUT TO CONSOLE  
    case ("/verbose")
      VERBO_EXPORT = 1
      nfound = 1
      
    ! ACTIVATE SILENT MODE
    case ("/silent")
      nsil = 1
      nfound = 1
    
    ! ACTIVATE COHERENT TEM WAVE FUNCTION CALCULATION  
    case ("/ctem")
      nfound = 1
      MSP_ctemmode = 1
      
    case ("/pimg")
      nfound = 1
      MSP_pimgmode = 1
      
    case ("/pdif")
      nfound = 1
      MSP_pdifmode = 1
      
    case ("/padif")
      nfound = 1
      MSP_padifmode = 1
      
    case ("/wave")
      nfound = 1
      MS_wave_export = 1
      nwef = 0
      
    case ("/waveft")
      nfound = 1
      MS_wave_export = 1
      nwef = 1
      
    case ("/avwave")
      nfound = 1
      MS_wave_avg_export = 1
      nawef = 0
      
    case ("/avwaveft")
      nfound = 1
      MS_wave_avg_export = 1
      nawef = 1
    
    case ("/silavwave")
      nfound = 1
      MS_wave_avg_export = 2
      nawef = 0
      
    case ("/silavwaveft")
      nfound = 1
      MS_wave_avg_export = 2
      nawef = 1
      
    case ("/gaussap")
      nfound = 1
      STF_cap_type = 1
      
    case ("/detimg")
      nfound = 1
      MSP_detimg_output = 1
      
    case ("/fre")
      nfound = 1
      MSP_use_fre = 1
      
    case ("/lapro")
      nfound = 1
      MSP_use_fre = 0
      
    !case ("/dftest")
    !  nfound = 1
    !  MSP_FFTW_FLAG = 64
      
    case ("/rti")
      nfound = 1
      MSP_runtimes = 1
      
    case ("/epc")
      nfound = 1
      MSP_ExplicitPSC = 1
      
    case ("/slod")
      nfound = 1
      MSP_SLC_lod = 1
    
    end select CHECK_COMMAND
    
    if (nfound == 0) goto 103
  
  end do
! ------------

! ------------
! post check output options for dominant silent mode
  if (nsil==1) then ! most significant flag, turns other settings off
    VERBO_EXPORT = -1
    DEBUG_EXPORT = 0
  end if
! ------------

! ------------
  if (nprm==0) then
    call ExplainUsage()
    call CriticalError("Command line error: missing parameter file, option -prm")
  end if
  if (nout==0) then
    call PostWarning("No output file specified, using default output file name ["//trim(MSP_outfile)//"]")
  end if
  if (MSP_ctemmode==1 .and. MSP_Kmomout>0) then
    call PostWarning("K-moment analysis not supported in CTEM mode. Option -kmom is ignored.")
    MSP_Kmomout = 0
  end if
  if (MSP_Kmomout>0) then
    if (MSP_KmomMmax < 0) then
      call PostWarning("K-moment maximum order is smaller than zero. Option -kmom is ignored.")
      MSP_Kmomout = 0
    else
      MSP_KmomNum = (MSP_KmomMmax+2) * (MSP_KmomMmax+1) / 2
    end if
    if (MSP_KmomRange <= 0.0) then
      call PostWarning("K-moment integration range is smaller than zero. Option -kmom is ignored.")
      MSP_Kmomout = 0
    end if
  end if
! ------------

! ------------
! check wave export
  if (MS_wave_export > 0 .or. MS_wave_avg_export > 0 .or. &
    & MSP_pimgmode > 0 .or. MSP_pdifmode > 0 .or. MSP_padifmode > 0) then
  ! preset file names and backup names
    MS_wave_filenm = MSP_outfile
    MS_wave_filenm_bk = MS_wave_filenm
    MS_wave_filenm_avg = MS_wave_filenm
  ! determine export form from input parameters before manipulating them
    if (MS_wave_export > 0) then ! if wave export requested ...
      MS_wave_export_form = nwef ! ... this preferentially decides the export form
    else ! if no direct wave export is requested, ...
      if (MS_wave_avg_export > 0) then ! ... but avg. wave export is requested
        MS_wave_export_form = nawef ! ... the avg. wave export form is used
      end if
    end if
    MS_incwave_export = 1
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > ParseCommandLine: EXIT."
  return
  
! error handling
100 continue
  call ExplainUsage()
  call CriticalError("Command line error: missing required program arguments.")
  return
101 continue ! missing parameter behind -option
  call ExplainUsage()
  call CriticalError("Command line error: missing data for option "// &
    & cmd(1:clen))
  return
102 continue ! failed reading data behind -option
  call ExplainUsage()
  call CriticalError("Command line error: failed reading data for option "// &
    & cmd(1:clen))
  return
103 continue ! unknown option
  call ExplainUsage()
  call CriticalError("Command line error. Unknown command ["// &
    & cmd(1:clen)//"].")
  return
  

END SUBROUTINE ParseCommandLine
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE LoadParameters(sprmfile)
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams
  use MultiSlice
  
  implicit none

! ------------
! DECLARATION
  character*(*), intent(in) :: sprmfile
  logical :: ex
  integer*4 :: lfu, lfu0, lfumax, nerr, nalloc
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > LoadParameters: INIT."
  lfu0 = 20
  lfumax = 500
  nalloc = 0
! ------------

! ------------
! check file
  inquire(file=trim(sprmfile),EXIST=ex)
  if (.not.ex) then
    call CriticalError("LoadParameters: Missing file ["//trim(sprmfile)//"].")
  end if
! ------------


! ------------
! open file for reading
  call GetFreeLFU(lfu,lfu0,lfumax)
  
  open(unit=lfu,file=trim(sprmfile),iostat=nerr,action="read",share='DENYNONE')
  if (nerr/=0) then
    call CriticalError("LoadParameters: Failed to open file ["//trim(sprmfile)//"].")
  end if
  
  call PostMessage("Opened parameter file ["//trim(sprmfile)//"].")
  write(unit=MSP_stmp,fmt='(A,I4)') "- using logical unit ", lfu
  call PostDebugMessage(trim(MSP_stmp))
! ------------

! ------------
  call MSP_READparams(lfu)
! ------------

! ------------
! close file
  close(lfu)
  write(unit=MSP_stmp,fmt='(A,I4,A)') "- logical unit ", lfu," closed."
  call PostDebugMessage(trim(MSP_stmp))
! ------------


! ------------
! read detector file if switched ON
  if (MSP_usedetdef/=0) then
    ! check file
    inquire(file=trim(MSP_detfile),EXIST=ex)
    if (.not.ex) then
      call CriticalError("LoadParameters: Missing file ["//trim(MSP_detfile)//"].")
    end if
    ! get free lfu
    call GetFreeLFU(lfu,lfu0,lfumax)
  
    open(unit=lfu,file=trim(MSP_detfile),iostat=nerr,action="read",share='DENYNONE')
    if (nerr/=0) then
      call CriticalError("LoadParameters: Failed to open file ["//trim(MSP_detfile)//"].")
    end if
  
    call PostMessage("Opened detector definition file ["//trim(MSP_detfile)//"].")
    write(unit=MSP_stmp,fmt='(A,I4)') "- using logical unit ", lfu
    call PostDebugMessage(trim(MSP_stmp))
    
    call MSP_READdetdef(lfu)
    
    close(lfu)
    write(unit=MSP_stmp,fmt='(A,I4,A)') "- logical unit ", lfu," closed."
    call PostDebugMessage(trim(MSP_stmp))
  
  end if
  
  
  if (MSP_usedetdef/=1 .or. MSP_detnum<=0) then ! no special detector parameters, init standard detector arrays
  
    MSP_detnum = 1
    
    ! allocate detector parameter array
    if (allocated(MSP_detdef)) then
      deallocate(MSP_detdef,stat=nalloc)
      if (nalloc/=0) then
        call CriticalError("Failed to deallocate memory of previous detector definitions.")
      end if
    end if
    allocate(MSP_detdef(0:10,MSP_detnum), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("Failed to allocate memory for new detector definitions.")
    end if
    MSP_detdef(:,:) = 0.0_fpp
    if (allocated(MSP_detname)) then
      deallocate(MSP_detname,stat=nalloc)
      if (nalloc/=0) then
        call CriticalError("Failed to deallocate memory of previous detector names.")
      end if
    end if
    allocate(MSP_detname(MSP_detnum), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("Failed to allocate memory for new detector names.")
    end if
    MSP_detname = ""
    if (allocated(MSP_detrspfile)) then
      deallocate(MSP_detrspfile,stat=nalloc)
      if (nalloc/=0) then
        call CriticalError("Failed to deallocate memory of previous detector sensitivity profile file names.")
      end if
    end if
    allocate(MSP_detrspfile(MSP_detnum), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("Failed to allocate memory for new detector detector sensitivity profile file names.")
    end if
    MSP_detrspfile = ""
    if (allocated(MSP_detrsphdr)) then
      deallocate(MSP_detrsphdr,stat=nalloc)
      if (nalloc/=0) then
        call CriticalError("Failed to deallocate memory of previous detector sensitivity profile headers.")
      end if
    end if
    allocate(MSP_detrsphdr(MSP_detrspnhdr,MSP_detnum), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("Failed to allocate memory for new detector detector sensitivity profile headers.")
    end if
    MSP_detrsphdr = 0.0_fpp
    
    MSP_detdef(0,1) = 1.0_fpp
    MSP_detdef(1,1) = MS_detminang
    MSP_detdef(2,1) = MS_detmaxang
    MSP_detname(1) = "StdDet"
  
  end if
  
  ! allocate the detection result array
  if (allocated(MSP_detresult)) deallocate(MSP_detresult,stat=nalloc)
  if (allocated(MSP_detresult_ela)) deallocate(MSP_detresult_ela,stat=nalloc)
  allocate(MSP_detresult(MSP_detnum,0:MS_stacksize), stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("Failed to allocate detector array.")
  end if
  MSP_detresult = 0.0_fpp
  if (MS_wave_avg_export>0) then
    allocate(MSP_detresult_ela(MSP_detnum,0:MS_stacksize), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("Failed to allocate detector array for elastic channel.")
    end if
    MSP_detresult_ela = 0.0_fpp
  end if
  
  ! allocate the k-momentum result array
  if (allocated(MSP_Kmomresult)) deallocate(MSP_Kmomresult,stat=nalloc)
  allocate(MSP_Kmomresult(MSP_KmomNum,0:MS_stacksize), stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("Failed to allocate k-moment array.")
  end if
  MSP_Kmomresult = 0.0_fpp
  if (MS_wave_avg_export>0) then
    allocate(MSP_Kmomresult_ela(MSP_KmomNum,0:MS_stacksize), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("Failed to allocate k-moment array for elastic channel.")
    end if
    MSP_Kmomresult_ela = 0.0_fpp
  end if
! ------------


! ------------
!  write(unit=*,fmt=*) " > LoadParameters: EXIT."
  return

END SUBROUTINE LoadParameters
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE GetSliceFileName(nslc,nvar,sfname,nerr) !,ndslc,ndvar)
! function: creates a slice file name from gloabl parameters and indices
! -------------------------------------------------------------------- !
! parameter: 
!   INPUT:
!     integer*4 :: nslc             ! slice index
!     integer*4 :: nvar             ! slice variant
!   IN/OUTPUT:
!     character(len=*) :: sfname    ! created file name
!     integer*4 :: nerr             ! error code
! -------------------------------------------------------------------- !

  use MSAparams
  use MultiSlice

  implicit none

! ------------
! DECLARATION
  integer*4, intent(in) :: nslc, nvar
  !integer*4, intent(in), optional :: ndslc, ndvar
  integer*4, intent(inout) :: nerr
  character(len=*), intent(inout) :: sfname
  
  integer*4 :: nvd, nsd
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > GetSliceFileName: INIT."
  nerr = 0
  nsd = MSP_nslcd
  nvd = MSP_nvard
  !if (present(ndslc)) nsd = ndslc
  !if (present(ndvar)) nvd = ndvar
! ------------

! ------------
  if (MS_slicenum<=0) then
    nerr = 1
    call CriticalError("No slice definitions specified.")
  end if
  if (MSP_FL_varnum<=0) then
    nerr = 2
    call CriticalError("No slice variants defined.")
  end if
  if (nslc<1 .or. nslc>MS_slicenum) then
    nerr = 3
    write(unit=MSP_stmp,fmt=*) nslc
    call CriticalError("Invalid parameter: slice index ("// &
      & trim(adjustl(MSP_stmp))//").")
  end if
  if (nvar<1 .or. nvar>MSP_FL_varnum) then
    nerr = 4
    write(unit=MSP_stmp,fmt=*) nvar
    call CriticalError("Invalid parameter: variant index ("// &
      & trim(adjustl(MSP_stmp))//").")
  end if
! ------------

! ------------
  if (MSP_FL_varnum<=1 .or. MSP_SLI_filenamestruct==0) then
    ! not more than one variant expected OR multi-variant file structure expected
    ! generate slice file name without variant index
    write(unit=sfname,fmt='(A,I<nsd>.<nsd>,A)') &
     &  trim(MSP_SLC_filenames)//"_",nslc,".sli"
  else
    ! more than one variant expected AND NOT multi-variant file structure expected
    ! generate slice file name with variant index
    write(unit=sfname,fmt='(A,I<nvd>.<nvd>,A,I<nsd>.<nsd>,A)') &
     &  trim(MSP_SLC_filenames)//"_",nvar,"_",nslc,".sli"
  end if

! ------------
!  write(unit=*,fmt=*) " > GetSliceFileName: EXIT."
  return

END SUBROUTINE GetSliceFileName
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE SetGlobalCellParams()
! function: get global supercell data and calculation dimension
!           and global sampling
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MSAparams
  use EMSdata
  use MultiSlice

  implicit none

! ------------
! DECLARATION
  integer*4 :: i, j, i1, nslipres
  integer*4 :: nx, ny, nerr, ierr
  integer*4 :: nx1, ny1, nv, nvu, nalloc, nslpot
  integer*4, allocatable :: slipresent(:,:)
  real*4 :: szx, szy, szz, szx1, szy1, samptest, detmax, htf, ht
  character(len=MSP_ll) :: sfilename
  logical :: fexists
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > SetGlobalCellParams: INIT."
  nerr = 0
  ierr = 0
  nslpot = 0
! ------------

! ------------
  if (MS_slicenum<=0) then
    ierr = 1
    goto 99
  end if
! ------------

! ------------
  call PostMessage("Checking presence of slice files:")
  write(unit=MSP_stmp,fmt='(A,I5)') "  Expected number of slices: ", MS_slicenum
  call PostMessage(trim(MSP_stmp))
  write(unit=MSP_stmp,fmt='(A,I5)') "  Maximum number of variants per slice: ", MSP_FL_varnum
  call PostMessage(trim(MSP_stmp))
  allocate( slipresent(0:MSP_FL_varnum, MS_slicenum), stat=nalloc )
  if (nalloc/=0) then
    ierr = 2
    goto 99
  end if
  slipresent = 0
  nslipres = MS_slicenum
  ! distinguish the file structure modes
  ! 0 : multi-variant file structure
  ! 1 : single variant file structure
  if (MSP_SLI_filenamestruct==0) then ! multi-variant file structure or no variants used
    ! this should be the default way to go ! this takes longer as parameters are load from each file
    i1 = 1
    do i=1, MS_slicenum
      call GetSliceFileName(i,1,sfilename,nerr)
      inquire(file=trim(sfilename),exist=fexists)
      if (.not.fexists) cycle
      call EMS_SLI_loadparams(trim(sfilename),nx1,ny1,nv,ht,szx1,szy1,szz,nerr)
      if (nerr/=0) then
        ierr = 3
        goto 99
      end if
      if (i1==1) then ! remember parameters of the first slice
        htf = ht
        nx = nx1
        ny = ny1
        szx = szx1
        szy = szy1
        i1 = 0
      end if
      ! the number of variants used can be limited by the parameter MSP_FL_varnum
      nvu = min(nv, MSP_FL_varnum)
      !
      slipresent(0,i) = nvu ! save used variants for slice i
      slipresent(1:nvu,i) = 1 ! flag the used variants
      nslipres = nslipres - 1 ! indicate that slice data is present
      ! store parameters relevant for loading
      MSP_SLC_iprm(1,i) = EMS_SLI_data_alt_offset ! data offset byte
      MSP_SLC_iprm(2,i) = EMS_SLI_data_swap ! byte swap flag
      MSP_SLC_iprm(3,i) = nx ! number of x samples
      MSP_SLC_iprm(4,i) = ny ! number of y samples
      MSP_SLC_iprm(5,i) = slipresent(0,i) ! number of used variants
      MSP_SLC_iprm(6,i) = EMS_SLI_data_ctype ! complex data type
      ! pre-store slice thickness
      MSP_SLC_fprm(1,i) = szz
      if (EMS_SLI_data_ctype>0) nslpot = nslpot + 1
      !
    end do
  else if (MSP_SLI_filenamestruct==1) then ! single variant file structure
    ! this should usually not happen
    i1 = 1
    do i=1, MS_slicenum
      nvu = 0
      do j=1, MSP_FL_varnum
        call GetSliceFileName(i,j,sfilename,nerr)
        inquire(file=trim(sfilename),exist=fexists)
        if (.not.fexists) cycle
        if (i1==1) then ! only load parameters from the first slice
          call EMS_SLI_loadparams(trim(sfilename),nx,ny,nv,htf,szx,szy,szz,nerr)
          if (nerr/=0) then
            ierr = 4
            goto 99
          end if
          i1 = 0
        end if
        nvu = nvu + 1
        ! the number of variants used can be limited by the parameter MSP_FL_varnum
        if (nvu>MSP_FL_varnum) cycle
        slipresent(0,i) = nvu ! save used variants for slice i
        slipresent(j,i) = 1 ! flag the used variants
      end do
      if (slipresent(0,i)>0) nslipres = nslipres - 1 ! indicate that slice data is present
      ! store parameters relevant for loading
      MSP_SLC_iprm(1,i) = EMS_SLI_data_alt_offset ! data offset byte
      MSP_SLC_iprm(2,i) = EMS_SLI_data_swap ! byte swap flag
      MSP_SLC_iprm(3,i) = nx ! number of x samples
      MSP_SLC_iprm(4,i) = ny ! number of y samples
      MSP_SLC_iprm(5,i) = slipresent(0,i) ! number of used variants
      MSP_SLC_iprm(6,i) = EMS_SLI_data_ctype ! complex data type
      ! pre-store slice thickness
      MSP_SLC_fprm(1,i) = szz
      if (EMS_SLI_data_ctype>0) nslpot = nslpot + 1
    end do
  end if
  if (nslipres>0) then ! still not all slice files present,
                       ! this case is treated as critical error
    deallocate( slipresent, stat=nalloc )
    ierr = 5
    goto 99
    return
  end if
  if (nslpot>0 .and. MSP_SLC_lod>0) then
    call PostWarning("Performance loss due to load-on-demad with input potential data.")
    call PostMessage("- good performance only with phase-grating data.")
  end if
  !
  
  

  ! ------------
  ! At this point we have valid slice data. 
  ! Now we need to index the slice data in MSP_SLC_setup
  ! This is an important step.
  ! ----> The slice index in the phase grating buffer is setup in MSP_SLC_setup
  !       This array will be used by MSP_ALLOCPGR to determine the total number
  !       of phase gratings. It will also serve as index table for later access
  !       of the central phase grating buffer MSP_phasegrt.
  nslipres = 1
  MSP_SLC_num = 0
  do i=1, MS_slicenum
    MSP_SLC_setup(0,i) = slipresent(0,i) ! store the number of used variants
    ! check variants
    do j=1, MSP_FL_varnum
      if (slipresent(j,i)==0) cycle ! variant not used, cycle
      MSP_SLC_setup(j,i) = nslipres ! set hash index of the used variant
      nslipres = nslipres + 1 ! increase hash index
      MSP_SLC_num = MSP_SLC_num + 1 ! increase total variant count
    end do
  end do 
  write(unit=MSP_stmp,fmt='(A,I3)') &
     &    "   SLC setup: total number of phase gratings: ",MSP_SLC_num
  call PostDebugMessage(trim(MSP_stmp))


  !if ( nx>FFT_BOUND_MAX .or. nx<=0 .or. ny>FFT_BOUND_MAX .or. ny<=0 ) then
  !  ierr = 6
  !  goto 99
  !end if
  if ( szx<=0.0 .or. szy<=0.0 ) then
    ierr = 7
    goto 99
  end if
  ht = STF_WL2HT(STF_lamb)
  if (abs(ht-htf)>1.0) then
    call PostSureMessage("")
    write(unit=MSP_stmp,fmt='(A,G11.4,A,G11.4)') "prm: ",ht," - sli: ",htf
    call PostSureMessage("Inconsistent electron energies [keV] - "//trim(adjustl(MSP_stmp)))
    ierr = 8
    goto 99
  end if
! ------------

! ------------
! set EMS-internal data size
  EMS_SLI_data_dimx = nx
  EMS_SLI_data_dimy = ny
  write(unit=MSP_stmp,fmt=*) "loading with data dim x:",EMS_SLI_data_dimx
  call PostMessage(trim(MSP_stmp))
  write(unit=MSP_stmp,fmt=*) "loading with data dim y:",EMS_SLI_data_dimy
  call PostMessage(trim(MSP_stmp))
  
! set supercell data size
  MSP_dimcellx = nx
  MSP_dimcelly = ny
  write(unit=MSP_stmp,fmt=*) "using global cell data dim x:",MSP_dimcellx
  call PostMessage(trim(MSP_stmp))
  write(unit=MSP_stmp,fmt=*) "using global cell data dim y:",MSP_dimcelly
  call PostMessage(trim(MSP_stmp))
  
! set MS-internal data size
  MS_dimx = nx*MSP_SC_repeatx
  MS_dimy = ny*MSP_SC_repeaty
  write(unit=MSP_stmp,fmt=*) "using global wave data dim x:",MS_dimx
  call PostMessage(trim(MSP_stmp))
  write(unit=MSP_stmp,fmt=*) "using global wave data dim y:",MS_dimy
  call PostMessage(trim(MSP_stmp))
  
! set MS-internal data sampling
  MS_samplingx = szx/real(EMS_SLI_data_dimx)
  MS_samplingy = szy/real(EMS_SLI_data_dimy)
  write(unit=MSP_stmp,fmt=*) "using global cell sampling rate x [nm/pix]:",MS_samplingx
  call PostMessage(trim(MSP_stmp))
  write(unit=MSP_stmp,fmt=*) "using global cell sampling rate y [nm/pix]:",MS_samplingy
  call PostMessage(trim(MSP_stmp))
  ! check detector maximum radius against simulation maximum spatial frequency
  detmax = MS_detmaxang
!  if (MSP_usedetdef==1 .and. MSP_detnum>0) then
  detmax = 0.0
  do i=1, MSP_detnum
    detmax = max(MSP_detdef(10,i),detmax)
  end do
!  end if
  if (MSP_ctemmode==0) then
    samptest = 0.5/MS_samplingx*MS_lamb*0.667*1000.0
    if (samptest<detmax) then
      call PostWarning("Maximum detection angle is larger than 2/3-gmax multislice aperture (X).")
      write(unit=MSP_stmp,fmt=*) "  max. detection angle [mrad]:",detmax
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt=*) "  2/3-gmax multislice aperture (X) [mrad]:",samptest
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt=*) "  Use a finer horizontal supercell sampling or smaller detection angles!"
      call PostMessage(trim(MSP_stmp))
    end if
    samptest = 0.5/MS_samplingy*MS_lamb*0.667*1000.0
    if (samptest<detmax) then
      call PostWarning("Maximum detection angle is larger than 2/3-gmax multislice aperture (Y).")
      write(unit=MSP_stmp,fmt=*) "  max. detection angle [mrad]:",detmax
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt=*) "  2/3-gmax multislice aperture (Y) [mrad]:",samptest
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt=*) "  Use finer vertical supercell sampling or smaller detector semi angles!"
      call PostMessage(trim(MSP_stmp))
    end if
  end if
! ------------

! ------------
! exit part
99 if (allocated(slipresent)) deallocate( slipresent, stat=nalloc )
  select case (ierr) ! error handling
  case (0)
    return
  case (1)
    call CriticalError("No slice definitions specified.")
  case (2)
    call CriticalError("Memory allocation failed.")
  case (3)
    call CriticalError("Failed to read parameters from slice file.")
  case (4)
    call CriticalError("Failed to read size parameters from first slice file.")
  case (5)
    call CriticalError("Required slice data is not completely present.")
  case (6)
    write(unit=MSP_stmp,fmt='(I5," x ",I5)') nx, ny
    call CriticalError("Invalid slice file sampling ("// &
       & trim(adjustl(MSP_stmp))//").")
  case (7)
    write(unit=MSP_stmp,fmt='(G11.4," x ",G11.4)') szx, szy
    call CriticalError("Invalid slice physical frame size ("// &
       & trim(adjustl(MSP_stmp))//") nm^2.")
  case (8)
    call CriticalError("Electron energy input parameter and slice files are inconsistent.")
  end select ! case(nerr)
!  write(unit=*,fmt=*) " > SetGlobalCellParams: EXIT."
  return

END SUBROUTINE SetGlobalCellParams
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE PrepareSupercells()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! REMARKS:
! - requires a previous call of MSP_PREALLOCPGR, SetGlobalCellParams,
!   and MSP_ALLOCPGR
! - requires MSP_SLC_setup to be set up completely
! -------------------------------------------------------------------- !

  use MSAparams
  use EMSdata
  use MultiSlice

  implicit none

! ------------
! DECLARATION
  integer*4 :: i, j, nerr, nv, ni1, ni2, nx, ny
  real*4 :: szz
  real(fpp) :: ht, dz, fabs, buni
  character(len=MSP_ll) :: sfilename
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PrepareSupercells: INIT."
  if (.not.allocated(MSP_phasegrt)) then
    call CriticalError("Slice data memory not allocated.")
  end if
! ------------

! ------------
  if (MS_slicenum<=0) then
    call CriticalError("No slice definitions specified.")
  end if
  if (MS_stacksize<=0) then
    call CriticalError("No object slices specified.")
  end if
  if (EMS_SLI_data_dimx*EMS_SLI_data_dimy<=0) then
    call CriticalError("Invalid data size.")
  end if
  if (MSP_FL_varnum<1) then
    call CriticalError("Invalid number of slice variants.")
  end if
! ------------


! ------------
  MS_slicestack(1:MS_stacksize) = MSP_SLC_object(1:MS_stacksize)
  call PostMessage("Object slice sequence set.")
! ------------


! ------------
! load data
  if (MSP_SLC_lod==1) then
    
    call PostMessage("Load-on-demand: slice data loading postponed.")
    ! set slice thickness in MS_Multislice
    do i=1, MS_slicenum ! loop through slices
      MS_slicethick(i) = MSP_SLC_fprm(1,i)
    end do
    goto 50 ! proceed to propagator loading
    
  end if
  call PostMessage("Start loading slice data.")
  nx = MSP_dimcellx
  ny = MSP_dimcelly
  do i=1, MS_slicenum ! loop through slices
  
    if (0==MSP_SLI_filenamestruct) then
      
      ! generate slice file name
      call GetSliceFileName(i,1,sfilename,nerr)
      if (nerr/=0) call CriticalError("Failed to generate slice file name from parameters.")
      ! this loading version expects also that the variants of one slice appear
      ! in sequence in MSP_phasegrt
      nv = MSP_SLC_setup(0,i) ! number of variants present for this slice
      ni1 = MSP_GetPGRIndex(1,i,nerr) ! get index of the first variant
      ni2 = MSP_GetPGRIndex(nv,i,nerr) ! get index of the last variant
      
      call EMS_SLI_loaddata(trim(sfilename), nx, ny, nv, &
     &       MSP_phasegrt(1:nx, 1:ny, ni1:ni2), szz, nerr)
      if (nerr/=0) then
        call CriticalError("Failed to load from EMS SLI file ["//trim(sfilename)//"].")
      end if
      
      ! file info
      write(unit=MSP_stmp,fmt='(A)') "   file : "//trim(sfilename)
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt='(A)') "   title: "//trim(EMS_SLI_data_title)
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt='(A,G11.4)') "   slice thickness [nm]: ",szz
      call PostMessage(trim(MSP_stmp))
      write(unit=MSP_stmp,fmt='(A,G11.4)') "   slice electron energy [keV]: ",EMS_SLI_data_ht
      call PostMessage(trim(MSP_stmp))
      !
      MSP_SLC_title(i) = EMS_SLI_data_title
      !
      write(unit=MSP_stmp,fmt='(A,I4,A,I4,A,I4,A)') &
     &      "   data : ",nv," variants loaded to indices ",ni1," .. ",ni2,"."
      call PostDebugMessage(trim(MSP_stmp))
      
    else ! (0/=MSP_SLI_filenamestruct)
      
      ! this loading version expects only one variant per file
      nv = MSP_SLC_setup(0,i) ! number of variants present for this slice
      
      do j=1, MSP_FL_varnum ! loop through all possible slice variants  
        ni1 = MSP_GetPGRIndex(j,i,nerr) ! get index of the variant, also indicates file presence, which is not checked again here
        if (ni1<=0) cycle ! skip the variant index in case of absent data
        ! generate slice file name
        call GetSliceFileName(i,j,sfilename,nerr)
        if (nerr/=0) call CriticalError("Failed to generate slice file name from parameters.")
        
        call EMS_SLI_loaddata(trim(sfilename), nx, ny, 1, &
     &       MSP_phasegrt(1:nx ,1:ny ,ni1:ni1), szz, nerr)
        if (nerr/=0) then
          call CriticalError("Failed to load from EMS SLI file ["//trim(sfilename)//"].")
        end if
        ! file info
        write(unit=MSP_stmp,fmt='(A)') "   file : "//trim(sfilename)
        call PostMessage(trim(MSP_stmp))
        write(unit=MSP_stmp,fmt='(A)') "   title: "//trim(EMS_SLI_data_title)
        call PostMessage(trim(MSP_stmp))
        write(unit=MSP_stmp,fmt='(A,G11.4)') "   slice thickness [nm]:",szz
        call PostMessage(trim(MSP_stmp))
        write(unit=MSP_stmp,fmt='(A,G11.4)') "   slice electron energy [keV]: ",EMS_SLI_data_ht
        call PostMessage(trim(MSP_stmp))
        
        write(unit=MSP_stmp,fmt='(A,I4,A,I4,A)') &
     &      "   data : variant #",j," loaded to index ",ni1,"."
        call PostDebugMessage(trim(MSP_stmp))
      
        if (j==1) then 
          MSP_SLC_title(i) = EMS_SLI_data_title
        end if
        
      end do ! variant loop
    
    end if
    
    call PostMessage("Preparing slice data for usage.")
    nerr = MS_err_num
    ! set slice thickness in MS_Multislice
    MS_slicethick(i) = MSP_SLC_fprm(1,i)
    !
    if (EMS_SLI_data_ctype==1) then ! loaded potential data
      call PostMessage("Transforming loaded potential data to phase gratings")
      nv = MSP_SLC_setup(0,i) ! number of variants present for this slice
      ni1 = MSP_GetPGRIndex(1,i,nerr) ! get index of the first variant
      ni2 = MSP_GetPGRIndex(nv,i,nerr) ! get index of the last variant
      ht = real(EMS_SLI_data_ht, kind=fpp)
      dz = real(MSP_SLC_fprm(1,i), kind=fpp)
      fabs = real(MSP_Absorption, kind=fpp)
      buni = real(MSP_Buni, kind=fpp)
      call MS_SlicePot2Pgr( ht, dz, MSP_nabf, fabs, &
                          & MSP_nbuni, buni, &
                          & nx, ny, ni2-ni1+1, &
                          & MSP_phasegrt(1:nx ,1:ny ,ni1:ni2), nerr )
      if (nerr/=0) then
        call CriticalError("Failed to load phase gratings from EMS SLI file ["//trim(sfilename)//"].")
      end if
    end if
  
  end do ! slice loop (i)
  call PostMessage("Finished loading slice data.")
! ------------

! ------------
! Prepare propagators
50 continue ! jump label for skipping slice loading
  call PostMessage("Preparing propagators.")
  if (MSP_use_fre==1) then
    call PostMessage("Preparing Fresnel propagators.")
    call MS_PreparePropagators() ! this creates Fresnel propagators (paraxial approximation)
  else
    call PostMessage("Preparing large-angle propagators.")
    call MS_PreparePropagators2() ! this creates large-angle propagators (no paraxial approximation)
  end if
! ------------

! ------------
  call PostMessage("Supercell data prepared.")
  !  write(unit=*,fmt=*) " > PrepareSupercells: EXIT."
  return

END SUBROUTINE PrepareSupercells
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE PrepareWavefunction()
! function: 
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MultiSlice
  use STEMfunctions
  use MSAparams

  implicit none

! ------------
! DECLARATION
  integer*4 :: nerr, nalloc
  complex(fpp), allocatable :: cdata(:,:)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > PrepareWavefunction: INIT."
  allocate(cdata(MS_dimx,MS_dimy),stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("Failed to allocate memory for incident wavefunction.")
  end if
  cdata = cmplx(0.0,0.0)
! ------------


! --------------
  nerr = STF_err_num
  STF_beam_tiltx = MSP_BeamTiltX
  STF_beam_tilty = MSP_BeamTiltY
  if (MSP_ctemmode==1) then ! PLANE WAVE FOR TEM
    call STF_PreparePlaneWaveFourier(cdata,MS_dimx,MS_dimy,MS_samplingx,MS_samplingy)
  else ! STEM
    if (MSP_Vortex==0) then ! DEFAULT STEM PROBE
      call STF_PrepareProbeWaveFourier(cdata,MS_dimx,MS_dimy,MS_samplingx,MS_samplingy)
    else ! STEM VORTEX PROBE (added 2017-10-17 by JB)
      call STF_PrepareVortexProbeWaveFourier(cdata,MSP_Vortex,MS_dimx,MS_dimy,MS_samplingx,MS_samplingy)
    end if
  end if
  if (nerr/=STF_err_num) then
    if (allocated(cdata)) deallocate(cdata, stat=nalloc)
    call CriticalError("Failed to create wavefunction.")
  end if
! --------------

! ------------
! transfer to multislice module
  nerr = MS_err_num
  call MS_SetIncomingWave(cdata)
  if (nerr/=MS_err_num) then
    if (allocated(cdata)) deallocate(cdata, stat=nalloc)
    call CriticalError("Failed to save wavefunction in multislice module.")
  end if
! ------------

! ------------
  if (allocated(cdata)) deallocate(cdata, stat=nalloc)
!  write(unit=*,fmt=*) " > PrepareWavefunction: EXIT."
  return

END SUBROUTINE PrepareWavefunction
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE InsertExternalWavefunction()
! function: 
!   Loads wave function data (real-space or Fourier space) from file
!   and transfers the data to the module memory
!
!   This function was inserted with version 0.60b - 14.06.2012 (JB)
!   * modified with version 0.78b - 09.11.2017 (JB)
!     now supporting the insertion of Fourier space wave functions
!
!   - loads the input wave function data from file
!   - Fourier-Transforms the data (depends on MSP_extinwform)
!   - copies the FS data to the MS-module backup wave
!   - Fourier-space data is expected to be in scrambled form
!
!   Use MS_OffsetIncomingWave(0.0,0.0,0.0) ONCE before running the
!   modified multislice
!
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none
  
! ------------
! DECLARATION
  integer*4 :: nerr, nx, ny, nlfu, nalloc
  complex(fpp), allocatable :: cdata(:,:)
  real(fpp) :: rsca
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > InsertExternalWavefunction: INIT."
  if (MS_status<1) goto 101
  nx = MS_dimx
  ny = MS_dimy
  if (nx <= 0 .or. ny <= 0) goto 106
! allocate memory for loading the wave function
  allocate( cdata(nx,ny), stat=nalloc)
  if (nalloc/=0) goto 102
  cdata = cmplx(0.0,0.0,kind=fpp)
  rsca = 1.0 / sqrt(real(nx*ny,kind=fpp))
! ------------


! ------------
! Load the wave data from file
  call PostMessage("- Loading input wave function from file ["//trim(MSP_inwfile)//"].")
  call GetFreeLFU(nlfu,20,100)
! open the input file
  open(unit=nlfu, file=trim(MSP_inwfile), form="binary", access="sequential", &
     & iostat=nerr, status="old", action="read", share='DENYNONE' )
  if (nerr/=0) goto 103
! load wave function
  read(unit=nlfu,iostat=nerr) cdata
! handle loading error
  if (nerr/=0) goto 104
! close the input file  
  close(unit=nlfu)
  call PostMessage("- Successful read of insert wave function.")
! ------------

! ------------
! Apply FT in case of a real-space input
  if (MSP_extinwform==0) then
    MS_work(1:nx,1:ny) = cdata(1:nx,1:ny)
  ! transform the data to Fourier space
    call MS_FFT_WORK(1)
  ! transfer the loaded data to the result array
    cdata(1:nx,1:ny) =  MS_work(1:nx,1:ny) * rsca
  end if
! ------------

! ------------
! transfer the wave data to multislice module (backup)
  nerr = MS_err_num
  call MS_SetIncomingWave(cdata)
  if (nerr/=MS_err_num) goto 105
! ------------

! ------------
!  write(unit=*,fmt=*) " > InsertExternalWavefunction: EXIT."
  goto 1000

! error handling
101 continue
  call CriticalError("Multislice module is not initialized.")
  goto 1000
102 continue
  call CriticalError("Allocation of memory for input wave failed.")
  goto 1000
103 continue
  write(unit=MSP_stmp,fmt='(I)') nerr
  call CriticalError("Failed to connect to file ["// &
     & trim(MSP_inwfile)//"] - Code ("//trim(adjustl(MSP_stmp))//").")
  goto 1000
104 continue
  write(unit=MSP_stmp2,fmt='(I)') nerr
  write(unit=MSP_stmp,fmt='(A,I4,A,I4,A)') &
     &    "Failed to read all ",nx," x ",ny, &
     &    " 64-bit complex*8 data values - Code("//trim(adjustl(MSP_stmp2))//")."
  call CriticalError(trim(MSP_stmp))
  goto 1000
105 continue
  call CriticalError("Failed to save wavefunction in multislice module.")
  goto 1000
106 continue
  write(unit=MSP_stmp,fmt='(A,I4,A,I4,A)') &
     &    "Invalid size of the wavefunction array (",nx," x ",ny,")."
  call CriticalError(trim(MSP_stmp))
  goto 1000
  
! final exit
1000 continue
  if (allocated(cdata)) deallocate(cdata,stat=nalloc) 
  return

END SUBROUTINE InsertExternalWavefunction
!**********************************************************************!

  
!**********************************************************************!
! straight summation on single precision accumulator
!**********************************************************************!
SUBROUTINE FSTRSUM(a, n, s)
  IMPLICIT NONE
  ! interface
  real*4, intent(inout) :: a(n) ! reference to data array 
  integer*4, intent(in) :: n ! number of items (32 bit should be sufficient)
  real*4, intent(out) :: s ! result of summation
  ! locals
  integer*4 :: i
  ! init
  s = 0.0
  ! summation
  do i=1, n
    s = s + a(i)
  end do
  return
END SUBROUTINE FSTRSUM
  
!**********************************************************************!
! straight summation on double precision accumulator
!**********************************************************************!
SUBROUTINE DSTRSUM(a, n, s)
  IMPLICIT NONE
  ! interface
  real*4, intent(inout) :: a(n) ! reference to data array 
  integer*4, intent(in) :: n ! number of items (32 bit should be sufficient)
  real*4, intent(out) :: s ! result of summation
  ! locals
  integer*4 :: i
  real*8 :: stmp
  ! init
  stmp = 0.0
  ! summation with type cast
  do i=1, n
    stmp = stmp + dble(a(i))
  end do
  ! copy result with type cast
  s = real(stmp, kind=4)
  return
END SUBROUTINE DSTRSUM
  
!**********************************************************************!
! summation on program precision i/o with double precision accumulator
!**********************************************************************!
SUBROUTINE FPPSTRSUM(a, n, s)
  use precision
  IMPLICIT NONE
  ! interface
  real(fpp), intent(inout) :: a(n) ! reference to data array 
  integer*4, intent(in) :: n ! number of items (32 bit should be sufficient)
  real(fpp), intent(out) :: s ! result of summation
  ! locals
  integer*4 :: i
  real*8 :: stmp
  ! init
  stmp = 0.0
  ! summation with type cast
  do i=1, n
    stmp = stmp + real(a(i),kind=8)
  end do
  ! copy result with type cast
  s = real(stmp, kind=fpp)
  return
END SUBROUTINE FPPSTRSUM
  

!**********************************************************************!
! strided 2-fold butterfly summation
! ! Warning: This routine may modify the input a.
! ! Make a backup of the data or provide a copy if you still need it!
!**********************************************************************!
SUBROUTINE FDNCS2M(a, n, s)
  IMPLICIT NONE
  ! parameters
  integer*4, parameter :: nbuf = 4096
  integer*4, parameter :: nlim = 32
  ! interface
  real*4, intent(inout) :: a(n) ! reference to data array 
  integer*4, intent(in) :: n ! number of items (32 bit should be sufficient)
  real*4, intent(out) :: s ! result of summation
  ! workers
  integer*4 :: itmp, idx, idxn0, nalloc ! indices and helpers
  integer*4 :: nc, n2, n1, n0 ! stride sizes
  integer*4 :: ntmp ! number of strides through buffer
  real*4 :: r ! rest value
  real*4, allocatable :: dst(:) ! allocatble temp data
  ! init
	s = 0.0
	if (n <= 0) return ! handle invalid n
  if (n == 1) then
    s = a(1)
    return
  end if
  if (n == 2) then
    s = a(1) + a(2)
    return
  end if
  if (n < nlim) then ! small array size, do simple sim
    call FSTRSUM(a,n,s)
    return
  end if
  ! - calculate number of strides through the buffer
	ntmp = int(ceiling((dble(n) / dble(nbuf))))
	if (ntmp > 1) then ! there will be more than one stride -> array larger than buffer length
		allocate(dst(ntmp), stat=nalloc) ! allocate new destination buffer
    dst = 0.0 ! preset dst with zeroes
    if (0/=nalloc) return
  else ! array smaller than or equal to the stride buffer, butterfly on size n
		r = 0.0
    n2 = 2
    n1 = 1
    do while (n2 <= n) ! butterfly on size n
			do idx=n2, n, n2
				a(idx) = a(idx) + a(idx - n1)
			end do
			if (n1 <= modulo(n,n2)) then ! handle left-over component
				r = r + a(idx - n1)
			end if
			n1 = n2
			n2 = n2 + n2
		end do
		s = a(n1) + r
		return ! done
  end if
  n0 = 0
  itmp = 1
	do while (n0 < n) ! loop over strides (n > nbuf)
		nc = min(nbuf, n - n0) ! number of items to take for in this stride
		if (nc == nbuf) then ! working on full buffer length. This is repeated ntmp-1 times
			n2 = 2 ! init
      n1 = 1
			do while (n2 <= nbuf) ! butterfly on size nbuf
        do idx=n2, nbuf, n2
          idxn0 = n0+idx
					a(idxn0) = a(idxn0) + a(idxn0 - n1)
				end do
				n1 = n2
				n2 = n2 + n2
			end do
			dst(itmp) = a(n0 + n1) ! store intermediate result in stride slot of destination
    else ! working on reduced buffer length (not power of two), this happens only once!
			if (nc > 0) then
        if (nc == 1) then ! roll-out for 1
          dst(itmp) = a(n0+1)
        else if (nc == 2) then ! roll-out for 2
          dst(itmp) = a(n0+1) + a(n0+2)
        else if (nc < nlim) then ! small array case
          call FSTRSUM(a(n0+1:n0+nc),nc,s) ! straight sum on size nc
        else ! larger array case (nc)
          r = 0.0
          n2 = 2
          n1 = 1
          do while (n2 <= nc) ! butterfly on size nc
			      do idx=n2, nc, n2
              idxn0 = n0+idx
				      a(idxn0) = a(idxn0) + a(idxn0 - n1)
			      end do
			      if (n1 <= modulo(nc,n2)) then ! handle left-over component
				      r = r + a(n0 + idx - n1)
			      end if
			      n1 = n2
			      n2 = n2 + n2
		      end do
		      dst(itmp) = dst(itmp) + a(n0 + n1) + r
        end if
      end if
    end if
		n0 = n0 + nbuf
		itmp = itmp + 1
	end do
	if (ntmp > 1) then ! recurse if more than one buffer stride happened, this should always happen here
		call FDNCS2M(dst, ntmp, s) ! sum on dst buffer
		deallocate(dst, stat=nalloc) ! release dst buffer memory
	end if
	return
END SUBROUTINE FDNCS2M
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE DetectorReadout(rdata, ndat, nret)
! SUBROUTINE: Performs readout of all detectors and k-moment
!             integration.
!             The output array rdata is a list of integrated detector
!             values. There is no integration done here. The values
!             are set and not summed to rdata.
! -------------------------------------------------------------------- !
! parameter:  rdata (real*4) expect one item per detector and k-moment
!             component, output
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams
  
  implicit none

! ------------
! DECLARATION
  real(fpp), intent(out) :: rdata(1:ndat)
  integer*4, intent(in) :: ndat
  integer*4, intent(out) :: nret
  integer*4 :: ndet, nkmom, nx, ny, i, j, k, l, m, nmlen, idy, idx, nalloc
  real(fpp) :: rval, gxk, gyl
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > DetectorReadout: INIT."
  nret = 1
  ndet = MSP_detnum
  nkmom = MSP_KmomNum
  nalloc = 0
  nx = MS_dimx
  ny = MS_dimy
  if (0==MS_status) then
    call CriticalError("DetectorReadout: module not initialized.")
    return ! skip if data array not ready
  end if
  if (0>=nx .or. 0>=ny) then ! invalid size of diffraction grid
    call CriticalError("DetectorReadout: invalid grid size.")
    return
  end if
  if (0>=ndet .and. 0>=nkmom) then ! no detectors and no other integrating readout
    nret = 0
    return ! nothing to readout, leave rdata as is
  end if
  if (ndat < ndet+nkmom) then
    call CriticalError("DetectorReadout: conflict between number of "// &
       & "detectors and size of provided output array.")
    return
  end if
  if (.not.allocated(MSP_pdiftmp)) then ! allocate temp readout buffer
    allocate(MSP_pdiftmp(MS_dimx*MS_dimy), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("DetectorReadout: allocation failed.")
      return ! skip if data array not ready
    end if
  end if
  if (.not.allocated(MSP_pdettmp)) then ! allocate temp readout buffer
    allocate(MSP_pdettmp(MS_dimx*MS_dimy), stat=nalloc)
    if (nalloc/=0) then
      call CriticalError("DetectorReadout: allocation failed.")
      return ! skip if data array not ready
    end if
  end if
! ------------
  
  
! ------------
! calculate the diffraction data stream used by all detectors and integrators below
  do j=1, ny
    idy = (j-1)*nx
    do i=1, nx
      idx = i + idy
      MSP_pdiftmp(idx) = real(MS_wave(i,j)*conjg(MS_wave(i,j)),kind=fpp) ! get power of wave function as stream of data
    end do
  end do
! reset output array
  rdata(1:ndat) = 0.0_fpp
! default STEM detectors
  if (ndet>0) then
    ! readout all integrating detectors
    do k=1, ndet 
      nmlen = MSP_detmasklen(k)
      if ( 0 == nint(MSP_detdef(0,k)) .or. 0 >= nmlen) cycle ! skip empty or undefined detector
      do i=1, nmlen ! readout each masked Fourier pixel
        idx = MSP_detmask(i,k) ! get stream index
        MSP_pdettmp(i) = MSP_pdiftmp(idx) * MSP_detarea(idx,k) ! store detected intensity
      end do
      !call FDNCS2M(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a 2-fold butterfly
      call FPPSTRSUM(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a double precision accumulator
      rdata(k) = rval 
    end do
  end if
! k-space moments
  if (nkmom>0) then
    ! readout k-moments
    idy = 0
    nmlen = MSP_Kmommasklen ! assuming everything is OK with MSP_Kmommasklen
    ! accumulate sums
    do m=0, MSP_KmomMmax ! loop from 0 to maximum order of moments
      do l=0, m ! loop through all 2d components of moment k
        idy = idy + 1
        k = m - l ! second parameter power
        do i=1, nmlen ! readout each masked Fourier pixel
          idx = MSP_Kmommask(i) ! get stream index
          gxk = MSP_Kmomgx(MSP_kmomhash(1,i),k) ! get gx^k
          gyl = MSP_Kmomgy(MSP_kmomhash(2,i),l) ! get gy^l
          MSP_pdettmp(i) = gxk * gyl * MSP_pdiftmp(idx) * MSP_Kmomwgt(idx) ! store detected intensity
        end do
        !call FDNCS2M(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a 2-fold butterfly
        call FPPSTRSUM(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a double precision accumulator
        rdata(ndet+idy) = rval ! store result in the channels behind the default detectors
        if (m>0) then ! normalize by the 0-th moment rdata(ndet+1)
          rdata(ndet+idy) = rdata(ndet+idy) / rdata(ndet+1)
        end if
      end do
    end do
  end if
! set success code
  nret = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > DetectorReadout: EXIT."
  return

  END SUBROUTINE DetectorReadout
!**********************************************************************!
  
  
  
  
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE DetectorReadoutElastic(nerr)
! SUBROUTINE: Performs readout of all detectors and k-moment
!             integration for the elastic channel.
!             The output is written to the module arrays
!             MSP_detresult_ela and MSP_Kmomresult_ela, which are
!             only allocated if MS_wave_avg_export > 0. This is the
!             condition by which this routine should be triggered.
!             A respective if statement is performed at the beginning.
!             The readout is done here for all slices, assuming that
!             the elastic wave functions exist in MS_wave_avg.
! -------------------------------------------------------------------- !
! parameter:  nerr, integer*4 = error code
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams
  
  implicit none

! ------------
! DECLARATION
  integer*4, intent(out) :: nerr
  integer*4 :: ncrit, ndet, nkmom, nx, ny, nalloc, nmlen, npln
  integer*4 :: i, j, k, l, m, ipln, islc, idy, idx
  real(fpp) :: rval, rsca, renorm, gxk, gyl
  complex(fpp), allocatable :: wave(:,:)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > DetectorReadoutElastic: INIT."
  nerr = 0
  ncrit = 0
  nalloc = 0
  nx = MS_dimx
  ny = MS_dimy
  if (MS_wave_avg_export <= 0) goto 100 ! just quit, no output of elastic channel requested
  ! We assume now, that MSP_detresult_ela is allocated
  !    and that MSP_Kmomresult_ela is allocated
  !    if the "-kmom" option is used and MSP_Kmomresult is allocated.
  ndet = MSP_detnum
  nkmom = MSP_KmomNum
  if (0>=ndet .and. 0>=nkmom) goto 100 ! no detectors and no other integrating readout
  if (0>=MSP_detpln) goto 100 ! no planes with readout
  if (0==MS_status) then
    nerr = 1
    goto 200
  end if
  if (0>=nx .or. 0>=ny) then ! invalid size of diffraction grid
    nerr = 2
    goto 200
  end if
  if (.not.allocated(MSP_pdiftmp)) then ! allocate temp readout buffer
    allocate(MSP_pdiftmp(nx*ny), stat=nalloc)
    if (nalloc/=0) then
      nerr = 3
      goto 200
    end if
  end if
  if (.not.allocated(MSP_pdettmp)) then ! allocate temp readout buffer
    allocate(MSP_pdettmp(nx*ny), stat=nalloc)
    if (nalloc/=0) then
      nerr = 3
      goto 200
    end if
  end if
  ! determine total number of planes
  npln = MS_wave_avg_num
  ! allocate memory for a temporary wave function
  allocate(wave(nx,ny), stat=nalloc)
  if (nalloc/=0) then
    nerr = 3
    goto 200
  end if
  wave = cmplx(0.,0.,kind=fpp)
  rsca = 1.0_fpp / real(nx*ny,kind=fpp) ! FFT transform rescale
! ------------
  
  
! ------------
! loop over all detection planes for which we have averaged wave functions
  do ipln=0, MSP_detpln-1
    islc = MSP_hdetpln(ipln) ! slice index
    if (ndet>0) MSP_detresult_ela(1:ndet, islc) = 0.0_fpp ! clear data
    if (nkmom>0) MSP_Kmomresult_ela(1:nkmom, islc) = 0.0_fpp
    if (MS_wave_avg_nac(ipln) > 0) then
      renorm = 1.0_fpp / real(MS_wave_avg_nac(ipln),kind=fpp) ! renormalization factor due to averaging
    else
      renorm = 1.0_fpp
    end if
    if (MS_wave_export_form==0) then ! real-space average -> need to transform
      ! real space export (inverse FT)
      ! - transfer data
      !   for an unknown reason, the following line causes a stack overflow and access violation
      !   MS_work(1:nx,1:ny) = complex(MS_wave_avg(1:nx,1:ny, ipln), kind=4)
      !   Though, the explicit assignement below element by element works.
      do j=1, ny
        do i=1, nx
          MS_work(i,j) = MS_wave_avg(i,j, ipln)
        end do
      end do
      ! call MS_FFT(work,MS_dimx,MS_dimy,'forwards') ! RS -> FS
      call MS_FFT_WORK(1)
      wave(1:nx,1:ny) = MS_work(1:nx,1:ny) * sqrt(rsca) * renorm ! renormalize after iDFT and averaging
    else ! Fourier-space averages, can take data as is, but renormalize
      wave(1:nx,1:ny) = MS_wave_avg(1:nx,1:ny, ipln) * renorm
    end if
    ! calculate the diffraction data stream used by all detectors and integrators below
    do j=1, ny
      idy = (j-1)*nx
      do i=1, nx
        idx = i + idy
        MSP_pdiftmp(idx) = real(wave(i,j)*conjg(wave(i,j)),kind=fpp) ! get power of wave function as stream of data
      end do
    end do
    ! default STEM detectors
    if (ndet>0) then
      ! readout all integrating detectors
      do k=1, ndet 
        nmlen = MSP_detmasklen(k)
        if ( 0 == nint(MSP_detdef(0,k)) .or. 0 >= nmlen) cycle ! skip empty or undefined detector
        do i=1, nmlen ! readout each masked Fourier pixel
          idx = MSP_detmask(i,k) ! get stream index
          MSP_pdettmp(i) = MSP_pdiftmp(idx) * MSP_detarea(idx,k) ! store detected intensity
        end do
        !call FDNCS2M(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a 2-fold butterfly
        call FPPSTRSUM(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a double precision accumulator
        MSP_detresult_ela(k, islc) = rval ! store the elastic channel intensity
      end do
    end if
  ! k-space moments
    if (nkmom>0) then
      ! readout k-moments
      idy = 0
      nmlen = MSP_Kmommasklen ! assuming everything is OK with MSP_Kmommasklen
      ! accumulate sums
      do m=0, MSP_KmomMmax ! loop from 0 to maximum order of moments
        do l=0, m ! loop through all 2d components of moment k
          idy = idy + 1
          k = m - l ! second parameter power
          do i=1, nmlen ! readout each masked Fourier pixel
            idx = MSP_Kmommask(i) ! get stream index
            gxk = MSP_Kmomgx(MSP_kmomhash(1,i),k) ! get gx^k
            gyl = MSP_Kmomgy(MSP_kmomhash(2,i),l) ! get gy^l
            MSP_pdettmp(i) = gxk * gyl * MSP_pdiftmp(idx) * MSP_Kmomwgt(idx) ! store detected intensity
          end do
          !call FDNCS2M(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a 2-fold butterfly
          call FPPSTRSUM(MSP_pdettmp(1:nmlen), nmlen, rval) ! sum up the intensities with a double precision accumulator
          MSP_Kmomresult_ela(idy, islc) = rval ! store elastic channel result
          if (m>0) then ! normalize by the 0-th moment of the total (!) MSP_Kmomresult(1, islc)
            MSP_Kmomresult_ela(idy, islc) = MSP_Kmomresult_ela(idy, islc) / MSP_Kmomresult(1, islc)
          end if
        end do
      end do
    end if
  end do
! ------------
  
! ------------
! normal routine exit
100 continue
  nerr = 0
  goto 1000 ! now really exit
! error exits (handle different error cases)
200 continue
  select case (nerr)
  case (1)
    call PostMessage("DetectorReadoutElastic failed: MultiSlice module not initialized.")
    ncrit = 1
  case (2)
    call PostMessage("DetectorReadoutElastic failed: Invalid grid size.")
    ncrit = 1
  case (3)
    call PostMessage("DetectorReadoutElastic failed: Memory allocation failed.")
    ncrit = 1
  case default
    call PostMessage("DetectorReadoutElastic failed due to unknown reason.")
  end select ! case (nerr)
  goto 1000
! final exit
1000 continue
  ! deallocate
  if (allocated(wave)) deallocate(wave, stat=nalloc)
  if (ncrit==1) then  ! trigger critical error stopping the program
    call CriticalError("DetectorReadoutElastic failed.")
  end if
  return ! exit

END SUBROUTINE DetectorReadoutElastic
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE ApplySpatialCoherence()
! purpose: Does source function convolution, expects real data
!          with given stem image size.
!          Will try to load many images sequentially from the input
!          file, convolute each, and stores them back to the output.
! -------------------------------------------------------------------- !
! parameter: 
!          None.
! -------------------------------------------------------------------- !

  use MultiSlice
  use STEMfunctions
  use MSAparams
  
  implicit none

! ------------
! DECLARATION
  integer*4 :: nlfu, mlfu, nerr, nalloc, npln, ioerr
  real(fpp) :: sx, sy
  real(fpp), allocatable :: rdata(:,:)
  external :: createfilefolder
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ApplySpatialCoherence: INIT."
  nerr = 0
  ioerr = 0
  nalloc = 0
  if (MSP_ctemmode==1) then
    call PostWarning("Inconsistent parameters: Deactivate either"//&
     & " CTEM mode or application of partial spatial coherence.")
    call PostMessage("Aborting application of partial spatial coherence")
    return
  end if
  ! calculate sampling of multislice image
  sx = real(MSP_SF_sizex/MSP_SF_ndimx,kind=fpp)
  sy = real(MSP_SF_sizey/MSP_SF_ndimy,kind=fpp)
  ! allocate data array for loading images
  allocate(rdata(MSP_SF_ndimx,MSP_SF_ndimy), stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("Failed to allocate memory for reading data.")
    goto 99
  end if
  rdata = 0.0
! ------------

! ------------
! - connect input file
  call GetFreeLFU(nlfu,20,100)
  open(unit=nlfu, file=trim(MSP_infile), form='binary', &
     & access='sequential', iostat=ioerr, status="old", &
     & action="read", share='DENYNONE' )
  if (ioerr/=0) then
    nerr = 1
    goto 99
  end if
  call PostMessage("Reading STEM image from file ["//trim(MSP_infile)//"].")
! - connect output file
  call GetFreeLFU(mlfu,20,100)
  call createfilefolder(trim(MSP_outfile),nerr)
  open(unit=mlfu, file=trim(MSP_outfile), form='binary', access='sequential', &
     & iostat=nerr, status="replace", action="write", share='DENYRW' )
  if (ioerr/=0) then
    nerr = 2
    goto 99
  end if
  call PostMessage("Writing convoluted STEM image to file ["//trim(MSP_outfile)//"].")
! ------------

! ------------
! - read - convolution - write loop
  npln = 0
  do ! repeated reading trials
    !
    ! - read data
    !
    read(unit=nlfu, iostat=ioerr) rdata
    if (ioerr/=0) exit ! read failed: no (more) data to handle
    !
    ! - convolute data
    !
    npln = npln + 1
    write(unit=MSP_stmp, fmt=*) npln
    call PostMessage("- convoluting image "//trim(adjustl(MSP_stmp)))
    ioerr = MS_err_num
    call MS_ApplyPSpatialCoh(rdata,MSP_SF_ndimx,MSP_SF_ndimy,&
     &                       sx,sy,STF_srcradius,MSP_PC_spatial)
    if (ioerr/=MS_err_num) then
      nerr = 3
      goto 99
    end if
    !
    ! - write data
    !
    write(unit=mlfu, iostat=ioerr) rdata
    if (ioerr/=0) then
      nerr = 4
      goto 99
    end if
    !
  end do ! loop
  !
  ! disconnect from files
  close(unit=mlfu, iostat=ioerr)
  close(unit=nlfu, iostat=ioerr)
  !
  call PostMessage("STEM image convolution finished.")
  !
! ------------


! ------------
! exit point (in any case jump here finally)
99 if (allocated(rdata)) deallocate(rdata, stat=nalloc)
! finishing and critical error handling 
!  write(unit=*,fmt=*) " > ApplySpatialCoherence: EXIT."
  select case (nerr)
    case (1)
      call CriticalError("Failed to connect to file ["// &
         & trim(MSP_infile)//"].")
    case (2)
      close(unit=nlfu, iostat=ioerr)
      call CriticalError("Failed to connect to file ["// &
         & trim(MSP_outfile)//"].")
    case (3)
      close(unit=nlfu, iostat=ioerr)
      close(unit=mlfu, iostat=ioerr)
      call CriticalError("Image convolution failed.")
    case (4)
      close(unit=nlfu, iostat=ioerr)
      close(unit=mlfu, iostat=ioerr)
      call CriticalError("Failed to write output.")
  end select ! case (nerr)
  
  return

END SUBROUTINE ApplySpatialCoherence
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE MSACalculate()

  use MultiSlice
  use MSAparams

  implicit none
  
  integer*4 :: nerr, nx, ny, nslc, nvc, nalloc, nslcidx
  integer*4, allocatable :: lvar(:)
  
  real*4, external :: UniRand
  
  nerr = 0
  nalloc = 0
  allocate(lvar(MS_stacksize),stat=nalloc)
  if (nalloc/=0) goto 99
  nx = MSP_dimcellx
  ny = MSP_dimcelly

! setup new multislice
  ! calculate variance Monte-Carlo for the next run (draw without returning)
  call MSP_GetVarList(MS_slicenum, MS_stacksize, MS_slicestack, lvar, nerr)
  if (nerr/=0) goto 99
  nerr = MS_err_num
  if (MSP_use_extinwave==1) then
    call MS_Start(MSP_extinwslc)
  else
    call MS_Start()
  end if
  if (nerr/=MS_err_num) goto 99
  
  ! stack load on demand is not clever, removed 2019-11-06 v.1.0.6 in favor of single slice load-on-demand
  !if (MSP_SLC_lod == 1) then ! load slice data on demand
  !  call MSP_LoadStack(lvar, nerr)
  !  if (nerr/=0) goto 99
  !end if

! loop over all slices the multislice
  do while (MS_slicecur >= 0.and. MS_slicecur<MS_stacksize)
    nerr = MS_err_num
    nslc = MS_slicestack(MS_slicecur+1)+1 ! get current slice index
    !nvc = 1 + int( UniRand()*real(MSP_SLC_setup(0,nslc))*0.9999 ) ! get current variant index
    nvc = 1 + lvar(1+MS_slicecur)
        
    ! slice
    if (MSP_SLC_lod==0) then ! get slice index in pre-loaded set
      nslcidx = MSP_GetPGRIndex(nvc,nslc,nerr)
      if (nerr/=0) goto 99
    else ! load-on-demand
      nslcidx = 1
      call MSP_LoadPGR(nslc, lvar(1+MS_slicecur), nslcidx, nerr)
      if (nerr/=0) goto 99
    end if
    
    call MS_CalculateNextSlice(MSP_phasegrt(1:nx, 1:ny, nslcidx), nx, ny)
    if (nerr/=MS_err_num) goto 99
  end do ! while (MS_slicecur >= 0)

! Stop the multislice
  nerr = MS_err_num
  call MS_Stop()
  if (nerr/=MS_err_num) goto 99
  
99 if (allocated(lvar)) deallocate(lvar, stat=nalloc)
  return

END SUBROUTINE MSACalculate
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
SUBROUTINE STEMMultiSlice()
! function: + variants
!           + focus spread
!           + source spread
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  integer*4 :: nx, ny, nz, nerr, nznum, ndet, nkmom, ndat, nalloc, nslcidx
  integer*4 :: nv, nvc, nvar, nvarnum, nvartot, nvdigits
  integer*4 :: nslc, ncalcslc
  real(fpp) :: scansampx, scansampy, scanposx, scanposy
  real(fpp) :: dxcurr, dycurr, dzcurr
  real(fpp) :: zstep, zoffset, zpow, zrescale, fafac, fascal, vrescale
  real(fpp) :: ffac
  real(fpp), allocatable :: rtmpresult(:)
  integer*4, allocatable :: lvar(:)
  character(len=1000) :: swavfile, stmp
  real*4 :: ffs, fdz, fsc, fdx, fdy
  real*4, external :: UniRand
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > STEMMultiSlice: INIT."
  nx = MSP_dimcellx
  ny = MSP_dimcelly
  ndet = MSP_detnum
  nkmom = 0
  if (MSP_Kmomout>0) nkmom = MSP_KmomNum
  ndat = ndet + nkmom
  allocate(rtmpresult(ndat), lvar(MS_stacksize),stat=nalloc)
  scansampx = 0.0_fpp
  if (MSP_SF_sizex>0.0.and.MSP_SF_ndimx>1) then
    scansampx = MSP_SF_sizex/real(MSP_SF_ndimx,kind=fpp)
  end if
  scansampy = 0.0_fpp
  if (MSP_SF_sizey>0.0.and.MSP_SF_ndimy>1) then
    scansampy = MSP_SF_sizey/real(MSP_SF_ndimy,kind=fpp)
  end if
  scanposx =   MSP_SF_offsetx &
     &       + MSP_SF_rotcos*MSP_ScanPixelX*scansampx &
     &       - MSP_SF_rotsin*MSP_ScanPixelY*scansampy
  scanposy =   MSP_SF_offsety &
     &       + MSP_SF_rotcos*MSP_ScanPixelY*scansampy &
     &       + MSP_SF_rotsin*MSP_ScanPixelX*scansampx
  ffs = real(STF_Defocusspread, kind=4)
  fsc = real(STF_srcradius, kind=4)
  fdz = 0.0
  fdx = 0.0
  fdy = 0.0
! ------------

! ------------
! init explicit focal averaging
  nznum = 1
  zstep = 0.0_fpp
  zoffset = 0.0_fpp
  zrescale = 1.0_fpp
  zpow = 0.0_fpp
  fafac = 0.0_fpp
  fascal = 0.0_fpp
  if (MSP_PC_temporal/=0 .and. MSP_ExplicitPSC==0) then
    nznum = STF_DEFOCUS_KERNEL_STEPS
    zoffset = -STF_DEFOCUS_KERNEL_SPREAD*STF_defocusspread
    zstep = -2.0_fpp*zoffset / real( nznum-1 , kind=fpp)
    fafac = -1.0_fpp/STF_defocusspread/STF_defocusspread
    write(unit=MSP_stmp,fmt=*) "Performing explicit focal averaging:"
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt=*) "  focal range [nm]: +/-",STF_DEFOCUS_KERNEL_SPREAD*STF_defocusspread
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt=*) "  number of steps:",nznum
    call PostMessage(trim(MSP_stmp))
    write(unit=MSP_stmp,fmt=*) "  focal change per step [nm]:",zstep
    call PostMessage(trim(MSP_stmp))
  end if
  if (MSP_ExplicitPSC/=0) then
    call PostMessage("Performing explicit averaging for:")
    if (MSP_PC_temporal/=0) then
      call PostMessage("- focus spread")
    end if
    if (MSP_PC_spatial/=0) then
      call PostMessage("- source size")
    end if
    nznum = 1 ! no fix focus loop
    zstep = 0.0_fpp
    zoffset = 0.0_fpp
    zrescale = 1.0_fpp
    zpow = 0.0_fpp
    fafac = 0.0_fpp
    fascal = 0.0_fpp
  end if
! ------------


! ------------
! init multi-slice core
  MSP_detresult = 0.0_fpp
  MSP_Kmomresult = 0.0_fpp
  rtmpresult = 0.0_fpp
  MSP_TheResult = 0.0_fpp
  nvarnum = max(1,ceiling(real(MSP_FL_varcalc)/real(nznum))) ! number of variant calculations per focus spread loop
  write(unit=stmp,fmt='(I)') nvarnum
  nvdigits = MAX( 3, LEN_TRIM(adjustl(stmp)) )
  vrescale = 1.0_fpp/real(nvarnum, kind=fpp)
  nvar = 0
  nvartot = nvarnum*nznum
  swavfile = trim(MS_wave_filenm) ! wave file name backup
  
  ! The following code does nznum * nvarnum multislice passes
  
  ! loop through focal variants
  do nz=1, nznum
    
    ! PROBE WAVEFUNCTION SHIFTS FOR PURE FOCAL KERNEL MODE
    dzcurr = zoffset + zstep*real(nz-1, kind=fpp)
    nerr = MS_err_num
    if (MSP_use_extinwave==1) then ! use external wavefunction
      call MS_OffsetIncomingWave(0.0_fpp,0.0_fpp,0.0_fpp) ! the inserted wave is inside the crystal, no shift and focus is needed
    else
      if (MSP_ExplicitPSC==0) then ! explicit focal kernel is used
        call MS_OffsetIncomingWave(scanposx, scanposy, dzcurr)
        if (nerr/=MS_err_num) goto 11
        if (DEBUG_EXPORT>0) then
          if (MSP_use_extinwave==1) then
            call PostMessage("Inserted wavefunction was applied.")
          else
            call PostMessage("Offset applied to probe wavefunction.")
            write(unit=MSP_stmp,fmt='(A,2G13.4)') "- scan shift: ", scanposx, scanposy
            call PostMessage(trim(MSP_stmp))
            if (MSP_PC_temporal/=0) then
              write(unit=MSP_stmp,fmt='(A, G13.4)') "- extra defocus  : ", dzcurr
              call PostMessage(trim(MSP_stmp))
            end if
          end if
        end if
      end if
    end if
    
    ! loop through frozen lattice variants
    do nv = 1, nvarnum ! variants loop per focus spread cycle
    
      if (MSP_use_extinwave==0 .and. MSP_ExplicitPSC==1) then
        ! PROBE WAVEFUNCTION SHIFTS FOR EXPLICIT PARTIAL COHERENCE MODE
        nerr = MS_err_num
        dzcurr = 0.0_fpp
        if (MSP_PC_temporal/=0) then
          call MSP_GetRandomProbeDefocus(ffs, fdz)
          dzcurr = real(fdz, kind=fpp)
        end if
        dxcurr = 0.0_fpp
        dycurr = 0.0_fpp
        if (MSP_PC_spatial/=0) then
          call MSP_GetRandomProbeShift(fsc, fdx, fdy)
          dxcurr = real(fdx, kind=fpp)
          dycurr = real(fdy, kind=fpp)
        end if
        call MS_OffsetIncomingWave(scanposx+dxcurr,scanposy+dycurr,dzcurr)
        if (nerr/=MS_err_num) goto 11
        if (DEBUG_EXPORT>0) then
          call PostMessage("Offset applied to probe wavefunction.")
          write(unit=MSP_stmp,fmt='(A,2G13.4)') "- scan shift: ",scanposx, scanposy
          call PostMessage(trim(MSP_stmp))
          if (MSP_PC_spatial/=0) then
            write(unit=MSP_stmp,fmt='(A, 2G13.4)') "- partial coherence shift  : ",dxcurr, dycurr
            call PostMessage(trim(MSP_stmp))
          end if
          if (MSP_PC_temporal/=0) then
            write(unit=MSP_stmp,fmt='(A, G13.4)') "- partial coherence defocus  : ",dzcurr
            call PostMessage(trim(MSP_stmp))
          end if
        end if
      end if
    
      ! PERFORM FULL MULTISLICE WITH CURRENT WAVE
      !
      ! <-- SINGLE MULTISLICE STARTS HERE
      !
      nerr = MS_err_num
      write(unit=MSP_stmp,fmt='(A,I6,A,I6)') "- repeat: ",nv + (nz-1)*nvarnum, " / ", nznum*nvarnum
      call PostMessage(trim(MSP_stmp))
      
      ! update wave file name (insert variation number) when more than one variant per multislice
      if (MS_wave_export>=1 .and. nvartot>1) then
        call saddsuffix(trim(swavfile), "_vr", nvar+1, nvdigits, MSP_stmp)
        MS_wave_filenm = trim(MSP_stmp)
      end if
      
      ! calculate variance Monte-Carlo for the next run (draw without returning)
      call MSP_GetVarList(MS_slicenum, MS_stacksize, MS_slicestack, lvar, nerr)
      if (nerr/=0) goto 14
      
      if (MSP_use_extinwave==1) then
        call MS_Start(MSP_extinwslc)
      else
        call MS_Start()
      end if
      if (nerr/=MS_err_num) goto 16
      
      ! stack load on demand is not clever, removed 2019-11-06 v.1.0.6 in favor of single slice load-on-demand
      !if (MSP_SLC_lod == 1) then ! load slice data on demand
      !  call MSP_LoadStack(lvar, nerr)
      !  if (nerr/=0) goto 17
      !end if
      
      ncalcslc = MS_slicecur ! reset number of calculated slices
      
      ! determine current scaling factor for summation over foci and over variants
      ffac = exp(fafac*dzcurr*dzcurr)*vrescale
      ! write(*,*) ffac, exp(fafac*dzcurr*dzcurr), vrescale
	    ! .. and sum up the applied weights for later rescaling to 1.0
	    fascal = fascal + ffac;
      
      if (MSP_ldetpln(ncalcslc) >= 0) then ! readout at incident plane
        ! COLLECT DATA FROM DETECTOR
        call DetectorReadout(rtmpresult, ndat, nerr)
        if (nerr/=0) goto 12
        ! sum up weighted data depending on slice number and detector number
	      MSP_detresult(1:ndet,ncalcslc) = MSP_detresult(1:ndet,ncalcslc) + rtmpresult(1:ndet)*ffac
        if (nkmom>0) then
          MSP_Kmomresult(1:nkmom,ncalcslc) = MSP_Kmomresult(1:nkmom,ncalcslc) + rtmpresult(ndet+1:ndet+nkmom)*ffac
        end if
        !
	      if (DEBUG_EXPORT>0) then
	        write(unit=MSP_stmp,fmt='(A,I4,A)') "- Detection after ",ncalcslc," slices."
          call PostMessage(trim(MSP_stmp))
          write(unit=MSP_stmp,fmt='(A,<ndet>G13.5)') "- Current detector readout:",rtmpresult(1:ndet)
          call PostMessage(trim(MSP_stmp))
          write(unit=MSP_stmp,fmt='(A,<ndet>G13.5)') "- Current total result    :",MSP_detresult(1:ndet,ncalcslc)
          call PostMessage(trim(MSP_stmp))
        end if
      end if
      
      !
      ! loop over all slices and do the multislice
      !    
      do while (MS_slicecur >= 0.and. MS_slicecur<MS_stacksize)
        nerr = MS_err_num ! backup error
        nslc = MS_slicestack(MS_slicecur+1)+1 ! get current slice index
        !nvc = 1 + int( UniRand()*real(MSP_SLC_setup(0,nslc))*0.9999 ) ! get current variant index
        nvc = 1 + lvar(1+MS_slicecur)
        
        ! slice
        if (MSP_SLC_lod==0) then ! get slice index in pre-loaded set
          nslcidx = MSP_GetPGRIndex(nvc,nslc,nerr)
          if (nerr/=0) goto 15
        else ! load-on-demand
          nslcidx = 1
          call MSP_LoadPGR(nslc, lvar(1+MS_slicecur), nslcidx, nerr)
          if (nerr/=0) goto 17
        end if
        call MS_CalculateNextSlice(MSP_phasegrt(1:nx, 1:ny, nslcidx), nx, ny)
        ncalcslc = ncalcslc + 1 ! increase number of calculated slices
        
        if (nerr/=MS_err_num) goto 16 ! error check
        
        if (MSP_ldetpln(ncalcslc) >= 0) then ! detection?
          ! COLLECT DATA FROM DETECTOR
          call DetectorReadout(rtmpresult, ndat, nerr)
          if (nerr/=0) goto 12
          ! sum up weighted data depending on slice number and detector number
	        MSP_detresult(1:ndet,ncalcslc) = MSP_detresult(1:ndet,ncalcslc) + rtmpresult(1:ndet)*ffac
          if (nkmom>0) then
            MSP_Kmomresult(1:nkmom,ncalcslc) = MSP_Kmomresult(1:nkmom,ncalcslc) + rtmpresult(ndet+1:ndet+nkmom)*ffac
          end if
          !
	        if (DEBUG_EXPORT>0) then
	          write(unit=MSP_stmp,fmt='(A,I4,A)') "- Detection after ",ncalcslc," slices."
            call PostMessage(trim(MSP_stmp))
            write(unit=MSP_stmp,fmt='(A,<ndet>G13.5)') "- Current detector readout:",rtmpresult(1:ndet)
            call PostMessage(trim(MSP_stmp))
            write(unit=MSP_stmp,fmt='(A,<ndet>G13.5)') "- Current total result    :",MSP_detresult(1:ndet,ncalcslc)
            call PostMessage(trim(MSP_stmp))
          end if
          !
        end if ! detector readout
        
      end do ! while (MS_slicecur >= 0.and. MS_slicecur<MS_stacksize)

      ! Stop the multislice
      nerr = MS_err_num
      call MS_Stop()
      if (nerr/=MS_err_num) goto 16
      !
      ! <-- SINGLE MULTISLICE ENDS HERE
      !
      
      nvar = nvar + 1 ! increase number of applied variants
    
    end do ! variants loop

  end do ! focus spread loop
! ------------

! ------------
! rescale result
  zrescale = 1.0_fpp/fascal
  MSP_detresult = MSP_detresult * zrescale
  if (nkmom>0) then
    MSP_Kmomresult = MSP_Kmomresult * zrescale
  end if
! ------------
  
  
! ------------
! elastic channel output
  if (MS_wave_avg_export>0) then
    call DetectorReadoutElastic(nerr)
    if (nerr/=0) goto 13
  end if
! ------------


! ------------
  if (MSP_txtout==0) call ExportSTEMData(trim(MSP_outfile))
  call MSP_WriteTextOutput(nerr)
! ------------


! ------------
! normal/successful exit
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  return
  
! ------------
! failure exits
11 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to shift incoming wavefunction.")
  return

12 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to readout detectors.")
  return
  
13 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to readout detectors for elastic channel.")
  return
  
14 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to determine random variant sequence.")
  return

15 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to connect to phase-grating buffer.")
  return

16 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to perform multislice algorithm.")
  return
  
17 continue
  if (allocated(rtmpresult)) deallocate(rtmpresult,stat=nalloc)
  if (allocated(lvar)) deallocate(lvar,stat=nalloc)
  call CriticalError("Failed to load slice data on demand.")
  return

END SUBROUTINE STEMMultiSlice
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE CTEMMultiSlice()
! function: + variants
!           - focus spread
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  integer*4 :: nerr, nv, nvarnum, nvd
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > CTEMMultiSlice: INIT."
  nerr = 0
  MSP_PC_temporal = 0
! ------------


! ------------
  nvarnum = MSP_FL_varcalc ! number of variant calculations
  
  nvd = MSP_nvard
    
  ! OFFSET THE BACKUP WAVE AND GENERATE ACTIVE INCOMING WAVE
  call MS_OffsetIncomingWave(0.0_fpp,0.0_fpp,0.0_fpp)
    
    
  do nv = 1, nvarnum ! variants loop per focus spread cycle
  
    ! update wave file name (insert variation number) when more than one variant per multislice
    if (MS_wave_export > 0 .and. nvarnum > 1) then
      MS_wave_filenm = MS_wave_filenm_bk
      call saddsuffix(trim(MS_wave_filenm), "_vr", nv, nvd, MSP_stmp)
      MS_wave_filenm = trim(MSP_stmp)
    end if
    
    ! PERFORM FULL MULTISLICE WITH CURRENT WAVE
    nerr = MS_err_num
    call MSACalculate()
    if (nerr/=MS_err_num) then
      call CriticalError("Failed to perform multislice algorithm with current wave.")
    end if

  end do ! variants loop
  MS_wave_filenm = MS_wave_filenm_bk ! reset wave file name
! ------------


! ------------
!  write(unit=*,fmt=*) " > CTEMMultiSlice: EXIT."
  return

END SUBROUTINE CTEMMultiSlice
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE CreateSTEMFile(sfile,ndata,nerr)
! functions: Creates a file for stem data output with sufficient size
!            to write all data into it later.
! -------------------------------------------------------------------- !
! parameter:
!   character(len=*) :: sfile = file name
!   integer*4 :: ndata = number of data (assuming 4 bytes per data)
!   integer*4 :: nerr = error code in case the routine fails
! -------------------------------------------------------------------- !

  use MSAparams

  implicit none
  
  character(len=*), intent(in) :: sfile ! input file name
  integer*4, intent(in) :: ndata ! total number of 4-byte data items in the file
  integer*4, intent(inout) :: nerr ! error code
  integer*4 :: lfu, status
  real(fpp), allocatable :: tmpdat(:)
  
  ! init
  nerr = 0

  call PostMessage("- Creating new output file ["//trim(sfile)//"].")
  ! get new logical fle unit
  call GetFreeLFU(lfu,40,1000)
  ! open the file and replace the old file / create new file
  call createfilefolder(trim(sfile),status)
  ! record length is equal to total data size of (4 byte)*datanum
  open (unit=lfu, file=trim(sfile), form="binary",iostat=status,&
     &  access="direct", recl=fpp*ndata, status="replace", &
     &  action="write", share='DENYRW')
  if (status/=0) then
    call CriticalError("Output file creation failed.")
    nerr = 1
    return
  end if
  
  ! write zeroes to create the file of complete size with zero values preset
  allocate(tmpdat(ndata),stat=status)
  if (status/=0) then
    write(unit=MSP_stmp,fmt='(A,I8,A)') &
     &    "Failed to allocate memory (",int(ndata*4/1024),"kB)."
    call CriticalError(trim(MSP_stmp))
    nerr = 2
    return
  end if
  tmpdat = 0.0
  write(unit=lfu,rec=1,iostat=status) tmpdat
  deallocate(tmpdat,stat=status)
  close (unit=lfu)
  return
   
END SUBROUTINE CreateSTEMFile
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE ExportSTEMData(sfile)
! function: writes calculated signal of one scan pixel to raw binary
!           output file(s).
! -------------------------------------------------------------------- !
! parameter: file name
! -------------------------------------------------------------------- !

  use MSAparams
  use MultiSlice

  implicit none

! ------------
! DECLARATION
  character*(*), intent(in) :: sfile
  integer*4 :: nfil, lfu(3), nerr, i, j, k, l, l1, iff
  integer*4 :: datapos, ndet, ipos
  integer*4 :: ndatanum
  real(fpp) :: rsignal(3)
  logical :: fex(3)
  character(len=12) :: snumk, snuml
  character(len=1024) :: stmp, ssufdet, ssufsep(3), spfile(3), sdfile(3)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ExportSTEMData: INIT."
  write(unit=stmp,fmt='(A,I4,A,I4)') &
     &  "Saving detector signals for scan pixel ", &
     &  MSP_ScanPixelX,", ",MSP_ScanPixelY
  call PostMessage(trim(stmp))
  ! check data position in the image
  datapos = 1 + MSP_ScanPixelX + MSP_ScanPixelY*MSP_SF_ndimx
  ! determine number of scan points
  ndatanum = MSP_SF_ndimx*MSP_SF_ndimy
  !
  if (datapos<=0.or.datapos>ndatanum) then
    call CriticalError("ExportSTEMData: Invalid datafile position.")
  end if
  nfil = 1
  ssufsep = ""
  if (MS_wave_avg_export > 0) then ! separate signal scattering channels
    nfil = 3
    call PostMessage("Saving also separated signal of elastic and thermal-diffuse scattering.")
    ssufsep(1) = "_tot"
    ssufsep(2) = "_ela"
    ssufsep(3) = "_tds"
  end if
  ssufdet = ""
  rsignal = 0.0_fpp
! ------------


! ------------
!  STEM DETECTOR OUTPUT
! ------------
  ! get number of detectors
  ndet = max(1,MSP_detnum)
  ! loop over detectors
  do k=1, ndet
    ! prepare the output file name (append detector name and output channel)
    sdfile = trim(sfile) ! default preset
    if (MSP_usedetdef/=0) then ! update output file name with detector name
      ssufdet = "_"//trim(MSP_detname(k))
    end if
    do iff=1, nfil ! add detector suffic and separation suffix for each file to open
      call saddsuffix(trim(sfile), trim(ssufdet)//trim(ssufsep(iff)), 0, 0, sdfile(iff))
    end do
    ! write data to file
    if (MSP_3dout > 0) then ! prepare for output to 3d data file
      !
      ! BEGIN OF 3D FILE OUTPUT
      !
      do iff=1, nfil ! open all files
        ! - check existence of current output file
        inquire(file=trim(sdfile(iff)),exist=fex(iff))
        if (.not.fex(iff)) then ! doesn't exist, create new
          call CreateSTEMFile(trim(sdfile(iff)),ndatanum*MSP_detpln,nerr)
          if (nerr/=0) call CriticalError("Output file creation failed.")
        end if
        ! - get logical unit
        call GetFreeLFU(lfu(iff),20,100)
        ! - open file shared access
        open (unit=lfu(iff), file=trim(sdfile(iff)), form="binary", access="direct", &
            & iostat=nerr, status="old", action="write", recl=fpp, share='DENYNONE' )
        if (nerr/=0) then
          call CriticalError("ExportSTEMData: Failed to open file ["//trim(sdfile(iff))//"].")
        end if
      end do
      ! - write to the files ... 
      ! preset plane offset for 3d data stacks
      do i=0, MS_stacksize ! loop over all slices beginning with incident probe plane (0)
        if (MSP_ldetpln(i) < 0) cycle ! skip this slice
        ipos = datapos + MSP_ldetpln(i)*ndatanum
        rsignal(1) = MSP_detresult(k,i)
        if (nfil>1) then
          rsignal(2) = MSP_detresult_ela(k,i) ! ela
          rsignal(3) = rsignal(1) - rsignal(2) ! tds = tot - ela
        end if
        ! - write the data at correct positions
        do iff=1, nfil
          write(unit=lfu(iff), rec=ipos, iostat=nerr) rsignal(iff)
        end do
        ! - report stored signal
        write (unit=MSP_stmp,fmt='(A,<nfil>G13.5,A)') "- saved "// &
          &   trim(MSP_detname(k))//" signal (",rsignal(1:nfil),")"
        call PostMessage(trim(MSP_stmp))
      end do ! i-loop over slices
      !
      ! - close the files ...
      do iff=1, nfil
        close(unit=lfu(iff), iostat=nerr)
      end do
      !
      ! END OF 3D FILE OUTPUT
      !
    else ! perpare for output to 2D files (one per export plane)
      !
      ! BEGIN OF SINGLE PLANE FILE OUTPUT
      !
      ! loop over all slices
      do i=0, MS_stacksize ! loop over all slices beginning with incident probe plane (0)
        if (MSP_ldetpln(i) < 0) cycle ! skip this slice
        spfile = sdfile ! default preset
        ! - modify file names (append slice index)
        if (MSP_detpln > 1) then ! append slice index to file name
          do iff=1, nfil
            call saddsuffix(trim(sdfile(iff)), "_sl", i, MSP_nslid, spfile(iff))
          end do
        end if
        ! - open all files
        do iff=1, nfil
          ! - check existence of current output file
          inquire(file=trim(spfile(iff)),exist=fex(iff))
          if (.not.fex(iff)) then ! doesn't exist, create new (single plane file)
            call CreateSTEMFile(trim(spfile(iff)),ndatanum,nerr)
            if (nerr/=0) call CriticalError("Output file creation failed.")
          end if
          ! - get logical unit
          call GetFreeLFU(lfu(iff),20,100)
          ! - open file shared access
          open (unit=lfu(iff), file=trim(spfile(iff)), form="binary", access="direct", &
           &    iostat=nerr, status="old", action="write", recl=fpp, share='DENYNONE' )
          if (nerr/=0) then
            call CriticalError("ExportSTEMData: Failed to open file ["//trim(spfile(iff))//"].")
          end if
        end do
        ! - set signal to store
        rsignal(1) = MSP_detresult(k,i) ! tot
        if (nfil>1) then
          rsignal(2) = MSP_detresult_ela(k,i) ! ela
          rsignal(3) = rsignal(1) - rsignal(2) ! tds = tot - ela
        end if
        ! - write the signal at correct positions ... and close the files
        do iff=1, nfil
          write(unit=lfu(iff), rec=datapos, iostat=nerr) rsignal(iff)
          close(unit=lfu(iff), iostat=nerr)
        end do
        ! - report stored signal
        write (unit=MSP_stmp,fmt='(A,<nfil>G13.5,A)') "- saved "// &
          &   trim(MSP_detname(k))//" signal (",rsignal(1:nfil),")"
        call PostMessage(trim(MSP_stmp))
      end do ! i-loop over slices
      !
      ! END OF SINGLE PLANE FILE OUTPUT
      !
    end if ! SWITCH 3d output or single plane files

  end do ! k-loop over detectors
! ------------
  
! ------------
!  k-MOMENT DETECTOR OUTPUT
!  This is an experimental mode
!  - datapos and ndatanum are modified
!  - for each momentum order, ncomp = order+1 components are stored as ncomp-tuples per scan position
! ------------
  if (MSP_Kmomout>0) then
    call PostMessage("Integrated k-space momentum output.")
    ! loop over all moment orders
    j = 0 ! offset of current order in result arrays MSP_Kmomresult etc.
    do k=0, MSP_KmomMmax ! ... from zero to m_max
      write(unit=snumk,fmt=*) k ! order to string
      do l=0, k ! component in current order
        l1 = l + 1 ! shifted index
        ! prepare the output file name
        write(unit=snuml,fmt=*) l ! component to string
        do iff=1, nfil ! construct file names ... 
          call saddsuffix(trim(sfile), "_kmom"//trim(adjustl(snumk))// &
            & "-"//trim(adjustl(snuml))//trim(ssufsep(iff)), 0, 0, &
            & sdfile(iff))
        end do
        ! write data to file
        if (MSP_3dout > 0) then ! prepare for output to 3d data file
          !
          ! BEGIN OF 3D FILE OUTPUT
          !
          do iff=1, nfil ! open all files
            ! - check existence of current output file
            inquire(file=trim(sdfile(iff)),exist=fex(iff))
            if (.not.fex(iff)) then ! doesn't exist, create new
              call CreateSTEMFile(trim(sdfile(iff)),ndatanum*MSP_detpln,nerr)
              if (nerr/=0) call CriticalError("Output file creation failed.")
            end if
            ! - get logical unit
            call GetFreeLFU(lfu(iff),20,100)
            ! - open file shared access
            open (unit=lfu(iff), file=trim(sdfile(iff)), form="binary", &
              &   access="direct", iostat=nerr, status="old", &
              &   action="write", recl=fpp, share='DENYNONE' )
            if (nerr/=0) then
              call CriticalError("ExportSTEMData: Failed to open file [" &
                &  //trim(sdfile(iff))//"].")
            end if
          end do
          ! - write to the file ... 
          do i=0, MS_stacksize ! loop over all slices beginning with incident probe plane (0)
            if (MSP_ldetpln(i) < 0) cycle ! skip this slice
            ! - write the total intensity data at correct position to files
            ipos = datapos + MSP_ldetpln(i)*ndatanum
            ! - get the signal of current plane and component
            rsignal(1) = MSP_Kmomresult(l1+j, i) ! total signal
            if (nfil>1) then ! store also ela and tds
              rsignal(2) = MSP_Kmomresult_ela(l1+j, i) ! ela
              rsignal(3) = rsignal(1) - rsignal(2) ! tds = tot - ela
            end if
            ! - store
            do iff=1, nfil
              write(unit=lfu(iff), rec=ipos, iostat=nerr) rsignal(iff)
            end do
            ! - report stored signal values
            write (unit=MSP_stmp,fmt='(A,<nfil>G13.5,A)') &
              &   "- saved moment "//trim(adjustl(snumk))// &
              &   "("//trim(adjustl(snuml))//") signal (", &
              &   rsignal(1:nfil),")"
            call PostMessage(trim(MSP_stmp))
          end do ! i-loop over slices
          ! - close the files ...
          do iff=1, nfil
            close(unit=lfu(iff), iostat=nerr)
          end do
          !
          ! END OF 3D FILE OUTPUT
          !
        else ! perpare for output to 2D files (one per export plane)
          !
          ! BEGIN OF SINGLE PLANE FILE OUTPUT
          !
          ! loop over all slices
          do i=0, MS_stacksize ! loop over all slices beginning with incident probe plane (0)
            if (MSP_ldetpln(i) < 0) cycle ! skip this slice
            spfile = sdfile ! default preset
            ! - modify file name with slice index
            if (MSP_detpln > 1) then ! append slice index to file name
              ! - update output file names ...
              do iff=1, nfil
                call saddsuffix(trim(sdfile(iff)), "_sl", i, MSP_nslid, spfile(iff))
              end do
            end if
            ! - open the files ...
            do iff=1, nfil
              ! - check existence of current output files
              inquire(file=trim(spfile(iff)),exist=fex(iff))
              if (.not.fex(iff)) then ! doesn't exist, create new (single plane files)
                call CreateSTEMFile(trim(spfile(iff)),ndatanum,nerr)
                if (nerr/=0) call CriticalError("Output file creation failed.")
              end if
              ! - get logical unit
              call GetFreeLFU(lfu(iff),20,100)
              ! - open file shared access
              open (unit=lfu(iff), file=trim(spfile(iff)), form="binary",&
                &   access="direct", iostat=nerr, status="old", &
                &   action="write", recl=fpp, share='DENYNONE' )
              if (nerr/=0) then
                call CriticalError("ExportSTEMData: Failed to open file ["&
                  &  //trim(spfile(iff))//"].")
              end if
            end do
            ! - get the signal of current plane and component
            rsignal(1) = MSP_Kmomresult(j+l1, i) ! total signal
            if (nfil>1) then ! store also ela and tds
              rsignal(2) = MSP_Kmomresult_ela(j+l1, i) ! ela
              rsignal(3) = rsignal(1) - rsignal(2) ! tds = tot - ela
            end if
            ! - store
            do iff=1, nfil 
              write(unit=lfu(iff), rec=datapos, iostat=nerr) rsignal(iff)
            end do
            ! - report stored signal values
            write (unit=MSP_stmp,fmt='(A,<nfil>G13.5,A)') &
              &   "- saved moment "//trim(adjustl(snumk))// &
              &   "("//trim(adjustl(snuml))//") signal (", &
              &   rsignal(1:nfil),")"
            call PostMessage(trim(MSP_stmp))
            ! - close files ...
            do iff=1, nfil
              close(unit=lfu(iff), iostat=nerr)
            end do
          end do ! i-loop over slices
          !
          ! END OF SINGLE PLANE FILE OUTPUT
          !
        end if ! SWITCH 3d output or single plane files 
        !
      end do ! loop l over components
      !
      ! update offset
      j = j + 1 + k
      !
    end do ! loop k over momentum orders
  end if
! ------------

! ------------

!  write(unit=*,fmt=*) " > ExportSTEMData: EXIT."
  return

END SUBROUTINE ExportSTEMData
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE InitProbeIntegration()
! function: Initializes all variables needed to handle probe intensity
!           integration and output
!           - determines the number of exit planes
!           - allocates the array holding the integrated intensities
!           - resets all variable handling the access to the arrays
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function after all multislice parameters are set up.
! Call this function if MS_pint_export > 0.
! Call this function before starting a new multislice calculation
!   for STEM mode for each scan pixel.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams
  use Plasmon

  implicit none

! ------------
! DECLARATION
  integer*4 :: nerr, nx, ny
  integer*4 :: nepw
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > InitWaveAvg: INIT."
  call PostMessage("Initializing probe intensity integration.")
  nepw = MSP_detpln
  nx = MS_dimx
  ny = MS_dimy
! ------------

! ------------
! Allocate the array holding the integrated probe intensities
  if (MSP_pimgmode/=0) then ! array is already allocated
    if (allocated(MSP_pimg)) deallocate(MSP_pimg,stat=nerr) ! deallocate
    allocate(MSP_pimg(1:nx, 1:ny, 0:nepw-1),stat=nerr) ! allocate
    if (nerr/=0) goto 101
    MSP_pimg = 0.0
  end if
  if (MSP_pdifmode/=0 .or. MSP_padifmode/=0) then ! array is already allocated
    if (allocated(MSP_pdif)) deallocate(MSP_pdif,stat=nerr) ! deallocate
    allocate(MSP_pdif(1:nx, 1:ny, 0:nepw-1),stat=nerr) ! allocate
    if (nerr/=0) goto 101
    MSP_pdif = 0.0
    if (MSP_padifmode/=0) then
      if (allocated(MSP_padif)) deallocate(MSP_padif,stat=nerr) ! deallocate
      allocate(MSP_padif(1:nx, 1:ny, 0:nepw-1),stat=nerr) ! allocate
      if (nerr/=0) goto 101
      MSP_padif = 0.0
      if (MS_wave_avg_export>0) then
        if (allocated(MSP_padif_ela)) deallocate(MSP_padif_ela,stat=nerr) ! deallocate
        allocate(MSP_padif_ela(1:nx, 1:ny, 0:nepw-1),stat=nerr) ! allocate
        if (nerr/=0) goto 101
        MSP_padif_ela = 0.0
      end if
    end if
  end if
  if (MSP_pimgmode/=0 .or. MSP_pdifmode/=0 .or. MSP_padifmode/=0) then ! arrays are allocated
    if (allocated(MSP_pint_nac)) deallocate(MSP_pint_nac,stat=nerr) ! deallocate
    MSP_pint_num = 0
    allocate(MSP_pint_nac(0:nepw-1, 0:PL_npemax),stat=nerr)
    if (nerr/=0) goto 101
    MSP_pint_nac = 0
    MSP_pint_num = nepw ! store the number of exit-planes
  end if
  MS_pint_idx = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > InitProbeIntegration: EXIT."
  return
  
! error handling
101 continue
  call CriticalError("InitProbeIntegration: Failed to allocate memory.")
  return

END SUBROUTINE InitProbeIntegration
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE ResetProbeIntegration()
! function: Resets all variables needed to handle probe intensity
!           integration and output
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function after all multislice parameters are set up.
! Call this function if MS_pint_export > 0.
! Call this function before starting a new multislice calculation
!   for STEM mode for each scan pixel.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ResetProbeIntegration: INIT."
! ------------

! ------------
! Reset the arrays holding the integrated probe intensities
  if (MSP_pimgmode/=0 .and. allocated(MSP_pimg)) then ! array is already allocated
    MSP_pimg = 0.0
  end if
  if ((MSP_pdifmode/=0 .or. MSP_padifmode/=0) .and. allocated(MSP_pdif)) then ! array is already allocated
    MSP_pdif = 0.0
    !! Do not reset MSP_padif and MSP_padif_ela here! These will accumulate over probe positions.
  end if
  if ((MSP_pimgmode/=0 .or. MSP_pdifmode/=0 .or. MSP_padifmode/=0) &
    & .and. allocated(MSP_pint_nac)) then ! array is already allocated
    MSP_pint_nac = 0
  end if
  MS_pint_idx = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > ResetProbeIntegration: EXIT."
  return

END SUBROUTINE ResetProbeIntegration
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE UnInitProbeIntegration()
! function: Uninitializes all variables needed to handle probe
!           intensity integrations
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function at the end of the program or a multislice run.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  integer*4 :: nerr
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > UnInitProbeIntegration: INIT."
  call PostMessage("Uninitializing probe intensity integration.")
! ------------

! ------------
! Deallocate the array holding the average wavefunctions
  if (allocated(MSP_pimg)) deallocate(MSP_pimg, stat=nerr)
  if (allocated(MSP_pdif)) deallocate(MSP_pdif, stat=nerr)
  if (allocated(MSP_padif)) deallocate(MSP_padif, stat=nerr)
  if (allocated(MSP_padif_ela)) deallocate(MSP_padif_ela, stat=nerr)
  if (allocated(MSP_pint_nac)) deallocate(MSP_pint_nac, stat=nerr)
! reset the access and accumulation indices
  MSP_pint_num = 0
  MS_pint_idx = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > UnInitProbeIntegration: EXIT."
  return

END SUBROUTINE UnInitProbeIntegration
!**********************************************************************!
  
!**********************************************************************!
! Warning: No controls are done here to check whether the file is OK.
SUBROUTINE AddImageToOpenFile(lun, lunpos, img, imgbuf, nx, ny, nerr)
  use IFPORT
  integer*4, intent(in) :: lun, lunpos, nx, ny
  real*4, intent(in) :: img(1:nx,1:ny)
  real*4, intent(inout) :: imgbuf(1:nx,1:ny)
  integer*4, intent(inout) :: nerr
  integer*4 :: ioerr
  ioerr = 0
  nerr = 0
  ! read old data from file
  ioerr = fseek(lun, lunpos, 0)
  if (ioerr/=0) then
    nerr = 1
    return
  end if
  read(unit=lun,iostat=ioerr) imgbuf(1:nx,1:ny)
  if (ioerr/=0) then
    nerr = 2
    return
  end if
  ! add the new data to the old data
  imgbuf(1:nx,1:ny) = imgbuf(1:nx,1:ny) + img(1:nx,1:ny)
  ! write the data for slice islc in slot k
  ioerr = fseek(lun, lunpos, 0)
  if (ioerr/=0) then
    nerr = 1
    return
  end if
  write(unit=lun,iostat=ioerr) imgbuf(1:nx,1:ny)
  if (nerr/=0) then
    nerr = 3
    return
  end if
  return
end SUBROUTINE AddImageToOpenFile
!**********************************************************************!
  
  
!**********************************************************************!
!**********************************************************************!
SUBROUTINE ExportProbeAvgDif(sfile)
! function: Exports the averaged diffraction pattern accumulated
!           in MSP_padif* to files.
!           - loops through all exit planes used for recording
!
!           This function will open with exclusive access since it
!           may need to read existing data before writing new data.
!           This will allow multiple-processes to write to the same
!           file. A repeat and waiting scheme is implemented to
!           queue up file acces from multiple processes.
!

  use IFPORT
  use MultiSlice
  use MSAparams

  implicit none
  
! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  character(len=MSP_ll) :: isfile, sexpfile(3)
  integer*4 :: nerr, nexist, nela, nslc, lun, islc, iout, k
  integer*4 :: nx, ny, nslotbytes
  real*4, allocatable :: ptmp(:,:), patmp(:,:), patds(:,:)
  external :: fileopenexclrw ! (sfile, lun, nexist, nerr)
  external :: createfilefolder ! (sfile,nerr) this file
  external :: sinsertslcidx ! (idx,idxlen,sfnin,sfnadd,sfnext,sfnout) this file
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ExportProbeAvgDif: INIT."
  nerr = 0
  if (MSP_padifmode==0 .or. (.not.allocated(MSP_padif))) goto 99
  nx = MS_dimx
  ny = MS_dimy
  if (nx <= 0 .or. ny <= 0) goto 99 ! no valid setup
  nslotbytes = nx*ny*4 ! set number of bytes per slot
  iout = 0
  lun = 0
  nslc = MS_stacksize ! max. number of object slices
  nela = 0 ! init without elastic pattern data (0)
  if (MS_wave_avg_export>0 .and. allocated(MSP_padif)) nela = 1 ! there should be elastic data -> trigger tds output
  if (LEN_TRIM(sfile)==0) then ! set a default output file name in case of invalid input
    isfile = "probe.dat"
  else ! use the input file name
    isfile = trim(sfile)
  end if
  allocate(ptmp(1:nx,1:ny), patmp(1:nx,1:ny), stat=nerr)
  if (nerr/=0) then
    call CriticalError("Failed to allocate memory for storing <padif> data.")
    goto 99
  end if
  ptmp = 0
  patmp = 0
  if (nela>0) then
    allocate(patds(1:nx,1:ny), stat=nerr)
    if (nerr/=0) then
      call CriticalError("Failed to allocate memory for storing <padif_tds> data.")
      goto 99
    end if
    patds = 0
  end if
! ------------

! ------------
  if (MSP_3dout > 0) then ! full stack output (open file only once)
    ! prepare file names
    call sinsertslcidx(0,0,trim(isfile),"_padif_tot",".dat",sexpfile(1))
    call sinsertslcidx(0,0,trim(isfile),"_padif_ela",".dat",sexpfile(2))
    call sinsertslcidx(0,0,trim(isfile),"_padif_tds",".dat",sexpfile(3))
    if (iout==0) call createfilefolder(trim(sexpfile(1)),nerr) ! try to create the folder
    iout = 1
    !
    ! * TOTAL INTENSITY DATA 
    !
    ! file opening
    call fileopenexclrw(sexpfile(1), lun, nexist, nerr)
    if (nerr/=0) goto 100
    ! work on the open file
    if (nexist>0) then ! existing file -> add the data
      do islc=0, nslc ! Loop over all slices
        k = MSP_ldetpln(islc) ! get storage slot
        if (k < 0) cycle ! nothing stored for islc
        ! add data to slot k in the file
        ptmp(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k), kind=4)
        call AddImageToOpenFile(lun, k*nslotbytes, &
           & ptmp(1:nx,1:ny), patmp(1:nx,1:ny), nx, ny , nerr)
        if (nerr/=0) goto 101
      end do
    else ! new file -> just store all data
      do islc=0, nslc ! Loop over all slices
        k = MSP_ldetpln(islc) ! get storage slot
        if (k < 0) cycle ! nothing stored for islc
        ! write the data for slice islc in slot k
        ptmp(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k), kind=4)
        write(unit=lun,iostat=nerr) ptmp(1:nx,1:ny)
        if (nerr/=0) then
          nerr = 3
          goto 101
        end if
      end do
    end if
    close(unit=lun) ! close the file
    !
    if (nela>0) then ! elastic and tds data is present
      !
      ! * ELASTIC INTENSITY DATA 
      !
      ! file opening
      call fileopenexclrw(sexpfile(2), lun, nexist, nerr)
      if (nerr/=0) goto 100
      ! work on the open file
      if (nexist>0) then ! existing file -> add the data
        do islc=0, nslc ! Loop over all slices
          k = MSP_ldetpln(islc) ! get storage slot
          if (k < 0) cycle ! nothing stored for islc
          ! add data to slot k in the file
          ptmp(1:nx,1:ny) = real(MSP_padif_ela(1:nx,1:ny,k), kind=4)
          call AddImageToOpenFile(lun, k*nslotbytes, &
             & ptmp(1:nx,1:ny), patmp(1:nx,1:ny), nx, ny , nerr)
          if (nerr/=0) goto 101
        end do
      else ! new file -> just store all data
        do islc=0, nslc ! Loop over all slices
          k = MSP_ldetpln(islc) ! get storage slot
          if (k < 0) cycle ! nothing stored for islc
          ! write the data for slice islc in slot k
          ptmp(1:nx,1:ny) = real(MSP_padif_ela(1:nx,1:ny,k), kind=4)
          write(unit=lun,iostat=nerr) ptmp(1:nx,1:ny)
          if (nerr/=0) then
            nerr = 3
            goto 101
          end if
        end do
      end if
      close(unit=lun) ! close the file
      !
      ! * TDS INTENSITY DATA 
      !
      ! file opening
      call fileopenexclrw(sexpfile(3), lun, nexist, nerr)
      if (nerr/=0) goto 100
      ! work on the open file
      if (nexist>0) then ! existing file -> add the data
        do islc=0, nslc ! Loop over all slices
          k = MSP_ldetpln(islc) ! get storage slot
          if (k < 0) cycle ! nothing stored for islc
          ! get the tds data
          patds(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k) - MSP_padif_ela(1:nx,1:ny,k), kind=4)
          ! add data to slot k in the file
          call AddImageToOpenFile(lun, k*nslotbytes, &
             & patds(1:nx,1:ny), patmp(1:nx,1:ny), nx, ny , nerr)
          if (nerr/=0) goto 101
        end do
      else ! new file -> just store all data
        do islc=0, nslc ! Loop over all slices
          k = MSP_ldetpln(islc) ! get storage slot
          if (k < 0) cycle ! nothing stored for islc
          ! get the tds data
          patds(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k) - MSP_padif_ela(1:nx,1:ny,k), kind=4)
          ! write the data for slice islc in slot k
          write(unit=lun,iostat=nerr) patds(1:nx,1:ny)
          if (nerr/=0) then
            nerr = 3
            goto 101
          end if
        end do
      end if
      close(unit=lun) ! close the file
    end if
    !
  else ! single slice output (need to open many files)
    do islc=0, nslc ! Loop over all slices
      k = MSP_ldetpln(islc) ! get storage slot
      if (k < 0) cycle ! nothing stored for islc
      ! prepare file names
      call sinsertslcidx(islc,MS_nslid,trim(isfile),"_padif_tot",".dat",sexpfile(1))
      call sinsertslcidx(islc,MS_nslid,trim(isfile),"_padif_ela",".dat",sexpfile(2))
      call sinsertslcidx(islc,MS_nslid,trim(isfile),"_padif_tds",".dat",sexpfile(3))
      if (iout==0) call createfilefolder(trim(sexpfile(1)),nerr) ! try to create the folder
      iout = iout + 1
      !
      ! * TOTAL INTENSITY DATA 
      !
      ! open the file 
      call fileopenexclrw(sexpfile(1), lun, nexist, nerr)
      if (nerr/=0) goto 100
      ! work on the open file
      if (nexist>0) then ! existing file -> add the data
        ! add data to the file
        ptmp(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k), kind=4)
        call AddImageToOpenFile(lun, 0,ptmp(1:nx,1:ny), &
          & patmp(1:nx,1:ny), nx, ny , nerr)
        if (nerr/=0) goto 101
      else ! new file -> just store all data
        ! write the data
        ptmp(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k), kind=4)
        write(unit=lun,iostat=nerr) ptmp(1:nx,1:ny)
        if (nerr/=0) then
          nerr = 3
          goto 101
        end if
      end if
      close(unit=lun) ! close the file
      !
      if (nela>0) then ! elastic and tds data
        !
        ! * ELASTIC INTENSITY DATA 
        !
        ! open the file 
        call fileopenexclrw(sexpfile(2), lun, nexist, nerr)
        if (nerr/=0) goto 100
        ! work on the open file
        if (nexist>0) then ! existing file -> add the data
          ! add data the file
          ptmp(1:nx,1:ny) = real(MSP_padif_ela(1:nx,1:ny,k), kind=4)
          call AddImageToOpenFile(lun, 0, ptmp(1:nx,1:ny), &
            & patmp(1:nx,1:ny), nx, ny , nerr)
          if (nerr/=0) goto 101
        else ! new file -> just store all data
          ! write the data
          ptmp(1:nx,1:ny) = real(MSP_padif_ela(1:nx,1:ny,k), kind=4)
          write(unit=lun,iostat=nerr) ptmp(1:nx,1:ny)
          if (nerr/=0) then
            nerr = 3
            goto 101
          end if
        end if
        close(unit=lun) ! close the file
        !
        ! * TDS INTENSITY DATA 
        !
        ! get the tds data
        patds(1:nx,1:ny) = real(MSP_padif(1:nx,1:ny,k) - MSP_padif_ela(1:nx,1:ny,k), kind=4)
        ! open the file 
        call fileopenexclrw(sexpfile(1), lun, nexist, nerr)
        if (nerr/=0) goto 100
        ! work on the open file
        if (nexist>0) then ! existing file -> add the data
          ! add data to slot k in the file
          call AddImageToOpenFile(lun, 0, patds(1:nx,1:ny), &
            & patmp(1:nx,1:ny), nx, ny , nerr)
          if (nerr/=0) goto 101
        else ! new file -> just store all data
          ! write the data for slice islc in slot k
          write(unit=lun,iostat=nerr) patds(1:nx,1:ny)
          if (nerr/=0) then
            nerr = 3
            goto 101
          end if
        end if
        close(unit=lun) ! close the file
        !
      end if
      !
    end do ! islc=0, nslc ! Loop over all slices
  end if
! ------------
  
! ------------
99 continue
  if (allocated(patmp)) deallocate(patmp,stat=nerr)
  if (allocated(patds)) deallocate(patds,stat=nerr)
  return
100 continue
  close(unit=lun) ! close the file
  call CriticalError("Failed to open file for <padif> output.")
  goto 99
101 continue
  close(unit=lun) ! close the file
  if (nerr==1) call CriticalError("Failed to position in file for <padif> output.")
  if (nerr==2) call CriticalError("Failed to read old data from file for adding <padif> data.")
  if (nerr==3) call CriticalError("Failed to write <padif> data to output file.")
  goto 99

  return
  
END SUBROUTINE ExportProbeAvgDif
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE ExportProbeIntensity(sfile)
! function: Exports the accumulated probe intensity to files.
!           - loops through all exit planes used for recording
!           - normalizes the accumulated intensities
!           - stores the normalized intensities to disk
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function if MS_pint_export>0.
! Call this function after a multislice calculation
!   for STEM mode only, for each scan pixel.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  character(len=MSP_ll) :: isfile, sexpfile(3)
  integer*4 :: nintout, nwavavg, ntransform, nuidx
  integer*4 :: nx, ny, nerr, nalloc, i, j, k
  integer*4 :: islc, iout, nslc
  real(fpp) :: rnorm, pint, rsca, rscas
  real(fpp), dimension(:,:), allocatable :: pimg, pela, ptds
  complex(fpp), dimension(:,:), allocatable :: wave !, work
  external :: SaveDataC, SaveDataR ! (sfile,dat,n,nerr) this file
  external :: AppendDataC, AppendDataR ! (sfile,dat,n,nerr) this file
  !external :: AddDataR4 ! (sfile,dat,n,nerr) this file
  external :: sinsertslcidx ! (idx,idxlen,sfnin,sfnadd,sfnext,sfnout) this file
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ExportProbeIntensity: INIT."
  nwavavg = 0 ! init without average wave function data: 0 -> none, 1 -> Fourier space, 2 -> real space
  nintout = 0 ! intensity output request strength: 0 -> none, 1 -> Fourier space, 2 -> real space, 3 -> both
  ntransform = 0 ! init without the need to transform data: 0 -> no, 1 -> yes
  nx = MS_dimx
  ny = MS_dimy
  if (nx <= 0 .or. ny <= 0) return ! no valid setup
  if (MS_pint_export<=0) return ! no valid setup
  if (MSP_pdifmode>0 .or. MSP_padifmode>0) nintout = 1
  if (MSP_pimgmode>0) nintout = nintout + 2
  rsca = 1.0_fpp / real(nx*ny,kind=fpp) ! for DFT renormalizations
  rscas = 1.0_fpp / real(MSP_SF_ndimx*MSP_SF_ndimy, kind=fpp) ! for average diffraction pattern normalization
  !
  ! Handle the case of present average wave function data.
  ! In this case, we want to export the elastic images as well as 
  ! also the difference (TDS) images to the total intensity.
  ! We then need extra data arrays
  if (MS_wave_avg_export>0 .and. allocated(MS_wave_avg)) then ! there is average wave data
    nwavavg = 2 - MS_wave_export_form ! 1 -> Fourier space, 2 -> real space wave function
    if ( nwavavg /= nintout ) then
      ! this means, average wave and at least one intended output are not in the same space
      ! and a fourier transform has to be done.
      ntransform = 1
    end if
  end if
  !
  if (LEN_TRIM(sfile)==0) then ! set a default output file name in case of invalid input
    isfile = "probe.dat"
  else ! use the input file name
    isfile = trim(sfile)
  end if
  !
  nuidx = 1
  if (MSP_3dout > 0) nuidx = 0 ! switch off using the index for file names
  nslc = MS_stacksize ! max. number of object slices
! ------------

! ------------
! OUTPUT OF PROBE INTENSITIES IN REAL SPACE
  if (MSP_pimgmode==1) then
    ! allocations
    allocate(pimg(nx,ny), stat=nalloc)
    if (nwavavg>0) then
      allocate(wave(nx,ny), stat=nalloc)
      allocate(pela(nx,ny), stat=nalloc)
      allocate(ptds(nx,ny), stat=nalloc)
    end if
    iout = 0 ! initialize output plane counter
    do islc=0, nslc ! Loop over all slices
      k = MSP_ldetpln(islc) ! get storage slot
      if (k < 0) cycle ! nothing stored for islc
      ! normalize (we assume that wave and images have the same number of contributions)
      if (MSP_pint_nac(k,0)>0) then
        rnorm = 1.0_fpp/real(MSP_pint_nac(k,0),kind=fpp)
      else
        rnorm = 1.0_fpp
      end if
      !
      ! get total intensity
      pimg(1:nx,1:ny) = real(MSP_pimg(1:nx,1:ny,k), kind=fpp) * rnorm
      ! prepare file names
      call sinsertslcidx(nuidx*islc,nuidx*MS_nslid,trim(isfile),"_pimg_tot",".dat",sexpfile(1))
      call sinsertslcidx(nuidx*islc,nuidx*MS_nslid,trim(isfile),"_pimg_ela",".dat",sexpfile(2))
      call sinsertslcidx(nuidx*islc,nuidx*MS_nslid,trim(isfile),"_pimg_tds",".dat",sexpfile(3))
      if (MSP_3dout > 0 .and. iout > 0) then ! /3dout append
        call AppendDataR(trim(sexpfile(1)), pimg, nx*ny, nerr) ! append tot
      else ! writing to new new file 
        call PostMessage("  Writing probe image total intensity to file ["//trim(sexpfile(1))//"].")
        call SaveDataR(trim(sexpfile(1)), pimg, nx*ny, nerr) ! save tot
      end if
      ! 
      if (nwavavg>0) then
        ! get elastic and tds images
        if (nwavavg==1) then ! wave data is in Fourier space, need to transform
          MS_work(1:nx,1:ny) = cmplx(MS_wave_avg(1:nx,1:ny,k), kind=fpp) * rnorm
          ! call MS_FFT(work,MS_dimx,MS_dimy,'backwards')
          call MS_FFT_WORK(-1)
          wave(1:nx,1:ny) = MS_work(1:nx,1:ny) * sqrt(rsca) ! renormalize after iDFT
        else ! wave data is in real space, just copy
          wave(1:nx,1:ny) = cmplx(MS_wave_avg(1:nx,1:ny,k), kind=fpp) * rnorm
        end if
        ! calculate elastic image
        do j=1, ny
          do i=1, nx
            pint = real( wave(i,j)*conjg(wave(i,j)) )
            pela(i,j) = pint
          end do
        end do
        ptds = pimg - pela
        if (MSP_3dout > 0 .and. iout > 0) then ! /3dout -> single files append
          call AppendDataR(trim(sexpfile(2)), pela, nx*ny, nerr) ! append ela
          call AppendDataR(trim(sexpfile(3)), ptds, nx*ny, nerr) ! append tds
        else ! individual files per plane
          call PostMessage("  Writing probe image elastic intensity to file ["//trim(sexpfile(2))//"].")
          call SaveDataR(trim(sexpfile(2)), pela, nx*ny, nerr) ! save ela
          call PostMessage("  Writing probe image TDS intensity to file ["//trim(sexpfile(3))//"].")
          call SaveDataR(trim(sexpfile(3)), ptds, nx*ny, nerr) ! save tds
        end if
      end if
      !
      iout = iout + 1 ! increment output plane counter
      !
    end do ! k=0, MS_pint_num
    !
    ! deallocate
    if (allocated(pimg)) deallocate(pimg, stat=nalloc)
    if (allocated(wave)) deallocate(wave, stat=nalloc)
    if (allocated(pela)) deallocate(pela, stat=nalloc)
    if (allocated(ptds)) deallocate(ptds, stat=nalloc)
    !
  end if ! (MSP_pimgmode==1)
! ------------


! ------------
! OUTPUT OF PROBE INTENSITIES IN FOURIER SPACE
  if (MSP_pdifmode==1 .or. MSP_padifmode==1) then
    ! allocations
    allocate(pimg(nx,ny), stat=nalloc)
    if (nwavavg>0) then
      allocate(wave(nx,ny), stat=nalloc)
      allocate(pela(nx,ny), stat=nalloc)
      allocate(ptds(nx,ny), stat=nalloc)
    end if
    iout = 0
    do islc=0, nslc ! Loop over all slices
      k = MSP_ldetpln(islc) ! get storage slot
      if (k < 0) cycle ! nothing stored for islc
      ! normalize (we assume that wave and images have the same number of contributions)
      if (MSP_pint_nac(k,0)>0) then
        rnorm = 1.0_fpp/real(MSP_pint_nac(k,0),kind=fpp)
      else
        rnorm = 1.0_fpp
      end if
      !
      ! get total intensity
      pimg(1:nx,1:ny) = real(MSP_pdif(1:nx,1:ny,k), kind=fpp) * rnorm
      !
      if (MSP_padifmode==1) then ! add to average diffraction pattern
        MSP_padif(1:nx,1:ny,k) = MSP_padif(1:nx,1:ny,k) + MSP_pdif(1:nx,1:ny,k)*rscas ! ... normalized to number of scan points
      end if
      if (MSP_pdifmode==1) then ! ouput of diffraction pattern per scan position
        ! prepare file names and save
        call sinsertslcidx(nuidx*islc,nuidx*MS_nslid,trim(isfile),"_pdif_tot",".dat",sexpfile(1))
        call sinsertslcidx(nuidx*islc,nuidx*MS_nslid,trim(isfile),"_pdif_ela",".dat",sexpfile(2))
        call sinsertslcidx(nuidx*islc,nuidx*MS_nslid,trim(isfile),"_pdif_tds",".dat",sexpfile(3))
        if (MSP_3dout > 0 .and. iout >0) then ! /3dout append
          call AppendDataR(trim(sexpfile(1)), pimg, nx*ny, nerr) ! append to old file
        else ! single file per plane
          call PostMessage("  Writing probe diffraction total intensity to file ["//trim(sexpfile(1))//"].")
          call SaveDataR(trim(sexpfile(1)), pimg, nx*ny, nerr) ! save to new file
        end if
      end if
      ! 
      if (nwavavg>0) then
        ! get elastic and tds images
        if (nwavavg==2) then ! wave data is in real space, need to transform
          MS_work(1:nx,1:ny) = cmplx(MS_wave_avg(1:nx,1:ny,k), kind=fpp) * rnorm
          ! call MS_FFT(work,MS_dimx,MS_dimy,'forwards')
          call MS_FFT_WORK(1)
          wave(1:nx,1:ny) = MS_work(1:nx,1:ny) * sqrt(rsca) ! renormalize after DFT
        else ! wave data is in Fourier space, just copy
          wave(1:nx,1:ny) = cmplx(MS_wave_avg(1:nx,1:ny,k), kind=fpp) * rnorm
        end if
        ! calculate elastic image
        do j=1, ny
          do i=1, nx
            pint = real( wave(i,j)*conjg(wave(i,j)) )
            pela(i,j) = pint
          end do
        end do
        ptds = pimg - pela
        if (MSP_padifmode==1) then ! add elastic diffraction pattern to average elastic diffraction pattern
          MSP_padif_ela(1:nx,1:ny,k) = MSP_padif_ela(1:nx,1:ny,k) + real(pela(1:nx,1:ny)*rscas, kind=fpp_ac) ! ... normalized to number of scan points
        end if
        if (MSP_pdifmode==1) then ! output of elastic diffraction pattern per scan position
          if (MSP_3dout > 0 .and. iout > 0) then ! /3dout -> single files append
            call AppendDataR(trim(sexpfile(2)), pela, nx*ny, nerr) ! append ela
            call AppendDataR(trim(sexpfile(3)), ptds, nx*ny, nerr) ! append tds
          else ! individual files per plane
            call PostMessage("  Writing probe diffraction elastic intensity to file ["//trim(sexpfile(2))//"].")
            call SaveDataR(trim(sexpfile(2)), pela, nx*ny, nerr) ! save ela
            call PostMessage("  Writing probe diffraction TDS intensity to file ["//trim(sexpfile(3))//"].")
            call SaveDataR(trim(sexpfile(3)), ptds, nx*ny, nerr) ! save tds
          end if
        end if
      end if
      !
      iout = iout + 1 ! increment probe diffraction plane output counter
      !
    end do ! k=0, MS_pint_num
    !
    ! deallocate
    if (allocated(pimg)) deallocate(pimg, stat=nalloc)
    if (allocated(wave)) deallocate(wave, stat=nalloc)
    if (allocated(pela)) deallocate(pela, stat=nalloc)
    if (allocated(ptds)) deallocate(ptds, stat=nalloc)
    !
  end if ! (MSP_pdifmode==1 .or. MSP_padifmode==1)
! ------------

! ------------
  if (allocated(wave)) deallocate(wave, stat=nerr)
!  write(unit=*,fmt=*) " > ExportProbeIntensity: EXIT."
  return

END SUBROUTINE ExportProbeIntensity
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE InitWaveAvg()
! function: Initializes all variables needed to handle avg. wavefunct.
!           - determines the number of exit planes
!           - allocates the array holding the avg. wavefunctions
!           - resets the avg. wavefunctions
!           - resets all variable handling the access to the avg. wf
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function after all multislice parameters are set up.
! Call this function if MS_wave_avg_export>0.
! Call this function before starting a new multislice calculation
!   for CTEM mode and for STEM mode.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  integer*4 :: nerr, nx, ny
  integer*4 :: nepw
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > InitWaveAvg: INIT."
  call PostMessage("Initializing wavefunction averaging.")
  nepw = MSP_detpln
  nx = MS_dimx
  ny = MS_dimy
! ------------

! ------------
! Allocate the array holding the average wavefunctions
  if (allocated(MS_wave_avg)) deallocate(MS_wave_avg,stat=nerr) ! deallocate on wrong size
  if (allocated(MS_wave_avg_nac)) deallocate(MS_wave_avg_nac,stat=nerr) ! 
  MS_wave_avg_num = 0
  ! alloation with an extra array for the incoming wavefunction, index 0 at dimension 3
  allocate(MS_wave_avg(1:nx, 1:ny, 0:nepw-1),stat=nerr)
  allocate(MS_wave_avg_nac(0:nepw-1),stat=nerr)
  if (nerr/=0) then
    call CriticalError("InitWaveAvg: Failed to allocate memory.")
  end if
  MS_wave_avg_num = nepw ! store the number of exit-planes
! reset the wavefunctions
  MS_wave_avg = cmplx(0.,0.,kind=fpp_ac)
  MS_wave_avg_nac = 0
! ------------

! ------------
! reset the average wf access and accumulation indices
  MS_wave_avg_idx = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > InitWaveAvg: EXIT."
  return

END SUBROUTINE InitWaveAvg
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE ResetWaveAvg()
! function: Resets all variables needed to handle avg. wavefunct.
!           - resets the avg. wavefunctions
!           - resets all variable handling the access to the avg. wf
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function if MS_wave_avg_export>0.
! Call this function before starting a new multislice calculation
!   in STEM mode for each pixel.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams
  use precision

  implicit none

! ------------
! DECLARATION
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ResetWaveAvg: INIT."
! ------------

! ------------
  if ((.not.allocated(MS_wave_avg)).or.MS_wave_avg_num<1) then ! array not allocated
    ! do nothing
    return
  end if
! reset the wavefunction data
  MS_wave_avg = cmplx(0.,0., kind=fpp_ac)
  MS_wave_avg_nac = 0
  MS_wave_avg_idx = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > ResetWaveAvg: EXIT."
  return

END SUBROUTINE ResetWaveAvg
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE UnInitWaveAvg()
! function: Uninitializes all variables needed to handle avg. wavefunct.
!           - deallocates the array holding the avg. wavefunctions
!           - resets all variable handling the access to the avg. wf
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function at the end of the program or a multislice run.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  integer*4 :: nerr
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > UnInitWaveAvg: INIT."
  call PostMessage("Uninitializing wavefunction averaging.")
! ------------

! ------------
! Deallocate the array holding the average wavefunctions
  if (allocated(MS_wave_avg)) deallocate(MS_wave_avg,stat=nerr) ! deallocate on wrong size
  if (allocated(MS_wave_avg_nac)) deallocate(MS_wave_avg_nac,stat=nerr) ! 
! reset the average wf access and accumulation indices
  MS_wave_avg_num = 0
  MS_wave_avg_idx = 0
! ------------

! ------------
!  write(unit=*,fmt=*) " > UnInitWaveAvg: EXIT."
  return

END SUBROUTINE UnInitWaveAvg
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE ExportWaveAvg(sfile)
! function: Exports the accumulated avg. wavefunctions to files.
!           - loops through all exit planes used for recording
!           - normalizes the accumulated wavefunctions
!           - stores the normalized wavefunctions to disk
! -------------------------------------------------------------------- !
! parameters: none
! -------------------------------------------------------------------- !
! remarks:
! Call this function if MS_wave_avg_export>0.
! Call this function after a multislice calculation
!   for CTEM mode (once) and for STEM mode for each scan pixel.
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  character(len=MSP_ll) :: isfile, sexpfile
  integer*4 :: nerr, k, nx, ny, nuidx
  integer*4 :: nslc, islc, iout
  real(fpp) :: rnorm
  complex(fpp), dimension(:,:), allocatable :: wave
  external :: SaveDataC
  !external :: AppendDataC
  external :: sinsertslcidx ! (idx,idxlen,sfnin,sfnadd,sfnext,sfnout)
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ExportWaveAvg: INIT."
  if (MS_wave_avg_export /= 1) return ! only export for MS_wave_avg_export == 1
  if (.not.allocated(MS_wave_avg)) return ! no data
  !
  if (MS_wave_export_form==0) then
    call PostMessage("Average wavefunction output (real space).")
  else
    call PostMessage("Average wavefunction output (Fourier space).")
  end if
  nx = MS_dimx
  ny = MS_dimy
  allocate(wave(1:nx,1:ny), stat=nerr)
  if (nerr/=0) then
    call CriticalError("ExportWaveAvg: Failed to allocate memory.")
    return
  end if
  if (LEN_TRIM(sfile)==0) then ! set a default output file name in case of invalid input
    isfile = "epw_avg.wav"
  else ! use the input file name
    isfile = trim(sfile)
  end if
  nuidx = 1
  if (MSP_3dout) nuidx = 0 ! slice index usage flag in file names
  nslc = MS_stacksize
! ------------

! ------------
! OUTPUT OF AVERAGE WAVE FUNCTION
  iout=0
  do islc = 0, nslc ! Loop over all sample slices
    k = MSP_ldetpln(islc) ! get storage slot
    if (k < 0) cycle ! skip, nothing stored for islc
    ! normalization factor
    if (MS_wave_avg_nac(k)>0) then
      rnorm = 1.0_fpp/real(MS_wave_avg_nac(k), kind=fpp)
    else
      rnorm = 1.0_fpp
    end if
    wave(1:nx, 1:ny) = cmplx(MS_wave_avg(1:nx, 1:ny, k), kind=fpp) * rnorm
    call sinsertslcidx(islc*nuidx,MS_nslid*nuidx,trim(isfile),"_avg",".wav",sexpfile)
    ! export to file
    if (MSP_3dout>0 .and. iout >0) then ! /3dout append
      call AppendDataC(trim(sexpfile), wave, nx*ny, nerr) ! append
    else ! first or single slice plane output
      call PostMessage("  Writing average wave function to file ["//trim(sexpfile)//"].")
      call SaveDataC(trim(sexpfile), wave, nx*ny, nerr) ! save
    end if
    iout = iout + 1
    ! ...
  end do ! k=0, MS_wave_avg_num
! ------------

! ------------
  if (allocated(wave)) deallocate(wave, stat=nerr)
!  write(unit=*,fmt=*) " > ExportWaveAvg: EXIT."
  return

END SUBROUTINE ExportWaveAvg
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE ExportWave(sfile, islice)
! function: saves the wave function to file
!           - expects the wave function in MS_wave in Fourier form
!           - does inverse FT
!           - writes the binary data to a file with given file name
!           - accumulates probe images, diffraction patterns and
!             average wave functions
! -------------------------------------------------------------------- !
! parameter: character(len=*) :: sfile    ! the basic output file name
!            integer*4 :: islice ! slice index for extra file naming
!                                ! 0 = incident probe plane
!                                ! >0 = planes in the sample
! -------------------------------------------------------------------- !

  use MultiSlice
  use MSAparams
  use precision

  implicit none

! ------------
! DECLARATION
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: islice
  integer*4 :: nerr, i, j, nx, ny, nwavrs, nuidx, iwav, iimg
  real(fpp) :: pint, rsca
  complex(fpp), allocatable :: wave(:,:) !(MS_dimx,MS_dimy)
  character(len=1024) :: sexpfile
  external :: sinsertslcidx
  external :: SaveDataC, AppendDataC
! ------------

! ------------
! INIT
!  write(unit=*,fmt=*) " > ExportWave: INIT."
  nerr = 0
  nwavrs = 0
  nx = MS_dimx
  ny = MS_dimy
  rsca = 1.0_fpp / real(nx*ny, kind=fpp) ! /!\ we do not check for nx*ny > 0
  allocate(wave(nx,ny), stat=nerr)
  if (nerr/=0) then
    call CriticalError("ExportWave: Failed to allocate memory.")
  end if
  nuidx = 1
  if (MSP_3dout>0) nuidx = 0 ! /3dout -> add no index suffix to file name
  ! create output file name for the wave function
  call sinsertslcidx(islice*nuidx,MS_nslid*nuidx,trim(sfile),"",".wav",sexpfile)
  iwav = MS_wave_avg_idx
  iimg = MS_pint_idx
! ------------

! ------------
! wave function export and averaging
  if (MS_wave_export>0 .or. MS_wave_avg_export>0) then
    ! from calculation frame (MS_wave) to local frame (wave)
    if (MS_wave_export_form==1) then
      ! Fourier space export (direct copy)
      wave(1:nx,1:ny) = MS_wave(1:nx,1:ny)
    else
      ! real space export (inverse FT)
      ! - transfer data
      !   for an unknown reason, the following line causes a stack overflow and access violation
      !   MS_work(1:nx,1:ny) = MS_wave(1:nx,1:ny)
      !   Though, the explicit assignement below element by element works.
      do j=1, ny
        do i=1, nx
          MS_work(i,j) = MS_wave(i,j)
        end do
      end do
      ! call MS_FFT(work,MS_dimx,MS_dimy,'backwards')
      call MS_FFT_WORK(-1)
      wave(1:nx,1:ny) = MS_work(1:nx,1:ny) * sqrt(rsca) ! renormalize after iDFT
      nwavrs = 1
    end if
    if (MS_wave_export>0) then ! individual wave export to disk
      
      if (MSP_3dout > 0 .and. iwav > 0) then ! 3d and not the first slot
        call AppendDataC(trim(sexpfile), wave, nx*ny, nerr) ! append file
      else ! 2d or first 3d slot
        call PostMessage("  Writing wave function to file ["// &
          &  trim(sexpfile)//"].")
        call SaveDataC(trim(sexpfile), wave, nx*ny, nerr) ! store new file
      end if
    end if
    if (MS_wave_avg_export>0) then ! accumulation of the elastic wave
      MS_wave_avg(1:nx,1:ny,iwav) = MS_wave_avg(1:nx,1:ny,iwav) + cmplx(wave(1:nx,1:ny), kind=fpp_ac)
      MS_wave_avg_nac(iwav) = MS_wave_avg_nac(iwav) + 1
    end if
  end if
! ------------

! ------------
! probe intensity accumulations
  if (MS_pint_export>0) then ! accumulation of probe intensities
    if (MSP_pimgmode>0) then ! real-space intensity accumulation
      if (nwavrs==0) then ! no real-space wave function generated here
        ! transform current wave function (MS_wave) to real space
        ! - transfer data
        !   for an unknown reason, the following line causes an stack overflow and access violation
        !   MS_work(1:nx,1:ny) = MS_wave(1:nx,1:ny)
        !   Though, the explicit assignement below element by element works.
        do j=1, ny
          do i=1, nx
            MS_work(i,j) = MS_wave(i,j)
          end do
        end do
        ! call MS_FFT(work,MS_dimx,MS_dimy,'backwards')
        call MS_FFT_WORK(-1)
        ! accumulate probe image from MS_work
        do j=1, MS_dimy
          do i=1, MS_dimx
            pint = real( MS_work(i,j)*conjg(MS_work(i,j)) ) * rsca ! renormalize after iDFT 
            MSP_pimg(i,j,iimg) = MSP_pimg(i,j,iimg) + real(pint, kind=fpp_ac)
          end do
        end do
      else ! real-space wave function was already created above
        ! accumulate probe image from wave
        do j=1, MS_dimy
          do i=1, MS_dimx
            pint = real( wave(i,j)*conjg(wave(i,j)) ) ! wave should already be normalized in this case
            MSP_pimg(i,j,iimg) = MSP_pimg(i,j,iimg) + real(pint, kind=fpp_ac)
          end do
        end do
      end if
    end if
    if (MSP_pdifmode>0 .or. MSP_padifmode>0) then ! Fourier-space intensity accumulation
      ! accumulate probe diffraction pattern from MS_wave
      do j=1, MS_dimy
        do i=1, MS_dimx
          pint = real( MS_wave(i,j)*conjg(MS_wave(i,j)) )
          MSP_pdif(i,j,iimg) = MSP_pdif(i,j,iimg) + real(pint, kind=fpp_ac)
        end do
      end do
    end if
    MSP_pint_nac(iimg,0) = MSP_pint_nac(iimg,0) + 1
  end if
! ------------

! ------------
! dealloc
  if (allocated(wave)) deallocate(wave, stat=nerr)
!  write(unit=*,fmt=*) " > ExportWave: EXIT."
  return

END SUBROUTINE ExportWave
!**********************************************************************!























!**********************************************************************!
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
!SUBROUTINE COMMENT TEMPLATE
!**********************************************************************!
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