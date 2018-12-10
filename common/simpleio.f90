!**********************************************************************!
!**********************************************************************!
!
! F90 Source - simpleio
!
! Author: Dr. J. Barthel
!         Forschungszentrum Juelich GmbH
!         Jülich, Germany
!         ju.barthel@fz-juelich.de
!         first version: 12.01.2016
!         last version: 12.01.2016
!
! Purpose: Provides simple input and output routines used by many
!          programs.
!
! Linkage: none
!
!**********************************************************************!
!**********************************************************************!



!**********************************************************************!
!
! subroutine PostWarning
!
! posts a warning message to stdio.
!
! INPUT:
!   character(len=*) :: smessage = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostWarning(smessage)
  implicit none
  character*(*) :: smessage
  write (unit=*,fmt='(A)') " > Warning: "//trim(smessage)
  return
end subroutine PostWarning


!**********************************************************************!
!
! subroutine PostMessage
!
! posts a normal message to stdio.
!
! INPUT:
!   character(len=*) :: smessage = the meassage as string
!
! IN/OUTPUT: none
!
subroutine PostMessage(smessage)
  implicit none
  character*(*) :: smessage
  write (unit=*,fmt='(A)') " > "//trim(smessage)
  return
end subroutine PostMessage


!**********************************************************************!
!
! function getfreelun
!
! returns a free logical unit number for file input and output
! between 20 and 99
!
! INPUT: none
!
! OUTPUT: integer*4 >= 20: logical unit number 
!                   <  0 : error code
!
integer*4 function getfreelun()
  implicit none
  integer*4, parameter :: minlun = 20
  integer*4, parameter :: maxlun = 99
  integer*4 :: lun, i           ! iterators and lun temps
  integer*4 :: nerr             ! error code
  logical :: isopen             ! available lun flag
  lun = -1
  do i=minlun, maxlun
    inquire(unit=i,opened=isopen,iostat=nerr)
    if (.not.isopen) then
      lun = i
      exit
    end if
  end do
  getfreelun = lun
  return
end function getfreelun


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