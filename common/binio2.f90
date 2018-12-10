!**********************************************************************!
!**********************************************************************!
!
! FILE: "binio2.f90"
!
! AUTHOR: Dr. J. Barthel
!         Forschungszentrum Jülich
!         Jülich, Germany
!
! PURPOSE: Implementation for simple binary data input and output
!
! VERSION: 1.0, J.B., 18.12.2007
!          1.1, J.B., 27.08.2010
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


!**********************************************************************!
!
! subroutine saver4datatable
!   saves real*4 data as text table to file
!
! INPUT:
!   integer*4 :: ncol               ! number of data columns
!   integer*4 :: nrow               ! number of data rows
!   real*4 :: a(ncol,nrow)          ! data array
!   character(len=*) :: sf          ! file name
!   character(len=*) :: sh          ! header string
!   character(len=*) :: sep         ! separator string
!
! IN/OUTPUT:
!   integer*4 :: nerr               ! error code
!                                      0 = success
!                                      1 = fail, invalid size parameter
!                                      2 = fail, invalid separator
!                                      3 = fail, no free logical file unit
!                                      4 = fail, opening file
!                                      5 = fail, writing data
!
subroutine saver4datatable(a,ncol,nrow,sf,sh,sep,nerr)

  implicit none
  
  integer*4, intent(in) :: ncol, nrow
  real*4, intent(in) :: a(ncol,nrow)
  character(len=*), intent(in) :: sf, sh, sep
  integer*4, intent(inout) :: nerr
  
  integer*4 :: i, j                 ! iterators
  integer*4 :: lun                  ! logical file unit used
  integer*4 :: ierr                 ! internal error code
  integer*4 :: lensep               ! length of the separator string
  
  integer*4, external :: getfreelun
  external :: createfilefolder
  
  ierr = 0
  nerr = 0
  
  if (nrow<=0 .or. ncol<=0) then
    nerr = 1
    return
  end if
  
  lensep = len(sep)
  if (lensep<=0) then
    nerr = 2
    return
  end if
  
  lun = getfreelun()
  if (lun<=0) then
    nerr = 3
    return
  end if
  
  call createfilefolder(trim(sf),ierr)
  if (ierr/=0) then
    nerr = 40
    return
  end if
  open(unit=lun, file=trim(sf), action="write", status="replace", iostat=ierr)
  if (ierr/=0) then
    nerr = 4
    return
  end if
  
  write(unit=lun, fmt='(A)', iostat=ierr) sh
  if (ierr/=0) then
    nerr = 5
    goto 99
  end if
  
  do j=1, nrow
    do i=1, ncol
      write(unit=lun,fmt='(G13.5,$)') a(i,j)
      if (i<ncol) write(unit=lun,fmt='(A,$)') sep
    end do
    write(unit=lun, fmt='(A)', iostat=ierr) ""
    if (ierr/=0) then
      nerr = 5
      goto 99
    end if
  end do
  
99 close(unit=lun,iostat=ierr)
  
  return
  
end subroutine saver4datatable


!**********************************************************************!
!
! subroutine loaddata
!
! Loads data from file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   integer*4 :: noff           = data offset in bytes
!   integer*4 :: ntype          = data type
!                                   0: real*4
!                                   1: real*8
!                                   2: integer*1
!                                   3: integer*2
!                                   4: integer*4
!
! IN/OUTPUT
!   real*4 :: a(n)              = data array
!   integer*4 :: nerr           = error code
!
! WORKINGMEM:
!   templates for different data formats
!   integer*1, allocatable, dimension(:) :: odi1
!   integer*2, allocatable, dimension(:) :: odi2
!   integer*4, allocatable, dimension(:) :: odi4
!   real*4, allocatable, dimension(:) :: odr4
!   real*8, allocatable, dimension(:) :: odr8
!
subroutine loaddata(sfile, ntype, n, noff, a, nerr)
  USE IFPORT
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: ntype, n, noff
  real*4, intent(inout) :: a(1:n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4 :: ipos
  integer*4 :: i
  integer*4 :: alloc, nstat
  logical :: fexists
  integer*4, external :: getfreelun
  integer*1, allocatable, dimension(:) :: odi1
  integer*2, allocatable, dimension(:) :: odi2
  integer*4, allocatable, dimension(:) :: odi4
  real*4, allocatable, dimension(:) :: odr4
  real*8, allocatable, dimension(:) :: odr8
  !
  ! check data type and allocate working data
  !
  if (ntype>=0.and.ntype<=4) then
    alloc = 0
    select case (ntype)
      case (0)
        allocate(odr4(n),stat=alloc)
      case (1)
        allocate(odr8(n),stat=alloc)
      case (2)
        allocate(odi1(n),stat=alloc)
      case (3)
        allocate(odi2(n),stat=alloc)
      case (4)
        allocate(odi4(n),stat=alloc)
    end select
    if (alloc/=0) then
      nerr = 1
      return
    end if
  else
    nerr = -6
    return
  end if
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! check if file exists
  !
  inquire(file=trim(sfile),exist=fexists,iostat=nstat)
  if (.not.fexists) then
    nerr = -2
    return
  end if
  !
  ! open file, connect to lun
  !
  open(unit=lun,file=trim(sfile),iostat=nstat,&
     & form='binary',action='read',status='old')
  if (nstat/=0) then
    nerr = -3
    return
  end if
  !
  ! move to offset position
  !   be shure that RECL is 1 byte, check compiler options
  !   or that noff is set appropriately
  !
  if (noff>0) then
    ipos = 0
    nerr = fseek(lun,noff,ipos)
    if (nerr/=0) then
      nerr = -4
      return
    end if
  end if
  !
  ! read data from file, sequential binary
  !
  select case (ntype)
    case (0)
      read(unit=lun,iostat=nstat) odr4
    case (1)
      read(unit=lun,iostat=nstat) odr8
    case (2)
      read(unit=lun,iostat=nstat) odi1
    case (3)
      read(unit=lun,iostat=nstat) odi2
    case (4)
      read(unit=lun,iostat=nstat) odi4
  end select
  if (nstat/=0) then
    nerr = -5
    return
  end if
  !
  ! transfer data to real*4
  !
  select case (ntype)
    case (0)
      do i=1, n
        a(i) = odr4(i)
      end do
    case (1)
      do i=1, n
        a(i) = real(odr8(i),kind(4))
      end do
    case (2)
      do i=1, n
        a(i) = real(odi1(i),kind(4))
      end do
    case (3)
      do i=1, n
        a(i) = real(odi2(i),kind(4))
      end do
    case (4)
      do i=1, n
        a(i) = real(odi4(i),kind(4))
      end do
  end select
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! deallocate arrays
  !
  select case (ntype)
    case (0)
      deallocate(odr4,stat=alloc)
    case (1)
      deallocate(odr8,stat=alloc)
    case (2)
      deallocate(odi1,stat=alloc)
    case (3)
      deallocate(odi2,stat=alloc)
    case (4)
      deallocate(odi4,stat=alloc)
  end select
  if (alloc/=0) then
    nerr = 2
    return
  end if
  !
  ! done
  !
  nerr = 0
  return
end subroutine loaddata

!**********************************************************************!
!
! subroutine savedata
!
! Saves data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   real*4 :: a(n)              = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedata(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real*4, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedata


!**********************************************************************!
!
! subroutine savedata1
!
! Saves data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   real*4 :: a(n)              = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedata1(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real*4, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedata1


!**********************************************************************!
!
! subroutine savedatar8
!
! Saves real*8 data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   real*8 :: a(n)              = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedatar8(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  real*8, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedatar8

!**********************************************************************!
!
! subroutine savedatac8
!
! Saves complex*8 data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   complex*8 :: a(n)           = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedatac8(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  complex*8, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedatac8


!**********************************************************************!
!
! subroutine savedatai4
!
! Saves integer*4 data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   integer*4 :: a(n)           = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedatai4(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  integer*4, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedatai4


!**********************************************************************!
!
! subroutine savedatai2
!
! Saves integer*2 data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   integer*2 :: a(n)           = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedatai2(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  integer*2, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedatai2


!**********************************************************************!
!
! subroutine savedatai1
!
! Saves integer*1 data to file
!
! INPUT:
!   character(len=*) :: sfile   = file name string
!   integer*4 :: n              = number of samples
!   integer*1 :: a(n)           = data array
!
! IN/OUTPUT
!   integer*4 :: nerr           = error code
!
subroutine savedatai1(sfile, n, a, nerr)
  implicit none
  character(len=*), intent(in) :: sfile
  integer*4, intent(in) :: n
  integer*1, intent(inout) :: a(n)
  integer*4, intent(inout) :: nerr
  integer*4 :: lun
  integer*4, external :: getfreelun
  external :: createfilefolder
  !
  ! get logical unit number
  !
  lun = getfreelun()
  if (lun<0) then
    nerr = -1
    return
  end if
  !
  ! open file, connect to lun
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    nerr = -20
    return
  end if
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='write',status='replace')
  if (nerr/=0) then
    nerr = -2
    return
  end if
  !
  ! write data to file, sequential binary
  !
  write(unit=lun,iostat=nerr) a
  if (nerr/=0) then
    nerr = -3
    return
  end if
  !
  ! close and disconnect
  !
  close(unit=lun)
  !
  ! done
  !
  nerr = 0
  return
end subroutine savedatai1