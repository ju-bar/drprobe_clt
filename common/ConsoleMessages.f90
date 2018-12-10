! ---------------------------------------------------------------------
!
! File: 'ConsoleMessaged.f90'
! Author: J. Barthel, juribarthel@gmx.de
! Date: 05.10.2012
!
! Subroutines for display of a simple messages in console apps
!
! - use ConsoleMessages
!
! - call CSM_VLEVEL to initialize the verbosity level
! - call CSM_DLEVEL to initialize the current debug level
! - call CSM_SILENT to suppress all non-critical messaging
! - call CSM_TALKATIVE to activate the uppermost messaging levels
!
! - call CSM_post to post a normal message
! - call CSM_mustpost to post a message regardless of the VERBOSE mode
! - call CSM_debug to post a debug message
! - call CSM_error to post an error message
! - call CSM_criterr to post a critical error message + STOP !!!!
!
! - call CSM_GETTIMESTAMP to get a 17 character time stamp string
!
! - call CSM_SETLOG to activate message logging to a text file
! - call CSM_SETDBGLOG to activate debug logging to a text file
! - call CSM_STOPLOG to stop message logging to a text file
! - call CSM_STOPDBGLOG to stop debug logging to a text file
!
! ---------------------------------------------------------------------


!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
!
MODULE ConsoleMessages

  IMPLICIT NONE

  SAVE

!
!**********************************************************************!
!**********************************************************************!
!  
! ---------------------------------------------------------------------
! routine names declaration
!
  public :: CSM_VLEVEL
  public :: CSM_DLEVEL
  public :: CSM_SILENT
  public :: CSM_TALKATIVE
  public :: CSM_SETOUTUNIT
  public :: CSM_GETOUTUNIT
  
  public :: CSM_post
  public :: CSM_mustpost
  public :: CSM_debug
  public :: CSM_error
  public :: CSM_criterr
  
  public :: CSM_SETLOG
  public :: CSM_SETDBGLOG
  
  private :: CSM_GETLUN
  public :: CSM_GETTIMESTAMP
  
  private :: CSM_LOG
  private :: CSM_DBGLOG
  
  public :: CSM_STOPLOG
  public :: CSM_STOPDBGLOG
  
!
!
! ---------------------------------------------------------------------
! parameter declaration and definition
!
  integer*4, private, parameter :: CSM_OUTUNIT_DEFAULT = 6
  integer*4, private, parameter :: CSM_LEVEL_MAX = 5
  integer*4, private, parameter :: CSM_VERBOSE_LEVEL_DEFAULT = 1
  integer*4, private, parameter :: CSM_DEBUG_LEVEL_DEFAULT = 0
!
! variable declaration and default definition
  integer*4, private :: CSM_OUTUNIT
  DATA CSM_OUTUNIT /CSM_OUTUNIT_DEFAULT/
  integer*4, private :: CSM_DOLOG
  DATA CSM_DOLOG /0/
  integer*4, private :: CSM_DODBGLOG
  DATA CSM_DODBGLOG /0/
  integer*4, private :: CSM_VERBOSE_LEVEL
  DATA CSM_VERBOSE_LEVEL /CSM_VERBOSE_LEVEL_DEFAULT/
  integer*4, private :: CSM_DEBUG_LEVEL
  DATA CSM_DEBUG_LEVEL /CSM_DEBUG_LEVEL_DEFAULT/
  character(len=4096), private :: CSM_LOGFILE
  DATA CSM_LOGFILE /"message.log"/
  character(len=4096), private :: CSM_DBGLOGFILE
  DATA CSM_DBGLOGFILE /"debug.log"/
  integer*4, public :: CSM_NERR
  DATA CSM_NERR /0/
!
!
!
! ---------------------------------------------------------------------
!  
  CONTAINS
!
!**********************************************************************!
!**********************************************************************!
!


!**********************************************************************!
!
SUBROUTINE CSM_VLEVEL(verbose_level)
!
! purpose:
! - set a new verbose level
!   0 = off, no message posting
!
! interface:
! - integer*4, intent(in) :: verbose_level = the new verbose level
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 100
!
  integer*4, intent(in) :: verbose_level
!
! ---------------------------------------------------------------------
!
  CSM_VERBOSE_LEVEL = MIN(CSM_LEVEL_MAX,MAX(0,verbose_level))
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_VLEVEL
!
!**********************************************************************!




!**********************************************************************!
!
SUBROUTINE CSM_DLEVEL(debug_level)
!
! purpose:
! - set a new debug level
!   0 = off, no debug message posting
!
! interface:
! - integer*4, intent(in) :: verbose_level = the new verbose level
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 200
!
  integer*4, intent(in) :: debug_level
!
! ---------------------------------------------------------------------
!
  CSM_DEBUG_LEVEL = MIN(CSM_LEVEL_MAX,MAX(0,debug_level))
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_DLEVEL
!
!**********************************************************************!



!**********************************************************************!
!
SUBROUTINE CSM_SILENT()
!
! purpose:
! - sets the message posting and debug levels to zero
!
! interface:
! - no params
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 300
!
! ---------------------------------------------------------------------
!
  CSM_DEBUG_LEVEL = 0
  CSM_VERBOSE_LEVEL = 0
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_SILENT
!
!**********************************************************************!



!**********************************************************************!
!
SUBROUTINE CSM_TALKATIVE()
!
! purpose:
! - sets the message posting and debug levels to maximum
!
! interface:
! - no params
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 400
!
! ---------------------------------------------------------------------
!
  CSM_DEBUG_LEVEL = CSM_LEVEL_MAX
  CSM_VERBOSE_LEVEL = CSM_LEVEL_MAX
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_TALKATIVE
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_SETOUTUNIT(nou)
!
! purpose:
! - sets a new logical output unit (default: 6)
!
! interface:
! - integer*4, intent(in) :: nou = the new output unit
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 500
  integer*4, intent(in) :: nou
!
! ---------------------------------------------------------------------
!
  CSM_OUTUNIT = nou
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_SETOUTUNIT
!
!**********************************************************************!



!**********************************************************************!
!
SUBROUTINE CSM_GETOUTUNIT(nou)
!
! purpose:
! - retrieves the current output unit
!
! interface:
! - integer*4, intent(out) :: nou = the current output unit
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, intent(inout) :: nou
  integer*4, parameter :: subnum = 600
!
! ---------------------------------------------------------------------
!
  nou = CSM_OUTUNIT
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_GETOUTUNIT
!
!**********************************************************************!



!**********************************************************************!
!
SUBROUTINE CSM_post(smsg, nlev)
!
! purpose:
! - posts a message string to the console output
!
! interface:
! - character(len=*), intent(in) :: smsg
! - integer*4, intent(in), optional :: nlev
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 700
!
  character(len=*), intent(in) :: smsg
  integer*4, intent(in), optional :: nlev
!
! ---------------------------------------------------------------------
!
  if (PRESENT(nlev)) then
    if (nlev>CSM_VERBOSE_LEVEL) return
  end if
!
  if (0==CSM_VERBOSE_LEVEL) return
!
  if (CSM_OUTUNIT>0) write(unit=CSM_OUTUNIT,fmt='(A)') trim(smsg)
  if (CSM_DOLOG>0) call CSM_LOG(trim(smsg))
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_post
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_mustpost(smsg)
!
! purpose:
! - posts a message string to the console output
!   independent of the verbose level
!
! interface:
! - character(len=*), intent(in) :: smsg
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 800
!
  character(len=*), intent(in) :: smsg
!
! ---------------------------------------------------------------------
!
  if (CSM_OUTUNIT>0) write(unit=CSM_OUTUNIT,fmt='(A)') trim(smsg)
  if (CSM_DOLOG>0) call CSM_LOG(trim(smsg))
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_mustpost
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_debug(smsg,nlev)
!
! purpose:
! - posts a debug message string to the console output
!
! interface:
! - character(len=*), intent(in) :: smsg
! - integer*4, intent(in), optional :: nlev
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 900
!
  character(len=*), intent(in) :: smsg
  integer*4, intent(in), optional :: nlev
!
! ---------------------------------------------------------------------
!
  if (PRESENT(nlev)) then
    if (nlev>CSM_DEBUG_LEVEL) return
  end if
!
  if (0==CSM_DEBUG_LEVEL) return
!
  if (CSM_OUTUNIT>0) write(unit=CSM_OUTUNIT,fmt='(A)') " dbg  > "//trim(smsg)
  if (CSM_DODBGLOG>0) call CSM_DBGLOG(trim(smsg))
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_debug
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_error(smsg,ncode)
!
! purpose:
! - posts an error message string to the console output
!
! interface:
! - character(len=*), intent(in) :: smsg
! - integer*4, intent(in), optional :: ncode
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1000
!
  character(len=*), intent(in) :: smsg
  integer*4, intent(in), optional :: ncode
!  
  character(len=4096) :: sumsg
!
! ---------------------------------------------------------------------
!
  if (PRESENT(ncode)) then
    write(unit=sumsg, fmt='(A,I)') " err  > "//trim(smsg)//" - code: ", ncode
  else
    write(unit=sumsg, fmt='(A)') " err  > "//trim(smsg)
  end if
  CSM_NERR = CSM_NERR + 1
  if (CSM_OUTUNIT>0) write(unit=CSM_OUTUNIT,fmt='(A)') trim(sumsg)
  if (CSM_DOLOG/=0) call CSM_LOG(trim(sumsg))
  if (CSM_DODBGLOG/=0) call CSM_DBGLOG(trim(sumsg))
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_error
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_criterr(smsg, ncode)
!
! purpose:
! - posts a critical error message string to the console output
!   AND STOPS THE PROGRAM
!
! interface:
! - character(len=*), intent(in) :: smsg
! - integer*4, intent(in), optional :: ncode
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1100
!
  character(len=*), intent(in) :: smsg
  integer*4, intent(in), optional :: ncode
!  
  character(len=4096) :: sumsg
!
! ---------------------------------------------------------------------
!
  if (PRESENT(ncode)) then
    write(unit=sumsg, fmt='(A,I)') " crit > "//trim(smsg)//" - code: ", ncode
  else
    write(unit=sumsg, fmt='(A)') " crit > "//trim(smsg)
  end if
  CSM_NERR = CSM_NERR + 1
  
  if (CSM_OUTUNIT>0) write(unit=CSM_OUTUNIT,fmt='(A)') trim(sumsg)
  if (CSM_DOLOG/=0) call CSM_LOG(trim(sumsg))
  if (CSM_DODBGLOG/=0) call CSM_DBGLOG(trim(sumsg))
!
  sumsg = " crit > Halting program."
  if (CSM_OUTUNIT>0) write(unit=CSM_OUTUNIT,fmt='(A)') trim(sumsg)
  if (CSM_DOLOG/=0) call CSM_LOG(trim(sumsg))
  if (CSM_DODBGLOG/=0) call CSM_DBGLOG(trim(sumsg))
!
  STOP
!
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_criterr
!
!**********************************************************************!



!**********************************************************************!
!
SUBROUTINE CSM_SETLOG(sfname, nerase, nerr)
!
! purpose:
! - sets log file name and tests writing to the log file
!
! interface:
! - character(len=*) :: sfname :: log file name (input)
! - integer*4 :: nerase :: flag: 1=erase existing log (input)
! - integer*4 :: nerr :: error code
!                0      = success (input&output)
!                1401   = failed to get free logical file unit
!                1402   = failed to open / create / erase file for writing
!                1403   = failed to open / create file for writing
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1400
  
  character(len=*), intent(in) :: sfname
  integer*4, intent(in) :: nerase
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lun, ioerr
  character(len=30) :: st
!
! ---------------------------------------------------------------------
!
  CSM_DOLOG = 0
  call CSM_GETLUN(lun)
  if (lun<0) goto 13
  
  if (nerase/=0) then
    open (unit=lun, file=trim(sfname), STATUS='REPLACE', &
     &    ACTION='WRITE', FORM='FORMATTED', iostat=ioerr )
    if (ioerr/=0) goto 14
  else
    open (unit=lun, file=trim(sfname), STATUS='UNKNOWN', &
     &    ACTION='WRITE', FORM='FORMATTED', &
     &    POSITION='APPEND', iostat=ioerr )
    if (ioerr/=0) goto 15
  end if
  
  write (unit=lun, iostat=ioerr, fmt='(A)') ""
  call CSM_GETTIMESTAMP(st)
  write (unit=lun, iostat=ioerr, fmt='(A)') "Log started - "//trim(st)
  close (unit=lun, iostat=ioerr)
  
  CSM_LOGFILE = trim(sfname) ! store the log file name
  CSM_DOLOG = 1 ! activate logging
  
  return
!
! ---------------------------------------------------------------------
!
13 nerr = subnum + 1
  return
14 nerr = subnum + 2
  return
15 nerr = subnum + 3
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_SETLOG
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_SETDBGLOG(sfname, nerase, nerr)
!
! purpose:
! - sets debug log file name and tests writing to the file
!
! interface:
! - character(len=*) :: sfname :: file name (input)
! - integer*4 :: nerase :: flag: 1=erase existing log (input)
! - integer*4 :: nerr :: error code
!                0      = success (input&output)
!                1501   = failed to get free logical file unit
!                1502   = failed to open / create / erase file for writing
!                1503   = failed to open / create file for writing
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1500
  
  character(len=*), intent(in) :: sfname
  integer*4, intent(in) :: nerase
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lun, ioerr
  character(len=30) :: st
!
! ---------------------------------------------------------------------
!
  CSM_DODBGLOG = 0
  call CSM_GETLUN(lun)
  if (lun<0) goto 13
  
  if (nerase/=0) then
    open (unit=lun, file=trim(sfname), STATUS='REPLACE', &
     &    ACTION='WRITE', FORM='FORMATTED', iostat=ioerr )
    if (ioerr/=0) goto 14
  else
    open (unit=lun, file=trim(sfname), STATUS='UNKNOWN', &
     &    ACTION='WRITE', FORM='FORMATTED', &
     &    POSITION='APPEND', iostat=ioerr )
    if (ioerr/=0) goto 15
  end if
  
  write (unit=lun, iostat=ioerr, fmt='(A)') ""
  call CSM_GETTIMESTAMP(st)
  write (unit=lun, iostat=ioerr, fmt='(A)') "Debug log started - "//trim(st)
  close (unit=lun, iostat=ioerr)
  
  CSM_DBGLOGFILE = trim(sfname) ! store the log file name
  CSM_DODBGLOG = 1 ! activate logging
  
  return
!
! ---------------------------------------------------------------------
!
13 nerr = subnum + 1
  return
14 nerr = subnum + 2
  return
15 nerr = subnum + 3
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_SETDBGLOG
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_GETLUN(lun)
!
! purpose:
! - returns a free logical unit number for file input and output
!   between 20 and 99
!
! interface:
! - integer*4 :: lun (inout) ( < 0 in case of error )
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1600
  integer*4, parameter :: minlun = 20
  integer*4, parameter :: maxlun = 99
  
  integer*4, intent(inout) :: lun
  
  integer*4 :: l, i             ! iterators and lun temps
  integer*4 :: nerr             ! error code
  logical :: isopen             ! available lun flag
!
! ---------------------------------------------------------------------
!
  l = -1
  do i=minlun, maxlun
    inquire(unit=i,opened=isopen,iostat=nerr)
    if (.not.isopen) then
      l = i
      exit
    end if
  end do
  lun = l
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_GETLUN
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_GETTIMESTAMP(stime)
!
! purpose:
! - returns current time, format "YY/MM/DD hh:mm:ss"
!
! interface:
! - character(len=*) :: stime :: output (min length 17)
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1700
!
  character(len=*), intent(inout) :: stime
!
  character*8 :: sysdate
  character*10 :: systime
  character*5 :: syszone
  integer, dimension(8) :: values
!
! ---------------------------------------------------------------------
!
  call date_and_time(sysdate,systime,syszone,values)
!
  stime=sysdate(3:4)//"/"//sysdate(5:6)//"/"//sysdate(7:8)//" "&
     & //systime(1:2)//":"//systime(3:4)//":"//systime(5:6)
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_GETTIMESTAMP
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_LOG(sline)
!
! purpose:
! - appends a line to the current log file
!
! interface:
! - character(len=*) :: sline :: line to write (input)
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1800
!
  character(len=*), intent(in) :: sline
!
  integer*4 :: lun, ioerr
!
! ---------------------------------------------------------------------
!
!
  call CSM_GETLUN(lun)
  if (lun<0) return
  
  open  (unit=lun, file=trim(CSM_LOGFILE), STATUS='UNKNOWN', &
     &   ACTION='WRITE', FORM='FORMATTED', &
     &   POSITION='APPEND', iostat=ioerr )
  if (ioerr/=0) return
  
  write (unit=lun, iostat=ioerr, fmt='(A)') trim(sline)
  
  close (unit=lun, iostat=ioerr)

  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_LOG
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_DBGLOG(sline)
!
! purpose:
! - appends a line to the current debug log file
!
! interface:
! - character(len=*) :: sline :: line to write (input)
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 1900
!
  character(len=*), intent(in) :: sline
!
  integer*4 :: lun, ioerr
!
! ---------------------------------------------------------------------
!
!
  call CSM_GETLUN(lun)
  if (lun<0) return
  
  open  (unit=lun, file=trim(CSM_DBGLOGFILE), STATUS='UNKNOWN', &
     &   ACTION='WRITE', FORM='FORMATTED', &
     &   POSITION='APPEND', iostat=ioerr )
  if (ioerr/=0) return
  
  write (unit=lun, iostat=ioerr, fmt='(A)') trim(sline)
  
  close (unit=lun, iostat=ioerr)

  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_DBGLOG
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_STOPLOG()
!
! purpose:
! - stops message logging
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 2000
!
! ---------------------------------------------------------------------
!
  CSM_DOLOG = 0
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_STOPLOG
!
!**********************************************************************!


!**********************************************************************!
!
SUBROUTINE CSM_STOPDBGLOG()
!
! purpose:
! - stops debug logging
!
! ---------------------------------------------------------------------
!
  implicit none
!
  integer*4, parameter :: subnum = 2100
!
! ---------------------------------------------------------------------
!
  CSM_DODBGLOG = 0
!
  return
!
! ---------------------------------------------------------------------
!
END SUBROUTINE CSM_STOPDBGLOG
!
!**********************************************************************!


!  
!**********************************************************************!
!**********************************************************************!
!
END MODULE ConsoleMessages
!
!**********************************************************************!
!**********************************************************************!
!**********************************************************************!














!**********************************************************************!
!
! SUBROUTINE TEMPLATE
!
!**********************************************************************!




!**********************************************************************!
!
! SUBROUTINE TEMPLATE()
!
! purpose:
! - 
!
! interface:
! -
!
! ---------------------------------------------------------------------
!
!  implicit none
!
!  integer*4, parameter :: subnum = 2200
!
! ---------------------------------------------------------------------
!
!  return
!
! ---------------------------------------------------------------------
!
! END SUBROUTINE TEMPLATE
!
!**********************************************************************!