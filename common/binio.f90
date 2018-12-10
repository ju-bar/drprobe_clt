!**********************************************************************!
!**********************************************************************!
!                                                                      !
!                    file   "binio".f95                                !
!                                                                      !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                  !
!    Version  :  1.0.0, May 09, 2005                                   !
!                                                                      !
!                                                                      !
!**********************************************************************!


!**********************************************************************!
!                                                                      !
! Purpose: implementing subroutines for buffered binary data in- and   !
!          output to fortran file units                                !
!                                                                      !
!**********************************************************************!
!                                                                      !
!   External  sources to link:                                         !
!     F90/95:                                                          !
!     F77:                                                             !
!     LIBS:                                                            !
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
MODULE binio

  implicit none

  SAVE

! internal long string length and buffer length in bytes
  integer*4, private, parameter :: BIO_ll = 1024


! declare internal data types
  type, public :: FileHandle ! Can be used outside this unit (is public)
    private                  ! but all fields are unknown outside (private)
    character(len=BIO_ll), pointer  :: sbuf
    character(len=BIO_ll)           :: sname
    integer*4                       :: nFilePos, nBufPos
    logical                         :: bForWriting
    logical                         :: bConnected
  end type FileHandle


! accessibility of subroutines or functions
!  private :: 
  private :: BIO_ERROR
  private :: BIO_createhandle, BIO_deletehandle
  private :: BIO_flush, BIO_endread
  private :: BIO_bufread, BIO_bufwrite
  
!  public :: 
  public :: BIO_INIT, BIO_UNINIT
  public :: BIO_READ_connect, BIO_WRITE_connect, BIO_disconnect
  public :: BIO_GetNextFreeLU
  public :: BIO_READ_skip, BIO_WRITE_skip
  public :: BIO_READ_char, BIO_WRITE_char
  public :: BIO_READ_i8, BIO_WRITE_i8
  public :: BIO_READ_i16, BIO_WRITE_i16
  public :: BIO_READ_i32, BIO_WRITE_i32
  public :: BIO_READ_r32, BIO_WRITE_r32
  public :: BIO_READ_r64, BIO_WRITE_r64
    
! declare module global (but private) variables, params and arrays

! max. number of simultaneously connectable files
  integer*4, private, parameter :: BIO_files_max = 20

! file unit limits
  integer*4, private, parameter :: BIO_minunit = 10
  integer*4, private, parameter :: BIO_maxunit = 1000
 
! empty line memory
  character(len=BIO_ll), private :: BIO_el

! file info handles
  type(FileHandle), private, dimension(BIO_files_max) :: BIO_fhInfo
! file unit to handle hasher
  integer*4, private, dimension(BIO_minunit:BIO_maxunit) :: BIO_nfh


! declare module global (but public) variables, params and arrays
  integer*4, public :: BIO_endian_swap
  

  CONTAINS





























!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!
!* >>
!* >>
!* >> SETUP & INIT ROUTINES

!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_INIT()
! function: init module
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
!  write(unit=*,fmt=*) " > BIO_INIT: INIT."
! ------------


! ------------
  BIO_el = REPEAT(" ",BIO_ll)
  BIO_nfh(:) = 0
  BIO_endian_swap = 0
! ------------


! ------------
!  write(unit=*,fmt=*) " > BIO_INIT: EXIT."
  return

END SUBROUTINE BIO_INIT
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_UNINIT()
! function: init module
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 200
  integer*4 :: i, nstat
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > BIO_UNINIT: INIT."
! ------------


! ------------
! disconnect all files
  do i=BIO_minunit, BIO_maxunit

    nstat = 0

    call BIO_disconnect(i,nstat)
    
    if (nstat/=0) then
      call BIO_ERROR("Logical unit disconnection problem.", subnum+1)
    end if
    
  end do
! ------------


! ------------
!  write(unit=*,fmt=*) " > BIO_UNINIT: EXIT."
  return

END SUBROUTINE BIO_UNINIT
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
SUBROUTINE BIO_createhandle(nunit, sfname, bopenforwrite, nstat)
! function: creates a file handle for a file to open
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            sfname : character(len=*) : file name
!            bopenforwrite : logical : write-flag
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 300
  
  integer*4, intent(in) :: nunit
  character(len=*), intent(in) :: sfname
  logical, intent(in) :: bopenforwrite
  integer*4, intent(out) :: nstat
  
  integer*4 :: i, nalloc, nfree
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_createhandle: INIT."
  nstat = 0
! ------------


! ------------
! find free handle
!   check unit connection
  if (nunit<BIO_minunit.or.nunit>BIO_maxunit) then
    nstat = -1
    call BIO_ERROR("Invalid unit number.",subnum+1)
    write(unit=*,fmt=*) nunit
    return
  end if
  if (BIO_nfh(nunit)/=0) then
    nstat = -2
    call BIO_ERROR("Unit already connected.",subnum+2)
    return
  end if
  
!   search all handles
  nfree = 0
  do i=1, BIO_files_max
    if (.not.BIO_fhInfo(i)%bConnected) then
      nfree = i
      exit ! finder loop
    end if
  end do
  
!   catch no handle
  if (nfree==0) then
    nstat = -3
    call BIO_ERROR("No free file handle.",subnum+3)
    return
  end if
! ------------


! ------------
! create unit and reserve file handle
  
  allocate( BIO_fhInfo(nfree)%sbuf, stat=nalloc)
  if (nalloc/=0) then
    nstat = -4
    call BIO_ERROR("Buffer allocation failed.", &
     & subnum+4)
    return
  end if
  BIO_nfh(nunit) = nfree
  BIO_fhInfo(nfree)%sname = BIO_el
  BIO_fhInfo(nfree)%sname = TRIM(sfname)
  BIO_fhInfo(nfree)%nFilePos = 0
  BIO_fhInfo(nfree)%nBufPos = 0
  BIO_fhInfo(nfree)%bConnected = .false.
  BIO_fhInfo(nfree)%bForWriting = bopenforwrite
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_createhandle: EXIT."
  return

END SUBROUTINE BIO_createhandle
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_deletehandle(nunit, nstat)
! function: creates a file handle for a file to open
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            sfname : character(len=*) : file name
!            bopenforwrite : logical : write-flag
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 700
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat
  
  integer*4 :: i, nalloc, nfree
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_deletehandle: INIT."
  nstat = 0
! ------------


! ------------
! find free handle
!   check unit connection
  if (nunit<BIO_minunit.or.nunit>BIO_maxunit) then
    nstat = -1
    call BIO_ERROR("Invalid unit number.",subnum+1)
    return
  end if
  nfree = BIO_nfh(nunit)
  if (nfree==0) then
!    nstat = -2
!    call BIO_ERROR("BIO_createhandle: handle already deleted.",subnum+2)
!   **** NO ERROR, JUST NOTHING TO DO, RETURN TO CALLING SUB
    return
  end if
! ------------


! ------------
! delete unit and un-reserve file handle
  
  deallocate( BIO_fhInfo(nfree)%sbuf, stat=nalloc)
  if (nalloc/=0) then
    nstat = -3
    call BIO_ERROR("Buffer deallocation failed.", &
     & subnum+3)
    return
  end if
  BIO_nfh(nunit) = 0
  BIO_fhInfo(nfree)%sname = BIO_el
  BIO_fhInfo(nfree)%nFilePos = 0
  BIO_fhInfo(nfree)%nBufPos = 0
  BIO_fhInfo(nfree)%bForWriting = .false.
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_deletehandle: EXIT."
  return

END SUBROUTINE BIO_deletehandle
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_flush(nunit, nstat)
! function: writes the last records bytewise to the file
! **** supposing that nunit is valid
!      DO ONLY CALL THIS UNIT DIRECTLY BEFORE CLOSING A WRITE-TO-FILE
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 800
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*4 :: i, nh, fp
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_flush: INIT."
  nstat = 0
! ------------


! ------------
! get handle
  nh = BIO_nfh(nunit)
! ------------


! ------------
! close and reopen
  close(unit=nunit, iostat=nstat)
  if (nstat/=0) then
    call BIO_ERROR("Closing file failed.",subnum+1)
    return
  end if
  open (unit=nunit, file=trim(BIO_fhInfo(nh)%sname),form="unformatted",&
     & iostat=nstat, access="direct", recl=1, status="old", &
     & action="write")
  if (nstat/=0) then
    call BIO_ERROR( &
     & "Connection to file unit failed.", subnum+2)
    return
  end if
! ------------


! ------------
! write last buffer to file
  fp = BIO_fhInfo(nh)%nFilepos*BIO_ll
  do i=1, BIO_fhInfo(nh)%nBufPos
    write(unit=nunit, rec=fp+i, iostat=nstat) BIO_fhInfo(nh)%sbuf(i:i)
    if (nstat/=0) then
      call BIO_ERROR( &
     & "Writing to file failed.", subnum+3)
      return
    end if
  end do
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_flush: EXIT."
  return

END SUBROUTINE BIO_flush
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_endread(nunit, newpos, nstat)
! function: reads from file to buffer unitl end of file or end of buffer
! **** supposing that nunit is valid
!      DO ONLY CALL THIS UNIT IF THe STD BLOCK READ FUNCTION FAILED DUE
!      TO EOF, CLOSE THE FILE AFTERWARDS
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 900

  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: newpos, nstat

  integer*4 :: i, nh, fp, nread
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_endread: INIT."
  nstat = 0
! ------------


! ------------
! get handle
  nh = BIO_nfh(nunit)
! ------------


! ------------
! close and reopen
  close(unit=nunit, iostat=nstat)
  if (nstat/=0) then
    call BIO_ERROR("Closing file failed.",subnum+1)
    return
  end if
  open (unit=nunit, file=trim(BIO_fhInfo(nh)%sname),form="unformatted",&
     & iostat=nstat, access="direct", recl=1, status="old", &
     & action="read", share='DENYNONE' )
  if (nstat/=0) then
    call BIO_disconnect(nunit, nstat)
    call BIO_ERROR( &
     & "Connection to file unit failed.", subnum+2)
    return
  end if
! ------------

! ------------
! read last bytes to buffer
  nread = 0
  fp = BIO_fhInfo(nh)%nFilepos*BIO_ll
  do i=1, BIO_ll
    read(unit=nunit, rec=fp+i, iostat=nstat) BIO_fhInfo(nh)%sbuf(i:i)
    if (nstat/=0) then
!     stop reading
      nread = i-1
      exit
    end if
  end do
! ------------


! ------------
! new position in buffer
  newpos = BIO_ll - nread + 1
! check number of bytes read
  if (nread>0) then
!   some bytes were read
    nstat = 0
!   copy read data to end of buffer, so that buffer reading ends
!   correctly then
    BIO_fhInfo(nh)%sbuf(newpos:BIO_ll) = BIO_fhInfo(nh)%sbuf(1:nread)
  end if
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_endread: EXIT."
  return

END SUBROUTINE BIO_endread
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_bufread(nunit, nstat)
! function: reads the next buffer and updates file handle infos
!  **** supposing that nunit is valid
!  **** WILL READ NEW DATA TO BUFFER OVERWRITING OLD BUFFER CONTENT
!       BUFFER POSITION WILL BE RESET
!  **** returns /=0 if EOF or other reading problem occured
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1000
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat
  
  integer*4 :: nh, fp, bp
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_bufread: INIT."
  nstat = 0
  nh = BIO_nfh(nunit)
  fp = BIO_fhInfo(nh)%nFilePos
! ------------


! ------------
! try read
  read(unit=nunit,rec=fp+1,iostat=nstat) BIO_fhInfo(nh)%sbuf
  if (nstat==0) then
    BIO_fhInfo(nh)%nBufPos = 1
  else
    call BIO_endread(nunit, bp, nstat)
    BIO_fhInfo(nh)%nBufPos = bp
  end if
  BIO_fhInfo(nh)%nFilePos = fp + 1
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_bufread: EXIT."
  return

END SUBROUTINE BIO_bufread
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_bufwrite(nunit, nstat)
! function: writes the current buffer to file
!  **** supposing that nunit is valid
!  **** BUFFER POSITION WILL BE RESETTED
!  **** returns /=0 if writing problem occured
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1100
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat
  
  integer*4 :: nh, fp, bp
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_bufwrite: INIT."
  nstat = 0
  nh = BIO_nfh(nunit)
  fp = BIO_fhInfo(nh)%nFilePos
! ------------


! ------------
! try write
  write(unit=nunit,rec=fp+1,iostat=nstat) BIO_fhInfo(nh)%sbuf
  if (nstat/=0) then
    call BIO_ERROR("Writing failed.",subnum+1)
    return
  end if
  BIO_fhInfo(nh)%nFilePos = fp + 1
  BIO_fhInfo(nh)%nBufPos = 0
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_bufwrite: EXIT."
  return

END SUBROUTINE BIO_bufwrite
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
!* >> CALCULATIONS


!* << CALCULATIONS
!* <<
!* <<
!* >> ------------------------------------------------------------ << *!
!* >> ------------------------------------------------------------ << *!



















!**********************************************************************!
!********************* INTERFACE ROUTINES *****************************!
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_connect(nunit, sfname, nstat)
! function: connects a file to a unit and creates file handle
!           returns error code in nstat == 0 if successful
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4, : file unit to connect to
!            sfname : character(len=*) : file name
!            nstat : integer*4 : error code reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 400
  
  integer*4, intent(in) :: nunit
  character(len=*), intent(in) :: sfname
  integer*4, intent(out) :: nstat
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_connect: INIT."
  nstat = 0
! ------------


! ------------
! create a file handle / also checks validity of unit and handles
  call BIO_createhandle(nunit, trim(sfname), .false., nstat)
  if (nstat/=0) then
    call BIO_ERROR("Creation of handle failed.", &
     & subnum+1)
    return
  end if
! ------------


! ------------
! open the file
  open (unit=nunit, file=trim(sfname), form="unformatted",iostat=nstat,&
     &  access="direct", recl=BIO_ll, status="old", action="read", share='DENYNONE' )
  if (nstat/=0) then
    call BIO_ERROR("Connection to file unit failed.",&
     & subnum+2)
    return
  end if
  BIO_fhInfo(BIO_nfh(nunit))%bConnected = .true.
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_connect: EXIT."
  return

END SUBROUTINE BIO_READ_connect
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_connect(nunit, sfname, nstat)
! function: connects a file to a unit and creates file handle
!           returns error code in nstat == 0 if successful
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4, : file unit to connect to
!            sfname : character(len=*) : file name
!            nstat : integer*4 : error code reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 500
  
  integer*4, intent(in) :: nunit
  character(len=*), intent(in) :: sfname
  integer*4, intent(out) :: nstat
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_WRITE_connect: INIT."
  nstat = 0
! ------------


! ------------
! create a file handle / also checks validity of unit and handles
  call BIO_createhandle(nunit, trim(sfname), .true., nstat)
  if (nstat/=0) then
    call BIO_ERROR("Creation of handle failed.", &
     & subnum+1)
    return
  end if
! ------------


! ------------
! open the file
  open (unit=nunit, file=trim(sfname), form="unformatted",iostat=nstat,&
     &  access="direct", recl=BIO_ll, status="replace", action="write")
  if (nstat/=0) then
    call BIO_ERROR( &
     & "Connection to file unit failed.", subnum+2)
    return
  end if
  BIO_fhInfo(BIO_nfh(nunit))%bConnected = .true.
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_connect: EXIT."
  return

END SUBROUTINE BIO_WRITE_connect
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_disconnect(nunit, nstat)
! function: disconnects a file from a unit and calls last writing ops.
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : unit number
!            nstat : integer*4 : error return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 600
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat
  
  integer*4 :: nh, nalloc
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_disconnect: INIT."
  nstat = 0
! ------------


! ------------
! check unit connection state
  if (nunit<BIO_minunit.or.nunit>BIO_maxunit) then
    nstat = -1
    call BIO_ERROR("Invalid unit number.",subnum+1)
    return
  end if

! closing the unit
  nh = BIO_nfh(nunit)
  if (nh/=0) then
  
    if (BIO_fhInfo(nh)%bConnected) then

!     flushing last record if the unit was connected for writing
      if (BIO_fhInfo(nh)%bForWriting .and. BIO_fhInfo(nh)%nBufPos>0) then
        call BIO_flush(nunit, nstat)
      end if

      close(unit=nunit,iostat=nstat)
      if (nstat/=0) then
        call BIO_ERROR("Closing file failed.",subnum+2)
!        return
      end if
      BIO_fhInfo(nh)%bConnected=.false.
    end if
    
!   **** OTHERWISE NOTHING TO DISCONNECT, HANDLE WILL BE DELETED BELOW

  end if  
  
! delete the handle
  call BIO_deletehandle(nunit, nstat)
  if (nstat/=0) then
    call BIO_ERROR("Deleting handle failed.",subnum+3)
    return
  end if
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_disconnect: EXIT."
  return

END SUBROUTINE BIO_disconnect
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_skip(nunit, nskip, nstat)
! function: skips by read (full length)
!  **** supposing that nunit is valid
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nskip : integer*4 : number of bytes to skip
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10000

  integer*4, intent(in) :: nunit, nskip

  integer*4, intent(out) :: nstat
  
  integer*4 :: schlen, nh, bp, sl1, sl2, sl, bp1, bp2
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_skip: INIT."
  nstat = 0
  schlen = nskip
  if (schlen<1) return
!  if (nunit<BIO_minunit .or. nunit>BIO_maxunit) then
!    nstat = - 1
!    call BIO_ERROR("BIO_READ_skip: invalid unit number.",subnum+1)
!    return
!  end if
  nh = BIO_nfh(nunit)
!  if (nh<1.or.nh>BIO_files_max) then
!    nstat = - 2
!    call BIO_ERROR("BIO_READ_skip: invalid file handle.",subnum+2)
!    return
!  end if
!  if (.not.BIO_fhInfo(nh)%bConnected.or.BIO_fhInfo(nh)%bForWriting) then
!    nstat = - 3
!    call BIO_ERROR("BIO_READ_skip: file not connected for reading.", &
!     & subnum+3)
!    return
!  end if
! ------------


! ------------
! try reading
  sl1 = 1
  sl2 = 1
  sl = 0
  do ! continuous reader loop

    bp = BIO_fhInfo(nh)%nBufPos ! get current read-buffer position
    if (bp<1.or.bp>BIO_ll) then ! update read buffer if invalid position
      call BIO_bufread(nunit, nstat)
      if (nstat/=0) then
        ! return, EOF
        return
      end if

      bp = BIO_fhInfo(nh)%nBufPos ! get new read-buffer position
    end if

    if (sl>=schlen) exit ! exit condition: number of bytes read equals string length
    
    bp1 = bp ! set buffer start read-pos
    bp2 = MIN(bp+(schlen-sl)-1,BIO_ll) ! set buffer end read pos
    sl2 = sl1 + (bp2-bp1) ! new string end length
!    write(unit=*,fmt=*) " > BIO_READ_skip:",sl1,sl2,bp1,bp2
!    sch(sl1:sl2) = BIO_fhInfo(nh)%sbuf(bp1:bp2) ! copy from buffer to transfer var.
    sl1 = sl2 + 1 ! update string locs
    sl = sl2
!!! update buffer pos
    BIO_fhInfo(nh)%nBufPos = bp2+1
  end do
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_skip: EXIT."
  return

END SUBROUTINE BIO_READ_skip
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_skip(nunit, nskip, nstat)
! function: skips by write (full length)
!  **** supposing that nunit is valid
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nskip : integer*4 : number of bbytes to skip
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10100

  integer*4, intent(in) :: nunit, nskip
  
  integer*4, intent(out) :: nstat
  
  integer*4 :: schlen, nh, bp, sl, sl1, sl2, bp1, bp2, i
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > BIO_WRITE_skip: INIT ... |",nskip,"|"
  nstat = 0
  schlen = nskip
  if (schlen<1) return
!  if (nunit<BIO_minunit .or. nunit>BIO_maxunit) then
!    nstat = - 1
!    call BIO_ERROR("BIO_WRITE_skip: invalid unit number.",subnum+1)
!    return
!  end if
  nh = BIO_nfh(nunit)
!  if (nh<1.or.nh>BIO_files_max) then
!    nstat = - 2
!    call BIO_ERROR("BIO_WRITE_skip: invalid file handle.",subnum+2)
!    return
!  end if
!  if (.not.(BIO_fhInfo(nh)%bConnected.and.BIO_fhInfo(nh)%bForWriting)) &
!     & then
!    nstat = - 3
!    call BIO_ERROR("BIO_WRITE_skip: file not connected for writing.", &
!     & subnum+3)
!    return
!  end if
! ------------


! ------------
! try writing
  bp = BIO_fhInfo(nh)%nBufPos
!  write(unit=*,fmt=*) bp, schlen
  sl = 0
  do

    if (bp>=BIO_ll) then
      call BIO_bufwrite(nunit, nstat)
      if (nstat/=0) then
        call BIO_ERROR("Writing failed.", subnum+4)
        return
      end if
      bp = BIO_fhInfo(nh)%nBufPos
    end if

    if (sl>=schlen) exit

    bp1 = bp + 1
    bp2 = MIN(bp + (schlen-sl), BIO_ll)
    sl1 = sl + 1
    sl2 = sl1 + (bp2-bp1)
!    write(unit=*,fmt=*) sl1,sl2,"->",bp1,bp2,(ichar(sch(i:i)),i=sl1,sl2),BIO_fhInfo(nh)%nBufPos,BIO_fhInfo(nh)%nFilePos
!    BIO_fhInfo(nh)%sbuf(bp1:bp2) = sch(sl1:sl2)
    sl = sl2
    bp = bp2
    BIO_fhInfo(nh)%nBufPos = bp
  end do
!  write(unit=*,fmt=*) bp
  
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_skip: EXIT."
  return

END SUBROUTINE BIO_WRITE_skip
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_char(nunit, sch, nstat)
! function: reads from file buffer to char variable (full length)
!  **** supposing that nunit is valid
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            sch : character(len=*) : char variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10000

  integer*4, intent(in) :: nunit
  character(len=*), intent(out) :: sch
  integer*4, intent(out) :: nstat
  
  integer*4 :: schlen, nh, bp, sl1, sl2, sl, bp1, bp2
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_char: INIT."
  nstat = 0
  schlen = LEN(sch)
  if (schlen<1) return
!  if (nunit<BIO_minunit .or. nunit>BIO_maxunit) then
!    nstat = - 1
!    call BIO_ERROR("BIO_READ_char: invalid unit number.",subnum+1)
!    return
!  end if
  nh = BIO_nfh(nunit)
!  if (nh<1.or.nh>BIO_files_max) then
!    nstat = - 2
!    call BIO_ERROR("BIO_READ_char: invalid file handle.",subnum+2)
!    return
!  end if
!  if (.not.BIO_fhInfo(nh)%bConnected.or.BIO_fhInfo(nh)%bForWriting) then
!    nstat = - 3
!    call BIO_ERROR("BIO_READ_char: file not connected for reading.", &
!     & subnum+3)
!    return
!  end if
! ------------


! ------------
! try reading
  sl1 = 1
  sl2 = 1
  sl = 0
  do ! continuous reader loop

    bp = BIO_fhInfo(nh)%nBufPos ! get current read-buffer position
    if (bp<1.or.bp>BIO_ll) then ! update read buffer if invalid position
      call BIO_bufread(nunit, nstat)
      if (nstat/=0) then
        ! return, EOF
        return
      end if

      bp = BIO_fhInfo(nh)%nBufPos ! get new read-buffer position
    end if

    if (sl>=schlen) exit ! exit condition: number of bytes read equals string length
    
    bp1 = bp ! set buffer start read-pos
    bp2 = MIN(bp+(schlen-sl)-1,BIO_ll) ! set buffer end read pos
    sl2 = sl1 + (bp2-bp1) ! new string end length
    sch(sl1:sl2) = BIO_fhInfo(nh)%sbuf(bp1:bp2) ! copy from buffer to transfer var.
    sl1 = sl2 + 1 ! update string locs
    sl = sl2
!!! update buffer pos
    BIO_fhInfo(nh)%nBufPos = bp2+1
  end do
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_char: EXIT."
  return

END SUBROUTINE BIO_READ_char
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_char(nunit, sch, nstat)
! function: writes from char vaiable to file buffer (full length)
!  **** supposing that nunit is valid
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            sch : character(len=*) : char variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10100

  integer*4, intent(in) :: nunit
  character(len=*), intent(in) :: sch
  integer*4, intent(out) :: nstat
  
  integer*4 :: schlen, nh, bp, sl, sl1, sl2, bp1, bp2, i
! ------------


! ------------
! INIT
!  write(unit=*,fmt=*) " > BIO_WRITE_char: INIT ... |",sch,"|"
  nstat = 0
  schlen = LEN(sch)
  if (schlen<1) return
!  if (nunit<BIO_minunit .or. nunit>BIO_maxunit) then
!    nstat = - 1
!    call BIO_ERROR("BIO_WRITE_char: invalid unit number.",subnum+1)
!    return
!  end if
  nh = BIO_nfh(nunit)
!  if (nh<1.or.nh>BIO_files_max) then
!    nstat = - 2
!    call BIO_ERROR("BIO_WRITE_char: invalid file handle.",subnum+2)
!    return
!  end if
!  if (.not.(BIO_fhInfo(nh)%bConnected.and.BIO_fhInfo(nh)%bForWriting)) &
!     & then
!    nstat = - 3
!    call BIO_ERROR("BIO_WRITE_char: file not connected for writing.", &
!     & subnum+3)
!    return
!  end if
! ------------


! ------------
! try writing
  bp = BIO_fhInfo(nh)%nBufPos
!  write(unit=*,fmt=*) bp, schlen
  sl = 0
  do

    if (bp>=BIO_ll) then
      call BIO_bufwrite(nunit, nstat)
      if (nstat/=0) then
        call BIO_ERROR("Writing failed.", subnum+4)
        return
      end if
      bp = BIO_fhInfo(nh)%nBufPos
    end if

    if (sl>=schlen) exit

    bp1 = bp + 1
    bp2 = MIN(bp + (schlen-sl), BIO_ll)
    sl1 = sl + 1
    sl2 = sl1 + (bp2-bp1)
!    write(unit=*,fmt=*) sl1,sl2,"->",bp1,bp2,(ichar(sch(i:i)),i=sl1,sl2),BIO_fhInfo(nh)%nBufPos,BIO_fhInfo(nh)%nFilePos
    BIO_fhInfo(nh)%sbuf(bp1:bp2) = sch(sl1:sl2)
    sl = sl2
    bp = bp2
    BIO_fhInfo(nh)%nBufPos = bp
  end do
!  write(unit=*,fmt=*) bp
  
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_char: EXIT."
  return

END SUBROUTINE BIO_WRITE_char
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_i8(nunit,ndest,nstat)
! function: reads 1 byte integer
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            ndest : integer*1 : destination variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10200
  integer*4, parameter :: tbyte = 1
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*1, intent(out) :: ndest

  integer*4 :: i
  character(len=tbyte) :: ch
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_i8: INIT."
  nstat = 0
  ndest = 0
! ------------


! ------------
  call BIO_READ_char(nunit, ch, nstat)
  ndest = TRANSFER(ch,ndest)
!  ndest = ichar(ch)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_i8: EXIT."
  return

END SUBROUTINE BIO_READ_i8
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_i8(nunit,nsrc,nstat)
! function: writes 1 byte integer
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nsrc : integer*1 : source variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10300
  integer*4, parameter :: tbyte = 1
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*1, intent(in) :: nsrc

  integer*4 :: i
  character(len=tbyte) :: ch
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_WRITE_i8: INIT."
  nstat = 0
! ------------


! ------------
!  ch = char( modulo(nsrc,256) )
  ch = TRANSFER(nsrc,ch)
  call BIO_WRITE_char(nunit, ch, nstat)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_i8: EXIT."
  return

END SUBROUTINE BIO_WRITE_i8
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_i16(nunit,ndest,nstat)
! function: reads 2 byte integer
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            ndest : integer*2 : destination variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10400
  integer*4, parameter :: tbyte = 2
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*2, intent(out) :: ndest

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_i16: INIT."
  nstat = 0
  ndest = 0
! ------------


! ------------
  call BIO_READ_char(nunit, ch, nstat)
!  write(unit=*,fmt=*) ichar(ch(1:1)),ichar(ch(2:2))

  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if

  ndest = TRANSFER(chs,ndest)
! ------------


! ------------
!  write(unit=*,fmt=*) " > BIO_READ_i16: EXIT.",ndest
  return

END SUBROUTINE BIO_READ_i16
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_i16(nunit,nsrc,nstat)
! function: writes 2 byte integer
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nsrc : integer*2 : source variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10500
  integer*4, parameter :: tbyte = 2
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*2, intent(in) :: nsrc

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_WRITE_i16: INIT.",nsrc
  nstat = 0
! ------------


! ------------
  ch = TRANSFER(nsrc,ch)

  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  
  call BIO_WRITE_char(nunit, chs, nstat)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_i16: EXIT. |",chs,"|",ichar(chs(1:1)),ichar(chs(2:2))
  return

END SUBROUTINE BIO_WRITE_i16
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_i32(nunit,ndest,nstat)
! function: reads 4 byte integer
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            ndest : integer*4 : destination variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10600
  integer*4, parameter :: tbyte = 4
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*4, intent(out) :: ndest

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_i32: INIT."
  nstat = 0
  ndest = 0
! ------------


! ------------
  call BIO_READ_char(nunit, ch, nstat)

  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  

  ndest = TRANSFER(chs,ndest)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_i32: EXIT."
  return

END SUBROUTINE BIO_READ_i32
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_i32(nunit,nsrc,nstat)
! function: writes 4 byte integer
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            nsrc : integer*4 : source variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10700
  integer*4, parameter :: tbyte = 4
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  integer*4, intent(in) :: nsrc

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_WRITE_i32: INIT."
  nstat = 0
! ------------


! ------------
  ch = TRANSFER(nsrc,ch)
  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  
  call BIO_WRITE_char(nunit, chs, nstat)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_i32: EXIT."
  return

END SUBROUTINE BIO_WRITE_i32
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_r32(nunit,rdest,nstat)
! function: reads 4 byte real
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            rdest : real*4 : destination variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10800
  integer*4, parameter :: tbyte = 4
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  real*4, intent(out) :: rdest

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_r32: INIT."
  nstat = 0
  rdest = 0.0
! ------------


! ------------
  call BIO_READ_char(nunit, ch, nstat)

  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  

  rdest = TRANSFER(chs,rdest)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_r32: EXIT."
  return

END SUBROUTINE BIO_READ_r32
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_r32(nunit,rsrc,nstat)
! function: writes 4 byte real
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            rsrc : real*4 : source variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 10900
  integer*4, parameter :: tbyte = 4
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  real*4, intent(in) :: rsrc

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_WRITE_r32: INIT."
  nstat = 0
! ------------


! ------------
  ch = TRANSFER(rsrc,ch)
  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  
  call BIO_WRITE_char(nunit, chs, nstat)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_r32: EXIT."
  return

END SUBROUTINE BIO_WRITE_r32
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_READ_r64(nunit,rdest,nstat)
! function: reads 8 byte real
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            rdest : real*8 : destination variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 11000
  integer*4, parameter :: tbyte = 8
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  real*8, intent(out) :: rdest

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_READ_r64: INIT."
  nstat = 0
  rdest = 0.0
! ------------


! ------------
  call BIO_READ_char(nunit, ch, nstat)

  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  

  rdest = TRANSFER(chs,rdest)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_READ_r64: EXIT."
  return

END SUBROUTINE BIO_READ_r64
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_WRITE_r64(nunit,rsrc,nstat)
! function: writes 8 byte real
! -------------------------------------------------------------------- !
! parameter: nunit : integer*4 : file unit number
!            rsrc : real*8 : source variable
!            nstat : integer*4 : error status return reference
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 11100
  integer*4, parameter :: tbyte = 8
  
  integer*4, intent(in) :: nunit
  integer*4, intent(out) :: nstat

  real*8, intent(in) :: rsrc

  integer*4 :: i, j
  character(len=tbyte) :: ch, chs
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_WRITE_r64: INIT."
  nstat = 0
! ------------


! ------------
  ch = TRANSFER(rsrc,ch)
  chs = ch
  if (BIO_endian_swap/=0) then
    do i=tbyte,1,-1
      j = 1+tbyte-i
      chs(j:j) = ch(i:i)
    end do
  end if
  
  call BIO_WRITE_char(nunit, chs, nstat)
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_WRITE_r64: EXIT."
  return

END SUBROUTINE BIO_WRITE_r64
!**********************************************************************!
















!**********************************************************************!
!********************* SPECIAL STD ROUTINES ***************************!
!**********************************************************************!

!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_GetNextFreeLU(nunit)
! function: finds next free logical external file unit and returns this
!           unit in nunit, returns -1 if no free unit was found
! -------------------------------------------------------------------- !
! parameter: 
! -------------------------------------------------------------------- !

  implicit none

! ------------
! DECLARATION
  integer*4, parameter :: subnum = 1200
  
  integer*4, intent(inout) :: nunit
  integer*4 :: nu
  logical :: isopen
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > BIO_GetNextFreeLU: INIT."
  nu = MAX(nunit,BIO_minunit)
! ------------


! ------------
  do
  
    inquire(unit=nu,OPENED=isopen)
  
    if (.not.isopen) exit
        
    nu = nu + 1
    
    if (nu>BIO_maxunit) then
      nunit = -1
      return
    end if
    
  end do
! ------------


! ------------
! write(unit=*,fmt=*) " > BIO_GetNextFreeLU: EXIT."
  return

END SUBROUTINE BIO_GetNextFreeLU
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
SUBROUTINE BIO_ERROR(sTxt,nErr)
! function: print error message and stop program
! -------------------------------------------------------------------- !
! parameter: CHARACTER(LEN=*) :: sTxt : error message text
!            INTEGER*4 :: nErr : error code, displayed when /= 0
! -------------------------------------------------------------------- !

  implicit none

! ----------  
  character(len=*), intent(in) :: sTxt
  integer*4, intent(in) :: nErr
! ----------

! ----------
! init
! ----------

! ----------
! messaging (to screen)
  write(unit=*,fmt=*) " > ERROR: ",trim(sTxt)," Error code:",nErr
! ----------

! ----------
  return

END SUBROUTINE BIO_ERROR
!**********************************************************************!











END MODULE binio
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
!  integer*4, parameter :: subnum = 1400
! ------------


! ------------
! INIT
! write(unit=*,fmt=*) " > <NAME>: INIT."
! ------------


! ------------
! 
! ------------


! ------------
! write(unit=*,fmt=*) " > <NAME>: EXIT."
!  return

!END SUBROUTINE <NAME>
!**********************************************************************!

