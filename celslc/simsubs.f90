!**********************************************************************!
!**********************************************************************!
!
! FILE: "simsubs.f90"
!
! AUTHOR: Dr. J. Barthel
!         Ernst Ruska-Centre
!         Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
!           and
!         RWTH Aachen University, 52074 Aachen, Germany
!
! PURPOSE: Implementation of subroutines for CELSLC
!
! VERSION: 0.41b, J.B., 21.09.2015
!
!**********************************************************************!
!**********************************************************************!

!**********************************************************************!
!
! subroutine Introduce
!
! writes introduction info to console
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine Introduce
  
  implicit none
  
  include 'global.fi'
  
  character(len=5) :: snft
  
  write(unit=snft,fmt='(I5)') fft_dmax
  
  write(unit=stdout,fmt='(A)') ""
  write(unit=stdout,fmt='(A)') "  +---------------------------------------------------+"
  write(unit=stdout,fmt='(A)') "  | Program [celslc]                                  |" 
  write(unit=stdout,fmt='(A)') "  | Version: 0.41b (20150921) 64-bit                  |"
  write(unit=stdout,fmt='(A)') "  | Author : Dr. J. Barthel, FZ-Juelich, Germany      |"
  write(unit=stdout,fmt='(A)') "  +---------------------------------------------------+"
  write(unit=stdout,fmt='(A)') "  |          max. slice size: "//snft//"                   |"
  write(unit=stdout,fmt='(A)') "  +---------------------------------------------------+"
  write(unit=stdout,fmt='(A)') ""
  write(unit=stdout,fmt='(A)') ""
  
  return

end subroutine Introduce

!**********************************************************************!
!
! subroutine Outroduce
!
! writes introduction info to console
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine Outroduce
  
  implicit none
  
  include 'global.fi'
  
  write(unit=stdout,fmt='(A)') ""
  write(unit=stdout,fmt='(A)') ""
  
  return

end subroutine Outroduce


!**********************************************************************!
!
! subroutine CriticalError
!
! posts an error message to console and halts the program.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine CriticalError(smessage)

  implicit none
  
  include 'global.fi'

  character*(*) :: smessage

  write (unit=stdout,fmt='(A)') " > "
  write (unit=stdout,fmt='(A)') " > "//trim(smessage)
  write (unit=stdout,fmt='(A)') " > Critical error. Halting program."
  write (unit=stdout,fmt='(A)') " > "
  stop

  return

end subroutine CriticalError
!**********************************************************************!



!**********************************************************************!
!
! subroutine PostWarning
!
! posts a warning message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostWarning(smessage)

  implicit none
  
  include 'global.fi'

  character*(*) :: smessage

  write (unit=stdout,fmt='(A)') " > Warning: "//trim(smessage)

  return

end subroutine PostWarning


!**********************************************************************!
!
! subroutine PostMessage
!
! posts a normal message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostMessage(smessage)

  implicit none
  
  include 'global.fi'

  character*(*) :: smessage

  write (unit=stdout,fmt='(A)') " > "//trim(smessage)

  return

end subroutine PostMessage




!**********************************************************************!
!
! subroutine PostCellInfo
!
! posts a normal message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostCellInfo()

  use CellSlicer

  implicit none
  
  include "global.fi"

  character(len=400) :: smsg
  
  ! start
  call PostMessage("Info on loaded super-cell data:")
  
  sdx = CS_scsx / real(nx)
  sdy = CS_scsy / real(ny)
  sdz = CS_scsz / real(nz)
  
  write(unit=smsg,fmt='("Super-cell size (x,y,z) in nm: ",3F12.8)') CS_scsx, CS_scsy, CS_scsz 
  call PostMessage(trim(smsg))
  
  write(unit=smsg,fmt='("Super-cell sampling (dx,dy,dz) in nm: ",3F12.8)') sdx, sdy, sdz 
  call PostMessage(trim(smsg))
  
  write(unit=smsg,fmt='("Total number of atoms in super-cell: ",I5)') CS_numat
  call PostMessage(trim(smsg))
    
  return

end subroutine PostCellInfo





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
!  write(unit=stdout,fmt=*) " > factorial: INIT."
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
!  write(unit=stdout,fmt=*) " > factorial: EXIT."
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
!  write(unit=stdout,fmt=*) " > binomial: INIT."
  binomial = 0 ! precheck default -> this means ERROR!
  if (n<0.or.k<0.or.n<k) return
  binomial = factorial(n)/( factorial(n-k) * factorial(k) )
! ------------

! ------------
!  write(unit=stdout,fmt=*) " > binomial: EXIT."
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

  implicit none

! ------------
! declaration
  real*4 :: sigmoid
  real*4, intent(in) :: x, x0, dx
! ------------

! ------------
! init
!  write(unit=stdout,fmt=*) " > sigmoid: INIT."
  sigmoid = 0.5*(tanh((x-x0)/dx)+1.0)
! ------------



! ------------
!  write(unit=stdout,fmt=*) " > sigmoid: EXIT."
  return

END FUNCTION sigmoid
!**********************************************************************!




!**********************************************************************!
SUBROUTINE UPPERCASE(strin, strout)
  
  implicit none
  
  character(len=*), intent(in) :: strin
  character(len=len(strin)), intent(out) :: strout
  integer :: i,j
  
  ! define captial and lower-case characters to recognize for transform
  character(29), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZÄÖÜ'
  character(29), parameter :: low = 'abcdefghijklmnopqrstuvwxyzäöü'
  
  strout = strin
  
  do i=1, len_trim(strin)
    j = index(low,strin(i:i))
    if (j>0) strout(i:i) = cap(j:j)
  end do
  
  return
  
END SUBROUTINE UPPERCASE
!**********************************************************************!



!**********************************************************************!
SUBROUTINE LOWERCASE(strin, strout)
  
  implicit none
  
  character(len=*), intent(in) :: strin
  character(len=len(strin)), intent(out) :: strout
  integer :: i,j
  
  ! define captial and lower-case characters to recognize for transform
  character(29), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZÄÖÜ'
  character(29), parameter :: low = 'abcdefghijklmnopqrstuvwxyzäöü'
  
  strout = strin
  
  do i=1, len_trim(strin)
    j = index(cap,strin(i:i))
    if (j>0) strout(i:i) = low(j:j)
  end do
  
  return
  
END SUBROUTINE LOWERCASE
!**********************************************************************!





!**********************************************************************!
!
! subroutine ExplainUsage
!
! posts usage info to console.
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine ExplainUsage()

  implicit none
  
  include 'global.fi'
  character(len=600) :: stmp

  write(unit=stdout,fmt='(A)') ' > '
  write(unit=stdout,fmt='(A)') ' > CELSLC parameters'
  write(unit=stdout,fmt='(A)') ' >  -cel <string> = super-cell file (e.g. "atoms.cel")'
  write(unit=stdout,fmt='(A)') ' >  -slc <string> = slice file name (e.g. "scl")'
  write(unit=stdout,fmt='(A)') ' >  -nx <number>  = horizontal cell sampling (32 ... 2048)'
  write(unit=stdout,fmt='(A)') ' >  -ny <number>  = vertical cell sampling (32 ... 2048)'
  write(unit=stdout,fmt='(A)') ' >  -nz <number>  = number of slices (1 ... 2048)'
  write(unit=stdout,fmt='(A)') ' >  -ht <number>  = electron energy in keV (10 ... 1300)'
  write(unit=stdout,fmt='(A)') ' > [-rev         -> slicing in reversed sequence]'
  write(unit=stdout,fmt='(A)') ' > [-fl          -> frozen-lattice simulation]'
  write(unit=stdout,fmt='(A)') ' > [-nv <number>  = number of FL variants per slice (1 ... 2048)]'
  write(unit=stdout,fmt='(A)') ' > [-dwf         -> apply Debye-Waller factors]'
  write(unit=stdout,fmt='(A)') ' > [-abs         -> apply built-in absorption factors]'
  write(unit=stdout,fmt='(A)') ' > [-abf <number> = apply fix rel. absorption factor]'
  write(unit=stdout,fmt='(A)') ' > [-pot         -> export potentials to files *.pot]'
  write(unit=stdout,fmt='(A)') ' > [-ssc <number> = single slice calculation]'
  write(unit=stdout,fmt='(A)') ' > '
  write(unit=stdout,fmt='(A)') ' > [] = optional parameters'
  write(unit=stdout,fmt='(A)') ' > '
  
  return

end subroutine ExplainUsage


!**********************************************************************!
!
! subroutine ParseCommandLine
!
! parses the command line options and sets global variables.
! also performs parameter checks
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine ParseCommandLine()

  implicit none
  
  include "global.fi"

  character*512 :: buffer, cmd
  logical :: fex
  integer*4 :: lfu
  integer*4 :: i, cnt, status, len, nfound
  integer*4 :: nprm, nout, nht, nnx, nny, nnz
  integer*4, external :: getfreelun
  real*4, external :: HT2WL
  character(len=512) :: smsg

! ------------
! initialize
!  write(unit=stdout,fmt=*) " > ParseCommandLine: INIT."
  i = 0
  cnt = command_argument_count()
  if (cnt==0) then
    call ExplainUsage()
    call CriticalError("No arguments found.")
  end if
  nprm = 0
  nout = 0
  nnx = 0
  nny = 0
  nnz = 0
  nht = 0
  nfl = 0
  ndwf = 0
  nabs = 0
  nabf = 0
  npot = 0
  npot = 0
  nzd = 0
  nfx = 0
  nfe = 0
  nrev = 0
  nx = 0
  ny = 0
  nz = 0
  nv = 1
  nfin = 0 ! default input format = 0 = cel structure file
  abf = 0.0
  ssc = 0
  buni = 0
  buniv = 0.0

! ------------
! LOOP OVER ALL GIVEN ARGUMENTS
  do
    i = i + 1
    if (i>cnt) exit
    
    call get_command_argument (i, buffer, len, status)
    if (status/=0) then
      call ExplainUsage()
      call CriticalError("Command line parsing error.")
    end if
    
    ! CHECK COMMAND
    nfound = 0
    cmd = buffer(1:len)
    CHECK_COMMAND: select case (cmd(1:len))
    
    ! THE STRUCTURE PARAMETER FILE
    case ("-cel")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status)
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = scellfile, fmt='(A)') buffer(1:len)
      inquire(file=trim(scellfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid argument: Specified super-cell file ["//trim(scellfile)//"] does not exist.")
      end if
      nprm = 1
    
    ! THE OUTPUT FILE
    case ("-slc")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = sslcfile, fmt='(A)') buffer(1:len)
      nout = 1
    
    ! THE TEM HIGH-TENSION VALUE IN KILOVOLTS  
    case ("-ht")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) ht
      if (status/=0 .or. ht<10.0 .or. ht>1300.0) then
        call ExplainUsage()
        call CriticalError("Failed to read HT parameter.")
      end if
      wl = HT2WL(ht)
      nht = 1
      
    
      
    ! THE HOR DISRETIZATION
    case ("-nx")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) nx
      if (status/=0 .or. nx<fft_dmin .or. nx>fft_dmax) then
        call ExplainUsage()
        call CriticalError("Failed to read x-discretization.")
      end if
      nnx = 1
      
    ! THE VER DISRETIZATION
    case ("-ny")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) ny
      if (status/=0 .or. nx<fft_dmin .or. nx>fft_dmax) then
        call ExplainUsage()
        call CriticalError("Failed to read y-discretization.")
      end if
      nny = 1
      
    ! THE SLICE DISRETIZATION
    case ("-nz")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) nz
      if (status/=0 .or. nz<1 .or. nz>2048) then
        call ExplainUsage()
        call CriticalError("Failed to read number of slices.")
      end if
      nnz = 1
    
    ! THE NUMBER OF VARIANTS
    case ("-nv")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) nv
      if (status/=0 .or. nv<1 .or. nv>2048) then
        call ExplainUsage()
        call CriticalError("Failed to read number of variants per slice.")
      end if
      
    ! FROZEN LATTICE USAGE
    case ("-fl")
      nfound = 1
      nfl = 1 ! switch on frozen lattice calculation
      
    ! DEBYE-WALLER FACTOR USAGE
    case ("-dwf")
      nfound = 1
      ndwf = 1 ! switch on Debye-Waller factors
      
    ! ABSORPTION USAGE FROM WEICKENMEYER&KOHL, Acta Cryst. A 47 (1991) p. 597
    case ("-abs")
      nfound = 1
      nabs = 1 ! switch on absorption
      
    ! ABSORPTION USAGE FROM HASHIMOTO, HOWIE, WHELAN, Proc. R. Soc. London Ser. A, 269, 80-103 (1962)   
    case ("-abf")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! outputfile
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) abf
      if (status/=0 .or. abf<0.0 .or. abf>1.0) then
        call ExplainUsage()
        call CriticalError("Unreasonable absorption factor.")
      end if
      nabf = 1
      
    ! EXPORT POTENTIAL FILES
    case ("-pot")
      nfound = 1
      npot = 1 ! switch on potential export
      
    ! CREATE 3D POTENTIAL
    case ("-3dp")
      nfound = 1
      n3dp = 1 ! switch on 3d potential creation
      
    ! USE EXTRA SCATTERING X-RAY FACTORS
    case ("-fx")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! input file
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = sfxfile, fmt='(A)') buffer(1:len)
      nfx = 1
      
    ! USE EXTRA SCATTERING X-RAY FACTORS
    case ("-fe")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! input file
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = sfefile, fmt='(A)') buffer(1:len)
      nfe = 1
      
    case ("-rev")
      nfound = 1
      nrev = 1
      
    case ("-inf")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! input file format
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) nfin
      if (status/=0 .or. nfin<0) then
        call ExplainUsage()
        call CriticalError("Failed to read input file format selector.")
      end if
      
    case ("-ssc")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! single slice caclualtion index
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) ssc
      if (status/=0 .or. nfin<0) then
        call ExplainUsage()
        call CriticalError("Failed to read input file format selector.")
      end if
      
    case ("-buni")
      nfound = 1
      buni = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! single slice caclualtion index
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) buniv
      if (status/=0 .or. buniv<0.0) then
        call ExplainUsage()
        call CriticalError("Failed to read universal B_ISO parameter.")
      end if
    
    end select CHECK_COMMAND
    
    if (nfound == 0) then
      call ExplainUsage()
      call CriticalError("Command line parsing error. Unknown command ["//cmd(1:len)//"].")
    end if
    
    
  
  end do
  
! ------------
! treat unwanted combination of optional parameters -fl -dwf -abs
  if (nfl<0) nfl = 0
  if (nfl>1) nfl = 1
  if (nabs<0) nabs = 0
  if (nabs>1) nabs = 1
  if (nabf<0) nabf = 0
  if (nabf>1) nabf = 1
  if (ndwf<0) ndwf = 0
  if (ndwf>1) ndwf = 1
  if (npot<0) npot = 0
  if (npot>1) npot = 1
  if (n3dp<0) n3dp = 0
  if (n3dp>1) n3dp = 1
  if (nrev<0) nrev = 0
  if (nrev>1) nrev = 1
  if (nfl==1) then ! turn OFF dwf and abs
    ndwf = 0
    nabs = 0
  else
    nv = 1 ! set number of variants back to default 1.
    if (ndwf==0) then ! turn OFF abs
      nabs = 0
    end if
  end if
  if (nabs==1) then ! turn OFF abf
    nabf = 0
  end if
  if (nfe==1.and.nfx==1) then ! prefer fe-prm
    nfx = 0
  end if
  if (nfe==1.or.nfx==1) then ! prefer fe-prm
    if (nabs==1) then ! do not allow automatic absorption factors
      nabs = 0
      call PostWarning("Automatic absorption factor calculation is "// &
                     & "not supported with external scattering "// &
                     & "tables. The option -abs is ignored.")
    end if
  end if
  if (ssc<0) ssc = 0 ! turn single slice calculation off.
  
  
  if (ndwf==1) call PostMessage("Using Debye-Waller factors.")
  if (nabs==1) call PostMessage("Using absorption potentials.")
  if (nabf==1) then
    write(unit=smsg,fmt='(F8.3)') abf
    call PostMessage("Using absorption potentials with fix absorption factor "&
     &   //trim(adjustl(smsg))//".")
  else
    abf = 0.0
  end if
  if (nfl==1) call PostMessage("Generating frozen lattice configuration.")
  if (nfl==1.and.nv==1) call PostWarning(&
     & "Only one frozen lattice configuration is generated per slice.")
  if (nv>1) call PostMessage("Generating several variants per slice.")
  if (npot==1) call PostMessage("Exporting potentials to files *.pot.")
  if (n3dp==1) call PostMessage("Creating a 3D potential from the atomic structure model.")
  if (nfx==1) call PostMessage("Using external x-ray scattering factor tables.")
  if (nfe==1) call PostMessage("Using external electron scattering factor parameters.")
  if (nrev==0) call PostMessage("Sorting slices in default order from z/c=0 to z/c=1.")
  if (nrev==1) call PostMessage("Sorting slices in reverse order from z/c=1 to z/c=0.")
  if (nfin==0) call PostMessage("Input structure parameter file: default (cel).")
  if (nfin>=10.and.nfin<=19) call PostMessage("Input structure parameter file: 3D potential grid (asc).")
  if (ssc>0) call PostMessage("Single slice calculation mode.")
  write(unit=smsg,fmt='(F8.3)') buniv
  if (buni==1 .and. ndwf==1) call PostMessage("Using universal B_ISO parameter: "// &
     & trim(adjustl(smsg))//" nm^2.")

! ------------
! final option existence checks
  if (nprm==0) then
    call ExplainUsage()
    call CriticalError("Command line error. Structure data file not specified")
  end if
  if (nout==0) then
    call ExplainUsage()
    call CriticalError("Command line error. Output file not specified")
  end if
  if (nht==0) then
    call ExplainUsage()
    call CriticalError("Command line error. TEM high tension not specified")
  end if
  if (nfin==0) then
    if (nnx==0) then
      call ExplainUsage()
      call CriticalError("Command line error. x-discretization not specified")
    end if
    if (nny==0) then
      call ExplainUsage()
      call CriticalError("Command line error. y-discretization not specified")
    end if
    if (nnz==0) then
      call ExplainUsage()
      call CriticalError("Command line error. Number of slices not specified")
    end if
  end if

! ------------
  return

END SUBROUTINE ParseCommandLine


!**********************************************************************!
!
! function ht2wl
!
! returns electron wavelength for given high tension
!
real*4 function ht2wl(ht)
  implicit none
  real*4 :: ht
  ht2wl = 1.239842447 / sqrt( ht * ( 1022.0 + ht ) )
  return
end function ht2wl


!**********************************************************************!
!
! function writeslcprm
!
! writes slice file parameters to a text file
!
! parameter format:
! line 001: nz                          = number of slices defined below
! line 002: '[Slice Parameters]'        = slice parameter intro
! line 003: slc001                      = slice name string
! line 004: 'C:\slices\slc_001.sli'     = slice file name
! line 005: nx                          = slice x-discretization 
! line 006: ny                          = slice y-discretization 
! line 007: nv                          = number of structure variants
! line 008: sdx      	                = horizontal slice sampling [nm/pix]
! line 009: sdy    				        = vertical slice sampling [nm/pix]
! line 010: sdz                         = slice thickness [nm]
! line 011: parameter block for the next slice, analogous as from line 002
subroutine writeslcprm(sfile)

  implicit none
  
  include "global.fi"
  
  character(len=*) :: sfile
  
  integer*4 :: i, lun, slclun
  character(len=400) :: smsg, sslc, sslcn
  integer*4 :: getfreelun
  
  nerr = 0
  
  lun = getfreelun()
  if (lun<=0) call CriticalError("Failed to acquire free logical file unit.")
  
  open( file=trim(sfile),unit=lun,iostat=nerr, &
     &  action='WRITE',status='REPLACE')
  if (nerr/=0) call CriticalError("Failed to open file.")
  
  ! write number of slices
  write(unit=lun,fmt='(I6)') nz
  
  ! loop over all slices
  do i=1, nz
  
    ! write slice parameter intro
    write(unit=lun,fmt='(A)') "'[Slice Parameters]'"
    
    ! write slice name
    write(unit=lun,fmt='("slice",I<ndigsl>.<ndigsl>)') i
    
    ! get full name of slice file, and write it
    slclun = getfreelun()
    if (slclun<=0) call CriticalError("Failed to acquire free logical file unit.")
    write(unit=sslc,fmt='(A,"_",I<ndigsl>.<ndigsl>,".sli")') trim(sslcfile),i
    open( file=trim(sslc),unit=slclun,iostat=nerr,action='READ',status='OLD')
    if (nerr/=0) call CriticalError("Failed to open file.")
    inquire(unit=slclun, name=sslcn)
    close( unit=slclun, iostat=nerr )
    if (nerr/=0) call CriticalError("Failed to close file.")
    
    write(unit=lun,fmt='(A,A,A)') "'",trim(sslcn),"'"
    
    ! write slice dimensions x and y
    write(unit=lun,fmt='(I6)') nx
    write(unit=lun,fmt='(I6)') ny
    
    ! write slice number of variants
    write(unit=lun,fmt='(I6)') nv
    
    ! write slice sampling x and y
    write(unit=lun,fmt='(F12.8)') sdx
    write(unit=lun,fmt='(F12.8)') sdy
    
    ! write slice thickness
    write(unit=lun,fmt='(F12.8)') sdz
  
  end do
  
  close( unit=lun, iostat=nerr )
  if (nerr/=0) call CriticalError("Failed to close file.")
  
  return

end subroutine writeslcprm






!**********************************************************************!
!
! MAJOR SUBROUTINE
!
! Creates and saves phase gratings of cell slices
!
subroutine CEL2SLC()

  use CellSlicer
  use EMSdata

  implicit none
  
  include "global.fi"
  
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  character(len=16) :: sslnum, svrnum, svrmax
  real*4 :: fz0, fz1
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  
    
    ! generate scattering data
    call PostMessage("Preparing scattering amplitudes for all atom types.")
    CS_absorptionprm = abf
    call CS_PREPARE_SCATTAMPS(nx,ny,1,ndwf,nabs+2*nabf,wl,nerr)
    if (nerr/=0) call CriticalError("Failed to prepare scattering amplitudes")
 
    ! prepare slicing
    call PostMessage("Distributing atoms in slices.")
    call CS_SETSLICE_EQUIDIST(nz,nrev,nerr)
    if (nerr/=0) call CriticalError("Failed to distribute atoms in slices.")
    ! set absorption parameters
    CS_useabsorption = 0
    if (nabs/=0) CS_useabsorption = 1
  
    ! allocate slice memory
    call PostMessage("Allocating slice memory.")
    allocate(slcdat(nx,ny,nv),stat=nerr)
    if (nerr/=0) call CriticalError("Failed to allocate slice memory.")
  
    ! set potential backup option
    CS_backup_pot_flg = npot
  
    ! creating phase gratings for all slices
    do i=1, nz
    
      ! handle single slice calculation mode, added 2015-06-03, JB
      if (ssc>0 .and. i/=ssc) cycle ! skip all other slices in single slice calculation mode
    
      write(unit=sslnum,fmt='(I<ndigsl>.<ndigsl>)') i
      sslnum = trim(adjustl(sslnum))
    
      call PostMessage("Preparing slice "//trim(sslnum)//".")
      write(unit=smsg,fmt='("  ",I5," atoms.")') CS_slcatnum(i)
      call PostMessage(trim(smsg))
      !
      ! setup slice title using atomic content
      call CS_GETSLICETITLE(i,EMS_SLI_data_title,nerr)
      if (nerr/=0) call CriticalError("Failed analyse slice data.")
      call PostMessage("  "//trim(EMS_SLI_data_title))
      !
      ! calculate phase gratings
      CS_warn_num = 0
      do j=1, nv ! loop over numbr of variants (FL)
        
        write(unit=svrnum,fmt='(I<ndigvr>.<ndigvr>)') j
        svrnum = trim(adjustl(svrnum))
        write(unit=svrmax,fmt='(I<ndigvr>.<ndigvr>)') nv
        svrmax = trim(adjustl(svrmax))
        
        if (nv>1) then
          call PostMessage("- preparing variant "//trim(svrnum)// &
          &    " of "//trim(svrmax)//".")
        end if
        
        call CS_GETSLICE_PGR(i, nx, ny, 1, 1, nabs, nfl, ndwf, wl, &
           & slcdat(:,:,j), nerr)
        if (nerr/=0) call CriticalError("Slice preparation failed.")
      
        if (npot==1) then ! save potential backup
          if (nv>1) then
            sfil = trim(sslcfile)//"_"//trim(sslnum)//"_"// &
            &      trim(svrnum)//".pot"
          else
            sfil = trim(sslcfile)//"_"//trim(sslnum)//".pot"
          end if
          call savedatac8(trim(sfil), nx, ny, CS_backup_pot, nerr)
          if (nerr/=0) then
            call PostWarning("Failed to write potential to file.")
          else
            call PostMessage("Proj. slice potential saved to file ["// &
               & trim(sfil)//"].")
          end if
        end if
      end do
      !
      ! gather atom data to a table before saving (elastic data only, inelastic: TODO)
      nat = CS_slcatnum(i) ! number of atoms in this slice
      fz0 = CS_slczlim(1,i)/CS_scsz ! fractional slice z-offset
      fz1 = CS_slczlim(2,i)/CS_scsz ! fractional slice z-termination
      call EMS_SLI_settab1(fz0, fz1, nat, CS_slcatacc(1:nat,i), &
     &     CS_numat, CS_atnum(1:CS_numat), CS_atcrg(1:CS_numat), &
     &     CS_atdwf(1:CS_numat), CS_atocc(1:CS_numat), &
     &     CS_atpos(1:3,1:CS_numat) )
      !
      ! write data to a slice file
      sfil = trim(sslcfile)//"_"//trim(sslnum)//".sli"
      call EMS_SLI_save(trim(sfil),nx,ny,nv,CS_scsx,CS_scsy, &
                      & sdz,ht,slcdat,nerr)
      if (nerr/=0) then
        call CriticalError("Failed to write slice file.")
      else
        call PostMessage("Slice phase grating data saved to file ["// &
           & trim(sfil)//"].")
      end if
      !
    end do
    
  ! deallocate slice memory
  call PostMessage("Deallocating slice memory.")
  if (allocated(slcdat)) deallocate(slcdat,stat=nerr)
  if (nerr/=0) call CriticalError("Failed to deallocate slice memory.")
  
  ! write slicing parameters to file
  ! * modified 2015-06-03 for single slice calculation mode
  !   now the parameters are saved only in all-slice mode or in single slice mode, when the first slice is calculated
  if (ssc==0.or.ssc==1) then
    write(unit=sfil,fmt='(A,".prm")') trim(sslcfile)
    call PostMessage("Saving slice parameters to file ["//trim(sfil)//"].")
    call writeslcprm(trim(sfil))
    if (nerr/=0) call CriticalError("Failed to save slice parameter file.")
  end if
  
  return
  
end subroutine CEL2SLC





!**********************************************************************!
!
! MAJOR SUBROUTINE
!
! Creates 3D potentail and saves phase gratings of potential slices
!
subroutine CEL2POT3D2SLC()

  use CellSlicer
  use EMSdata
  use m3dpot

  implicit none
  
  include "global.fi"
  
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  real*4 :: fz0, fz1, relcor, fcorr
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  
    
    ! generate scattering data
    call PostMessage("Preparing scattering amplitudes for all atom types.")
    CS_absorptionprm = abf
    call CS_PREPARE_SCATTAMPS(nx,ny,nz,ndwf,nabs+2*nabf,wl,nerr)
    if (nerr/=0) call CriticalError("Failed to prepare scattering amplitudes")
    
    ! prepare slicing
    call PostMessage("Distributing atoms in slices.")
    call CS_SETSLICE_EQUIDIST(nz,nrev,nerr)
    if (nerr/=0) call CriticalError("Failed to distribute atoms in slices.")
 
    ! set absorption parameters
    CS_useabsorption = 0
    if (nabs/=0) CS_useabsorption = 1
    
    ! allocate potential memory
    if (allocated(M3D_pot)) deallocate(M3D_pot, stat=nerr)
    allocate(M3D_pot(nx,ny,nz), stat=nerr)
    if (nerr/=0) call CriticalError("Failed to allocate slice memory.")
    M3D_pot = cmplx(0.0,0.0)
    M3D_n1 = nx
    M3D_n2 = ny
    M3D_n3 = nz
    ! 
    M3D_b1 = (/ sdx , 0.0 , 0.0 /)
    M3D_b2 = (/ 0.0 , sdy , 0.0 /)
    M3D_b3 = (/ 0.0 , 0.0 , sdz /)
    !
    call CS_GETCELL_POT(nx, ny, nz, nfl, ndwf, wl, M3D_pot, nerr)
    if (nerr/=0) call CriticalError("Failed to create 3D potential.")
    ! dump the potential back to HD, remove the relativistic correction just for this
    fcorr = 1.0 / ( 1.0 + ht / 510.9989 )
    call savedata1( "pot3d-2.dat", 2*M3D_n1*M3D_n2*M3D_n3, M3D_pot*fcorr, nerr)
    
    
  !
  !
  ! --- make slices
  !
  call POT3D2SLC()
  !
  !
  
  
  return
  
end subroutine CEL2POT3D2SLC




!**********************************************************************!
!
! MAJOR SUBROUTINE
!
! Creates and saves phase gratings of 3d potential data
!
subroutine POT3D2SLC()

  use CellSlicer
  use EMSdata
  use m3dpot

  implicit none
  
  include "global.fi"
  
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  character(len=16) :: sslnum, svrnum, svrmax
  real*4 :: fz0, fz1
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  
    
    ! allocate slice memory
    call PostMessage("Allocating slice memory.")
    allocate(slcdat(nx,ny,1),stat=nerr)
    if (nerr/=0) call CriticalError("Failed to allocate slice memory.")
    
    ! set potential backup option
    M3D_backup_slcpot = npot
    
    ! creating phase gratings for all slices
    do i=1, nz
    
      ! handle single slice calculation mode, added 2015-06-03, JB
      if (ssc>0 .and. i/=ssc) cycle ! skip all other slices in single slice calculation mode
    
      write(unit=sslnum,fmt='(I<ndigsl>.<ndigsl>)') i
      sslnum = trim(adjustl(sslnum))
    
      call PostMessage("Preparing slice "//trim(sslnum)//".")
      !
      ! setup slice title from 3D-potential filename
      j = index(trim(scellfile),".",BACK=.TRUE.)
      if (j<1) j = len_trim(scellfile)+1
      k = index(trim(scellfile),"\",BACK=.TRUE.)
      if (k<1) k = index(trim(scellfile),"/",BACK=.TRUE.)
      if (k>j) then
        j = len_trim(scellfile)+1
        k = 0
      end if
      write(unit=stmp1, fmt='(I)') i
      write(unit=stmp2, fmt='(I)') nz
      EMS_SLI_data_title = scellfile(k+1:j-1)//" "// &
     &     trim(adjustl(stmp1))//"/"//trim(adjustl(stmp2))
      call PostMessage("  "//trim(EMS_SLI_data_title))
      !
      ! calculate phase gratings
      call M3D_getslice_pgr(i, nz, 1, 1, abf, ht, slcdat(:,:,1), nerr) 
      if (nerr/=0) call CriticalError("Slice preparation failed.")
    
      if (npot==1) then ! save potential backup
        sfil = trim(sslcfile)//"_"//trim(sslnum)//".pot"
        call savedatac8(trim(sfil), nx, ny, M3D_slcpot, nerr)
        if (nerr/=0) then
          call PostWarning("Failed to write potential to file.")
        else
          call PostMessage("Proj. slice potential saved to file ["// &
     &         trim(sfil)//"].")
        end if
      end if
      !
      !
      nat = 0
      if (allocated(CS_slcatnum)) then
        ! gather atom data to a table before saving (elastic data only, inelastic: TODO)
        nat = CS_slcatnum(i) ! number of atoms in this slice
      end if
      ! is there atom data?
      if (nat>0) then ! yes
        !
        ! Prepare the atom table
        fz0 = CS_slczlim(1,i)/CS_scsz ! fractional slice z-offset
        fz1 = CS_slczlim(2,i)/CS_scsz ! fractional slice z-termination
        call EMS_SLI_settab1(fz0, fz1, nat, CS_slcatacc(1:nat,i), &
     &     CS_numat, CS_atnum(1:CS_numat), CS_atcrg(1:CS_numat), &
     &     CS_atdwf(1:CS_numat), CS_atocc(1:CS_numat), &
     &     CS_atpos(1:3,1:CS_numat) )
        !
        ! Store a slice title made from the slice composition
        call CS_GETSLICETITLE(i,EMS_SLI_data_title,nerr)
        if (nerr/=0) call CriticalError("Failed analyse slice data.")
        call PostMessage("  "//trim(EMS_SLI_data_title))
        !
      else ! no
        ! clear atom data table, there is non with 3d potentials without loading a cell file
        call EMS_SLI_settab0()
      end if
      !
      ! write data to a slice file
      sfil = trim(sslcfile)//"_"//trim(sslnum)//".sli"
      call EMS_SLI_save(trim(sfil),nx,ny,nv,CS_scsx,CS_scsy, &
     &                  sdz,ht,slcdat,nerr)
      if (nerr/=0) then
        call CriticalError("Failed to write slice file.")
      else
        call PostMessage("Slice phase grating data saved to file ["// &
     &       trim(sfil)//"].")
      end if
      !
    end do
    
  ! deallocate slice memory
  call PostMessage("Deallocating slice memory.")
  if (allocated(slcdat)) deallocate(slcdat,stat=nerr)
  if (nerr/=0) call CriticalError("Failed to deallocate slice memory.")
  
  ! write slicing parameters to file
  ! * modified 2015-06-03 for single slice calculation mode
  !   now the parameters are saved only in all-slice mode or in single slice mode, when the first slice is calculated
  if (ssc==0.or.ssc==1) then
    write(unit=sfil,fmt='(A,".prm")') trim(sslcfile)
    call PostMessage("Saving slice parameters to file ["//trim(sfil)//"].")
    call writeslcprm(trim(sfil))
    if (nerr/=0) call CriticalError("Failed to save slice parameter file.")
  end if
  
  return
  
end subroutine POT3D2SLC