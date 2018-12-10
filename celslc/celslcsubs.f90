!**********************************************************************!
!**********************************************************************!
!
! FILE: "celslcsubs.f90"
!
! AUTHOR: Dr. J. Barthel
!         Ernst Ruska-Centre
!         Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
!
! PURPOSE: Implementation of subroutines for CELSLC
!
! VERSION: 0.70b, J.B., 28.11.2018
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
  use celslcprm  
  implicit none
  
  call PostMessage("")
  call PostMessage(" +---------------------------------------------------+")
  call PostMessage(" | Program [celslc]                                  |")
  call PostMessage(" | Version: 0.70b 64-bit  -  2018 Nov  28            |")
  call PostMessage(" | Author : Dr. J. Barthel, ju.barthel@fz-juelich.de |")
  call PostMessage(" |          Forschungszentrum Juelich GmbH, GERMANY  |")
  call PostMessage(" | License: GNU GPL 3 <http://www.gnu.org/licenses/> |")
  call PostMessage(" +---------------------------------------------------+")
  call PostMessage("")
  call PostMessage("")
  
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
  use celslcprm  
  implicit none
  
  call PostMessage("")
  call PostMessage("")
  
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
  use celslcprm
  implicit none
  
  character*(*) :: smessage

  write (unit=stdout,fmt='(A)') " "
  write (unit=stdout,fmt='(A)') trim(smessage)
  write (unit=stdout,fmt='(A)') "Critical error. Halting program."
  write (unit=stdout,fmt='(A)') " "
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
  use celslcprm
  implicit none
  
  character*(*) :: smessage

  if (nverbose>=0) write(unit=stdout,fmt='(A)') "Warning: "//trim(smessage)

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
  use celslcprm
  implicit none
  
  character*(*) :: smessage

  if (nverbose>=0) write(unit=stdout,fmt='(A)') trim(smessage)

  return

end subroutine PostMessage



!**********************************************************************!
!**********************************************************************!
subroutine PostRuntime(sinfo)
!
! Posts a run-time message with "sinfo" as add text and
!
! INPUT:
!   character(len=*) :: sinfo = add info

  use celslcprm
  
  implicit none

  character(len=*), intent(in) :: sinfo
  character(len=20) :: srt
  real*4 :: rt

  if (csprm_runtimes==1.and.nverbose>=0) then
    call csprm_getcls(rt)
    write(unit=srt,fmt='(G15.3)') rt
    call PostMessage("- "//trim(sinfo)//" ("//trim(adjustl(srt))//" s).")
  end if

  return

end subroutine PostRuntime




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
  use celslcprm
  use CellSlicer

  implicit none
  
  character(len=400) :: smsg
  
  ! start
  call PostMessage("Info on used super-cell:")
  
  ! determine sampling rates
  sdx = CS_scsx / real(nx)
  sdy = CS_scsy / real(ny)
  sdz = 0.0
  if (nz>0) sdz = CS_scsz / real(nz)
  
  write(unit=smsg,fmt='("Super-cell size (x,y,z) in nm: ",3F12.8)') CS_scsx, CS_scsy, CS_scsz 
  call PostMessage(trim(smsg))
  
  write(unit=smsg,fmt='("Super-cell sampling (dx,dy) in nm: ",2F12.8)') sdx, sdy
  call PostMessage(trim(smsg))
  
  if (nz>0) then
    write(unit=smsg,fmt=*) nz
    call PostMessage("User defined equidistant sampling along z with "// &
         & trim(adjustl(smsg))//" slices.")
  else if (nz==0) then
    call PostMessage("Automatic equidistant sampling along z.")
  else
    call PostMessage("Automatic variable sampling along z.")
  end if
  
  write(unit=smsg,fmt='("Total number of atoms in super-cell: ",I5)') CS_numat
  call PostMessage(trim(smsg))
    
  return

end subroutine PostCellInfo





!**********************************************************************!
!
! subroutine CheckCell
!
! Checks the current super cell for
! - strange volume
! - orthogonality
! - strange occupancies
! - strange thermal displacement parameter
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine CheckCell()

  use celslcprm
  use CellSlicer
  
  implicit none
  
  integer*4 :: i, j
  character(len=400) :: smsg
  real*4 :: tmp
  
  ! start
  call PostMessage("Checking input super-cell parameters:")
  
  ! volume
  call CS_GETCELL_VOL(tmp)
  write(unit=smsg,fmt=*) tmp
  call PostMessage("- super-cell volume = "//trim(adjustl(smsg))//" nm^3.")
  if (tmp<=0.0) then
    call CriticalError("Invalid super-cell volume.")
  end if
  
  ! orthogonality
  tmp = tmp / ( CS_scsx*CS_scsy*CS_scsz )
  write(unit=smsg,fmt=*) CS_scalpha
  call PostMessage("- super-cell angle alpha = "//trim(adjustl(smsg))//" deg.")
  write(unit=smsg,fmt=*) CS_scbeta
  call PostMessage("- super-cell angle beta  = "//trim(adjustl(smsg))//" deg.")
  write(unit=smsg,fmt=*) CS_scgamma
  call PostMessage("- super-cell angle gamma = "//trim(adjustl(smsg))//" deg.")
  CS_scorth_check = 0
  if (tmp<0.99) then
    call PostWarning("Non-Orthogonal input super-cell.")
    CS_scorth_check = -1
  else
    call PostMessage("- Input super-cell is orthogonal.")
    CS_scorth_check = 1
  end if
  
  ! strange occupancies
  if (CS_numat>0) then
    j = 0
    do i=1, CS_numat
      if (CS_atocc(i)<0.001) j = j + 1
    end do
    if (j>0) then
      write(unit=smsg,fmt=*) j
      call PostWarning("Found "//trim(adjustl(smsg))//" atoms with zero occupancy.")
    else
      call PostMessage("- occupancy test passed.")
    end if
  end if
  
  ! strange Biso
  if (CS_numat>0) then
    j = 0
    do i=1, CS_numat
      if (CS_atdwf(i)<0.00001) j = j + 1
    end do
    if (j>0) then
      write(unit=smsg,fmt=*) j
      call PostWarning("Found "//trim(adjustl(smsg))//" atoms with zero thermal displacement parameter.")
    else
      call PostMessage("- thermal displacement parameter test passed.")
    end if
  end if
  
end subroutine CheckCell








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
  use celslcprm
  implicit none
  
  call Introduce()
  call PostMessage(' CELSLC parameters')
  call PostMessage('  -cel <string> = CEL structure file (e.g. "atoms.cel")')
  call PostMessage('  -cif <string> = CIF structure file (e.g. "atoms.cif")')
  call PostMessage('  -ht <number>  = electron energy in keV (10 ... 1300)')
  call PostMessage('  -nx <number>  = horizontal cell sampling (32 ... 8192)')
  call PostMessage('  -ny <number>  = vertical cell sampling (32 ... 8192)')
  call PostMessage('  -slc <string> = slice file name (e.g. "scl")')
  call PostMessage(' [-abf <number> = apply fix rel. absorption factor]')
  call PostMessage(' [-abs         -> apply built-in absorption factors]')
  call PostMessage(' [-dwf         -> apply Debye-Waller factors]')
  call PostMessage(' [-fl          -> frozen-lattice simulation]')
  call PostMessage(' [-nv <number>  = number of FL variants per slice (1 ... 8192)]')
  call PostMessage(' [-nz <number>  = number of equidistant slices (0 ... 8192)]')
  call PostMessage(' [-pot         -> export potentials to files *.pot]')
  call PostMessage(' [-pps         -> output projected potentials to slice files]')
  call PostMessage(' [-prj <u,v,w,u,v,w,a,b,c> = super-cell re-orientation]')
  call PostMessage(' [-rev         -> slicing in reversed sequence]')
  call PostMessage(' [-ssc <number> = single slice calculation]')
  call PostMessage(' [-tla <x,y,z>  = additional shift of atoms]')
  call PostMessage('')
  call PostMessage(' [] = optional parameters')
  call Outroduce()
  
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
  use celslcprm
  implicit none
  
  character*512 :: buffer, cmd
  logical :: fex
  integer*4 :: lfu
  integer*4 :: i, cnt, status, len, nfound
  integer*4 :: nprm, nout, nht, nnx, nny
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
  nht = 0
  nfl = 0
  ndwf = 0
  nabs = 0
  nabf = 0
  npot = 0
  npps = 0
  nfx = 0
  nfe = 0
  nsca = 0
  nrev = 0
  nx = 0
  ny = 0
  nz = -1
  nv = 1
  nfin = 0 ! default input format = 0 = cel structure file
  abf = 0.0
  ssc = 0
  buni = 0
  buniv = 0.0
  block = 0
  bloh = 0.0
  blok = 0.0
  blol = 0.0
  blyh = 0.0
  blyk = 0.0
  blyl = 0.0
  blsa = 0.0
  blsb = 0.0
  blsc = 0.0
  ntla = 0
  tlax = 0.0
  tlay = 0.0
  tlaz = 0.0
  csprm_runtimes=0
  nffdec = 0
  nf2dec = 0

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
    case ("-cel") ! CEL input
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
      
    case ("-cif") ! CIF input
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
      nfin = 1
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
      if (status/=0 .or. nz>2048) then
        call PostWarning("Failed to read number of slices, using automatic slicing.")
        nz = 0 ! set back to default equidistant slicing
      end if
      if (nz<0) nz = -1 ! set back to default auto slicing
    
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
      
    ! FORM-FACTOR DECAY DATA OUTPUT
    case ("-ffdec")
      nfound = 1
      nffdec = 1
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
      read(unit=buffer,fmt=*,iostat=status) vffdec
      if (status/=0 .or. vffdec<=0.0 .or. vffdec>=1.0) then
        call ExplainUsage()
        call CriticalError("Invalid form-factor target decay (output declined).")
        nffdec = 0
      end if
      
    ! SCATTERIN-POWER LOSS OUTPUT
    case ("-f2dec")
      nfound = 1
      nf2dec = 1
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
      read(unit=buffer,fmt=*,iostat=status) vf2dec
      if (status/=0 .or. vffdec<=0.0 .or. vffdec>=1.0) then
        call ExplainUsage()
        call CriticalError("Invalid scattering power loss target (output declined).")
        nffdec = 0
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
    
    case ("-pps")  
    ! OUTPUT POTENTIALS TO SLICE FILES INSTEAD OF PHASE GRATINGS
      nfound = 1
      npps = 1
      
    ! CREATE 3D POTENTIAL
    case ("-3dp")
      nfound = 1
      n3dp = 1 ! switch on 3d potential creation
      
    ! USE EXTRA X-RAY SCATTERING FACTORS
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
      
    ! USE EXTRA ELECTRON SCATTERING FACTORS
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
    
    ! SLICE IN REVERSED ORDER 
    case ("-rev")
      nfound = 1
      nrev = 1
      
    case ("-rti")
      nfound = 1
      csprm_runtimes = 1
      
    case ("-verbose")
      nfound = 1
      nverbose = 3
    
    case ("-silent")
      nfound = 1
      nverbose = -1
    
    ! ALTERNATIVE INPUT FORMAT
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
    
    ! CALCULATE ONLY ONE SLICE  
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
      if (status/=0 .or. ssc<0) then
        call ExplainUsage()
        call CriticalError("Failed to read input file format selector.")
      end if
    
    ! APPLY UNIVERSAL THERMAL DISPLACEMENT PARAMETERS  
    case ("-buni")
      nfound = 1
      buni = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! universal Biso parameter
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) buniv
      if (status/=0 .or. buniv<0.0) then
        call ExplainUsage()
        call CriticalError("Failed to read universal B_ISO parameter.")
      end if
    
    ! CREATE A RE-ORIENTED ORTHOGONAL SUPER-CELL FROM INPUT STRUCTURE
    case ("-prj")
      nfound = 1
      block = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! h,k,l,m,n,o,a,b,c
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) bloh,blok,blol,blyh,blyk,blyl,blsa,blsb,blsc
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to read 9 orientation and size parameters.")
      end if
      
    ! SHIFT OF ALL ATOMS BY A CERTAIN VECTOR ALONG THE THREE AXES OF THE SUPER-CELL
    case ("-tla")
      nfound = 1
      ntla = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! x,y,z
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) tlax,tlay,tlaz
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to read 3 shift vector components.")
      end if
      
    ! THE INTERNAL FORM-FACTORS
    case ("-nsca")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status) ! form factor table index
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      read(unit=buffer,fmt=*,iostat=status) nsca
      if (status/=0 .or. nsca<0 .or. nsca>2 ) then
        call PostWarning("Failed to read form-factor table index. Switching to default.")
        nsca = 0 ! set back to default Weickenmeier and Kohl
      end if
      
    end select CHECK_COMMAND
    
    if (nfound == 0) then
      call ExplainUsage()
      call CriticalError("Command line parsing error. Unknown command ["//cmd(1:len)//"].")
    end if
  
  end do
  
  
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
  if (nfin==0.or.nfin==1) then
    if (nnx==0) then
      call ExplainUsage()
      call CriticalError("Command line error. x-discretization not specified")
    end if
    if (nny==0) then
      call ExplainUsage()
      call CriticalError("Command line error. y-discretization not specified")
    end if
  end if

! ------------
  return

END SUBROUTINE ParseCommandLine


SUBROUTINE CheckCommandLine
  use celslcprm
  use CellSlicer
  implicit none
  
  character(len=1024) smsg

  ! ------------
! treat unwanted combination of optional parameters
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
  if (npps<0) npps = 0
  if (npps>1) npps = 1
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
! This switch has become obsolete with the implementation of numerical
! integrations of absorptive form factors.
!  if (nfe==1.or.nfx==1) then ! prefer fe-prm
!    if (nabs==1) then ! do not allow automatic absorption factors
!      nabs = 0
!      call PostWarning("Automatic absorption factor calculation is "// &
!                     & "not supported with external scattering "// &
!                     & "tables. The option -abs is ignored. "// &
!                     & "Use the -abf option instead.")
!    end if
!  end if
  if (ssc<0) ssc = 0 ! turn single slice calculation off.
  
  if (nz<=0) then
    if (n3dp>1) then
      call CriticalError("Automatic slicing is not supported in "// &
                       & "combination with the option -3dp. Use the -nz option "// &
                       & "to explicitly set the sampling along z manually!")
    end if
    if (nfin>=10) then
      call PostMessage("Using automatic slicing according to the input 3D potential.")
    end if
    if (nz<0) then
      call PostMessage("Using non-equidistant automatic slicing of the super-cell.")
    else
      call PostMessage("Using equidistant automatic slicing of the super-cell.")
    end if
    
  else
    write(unit=smsg,fmt=*) nz
    call PostMessage( "Creating "//trim(adjustl(smsg))// &
                    & " equidistant slices of the super-cell.")
  end if
  if (ndwf==1) call PostMessage("Using Debye-Waller factors.")
  if (nabs==1) call PostMessage("Using absorption potentials.")
  if (nabf==1) then
    write(unit=smsg,fmt='(F8.3)') abf
    call PostMessage("Using absorption potentials with fix absorption factor "&
     &   //trim(adjustl(smsg))//".")
  else
    abf = 0.0
  end if
  if (nfl==1) call PostMessage("Generating frozen lattice configurations.")
  if (nfl==1.and.nv==1) call PostWarning(&
     & "Only one frozen lattice configuration is generated per slice.")
  if (nv>1) call PostMessage("Generating several variants per slice.")
  if (npot==1) call PostMessage("Exporting potentials to files *.pot.")
  if (npps==0) call PostMessage("Output in slice files are phase gratings.")
  if (npps==1) call PostMessage("Output in slice files are projected potentials.")
  if (n3dp==1) call PostMessage("Creating a 3D potential from the atomic structure model.")
  if (nfx==1) call PostMessage("Using external x-ray scattering factor tables.")
  if (nfe==1) call PostMessage("Using external electron scattering factor parameters.")
  if (nrev==0) call PostMessage("Sorting slices in default order from z/c=0 to z/c=1.")
  if (nrev==1) call PostMessage("Sorting slices in reverse order from z/c=1 to z/c=0.")
  if (nfin==0) call PostMessage("Input structure parameter file: CEL format (default).")
  if (nfin==1) call PostMessage("Input structure parameter file: CIF format.")
  if (nfin>=10.and.nfin<=19) call PostMessage("Input structure parameter file: ASCII 3D potential grid.")
  if (ssc>0) call PostMessage("Single slice calculation mode.")
  if (buni==1) then
    write(unit=smsg,fmt='(F8.3)') buniv
    call PostMessage("Using B_ISO = "//trim(adjustl(smsg))//" nm^2 for all atoms.")
  end if
  if (block==1) then
    call PostMessage("Re-orienting input structure to an orthogonal super-cell")
    write(unit=smsg,fmt='(F8.5,", ",F8.5,", ",F8.5)') blsa,blsb,blsc
    call PostMessage("- size: "//trim(adjustl(smsg))//" (nm)")
    write(unit=smsg,fmt='(F6.2,", ",F6.2,", ",F6.2)') bloh,blok,blol
    call PostMessage("- projection  axis: "//trim(adjustl(smsg))//" [uvw]")
    write(unit=smsg,fmt='(F6.2,", ",F6.2,", ",F6.2)') blyh,blyk,blyl
    call PostMessage("- projected y-axis: "//trim(adjustl(smsg))//" [uvw]")
  end if
  if (ntla==1) then
    call PostMessage("Shifting all atoms along the axes of the final orthogonal super-cell")
    write(unit=smsg,fmt='(G12.5,", ",G12.5,", ",G12.5)') tlax,tlay,tlaz
    call PostMessage("- shift vector: "//trim(adjustl(smsg))//" (fractional)")
  end if

  write(unit=smsg,fmt='(F8.2)') ht
  call PostMessage("Input electron energy: "//trim(adjustl(smsg))//" keV")
  write(unit=smsg,fmt='(G12.5)') wl*1000.0
  call PostMessage("Input electron wavelength: "//trim(adjustl(smsg))//" pm")
  write(unit=smsg,fmt='(G12.5)') CS_sig * wl
  call PostMessage("Interaction constant: "//trim(adjustl(smsg))//" (eV nm)^(-1) (2pi m0 e / h^2 * lambda)")
  if (nffdec==1) then
    write(unit=smsg,fmt=*) vffdec
    call PostMessage("Output of k-values where form factors decay by "// &
       & trim(adjustl(smsg))//" to file ffdec.txt")
  end if
  if (nf2dec==1) then
    write(unit=smsg,fmt=*) vf2dec
    call PostMessage("Output of k-values where loss of scattering power exceeds "// &
       & trim(adjustl(smsg))//" to file f2dec.txt")
  end if

  return

END SUBROUTINE CheckCommandLine


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
  
  use celslcprm
  use CellSlicer
  
  implicit none
  
  character(len=*) :: sfile
  
  integer*4 :: i, lun, slclun
  character(len=400) :: smsg, sslc, sslcn
  integer*4, external :: getfreelun
  external :: createfilefolder
  
  nerr = 0
  
  lun = getfreelun()
  if (lun<=0) call CriticalError("Failed to acquire free logical file unit.")
  call createfilefolder(trim(sfile),nerr)
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
    if (nerr/=0) call CriticalError("Failed to open slice file.")
    inquire(unit=slclun, name=sslcn)
    close( unit=slclun, iostat=nerr )
    if (nerr/=0) call CriticalError("Failed to close slice file.")
    
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
    sdz = CS_slczlim(2,i) - CS_slczlim(1,i)
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
  use celslcprm
  use CellSlicer
  use EMSdata

  implicit none
  
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  character(len=16) :: sslnum, svrnum, svrmax, snatoms, sthick
  real*4 :: fz0, fz1, mpot
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  real*4, external :: HT2WL
  
    
  ! generate scattering data
  call PostMessage("Preparing scattering amplitudes for all atom types.")
  CS_absorptionprm = abf
  CS_useppot = 1 ! use projected potentials
  call CS_PREPARE_SCATTAMPS(nx,ny,1,ndwf,nabs+2*nabf,wl,nerr)
  if (nerr/=0) call CriticalError("Failed to prepare scattering amplitudes")
  call CS_GET_MEANINNERPOT(mpot, ht, nerr)
  write(unit=smsg,fmt='(F10.5)') mpot
  call PostMessage("- mean inner potential: "//trim(adjustl(smsg))//" V")
  wl1 = HT2WL(ht+0.001*mpot)
  write(unit=smsg,fmt='(F10.5)') wl1*1000.
  call PostMessage("- wavelength corrected for refraction: "//trim(adjustl(smsg))//" pm")
  write(unit=smsg,fmt='(F10.5)') wl*1000.
  call PostMessage("- wavelength in vacuum               : "//trim(adjustl(smsg))//" pm")
 
  ! prepare slicing
  call PostMessage("Distributing atoms in slices.")
  if (nz>0) then ! equidistant slicing
    call CS_SETSLICE_EQUIDIST(nz,nrev,nerr)
  else if (nz==0) then ! equidistant auto slicing
    call CS_SUGGEST_NSLCEQUI(0.02,nz,nerr)
    call CS_SETSLICE_EQUIDIST(nz,nrev,nerr)
  else ! non-equidistant auto slicing
    call CS_SETSLICE_AUTO(nrev,nerr)
    nz = CS_nspsc ! update nz to number of slices
  end if
  if (nerr/=0) then ! slicing routine failure
    call CriticalError("Failed to distribute atoms in slices.")
  end if
  if (nz>0) then ! report number of slices
    write(unit=smsg, fmt=*) nz
    call PostMessage("- generated "//trim(adjustl(smsg))//" partitions of the super-cell.")
  else ! no slices, this is not allowed ... critical error
    call CriticalError("Invalid number of slices (0).")
  end if
    
  ! check single slice calculation index
  ssc = min(ssc,nz) ! limit slice index to number of slices
  if (ssc>0) then
    write(unit=stmp1, fmt=*) ssc ! smsg still contains the string of nz
    call PostMessage("- calculation restricted to slice #"// &
          & trim(adjustl(stmp1))//" of "//trim(adjustl(smsg))//".")
  end if
    
  ! set absorption parameters in calculation module
  CS_useabsorption = 0
  if (nabs/=0) CS_useabsorption = 1
  
  ! allocate slice memory
  call PostMessage("Allocating slice memory.")
  allocate(slcdat(nx,ny,nv),stat=nerr)
  if (nerr/=0) call CriticalError("Failed to allocate slice memory.")
  
  ! set potential backup option
  CS_backup_pot_flg = npot ! store potentials, we want to export them
  if (npps==1) CS_backup_pot_flg = 1 ! store potentials, we want to output them
  
  ! creating phase gratings for all slices
  CS_maxphase = 0.0
  call PostMessage("Calculating phase gratings ...")
  do i=1, nz
    
    ! handle single slice calculation mode, added 2015-06-03, JB
    if (ssc>0 .and. i/=ssc) cycle ! skip all other slices in single slice calculation mode
    
    write(unit=sslnum,fmt='(I<ndigsl>.<ndigsl>)') i
    sslnum = trim(adjustl(sslnum))
      
    write(unit=snatoms,fmt=*) CS_slcatnum(i)
    write(unit=sthick,fmt='(F10.3)') CS_slczlim(i,2)-CS_slczlim(i,1)
    call PostMessage("- slice #"//trim(sslnum)//", "// &
        & trim(adjustl(snatoms))//" atoms, "//trim(adjustl(sthick))// &
        & " nm thick")
    !
    ! setup slice title using atomic content
    call CS_GETSLICETITLE(i,EMS_SLI_data_title,nerr)
    if (nerr/=0) call CriticalError("Failed analyse slice data.")
    call PostMessage("  content: "//trim(EMS_SLI_data_title))
    !
    ! calculate phase gratings
    CS_warn_num = 0
    if (nv>1) then
      ! Start non-advancing progress display for variant calculations.
      ! Make sure to only write error outputs in the enclosed routines
      OPEN (UNIT=stdout,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')
    end if
    do j=1, nv ! loop over numbr of variants (FL)
        
      write(unit=svrnum,fmt='(I<ndigvr>.<ndigvr>)') j
      svrnum = trim(adjustl(svrnum))
      write(unit=svrmax,fmt='(I<ndigvr>.<ndigvr>)') nv
      svrmax = trim(adjustl(svrmax))
        
      if (nv>1) then
        ! update non-advancing
        write(unit=stdout,fmt='("+",A)') "  variant # "//trim(svrnum)// &
        &    " / "//trim(svrmax)//" "
      end if
        
      call CS_GETSLICE_PGR(i, nx, ny, 1, 1, nabs, nfl, ndwf, wl, &
          & slcdat(:,:,j), nerr)
      if (nerr/=0) call CriticalError("Slice preparation failed.")
      EMS_SLI_data_ctype = 0 ! set default slice data type (phase grating)
        
      if (npps==1) then ! replace pgr output with potential data
        slcdat(1:nx,1:ny,j) = CS_backup_pot(1:nx,1:ny)
        ! mark the change in the slice export module
        EMS_SLI_data_ctype = 1 ! set to potential
      end if
      
      if (npot==1) then ! save potential backup
        if (nv>1) then
          sfil = trim(sslcfile)//"_"//trim(sslnum)//"_"// &
          &      trim(svrnum)//".pot"
        else
          sfil = trim(sslcfile)//"_"//trim(sslnum)//".pot"
        end if
        call savedatac8(trim(sfil), nx* ny, CS_backup_pot, nerr)
        !if (nerr/=0) then
        !  call PostWarning("Failed to write potential to file.")
        !else
        !  call PostMessage("Proj. slice potential saved to file ["// &
        !     & trim(sfil)//"].")
        !end if
      end if
        
    end do ! loop over variants
    if (nv>1) then
      ! Stop non-advancing output and return to normal
      OPEN (UNIT=stdout,FORM='FORMATTED',CARRIAGECONTROL='LIST')
    end if
    !
    ! gather atom data to a table before saving (elastic data only, inelastic: TODO)
    nat = CS_slcatnum(i) ! number of atoms in this slice
    fz0 = CS_slczlim(1,i)/CS_scsz ! fractional slice z-offset
    fz1 = CS_slczlim(2,i)/CS_scsz ! fractional slice z-termination
    sdz = CS_slczlim(2,i) - CS_slczlim(1,i) ! slice thickness in nm
    !
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
      if (npps==1) then
        call PostMessage("  projected potential data saved to file ["// &
          & trim(sfil)//"].")
      else
        call PostMessage("  slice phase-grating data saved to file ["// &
          & trim(sfil)//"].")
      end if
    end if
    !
  end do ! loop over slices
  
  if (CS_maxphase > CS_pi4) then
    ! post strong phase shift warning
    write(unit=smsg,fmt='(G13.4)') CS_maxphase
    call PostWarning("Max. phase shift exceeds pi/4 ("//trim(adjustl(smsg))//").")
  end if
  
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
! Creates 3D potential and saves phase gratings of potential slices
!
subroutine CEL2POT3D2SLC()
  use celslcprm
  use CellSlicer
  use EMSdata
  use m3dpot

  implicit none
  
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  real*4 :: fz0, fz1, relcor, fcorr, mpot
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  real*4, external :: HT2WL
  
    
  ! generate scattering data
  CS_useppot = 0 ! use 3d potentials and not projected potentials
  call PostMessage("Preparing scattering amplitudes for all atom types.")
  CS_absorptionprm = abf
  call CS_PREPARE_SCATTAMPS(nx,ny,nz,ndwf,nabs+2*nabf,wl,nerr)
  if (nerr/=0) call CriticalError("Failed to prepare scattering amplitudes")
  call CS_GET_MEANINNERPOT(mpot, ht, nerr)
  write(unit=smsg,fmt='(F10.5)') mpot
  call PostMessage("- mean inner potential: "//trim(adjustl(smsg))//" V")
  wl1 = HT2WL(ht+0.001*mpot)
  write(unit=smsg,fmt='(F10.5)') wl1*1000.
  call PostMessage("- wavelength corrected for refraction: "//trim(adjustl(smsg))//" pm")
  write(unit=smsg,fmt='(F10.5)') wl*1000.
  call PostMessage("- wavelength in vacuum               : "//trim(adjustl(smsg))//" pm")
    
  ! prepare slicing
  call PostMessage("Distributing atoms in slices.")
  call CS_SETSLICE_EQUIDIST(nz,nrev,nerr) ! only equidistant slicing for 3D potentials
  if (nerr/=0) then
    call CriticalError("Failed to distribute atoms in slices.")
  end if
  if (nz>0) then ! report number of slices
    write(unit=smsg, fmt=*) nz
    call PostMessage("- generated "//trim(adjustl(smsg))//" partitions of the super-cell.")
  else ! no slices, this is not allowed ... critical error
    call CriticalError("Invalid number of slices (0).")
  end if
    
  ! check single slice calculation index
  ssc = min(ssc,nz) ! limit slice index to number of slices
  if (ssc>0) then
    write(unit=stmp1, fmt=*) ssc ! smsg still contains the string of nz
    call PostMessage("- calculation restricted to slice #"// &
          & trim(adjustl(stmp1))//" of "//trim(adjustl(smsg))//".")
  end if
 
  ! set absorption parameters
  CS_useabsorption = 0
  if (nabs/=0) CS_useabsorption = 1
  
  
  call PostMessage("Calculating 3D potential ...")
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
        
  allocate(slcdat(M3D_n1,M3D_n2,M3D_n3),stat=nerr)
  if (nerr==0) then
    ! dump the potential back to HD, remove the relativistic correction just for this
    fcorr = 1.0 / ( 1.0 + ht / 510.9989 )
    slcdat = M3D_pot*fcorr
    call savedata1( "pot3d-2.dat", 2*M3D_n1*M3D_n2*M3D_n3, slcdat, nerr)
    deallocate(slcdat,stat=nerr)
  end if
    
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
  use celslcprm
  use CellSlicer
  use EMSdata
  use m3dpot

  implicit none
  
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  character(len=16) :: sslnum, svrnum, svrmax, sthick
  real*4 :: fz0, fz1, fcorr
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  
    
  ! allocate slice memory
  call PostMessage("Allocating slice memory.")
  allocate(slcdat(nx,ny,1),stat=nerr)
  if (nerr/=0) call CriticalError("Failed to allocate slice memory.")
  CS_useppot = 0 ! use 3d potentials and not projected potentials
                  ! (not really required here as we use external potentials)
    
  ! set potential backup option
  M3D_backup_slcpot = npot
  fcorr = 1.0 / ( 1.0 + ht / 510.998928 ) ! relativistic correction, removal factor
    
  ! creating phase gratings for all slices
  CS_maxphase = 0.0
  call PostMessage("Calculating phase gratings ...")
  do i=1, nz
    
    ! handle single slice calculation mode, added 2015-06-03, JB
    if (ssc>0 .and. i/=ssc) cycle ! skip all other slices in single slice calculation mode
    
    write(unit=sslnum,fmt='(I<ndigsl>.<ndigsl>)') i
    sslnum = trim(adjustl(sslnum))
    write(unit=sthick,fmt='(F10.3)') CS_slczlim(i,2)-CS_slczlim(i,1)
    
    call PostMessage("- slice #"//trim(sslnum)//", "// &
        & trim(adjustl(sthick))//" nm thick")
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
    call PostMessage("  name: "//trim(EMS_SLI_data_title))
    !
    ! calculate phase gratings
    call M3D_getslice_pgr(i, nz, 1, 1, abf, ht, slcdat(:,:,1), nerr) 
    if (nerr/=0) call CriticalError("Slice preparation failed.")
    EMS_SLI_data_ctype = 0 ! set default slice data type (phase grating)
      
    if (npps==1) then ! replace the phase grating data with the stored potential
      slcdat(1:nx,1:ny,1) = M3D_slcpot(1:nx,1:ny)*fcorr ! rel. corr. removed
      EMS_SLI_data_ctype = 1 ! update the slice data type to potential
    end if
    
    if (npot==1) then ! save potential backup
      sfil = trim(sslcfile)//"_"//trim(sslnum)//".pot"
      call savedatac8(trim(sfil), nx* ny, M3D_slcpot*fcorr, nerr)
      !if (nerr/=0) then
      !  call PostWarning("Failed to write potential to file.")
      !else
      !  call PostMessage("Proj. slice potential saved to file ["// &
      !     & trim(sfil)//"].")
      !end if
    end if
    !
    !
    nat = 0
    if (allocated(CS_slcatnum)) then
      ! gather atom data to a table before saving (elastic data only, inelastic: TODO)
      nat = CS_slcatnum(i) ! number of atoms in this slice
    end if
    fz0 = CS_slczlim(1,i)/CS_scsz ! fractional slice z-offset
    fz1 = CS_slczlim(2,i)/CS_scsz ! fractional slice z-termination
    sdz = CS_slczlim(2,i) - CS_slczlim(1,i) ! slice thickness
    ! is there atom data?
    if (nat>0) then ! yes
      !
      ! Prepare the atom table
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
      if (npps==1) then
        call PostMessage("  projected potential data saved to file ["// &
          & trim(sfil)//"].")
      else
        call PostMessage("  slice phase-grating data saved to file ["// &
          & trim(sfil)//"].")
      end if
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