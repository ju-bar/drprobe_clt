!**********************************************************************!
!**********************************************************************!
!
! FILE: "wavimgsubs.f90"
!
! AUTHOR: Dr. J. Barthel
!         Forschungszentrum Jülich
!         Jülich, Germany
!
! PURPOSE: Implementations for image simulations
!
! VERSION: 0.70, J.B., 14.09.2018
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

  use wavimgprm  
  implicit none
  
  call PostMessage("")
  call PostMessage(" +---------------------------------------------------+")
  call PostMessage(" | Program: [wavimg]                                 |")
  call PostMessage(" | Version: 0.70b 64-bit  -  2018 Sept 14  -         |")
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

  use wavimgprm  
  implicit none
  
  character(len=400) :: smsg
  
  write(unit=smsg,fmt='(A,I3,A,I3,A)') "Number of warnings: ",numw,",  number of errors: ",nume,"."
  call PostMessage("")
  call PostMessage(trim(smsg))
  call PostMessage("")
  
  return

end subroutine Outroduce


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

  call Introduce()
  call PostMessage("Usage of wavimg in command line:")
  call PostMessage("wavimg  -prm 'parameter file name', e.g. 'wavimg.prm', MUST BE GIVEN")
  call PostMessage("       [-out 'image file name', e.g. 'img\image.dat', optional]")
  call PostMessage("       [-foc <float> - image defocus in nm, optional]")
  call PostMessage("       [-btx <float> - beam tilt x in mrad, optional]")
  call PostMessage("       [-bty <float> - beam tilt y in mrad, optional]")
  call PostMessage("       [-oar <float> - objective aperture radius in mrad, optional]")
  call PostMessage("       [/sil deactivates console output, switch]")
  call PostMessage("       [/dbg activates extra debug console output, switch]")
  call PostMessage("       [/nli no loop image output in map output mode, switch]")
  call Outroduce()

  return

end subroutine ExplainUsage





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

  use wavimgprm
  implicit none
  
  character*(*) :: smessage

  nume = nume + 1
  
  write(unit=*,fmt='(A)') " "
  write(unit=*,fmt='(A)') "Error: "//trim(smessage)
  write(unit=*,fmt='(A)') "Critical error. Halting program."
  write(unit=*,fmt='(A)') " "
  call Outroduce()
  stop

  return

end subroutine CriticalError
!**********************************************************************!


!**********************************************************************!
!
! subroutine PostError
!
! posts an error message to console, continues program
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostError(smessage)

  use wavimgprm
  implicit none
  
  character*(*) :: smessage

  nume = nume + 1

  write(unit=*,fmt='(A)') "Error: "//trim(smessage)

  return

end subroutine PostError
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

  use wavimgprm
  implicit none
  
  character*(*) :: smessage
  
  numw = numw + 1

  write(unit=*,fmt='(A)') "Warning: "//trim(smessage)

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

  use wavimgprm
  implicit none
  
  character*(*) :: smessage

  if (nsil==0) then ! do not post if silence flag is set in global options
    write (unit=*,fmt='(A)') trim(smessage)
  end if

  return

end subroutine PostMessage


!**********************************************************************!
!
! subroutine PostDBGMessage
!
! posts a normal message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostDBGMessage(smessage)

  use wavimgprm
  implicit none
  
  character*(*) :: smessage

  if (ndbg/=0 .and. nsil==0) then ! post if debug flag is set in global options
    write (*,*) "debug: "//trim(smessage)
  end if

  return

end subroutine PostDBGMessage



!**********************************************************************!
!
! subroutine PostRuntime
!
! posts a normal message with runtime information to console.
!
! INPUT:
!   character(len=*) :: smessage    = the runtime stamp info as string
!
! IN/OUTPUT: none
!
subroutine PostRuntime(smessage)

  use wavimgprm
  implicit none
  
  character*(*) :: smessage
  real*4 :: rt

   ! do not post if silence flag is set in global options
   ! and if runtime info is requested
  if (nsil==0 .and. wicl_runtimes==1) then
    call wavimgprm_getcls(rt)
    write(unit=stmp,fmt='(G15.3)') rt
    call PostMessage("- "//trim(smessage)//" ("//trim(adjustl(stmp))//" s).")
  end if

  return

end subroutine PostRuntime


!**********************************************************************!
!
! subroutine CheckLicense
!
! quietly check license
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine CheckLicense()

  implicit none
  
  integer*4, parameter :: noff = 46
  integer*4, parameter :: n3 = 62
  integer*4, parameter :: n5 = 70
  integer*4, parameter :: n2 = 71
  integer*4, parameter :: n4 = 0
  integer*4, parameter :: n1 = 52
  integer*4, parameter :: n7 = 67 
  integer*4, parameter :: n6 = 69
  
  logical :: fex
  
  character(len=100) :: f0
  character(len=7) :: f1
  character(len=400) :: f2
  
  f0 = "C:\"
  f1 = "       "
  f0 = trim(f0)//"windows\"
  write(unit=f1,fmt='(7a1)') char(n1+noff),char(n2+noff),char(n3+noff),&
     &         char(n4+noff),char(n5+noff),char(n6+noff),char(n7+noff)
  f0 = trim(f0)//"system32\"
  f2 = trim(f0)//trim(f1)
  
  inquire(file=trim(f2),exist=fex)
  
  if (.not.fex) stop

  return
  
end subroutine CheckLicense










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

  implicit none

! ------------
! declaration
  real*4 :: sigmoid
  real*4, intent(in) :: x, x0, dx
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

  use wavimgprm
  implicit none
  
  character*512 :: buffer, cmd
  logical :: fex
  integer*4 :: i, cnt, status, len, nfound
  integer*4 :: nprm
  integer*4, external :: getfreelun
  real*4, external :: HT2WL

! ------------
! initialize
!  write(unit=*,fmt=*) " > ParseCommandLine: INIT."
  i = 0
  cnt = command_argument_count()
  if (cnt==0) then
    call ExplainUsage()
    call CriticalError("Missing required arguments.")
  end if
  nprm = 0
  nfoc_ex = 0
  oapr_ex = -1.0
  simgfile_ex = ""
  sbshx = 0.0
  sbshy = 0.0
  dornsb = 0
  wicl_runtimes = 0

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
    
    ! THE PARAMETER FILE
    case ("-prm")
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
      write(unit = sprmfile, fmt='(A)') buffer(1:len)
      inquire(file=trim(sprmfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid argument: Specified parameter file does not exist.")
      end if
      nprm = 1
    
    ! THE IMAGE OUTPUT FILE NAME
    case ("-out")
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
      write(unit = simgfile_ex, fmt='(A)') buffer(1:len)
      
    ! THE OBJECTIVE APERTURE RADIUS
    case ("-oar")
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
      read(unit=buffer, fmt=*, iostat=status) oapr_ex
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to recognize objective aperture radius.")
      end if
      
    case ("-btx")
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
      read(unit=buffer, fmt=*, iostat=status) btx
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to recognize beam tilt x.")
      end if
      
    case ("-bty")
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
      read(unit=buffer, fmt=*, iostat=status) bty
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to recognize beam tilt y.")
      end if
      
    case ("-sbshx")
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
      read(unit=buffer, fmt=*, iostat=status) sbshx
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to recognize sideband shift x.")
      end if
      
    case ("-sbshy")
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
      read(unit=buffer, fmt=*, iostat=status) sbshy
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to recognize sideband shift y.")
      end if
      
    case ("-foc")
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
      read(unit=buffer, fmt=*, iostat=status) foc_ex
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Failed to recognize focus value.")
      end if
      nfoc_ex = 1
      
    case ("/rnsb")
      nfound = 1
      dornsb = 1
      

    case ("/sil")
      nfound = 1
      nsil = 1
      
    case ("/dbg")
      nfound = 1
      ndbg = 1
      
    case ("/nli")
      nfound = 1
      nnli = 1
      
    case ("/rti")
      nfound = 1
      wicl_runtimes = 1
    
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
    call CriticalError("Command line error. Parameter file not specified")
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
! subroutine loadprm
!
! loads the program parameters from file
!
! Format of the parameter file: TEXT ASCII
! Version 0.55b+
!
! line 01: '<wave name>'                ! Wave function file name string used to locate existing wave functions. Use quotation marks to secure the input including space characters.
! line 02: <nx>, <ny>                   ! Dimension of the wave data in pixels, <nx> = number of horizontal wave pixels, <ny>  = number of vertical wave pixels.
! line 03: <sx>, <sy>                   ! Sampling rate of the wave data (<sx> = horizontal, <sy> = vertical) [nm/pix].
! line 04: <HT>                         ! TEM high-tension as used for wave function calculation [kV].
! line 05: <notype>                     ! Image output type option: 0 = TEM image, 1 = complex image plane wave, 2 = wave amplitude, 3 = wave phase, 4 = wave real part, 5 = wave imaginary part, 6 = TEM image map of 2 variables
! line 06: '<image name>'               ! Image output file name string. Use quotation marks to secure the input including space characters.
! line 07: <ix>, <iy>                   ! Image output size (<ix> = horizontal , <iy> = vertical) in number of pixels.
! line 08: <intflg>, <mean>, <conv>, <rnoise> ! Flag and parameters for creating integer images with optional noise. Flag <intflg> 0 = off (default), 1 = 32-bit, 2 = 16-bit, Parameter: <mean> = mean vacuum intensity, <conv> = electron to counts conversion rate, <rnoise> detector readout noise level.
! line 09: <extfgl>                     ! Flag activating the extraction of a special image frame (0=OFF, 1=ON). The frame parameters are defined in the lines below.
! line 10: <so>                         ! Image output sampling rate [nm/pix], isotropic. The parameter is used only if the Flag in line 09 is set to 1.
! line 11: <opx>, <opy>                 ! Image frame offset in pixels of the input wave. The parameter is used only if the Flag in line 09 is set to 1.
! line 12: <orot>                       ! Image frame rotation in [deg] with respect to the input wave horizontal axis. The parameter is used only if the Flag in line 09 is set to 1.
! line 13: <coh-model>                  ! Coherence calculation model switch: 1 = averaging of coherent sub images explicit focal variation but quasi-coherent spatial envelope, 2 = averaging of coherent sub images with explicit focal and angular variation, 3 = quasi-coherent linear envelopes, 4 = Fourier-space synthesis with partially coherent TCC, 5: averaging of coherent sub images with explicit focal, angular, and frozen lattice variation)
! line 14: <ptcflg>, <f-spread>         ! Flag and parameters for partial temporal coherence: <ptcflg> = flag (0=OFF, 1=ON), <f-spread> = focus spread (1/e) half width [nm]
! line 15: <pscflg>, <s-conv>           ! Flag and parameters for partial spatial coherence: <pscflg> = flag (0=OFF, 1=ON), <s-conv> = beam convergence (1/e) half width [mrad]
! line 16: <mtfflg>, <mtf-scale>, '<mtf-file>' ! Flag and parameters for applying the detector MTF: <mtfflag> = flag (0=OFF, 1=ON), <mtf-scale> = calculation scale of the mtf = (sampling rate experiment)/(sampling rate simulation), <mtf-file> = File name string to locate the MTF data. Use quotation marks to secure the input including space characters.
! line 17: <vibflg>, <vibprm1>, <vibprm2>, <vibprm3> ! Flag and parameters for a vibration envelope: <vibflg> = flag (0=OFF, 1=ON-ISO, 2=ON-ANISO), <vibprm1>, <vibprm1> = vibration RMS amplitudes [nm], <vibprm3> = orientation [deg] of the primary vibration amplitude w.r.t. the horizontal image axis.
! line 18: <numabr>                     ! Number of aberration definitions following this line
! line 19: <aidx>, <ax>, <ay>           ! aberration definition (index, ax, ay) [nm]
! line 20: <aidx>, <ax>, <ay>           ! aberration definition (index, ax, ay) [nm] .. more similar lines possible depending on value of line 18
! line 21: <oap>                        ! Objective aperture radius [mrad]. Set to very large values to deactivate.
! line 22: <ocx>, <oxy>                 ! Center of the objective aperture with respect to the zero beam [mrad].
! line 23: <numloop>                    ! Number variable of loop definitions following below.
! line 24: <loop-class>                 ! Loop variable class: 1 = aberration, 2 = coherence params, 3 = wave file
! line 25: <loop-id>                    ! Loop variable identifier (e.g. aberration index)
! line 26: <loop-form>                  ! Loop variation form: 1 = ramp, 2 = oscillate, 3 = listed
! line 27: <range-min>, <range-max>, <num-steps> ! Loop variation range: [<loop-form> == 1] <range-min> = minimum variable value, <range-max> = maximum vaiable value, <num-steps> = number of loop steps (samples), [<loop-form> == 2] <range-min> = variation amplitude, <range-max> = variation frequency, <num-steps> = number of loop steps (samples)
! line 28: <loop-string>                ! Loop identifier string used for file naming or to identify the suffix sub-string for the location of numbers in input file names.
! More loop definitions are possible according to the number of loops define din line 23
!
subroutine loadprm

  use wavimgprm
  use AberrationFunctions
  
  implicit none
  
  integer*4 :: lun, na, ai, i
  real*4 :: wax, way
  integer*4, external :: getfreelun
  character(len=1024) :: smsg, sprm1, sprm2, sprm3
  
  call PostMessage("Loading parameters from file ["//trim(sprmfile)//"].")
  
  lun = getfreelun()
  if (lun<=0) goto 13
  
  open( unit=lun, file=trim(sprmfile), iostat=nerr, &
     &  action='read', status='old', err=14)
  
  sdbgmsg ="(line 1) swavfile"
  read(unit=lun, fmt=*, err=16) swavfile
  sdbgmsg ="(line 2) nwx, nwy"
  read(unit=lun, fmt=*, err=16) nwx, nwy
  sdbgmsg ="(line 3) swx, swy"
  read(unit=lun, fmt=*, err=16) swx, swy
  sdbgmsg ="(line 4) ht"
  read(unit=lun, fmt=*, err=16) ht
  sdbgmsg ="(line 5) notype"
  read(unit=lun, fmt=*, err=16) notype
  sdbgmsg ="(line 6) simgfile"
  read(unit=lun, fmt=*, err=16) simgfile
  sdbgmsg ="(line 7) nx, ny"
  read(unit=lun, fmt=*, err=16) nx, ny
  sdbgmsg ="(line 8) doint, iimean, el_conv, dark_noise"
  read(unit=lun, fmt=*, err=16) doint, iimean, el_conv, dark_noise
  sdbgmsg ="(line 9) dofrm"
  read(unit=lun, fmt=*, err=16) dofrm
  sdbgmsg ="(line 10) simg"
  read(unit=lun, fmt=*, err=16) simg
  sdbgmsg ="(line 11) ox, oy"
  read(unit=lun, fmt=*, err=16) ox, oy
  sdbgmsg ="(line 12) ax"
  read(unit=lun, fmt=*, err=16) ax
  sdbgmsg ="(line 13) ncohm"
  read(unit=lun, fmt=*, err=16) ncohm
  sdbgmsg ="(line 14) doptc, fs"
  read(unit=lun, fmt=*, err=16) doptc, fs 
  sdbgmsg ="(line 15) dopsc, sc"
  read(unit=lun, fmt=*, err=16) dopsc, sc 
  sdbgmsg ="(line 16) domtf, mtfscal, smtffile"
  read(unit=lun, fmt=*, err=16) domtf, mtfscal, smtffile
  sdbgmsg ="(line 17) dovib, sig1 (, sig2, sig3)"
  read(unit=lun, fmt=*, err=16) dovib, sprm1, sprm2, sprm3
  dovib = min(2,max(0,dovib)) ! limit to valid range
  vibamp = 0.0
  vibamp2 = 0.0
  vibdir = 0.0
  select case (dovib)
    case (1)
      sdbgmsg ="isotropic vibration amplitude: vibamp"
      read(unit=sprm1, fmt=*, err=16) vibamp
      vibamp2 = vibamp
    case (2)
      sdbgmsg ="anisotropic vibration amplitude1: vibamp"
      read(unit=sprm1, fmt=*, err=16) vibamp
      sdbgmsg ="anisotropic vibration amplitude2: vibamp2"
      read(unit=sprm2, fmt=*, err=16) vibamp2
      sdbgmsg ="anisotropic vibration direction: vibdir"
      read(unit=sprm3, fmt=*, err=16) vibdir
  end select ! case (dovib)
  sdbgmsg ="(line 18) # aberrations below"
  read(unit=lun, fmt=*, err=16) na
  AF_wa = 0.0 ! zero all aberrations now
  if (na<=0) then
    na = 0
    nal = 0
  else
    nal = na
    do i=1, nal
      sdbgmsg ="ai, wax, way"
      read(unit=lun, fmt=*, err=16) ai, wax, way
      ai = ai + 1
      call AF_SetAberration(ai,wax,way)
    end do
  end if
  sdbgmsg ="oapr"
  read(unit=lun, fmt=*, err=16) oapr
  sdbgmsg ="oapx, oapy"
  read(unit=lun, fmt=*, err=16) oapx, oapy
  sdbgmsg ="nloop"
  read(unit=lun, fmt=*, err=16) nloop
  nloop = max(0, min( nloopmax, nloop ) )
  if (nloop>0) then
    do i=1, nloop
      write(unit=sdbgmsg,fmt='(A,I1,A)') "loop #",i,": class"
      read(unit=lun, fmt=*, err=16) lpcl(i)
      write(unit=sdbgmsg,fmt='(A,I1,A)') "loop #",i,": variable"
      read(unit=lun, fmt=*, err=16) lpvr(i)
      write(unit=sdbgmsg,fmt='(A,I1,A)') "loop #",i,": form"
      read(unit=lun, fmt=*, err=16) lpvf(i)
      write(unit=sdbgmsg,fmt='(A,I1,A)') "loop #",i,": v0, v1, #samp"
      read(unit=lun, fmt=*, err=16) lpv0(i), lpv1(i), lpsz(i)
      write(unit=sdbgmsg,fmt='(A,I1,A)') "loop #",i,": id string"
      read(unit=lun, fmt=*, err=16) lpid(i)
      ! post-read calculation depending on variation form
      if (lpvf(i)==1) lpvd(i) = ( lpv1(i) - lpv0(i) ) / real( lpsz(i)-1 )
      if (lpvf(i)==2) lpvd(i) = twopi / real(lpsz(i))
      ! the case lpvf(i)==3 will be handled later in subroutine checkprm
    end do
  end if
  
  close( unit=lun, iostat=nerr, err=15)
  
  
  
  call PostMessage("Finished loading parameters.")
  
  write(unit=sdbgmsg,fmt=*) "prm load: wave file: ", trim(swavfile)
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: wave dim.: ", nwx, nwy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: wave samp.: ", swx, swy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: HT: ", ht
  call PostDBGMessage(trim(sdbgmsg))
  select case (notype)
    case (0)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = TEM image"
    case (1)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = complex wave"
    case (2)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = wave amplitude"
    case (3)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = wave phase"
    case (4)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = wave real part"
    case (5)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = wave imag part"
    case (6)
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," = Thickness - Parameter map image"
    case default
      write(unit=sdbgmsg,fmt=*) "prm load: output type: ", notype," UNKNOWN! WARNING! Set back to 0 = TEM image"
      notype = 0
  end select ! case (notype)
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image file: ", trim(simgfile)
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image dim.: ", nx, ny
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image integer output: ", doint, iimean, el_conv, dark_noise
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image area selection: ", dofrm
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image sampl.: ", simg
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image frame offset: ", ox, oy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image frame rot.: ", ax
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image calculation model: ", ncohm
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image PTC opt.: ", doptc, fs
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image PSC opt.: ", dopsc, sc
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image MTF opt.: ", domtf, mtfscal, " ["//trim(smtffile)//"]"
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image vib opt.: ", dovib
  call PostDBGMessage(trim(sdbgmsg))
  select case (dovib)
    case (1)
      write(unit=sdbgmsg,fmt=*) "prm load: image vib amplitude [nm]: ", vibamp
      call PostDBGMessage(trim(sdbgmsg))
    case (2)
      write(unit=sdbgmsg,fmt=*) "prm load: image vib amplitudes [nm,nm,deg]: ", vibamp, vibamp2, vibdir
      call PostDBGMessage(trim(sdbgmsg))
  end select ! case (dovib)
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: image # aberrations: ", na
  call PostDBGMessage(trim(sdbgmsg))
  if (na>0) then
    do i=1, AF_maxaberration
      call AF_GetAberrationSName(i,1,smsg)
      write(unit=sdbgmsg,fmt='(A,A,A,2G12.4)') "prm load: ",trim(smsg)," ",AF_wa(1,i),AF_wa(2,i)
      call PostDBGMessage(trim(sdbgmsg))
    end do
  end if
  write(unit=sdbgmsg,fmt=*) "prm load: objective aperture semi opening angle: ", oapr
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: objective aperture decenter (x,y): ", oapx, oapy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: number of loops: ", nloop
  call PostDBGMessage(trim(sdbgmsg))
  if (nloop>0) then
    do i=1, nloop
      write(unit=sdbgmsg,fmt=*) "prm load: loop #",i,":  class:", lpcl(i)
      call PostDBGMessage(trim(sdbgmsg))
      write(unit=sdbgmsg,fmt=*) "prm load: loop #",i,":  variable id:", lpvr(i)
      call PostDBGMessage(trim(sdbgmsg))
      write(unit=sdbgmsg,fmt=*) "prm load: loop #",i,":  form:", lpvf(i)
      call PostDBGMessage(trim(sdbgmsg))
      write(unit=sdbgmsg,fmt=*) "prm load: loop #",i,":  v0, v1:", lpv0(i), lpv1(i)
      call PostDBGMessage(trim(sdbgmsg))
      write(unit=sdbgmsg,fmt=*) "prm load: loop #",i,":  # samples, step size:", lpsz(i), lpvd(i)
      call PostDBGMessage(trim(sdbgmsg))
      write(unit=sdbgmsg,fmt=*) "prm load: loop #",i,":  id string:", trim(lpid(i))
      call PostDBGMessage(trim(sdbgmsg))
    end do
  end if
  
  
  return
  
  ! handle errors
13 nerr = 1
  call CriticalError("Failed to acquire free logical file unit.")
14 nerr = 2
  call CriticalError("Failed to open parameter file.")
15 nerr = 3
  call CriticalError("Failed to close parameter file.")
16 nerr = 4
  call CriticalError("Failed to read data ("//trim(sdbgmsg)//") from parameter file.")

end subroutine loadprm



!**********************************************************************!
!
! subroutine readavlist
!
! reads an aberration parameter variation list from a file
!
! Expected list format in the file:
! line 01 : <list length> number defining the length of the list to read
! line 02 : <header length> number defining the length of the list header
! line 03 : <header> list of numbers defining the column content
! line 04+: <data set> list of floating point numbers defining the aberrations
!           Each aberration is defined by two values (ax, ay).
!           The length of <data set> is therefor twice the length of <header>
!
subroutine readavlist(sfile, nset, naval, lheader, lvalue, nierr)

  use wavimgprm
  implicit none
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(out) :: nset, naval, lheader(1:namax)
  real*4, intent(out) :: lvalue(1:namax,1:nalmax)
  integer*4, intent(inout) :: nierr
  
  integer*4 :: lun, i
  integer*4, external :: getfreelun
  character(len=600) :: smsg1, smsg2
  
  ! initialize
  nierr = 0
  nset = 0
  naval = 0
  lheader = 0
  lvalue = 0.0
  call PostMessage("Loading aberration variation list from file ["//trim(sfile)//"].")
  
  ! get a free logical file unit
  lun = getfreelun()
  if (lun<=0) goto 13
  
  ! open the file for reading
  open( unit=lun, file=trim(sfile), iostat=nerr, &
     &  action='read', status='old', err=14)
  
  ! read the length of the list (number of data sets)
  sdbgmsg ="list-length"
  read(unit=lun, fmt=*, err=16) nset
  if (nset>nalmax) then
    write(unit=smsg1, fmt=*) nset
    write(unit=smsg2, fmt=*) nalmax
    nset = nalmax
    call PostWarning("Invalid list length ("//trim(adjustl(smsg1))// &
     &               "), limited to "//trim(adjustl(smsg2))//".")
  end if
  
  ! read the length of the header (number of aberrations)
  sdbgmsg ="header-length"
  read(unit=lun, fmt=*, err=16) naval
  if (naval>namax) then
    write(unit=smsg1, fmt=*) naval
    write(unit=smsg2, fmt=*) namax
    naval = namax
    call PostWarning("Invalid header length ("//trim(adjustl(smsg1))// &
     &               "), limited to "//trim(adjustl(smsg2))//".")
  end if
  
  if (naval>0) then
    ! read the list header (list of aberration indices)
    sdbgmsg ="list-header"
    read(unit=lun, fmt=*, err=16) lheader(1:naval)
    
    if (nset>0) then
      ! read the data sets (lists of aberration values)
      do i=1, nset
        read(unit=lun, fmt=*, err=16) lvalue(1:2*naval,i)
      end do
    end if
  end if
  
  close( unit=lun, iostat=nerr, err=15)
  
  call PostMessage("Finished loading parameters.")
  
  return
  
  ! handle errors
13 nierr = 1
  call CriticalError("Failed to acquire free logical file unit.")
14 nierr = 2
  call CriticalError("Failed to open the file.")
15 nierr = 3
  call CriticalError("Failed to close the file.")
16 nierr = 4
  call CriticalError("Failed to read data ("//trim(sdbgmsg)//") from file.")
  
end subroutine readavlist     
!**********************************************************************!



!**********************************************************************!
!
! subroutine checkprm
!
! checks input parameters and transfers data to additional variables
!
subroutine checkprm

  use wavimgprm
  use AberrationFunctions

  implicit none
  
  logical :: fex
  real*4 :: dZ, gmax, dinfo
  real*4, external :: HT2WL
  integer*4 :: i, j, k
  character(len=5) :: snft
  
  write(unit=snft,fmt='(I5)') nft
  
  ! check ht
  if (ht<10.0 .or. ht>1300.0) call CriticalError("TEM high tension parameter out of range (10 - 1300 kV).")
  wl = HT2WL(ht)
  AF_lamb = wl
  write(unit=sdbgmsg,fmt=*) "parameter calc.: wavelength [nm]: ",wl
  call PostDBGMessage(trim(sdbgmsg))
  ! image parameters
  if (len_trim(simgfile_ex)>0) then
    write(unit=stmp,fmt='(A)') "Overriding image file name ("//trim(simgfile)//&
     &    ") with command-line parameter ("//trim(simgfile_ex)//")."
    call PostDBGMessage(trim(stmp))
    simgfile = simgfile_ex
  end if
  if (doint<0) doint = 0
  if (doint>2) doint = 2
  if (iimean<=0.0) then
    doint = 0
    call PostWarning("Invalid mean value for integer image, floating point image will be saved.")
  end if
  if (dark_noise<0.0) then
    dark_noise = 0.0
    call PostWarning("Invalid dark noise level for integer image, dark noise deactivated.")
  end if
  if (el_conv<=0.0) then
    el_conv = 1.0
    call PostWarning("Invalid counts per electron rate, reset to default 1.0.")
  end if
  ! temp. calculate info-limit from focus spread parameter
  dinfo = ( (pi*wl*fs)**2 / 8.0 )**0.25
  ginfo = 1/dinfo
  gmax = 0.5/min(swx,swy)
  write(unit=sdbgmsg,fmt=*) "parameter calc.: info limit [nm]: ",dinfo
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "parameter calc.: info limit [1/nm]: ",ginfo
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "parameter calc.: image Nyquist [1/nm]: ",gmax
  call PostDBGMessage(trim(sdbgmsg))
  
  ! calculate additional parameters for kernels
  nsre = 0
  if (sqrt(2.0)*gmax<3.0*ginfo) then
    call PostDBGMessage("  Choosing focal kernel sampling such that it is outsides the image spectrum.")
    gres1 = 1.1*sqrt(2.0)*gmax
  else
    call PostDBGMessage("  Choosing focal kernel sampling such that first linear resonance is equal to 3*ginfo.")
    gres1 = 3.0*ginfo
    nsre = 1
  end if
  dZ = 2.0/wl/(gres1**2)
  write(unit=sdbgmsg,fmt=*) "parameter calc.: focal kernel sampling [nm/pixel]: ",dZ
  call PostDBGMessage(trim(sdbgmsg))
  fkw = 2.0
  cbkw = 2.0
  NKFS = 1+int(fs*fkw/dZ)
  NKCB = 5
  write(unit=sdbgmsg,fmt=*) "parameter calc.: focal kernel size: ",NKFS
  call PostDBGMessage(trim(sdbgmsg))
  
  ! check file existence
  ! mtf
  if (domtf/=0) then
    domtf = 1
    inquire(file=trim(smtffile),exist=fex)
    if (.not.fex) call CriticalError("Detector MTF data file does not exist.")
    !if (mtfscal>1.0) call CriticalError("Detector MTF scaling parameter out of range (0.0 - 1.0).")
  end if
  
  ! check coherence parameters
  if (ncohm<1) ncohm = 1
  if (ncohm>6) ncohm = 6
  if (ncohm==5) then ! switch off resonance suppression
    gres1 = 0.0
    nsre = 0
  end if
  if (doptc/=0) then
    if (fs<=0.0) then
      fs = 0.0
      doptc = 0
      NKFS = 0
    end if
    if (fs>50.0) call PostWarning("Unusually large focus-spread parameter!")
    if (fkw<=0.0) then
      fs = 0.0
      doptc = 0
      NKFS = 0
      fkw = 0.0
    end if
    if (NKFS<=0) then
      fs = 0.0
      doptc = 0
      NKFS = 0
    end if
  else 
    fs = 0.0
    NKFS = 0
    fkw = 0.0
  end if
  if (dopsc/=0) then
    dopsc = 1
    if (sc<=0.0) then
      sc = 0.0
      dopsc = 0
      NKCB = 0
    end if
    if (sc>5.0) call PostWarning("Unusually large semi angle of convergence!")
    if (cbkw<=0.0) then
      sc = 0.0
      dopsc = 0
      NKCB = 0
      cbkw = 0.0
    end if
    if (NKCB<=0) then
      sc = 0.0
      dopsc = 0
      NKCB = 0
    end if
    sc = sc * 0.001 ! translate to radian
  else
    sc = 0.0
    NKCB = 0
    cbkw = 0.0
  end if
  if (dofrm/=0) then
    dofrm = 1
  else
    call PostDBGMessage("Image area selection turned OFF, using wave dimensions.")
    nx = nwx
    ny = nwy
    simg = swx
    ax = 0.0
  end if
  
  if (dovib<=0) then
    dovib = 0
    vibamp = 0.0
  else
    dovib = 1
    vibamp=abs(vibamp)
  end if
  
  if (nfoc_ex>0) then
    write(unit=stmp,fmt='(A,F8.3,A,F8.3,A)') &
     &    "Overriding defocus (",AF_wa(1,2),&
     &    " nm) with command-line parameter (",foc_ex," nm)."
    call PostDBGMessage(trim(stmp))
    AF_wa(1,2) = foc_ex
    AF_wa(2,2) = 0.0
  end if
  
  if (oapr_ex>0.0) then
    write(unit=stmp,fmt='(A,F8.3,A,F8.3,A)') &
     &    "Overriding objective aperture radius (",oapr,&
     &    " mrad) with command-line parameter (",oapr_ex," mrad)."
    call PostDBGMessage(trim(stmp))
    oapr = oapr_ex
  end if
  if (oapr<=0.0) then
    ! non-physical aperture value.
    call PostWarning("The objective aperture radius is not larger than 0. The aperture is removed.")
    ! - remove the aperture by setting the radius to min. 2*Nyquist(wave) for both directions
    ! - set the aperture de-center to zero
    oapr = max(1.0/swx,1.0/swy)
    oapx = 0.0
    oapy = 0.0
  end if
  
  ! check wave size
  if (nwx<32 .or. nwx>nft) call CriticalError("Wave function x-dim out of range (32 - "//trim(snft)//").")
  if (nwy<32 .or. nwy>nft) call CriticalError("Wave function y-dim out of range (32 - "//trim(snft)//").")
  if (swy<=0.0 .or. swy<=0.0) call CriticalError("Cannot use this wave function sampling rate.")
  
  ! check image size
  if (nx<32 .or. nx>nft) call CriticalError("image x-dim out of range (32 - "//trim(snft)//").")
  if (ny<32 .or. ny>nft) call CriticalError("image y-dim out of range (32 - "//trim(snft)//").")
  if (simg<=0.0) call CriticalError("Cannot use this image sampling rate.")
  
  ! check loop data
  if (nloop>0) then
    do i=nloop, 1, -1
      j = 0
      write(unit=stmp, fmt='(I1)') i
      if (lpsz(i)<=1.and.lpvf(i)<3) then
        ! the case lpsz(i)<=1.and.lpvf(i)==3 is allowed and should not cause problems
        call PostWarning("Loop #"//trim(stmp)//": loop size <= 1 makes no sense.")
        j = j + 1
      end if
      if (lpcl(i)<1.or.lpcl(i)>3) then
        call PostWarning("Loop #"//trim(stmp)//": variable class unknown, use 1, 2, or 3!")
        j = j + 1
      end if
      if (lpvf(i)<1.or.lpvf(i)>3) then
        call PostWarning("Loop #"//trim(stmp)//": variation form unknown, use 1, 2, or 3!")
        j = j + 1
      end if
      if (lpvf(i)==3) then
        if (lpcl(i)/=1) then
          call PostWarning("Loop #"//trim(stmp)//": unsupported listed parameter variation class, use 1!")
          j = j + 1
        end if
        inquire(file=trim(lpid(i)),exist=fex)
        if (.not.fex) then
          call PostWarning("Loop #"//trim(stmp)//": parameter variation list ["//trim(lpid(i))//"] not found!")
          j = j + 1
        end if
      end if
      if (j>0) then ! deactivate the loop
        if (i<nloop) then ! move down all the outer loops by one
          do k=i+1, nloop
            lpcl(i-1) = lpcl(i)
            lpvr(i-1) = lpvr(i)
            lpvf(i-1) = lpvf(i)
            lpv0(i-1) = lpv0(i)
            lpv1(i-1) = lpv1(i)
            lpsz(i-1) = lpsz(i)
            lpvd(i-1) = lpvd(i)
            lpid(i-1) = lpid(i)
          end do
        end if
        nloop = nloop - 1
        call PostWarning("Loop #"//trim(stmp)//" deactivated.")
      end if
    end do
  end if
  
  if (nloop>0) then ! determine and read list-directed parameter variation
    !
    ! The list input is done here because some loops may have been
    ! deactivated above and I am too lazy to manage a removal of lists. (JB)
    !
    ! clear the lists
    lpali = 0
    lpalv = 0.0
    do i=1, nloop
      if (lpvf(i)==3) then
        ! read the aberration variation list and fill the global array
        call readavlist(trim(lpid(i)), lpali(0,0), lpali(0,i), &
     &                  lpali(1:namax,i), lpalv(1:namax,1:nalmax,i),&
     &                  nerr)
        if (nerr/=0) then
          ! error while reading the list, try anyway
        end if
        lpsz(i) = lpali(0,0)
      end if
    end do
  end if
  
  if (notype==6) then ! sanity check for map images
    if (nloop<1) then
      ! no loop is defined, it makes no sense to calculate a map image
      ! therefore, the map image creation is switched off here and set to normal image creation
      call PostWarning("At least one loop definition is required to create a map image, switching to normal TEM image mode.")
      notype = 0
    end if
  end if
  
  return

end subroutine checkprm



!**********************************************************************!
!
! subroutine loadwave
!
! loads wave function data from file
!
! parameters:
!   character(len=*) :: sfile           = file name
!   integer*4 :: nx, ny                 = wave dimensions
!   complex*8 :: cdata(nx,ny)           = wave data (out)
!   integer*4 :: nerr                   = error code (out)
!
subroutine loadwave(sfile,nx,ny,cdata,nerr)

  implicit none
  
  integer*4, intent(in) :: nx, ny
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  complex*8, intent(inout) :: cdata(nx,ny)
  
  integer*4 :: lun
  integer*4, external :: getfreelun
  complex*8 :: c0
  
  nerr = 0
  c0 = cmplx(0.0,0.0)
  
  cdata = c0 ! aop
  
  lun = getfreelun()
  if (lun<=0) call CriticalError("Failed to acquire free logical file unit.")
  
  open(unit=lun,file=trim(sfile),iostat=nerr,&
     & form='binary',action='read',status='old')
  if (nerr/=0) call CriticalError("Failed to open the file ["//trim(sfile)//"].")
  
  read(unit=lun,iostat=nerr) cdata
  if (nerr/=0) call CriticalError("Failed while reading data from file.")
  
  close(unit=lun,iostat=nerr)
  if (nerr/=0) call CriticalError("Failed to close an open file ")
  
  return
  
end subroutine loadwave


!**********************************************************************!
!
! subroutine loadmtf
!
! loads the program mtf from file
!
! Format of the mtf file: TEXT ASCII
! line 01: <numpix>                     ! nuber of MTF pixels (1+CCD Nyquist)
! line 02: MTF(0)                       ! first MTF value (frequency i=0)
! line 03: MTF(1)                       ! 2nd MTF value (frequency i=1)
! ....
! line <numpix>: MTF(Nyquist)           ! last MTF value (frequency i=<numpix>-1)
!
subroutine loadmtf

  use wavimgprm
  implicit none
  
  integer*4 :: lun, i, i1, i2, nfit
  integer*4, external :: getfreelun
  real*4 :: ntfprm, ntfoff, x, y, mx, mx2, mxy, my
  
  nerr = 0
  nmtf = 0
  mtfdata = 0.0
  ntfdata = 0.0
  
  call PostMessage("Loading mtf from file ["//trim(smtffile)//"].")
  
  lun = getfreelun()
  if (lun<=0) goto 13
  
  open( unit=lun, file=trim(smtffile), iostat=nerr, &
     &  action='read', status='old', err=14)
  
  read(unit=lun, fmt=*, err=16) nmtf
  
  write(unit=sdbgmsg,fmt=*) "parameter load: mtf length: ", nmtf
  call PostDBGMessage(trim(sdbgmsg))
  
  if (nmtf<=0) goto 17
  nmtf = min(ntfmax, nmtf)
  
  do i=1, nmtf
    read(unit=lun, fmt=*, err=16) mtfdata(i)
  end do
  
  close( unit=lun, iostat=nerr, err=15)
  
  call PostMessage("Finished loading of mtf.")
  
  ! fill outer frequencies with last mtf value
  if (nmtf<ntfmax) then
    mtfdata(nmtf+1:ntfmax) = mtfdata(nmtf)
  end if
  
  call PostMessage("Extracting noise transfer function.")
  !  
  ! do a linear regression to the central quarter region
  !
  ! get the central quarter indices
  i1 = min(nmtf,max(1,nint( real(nmtf)*0.5 - real(nmtf)*0.125 )))
  i2 = min(nmtf,max(1,nint( real(nmtf)*0.5 + real(nmtf)*0.125 )))
  ! get the number of involved mtf data points
  nfit = i2-i1 + 1
  !
  ! get the log of the mtf data and prepare the means
  mx = 0.0
  my = 0.0
  mx2 = 0.0
  mxy = 0.0
  do i=i1, i2
    x = real(i-1)
    y = log(mtfdata(i))
    mx = mx + x
    my = my + y
    mx2 = mx2 + x*x
    mxy = mxy + x*y
  end do
  ! normalize
  mx = mx / real(nfit)
  my = my / real(nfit)
  mx2 = mx2 / real(nfit)
  mxy = mxy / real(nfit)
  !
  ! do the linear regression
  ! - slope parameter
  ntfprm = ( mxy - mx*my ) / ( mx2 - mx*mx )
  ! - offset parameter
  ntfoff = my - ntfprm * mx
  !
  write(unit=sdbgmsg,fmt=*) "parameter fit: ntf decay parameter: ", ntfprm
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "parameter fit: ntf offset parameter: ", exp( ntfoff )
  call PostDBGMessage(trim(sdbgmsg))
  !
  ! calculate the ntf
  do i=1, ntfmax
    ntfdata(i) = exp( ntfoff ) * ( exp( ntfprm * real(i-1) ) - 1.0 ) + 1.0
    ! ntfdata(i) = exp( ntfprm * real(i-1) )
  end do
  !
  
  return
  
  ! handle errors
13 nerr = 1
  call CriticalError("Failed to acquire free logical file unit.")
14 nerr = 2
  call CriticalError("Failed to open mtf file.")
15 nerr = 3
  call CriticalError("Failed to close mtf file.")
16 nerr = 4
  call CriticalError("Failed to read data from mtf file.")
17 nerr = 5
  call CriticalError("Invalid number of mtf pixels.")

end subroutine loadmtf


!**********************************************************************!
!
! FT
!
! does a fourier transform
!
! parameters:
!   integer*4 :: ndx, ndy               = real space data dim (in)
!   complex*8 :: crs(ndx, ndy)          = real space data (inout)
!   complex*8 :: cfs(ndy, ndx)          = fourier space data (inout)
!   character(len=3) :: dir             = transformation direction (in) ('for' or 'bac')
!
subroutine FT(ndx,ndy,crs,cfs,dir)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: ndx, ndy
  character(len=3), intent(in) :: dir
  complex*8, intent(inout) :: crs(ndx,ndy), cfs(ndy,ndx)
  
  integer*4 :: i, j
  
  if (ndx<32 .or. ndx>nft .or. ndy<32 .or. ndy>nft) call CriticalError("Invalid array size.")
  
  !cft = cmplx(0.0,0.0)
  
  if (dir=="for" .or. dir == "For" .or. dir=="FOR") then
    !call PostDBGMessage("Forward FFT in.")
    do j=1,ndy
      !do i=1,ndx
        cft(1:ndx,j) = crs(1:ndx,j)
      !end do
    end do
  end if
  
  if (dir=="bac" .or. dir == "Bac" .or. dir=="BAC") then
    !call PostDBGMessage("Backward FFT in.")
    do j=1,ndx
      !do i=1,ndy
        cft(1:ndy,j) = cfs(1:ndy,j)
      !end do
    end do
  end if
  
  call FSFFT(cft,ndx,ndy,dir)
  
  if (dir=="for" .or. dir == "For" .or. dir=="FOR") then
    !call PostDBGMessage("Forward FFT out.")
    do j=1,ndx
      !do i=1,ndy
        cfs(1:ndy,j) = cft(1:ndy,j)
      !end do
    end do
  end if
  
  if (dir=="bac" .or. dir == "Bac" .or. dir=="BAC") then
   !call PostDBGMessage("Backward FFT out.")
    do j=1,ndy
      do i=1,ndx
        crs(1:ndx,j) = cft(1:ndx,j)
      end do
    end do
  end if

end subroutine FT


!**********************************************************************!
!
! FSFFT
!
! does a fix size fourier transform
!
! parameters:
!   integer*4 :: ndx, ndy               = real space data dim (in)
!   complex*8 :: cdata(nft, nft)        = real space data (inout)
!   character(len=3) :: dir             = transformation direction (in) ('for' or 'bac')
!
subroutine FSFFT(cdata,ndx,ndy,dir)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: ndx, ndy
  character(len=*), intent(in) :: dir
  complex*8, intent(inout) :: cdata(nft,nft)
    
  character*40 :: direction
  external :: ODDCC128S, ODDCC256S, ODDCC512S, ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S ! link FFTs.f
  
  if (ndx<32 .or. ndx>nft .or. ndy<32 .or. ndy>nft) call CriticalError("Invalid array size.")
  direction = dir
  
  if (nft==128) call ODDCC128S(cdata,ndx,ndy,direction)
  if (nft==256) call ODDCC256S(cdata,ndx,ndy,direction)
  if (nft==512) call ODDCC512S(cdata,ndx,ndy,direction)
  if (nft==1024) call ODDCC1024S(cdata,ndx,ndy,direction)
  if (nft==2048) call ODDCC2048S(cdata,ndx,ndy,direction)
  if (nft==4096) call ODDCC4096S(cdata,ndx,ndy,direction)
  if (nft==8192) call ODDCC8192S(cdata,ndx,ndy,direction)
  
  return
  
end subroutine FSFFT

!**********************************************************************!
!
! subroutine AberrateWaveFourier
!
! function: applies aberrations to the given wave function
!
! parameter: complex*8 :: wave(ndimy,ndimx) (in fourier space)
!            integer*4 :: ndimx,ndimy : wave dimension
!            real*4 :: samplingx,samplingy : wave result sampling x & y
!
subroutine AberrateWaveFourier(wave,ndimx,ndimy,samplingx,samplingy)

  use AberrationFunctions
  use wavimgprm
  implicit none
  
  complex*8 :: wave(ndimy,ndimx)
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy
  
  integer*4 :: i, j, i1, j1, ndim2x, ndim2y, nd2x1, nd2y1
  real*4 :: wx, wy, wa, chi, sx, sy, tx, ty, apr, apx, apy
  real*4 :: ewx, ewy
  real*4 :: itogx, itogy, itowx, itowy
  real*4 :: aperture
  complex*8 :: cval, cphpl
  complex*8, allocatable :: phaseplate(:,:)
  external :: savedatac8

! ------------
! initialize
  if (AF_maxaberration<=0) call CriticalError("AF module not initialized.")
  
  ndim2x = int(ndimx/2)
  nd2x1  = ndim2x-1
  ndim2y = int(ndimy/2)
  nd2y1  = ndim2y-1
  sx = samplingx*real(ndimx)
  sy = samplingy*real(ndimy)
  itogx = 1.0/sx
  itowx = itogx*AF_lamb
  itogy = 1.0/sy
  itowy = itogy*AF_lamb
  tx = btx*0.001
  ty = bty*0.001
  apr = oapr*0.001
  apx = oapx*0.001
  apy = oapy*0.001
  
!  write(*,*) ndimx, ndim2x, samplingx, sx, itowx
!  write(*,*) ndimy, ndim2y, samplingy, sy, itowy
  
  allocate(phaseplate(ndimy,ndimx),stat=nerr)
  !write(*,*) AF_wa
  
! ------------
  do j=1,ndimx
    j1 = mod((j+nd2x1),ndimx)-ndim2x
    wx = real(j1)*itowx
    ewx = wx + tx
    
!    write(*,*) "j->wx: ",j,j1,wx
    
    do i=1,ndimy
      i1 = mod((i+nd2y1),ndimy)-ndim2y
      wy = real(i1)*itowy
      ewy = wy + ty
      
!      write(*,*) "i->wy: ",i,i1,wy
      
      chi = AF_AberrationFunction(ewx,ewy)
      
!      write(*,*) "chi=",chi
      
      cphpl = cmplx(0.0,-chi)
      
      cval = exp(cphpl)
      
      phaseplate(i,j) = cval
!      write(*,*) "phaseplate=",cval
      
      wa = sqrt((ewx-apx)**2 + (ewy-apy)**2)
      aperture = 0.5-0.5*tanh((wa-apr)/apr*100.0)
      
      wave(i,j) = cval*wave(i,j)*aperture
      
    end do
  end do
! ------------


! ------------
!  call PostDBGMessage("Saving current phase plate to [phplt.dat].")
!  call savedatac8("phplt.dat", ndimy*ndimx, phaseplate, nerr)
!  deallocate(phaseplate,stat=nerr)
  return

end subroutine AberrateWaveFourier
!**********************************************************************!



!**********************************************************************!
!
! subroutine DampenWaveFourier
!
! function: applies quasi-coherent (linear) dampening envelopes
!           and incoherent dampening envolopes (MTF and vibration)
!           to a wave function.
!           This should be the last modification of a wavefunction
!           before output. The dampening is applied to simulate the
!           resolution limiting effects in electron holography
!           and includes a sideband shift for the application of the
!           detector MTF.
!
! parameter: complex*8 :: wave(ndimy,ndimx) (in fourier space)
!            integer*4 :: ndimx,ndimy : wave dimension
!
subroutine DampenWaveFourier(wave,ndimx,ndimy)

  use AberrationFunctions
  use wavimgprm
  implicit none
  
  complex*8, intent(inout) :: wave(ndimy,ndimx)
  integer*4, intent(in) :: ndimx,ndimy
  
  integer*4 :: nalloc
  integer*4 :: i, j, i1, j1, k1, ndim2x, ndim2y, nd2x1, nd2y1
  real*4 :: six, siy, wx, wy, sx, sy, tx, ty
  real*4 :: pfkm, pfk, px, py, px2, py2, rx2, r2, pr, fk1
  real*4 :: ewx, ewy, dchiwx, dchiwy, dchitx, dchity
  real*4 :: v1, v12, v2, v22, cd, sd, vmean, vtmp
  real*4 :: itogx, itogy, itowx, itowy
  real*4 :: fspf, scpf, mtfv, Envt, Envs, rnmtf
  complex*8 :: Envv
  complex*8, allocatable :: cvib(:,:)

! ------------
! initialize
  if (AF_maxaberration<=0) call CriticalError("AF module not initialized.")
! - wave function nyquist numbers
  ndim2x = (ndimx-modulo(ndimx,2))/2 ! Nyquist X
  nd2x1  = ndim2x-1     ! ... minus 1
  ndim2y = (ndimy-modulo(ndimy,2))/2 ! Nyquist Y
  nd2y1  = ndim2y-1     ! ... minus 1
! - wave function output sampling rate
  if (dofrm/=0) then
    six = simg
    siy = simg
  else
    six = swx
    siy = swy
  end if
! - total wav function size
  sx = six*real(ndimx)  ! wave size X (nm)
  sy = siy*real(ndimy)  ! wave size Y (nm)
  itogx = 1.0/sx        ! Fourier-space sampling rate X (1/nm/pix)
  itowx = itogx*AF_lamb ! Fourier-space sampling rate X (mrad/pix)
  itogy = 1.0/sy        ! Fourier-space sampling rate Y (1/nm/pix)
  itowy = itogy*AF_lamb ! Fourier-space sampling rate Y (mrad/pix)
  tx = btx*0.001        ! optical transfer beam tilt X (rad)
  ty = bty*0.001        ! optical transfer beam tilt Y (rad)
  call AF_AberrationFunctionGradient(tx, ty, dchitx, dchity) ! zero beam offset
  v1 = vibamp           ! vibration amplitude (rms) (nm)
  v12 = v1*v1           ! ... square
  v2 = vibamp2          ! 2nd vibration amplitide (rms) (nm)
  v22 = v2*v2           ! ... square
  cd = COS(d2r*vibdir)  ! vibration main axis rotation cosine
  sd = SIN(d2r*vibdir)  ! vibration main axis rotation sine
  pfkm = real(nmtf-1)   ! number of MTF frequencies (Nyquist)
  mtfv = 1.0            ! preset mtf value
  Envv = cmplx(1.0, 0.0) ! preset vibration envelope value
  Envt = 1.0             ! preset temporal coherence envelope value
  Envs = 1.0             ! preset spatial coherence envelope value
  scpf = (sc/2)**2       ! spatial coherence parameter
  fspf = (pi*fs/wl/2)**2 ! temporal coherence parameter
  rnmtf = 1.0            ! MTF renormalization factor
  
  
  !
  ! preparing vibration envelope
  if (dovib/=0) then
    allocate(cvib(nft,nft), stat=nalloc)
    if (nalloc/=0) then
      call PostWarning("Allocation of memory for vibration envelope failed. Switching off vibration.")
      dovib=0
    else 
      ! setup the vibration kernel
      cvib = cmplx(0.0,0.0)
      ! The vibration kernel is set up in real-space and transform later to Fourier space.
      ! This is required since the kernel can be very narrow which requires a correct
      ! aliasing. Due to my laziness I decided to do this in real-space.
      !
      vmean = 0.0
      do j=1, ndimy
        ! calculate the y-coordinate
        if (j<ndim2y) then
          py = siy*real(j-1)
        else
          py = siy*real(j-1-ndimy)
        end if
        py2 = py * py
        do i=1, ndimx
        ! calculate the x-coordinate
          if (i<ndim2x) then
            px = six*real(i-1)
          else
            px = six*real(i-1-ndimx)
          end if
          px2 = px * px
          !
          ! calculate the convolution kernel
          ! using a gaussian of rms - half width
          vtmp = EXP( -((py2*v12 + px2*v22)*cd*cd +((px2*v12 + py2*v22)*sd &
                  & -2.*px*py*(v12-v22)*cd)*sd )/(2.*v12*v22) )
          cvib(i,j) = cmplx(vtmp ,0.0)
       end do
      end do
      ! tranform the kernel to Fourier space
      call FSFFT(cvib,ndimx,ndimy,'for')
      ! normalize
      vmean = cabs(cvib(1,1))
      cvib = cvib / vmean
      !
    end if
  end if
  !
  !
  ! side-band shifts are with respect to the input wave sampling
  if (sbshx/=0.and.ndbg/=0) then
    write(unit=sdbgmsg,fmt='(A,F)') "mtf: sideband shift x: ", sbshx
    call PostDBGMessage(trim(sdbgmsg))
  end if
  if (sbshy/=0.and.ndbg/=0) then
    write(unit=sdbgmsg,fmt='(A,F)') "mtf: sideband shift y: ", sbshy
    call PostDBGMessage(trim(sdbgmsg))
  end if
  if (dornsb/=0) then
    rx2 = ( sbshx / real(ndim2x) )**2
    r2 = rx2 + ( sbshy / real(ndim2y) )**2
    pr = sqrt(r2) ! relative (effective) frequency 
    pfk = pr * pfkm * mtfscal ! translate / re-scale to mtf pixel number
    ! access mtf data
    k1 = int(pfk)
    fk1 = pfk - real(k1)
    k1 = k1 + 1
    ! limit to defined mtf range
    if (k1<1) then
      k1 = 1
      fk1 = 0.0
    end if
    if (k1>=nmtf) then
      k1 = nmtf
      fk1 = 0.0
    end if
    ! interpolate mtf (linear)
    mtfv = mtfdata(k1)*(1.0-fk1) + mtfdata(k1+1)*fk1
    rnmtf = 1.0 / mtfv
    mtfv = 1.0
    if (ndbg/=0) then
      write(unit=sdbgmsg,fmt='(A,F,A)') "mtf: re-normalizing output wave (factor = ",rnmtf,")"
      call PostDBGMessage(trim(sdbgmsg))
    end if
  end if
  !
  do j=1,ndimx
    j1 = mod((j+nd2x1),ndimx)-ndim2x ! get X frequncy (pixel number)
    px = real(j1) ! ... same
    rx2 = ( (px + sbshx) / real(ndim2x) )**2 ! apply X side-band shift and rescale w.r.t Nyquist
    wx = px*itowx ! ... frequency as diffraction angle (rad)
    ewx = wx + tx ! ... plus beam tilt X
    
    do i=1,ndimy
      i1 = mod((i+nd2y1),ndimy)-ndim2y ! get y frequency (pixel number)
      py = real(i1) ! ... same
      r2 = rx2 + ( (py + sbshy) / real(ndim2y) )**2 ! apply Y side-band shift and rescale to Nyquist, and add up to rel. frequency^2
      wy = py*itowy ! ... frequency as diffraction angle (rad)
      ewy = wy + ty ! ... plus beam tilt Y
      
      if (dopsc/=0) then
        ! aberration function gradient
        call AF_AberrationFunctionGradient(ewx, ewy, dchiwx, dchiwy)
        ! spatial envelope value
        Envs = exp( -scpf*( (dchiwx-dchitx)**2 + (dchiwy-dchity)**2 ) )
      end if
      if (doptc/=0) then
        ! temporal envelope value
        Envt = exp( -fspf*( ( ewx**2 - tx**2 + ewy**2 - ty**2 )**2) )
      end if
      
      ! get the mtf value
      if (domtf/=0) then
        ! determine scale to mtf 
        pr = sqrt(r2) ! relative (effective) frequency 
        pfk = pr * pfkm * mtfscal ! translate / re-scale to mtf pixel number
        ! access mtf data
        k1 = int(pfk)
        fk1 = pfk - real(k1)
        k1 = k1 + 1
        ! limit to defined mtf range
        if (k1<1) then
          k1 = 1
          fk1 = 0.0
        end if
        if (k1>=nmtf) then
          k1 = nmtf
          fk1 = 0.0
        end if
        ! interpolate mtf (linear)
        mtfv = rnmtf * (mtfdata(k1)*(1.0-fk1) + mtfdata(k1+1)*fk1)
      end if
      
      if (dovib/=0) then
        Envv = conjg( cvib(i,j) )
      end if
      
      !
      ! apply all envelopes
      !
      wave(i,j) = wave(i,j) * mtfv * Envs * Envt * Envv
      
    end do
  end do
  
  if (allocated(cvib)) deallocate(cvib, stat=nalloc)

  return

end subroutine DampenWaveFourier
!**********************************************************************!




!**********************************************************************!
!
! subroutine ApplyPTCExpl
!
! function: applies partial temporal coherence
!           by explicit focal averaging
!
subroutine ApplyPTCExpl(nx,ny,swx,swy,fs,wl,tx,ty,waveft,wave)

  implicit none

  real*4, parameter :: pi = 3.1415927
  real*4, parameter :: dt = 9.0/16.0
  integer*4, parameter :: NFK = 21
  real*4, parameter :: FSZ = 3.0

  complex*8 :: waveft(ny,nx), wave(nx,ny)
  integer*4 :: nx,ny
  real*4 :: swx,swy,fs,wl,tx,ty
  
  integer*4 :: i, j, i1, j1, ndim2x, ndim2y, nd2x1, nd2y1, k
  integer*4 :: nalloc
  real*4 :: wx, wy, sx, sy, wx2, w2, foc, focfac, fks, fwt, dmp
  real*4 :: itogx, itogy, itowx, itowy, ph
  complex*8 :: cval, c0
  character(len=400) :: smsg
  real*4, allocatable :: phgr(:,:)
  complex*8, allocatable :: tmpft(:,:), tmprs(:,:)
  external :: FSFFT

! ------------
! initialize
  nalloc = 0
  c0 = cmplx(0.0,0.0)
  tmpft = c0
  wave = c0
  phgr = 0.0
  ndim2x = int(nx/2)
  nd2x1  = ndim2x-1
  ndim2y = int(ny/2)
  nd2y1  = ndim2y-1
  sx = swx*real(nx)
  sy = swy*real(ny)
  itogx = 1.0/sx
  itowx = itogx*wl
  itogy = 1.0/sy
  itowy = itogy*wl
  focfac = pi/wl

! ------------
! debug noise  
  call PostDBGMessage("Applying ptc explicitly.")
  write(unit=smsg,fmt=*) " defocus spread [nm]:",fs
  call PostDBGMessage(trim(smsg))
  write(unit=smsg,fmt=*) " input img size (RS):",nx,ny
  call PostDBGMessage(trim(smsg))
  write(unit=smsg,fmt=*) " input sampling (RS):",swx,swy
  call PostDBGMessage(trim(smsg))
  write(unit=smsg,fmt=*) " input sampling (FS):",itowx,itowy
  call PostDBGMessage(trim(smsg))
  write(unit=smsg,fmt=*) " input beam tilt [rad]:",tx,ty
  call PostDBGMessage(trim(smsg))
  write(unit=smsg,fmt=*) " size if focal kernel:",NFK
  call PostDBGMessage(trim(smsg))
  write(unit=smsg,fmt=*) " focal extent of kernel:",FSZ
  call PostDBGMessage(trim(smsg))
  
! ------------
! allocations
  allocate(tmpft(ny,nx), tmprs(nx,ny), stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("Memory allocation failed on init of (ApplyPTCExpl).")  
  end if
  tmpft = cmplx(0.0,0.0)
  tmprs = cmplx(0.0,0.0)
  allocate(phgr(ny,nx), stat=nalloc)
  if (nalloc/=0) then
    call CriticalError("Memory allocation failed on init of (ApplyPTCExpl).")
  end if
  phgr = 0.0
  
! ------------
! prepare focus phase shifts
  do j=1,nx
    j1 = mod((j+nd2x1),nx)-ndim2x
    wx = real(j1)*itowx+tx
    wx2 = wx*wx
    do i=1,ny
      i1 = mod((i+nd2y1),ny)-ndim2y
      wy = real(i1)*itowy+ty
      w2 = wy*wy+wx2
      phgr(i,j) = -focfac*w2
    end do
  end do
  
! ------------
  fks = 0.0
  do k=1,NFK
    foc = (2.0*real(k-1)/real(NFK-1)-1.0)*FSZ*fs
    fwt = exp(-foc*foc/fs/fs)
    fks = fks + fwt
    write(unit=smsg,fmt='(A,I2,A,F7.2,A,F7.4)') "  Kernel Step: ",k,"  defocus: ",foc," nm  weight: ",fwt
    call PostDBGMessage(trim(smsg))
    do j=1,nx
      do i=1,ny
        ph = foc*phgr(i,j)
        cval = cmplx(cos(ph),sin(ph))
        tmpft(i,j) = cval*waveft(i,j)
      end do
    end do
    call FT(nx,ny,tmprs,tmpft,"bac")
    do j=1,ny
      do i=1,nx
        cval = tmprs(i,j)
        wave(i,j) = wave(i,j) + cval*conjg(cval)*fwt
      end do
    end do
  end do
  wave = wave / fks
  !call savedatac8("wavexp.dat", nx* ny, wave, i)
  
  ! save apply exponental envelope at outer parts of the FT to
  ! suppress resonances of the focal kernel
  call PostMessage("  Suppressing explicit focal averaging resonances.")
  call FT(nx,ny,wave,waveft,"for")
  fwt = (0.5*pi*fs/wl)**2.0 ! damping constant
  ph = dt*wl*real(NFK-1)/fs/FSZ
  !write(*,*) sqrt(ph)
  do j=1,nx
    j1 = mod((j+nd2x1),nx)-ndim2x
    wx = real(j1)*itowx+tx
    wx2 = wx*wx
    do i=1,ny
      i1 = mod((i+nd2y1),ny)-ndim2y
      wy = real(i1)*itowy+ty
      w2 = wy*wy+wx2
      if (w2>ph) then
        dmp = exp(-fwt*w2*w2)
        waveft(i,j) = waveft(i,j) * dmp
      end if
      
    end do
  end do
  call FT(nx,ny,wave,waveft,"bac")
! ------------

  deallocate(phgr, tmpft, tmprs, stat=nalloc)
  return

END SUBROUTINE ApplyPTCExpl
!**********************************************************************!

!**********************************************************************!
!
! subroutine ApplyPTCEnv
!
! function: applies partial temporal coherence envelope to FS-image data
!
! parameter: complex*8 :: cdat(ndimy,ndimx) (in fourier space)
!            integer*4 :: ndimx,ndimy : wave dimension
!            real*4 :: samplingx,samplingy : RS sampling x & y
!            real*4 :: fs : focus-spread parameter in nm
!            real*4 :: wl : el. wavelength in nm
!            real*4 :: tx, ty : beam tilt in mrad
!
subroutine ApplyPTCEnv(ndimx,ndimy,samplingx,samplingy,cdat,fs,wl,tx,ty)

  implicit none

  real*4, parameter :: pi = 3.1415927

  complex*8 :: cdat(ndimy,ndimx)
  integer*4 :: ndimx,ndimy
  real*4 :: samplingx,samplingy,fs,wl,tx,ty
  
  integer*4 :: i, j, i1, j1, ndim2x, ndim2y, nd2x1, nd2y1
  real*4 :: wx, wy, sx, sy, rdmp, dmpfac,wx2,w2
  real*4 :: itogx, itogy, itowx, itowy
  complex*8 :: cval
  character(len=1024) :: sdbgmsg

! ------------
! initialize
  ndim2x = int(ndimx/2)
  nd2x1  = ndim2x-1
  ndim2y = int(ndimy/2)
  nd2y1  = ndim2y-1
  sx = samplingx*real(ndimx)
  sy = samplingy*real(ndimy)
  itogx = 1.0/sx
  itowx = itogx*wl
  itogy = 1.0/sy
  itowy = itogy*wl
  dmpfac = - (0.5*pi*fs/wl)**2.0
  
  
! ------------
  call PostDBGMessage("Applying ptc envelope.")
  write(unit=sdbgmsg,fmt=*) " defocus spread [nm]:",fs
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input img size (RS):",ndimx,ndimy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input sampling (RS):",samplingx,samplingy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input sampling (FS):",itowx,itowy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input beam tilt [rad]:",tx, ty
  call PostDBGMessage(trim(sdbgmsg))
  do j=1,ndimx
    j1 = mod((j+nd2x1),ndimx)-ndim2x
    wx = real(j1)*itowx+tx
    wx2 = wx*wx
    do i=1,ndimy
      i1 = mod((i+nd2y1),ndimy)-ndim2y
      wy = real(i1)*itowy+ty
      w2 = wy*wy+wx2
      rdmp = exp( dmpfac*w2*w2)
      cval = cdat(i,j)*rdmp
      cdat(i,j) = cval
    end do
  end do

! ------------
  return

end subroutine ApplyPTCEnv
!**********************************************************************!



!**********************************************************************!
!
! subroutine ApplyPSCEnv
!
! function: applies partial spatial coherence envelope to image intensity
!
! parameter: (global common parameters in wavimgprm.f90)
!            complex*8 :: cdat(ndimy,ndimx) (in real space) (in/out)
!            integer*4 :: nix, niy : wave dimension (in)
!
! requires use of module AberrationFunctions, must be initialzed outsides
!
subroutine ApplyPSCEnv(nix,niy,cdat)

  use AberrationFunctions
  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy
  complex*8, intent(inout) :: cdat(nix,niy)
    
  integer*4 :: i, j, i1, j1, n2x, n2y, n2x1, n2y1
  real*4 :: wx, wy, rdmp, dmpfac
  real*4 :: dcx, dcy
  real*4 :: itogx, itogy, itowx, itowy, tx, ty
!  real*4, allocatable :: pscdmp(:,:)
  complex*8 :: cval, cval0
  
  external :: FSFFT ! FFT, link FFTs.f
  external :: savedata

! ------------
! initialize
  nerr = 0
  cval0 = cmplx(0.0,0.0)
  n2x = int(nix/2)
  n2x1  = n2x-1
  n2y = int(niy/2)
  n2y1  = n2y-1
  itogx = 1.0/swx/real(nix)
  itowx = itogx*wl
  itogy = 1.0/swy/real(niy)
  itowy = itogy*wl
  dmpfac = - (0.5*sc)**2.0
  cft = cval0
  tx = btx*0.001
  ty = bty*0.001
!  allocate(pscdmp(niy,nix),stat=nerr)
  
  call PostDBGMessage("Applying psc envelope.")
  write(unit=sdbgmsg,fmt=*) " convergence semi angle [rad]:",sc
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input img size (RS):",nix,niy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input sampling (RS):",swx,swy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) " input sampling (FS):",itowx,itowy
  call PostDBGMessage(trim(sdbgmsg))

! ------------
! copy data
  do j=1, niy
    do i=1, nix
      cft(i,j) = cdat(i,j)
    end do
  end do
  call PostDBGMessage("input data copied.")

! ------------
! FFT, image fourier transform
  call PostDBGMessage("Forward FFT.")
  call FSFFT(cft,nix,niy,'for')
  

! ------------
! apply envelope to image FT
  do j=1,nix
    j1 = mod((j+n2x1),nix)-n2x
    wx = real(j1)*itowx
    do i=1,niy
      i1 = mod((i+n2y1),niy)-n2y
      wy = real(i1)*itowy
      call AF_AberrationFunctionGradient(wx+tx, wy+ty, dcx, dcy)
      rdmp = exp( dmpfac*(dcx*dcx+dcy*dcy))
      !write(*,*) i,j,wx,wy,dcx,dcy,dmpfac,rdmp
!      pscdmp(i,j) = rdmp
      cval = cft(i,j)*rdmp
      cft(i,j) = cval
    end do
  end do
  
!  call savedata("pscdmp.dat", niy* nix, pscdmp, nerr)
!  deallocate(pscdmp,stat=nerr)
  
! ------------
! FFT, image fourier transform
  call PostDBGMessage("Backward FFT.")
  call FSFFT(cft,nix,niy,'bac')
  
! ------------
! copy data
  do j=1, niy
    do i=1, nix
      cdat(i,j) = cft(i,j)
    end do
  end do
  call PostDBGMessage("result data copied.")

  return

end subroutine ApplyPSCEnv


!**********************************************************************!
!
! subroutine GetCohTF
!
! function: Calculates the value of the coherent (wave) transfer function
!           for a given g-vector and extra defocus dZ 
!
! parameter: (global common parameters in wavimgprm.f90)
!            complex*8 :: cohtf,  (in/out)
!            real*4 :: wx, wy, dZ : local parameters (in)
!
! requires use of module AberrationFunctions, must be initialzed outsides
!
subroutine GetCohTF(wx,wy,dZ,cohtf)

  use AberrationFunctions
  use wavimgprm
  implicit none
  
  real*4, intent(in) :: wx, wy, dZ
  complex*8, intent(inout) :: cohtf
  real*4 :: tx, ty ! beam tilts
  real*4 :: chi ! current phase shift due to local aberrations
  tx = btx*0.001
  ty = bty*0.001
  chi = AF_AberrationFunction(wx+tx,wy+ty) + pi*((wx+tx)**2+(wy+ty)**2)*dZ/wl
  cohtf = exp(cmplx(0.0,-chi))
  
  return

end subroutine GetCohTF

!**********************************************************************!
!
! function bilinr4
!
! returns bilinear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*4 function bilinr4(a,n,m,x,y)
  implicit none
  real*4, intent(in) :: a(n,m), x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rval                    ! temp values
  rval = 0.0
  ! get quasi pixel numbers left and below the target point
  i1 = FLOOR(x)
  j1 = FLOOR(y)
  ! get quasi pixel numbers right and above the target point
  i2 = i1+1
  j2 = j1+1
  ! get the distances from target point to the lower left data point
  f1 = x-real(j1)
  f2 = y-real(j2)
  ! wrap the quasi pixel numbers to the bounds of the data array
  i1 = 1+modulo(i1-1,n)
  j1 = 1+modulo(j1-1,m)
  i2 = 1+modulo(i2-1,n)
  j2 = 1+modulo(j2-1,m)
  ! sum pu the values
  rval = rval + (1.0-f1)*(1.0-f2)*a(i1,j1)
  rval = rval + f1      *(1.0-f2)*a(i2,j1)
  rval = rval + (1.0-f1)*f2      *a(i1,j2)
  rval = rval + f1      *f2      *a(i2,j2)
  bilinr4 = rval
end function bilinr4


!**********************************************************************!
!
! real*4 function cubickernelr4
!
! defines cubic interpolation kernel
!
! INPUT:
!   real*4 :: t             = offset of supporting point from
!                             interpolation position
!   real*4 :: a             = interpolation smoothness parameter
!
real*4 function cubickernelr4(t,a)
  implicit none
  real*4, intent(in) :: t, a
  real*4 :: tabs
  tabs = abs(t)
  cubickernelr4 = 0.0
  if (tabs<=1.0) then
    cubickernelr4 = 1.0+tabs*tabs*((a+2.0)*tabs-(a+3.0))
  else if (tabs<2.0 .and. tabs>1.0) then
    cubickernelr4 = a*(tabs*(8.0-tabs*(5.0-tabs))-4.0)
  end if
  return
end function cubickernelr4

!**********************************************************************!
!
! function bicubr4
!
! returns bicubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*4 function bicubr4(a,n,m,x,y)
  implicit none
  real*4, intent(in) :: a(n,m), x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rval, bcx, bcy          ! temp values
  real*4 :: bip
  real*4, external :: cubickernelr4
  bip = -0.5
  rval = 0.0
  ! get quasi pixel numbers i1 and j1 left and below of the fractional coordinate
  ! the pixel numbers may be out of bounds of the data array, we wrap them in later
  i1 = FLOOR(x)
  j1 = FLOOR(y)
  do j=j1-1, j1+2 ! kernel loop
    ! get the y-distance between target point y and kernel point j
    f2 = y-real(j)
    ! get the kernel weight for the y-direction
    bcy = cubickernelr4(f2,bip)
    ! wrap the pixel number to the data array y-bounds
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      ! get the x-distance between the target point x and kernel point i
      f1 = x-real(i)
      ! get the kernel weight for the x-direction
      bcx = cubickernelr4(f1,bip)
      ! wrap the pixel number to the data array x-bounds
      i2 = 1+modulo(i-1,n)
      ! sum up the kernel-weighted data
      rval = rval + a(i2,j2)*bcx*bcy
    end do
  end do
  bicubr4 = rval ! transfer to interface
end function bicubr4


!**********************************************************************!
!
! function bicubc8
!
! returns bicubic interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n,m)     = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
complex*8 function bicubc8(a,n,m,x,y)
  implicit none
  complex*8, intent(in) :: a(n,m)
  real*4, intent(in) :: x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: bcx, bcy                ! temp values
  real*4 :: bip
  complex*8 :: cval
  real*4, external :: cubickernelr4
  bip = -0.5
  cval = cmplx(0.0,0.0)
  ! get quasi pixel numbers i1 and j1 left and below of the fractional coordinate
  ! the pixel numbers may be out of bounds of the data array, we wrap them in later
  i1 = FLOOR(x)
  j1 = FLOOR(y)
  do j=j1-1, j1+2 ! kernel loop
    ! get the y-distance between target point y and kernel point j
    f2 = y-real(j)
    ! get the kernel weight for the y-direction
    bcy = cubickernelr4(f2,bip)
    ! wrap the pixel number to the data array y-bounds
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      ! get the x-distance between the target point x and kernel point i
      f1 = x-real(i)
      ! get the kernel weight for the x-direction
      bcx = cubickernelr4(f1,bip)
      ! wrap the pixel number to the data array x-bounds
      i2 = 1+modulo(i-1,n)
      ! sum up the kernel-weighted data
      cval = cval + a(i2,j2)*bcx*bcy
    end do
  end do
  bicubc8 = cval ! transfer to interface
end function bicubc8


!**********************************************************************!
!
! subroutine getimageframe
!
! extracts real*4 image data from given frame of wave field (complex*8)
! data, however the complex*8 data is already intensity data
! just do the typecast here
!
! parameters:
!   integer*4 :: nix, niy               = input data dim (in)
!   complex*8 :: cdat(nix, niy)         = input data (in)
!   integer*4 :: nox, noy               = output data dim (in)
!   complex*8 :: rdat(nox, noy)         = output data (out)
!
subroutine getimageframe(nix,niy,cdat,nox,noy,rdat)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy, nox, noy
  complex*8, intent(in) :: cdat(nix,niy)
  real*4, intent(inout) :: rdat(nox,noy)
  
  integer*4 :: i, j, i1, j1, nalloc
  real*4 :: di, dj, dwx, dwy, cs, ss, px, py, rval
  real*4, allocatable :: rin(:,:)
  real*4, external :: bilinr4, bicubr4
  
  write(unit=sdbgmsg,fmt=*) "wave data size:", nix, niy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "image data size:", nox, noy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "image area sel. mode:", dofrm
  call PostDBGMessage(trim(sdbgmsg))
  
  if (dofrm/=0) then
  
    write(unit=sdbgmsg,fmt=*) "wave offset:", ox, oy
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "wave rotation:", ax
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "wave sampling:", swx, swy
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "image sampling:", simg
    call PostDBGMessage(trim(sdbgmsg))
    
    ! allocate array
    nalloc = 0
    allocate(rin(nix,niy), stat=nalloc)
    rin = 0.0
  
    ! type-cast data to real*4
    do j=1, niy
      do i=1, nix
        rin(i,j) = abs(real(cdat(i,j)))
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "image data copied from wave sampling."
    call PostDBGMessage(trim(sdbgmsg))
    
    if (ndbg/=0) then ! speed issue (full array operation)
      rval = sum(rin)/real(nix*niy)
      write(unit=stmp,fmt='(G12.4)') rval
      call PostDBGMessage("wave size image mean value: "//trim(adjustl(stmp)))
    end if
  
    cs = cos(d2r*ax)
    ss = sin(d2r*ax)
  
    ! now to the frame op
    do j=1, noy
      dj = simg*real(j-1) ! image y distance from offset
      do i=1, nox
        di = simg*real(i-1) ! image x distance from offset
        dwx = di*cs + dj*ss ! wave x-coordinate
        dwy = -di*ss + dj*cs ! wave y-coordinate
        px = 1.0 + ox + dwx/swx ! translate from physical wave coords to pixels
        py = 1.0 + oy + dwy/swy ! translate from physical wave coords to pixels
        rdat(i,j) = bicubr4(rin,nix,niy,px,py) ! bicubic interpolation
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "image data interpolated."
    call PostDBGMessage(trim(sdbgmsg))
    
    if (ndbg/=0) then ! speed issue (full array operation)
      rval = sum(rdat)/real(nox*noy)
      write(unit=stmp,fmt='(G12.4)') rval
      call PostDBGMessage("interpolated image mean value: "//trim(adjustl(stmp)))
    end if
    
    deallocate(rin, stat=nalloc)
  
  else
  
    do j=1, noy
      j1 = 1 + modulo(j-1,niy)
      do i=1, nox
        i1 = 1 + modulo(i-1,nix)
        rdat(i,j) = abs(real(cdat(i1,j1)))
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "image data copied directly."
    call PostDBGMessage(trim(sdbgmsg))
    
    if (ndbg/=0) then ! speed issue (full array operation)
      rval = sum(rdat)/real(nox*noy)
      write(unit=stmp,fmt='(G12.4)') rval
      call PostDBGMessage("copied image mean value: "//trim(adjustl(stmp)))
    end if
  
  end if

end subroutine getimageframe



!**********************************************************************!
!
! subroutine getimageframe2
!
! applies area selection to 2d-data, real*4
!
! parameters:
!   integer*4 :: nix, niy               = input data dim (in)
!   real*4 :: rin(nix, niy)             = input data (in)
!   integer*4 :: nox, noy               = output data dim (in)
!   real*4 :: rout(nox, noy)            = output data (out)
!
subroutine getimageframe2(nix,niy,rin,nox,noy,rout)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy, nox, noy
  real*4, intent(in) :: rin(nix,niy)
  real*4, intent(inout) :: rout(nox,noy)
  
  integer*4 :: i, j
  real*4 :: di, dj, dwx, dwy, cs, ss, px, py
  real*4, external :: bilinr4, bicubr4
  
  write(unit=sdbgmsg,fmt=*) "input data size:", nix, niy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "output data size:", nox, noy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "image area sel. mode:", dofrm
  call PostDBGMessage(trim(sdbgmsg))
  
  if (dofrm/=0) then
  
    write(unit=sdbgmsg,fmt=*) "input offset:", ox, oy
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "input rotation:", ax
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "input sampling:", swx, swy
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "output sampling:", simg
    call PostDBGMessage(trim(sdbgmsg))
    
  
    cs = cos(d2r*ax)
    ss = sin(d2r*ax)
  
    ! now to the frame op
    do j=1, noy
      dj = simg*real(j-1) ! image y distance from offset
      do i=1, nox
        di = simg*real(i-1) ! image x distance from offset
        dwx = di*cs + dj*ss ! wave x-coordinate
        dwy = -di*ss + dj*cs ! wave y-coordinate
        px = 1.0 + ox + dwx/swx ! translate from physical wave coords to pixels
        py = 1.0 + oy + dwy/swy ! translate from physical wave coords to pixels
        rout(i,j) = bicubr4(rin,nix,niy,px,py) ! bilinear interpolation
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "data interpolated."
    call PostDBGMessage(trim(sdbgmsg))
  
  else
  
    do j=1, noy
      do i=1, nox
        rout(i,j) = rin(i,j)
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "data copied directly."
    call PostDBGMessage(trim(sdbgmsg))
  
  end if

end subroutine getimageframe2


!**********************************************************************!
!
! subroutine getimageframe3
!
! applies area selection to 2d-data, complex*8
!
! parameters:
!   integer*4 :: nix, niy               = input data dim (in)
!   complex*8 :: cin(nix, niy)          = input data (in)
!   integer*4 :: nox, noy               = output data dim (in)
!   complex*8 :: cout(nox, noy)         = output data (out)
!
subroutine getimageframe3(nix,niy,cin,nox,noy,cout)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy, nox, noy
  complex*8, intent(in) :: cin(nix,niy)
  complex*8, intent(inout) :: cout(nox,noy)
  
  integer*4 :: i, j
  real*4 :: di, dj, dwx, dwy, cs, ss, px, py
  complex*8, external :: bicubc8
  
  call PostDBGMessage("(getimageframe3)")
  write(unit=sdbgmsg,fmt=*) "input data size:", nix, niy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "output data size:", nox, noy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "image area sel. mode:", dofrm
  call PostDBGMessage(trim(sdbgmsg))
  
  if (dofrm/=0) then
  
    write(unit=sdbgmsg,fmt=*) "input offset:", ox, oy
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "input rotation:", ax
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "input sampling:", swx, swy
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt=*) "output sampling:", simg
    call PostDBGMessage(trim(sdbgmsg))
    
  
    cs = cos(d2r*ax)
    ss = sin(d2r*ax)
  
    ! now to the frame op
    do j=1, noy
      dj = simg*real(j-1) ! image y distance from offset
      do i=1, nox
        di = simg*real(i-1) ! image x distance from offset
        dwx = di*cs + dj*ss ! wave x-coordinate
        dwy = -di*ss + dj*cs ! wave y-coordinate
        px = 1.0 + ox + dwx/swx ! translate from physical wave coords to pixels
        py = 1.0 + oy + dwy/swy ! translate from physical wave coords to pixels
        cout(i,j) = bicubc8(cin,nix,niy,px,py) ! bilinear interpolation
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "data interpolated."
    call PostDBGMessage(trim(sdbgmsg))
  
  else
  
    do j=1, noy
      do i=1, nox
        cout(i,j) = cin(i,j)
      end do
    end do
    
    write(unit=sdbgmsg,fmt=*) "data copied directly."
    call PostDBGMessage(trim(sdbgmsg))
  
  end if

end subroutine getimageframe3

!! removed wrong code
!
!!**********************************************************************!
!!
!! subroutine ApplyMTF
!!
!! Applies mtf to given image
!! (only if domtf flag is set in global parameters)
!!
!! INPUT: (global parameters in wavimgprm.f90)
!!   integer*4 :: nix, niy           = image dimensions
!!
!! IN/OUTPUT:  (global parameters in wavmgprm.f90)
!!   real*4 :: img(nix, niy)         = image data
!!
!subroutine ApplyMTF(img,nix,niy)
!
!  use wavimgprm
!  implicit none
!  
!  integer*4, intent(in) :: nix, niy
!  real*4, intent(inout) :: img(nix,niy)
!  
!  integer*4 :: i, j, i1, j1         ! iterators
!  integer*4 :: nnx, nny, nnx1, nny1 ! nyquist numbers
!  integer*4 :: k1                   ! mtf pos
!  real*4 :: px, py                  ! positions
!  real*4 :: px2, pr2                ! position squared radius
!  real*4 :: pf, pfk, pfim, pfkm     ! frequency interpolation pos
!  real*4 :: fk1                     ! mtf pos fraction
!  real*4 :: vmtf                    ! mtf value to apply
!  
!  
!  external :: FSFFT            ! fourier transform, link FFTs.f
!  
!  ! initialize
!  nerr = 0
!  !   check mtf application request flag
!  if (domtf==0) return
!  cft = cmplx(0.0,0.0)
!  nnx = nix / 2
!  nnx1 = nnx - 1
!  nny = niy / 2
!  nny1 = nny - 1
!  pfim = real(max(nnx,nny))
!  pfkm = real(nmtf-1)
!  
!  ! call PostDBGMessage("Applying MTF forward FFT.")
!  ! copy data
!  do j=1, niy
!    do i=1, nix
!      cft(i,j) = cmplx(img(i,j),0.0)
!    end do
!  end do
!  
!  ! tranform to fourier space
!  call FSFFT(cft,nix,niy,'for')
!  
!  ! apply mtf to ft image data
!  do j=1, nix
!    ! get x-frequency
!    j1 = mod((j+nnx1),nix)-nnx
!    px = real(j1)
!    px2 = px*px
!    
!    do i=1, niy
!      ! get y-frequency
!      i1 = mod((i+nny1),niy)-nny
!      py = real(i1)
!      pr2 = px2 + py*py
!      
!      ! 2d-frequency modulus
!      pf = sqrt(pr2)
!      
!      ! determine scale to mtf 
!      pfk = pf / pfim * pfkm * mtfscal
!      k1 = int(pfk)
!      fk1 = pfk - real(k1)
!      k1 = k1 + 1
!      
!      ! limit to defined mtf range
!      if (k1<1) then
!        k1 = 1
!        fk1 = 0.0
!      end if
!      if (k1>=nmtf) then
!        k1 = nmtf
!        fk1 = 0.0
!      end if
!      
!      vmtf = mtfdata(k1)*(1.0-fk1) + mtfdata(k1+1)*fk1
!      
!      cft(i,j) = cft(i,j) * vmtf
!    
!    end do
!  end do
!  
!  ! call PostDBGMessage("Applying MTF backward FFT.")
!  ! tranform to real space
!  call FSFFT(cft,nix,niy,'bac')
!  
!  ! copy data
!  do j=1, niy
!    do i=1, nix
!      img(i,j) = real(cft(i,j))
!    end do
!  end do
!  
!  return
!  
!end subroutine ApplyMTF
!
!
!!**********************************************************************!
!!
!! subroutine ApplyNTF
!!
!! Applies noise-tf to given noise image
!! (only if domtf flag is set in global parameters)
!!
!! INPUT: (global parameters in wavimgprm.f90)
!!   integer*4 :: nix, niy           = image dimensions
!!
!! IN/OUTPUT:  (global parameters in wavimgprm.f90)
!!   real*4 :: img(nix, niy)         = image data
!!
!subroutine ApplyNTF(img,nix,niy)
!
!  use wavimgprm
!  implicit none
!  
!  integer*4, intent(in) :: nix, niy
!  real*4, intent(inout) :: img(nix,niy)
!  
!  integer*4 :: i, j, i1, j1         ! iterators
!  integer*4 :: nnx, nny, nnx1, nny1 ! nyquist numbers
!  integer*4 :: k1                   ! mtf pos
!  real*4 :: px, py                  ! positions
!  real*4 :: px2, pr2                ! position squared radius
!  real*4 :: pf, pfk, pfim, pfkm     ! frequency interpolation pos
!  real*4 :: fk1                     ! mtf pos fraction
!  real*4 :: vntf                    ! mtf value to apply
!  
!  
!  external :: FSFFT            ! fourier transform, link FFTs.f
!  
!  ! initialize
!  nerr = 0
!  !   check mtf application request flag
!  if (domtf==0) return
!  cft = cmplx(0.0,0.0)
!  nnx = nix / 2
!  nnx1 = nnx - 1
!  nny = niy / 2
!  nny1 = nny - 1
!  pfim = real(max(nnx,nny))
!  pfkm = real(nmtf-1)
!  
!  ! call PostDBGMessage("Applying NTF forward FFT.")
!  ! copy data
!  do j=1, niy
!    do i=1, nix
!      cft(i,j) = cmplx(img(i,j),0.0)
!    end do
!  end do
!  
!  ! tranform to fourier space
!  call FSFFT(cft,nix,niy,'for')
!  
!  ! apply mtf to ft image data
!  do j=1, nix
!    ! get x-frequency
!    j1 = mod((j+nnx1),nix)-nnx
!    px = real(j1)
!    px2 = px*px
!    
!    do i=1, niy
!      ! get y-frequency
!      i1 = mod((i+nny1),niy)-nny
!      py = real(i1)
!      pr2 = px2 + py*py
!      
!      ! 2d-frequency modulus
!      pf = sqrt(pr2)
!      
!      ! determine scale to mtf 
!      pfk = pf / pfim * pfkm * mtfscal
!      k1 = int(pfk)
!      fk1 = pfk - real(k1)
!      k1 = k1 + 1
!      
!      ! limit to defined mtf range
!      if (k1<1) then
!        k1 = 1
!        fk1 = 0.0
!      end if
!      if (k1>=nmtf) then
!        k1 = nmtf
!        fk1 = 0.0
!      end if
!      
!      vntf = ntfdata(k1)*(1.0-fk1) + ntfdata(k1+1)*fk1
!      
!      cft(i,j) = cft(i,j) * vntf
!    
!    end do
!  end do
!  
!  ! call PostDBGMessage("Applying MTF backward FFT.")
!  ! tranform to real space
!  call FSFFT(cft,nix,niy,'bac')
!  
!  ! copy data
!  do j=1, niy
!    do i=1, nix
!      img(i,j) = real(cft(i,j))
!    end do
!  end do
!  
!  return
!  
!end subroutine ApplyNTF


!**********************************************************************!
!
! subroutine ApplyImgEnvelopes
!
! Applies several envelopes to a given image
! - suppressing linear focal kernel resonances, only if nsre==1
! - vib, only if dovib flag is set in global parameters
! - convert to mean number of electrons in vacuum
! - apply the shot noise
! - apply the detector conversion rate
! - mtf, only if domtf flag is set in global parameters
! - apply dark/readout noise
!
! I0
! |
! | Amplify to e-                             1.1
! |
! Ia
! |
! | *RES                                      2.1
! | *VIB                                      2.2
! | *MTF                                      2.3
! |
! Ib -- If(doint) -+----------+
! |                |          |
! |                | Get SN   | Get RN (cts)  3.1, 3.2
! |                |          |
! |                | *NTF     |               4.1
! |                |          |
! Ib               Isn        Irn
! |                |          |
! Ic = Ib    +     Isn        |               5.1
! |                           |
! | Amplify to cts            |
! |                           |
! Ica                         |               6.1
! |                           |
! Ici                         |               7.1
! |                           |
! I = Ici                   + Irn             8.1
!
!
! INPUT: (global parameters in wavimgprm.f90)
!   integer*4 :: nix, niy           = image dimensions
!
! IN/OUTPUT:  (global parameters in wavimgprm.f90)
!   real*4 :: img(nix, niy)         = image data
!
subroutine ApplyImgEnvelopes(img,nix,niy)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy
  real*4, intent(inout) :: img(nix,niy)
  
  integer*4 :: i, j, i1, j1         ! iterators
  integer*4 :: nalloc               ! nalloc allocation status
  integer*4 :: nnx, nny, nnx1, nny1 ! nyquist numbers
  integer*4 :: k1                   ! mtf pos
  real*4 :: gmax, gmaxx, gmaxy      ! image nyquist frequency
  real*4 :: px, py                  ! positions
  real*4 :: px2, py2, pr2, sd, cd   ! position squared radius
  real*4 :: v1, v2, v12, v22        ! vibration parameters
  real*4 :: pf, pfk, pfim, pfkm     ! frequency interpolation pos
  real*4 :: gx, gy, gx2, g2, g      ! spatial frequency values [1/nm]
  real*4 :: itogx, itogy            ! fourier space samplings
  real*4 :: fk1                     ! mtf pos fraction
  real*4 :: gsr0                    ! resonance suppressor centre [1/nm]
  real*4 :: gsrd                    ! resonance suppressor width [1/nm]
  real*4 :: vmtf, vntf              ! mtf & ntf values to apply
  real*4 :: vvib                    ! vibration damping value to apply
  real*4 :: vres                    ! resonance suppressor value to apply
  real*4 :: el_mean                 ! mean number of electrons per pixel in vacuum
  real*4 :: cpel                    ! counts per electron
  real*4 :: vtmp, vmean             ! temp. calculation value
  real*4 :: tx, ty                  ! beam tilt  in rad
  real*4 :: six, siy                ! image sampling rates x, y
  real*4, allocatable :: snimg(:,:), rnimg(:,:) ! noise images
  complex*8, allocatable :: cvib(:,:) ! vibration envelope
  
  
  external :: FSFFT            ! fourier transform, link FFTs.f
  integer*4, external :: PoissonRand
  real*4, external :: GaussRand, UniRand
  
  ! --------------------------------------------------------------------
  !
  ! Initialize
  !
  nerr = 0
  cft = cmplx(0.0,0.0)
  nnx = nix / 2
  nnx1 = nnx - 1
  nny = niy / 2
  nny1 = nny - 1
  pfim = real(max(nnx,nny))
  pfkm = real(nmtf-1)
  vmtf = 1.0
  vntf = 1.0
  vres = 1.0
  vvib = 1.0
  if (dofrm/=0) then
    six = simg
    siy = simg
  else
    six = swx
    siy = swy
  end if
  gmaxx = 0.5/six
  gmaxy = 0.5/siy
  gmax = max(gmaxx,gmaxy)
  itogx = 1.0/(six*real(nix))
  itogy = 1.0/(six*real(niy))
  v1 = vibamp
  v12 = v1*v1
  v2 = vibamp2
  v22 = v2*v2
  cd = COS(d2r*vibdir)
  sd = SIN(d2r*vibdir)
  el_mean = 1.0
  cpel = 1.0
  if (abs(el_conv)>0.0) cpel = abs(el_conv)
  el_mean = abs(iimean/cpel)
  tx = btx*0.001
  ty = bty*0.001
  !
  ! allocations
  allocate(snimg(nix,niy), rnimg(nix,niy), stat=nalloc)
  if (nalloc/=0) then
    call PostWarning("Allocation of noise image memory failed, switching off integer image calculation and calculating a floating point image.")
    doint=0
  else 
    snimg = 0.0
    rnimg = 0.0
  end if
  !
  ! vibration envelope
  if (dovib/=0) then
    allocate(cvib(nft,nft), stat=nalloc)
    if (nalloc/=0) then
      call PostWarning("Allocation of memory for vibration envelope failed. Switching off vibration.")
      dovib=0
    else 
      ! setup the vibration kernel
      cvib = cmplx(0.0,0.0)
      ! The vibration kernel is set up in real-space and transform later to Fourier space.
      ! This is required since the kernel can be very narrow which requires a correct
      ! aliasing. Due to my laziness I decided to do this in real-space.
      !
      vmean = 0.0
      do j=1, niy
        ! calculate the y-coordinate
        if (j<nny) then
          py = siy*real(j-1)
        else
          py = siy*real(j-1-niy)
        end if
        py2 = py * py
        do i=1, nix
        ! calculate the x-coordinate
          if (i<nnx) then
            px = six*real(i-1)
          else
            px = six*real(i-1-nix)
          end if
          px2 = px * px
          !
          ! calculate the convolution kernel
          ! using a gaussian of rms - half width
          vtmp = EXP( -((py2*v12 + px2*v22)*cd*cd +((px2*v12 + py2*v22)*sd &
                  & -2.*px*py*(v12-v22)*cd)*sd )/(2.*v12*v22) )
          cvib(i,j) = cmplx(vtmp ,0.0)
       end do
      end do
      ! tranform the kernel to Fourier space
      call FSFFT(cvib,nix,niy,'for')
      ! normalize
      vmean = real(cvib(1,1))
      cvib = cvib / vmean
      !
    end if
  end if
  !
  !
  ! check resonance damping
  if (sqrt(2.0)*gmax<gres1 .or. gres1<=0.0) then
    nsre=0
    !call PostDBGMessage("1st linear focal kernel resonance is outsides the image spectrum, resonance suppression deactivated")
  end if
  !
  ! 
  ! pre-calculate focal kernel resonances
  if (nsre/=0) then
    if (gres1<2.1*ginfo) then
      call PostWarning("Resonance of focal kernel is less than twice the information limit. Some non-linear image terms will be lost.")
      gsr0 = gres1*0.9
      gsrd = gres1*0.05
    else
      gsr0 = 0.5*(gres1+2.0*ginfo)
      gsrd = 0.2*(gres1-2.0*ginfo)
      write(unit=sdbgmsg,fmt=*) "image post-proc: reson. damp. prm: ", gsr0, gsrd
      call PostDBGMessage(trim(sdbgmsg))
    end if
  end if
  
  ! --------------------------------------------------------------------
  !
  ! 1) SCALE THE INPUT IMAGE TO THE AVERAGE VACUUM ELECTRON NUMBER
  !
  write(unit=sdbgmsg,fmt='(G12.4)') el_mean
  call PostDBGMessage("Image post-processing: 1) scaling to avg. input e-: "&
     & //trim(adjustl(sdbgmsg)))
  img = img * el_mean
  !
  !
  
  
  ! --------------------------------------------------------------------
  !
  ! 2) CONTRAST DAMPING
  !
  if (domtf==0 .and. dovib==0 .and. nsre==0) then
    call PostDBGMessage("Image post-processing: 2) No image contrast damping applied")
  else
    call PostDBGMessage("Image post-processing: 2) Applying image contrast damping")
  end if
  !
  ! prepare work arrays
  ! - avg. image work array
  do j=1, niy
    do i=1, nix
      cft(i,j) = cmplx(img(i,j),0.0)
    end do
  end do
  ! transform to fourier space
  call FSFFT(cft,nix,niy,'for')
  
  if (nsre/=0 .or. dovib/=0 .or. domtf/=0 ) then
  
    if (nsre/=0) then
    !
    ! 2.1) DAMPING OF FOCUS KERNAL RESONANCES
    !
    call PostDBGMessage("Image post-processing: 2.1) Suppressing focal kernel resonances")
    end if
    !
    if (dovib/=0) then
    !
    ! 2.2) DAMPING BY VIBRATION ENVELOPE
    !
    call PostDBGMessage("Image post-processing: 2.2) Applying vibration damping")
    end if
    !
    if (domtf/=0) then
    !
    ! 2.3) DAMPING BY DETECTOR MTF
    !
    call PostDBGMessage("Image post-processing: 2.3) Applying detector mtf to the average image")
    end if
    !
    !
    vres = 1.0
    vmtf = 1.0
    ! apply envelopes to ft of avg. image data
    do j=1, nix
      ! get x-frequency
      j1 = mod((j+nnx1),nix)-nnx
      px = real(j1)
      gx = px*itogx
      gx2 = gx*gx
      
      px = px / real(nnx)
      px2 = px*px

      do i=1, niy
        ! get y-frequency
        i1 = mod((i+nny1),niy)-nny
        py = real(i1)
        gy = py*itogy
        g2 = gx2 + gy*gy
        
        py = py/real(nny)
        pr2 = px2 + py*py
      
        ! 2d-frequency modulus
        pf = sqrt(pr2)
        g = sqrt(g2)
      
        if (domtf/=0) then
          ! determine scale to mtf 
          pfk = pf * pfkm * mtfscal
          k1 = int(pfk)
          fk1 = pfk - real(k1)
          k1 = k1 + 1
          ! limit to defined mtf range
          if (k1<1) then
            k1 = 1
            fk1 = 0.0
          end if
          if (k1>=nmtf) then
            k1 = nmtf
            fk1 = 0.0
          end if
          ! interpolate mtf (linear)
          vmtf = mtfdata(k1)*(1.0-fk1) + mtfdata(k1+1)*fk1
        end if
      
        if (nsre/=0) then ! get resonance damping
          vres = 0.5-0.5*tanh((g-gsr0)/gsrd)
        end if ! (nsre/=0)
      
        ! apply damping terms to FT of image
        if (dovib/=0) then
          cft(i,j) = cft(i,j) * conjg(cvib(i,j)) * vres * vmtf
        else
          cft(i,j) = cft(i,j) * vres * vmtf
        end if
        ! cft2(i,j) = cft(i,j)  ! copy to cft2 for noise image generation

      end do
    end do

    ! tranform to real space
    call FSFFT(cft,nix,niy,'bac')
    ! get dampened image
    do j=1, niy
      do i=1, nix
        img(i,j) = real(cft(i,j))
      end do
    end do
  
  end if ! (nsre, dovib)
  
  
  ! --------------------------------------------------------------------
  !
  ! 3) NOISE IMAGE GENERATION
  !
  if (doint/=0) then ! prepare integer image noise generation
    !
    ! init noise image
    snimg = 0.0
    !
!    ! get the pre-dampened avg. image to create the noise images
!    call FSFFT(cft2,nix,niy,'bac')
!    ! copy data
!    do j=1, niy
!      do i=1, nix
!        snimg(i,j) = real(cft2(i,j))
!      end do
!    end do
    !
    ! 3.1) GENERATE SHOT NOISE IMAGE FROM AVG NUMBER OF ELECTRONS PER PIXEL
    !
    call PostDBGMessage("Image post-processing: 3.1) shot noise generation")
    do j=1, niy
      do i=1, nix
        vmean = img(i,j) ! copy value of the pre-dampened image
        if (vmean>20.0) then
          vtmp = sqrt(abs(vmean))*GaussRand()
        else
          vtmp = real(PoissonRand(vmean))-vmean
        end if
        snimg(i,j) = vtmp
      end do
    end do
    !
    ! prepare work array in RS
    cft(:,:) = cmplx(0.0,0.0)
    ! - noise image work array
    do j=1, niy
      do i=1, nix
        cft(i,j) = cmplx(snimg(i,j),0.0)
      end do
    end do
    ! transform to Fourier space
    call FSFFT(cft,nix,niy,'for')
    !
    !
    ! 3.2) GENERATE THE READ-OUT NOISE IMAGE
    !
    write(unit=sdbgmsg,fmt='(G12.4)') abs(dark_noise)
    call PostDBGMessage("Image post-processing: 3.2) readout noise generation, RMS amplitude = "&
     & //trim(adjustl(sdbgmsg)))
    do j=1, niy
      do i=1, nix
        vtmp = abs(dark_noise)*GaussRand()
        rnimg(i,j) = vtmp
      end do
    end do
    !
    !
  end if
  
  
  ! --------------------------------------------------------------------
  !
  ! 4) NTF
  !
  if ( domtf/=0 .and. doint/=0 ) then ! apply the ntf here
  
    call PostDBGMessage("Image post-processing: 4.1) applying NTF")
    ! apply envelopes to ft image data
    do j=1, nix
      ! get x-frequency
      j1 = mod((j+nnx1),nix)-nnx
      px = real(j1)
      gx = px*itogx
      gx2 = gx*gx
      
      px = px / real(nnx)
      px2 = px*px
    
      do i=1, niy
        ! get y-frequency
        i1 = mod((i+nny1),niy)-nny
        py = real(i1)
        gy = py*itogy
        
        py = py / real(nny)
        pr2 = px2 + py*py
        
        g2 = gx2 + gy*gy
      
        ! 2d-frequency modulus
        pf = sqrt(pr2)
        g = sqrt(g2)
      
        ! determine scale to mtf 
        pfk = pf * pfkm * mtfscal
        k1 = int(pfk)
        fk1 = pfk - real(k1)
        k1 = k1 + 1
      
        ! limit to defined mtf range
        if (k1<1) then
          k1 = 1
          fk1 = 0.0
        end if
        if (k1>=nmtf) then
          k1 = nmtf
          fk1 = 0.0
        end if
        
        ! interpolate ntf (linear)
        vntf = ntfdata(k1)*(1.0-fk1) + ntfdata(k1+1)*fk1
        ! apply ntf
        cft(i,j) = cft(i,j) * vntf
      
      end do
    end do
    
  
    
    ! ------------------------------------------------------------------
    !
    ! 5) Adding noise and average image after mtf and ntf convolution
    !
    call PostDBGMessage("Image post-processing: 5.1) adding shot noise")
    ! tranform to real space
    call FSFFT(cft,nix,niy,'bac')
    ! copy data
    do j=1, niy
      do i=1, nix
        img(i,j) = img(i,j) + real(cft(i,j))
      end do
    end do
    !
  
  end if !(domtf and doint)
  
  
  ! --------------------------------------------------------------------
  !
  ! 6) SCALE THE IMAGE TO THE AVERAGE VACUUM COUNT NUMBER
  !
  write(unit=sdbgmsg,fmt='(G12.4)') cpel
  call PostDBGMessage("Image post-processing: 6) scaling from e- to cts, factor: "&
     & //trim(adjustl(sdbgmsg)))
  img = img * cpel
  !
  !
  
  ! --------------------------------------------------------------------
  !
  ! 7) APPLY THE PIXEL INTEGRATION TRANSFER FUNCTION
  !
  if (domtf/=0) then ! trigger pixelation with the use of an MTF
    call PostDBGMessage("Image post-processing: 7.1) pixel area integration")
    do j=1, niy
      do i=1, nix
        cft(i,j) = cmplx(img(i,j),0.0)
      end do
    end do
    ! transform to Fourier space
    call FSFFT(cft,nix,niy,'for')
    ! apply sync damping
    vtmp = 1.5707963 ! pi / 2
    do j=1, nix
      ! get x-frequency
      j1 = mod((j+nnx1),nix)-nnx
      px = real(j1) / real(nnx)
      v1 = sin(vtmp*px+0.0001) / (vtmp*px + 0.0001)
      do i=1, niy
        ! get y-frequency
        i1 = mod((i+nny1),niy)-nny
        py = real(i1) / real(nny)
        v2 = sin(vtmp*py+0.0001) / (vtmp*py + 0.0001)
        cft(i,j) = cft(i,j) * v1 * v2
      end do
    end do
    ! transform to real space
    call FSFFT(cft,nix,niy,'bac')
    ! copy data
    do j=1, niy
      do i=1, nix
        img(i,j) = real(cft(i,j))
      end do
    end do
  end if
  !
  !
  
  
  ! --------------------------------------------------------------------
  !
  ! 8) ADD READ-OUT NOISE
  !
  if (doint) then ! apply readout noise
    call PostDBGMessage("Image post-processing: 8.1) adding readout noise")
    do j=1, niy
      do i=1, nix
        img(i,j) = img(i,j) + rnimg(i,j)
      end do
    end do
  end if
  
  
  !
  ! deallocations
  if (allocated(snimg)) deallocate(snimg, stat=nalloc)
  if (allocated(rnimg)) deallocate(rnimg, stat=nalloc)
  if (allocated(cvib)) deallocate(cvib, stat=nalloc)
  
  
  call PostDBGMessage("Image post-processing: finished.")
  
  return
  
end subroutine ApplyImgEnvelopes




!**********************************************************************!
!
! subroutine SaveIntegerImage32
!
! creates a 32-bit integer image 
!
! INPUT: (global parameters in wavimgprm.f90)
!   integer*4 :: nix, niy           = image dimensions
!   real*4 :: rimg(nix, niy)        = image data
!   character(len=*) :: sfilename   = image file name string
!
subroutine SaveIntegerImage32(sfilename,nix,niy,rimg)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy
  real*4, intent(in) :: rimg(nix,niy)
  character(len=*), intent(in) :: sfilename
  
  real*4, external :: GaussRand
  external :: savedatai4 ! (sfile, n, a, nerr)
  
  integer*4 :: i, j                             ! iterators
  integer*4, allocatable :: i32img(:,:)         ! output
  
  allocate(i32img(nix,niy),stat=nerr)
  if (nerr/=0) return
  
  do j=1, niy
    do i=1, nix
      i32img(i,j) = nint(abs(rimg(i,j)),kind=4)
    end do
  end do
  
  call savedatai4(trim(sfilename),nix*niy,i32img,nerr)
  
  deallocate(i32img,stat=nerr)
  
  return

end subroutine SaveIntegerImage32


!**********************************************************************!
!
! subroutine SaveIntegerImage16
!
! creates a 16-bit integer image 
!
! INPUT: (global parameters in wavimgprm.f90)
!   integer*4 :: nix, niy           = image dimensions
!   real*4 :: rimg(nix, niy)        = image data
!   character(len=*) :: sfilename   = image file name string
!
subroutine SaveIntegerImage16(sfilename,nix,niy,rimg)

  use wavimgprm
  implicit none
  
  integer*4, intent(in) :: nix, niy
  real*4, intent(in) :: rimg(nix,niy)
  character(len=*), intent(in) :: sfilename
  
  real*4, external :: GaussRand
  external :: savedatai2 ! (sfile, n, a, nerr)
  
  integer*4 :: i, j                             ! iterators
  integer*2, allocatable :: i16img(:,:)         ! output
  
  allocate(i16img(nix,niy),stat=nerr)
  if (nerr/=0) return
  
  do j=1, niy
    do i=1, nix
      i16img(i,j) = nint(abs(rimg(i,j)),kind=2)
    end do
  end do
  
  call savedatai2(trim(sfilename),nix*niy,i16img,nerr)
  
  deallocate(i16img,stat=nerr)
  
  return

end subroutine SaveIntegerImage16


!**********************************************************************!
!
! subroutine getwaveamp
!
! caluclates the amplitude of a given wave
!
! INPUT:
!   complex*8 :: w(nx,ny)           = wave data
!   integer*4 :: nx, ny             = wave dimensions
!
! IN/OUTPUT:
!   real*4 :: a(nx, ny)             = amplitude data
!
subroutine getwaveamp(w,a,nx,ny)
  implicit none
  integer*4, intent(in) :: nx, ny
  complex*8, intent(in) :: w(nx,ny)
  real*4, intent(inout) :: a(nx,ny)
  if (nx*ny>0) a = cabs(w)
  return
end subroutine getwaveamp

!**********************************************************************!
!
! subroutine getwavephase
!
! caluclates the phase of a given wave
!
! INPUT:
!   complex*8 :: w(nx,ny)           = wave data
!   integer*4 :: nx, ny             = wave dimensions
!
! IN/OUTPUT:
!   real*4 :: a(nx, ny)             = phase data
!
subroutine getwavephase(w,a,nx,ny)
  implicit none
  integer*4, intent(in) :: nx, ny
  complex*8, intent(in) :: w(nx,ny)
  real*4, intent(inout) :: a(nx,ny)
  integer*4 :: i,j
  complex*8 :: cval
  if (nx*ny<=0) return
  do j=1,ny
    do i=1,nx
      cval = w(i,j)
      a(i,j) = atan2(imag(cval),real(cval))
    end do
  end do
  return
end subroutine getwavephase



!**********************************************************************!
!
! subroutine setloopvariable
!
! caluclates value of loop variables and applies change
!
! INPUT:
!   none
!
! IN/OUTPUT:
!   none
!
subroutine setloopvariable()
  use AberrationFunctions
  use wavimgprm
  implicit none
  integer*4 :: i, j, k, looptype
  real*4 :: loopval, tmp, vazi, vamp, v1, v2, vv1, vv2
  
  ! check if a loop is active
  if (nloop<=0) then
    call PostDBGMessage("No parameter loop defined.") 
    return
  end if
  
  ! Reset all aberrations to their backup value.
  AF_wa = AF_wabk
  ! The loops will apply changes to AF_wa
  siwav = trim(swavfile) ! re-init the modified wave file name
  
  ! loop through the possible loops
  do i=1, nloop
    write(unit=sdbgmsg,fmt='(A,I1)') "Setting loop variable #",i
    call PostDBGMessage(trim(sdbgmsg))
    write(unit=sdbgmsg,fmt='(A,I4,A,I4)') "- loop index = ",lpidx(i)," of ", lpsz(i)
    call PostDBGMessage(trim(sdbgmsg))
    ! get fractional loop position
    tmp = real(lpidx(i)-1)/real(lpsz(i)-1)
    ! determine the loop variation form lpvf(i)
    ! meaning of the loop variable parameters depending on the variation form:
    ! - form = RAMP (1):
    !   - lpv0 = start value
    !   - lpvd = step size
    !   - lpv1 = stop value or upper limit (loop may stop earlier to keep up the step size)
    ! - form = OSZI (2):
    !   - lpv0 = amplitude
    !   - lpvd = total number of samples per loop = lpsz
    !   - lpv1 = frequency ( number of full oscillations during the loop)
    ! - form = LIST (3):
    !   - parameters are taken from a pre-loaded list
    looptype = lpvf(i)
    ! - handle the case of a wave-function loop, where the variation form MUST BE = 1 (ramp)
    if (lpcl(i)==3) looptype = 1
    select case (looptype)
      case (1) ! ramp
        loopval = lpv0(i) + lpvd(i)*tmp*real(lpsz(i)-1)
        vamp = loopval ! amplitude, azimuth phase must be defined later
        write(unit=sdbgmsg,fmt='(A,G11.3)') "- current loop variable value = ",loopval
        call PostDBGMessage(trim(sdbgmsg))
      case (2) ! oscil phase
        loopval = twopi*tmp*lpv1(i)
        vazi = loopval ! azimuth phase
        vamp = lpv0(i) ! amplitude
        write(unit=sdbgmsg,fmt='(A,G11.3)') "- current loop variable value = ",loopval
        call PostDBGMessage(trim(sdbgmsg))
    end select
    !
    !
    ! switch between supported loop variable classes lpcl(i)
    ! and apply the current loop variable value
    SEL_SPECIES: select case (lpcl(i))
      !
      case (1) ! aberration, identified by index lpvr(:)+1
        if (looptype==1.or.looptype==2) then
          call AF_GetAberrationBackup(lpvr(i)+1,v1,v2) ! get offset value from backup memory (input data)
          vazi = atan2(v2,v1) ! determine offset azimuth
          v1 = vamp * cos( vazi )
          v2 = vamp * sin( vazi )
          ! set new aberration values
          ! this way of variation will override the input completely
          call AF_SetAberration(lpvr(i)+1,v1,v2)
          write(unit=sdbgmsg,fmt='(A,I2,A,G11.3,A,G11.3,A)') &
     &          "- set aberration #",lpvr(i)," to value (",v1,",",v2,") nm"
          call PostDBGMessage(trim(sdbgmsg))
        end if
        if (looptype==3) then
          do j=1, lpali(0,i)
            ! get the current aberration value
            call AF_GetAberration(lpali(j,i)+1,v1,v2) ! get offset value from current value (could be changed by other loops)
            ! get the aberration change for this loop
            k = 2*(j-1) + 1 ! list index
            vv1 = lpalv(k  , lpidx(i), i)  ! aberration x value of set # lpidx(i) in loop # i
            vv2 = lpalv(k+1, lpidx(i), i)  ! aberration y value of set # lpidx(i) in loop # i
            call AF_SetAberration( lpali(j,i)+1, & ! index of aberration # j in loop # i
     &                             v1+vv1, v2+vv2 ) ! aberration values (offset + variation)
            write(unit=sdbgmsg,fmt='(A,I2,A,G11.3,A,G11.3,A)') &
     &          "- set aberration #",lpali(j,i)," to value (",v1+vv1,",",v2+vv2,") nm"
            call PostDBGMessage(trim(sdbgmsg))
          end do
        end if
        
      !
      case (2) ! coherence parameter
        SEL_SPECIES_ID: select case (lpvr(i)) ! the coherence variable is identified by lpvr(i)
          case (1) ! temporal coherence
            if (looptype==1) vazi = 0.0 ! determine offset azimuth
            fs = vamp * cos(vazi)
            write(unit=sdbgmsg,fmt='(A,G11.3,A)') "- set focus-spread to ",fs," nm"
            call PostDBGMessage(trim(sdbgmsg))
          case (2) ! spatial coherence
            if (looptype==1) vazi = 0.0 ! determine offset azimuth
            sc = vamp*0.001 * cos(vazi) ! sc is in radian, while the loop variable is defined in milliradian
            write(unit=sdbgmsg,fmt='(A,G11.3,A)') "- set convergence angle to ",sc," rad"
            call PostDBGMessage(trim(sdbgmsg))
        end select SEL_SPECIES_ID
      !
      case (3) ! wave function
        ! the wave function to load is defined by a numerical index and a prefix sub-string
        ! the prefix sub-string MUST BE part of the wave function filename swavfile and is saved in lpid(:)
        ! the numerical index is determined from the loop parameters and used with leading zeros up to 3 chars
        ! the numerical index string is inserted directly after the sub-string
        ! e.g. swavfile = "Material_px001_py008_sl.wav", sub-string = "_sl", wave number = 12 -> "Material_px001_py008_sl012.wav"
        !      swavfile = "Material_px_py008_sl001.wav", sub-string = "_px", wave number = 3  -> "Material_px003_py008_sl001.wav"
        if (looptype==1) vazi = 0.0 ! determine offset azimuth
        j = nint( vamp * cos(vazi) ) ! identify the wave function index using the loopvariable without further math
        
        ! determine max. wave file index
        if (looptype==1) k = nint(lpv1(i))
        if (looptype==2) k = nint(lpv0(i))
        call addloopwavename(j,k,lpid(i),siwav)
    end select SEL_SPECIES
    !
  end do ! loop over active variable loops
  
  return
end subroutine setloopvariable



!**********************************************************************!
!
! subroutine addloopimagename
!
! adds indices to the current file name
! inserted before the file extension
! the inner loop (lowest loop index) is added last
!
! INPUT:
!
! IN/OUTPUT:
!   character(len=*) :: sname           = input and output name
!
subroutine addloopimagename(sname)
  use wavimgprm
  implicit none
  character(len=*), intent(inout) :: sname
  integer*4 :: ndot, i, nilen, idx
  character(len=600) :: snum
  if (nloop<=0) return ! no loop
  do i=nloop, 1, -1
    idx = lpidx(i)
    stmp = trim(sname)
    ndot = INDEX(trim(stmp),'.', BACK = .TRUE.)
    write(unit=snum,fmt='(I)') lpsz(i)
    nilen = max( 3, len_trim(adjustl(snum)) )
    if (ndot==0) then
      write(unit=sname,fmt='(A,"_",I<nilen>.<nilen>)') trim(stmp),idx
    else
      write(unit=sname,fmt='(A,"_",I<nilen>.<nilen>,A)') stmp(1:ndot-1),idx,stmp(ndot:len_trim(stmp))
    end if
  end do
  return
end subroutine addloopimagename




!**********************************************************************!
!
! subroutine addloopmapname
!
! adds indices to the current file name
! inserted before the file extension
! the inner loop (lowest loop index) is added last
!
! this special implementation is meant for map images
! the first two loops will be ignored for the file name creation
!
! INPUT:
!
! IN/OUTPUT:
!   character(len=*) :: sname           = input and output name
!
subroutine addloopmapname(sname)
  use wavimgprm
  implicit none
  character(len=*), intent(inout) :: sname
  integer*4 :: ndot, i, nilen, idx
  character(len=600) :: snum
  if (nloop<=2) return ! no loop around the map
  do i=nloop, 3, -1
    idx = lpidx(i)
    stmp = trim(sname)
    ndot = INDEX(trim(stmp),'.', BACK = .TRUE.)
    write(unit=snum,fmt='(I)') lpsz(i)
    nilen = max( 3, len_trim(adjustl(snum)) )
    if (ndot==0) then
      write(unit=sname,fmt='(A,"_",I<nilen>.<nilen>)') trim(stmp),idx
    else
      write(unit=sname,fmt='(A,"_",I<nilen>.<nilen>,A)') stmp(1:ndot-1),idx,stmp(ndot:len_trim(stmp))
    end if
  end do
  return
end subroutine addloopmapname




!**********************************************************************!
!
! subroutine addloopavgname
!
! adds indices to the current file name
! inserted before the file extension
! the inner loop (lowest loop index) is added last
!
! this special implementation is meant for average images
! the first loop will be ignored for the file name creation
!
! INPUT:
!
! IN/OUTPUT:
!   character(len=*) :: sname           = input and output name
!
subroutine addloopavgname(sname)
  use wavimgprm
  implicit none
  character(len=*), intent(inout) :: sname
  integer*4 :: ndot, i, nilen, idx
  character(len=600) :: snum
  if (nloop<=1) return ! no loop around the map
  do i=nloop, 2, -1
    idx = lpidx(i)
    stmp = trim(sname)
    ndot = INDEX(trim(stmp),'.', BACK = .TRUE.)
    write(unit=snum,fmt='(I)') lpsz(i)
    nilen = max( 3, len_trim(adjustl(snum)) )
    if (ndot==0) then
      write(unit=sname,fmt='(A,"_",I<nilen>.<nilen>)') trim(stmp),idx
    else
      write(unit=sname,fmt='(A,"_",I<nilen>.<nilen>,A)') stmp(1:ndot-1),idx,stmp(ndot:len_trim(stmp))
    end if
  end do
  return
end subroutine addloopavgname




!**********************************************************************!
!
! subroutine addloopwavename
!
! adds a index to the current file name
! inserted after the substring
!
! INPUT:
!   integer*4 :: idx                = loop index
!   integer*4 :: idxmax             = loop index maximum
!   character(len=*) :: ssub        = sub string preceeding index location
!
! IN/OUTPUT:
!   character(lem=*) :: sname       = input and output string
!
subroutine addloopwavename(idx,idxmax,ssub,sname)
  use wavimgprm
  implicit none
  integer*4, intent(in) :: idx, idxmax
  character(len=*), intent(inout) :: ssub
  character(len=*), intent(inout) :: sname
  integer*4 :: ninpos, ninlen, nsublen, j
  ! determine the length of the inserted index string (minimum is 3)
  write(unit=stmp,fmt='(I)') idxmax
  ninlen = max( 3, len_trim(adjustl(stmp)) )
  nsublen = len_trim(ssub)
  stmp = trim(sname)
  ninpos = 0
  if (nsublen>0) then
    ninpos = INDEX(trim(stmp),trim(ssub), BACK = .TRUE.)
  end if
  if (ninpos==0) then
    write(unit=sname,fmt='(A,I<ninlen>.<ninlen>)') trim(stmp),idx
  else
    j = ninpos+nsublen
    write(unit=sname,fmt='(A,I<ninlen>.<ninlen>,A)') stmp(1:j-1),idx,stmp(j:len_trim(stmp))
  end if
  return
end subroutine addloopwavename



!**********************************************************************!
!
! Increases the loop indices for the next loop run.
! The routine will manipulate lpidx by raising the innermost loop index
!   lpidx(1) by one, wrapping inner loops when the loop size is reached
!
subroutine nextloopindex
  use wavimgprm
  implicit none
  integer*4 :: i
  !
  if (nloop>0) then
    ! there are active loops
    lpidx(1) = lpidx(1) + 1 ! raise innermost loop index by 1
    iloopstep = iloopstep + 1 ! increment loop step counter
    if (nloop>1) then ! check for outer loops
      do i=1, nloop-1 ! loop through all inner loops and check index range
        if (lpidx(i)>lpsz(i)) then ! check for index range exceed
          lpidx(i) = 1 ! wrap back to one ...
          lpidx(i+1) = lpidx(i+1) + 1 ! and increase next loop
        end if
      end do
    end if
  end if
  return
end subroutine nextloopindex


logical function exitprmloop
  use wavimgprm
  implicit none
  integer*4 :: i
  logical :: doexit
  doexit = .TRUE. ! exit by default / also if no loop is active
  if (nloop>0) then ! check specified loops
    do i=1, nloop ! let i go over all loops
      if (lpidx(i)<lpsz(i)) then ! check if one of the loop indices is less than the loop size
        doexit = .FALSE. ! there is still something to loop through
        exit ! no need for further checks
      end if
    end do
  end if
  exitprmloop = doexit ! return the function result
end function exitprmloop
