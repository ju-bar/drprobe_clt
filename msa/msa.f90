!**********************************************************************
!
!  PROGRAM: msa, version 0.85
!  FILE: msa.f90
!  PURPOSE:  Entry point for the console application MSA.
!            Multislice calculation for electron diffraction
!
!**********************************************************************
!                                                                      
!   Date: 2018-12-05                                                   
!                                                                      
!   Author: Juri Barthel                                               
!           Ernst Ruska-Centre                                         
!           Forschungszentrum Jülich GmbH, 52425 Jülich, Germany       
!           RWTH Aachen University, 52074 Aachen, Germany              
!                                                                      
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
!
! When using this program, please consider to cite the following
! article as an appreciation of the many hours of work put into
! this project:
!
! J. Barthel, "Dr. Probe - A software for high-resolution STEM
! image simulation", Ultramicroscopy 193 (2018) 1-11.
! doi: 10.1016/j.ultramic.2018.06.003
!
!----------------------------------------------------------------------



!**********************************************************************
!* MAIN ***************************************************************
!**********************************************************************




!**********************************************************************
!**********************************************************************
program msa

! ------------
! link modules    
  use IFPORT
  use EMSdata
  use MSAparams
  use STEMfunctions
  use MultiSlice
    
  implicit none
  
! ------------
! declare constants and variables
  integer :: nerr, prm_px, prm_py, lpx, lpy, mpx, mpy, i, j
  
  external :: InitRand
  
  call InitRand()
  
!  call CheckLicense()

! ------------
! INITIALIZE MSP MODULE
  call MSP_INIT()
  call EMS_INIT()
  MSP_callApp =                   "[msa] MultiSlice Algorithm"
  MSP_verApp  =                   "0.85  64-bit  -  2018 Dec   5  -"
  MSP_authApp =                   "Dr. J. Barthel, ju.barthel@fz-juelich.de"
! GET COMMAND LINE ARGUMENTS
  call parsecommandline()
! INTRODUCE YOURSELF
  call Introduce()
! ------------  

! ------------
  call MSP_INITCL()
  ! save pixel parameters
  prm_px = MSP_ScanPixelX
  prm_py = MSP_ScanPixelY
! ------------

! ------------  
  if (MSP_ctemmode==1) then
    call PostSureMessage("")
    call PostMessage( " +---------------------------------------------------+" )
    call PostMessage( " | Running in ctem mode.                             |" )
    call PostMessage( " +---------------------------------------------------+" )
  else
    call PostSureMessage("")
    call PostMessage( " +---------------------------------------------------+" )
    call PostMessage( " | Running in stem mode.                             |" )
    call PostMessage( " +---------------------------------------------------+" )
  end if
  call PostMessage( "" )
! ------------


! ------------
! VORTEX PROBE (added 2017-10-17 JB as hidden option)
  if (MSP_Vortex/=0) then
    if (0==MSP_ctemmode) then
      write(unit=MSP_stmp,fmt=*) MSP_Vortex
      call PostMessage("Using a vortex probe with an orbital angular momentum of "// &
         & trim(adjustl(MSP_stmp))//".")
    else
      call PostMessage("Vortex input is ignored in CTEM mode.")
      MSP_Vortex = 0
    end if
  end if
! ------------

  
! ------------
! LOAD PARAMETER FILE AND SET ALL RELEVANT PARAMETERS IN MODULES
  call STF_INIT() ! required for loading aberrations
! load parameters
  call PostSureMessage("Loading parameter file.")
  call LoadParameters(MSP_prmfile)
  call PostRuntime("main parameter input finished", 1)
! ------------

! ------------
  if (MSP_ExplicitPSC/=0 .and. MSP_ctemmode/=0) then
    MSP_ExplicitPSC = 0
    call PostWarning("Explicit partial coherence averaging is not available in CTEM mode.")
    call PostSureMessage("- Option /epc is ignored.")
  end if
  if (MSP_ExplicitPSC/=0 .and. MSP_ApplyPSC/=0) then
    MSP_ExplicitPSC = 0
    call PostWarning("Two conflicting partial coherence modes requested.")
    call PostSureMessage("- Explicit partial coherence averaging deactivated.")
    call PostSureMessage("- Input image convolution is used.")
  end if
! ------------

! ------------
! early shortcut to source convolution JB 2018-01-25
  if (MSP_use_extsrcprm) then ! update the source radius parameter using the external command line input
    STF_srcradius = MSP_extsrcrad
    write(unit=MSP_stmp,fmt=*) MSP_extsrcrad
    call PostMessage("- overriding source size by command-line parameter: "// &
       & trim(adjustl(MSP_stmp))//" nm.")
  end  if
  if (MSP_ApplyPSC/=0) then ! jump to convolution ?
    ! init STF module
    nerr = STF_err_num
    call PostMessage("Initialize module for wave modelling.")
    call STF_INIT_ALLOC()
    if (STF_err_num/=nerr) then
      call CriticalError("Module init failed.")
    end if
    call PostMessage("Initialize module for multislice calculations.")
    nerr = MS_err_num
    call MS_INIT()
    if (MS_err_num/=nerr) then
      call CriticalError("Module init failed.")
    end if
    goto 199
  end if
! ------------

! ------------
  if (MSP_ExplicitPSC/=0 .and. MSP_FL_varcalc<10) then 
    call PostWarning("Very low number of repeats for explicit partial coherence averaging.")
    call PostSureMessage("- Inaccurate results expected.")
    call PostSureMessage("- Increase the number of frozen-lattice averages per scan pixel.")
  end if
! ------------

! ------------
! probe image calculation warning
  if ((1==MSP_pimgmode .or. 1==MSP_pdifmode) .and. MSP_ctemmode/=0) then
    call PostWarning("Probe intensity integration is not available in CTEM mode.")
    if (1==MSP_pimgmode) call PostSureMessage("- Option /pimg is ignored.")
    if (1==MSP_pdifmode) call PostSureMessage("- Option /pdif is ignored.")
    MSP_pimgmode = 0
    MSP_pdifmode = 0
  end if
  if (1==MSP_pimgmode .or. 1==MSP_pdifmode) then
    call PostWarning("Performing probe intensity integrations, increased computation time!.")
    MS_pint_export = 1 ! This triggers the multislice into the wave export routine
  end if
! ------------
  
! ------------
  call PostMessage("Pre-Allocating memory for slice data.")
  nerr = MS_err_num
  call MSP_PREALLOCPGR(MS_slicenum,MSP_FL_varnum,nerr)
  if (MS_err_num/=nerr) call CriticalError("Memory allocation failed.")
  call PostMessage("Extracting size parameters from first sli-file.")
  call SetGlobalCellParams()
  FFT_BOUND = max(FFT_BOUND_MIN, 2**CEILING( LOG( real( max(MS_dimx,MS_dimy) ) )/LOG(2.0) )) ! next 2^N above max(nx,ny) or 128
  if (FFT_BOUND>FFT_BOUND_MAX) call CriticalError("FFT plan exceeds maximum size (8192).")
  FFT_NYQ = FFT_BOUND/2
  STF_FFT_BOUND = FFT_BOUND
  STF_FFT_NYQ = FFT_NYQ 
  write(unit=MSP_stmp,fmt='(A,I5,A,I5)') "- wave function size: ", MS_dimx, " x ", MS_dimy
  call PostMessage( trim(MSP_stmp) )
  write(unit=MSP_stmp,fmt='(A,I5)') "- FFT working field size: ", FFT_BOUND
  call PostMessage( trim(MSP_stmp) )
  ! allocate memory for phasegratings, size: MSP_dimcellx * MSP_dimcelly
  call PostMessage("Allocating memory for slice data.")
  call MSP_ALLOCPGR(MSP_dimcellx,MSP_dimcelly,nerr)
  if (MS_err_num/=nerr) then
    call CriticalError("Memory allocation failed.")
  end if
  call PostRuntime("parameter input finished", 1)
! ------------


! ------------
! INITIALIZE WAVE GENERATION MODULE
   nerr = STF_err_num
  call PostMessage("Initialize module for wave modelling.")
  call STF_INIT_ALLOC()
  if (STF_err_num/=nerr) then
    call CriticalError("Module init failed.")
  end if
  MS_nslid = MSP_nslid
  if (MSP_use_extdefocus/=0) then ! update the defocus using the external command line parameter
    call STF_SetAberration(2,MSP_extdefocus,0.0)
    write(unit=MSP_stmp,fmt=*) MSP_extdefocus
    call PostMessage("- overriding defocus by command-line parameter: "// &
       & trim(adjustl(MSP_stmp))//" nm.")
  end if
  if (MSP_use_extot/=0) then ! update object tilt parameters from command line input
    MS_objtiltx = MSP_OTExX
    MS_objtilty = MSP_OTExY
    write(unit=MSP_stmp,fmt=*) MSP_OTExX
    write(unit=MSP_stmp2,fmt=*) MSP_OTExY
    call PostMessage("- overriding object tilt by command-line parameter: otx="// &
       & trim(adjustl(MSP_stmp))//", oty="//trim(adjustl(MSP_stmp2))//" deg.")
  end if
  if (MSP_ctemmode/=0) then ! by default switch the wave export ON in CTEM mode
    if (MS_wave_export<1) MS_wave_export = 1 ! set to rea-space export by default
    MS_wave_filenm = MSP_outfile
    MS_wave_filenm_bk = MS_wave_filenm
    MS_wave_filenm_avg = MS_wave_filenm
    call PostDebugMessage("Switched wave export ON by default in CTEM mode.")
    if (MS_wave_avg_export>0) MS_wave_avg_export = 1
  end if
  if (MSP_detslc>0) then
    MS_wave_export_pzp = MSP_detslc
    j = min(MSP_detslc,MS_stacksize)
    do i=1, MS_stacksize ! determine number of export planes
      if (0/=modulo(i,j)) cycle ! this slice is skipped for output
      MSP_detpln = MSP_detpln + 1
    end do
  else
    MS_wave_export_pzp = MS_stacksize ! default wave readout position
    MSP_detpln = 1
  end if
  write(unit=MSP_stmp,fmt=*) MSP_detpln
  call PostMessage("- number of output planes: "//trim(adjustl(MSP_stmp)))  
  call PostRuntime("probe generation initialized", 1)
! ------------


! ------------
! INITIALIZE MULTISLICE MODULE
  call PostMessage("Initialize module for multislice calculations.")
  nerr = MS_err_num
  call MS_INIT()
  if (MS_err_num/=nerr) then
    call CriticalError("Module init failed.")
  end if
  MS_DEBUG_EXPORT = DEBUG_EXPORT
  MS_VERBO_EXPORT = VERBO_EXPORT
  call PostRuntime("multislice initialized", 1)
! ------------


! ------------
! TRANSFER SUPERCELL DATA TO MODULE AND GENERATE PROPAGATORS
  call PostMessage("Preparing supercell data.")
  call PrepareSupercells()
  call PostRuntime("supercell data prepared", 1)
! ------------


! ------------
! DETECTOR ARRAY SETUP
  if (MSP_ctemmode==0) then
    call PostMessage("Preparing detector data.")
    call MSP_ALLOCDET(nerr)
    if (nerr/=0) then
      call CriticalError("Failed to allocate arrays for detector data.")
    end if
    call MSP_SetAnnularDetectors(nerr)
    if (nerr/=0) then
      call CriticalError("Failed to setup annular segment detectors.")
    end if
    call PostRuntime("detector functions prepared", 1)
  end if
! ------------

! ------------
  if (0==MSP_ctemmode) then ! stem
    call MSP_InitTextOutput(nerr)
    if (nerr/=0) call CriticalError("Failed to initialize text output file.")
  end if
  if (MS_pint_export>0) then ! initialize the probe integration and export
    call InitProbeIntegration()
  end if
  if (MS_wave_avg_export>0) then ! initialize the average wavefunction accumulation and export
    call InitWaveAvg()
  end if
  call PostRuntime("initialization finished", 1)
! ------------
  
! ------------
  call PostSureMessage("Begin of multislice calculations.")
! ------------

! ------------
! SETUP VERTICAL LOOP RANGE // SCAN ROWS
! - default values
  lpy = 0
  mpy = MSP_SF_ndimy-1
! - special sub-frame scan modifications
  if (MSP_ctemmode==0) then ! scan mode setup
    if (prm_py>=0) then ! setup for an individual scan row
      lpy = modulo(prm_py,MSP_SF_ndimy)
      mpy = lpy
    end if
    if (MSP_LastScanPixelY>=0) then ! setup with final scan row
      mpy = max(lpy,modulo(MSP_LastScanPixelY,MSP_SF_ndimy))
    end if
    write(unit=MSP_stmp,fmt='(A,I5,A,I5)') "Vertical loop range: ", lpy, "...",mpy
    call PostMessage(trim(MSP_stmp))
  end if
! ------------
! LOOP VERTICAL THROUGH IMAGE // SCAN ROWS
LV: do
    MSP_ScanPixelY = lpy
! ------------


! ------------
! SETUP HORIZONTAL LOOP RANGE // SCAN COLUMNS
! - default values
  lpx = 0
  mpx = MSP_SF_ndimx-1
! - special sub-frame scan modifications
  if (MSP_ctemmode==0) then ! scan mode setup
    if (prm_px>=0) then ! setup for an individual scan column
      lpx = modulo(prm_px,MSP_SF_ndimx)
      mpx = lpx
    end if
    if (MSP_LastScanPixelX>=0) then ! setup with final scan column
      mpx = max(lpx,modulo(MSP_LastScanPixelX,MSP_SF_ndimx))
    end if
    write(unit=MSP_stmp,fmt='(A,I5,A,I5)') "Horizontal loop range: ", lpx, "...",mpx
    call PostMessage(trim(MSP_stmp))
  end if
! ------------
! START LOOP HORIZONTAL THROUGH IMAGE // SCAN COLUMNS
LH: do
    MSP_ScanPixelX = lpx
! ------------

! ------------
! STANDARD MESSAGE
  if (MSP_ctemmode==0) then ! STEM presets per probe position
  
    write(unit=MSP_stmp,fmt='(A,I4,A,I4)') "Calculating STEM pixel ", &
     &    MSP_ScanPixelX,", ",MSP_ScanPixelY
    if (DEBUG_EXPORT==0 .and. VERBO_EXPORT==0) then ! almost silent mode: switch std output CR/LF control
      OPEN (UNIT=MSP_stdout,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')
      write(unit=MSP_stdout,fmt='("+",A)') trim(MSP_stmp)//"          "
      OPEN (UNIT=MSP_stdout,FORM='FORMATTED',CARRIAGECONTROL='LIST')
    else
      call PostSureMessage(trim(MSP_stmp))
    end if
    
    ! update wave file name (insert pixel numbers)
    if (MS_wave_export>0 .or. MS_wave_avg_export>0 .or. MS_pint_export>0) then
      MS_wave_filenm = MS_wave_filenm_bk
      i = index(trim(MS_wave_filenm),".",BACK=.TRUE.)
      if (i>0) then
        j = len_trim(MS_wave_filenm)
        write(unit=MSP_stmp, fmt='(A,I<MSP_nn1d>.<MSP_nn1d>,A,I<MSP_nn1d>.<MSP_nn1d>,A)') &
     &    MS_wave_filenm(1:i-1)//"_px",MSP_ScanPixelX, "_py", MSP_ScanPixelY, MS_wave_filenm(i:j)
      else 
        write(unit=MSP_stmp, fmt='(A,I<MSP_nn1d>.<MSP_nn1d>,A,I<MSP_nn1d>.<MSP_nn1d>)') &
     &    trim(MS_wave_filenm)//"_px",MSP_ScanPixelX, "_py", MSP_ScanPixelY
      end if
      MS_wave_filenm = trim(MSP_stmp)
      MS_wave_filenm_avg = MS_wave_filenm
    end if
    
    ! ------------
    ! RESET PER PROBE POSITION AVERAGING DATA
    if (MS_pint_export>0) then ! reset the probe integration and export
      call ResetProbeIntegration()
    end if
    if (MS_wave_avg_export>0) then ! reset the average wavefunction accumulation and export
      call ResetWaveAvg()
    end if
  
  end if
! ------------

! ------------
! GENERATE BACKUP WAVE (THIS CREATES OR INSERTS THE INITIAL WAVE FUNCTION)
  if (MSP_use_extinwave==1) then
    call PostMessage("Preparing inserted wave function.")
    call InsertExternalWavefunction()
  else
    call PostMessage("Preparing incoming wave function.")
    call PrepareWavefunction()
  end if
! ------------

! <<<< CORE OF CALCULATION
!
! -------------------------------------------------------------------
! DO MULTISLICE
  call PostMessage("Starting multislice calculations.")
  if (MSP_ctemmode==0) then
    call STEMMultislice()
  else
    call CTEMMultislice()
  end if
! -------------------------------------------------------------------
!
! <<<< CORE OF CALCULATION

!! ------------
!! SAVE THE RESULT TO FILE *** MOVED TO CTEMMultiSlice and STEMMultislice
!!  call CheckLicense()
!  call PostMessage("Saving result to output file.")
!  call SaveResult(MSP_outfile)
!! ------------

  if ((DEBUG_EXPORT>0 .or. VERBO_EXPORT>0) .and. MSP_runtimes==1 ) then ! per loop run time info
    call PostRuntime("multislice calculation done", 0)
  end if


! ------------
! SAVE ACCUMULATED PROBE INTENSITIES
  if (MS_pint_export>=1) then
    call ExportProbeIntensity(trim(MS_wave_filenm_avg))
  end if
  
! ------------
! SAVE ACCUMULATED AVERAGE WAVEFUNCTIONS
  if (MS_wave_avg_export>=1) then
    call ExportWaveAvg(trim(MS_wave_filenm_avg))
  end if

! ------------
! END LOOP HORIZONTAL THROUGH IMAGE
    lpx = lpx + 1
    if (lpx > mpx) exit
    if (MSP_ctemmode==1) exit
  end do LH
! ------------

! ------------
! END LOOP VERTICAL THROUGH IMAGE
    lpy = lpy + 1
    if (lpy > mpy) exit
    if (MSP_ctemmode==1) exit
  end do LV
  call PostSureMessage("End of multislice calculations.")
! ------------

! ------------
! FINISH TEXT OUTPUT
  if (0==MSP_ctemmode) then ! stem
    call MSP_FinishTextOutput(nerr)
    if (nerr/=0) call CriticalError("Failed to finalize text output.")
  end if
! ------------

! ------------
! UNINIT AVERAGE WAVEFUNCTION EXPORT
  if (MS_wave_avg_export>0) then ! Uninitialize the average wavefunction accumulation and export
    call UnInitWaveAvg()
  end if
! UNINIT PROBE INTENSITY INTEGRATION EXPORT
  if (MS_pint_export>0) then
    call UnInitProbeIntegration()
  end if
! ------------

! ------------
! UNINIT MULTISLICE MODULE
  call PostMessage("Uninitialize module for multislice calculations.")
  nerr = MS_err_num
  call MS_UNINIT()
  if (MS_err_num/=nerr) then
    call CriticalError("Module uninit failed.")
  end if
! ------------

! ------------
! APPLY PARTIAL SPATIAL COHERENCE DUE TO SOURCE SHAPE
199  if ((MSP_ApplyPSC/=0).and.(MSP_PC_spatial/=0)) then
    call PostSureMessage("Applying spatial convolution "//&
     & "to consider finite source size.")
    call ApplySpatialCoherence()
  end if
  if ((MSP_ApplyPSC/=0).and.(MSP_PC_spatial==0)) then
    call PostWarning("Inconsistent parameter flags "//&
     & "concerning the consideration of a finite source size.")
    call PostMessage("Set the spatial coherence flag in the "// &
     & "parameter file to 1.")
    call PostMessage("Or do not specify an input file as "// &
     & "command line argument.")
    call PostMessage("Nothing done.")
  end if
! ------------

! ------------
  call PostRuntime("calculation finished", 1)
! ------------  

! ------------
! UNINIT WAVE GENERATION MODULE
  call PostMessage("Uninitialize module for wave modelling.")
  nerr = STF_err_num
  call STF_UNINIT()
  if (STF_err_num/=nerr) then
    call CriticalError("Module uninit failed.")
  end if
! ------------
! UNINIT OTHER MODULES
  call EMS_UNINIT()
  call MSP_UNINIT()
  call Outroduce()
! ------------

! ------------
  call PostRuntime("terminating", 1)
! ------------

end program msa
!****************************************************************************
