!**********************************************************************
!
!  PROGRAM: msa, version 1.1.4
!  FILE: msa.f90
!  PURPOSE:  Entry point for the console application MSA.
!            Multislice calculation for electron diffraction
!
!**********************************************************************
!                                                                      
!   Date: 2023-08-11
!                                                                      
!   Author: Juri Barthel                                               
!           Ernst Ruska-Centre                                         
!           Forschungszentrum J�lich GmbH, 52425 J�lich, Germany       
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
  use Plasmon
  use precision
    
  implicit none
  
! ------------
! declare constants and variables
  integer*4 :: i,j
  integer*4 :: nerr, prm_px, prm_py, lpx, lpy, mpx, mpy
  
  external :: InitRand
  
  call InitRand()
  
!  call CheckLicense()

! ------------
! INITIALIZE MSP MODULE
  call MSP_INIT()
  call EMS_INIT()
  MSP_callApp =                   "[msa] MultiSlice Algorithm"
  MSP_verApp  =                   "1.1.4 64-bit  -  2023 November 13  -"
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
  call PostSureMessage("")
  call PostMessage( " +---------------------------------------------------+" )
  if (MSP_ctemmode==1) then
    call PostMessage( " | Running in plane wave CTEM mode.                  |" )
  else
    call PostMessage( " | Running in convergent probe STEM mode.            |" )
  end if
  if (fpp==4) then
    call PostMessage( " | This is a single precision compile.               |" )
  else
    call PostMessage( " | This is a double precision compile.               |" )
  end if
  call PostMessage( " +---------------------------------------------------+" )
  call PostMessage( "" )
! ------------


! ------------
! VORTEX PROBE (added 2017-10-17 JB)
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
! K-MOMENT ANALYSIS (added 2019-01-11 JB)
  if (MSP_Kmomout > 0) then
    write(unit=MSP_stmp,fmt=*) MSP_KmomMmax
    write(unit=MSP_stmp2,fmt='(F10.1)') MSP_KmomRange
    call PostMessage("Extracting diffraction space integral moments within "// &
         & trim(adjustl(MSP_stmp2))//" mrad up to order "//&
         & trim(adjustl(MSP_stmp))//".")
  end if
! ------------

  
! ------------
! LOAD PARAMETER FILE AND SET ALL RELEVANT PARAMETERS IN MODULES
  call STF_INIT() ! required for loading aberrations
! load parameters
  call PostSureMessage("Loading parameter file.")
  call LoadParameters(MSP_prmfile)
  call PostRuntime("main parameter input finished", 1)
  !MS_fft_flags = MSP_FFTW_FLAG
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
  if ((1==MSP_pimgmode .or. 1==MSP_pdifmode .or. 1==MSP_padifmode) .and. MSP_ctemmode/=0) then
    call PostWarning("Probe intensity integration is not available in CTEM mode.")
    if (1==MSP_pimgmode) call PostSureMessage("- Option /pimg is ignored.")
    if (1==MSP_pdifmode) call PostSureMessage("- Option /pdif is ignored.")
    if (1==MSP_padifmode) call PostSureMessage("- Option /padif is ignored.")
    MSP_pimgmode = 0
    MSP_pdifmode = 0
    MSP_padifmode = 0
  end if
  if (1==MSP_pimgmode .or. 1==MSP_pdifmode .or. 1==MSP_padifmode) then
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
  !FFT_BOUND = max(FFT_BOUND_MIN, 2**CEILING( LOG( real( max(MS_dimx,MS_dimy) ) )/LOG(2.0) )) ! next 2^N above max(nx,ny) or 128
  !if (FFT_BOUND>FFT_BOUND_MAX) call CriticalError("FFT plan exceeds maximum size (8192).")
  !FFT_NYQ = FFT_BOUND/2
  !STF_FFT_BOUND = FFT_BOUND
  !STF_FFT_NYQ = FFT_NYQ 
  write(unit=MSP_stmp,fmt='(A,I5,A,I5)') "- wave function size: ", MS_dimx, " x ", MS_dimy
  call PostMessage( trim(MSP_stmp) )
  !write(unit=MSP_stmp,fmt='(A,I5)') "- FFT working field size: ", FFT_BOUND
  !call PostMessage( trim(MSP_stmp) )
  ! allocate memory for phasegratings, size: MSP_dimcellx * MSP_dimcelly
  call PostMessage("Allocating memory for slice data.")
  call MSP_ALLOCPGR(MSP_dimcellx,MSP_dimcelly,nerr)
  if (MS_err_num/=nerr) then
    call CriticalError("Memory allocation failed.")
  end if
  call PostRuntime("parameter input finished", 1)
! ------------
  
! ------------
! HANDLE COMMAND-LINE PARAMETER OVERRIDES
  MS_nslid = MSP_nslid
  if (MSP_use_extalpha/=0) then ! update the probe convergence using the external command line parameter
    STF_caperture(1) = MSP_extalpha
    write(unit=MSP_stmp,fmt=*) MSP_extalpha
    call PostMessage("- overriding probe convergence by command-line parameter: "// &
       & trim(adjustl(MSP_stmp))//" mrad.")
  end if
  if (MSP_use_extdefocus/=0) then ! update the defocus using the external command line parameter
    call STF_SetAberration(2,1._fpp * MSP_extdefocus,0.0_fpp)
    write(unit=MSP_stmp,fmt=*) MSP_extdefocus
    call PostMessage("- overriding defocus by command-line parameter: "// &
       & trim(adjustl(MSP_stmp))//" nm.")
  end if
  if (MSP_use_extbsh/=0) then ! update the beam shift using the external command line parameter
    call STF_SetAberration(1,1._fpp * MSP_ext_bsx,1.0_fpp * MSP_ext_bsy)
    write(unit=MSP_stmp,fmt=*) MSP_ext_bsx
    write(unit=MSP_stmp2,fmt=*) MSP_ext_bsy
    call PostMessage("- overriding probe offset by command-line parameter: bsx="// &
       & trim(adjustl(MSP_stmp))//", bsy="//trim(adjustl(MSP_stmp2))//" nm.")
  end if
  if (MSP_use_extot/=0) then ! update object tilt parameters from command line input
    MS_objtiltx = real(MSP_OTExX, kind=fpp)
    MS_objtilty = real(MSP_OTExY, kind=fpp)
    write(unit=MSP_stmp,fmt=*) MSP_OTExX
    write(unit=MSP_stmp2,fmt=*) MSP_OTExY
    call PostMessage("- overriding object tilt by command-line parameter: otx="// &
       & trim(adjustl(MSP_stmp))//", oty="//trim(adjustl(MSP_stmp2))//" deg.")
  end if
  if (MSP_FL_varcalc_ex>0) then ! just inform on external parameter use, the new value was already set when loading the parameter file
    write(unit=MSP_stmp,fmt=*) MSP_FL_varcalc_ex
    call PostMessage("- overriding number of passes by command-line parameter: npass="// &
       & trim(adjustl(MSP_stmp)))
  end if
  if (MSP_use_SLC_filenames_ex>0) then ! just inform on external parameter use, the new value was already set when loading the parameter file
    call PostMessage("- overriding slice file name by command-line parameter: "// &
       & trim(adjustl(MSP_SLC_filenames_ex)))
  end if
! ------------
  
! ------------
! INITIALIZE WAVE GENERATION MODULE
  if (MSP_ctemmode/=0) then ! by default switch the wave export ON in CTEM mode
    MS_wave_export = 1 ! activate wave function export by default in CTEM mode
    !MS_wave_export_form = 0 ! set to real-space export by default ! commented out 2020-05-18 (JB) this should allow to use /wavft as switch for saving Fourier-space wave functions
    MS_incwave_export = 1 ! also export the trivial incident wave function
    MS_wave_filenm = MSP_outfile
    MS_wave_filenm_bk = MS_wave_filenm
    MS_wave_filenm_avg = MS_wave_filenm
    call PostDebugMessage("Switched wave export ON by default in CTEM mode.")
  end if
  call MSP_SetDetectionPlanes() ! set detection planes from input parameters
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
! PREPARE PLASMON SCATTERING CALCULATIONS
  if (MSP_do_plasm > 0 .and. MSP_do_plasm <= 2) then ! plasmon excitation calculations used
    ! set remaining input parameters required by the plasmon module
    call PostMessage("Initializing plasmon excitation calculation:")
    PL_ek = STF_ht*1000. ! electron kinetic energy [eV]
    PL_wthr = 0.01 ! run always with 1% probability threshold for remaining higher excitations
    PL_tmax = 0.0 ! init max sample thickness
    do i=1, MS_stacksize ! calculate max. sample tickness
      j = MS_slicestack(i)
      PL_tmax = PL_tmax + real(MS_slicethick(j+1),kind=4)
    end do
    if (MSP_do_plasm == 1) then ! bulk plasmon init
      call PL_init(nerr)
    elseif (MSP_do_plasm == 2) then ! low loss transition init
      !! determine critical angle from diffraction limit of the simulation
      !PL_qc = 1./max(MS_samplingx,MS_samplingy)*MS_lamb/3. ! 2/3 Lambda * Qmax = 2/3 Lambda * 1/2 / sampling-rate
      call PL_init2(nerr)
    else ! no supported type of plasmon calculation (safety catch)
      MSP_do_plasm = 0
    end if
    if (nerr/=0 .and. MSP_do_plasm > 0) then ! initialization error
      call PostWarning("Plasmon module returned error: "//trim(PL_msg_err))
      call PL_deinit()
      MSP_do_plasm = 0
      call PostWarning("Failed to setup plasmon module, option -ple or -plp ignored.")
    elseif (nerr==0 .and. MSP_do_plasm > 0) then ! initialization worked
      if (PL_npemax>0) then ! there will be plasmon excitations in the multislice
        ! post some information on used parameters
        write(unit=MSP_stmp,fmt='(I3)') PL_npemax
        call PostMessage("- calculating plasmon excitation up to "// &
            & trim(adjustl(MSP_stmp))//" levels.")
        if (MSP_do_plasm == 1) then
          write(unit=MSP_stmp,fmt='(F8.1)') PL_ep
          call PostMessage("- plasmon energy: "// &
              & trim(adjustl(MSP_stmp))//" eV.")
        end if
        write(unit=MSP_stmp,fmt='(F8.2)') PL_lp
        call PostMessage("- single plasmon excitation mean-free path: "// &
            & trim(adjustl(MSP_stmp))//" nm.")
        write(unit=MSP_stmp,fmt='(F8.1)') PL_tol
        call PostMessage("- t / lambda: " // trim(adjustl(MSP_stmp)))
        write(unit=MSP_stmp,fmt='(F8.4)') PL_qe*1000.
        call PostMessage("- characteristic angle: " // &
            & trim(adjustl(MSP_stmp)) // " mrad.")
        write(unit=MSP_stmp,fmt='(F8.2)') PL_qc*1000.
        call PostMessage("- crictical angle: " // &
            & trim(adjustl(MSP_stmp)) // " mrad.")
        if (MSP_FL_varcalc<100) then
          write(unit=MSP_stmp,fmt=*) MSP_FL_varcalc
          call PostWarning("Small amount of Monte-Carlo scattering runs. "// &
            & "The final result may not be converged ("// &
            & trim(adjustl(MSP_stmp))//"<100).")
        end if
      else ! (PL_npemax<=0) no plasmon excitations
        call PL_deinit()
        MSP_do_plasm = 0
        if (PL_npemax<=0) then
          call PostWarning("No plasmon excitations expected.")
        end if
        call PostWarning("Plasmon calculation rejected, option -ple is ignored.")
      end if
    end if
    
    !
    ! IMPORTANT: transfer the plasmon code activation flag to module MultiSlice
    MS_do_plasm = min(1, MSP_do_plasm)
    !
  end if
! ------------
  
! ------------
! DETECTOR ARRAY SETUP
  if (0==MSP_ctemmode) then
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
    if (MSP_Kmomout>0) then
      call MSP_SetKmomentDetector(nerr)
      if (nerr/=0) then
        MSP_Kmomout = 0
        call PostWarning("Failed to setup k-moment integration, option -kmom is ignored.")
      end if
    end if
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
      call saddsuffix(trim(MS_wave_filenm), "_px", MSP_ScanPixelX, MSP_nn1d, MSP_stmp2)
      call saddsuffix(trim(MSP_stmp2), "_py", MSP_ScanPixelY, MSP_nn1d, MSP_stmp)
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
!
! <<<< CORE OF CALCULATION
!
! -------------------------------------------------------------------
! DO MULTISLICE
  call PostMessage("Starting multislice calculations.")
  if (MSP_ctemmode==0) then
    call STEMMultiSlice()
  else
    call CTEMMultiSlice()
  end if
! -------------------------------------------------------------------
!
! <<<< CORE OF CALCULATION
!
!
  if ((DEBUG_EXPORT>0 .or. VERBO_EXPORT>0) .and. MSP_runtimes==1 ) then ! per loop run time info
    call PostRuntime("multislice calculation done", 0)
  end if
  
! ------------
! PLASMON STATISTICS
  if (MSP_do_plasm/=0) then
    call PostMessage("Plasmon excitation statistics:")
    j = sum(PL_mc_exc_pop)
    do i=0, PL_npemax
      write(unit=MSP_stmp,fmt='(A,I3,A,F5.3)') "- ",i," x Ep: ", real(PL_mc_exc_pop(i))/real(j)
      call PostMessage(trim(MSP_stmp))
    end do
  end if


! ------------
! SAVE ACCUMULATED PROBE INTENSITIES
  if (MS_pint_export>=1) then
    call ExportProbeIntensity(trim(MS_wave_filenm_avg))
  end if
  
! ------------
! SAVE ACCUMULATED AVERAGE WAVEFUNCTIONS
  if (MS_wave_avg_export > 0) then
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
! OUTPUT POSITION AVERAGED DIFFRACTION PATTERN
  if (MSP_padifmode>0) then
    call ExportProbeAvgDif(trim(MS_wave_filenm_bk))
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
  call PL_deinit()
  call EMS_UNINIT()
  call MSP_UNINIT()
! ------------

! ------------
  call PostRuntime("terminating", 1)
! ------------

end program msa
!****************************************************************************
