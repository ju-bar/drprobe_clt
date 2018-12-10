!****************************************************************************
!
!  wavimg.f90 
!
!  J. Barthel, FZ-Jülich, 2008 - 2018
!
!****************************************************************************
!
!  PROGRAM: wavimg
!
!  PURPOSE: Entry point for the console application.
!           CTEM image calculation from wavefunction
!
!  LINK:    wavimgsubs.f90
!           wavimgprm.f90
!           AberrationFunctions.f90
!           atonum.f
!           atosym.f
!           binio2.f90
!           CellSlicer.f90
!           FFTs.f
!           fscatt.f
!           random.f90
!           TCCformalism.f90
!           csprog.fi
!           ConsoleProgressBar.f90
!
!  DOCUMENTATION: see "wavimg howto.txt"
!
!****************************************************************************

program wavimg

  use IFPORT
  use wavimgprm
  use AberrationFunctions

  implicit none
   
  integer*4 :: mapndx, mapndy, mi, mj, mi1, mi2, mj1, mj2, maxn
  integer*4 :: navg, nwritemap
  complex*8, dimension(:,:), allocatable :: wave, waveft, wave2, wave2ft
  complex*8, dimension(:,:), allocatable :: wave_bk, waveft_bk
  real*4, dimension(:,:), allocatable :: rimg,rimg2,rmap,ravg ! ,rwavre,rwavim
  real*4 :: rval, sampuse
  character(len=5) :: snft
  logical :: fex
  external :: InitRand
  external :: getimageframe, getimageframe2, getimageframe3
  logical, external :: exitprmloop
  
  ! initialization
!  call CheckLicense()
  nerr = 0
  call wavimgprm_init(nerr)
  if (nerr/=0) call CriticalError("Failed to initialize global variables.")
  call InitRand()
  call AF_INIT()
  if (nerr/=0) call CriticalError("Program initialization failed.")
  nerr = 0
  nume = 0
  numw = 0
  nsil = 0
  ndbg = 0
  nnli = 0
  btx = 0.0
  bty = 0.0
  oapr_ex = -1.0
  simgfile_ex = ""
  
  ! *** 
  ! input
  call ParseCommandLine()
  if (nerr/=0) call CriticalError("Failed to determine input parameters.")
  call Introduce()
  call wavimgprm_initcl()
  
  ! read parameter file
  call loadprm()
  
  ! init large ft fields
  maxn = max( max(nwx,nwy), max(nx,ny) )
  nft = max(nft_min, 2**ceiling( log( real( maxn ) )/log(2.0) )) ! next 2^N above max(nx,ny) or 128
  if (nft>nft_max) call CriticalError("FFT plan exceeds maximum size (8192).")
  write(unit=snft,fmt='(I5)') nft
  call PostMessage("- internal array size: "//snft)
  call wavimgprm_init2(nerr)
  if (nerr/=0) call CriticalError("Failed to initialize global calculation arrays.")
  
  ! check parameters
  call checkprm()
  call PostMessage("Parameters checked successfully.")
  
  ! load mtf
  if (domtf) then
    call loadmtf()
  end if
  call PostRuntime("finished reading input parameters")
  
  ! ***
  ! initialize calculations
  ! - prepare parameter loops
  if (nloop>=1) then
    nloopsteps = 1
    do mi=1, nloop
      nloopsteps = nloopsteps * lpsz(mi)
    end do
    write(unit=stmp,fmt='(A,I5,A)') "Calculating ", nloopsteps, &
     &   " images with parameter variation loops."
    call PostMessage(trim(stmp))
    if (nnli==1) then
      call PostMessage("Option /nli: Single image output suppressed during parameter loop.")
    end if
  end if
  ! - initialize parameter variation map
  mapndx = nx
  mapndy = ny
  if (notype==6) then ! -> map is generated
    ! the map parameters are defined by the innermost loops (index 1 and index 2)
    if (nloop>=1) then
      call PostMessage("Determining map parameters:")
      mapndx = int( real(nx)/real(lpsz(1)) )
      if (nloop>=2) then
        mapndy = int( real(ny)/real(lpsz(2)) )
      end if
    
      write(unit=stmp,fmt='(A,I5," x ",I5)') "  Total Map image size: ", nx, ny
      call PostMessage(trim(stmp))
      write(unit=stmp,fmt='(A,I5," x ",I5)') "  Map image patch size: ", mapndx, mapndy
      call PostMessage(trim(stmp))
      write(unit=stmp,fmt='(A,I5," x ",I5)') "  Map image number of patches:", lpsz(1), lpsz(2)
      call PostMessage(trim(stmp))
    else
      ! no loop is defined, it makes no sense to calculate a map image
      ! therefore, the map image creation is switched off here and set to normal image creation
      notype = 0
      ! this should already happen in checkprm, but is done for security reasons again here
    end if
  end if
  
  ! allocate the wave function and image memory
  allocate(wave(nwx,nwy),waveft(nwy,nwx),rimg(nx,ny),ravg(nx,ny),rmap(nx,ny),stat=nerr)
  if (nerr/=0) call CriticalError("Failed to allocate memory.")
  allocate(wave_bk(nwx,nwy),waveft_bk(nwy,nwx),stat=nerr)
  if (nerr/=0) call CriticalError("Failed to allocate memory.")
  
  ravg = 0.0
  navg = 0
  nwritemap = 0
  
  ! setup output image parameters
  if (dofrm==0) then
    sampuse = (swx+swy)*0.5
    rval = abs(0.5*(swx-swy)/sampuse)
    if (rval>0.001 .and. domtf/=0) then
      write(unit=stmp,fmt=*) rval*100.0
      call PostWarning("The target sampling rate is anisotropic by "// &
       & trim(adjustl(stmp))//"%. This will cause a systematic error"// &
       & " for the dampening by the MTF.")
    end if
  else
    sampuse = simg
  end if
  
  ! reset the loop indices to start
  lpidx(:) = 1
  iloopstep = 1
  
  ! backup current aberration status
  call AF_BackupAllAberrations()
  
  ! prepare static image calculation arrays
  call PrepareIC(1)
  
  call PostRuntime("finished initialization")
  
  !----------------------------------------------------!
  !
  ! ******> parameter loops start here <******
  !
  prmloop: do
  
    ! reset the wave file names ( secure in case of loop over wave functions )
    siwav = trim(swavfile)
    siout = trim(simgfile)
    ! backup image file name
    sioutbk = siout
    
    ! update the loop variables depending on the current loop indices lpidx
    call setloopvariable()
  
    ! check wave file existence
    inquire(file=trim(siwav),exist=fex)
    if (.not.fex) then
      call PostWarning("Wave function file ["//trim(siwav)//"] not found.")
      goto 80 ! go for the next loop index or exit the loop, the decision is made at label 80
    end if
  
    ! load the wavefunction
    call PostMessage("Loading wave function data from file ["//trim(siwav)//"].")
    call loadwave(trim(siwav),nwx,nwy,wave,nerr)
    if (nerr/=0) call CriticalError("Failed load the wave function from file ["//trim(siwav)//"].")
    call PostDBGMessage("Transforming wave to Fourier space.")
    call FT(nwx,nwy,wave,waveft,"for")
  
    ! save backup waves
    wave_bk = wave
    waveft_bk = waveft
    
    if (ndbg>0) call PostRuntime("(dbg) wave function loaded")
  
  
    ! set new internal output image name
    siout = sioutbk
    call addloopimagename(siout)
    call PostDBGMessage("Current output file name is ["//trim(siout)//"].")
  
    ! ******>
    ! now we apply imaging properties, thus select between image and wave output
    if (notype==0.or.notype==6) then ! do image calculations ...
      call PostMessage("Creating TEM image from wave.")
      ! select the coherence model
      select case (ncohm)
      case (1) ! partial temporal coherence integerated, spatial coherence quasi by wave envelopes, no cross-terms
        call ExplicitTCC2(wave)
      case (2) ! partial temporal and spatial coherence with full cross-terms !!! lengthy calculations
        call ExplicitTCC(wave)
      case (3) ! partial temporal and spatial coherence with envelopes, ONLY LINEAR TERMS!
        call LinearEnvelopeTCC(wave)
      case (4) ! partial temporal and spatial coherence with envelopes, non-linear, no cross-terms, very lengthy calculations!
        call EnvelopeTCC(wave)
      case (5) ! frozen lattice coherent sub-image averaging !!! focal and tilt variation scheme is combined & random !!!
        call CoherentSubImage(wave)
      case (6) ! quasi-coherent image calculation
        call QuasiCoherent(wave)
      end select ! case (ncohm)
      if (ndbg>0) call PostRuntime("(dbg) image calculated on wave function frame")
      
      ! extract image intensity from selected area
      call PostMessage("Extracting image frame data.")
      call getimageframe(nwx,nwy,wave,nx,ny,rimg)
      if (ndbg>0) call PostRuntime("(dbg) image extracted to image frame")
    
      ! Apply image spread and MTF, optional: convert to output type + noise
      if (.not.(notype==0.and.ncohm==5)) then
        ! Do not apply envelopes in case of coherent sub-image averaging!
        ! This will be done later, when the final image is calculated.
        call PostMessage("Applying envelopes to image data.")
        call ApplyImgEnvelopes(rimg,nx,ny) 
        if (nerr/=0) call CriticalError("Failed to apply envelopes to image data.")
        if (ndbg/=0) then ! speed issue (full array operation)
          rval = sum(rimg)/real(nx*ny)
          write(unit=stmp,fmt='(G12.4)') rval
          call PostDBGMessage("- mean value of current image: "//trim(adjustl(stmp)))
          call PostRuntime("(dbg) image envelopes applied")
        end if ! dbg
      end if ! envelopes applied
    
      ! Image output
      if (nloop>0.and.nnli==1) then
        ! call PostMessage("Option /nli: Single image output suppressed during parameter loop.")
      else
        call PostMessage("Saving image data to file ["//trim(siout)//"].")
        if (notype==0.or.notype==6) then ! image calculation
          if (doint/=0) then ! integer value image output?
            select case (doint)
            case(1) ! 32-bit integer
              call PostMessage("Data format: 32-bit integer")
              call SaveIntegerImage32(trim(siout),nx,ny,rimg)
            case(2) ! 16-bit integer
              call PostMessage("Data format: 16-bit integer")
              call SaveIntegerImage16(trim(siout),nx,ny,rimg)
            end select ! case (doint)
          else ! float 32-bit output
            call PostMessage("Data format: 32-bit float")
            rimg = abs(rimg)
            call savedata(trim(siout),nx*ny,rimg,nerr)
          end if
          if (nerr/=0) call CriticalError("Failed to save image data to file.")
        end if
      end if
    
      ! accumulate map patches to the map
      if (notype==6) then                     ! add component to map image here!
        call PostMessage("Adding image patch to map.")
        mi1 = 1+(lpidx(1)-1)*mapndx           ! get x-start index of patch in image
        mi2 = lpidx(1)*mapndx                 ! get x-stop index of patch in image
        if (lpidx(1)>=lpsz(1)) mi2 = nx       ! for last patch, make shure that index does not run out of the image
        mj1 = 1+(lpidx(2)-1)*mapndy           ! get y-start index of patch in image
        mj2 = lpidx(2)*mapndy                 ! get y-stop index of patch in image
        if (lpidx(2)>=lpsz(2)) mj2 = ny       ! for last patch, make shure that index does not run out of the image
        write(unit=sdbgmsg,fmt='(A,I5,A,I5,A,I5,A,I5,A)') &
     &        " Adding patch frame to map [(",mi1,",",mi2,"),(",mj1,",",mj2,")]."
        call PostDBGMessage(trim(sdbgmsg))
        do mj = mj1, mj2
          do mi = mi1, mi2
            rmap(mi,mj) = rimg(mi,mj)
          end do
        end do
      end if
      
      ! accumulate cherent sub-image average
      if (notype==0.and.ncohm==5) then ! coherent sub-image averaging
        ravg = ravg + rimg
        navg = navg + 1
        write(unit=stmp,fmt='(I)') navg
        call PostDBGMessage("Added image #"//trim(adjustl(stmp))//" to average image.")
      end if
      if (ndbg>0) call PostRuntime("(dbg) image output finished")
    
    else ! (notype/=0 or 6)
      ! the output is a wave function
      ! - allocate helper arrays
      allocate(wave2(nx,ny),wave2ft(ny,nx),rimg2(nx,ny),stat=nerr)
      if (nerr/=0) call CriticalError("Failed to allocate memory.")
      
      ! apply coherent aberrations to wave
      call PostMessage("Applying coherent wave aberrations.")
      call AberrateWaveFourier(waveft,nwx,nwy,swx,swy)
      
      ! extract the wavefunction for the specified frame by
      ! using a bi-cubic interpolation
      call FT(nwx,nwy,wave,waveft,"bac")
      call getimageframe3(nwx,nwy,wave,nx,ny,wave2)
      call FT(nx,ny,wave2,wave2ft,"for")
      !
      ! apply quasi-coherent (linear) wave dampening
      ! by temporal and spatial coherence, and MTF
      call PostMessage("Applying quasi-coherent wave dampening.")
      call DampenWaveFourier(wave2ft,nx,ny)
      ! update real-space wave
      call FT(nx,ny,wave2,wave2ft,"bac")
  
      ! export wave data to file
      select case (notype)
        case (1) ! pure wave data
          call PostMessage("Saving wave function to file ["//trim(siout)//"].")
          call savedatac8(trim(siout),nx*ny,wave2,nerr)
          if (nerr/=0) call CriticalError("Failed to save wave function to file.")
        case (2) ! wave amp
          rimg = cabs(wave2)
          call PostMessage("Saving wave amplitude to file ["//trim(siout)//"].")
          call savedata(trim(siout),nx*ny,rimg,nerr)
          if (nerr/=0) call CriticalError("Failed to save wave amplitude to file.")
        case (3) ! wave phase
          call getwavephase(wave2,rimg,nx,ny)
          call PostMessage("Saving wave phase to file ["//trim(siout)//"].")
          call savedata(trim(siout),nx*ny,rimg,nerr)
          if (nerr/=0) call CriticalError("Failed to save wave phase to file.")
        case (4) ! wave re
          rimg = real(wave2)
          call PostMessage("Saving wave real part to file ["//trim(siout)//"].")
          call savedata(trim(siout),nx*ny,rimg,nerr)
          if (nerr/=0) call CriticalError("Failed to save wave real part to file.")
        case (5) ! wave im
          rimg = imag(wave2)
          call PostMessage("Saving wave imaginary part to file ["//trim(siout)//"].")
          call savedata(trim(siout),nx*ny,rimg,nerr)
          if (nerr/=0) call CriticalError("Failed to save wave imaginary part to file.")
      end select ! case (notype)
    
      ! deallocate secondary image memory
      deallocate(wave2,wave2ft,rimg2,stat=nerr)
      if (nerr/=0) call CriticalError("Failed to deallocate memory.")
      if (ndbg>0) call PostRuntime("(dbg) wave function output finished")
  
    end if ! (notype==0)
    
  
    !-> MAP IMAGE EXPORT 
    !
    ! check if it is time to save the map image
    ! this is done if either the single set loop or the two innermost loops are both at max. index
    nwritemap = 0
    if (notype==6) then
      if (nloop==1) then ! only the first loop is used to create a map
        if (lpidx(1)==lpsz(1)) then ! current image is last of the loop
          nwritemap = 1
        end if
      else if (nloop>=2) then ! the first and the second loop are used to create a map
        if (lpidx(1)==lpsz(1) .and. lpidx(2)==lpsz(2)) then ! current image is last of the loop
          nwritemap = 1
        end if
      end if
    end if
    ! no try to go for the map export
    if (notype==6.and.nwritemap==1) then ! save map image
      ! get output file backup name
      siout = sioutbk
      ! modify output file name for map export
      call addloopmapname(siout)
      ! repport the name to console
      call PostMessage("Saving map data to file ["//trim(siout)//"].")
      ! check the output format
      if (doint/=0) then ! convert to selected integer format
        select case (doint)
        case(1)
          call PostMessage("Data format: 32-bit integer")
          call SaveIntegerImage32(trim(siout),nx,ny,rmap)
        case(2)
          call PostMessage("Data format: 16-bit integer")
          call SaveIntegerImage16(trim(siout),nx,ny,rmap)
        end select ! case (doint)
      else ! save float data
        call PostMessage("Data format: 32-bit float")
        call savedata(trim(siout),nx*ny,rmap,nerr)
      end if
      if (nerr/=0) call CriticalError("Failed to save map image data to file.")
      if (ndbg>0) call PostRuntime("(dbg) map image output finished")
    end if
    !
    !<- MAP IMAGE EXPORT
  
  
    !-> Coherent sub-image average export
    !
    ! * This type of image requires the first loop to be a loop over wavefunctions
    !   Other setups lead to crap results.
    !   But be warned! No sanity check is done here.
    ! * The export will be initiated when loop 1 is at max. index
    ! * The image envelopes and noise is applied here also
    !   This sequence (averaging first, envelopes second) is mandatory in case of
    !   the noise and saves computation time concerning the envelope application.
    !  
    if (notype==0.and.ncohm==5.and.lpidx(1)>=lpsz(1)) then ! save averaged image
      ! renormalize to the average
      write(unit=sdbgmsg,fmt=*) "Normalizing the averaged image, dividing by ",real(navg)
      call PostDBGMessage(trim(sdbgmsg))
      ravg = ravg / real(navg)
      !
      ! the average image is ready now, apply envelopes to the average images
      call PostMessage("Applying envelopes to image data.")
      call ApplyImgEnvelopes(ravg,nx,ny) 
      if (nerr/=0) call CriticalError("Failed to apply envelopes to image data.")
      !
      if (ndbg/=0) then ! speed issue (full array operation)
        rval = sum(ravg)/real(nx*ny)
        write(unit=stmp,fmt='(G12.4)') rval
        call PostDBGMessage("- mean value of current image: "//trim(adjustl(stmp)))
      end if
      !
      ! get output file backup name
      siout = sioutbk
      ! modify output file name for map export
      call addloopavgname(siout)
      ! report the name to console
      call PostMessage("Saving data to file ["//trim(siout)//"].")
      ! check the output format
      if (doint/=0) then
        select case (doint)
        case(1)
          call PostMessage("Data format: 32-bit integer")
          call SaveIntegerImage32(trim(siout),nx,ny,ravg)
        case(2)
          call PostMessage("Data format: 16-bit integer")
          call SaveIntegerImage16(trim(siout),nx,ny,ravg)
        end select ! case (doint)
      else
        call PostMessage("Data format: 32-bit float")
        call savedata(trim(siout),nx*ny,ravg,nerr)
      end if
      ! reset the averaging count & image in case multiple images will be created (multi-loop scenario)
      navg = 0
      ravg = 0.0
      ! any error?
      if (nerr/=0) call CriticalError("Failed to save average image data to file.")
      if (ndbg>0) call PostRuntime("(dbg) output of average of coherent sub image finished")
    end if
    !
    !<- Coherent sub-image average export
    
    if (nloop>0) then ! output loop step information
      write(unit=stmpx,fmt='(I)') iloopstep
      write(unit=stmpy,fmt='(I)') nloopsteps
      call PostMessage("Finished loop step "//trim(adjustl(stmpx))// &
     &   " of "//trim(adjustl(stmpy))//".")
    end if
  
    ! loop exit position
80  if (exitprmloop()) exit prmloop
     
    call nextloopindex()
  
  end do prmloop
  ! ******> parameter loop ends here <******
  !
  !----------------------------------------------------!
  
  call PostRuntime("finished all calculations")
  call PostMessage("Finished all calculations.")
  call Outroduce()
  
  ! deallocate the wave memory
  deallocate(wave,waveft,rimg,ravg,rmap,stat=nerr)
  if (nerr/=0) call CriticalError("Failed to deallocate wave function memory.")
  deallocate(wave_bk,waveft_bk,stat=nerr)
  if (nerr/=0) call CriticalError("Failed to deallocate wave function memory.")
  call PrepareIC(0)
  
  ! *** 
  ! quit
  call AF_UNINIT()
  call wavimgprm_uninit(nerr)
  
end program wavimg


































! ---------------------------------------------------------------------------------
!
! REMOVED CODE
!
! JB 15.01.2010, removed wrong image calculation
!
!    ! partial temporal coherence
!    if (doptc>0) then
!    case (2)
!        call PostMessage("Applying partial temporal coherence explicitly.")
!        call ApplyPTCExpl(nwx,nwy,swx,swy,fs,wl,waveft,wave)
!        ! result is now the image intensity in real space in array wave
!    case (1)
!        ! make image in real space
!        call FT(nwx,nwy,wave,waveft,"bac")
!        wave = wave * conjg(wave)
!        ! go to fs
!        call FT(nwx,nwy,wave,waveft,"for")
!        call PostMessage("Applying partial temporal coherence envelope.")
!        call ApplyPTCEnv(nwx,nwy,swx,swy,waveft,fs,wl)
!        ! go back to rs
!        call FT(nwx,nwy,wave,waveft,"bac")
!        ! result is now the image intensity in real space in array wave
!    case (0)
!        ! make image in real space
!        call FT(nwx,nwy,wave,waveft,"bac")
!        wave = wave * conjg(wave)
!        call PostMessage("No partial temporal coherence considered.")
!    end if ! (doptc>0)
!        
!    ! partial spatial coherence
!    if (dopsc) then
!        call PostMessage("Applying partial spatial coherence envelope.")
!        call ApplyPSCEnv(nwx,nwy,wave)
!    end if ! (doptc>0)

!
!
! program code backup 17.02.2012
!
!
!program wavimg
!
!  use AberrationFunctions
!
!  implicit none
!  
!  include "global.fi"
!  
!  integer*4 :: lpidx, lpisz
!  integer*4 :: wlidx, wlcnt
!  integer*4 :: mapndx, mapndy, mi, mj, mi1, mi2, mj1, mj2
!  integer*4 :: navg
!  complex*8, dimension(:,:), allocatable :: wave, waveft, wave2
!  complex*8, dimension(:,:), allocatable :: wave_bk, waveft_bk
!  real*4, dimension(:,:), allocatable :: rimg,rimg2,rwavre,rwavim,rmap,ravg
!  
!  logical :: fex
!  external :: InitRand
!  
!  ! initialization
!!  call CheckLicense()
!  call InitRand()
!  call PostMessage("")
!  call AF_INIT()
!  call Introduce()
!  nerr = 0
!  nume = 0
!  numw = 0
!  nsil = 0
!  ndbg = 0
!  nnli = 0
!  btx = 0.0
!  bty = 0.0
!  
!  ! *** 
!  ! input
!  call PostMessage("Parsing command line.")
!  call ParseCommandLine()
!  if (nerr/=0) call CriticalError("Failed to determine input parameters.")
!  
!  ! read parameter file
!  call loadprm()
!  
!  ! check parameters
!  call checkprm()
!  call PostMessage("Parameters checked successfully.")
!  
!  ! load mtf
!  if (domtf) then
!    call loadmtf()
!  end if
!  
!  ! ***
!  ! calculations
!  mapndx = nx
!  mapndy = ny
!  if (notype==6) then
!    call PostDBGMessage("Determining map loop parameters:")
!    mapndx = int( real(nx)/real(loopsize) )
!    mapndy = int( real(ny)/ real( 1+int(real(wlstop-wlstart)/real(wlstep)) ) )
!    write(unit=sdbgmsg,fmt=*) "  Total Map image size:", nx, ny
!    call PostDBGMessage(trim(sdbgmsg))
!    write(unit=sdbgmsg,fmt=*) "  Map image patch size:", mapndx, mapndy
!    call PostDBGMessage(trim(sdbgmsg))
!    write(unit=sdbgmsg,fmt=*) "  Map image number of patches:", loopsize, 1+int(real(wlstop-wlstart)/real(wlstep))
!    call PostDBGMessage(trim(sdbgmsg))
!  end if
!  
!  ! allocate the wave memory
!  allocate(wave(nwx,nwy),waveft(nwy,nwx),rimg(nx,ny),ravg(nx,ny),rmap(nx,ny),stat=nerr)
!  if (nerr/=0) call CriticalError("Failed to allocate memory.")
!  allocate(wave_bk(nwx,nwy),waveft_bk(nwy,nwx),stat=nerr)
!  if (nerr/=0) call CriticalError("Failed to allocate memory.")
!  
!  ravg = 0.0
!  navg = 0
!  
!  
!  
!  !----------------------------------------------------!
!  !
!  ! ******> parameter loops start here <******
!  !
!  
!  ! reset the wave file names ( secure in case of loop over wave functions )
!  siwav = trim(swavfile)
!  siout = trim(simgfile)
!  ! backup image file name
!  sioutbk = siout
!  ! backup current aberration status
!  call AF_BackupAllAberrations()
!  ! reset the loop indices to start
!  lpidx(:) = 1
!  
!  ! loop begin
!  prmloop: do
!  
!  ! update the loop variables depending on the current loop indices lpidx
!  call setloopvariable()
!  
!  
!  
!  
!  !!!!! TO BE WORKED OVER !!!!!
!  wlcnt = 0
!  ! start loop
!  wavloop: do wlidx = wlstart, wlstop, wlstep
!  wlcnt = wlcnt + 1
!  if (dowloop>0 .or. notype==6) then ! update wave and image file names due to wave file loop
!    write(unit=sdbgmsg,fmt='(A,I4,A,I4)') "Current wave file loop index: ",wlidx," / ",wlstop
!    call PostMessage(trim(sdbgmsg))
!    ! set new internal wave name
!    siwav = trim(swavfile)
!    call addloopimagename(wlidx,wlstop,siwav)
!    call PostDBGMessage("Current wave file name is ["//trim(siwav)//"].")
!    ! set new internal output image name
!    siout = trim(simgfile) ! reset
!    call addloopimagename(wlidx,wlstop,siout)
!    call PostDBGMessage("Current output file name is ["//trim(siout)//"].")
!    sioutbk = siout ! save backup
!  end if
!  ! check wave file existence
!  inquire(file=trim(siwav),exist=fex)
!  if (.not.fex) then
!    call PostWarning("Wave function file ["//trim(siwav)//"] not found.")
!    cycle wavloop
!  end if
!  
!  ! load the wavefunction
!  call PostMessage("Loading wave function data from file ["//trim(siwav)//"].")
!  call loadwave(trim(siwav),nwx,nwy,wave,nerr)
!  if (nerr/=0) call CriticalError("Failed load the wave function from file ["//trim(siwav)//"].")
!  call PostMessage("Transforming wave to Fourier space.")
!  call FT(nwx,nwy,wave,waveft,"for")
!  
!  ! save backups
!  wave_bk = wave
!  waveft_bk = waveft
!  call AF_BackupAllAberrations()
!  
!  
!  
!  !----------------------------------------------------!
!  !
!  ! ******> parameter loop starts here <******
!  lpisz = loopsize
!  ! start parameter loop
!  if (doloop==0.and.notype/=6) lpisz = 1 ! no parameter loop and no map image creation requested, set parameter loop size to 1
!  prmloop: do lpidx = 1, lpisz
!  
!  if (doloop>0.or.notype==6) then ! loop is done, so save backup wave data and calculate current loop parameter
!    wave = wave_bk
!    waveft = waveft_bk
!    call setloopvariable(lpidx)
!    ! set new internal output image name
!    siout = sioutbk ! reset to outer loop backup
!    call addloopimagename(lpidx,lpisz,siout)
!    call PostDBGMessage("Current output file name is ["//trim(siout)//"].")
!    write(unit=sdbgmsg,fmt='(A,I4,A,I4)') "Current parameter loop index: ",lpidx," of ",lpisz
!    call PostMessage(trim(sdbgmsg))
!  end if
!  
!  ! ******>
!  ! now we apply imaging properties, thus select between image and wave output
!  if (notype==0.or.notype==6) then ! do image creation here
!    call PostMessage("Creating TEM image from wave.")
!  
!    select case (ncohm)
!    case (1) ! partial temporal coherence integerated, spatial coherence quasi by wave envelopes, no cross-terms
!      call ExplicitTCC2(wave)
!    case (2) ! partial temporal and spatial coherence with full cross-terms !!! lengthy calculations
!      call ExplicitTCC(wave)
!    case (3) ! partial temporal and spatial coherence with envelopes, ONLY LINEAR TERMS!
!      call LinearEnvelopeTCC(wave)
!    case (4) ! partial temporal and spatial coherence with envelopes, non-linear, no cross-terms, very lengthy calculations!
!      call EnvelopeTCC(wave)
!    case (5) ! frozen lattice coherent sub-image averaging
!      call CoherentSubImage(wave)
!    end select ! case (ncohm)
!      
!    ! extract image intensity from selected area
!    call PostMessage("Extracting image frame data.")
!    call getimageframe(nwx,nwy,wave,nx,ny,rimg)
!      
!    call PostMessage("Applying envelopes to image data.")
!    call ApplyImgEnvelopes(rimg,nx,ny) 
!    if (nerr/=0) call CriticalError("Failed to apply envelopes to image data.")
!      
!    ! *** 
!    ! output
!!    call CheckLicense()
!    if (dowloop==1.and.nnli==1) then
!      call PostMessage("Option /nli: Single image output suppressed while looping over wavefunctions.")
!    else
!      call PostMessage("Saving image data to file ["//trim(siout)//"].")
!      if (notype==0.or.notype==6) then
!        if (doint/=0) then
!          select case (doint)
!          case(1)
!            call PostMessage("Data format: 32-bit integer")
!            call SaveIntegerImage32(trim(siout),nx,ny,rimg)
!          case(2)
!            call PostMessage("Data format: 16-bit integer")
!            call SaveIntegerImage16(trim(siout),nx,ny,rimg)
!          end select ! case (doint)
!        else
!          call PostMessage("Data format: 32-bit float")
!          call savedata(trim(siout),nx*ny,rimg,nerr)
!        end if
!        if (nerr/=0) call CriticalError("Failed to save image data to file.")
!      end if
!    end if
!    
!    if (notype==6) then                 ! add component to map image here!
!      call PostMessage("Adding image patch to map.")
!      mi1 = 1+(lpidx-1)*mapndx           ! get x-start index of patch in image
!      mi2 = lpidx*mapndx                 ! get x-stop index of patch in image
!      if (lpidx>=lpisz) mi2 = nx        ! for last patch, make shure that index does not run out of the image
!      mj1 = 1+(wlcnt-1)*mapndy           ! get y-start index of patch in image
!      mj2 = wlcnt*mapndy                 ! get y-stop index of patch in image
!      if (wlidx+wlstep>wlstop) mj2 = ny ! for last patch, make shure that index does not run out of the image
!      write(unit=sdbgmsg,fmt='(A,I5,A,I5,A,I5,A,I5,A)') &
!     &  " Adding patch frame to map [(",mi1,",",mi2,"),(",mj1,",",mj2,")]."
!      call PostDBGMessage(trim(sdbgmsg))
!      do mj = mj1, mj2
!        do mi = mi1, mi2
!          rmap(mi,mj) = rimg(mi,mj)
!        end do
!      end do
!    end if
!    
!    if (notype==0.and.ncohm==5) then ! coherent sub-image averaging
!      ravg = ravg + rimg
!      navg = navg + 1
!    end if
!    
!  else ! (notype/=0 or 6) do wave output creation here, only wave data, no coherence issues and no mtf issues
!  
!    ! applying coherent aberrations to wave
!    call PostMessage("Applying coherent wave aberrations.")
!    call AberrateWaveFourier(waveft,nwx,nwy,swx,swy)
!    ! update real-space wave
!    call FT(nwx,nwy,wave,waveft,"bac")
!  
!    ! allocate secondary image memory
!    allocate(rwavre(nwx,nwy),rwavim(nwx,nwy),stat=nerr)
!    if (nerr/=0) call CriticalError("Failed to allocate memory.")
!    allocate(wave2(nx,ny),rimg2(nx,ny),stat=nerr)
!    if (nerr/=0) call CriticalError("Failed to allocate memory.")
!    
!    ! transfer wave data
!    select case (notype)
!      case (1) ! pure wave data
!        rwavre = real(wave)
!        rwavim = imag(wave)
!        call PostMessage("Extracting frame data.")
!        call getimageframe2(nwx,nwy,rwavre,nx,ny,rimg)
!        call getimageframe2(nwx,nwy,rwavim,nx,ny,rimg2)
!        wave2 = cmplx(rimg, rimg2)
!        ! output
!!        call CheckLicense()
!        call PostMessage("Saving data to file ["//trim(siout)//"].")
!        call savedatac8(trim(siout),nx*ny,wave2,nerr)
!        if (nerr/=0) call CriticalError("Failed to save data to file.")
!      case (2) ! wave amp
!        rwavre = cabs(wave)
!        call PostMessage("Extracting frame data.")
!        call getimageframe2(nwx,nwy,rwavre,nx,ny,rimg)
!        ! output
!!        call CheckLicense()
!        call PostMessage("Saving data to file ["//trim(siout)//"].")
!        call savedata(trim(siout),nx*ny,rimg,nerr)
!        if (nerr/=0) call CriticalError("Failed to save data to file.")
!      case (3) ! wave phase
!        call getwavephase(wave,rwavre,nwx,nwy)
!        call PostMessage("Extracting frame data.")
!        call getimageframe2(nwx,nwy,rwavre,nx,ny,rimg)
!        ! output
!!        call CheckLicense()
!        call PostMessage("Saving data to file ["//trim(siout)//"].")
!        call savedata(trim(siout),nx*ny,rimg,nerr)
!        if (nerr/=0) call CriticalError("Failed to save data to file.")
!      case (4) ! wave re
!        rwavre = real(wave)
!        call PostMessage("Extracting frame data.")
!        call getimageframe2(nwx,nwy,rwavre,nx,ny,rimg)
!        ! output
!!        call CheckLicense()
!        call PostMessage("Saving data to file ["//trim(siout)//"].")
!        call savedata(trim(siout),nx*ny,rimg,nerr)
!        if (nerr/=0) call CriticalError("Failed to save data to file.")
!      case (5) ! wave im
!        rwavre = imag(wave)
!        call PostMessage("Extracting frame data.")
!        call getimageframe2(nwx,nwy,rwavre,nx,ny,rimg)
!        ! output
!!        call CheckLicense()
!        call PostMessage("Saving data to file ["//trim(siout)//"].")
!        call savedata(trim(siout),nx*ny,rimg,nerr)
!        if (nerr/=0) call CriticalError("Failed to save data to file.")
!    end select ! case (notype)
!    
!    ! deallocate secondary image memory
!    deallocate(rwavre,rwavim,wave2,rimg2,stat=nerr)
!    if (nerr/=0) call CriticalError("Failed to deallocate memory.")
!  
!  end if ! (notype==0)
!  ! <******
!  
!  end do prmloop
!  ! ******> parameter loop ends here <******
!  !
!  !----------------------------------------------------!
!  
!  
!  if (dowloop==0 .and. notype/=6) exit wavloop ! security exit
!  
!  end do wavloop
!  ! ******> wave function loop ends here <*******
!  !
!  !----------------------------------------------------!
!  
!  
!  if (notype==6) then ! save map image
!    call PostMessage("Saving map data to file ["//trim(simgfile)//"].")
!    if (doint/=0) then
!      select case (doint)
!      case(1)
!        call PostMessage("Data format: 32-bit integer")
!        call SaveIntegerImage32(trim(simgfile),nx,ny,rmap)
!      case(2)
!        call PostMessage("Data format: 16-bit integer")
!        call SaveIntegerImage16(trim(simgfile),nx,ny,rmap)
!      end select ! case (doint)
!    else
!      call PostMessage("Data format: 32-bit float")
!      call savedata(trim(simgfile),nx*ny,rmap,nerr)
!    end if
!    if (nerr/=0) call CriticalError("Failed to save map image data to file.")
!  end if
!  
!  if (notype==0.and.ncohm==5) then ! save averaged image
!    ravg = ravg / real(navg)
!    call PostMessage("Saving average image to file ["//trim(simgfile)//"].")
!    if (doint/=0) then
!      select case (doint)
!      case(1)
!        call PostMessage("Data format: 32-bit integer")
!        call SaveIntegerImage32(trim(simgfile),nx,ny,ravg)
!      case(2)
!        call PostMessage("Data format: 16-bit integer")
!        call SaveIntegerImage16(trim(simgfile),nx,ny,ravg)
!      end select ! case (doint)
!    else
!      call PostMessage("Data format: 32-bit float")
!      call savedata(trim(simgfile),nx*ny,ravg,nerr)
!    end if
!    if (nerr/=0) call CriticalError("Failed to save average image data to file.")
!  end if
!  
!  call PostMessage("Finished.")
!  call Outroduce()
!  
!  ! deallocate the wave memory
!  deallocate(wave,waveft,rimg,ravg,rmap,stat=nerr)
!  if (nerr/=0) call CriticalError("Failed to deallocate wave function memory.")
!  deallocate(wave_bk,waveft_bk,stat=nerr)
!  if (nerr/=0) call CriticalError("Failed to deallocate wave function memory.")
!  
!  ! *** 
!  ! quit
!  call AF_UNINIT()
!  
!end program wavimg
