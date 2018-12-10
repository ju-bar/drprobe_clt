! --------------------------------------------------------------------------
!
!  subroutine PrepareIC
!
! Prepares static arrays for coherent image calculation
! This needs to be done just once before calling the other
! routines in this file.
!
! - input: 1 = prepare, 0 = unprepare
!
! - prepare
!   allocates and initializes:
!   * ic_iwx, ic_iwy
!   * ic_twx, ic_twy
!   * ic_objaper
!
! - unprepare
!   deallocates the arrays listed above
!
subroutine PrepareIC(pmode)

  use wavimgprm
  use AberrationFunctions
  
  implicit none
  
  integer*4, intent(in) :: pmode
  integer*4 :: nalloc
  integer*4 :: i, j
  integer*4 :: nx2, ny2
  real*4 :: itowx, itowy
  real*4 :: tx, ty
  real*4 :: awxs, awys, wa
  real*4 :: apr, apx, apy
  
  nalloc = 0
  call PostDBGMessage("Prapring global image calculation arrays")
  
  ! by default always deallocate
  if (allocated(ic_iwx)) deallocate(ic_iwx,stat=nalloc)
  if (allocated(ic_iwy)) deallocate(ic_iwy,stat=nalloc)
  if (allocated(ic_twx)) deallocate(ic_twx,stat=nalloc)
  if (allocated(ic_twy)) deallocate(ic_twy,stat=nalloc)
  if (allocated(ic_objaper)) deallocate(ic_objaper,stat=nalloc)
  
  if (pmode/=0) then ! allocate and initialize
    allocate(ic_iwx(nwx),ic_iwy(nwy), stat=nalloc)
    allocate(ic_twx(nwx),ic_twy(nwy), stat=nalloc)
    allocate(ic_objaper(nwy,nwx), stat=nalloc)
    nx2 = nwx/2
    ny2 = nwy/2
    itowx = wl/(swx*real(nwx))
    itowy = wl/(swy*real(nwy)) 
    tx = 0.001*btx
    ty = 0.001*bty
    apr = 0.001*oapr
    apx = 0.001*oapx
    apy = 0.001*oapy
    ! setup scrambled frequencies x
    do i=1, nwx
      ic_iwx(i) = modulo(i-1+nx2,nwx)-nx2
    end do
    ic_twx = itowx*real(ic_iwx) + tx
    ! setup scrambled frequencies y
    do i=1, nwy
      ic_iwy(i) = modulo(i-1-ny2,nwy)-ny2
    end do
    ic_twy = itowy*real(ic_iwy) + ty
    ! aperture and chi on twx,y
    do i=1, nwx
      awxs = (itowx*real(ic_iwx(i))-apx)**2
      do j=1, nwy
        awys = (itowy*real(ic_iwy(j))-apy)**2
        wa = sqrt(awxs + awys)
        ic_objaper(j,i) = 0.5-0.5*tanh((wa-apr)/apr*100.0)
      end do
    end do
  end if
  
  return
  
end subroutine PrepareIC
! --------------------------------------------------------------------






! --------------------------------------------------------------------------
!
!  subroutine CoherentSubImage
!
! Calculation of a coherent sub-image with random transfer variation
! of defocus and transfer beam tilt
!
! This is a simple coherent image calculation. Different from a normal
! calculation, randomly chosen defocus and beam-tilts are applied in the
! coherent transfer function, which is multiplied to the Fourier-transform.
! By this way you obtain a coherent sub-image for later averaging towards
! a partial coherent image.
!
! !!! call PrepareIC(1) before using this routine
!
! - input: real-space wave (coherent)
! - 1: FFT -> fourier space wave
! - 2: dice random numbers of an additional defocus and beam tilt
! - 2: - a: modify fourier space wave by current transfer function
!      - b: iFFT -> modified real-space wave
!      - c: take absolute square -> modified coherent sub-image      
! - output: real-space image (coherent but modified)
!
subroutine CoherentSubImage(wave)

  use wavimgprm
  use CSPROG
  use AberrationFunctions
  
  implicit none
  
  complex*8, intent(inout) :: wave(nwx,nwy) ! real space wave function (input)
                                            ! real space image (output) (in real-part of wave)
  
  integer*4 :: nx2, ny2                     ! data half sizes
  integer*4 :: i, j                         ! std. iterators
  real*4 :: uFS, uSC                        ! applied focus spread and semi convergence parameters
  real*4 :: dZ, dQx, dQy                    ! iteration physical parameters
  real*4 :: chi                             ! aberration value
  real*4 :: Zbk                             ! backup defocus
  complex*8 :: cval                         ! some complex value
  complex*8, allocatable :: wtmp(:,:)       ! temporary wave data (real-space)
  
  external :: FSFFT
  external :: PostMessage
  real*4, external :: GaussRand, UniRand
  
  ! init
  call PostMessage("Coherent sub-image calculation.")
  
  ! allocate wave arrays
  allocate(wtmp(nft,nft), stat=nerr)
  wtmp = cmplx(0.0,0.0)
  
  ! transfer wave data to temp array
  do j=1, nwy
    do i=1, nwx
      wtmp(i,j) = wave(i,j)
    end do
  end do
  
  ! forward FFT
  call FSFFT(wtmp,nwx,nwy,'for')
  
  ! fourier-space wave is now in wtmp, scrambled and transposed
  
  ! clear up data in array wave, this array will recieve the output
  wave = cmplx(0.0,0.0)
  
  ! calculate parameters
  nx2 = nwx/2
  ny2 = nwy/2
  uFS = 0.0
  uSC = 0.0
  if (doptc>0) uFS = sqrt(0.5)*fs
  if (dopsc>0) uSC = 0.5*sc
  
  ! dice randomized transfer parameters
  ! defocus
  dZ = uFS*GaussRand()
  ! beam tilt (X)
  dQx = uSC*GaussRand()
  ! beam tilt (Y)
  dQy = uSC*GaussRand()
  
  write(unit=sdbgmsg,fmt='(A,F6.2,A,G10.3,A,G10.3,A)') &
     & "- variation of defocus: ",dZ," nm and of beam tilt: (",dQx,",",dQy,") rad."
  call PostDBGMessage(trim(sdbgmsg))
  
  ! backup mean defocus
  Zbk = AF_wa(1,2)
  ! set new defocus
  AF_wa(1,2) = AF_wa(1,2) + dZ
  
  ! apply phaseplate with tilted transfer term and aperture
  do i=1, nwx
    do j=1, nwy
      chi = AF_AberrationFunction(ic_twx(i)+dQx, ic_twy(j)+dQy)
      wtmp(j,i) = wtmp(j,i)*exp( cmplx(0.0,-chi) )*ic_objaper(j,i)
    end do
  end do
  
  ! iFFT to real-space for wsub
  call FSFFT(wtmp,nwx,nwy,'bac')
  
  ! calculate real-space image intensity (take absolute square here)
  call PostDBGMessage("- caluclating image intensity from abs. square of the real-space wave function.")
  
  do j=1, nwy
    do i=1, nwx
      cval = wtmp(i,j)
      wave(i,j) = cval*conjg(cval)
    end do
  end do
  
  ! set back the defocus to the backup value
  AF_wa(1,2) = Zbk
  
  ! finish and deallocate
  deallocate(wtmp, stat=nerr)
  
  return

    end subroutine CoherentSubImage
! --------------------------------------------------------------------











! --------------------------------------------------------------------------
!
!  subroutine ExplicitTCC
!
! explicit TCC with partial temporal and partial spatial coherence
!
! This explicit focal and convergent beam averaging takes
! also the cross-terms between the two sources for partial coherence
! into account.
! The image intensities are calculated in real-space. The formalism
! is therefore fully non-linear!
!
! The distribution functions for the semi-convergence angles and the
! focus spread are of Gaussian shape with 1/e-half width fs and sc.
! Other distribution functions can be implemented!
!
! !!! call PrepareIC(1) before using this routine
!
! - input: real-space wave (coherent)
! - 1: FFT -> fourier space wave
! - 2: 3D-loop - a: modify fourier space wave by current transfer function
!              - b: iFFT -> modified real-space wave
!              - c: take absolute square -> modified coherent sub-image      
!              - d: add to total image intensity with respective weight
! - output: real-space image (partially coherent)
!
! 3D-loop:
!   The size of the 3D loop is set by parameters NKCB and NKFS, the number
!   of kernel samples offside the central beam and focus
!   Set NKCB or NKFS to zero if you want to switch off the respective convolution
!   The 3D-loop is ordered in an outside 2D-loop for the convergence and an
!   inside 1D-loop for the focus spread. By this way a full aberrations phase
!   plate is calculated for every sample of the convergence kernel, and the
!   focus phase shift is applied extra only inside the inner 1D-loop.
!
! TCC and explicit averaging:
!   I(R) = Int_Q1 d2Q1 Int_Q2 d2Q2 Psi(Q1)*Psi^*(Q2) Exp[ i R.(Q1-Q2) ] * ( Int_Eps dEps Int_Theta d2Theta Int_R0 d2R0 f(Eps) s(Theta) h(R0) Exp[ -i Chi(Q1+k0*Theta,C1+Eps) ] Exp[ i Chi(Q2+k0*Theta,C1+Eps) ] Exp[ -i R0.(Q1-Q2) ] )
!   or
!   I'(R) = Int_Eps dEps Int_Theta d2Theta  f(Eps) s(Theta) | Int_Q d2Q Exp[ i Q.R ] * Psi(Q) * Exp[ -i Chi(Q+k0*Theta,C1+Eps)] |^2
!   I(R) = Int_Q d2Q Exp[ -i Q.R ] Eim(Q) ( Int_R d2R I'(R) Exp[ -i Q.R ] )
!
subroutine ExplicitTCC(wave)

  use wavimgprm
  use AberrationFunctions
  use CSPROG
  
  implicit none
  
  complex*8, intent(inout) :: wave(nwx,nwy) ! real space wave function (input)
                                            ! real space image (output) (in real-part of wave)
  
  integer*4 :: nx2, ny2                     ! data half sizes
  integer*4 :: i, j                         ! std. iterators
  integer*4 :: ifoc, iqx, iqy               ! iterators
  real*4 :: dZ, dQx, dQy, dQ2, sc2          ! iteration physical parameters
  real*4 :: itoz                            ! focus sampling (focal kernel)
  real*4 :: itoq                            ! Fourier space sampling (CB kernel)
  real*4 :: maxq2                           ! maximum convergent beam value squared
  real*4 :: fwt, cwt, totwt, sumwt          ! kernel weights
  real*4 :: chi, wx2, w2, wi, wj            ! aberration value
  real*4 :: ffoc
  complex*8 :: cval                         ! some complex value
  complex*8, allocatable :: wtmp(:,:)       ! temporary wave data (real-space)
  complex*8, allocatable :: wsub(:,:)       ! sub-image data
  complex*8, allocatable :: cqpp(:,:)       ! convergent beam phase plate
    
  external :: FSFFT
  external :: PostMessage
  
  ! init
  call PostMessage("Image calculation by explicit TCC integration.")
  ! get number of steps
  i = 2+(1+2*NKCB)*(1+2*NKCB)*(1+2*NKFS)*1
  call CSPROG_INIT(i,60)
  ! calculate parameters
  nx2 = nwx/2
  ny2 = nwy/2
  itoz = 0.0
  itoq = 0.0
  if (doptc>0) itoz = fkw*fs/real(NKFS)
  if (dopsc>0) itoq = cbkw*sc/real(NKCB)
  maxq2 = (cbkw*sc)**2
  sc2 = sc*sc
  sumwt = 0.0
  cwt = 1.0
  fwt = 1.0
  totwt = 1.0
  
  ! allocate wave arrays
  allocate(wtmp(nft,nft), stat=nerr)
  allocate(wsub(nft,nft), stat=nerr)
  allocate(cqpp(nwy,nwx), stat=nerr)
  wtmp = cmplx(0.0,0.0)
  wsub = cmplx(0.0,0.0)
  cqpp = cmplx(0.0,0.0)
  
  ! transfer wave data to temp array
  do j=1, nwy
    do i=1, nwx
      wtmp(i,j) = wave(i,j)
    end do
  end do
  
  ! forward FFT
  call FSFFT(wtmp,nwx,nwy,'for')
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! fourier-space wave is now in wtmp, scrambled and transposed
  
  ! clear up data in array wave
  wave = cmplx(0.0,0.0)
  
  ! pre-apply aperture to wave in FS (wtmp)
  do i=1, nwx
    do j=1, nwy
      wtmp(j,i) = wtmp(j,i) * ic_objaper(j,i)
    end do
  end do
  
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! loop through number of beam convergence steps ! 2D
  do iqy=-NKCB, NKCB
    
    ! get current convergent beam (Y)
    dQy = itoq*real(iqy)
      
    do iqx=-NKCB, NKCB
        
      ! get current convergent beam (X)
      dQx = itoq*real(iqx)
      dQ2 = dQx**2 + dQy**2
        
      ! skip beams outside max. convergence of kernel
      if (maxq2<dQ2) cycle
        
      if (dopsc>0) cwt = exp( -dQ2/sc2 ) ! SemiConvergenceDistribution(dQx,dQy) ! exp( -dQ2/sc2 )
        
      ! calculate a phaseplate for current direction of convergence
      ! from current aberration function
      do i=1, nwx
        wi = ic_twx(i) + dQx
        do j=1, nwy
          wj = ic_twy(j) + dQy
          chi = AF_AberrationFunction(wi,wj)
          cqpp(j,i) = cmplx(0.0,-chi)
        end do
      end do
        
      ! loop through number of focal kernel steps ! 1D
      do ifoc=-NKFS, NKFS
  
        ! get current kernel focus
        dZ = itoz*real(ifoc)
        
        ! get focal weight from focal distribution function
        if (doptc>0) fwt = exp(-(dZ/fs)**2) ! FocalDistribution(dZ) ! exp(-(dZ/fs)**2) 
        
        ! get focal phase constant
        ffoc = pi*dZ/wl
        
        ! apply phaseplate and focal phase shifts
        do i=1, nwx
          wi = ic_twx(i) + dQx
          wx2 = wi*wi
          do j=1, nwy
            wj = ic_twy(j) + dQy
            w2 = wx2 + wj*wj
            wsub(j,i) = wtmp(j,i)*exp( cqpp(j,i)+cmplx(0.0,-ffoc*w2) )
          end do
        end do
        
        ! iFFT to real-space for wsub
        call FSFFT(wsub,nwx,nwy,'bac')
        
        ! get total weight
        totwt = fwt*cwt
        
        ! sum up the coherent sub images (take absolute square here)
        do j=1, nwy
          do i=1, nwx
            cval = wsub(i,j)
            wave(i,j) = wave(i,j) + totwt*cval*conjg(cval)
          end do
        end do
        
        ! sum up the applied total weights
        sumwt = sumwt + totwt
        
        FSTEP_CUR = FSTEP_CUR+1
        call CSPROG_UPDATE()
        
      end do ! loop NKFS  
        
    end do ! loop NKCB 2
    
  end do ! loop NKCB 1
    
  ! divide result by total sum
  wave = wave / sumwt
  
  call CSPROG_FINISH()
  
  deallocate(wtmp, wsub, stat=nerr)
  deallocate(cqpp, stat=nerr)
  
  return

end subroutine ExplicitTCC
! --------------------------------------------------------------------










! --------------------------------------------------------------------------
!
!  subroutine ExplicitTCC2
!
! explicit TCC with partial temporal and partial spatial coherence
!
! In this routine, the temporal coherence is treated by explicit focal
! averaging, while the spatial coherence is treated by the quasi-coherent
! approach of multiplication of the beams by exp(-grad(chi)**2 * (sc/2)**2)
!
! The distribution functions for the semi-convergence angles and the
! focus spread are of Gaussian shape with 1/e-half width fs and sc.
!
! !!! call PrepareIC(1) before using this routine
!
! - input: real-space wave (coherent)
! - 1: FFT -> fourier space wave
! - 2: 1D-loop - a: modify fourier space wave by current transfer function
!              - b: iFFT -> modified real-space wave
!              - c: take absolute square -> modified coherent sub-image      
!              - d: add to total image intensity with respective weight
! - output: real-space image (partially coherent)
!
! 1D-loop:
!   The size of the 1D loop is set by parameters NKFS, the number
!   of kernel samples offside the central focus
!   Set NKFS to zero if you want to switch off the convolution
!   
!
subroutine ExplicitTCC2(wave)

  use wavimgprm
  use AberrationFunctions
  use CSPROG  
  
  implicit none
  
  complex*8, intent(inout) :: wave(nwx,nwy) ! real space wave function (input)
                                            ! real space image (output) (in real-part of wave)
  
  integer*4 :: nx2, ny2                     ! data half sizes
  integer*4 :: i, j                         ! std. iterators
  integer*4 :: ifoc                         ! iterators
  real*4 :: dZ                              ! iteration physical parameters
  real*4 :: itoz                            ! focus sampling (focal kernel)
  real*4 :: fwt, sumwt                      ! kernel weights
  real*4 :: chi, wi, wj, wx2, w2            ! aberration value
  real*4 :: ffoc                            ! prefactor of focal averaging phase
  real*4 :: scpf                            ! prefactor for partial spatial coherence
  real*4 :: dchix, dchiy                    ! gradient values of the aberration function
  real*4 :: dchix0, dchiy0                  ! gradient values of the aberration function at central beam
  real*4 :: Envs                            ! envelope for beams due to spatial coherence
  real*4, allocatable :: wxs(:), wys(:)     ! prepared scrambled frequencies
  complex*8 :: cval                         ! some complex value
  complex*8, allocatable :: wtmp(:,:)       ! temporary wave data (real-space)
  complex*8, allocatable :: wsub(:,:)       ! sub-image data
  complex*8, allocatable :: cqpp(:,:)       ! convergent beam phase plate
  
  external :: FSFFT
  external :: PostMessage
  
  ! init
  call PostMessage("Image calculation by explicit focal integration.")
  ! calculate parameters
  nx2 = nwx/2
  ny2 = nwy/2
  itoz = 0.0
  if (doptc>0) itoz = fkw*fs/real(NKFS)
  scpf = (sc/2.)**2
  sumwt = 0.0
  fwt = 1.0
  Envs = 1.0
  ! get number of steps
  i = 3+(1+2*NKFS)*1
  call CSPROG_INIT(i,60)
  
  ! allocate wave arrays
  allocate(wtmp(nft,nft), stat=nerr)
  allocate(wsub(nft,nft), stat=nerr)
  allocate(cqpp(nwy,nwx), stat=nerr)
  allocate(wxs(nwx), wys(nwy), stat=nerr)
  wtmp = cmplx(0.0,0.0)
  wsub = cmplx(0.0,0.0)
  cqpp = cmplx(0.0,0.0)
  wxs = 0.0
  wys = 0.0
  
  ! transfer wave data to temp array
  do j=1, nwy
    do i=1, nwx
      wtmp(i,j) = wave(i,j)
    end do
  end do
  
  ! forward FFT
  call FSFFT(wtmp,nwx,nwy,'for')
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! fourier-space wave is now in wtmp, scrambled and transposed
  
  ! clear up data in array wave
  wave = cmplx(0.0,0.0)
  
  ! setup scrambled frequencies x
  wxs = ic_twx*ic_twx
  wys = ic_twy*ic_twy
  !
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! pre-calculate the quasi-coherent transfer terms
  ! from current aberration function
  call AF_AberrationFunctionGradient(ic_twx(1), ic_twy(1), dchix0, dchiy0) ! get central beam gradients
  do i=1, nwx
    wi = ic_twx(i)
    do j=1, nwy
      wj = ic_twy(j)
      if (dopsc>0) then ! linear term envelope for psc
        call AF_AberrationFunctionGradient(wi, wj, dchix, dchiy)
        Envs = exp( -scpf*((dchix-dchix0)**2+(dchiy-dchiy0)**2) )
      end if
      chi = AF_AberrationFunction(wi, wj)
      cqpp(j,i) = ic_objaper(j,i) * exp( cmplx(0.0,-chi) ) * Envs
    end do
  end do
  
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! loop through number of focal kernel steps ! 1D
  do ifoc=-NKFS, NKFS
  
    ! get current kernel focus
    dZ = itoz*real(ifoc)
        
    ! get focal weight from focal distribution function
    if (doptc>0) fwt = exp(-(dZ/fs)**2) ! FocalDistribution(dZ) ! exp(-(dZ/fs)**2) 
        
    ! get focal phase constant
    ffoc = pi*dZ/wl
    
    ! apply phaseplate, spatial coherence damping and focal phase shifts
    do i=1, nwx
      wx2 = wxs(i)
      do j=1, nwy
        w2 = wx2 + wys(j)
        wsub(j,i) = wtmp(j,i)*cqpp(j,i)*exp( cmplx(0.0,-ffoc*w2) )
      end do
    end do
        
    ! iFFT to real-space for wsub
    call FSFFT(wsub,nwx,nwy,'bac')
        
    ! sum up the coherent sub images (take absolute square here)
    do j=1, nwy
      do i=1, nwx
        cval = wsub(i,j)
        wave(i,j) = wave(i,j) + fwt*cval*conjg(cval)
      end do
    end do
        
    ! sum up the applied total weights
    sumwt = sumwt + fwt
        
    FSTEP_CUR = FSTEP_CUR+1
    call CSPROG_UPDATE()
        
  end do ! loop NKFS  
  
  ! divide result by total sum
  wave = wave / sumwt
  
  call CSPROG_FINISH()
  
  deallocate(wtmp, wsub, stat=nerr)
  deallocate(cqpp, stat=nerr)
  deallocate(wxs, wys, stat=nerr)
  
  return

end subroutine ExplicitTCC2
! --------------------------------------------------------------------











! --------------------------------------------------------------------
!
!  subroutine EnvelopeTCC
!
! Envelope-based image calculation by transmission cross-coefficients
! 
! The image calculation is done in Fourier-space, fully non-linear!
! The envelope for the partial spatial coherence is calculated from 
! a linearized form of the complete transfer function, where the first
! cross-term between partial spatial and partial temporal coherence
! is considered as in Ishizuka (1980) Ultramicroscopy 5.
!
! !!! call PrepareIC(1) before using this routine
!
! - input: complex*8 real-space wave, 2d array of arbitrary size and sampling
! - 1: FFT -> Fourier-space wave
! - 2: Array preparation -> wx, wy, chi, dchix, dchiy
! - 3: 4D-loop
!       - a: outer 2D-loop over image frequencies g
!       - b: inner 2D-loop over beams k
!       - c: sum up on image frequency g
!            TCC(g+k,k)*Psi(g+k)*PsiStar(k)
!            where TCC(g+k,k) = Exp(-ii*Chi(g+k)+ii*Chi(k))*P(g+k,k)*Et(g+k,k)*Es(g+k,k)*Ex(g+k,k),
!                  Et = Exp( -( (g+k)**2 - k**2 )**2 *(pi*wl*fs/2)**2 / u )
!                  Es = Exp( -( GradChi(g+k) - GradChi(k) )**2 *(sc/2)**2 / u)
!                  P  = Exp( 2*ii/u)*(pi*sc*wl*fs/2)**2 *((g+k)**2 - k**2) *((g+k-k).(GradChi(g+k) - GradChi(k))) )
!                  Ex = 1/Sqrt(u)*Exp( -sc**2/u *(pi*sc*wl*fs/2)**2 *((g+k-k) x (GradChi(g+k) - GradChi(k)))**2 )
!                  u  = 1 + (pi*sc*wl*fs)**2 * (g+k-k)^2
!            -> Fourier space image
! - 4: iFFT -> real-space image
! - output: complex*8 real-space image, 2d array of same size and sampling as input
!           the image data is found in the real part of the complex data
!
! In order to speed up the computation of envelope values, the spatial coherence
! envelope is not calculated by an external routine. Instead, the required gradients
! of the aberration function are pre-calculated in extra arrays. The same
! pre-calculation is done for the aberration function.
!
subroutine EnvelopeTCC(wave)

  use wavimgprm
  use AberrationFunctions
  use CSPROG
  
  implicit none

  complex*8, intent(inout) :: wave(nwx,nwy) ! real space wave function (input)
                                            ! real space image (output) (in real-part of wave)
  
  integer*4 :: nx2, ny2                     ! data half sizes
  integer*4 :: i, j, i1, j1, id, jd, i2, j2 ! std. iterators
  real*4 :: chi1, chi2, wx2, w2             ! aberration value
  real*4 :: wi, wj                          ! temp. frequencies
  real*4 :: scpf, sc2                       ! prefactor in spatial coherence envelope
  real*4 :: fspf                            ! prefactor in temporal coherence envelope
  real*4 :: xtpf                            ! prefactor for cross-terms
  real*4 :: Envt, Envs, Envx, u, sqrtu      ! Envelope values
  real*4, allocatable :: wxs(:), wys(:)     ! prepared squared scrambled frequencies
  real*4, allocatable :: dchix(:,:)         ! prepared gradient of aberration function (X)
  real*4, allocatable :: dchiy(:,:)         ! prepared gradient of aberration function (Y)
  complex*8 :: ctf, xph                     ! transfer function and cross-term phase
  complex*8, allocatable :: cctf(:,:)       ! pre-calculated transfer function values
  complex*8, allocatable :: wtmp(:,:)       ! temporary wave data
  complex*8, allocatable :: wtmpc(:,:)      ! conjugate of temporary wave data
  complex*8, allocatable :: wtmp2(:,:)      ! temporary wave data
    
  external :: FSFFT
  external :: PostMessage
  
  ! init
  call PostMessage("Image calculation by envelope TCC integration.")
  ! calculate parameters
  nx2 = nwx/2 ! wave function nyquist index x
  ny2 = nwy/2 ! wave function nyquist index y
  scpf = (sc/2)**2 ! prefactor for partial spatial coherence (semi convergence angle)
  sc2 = sc**2 ! additional prefactor for cross-term envelope
  fspf = (pi*fs/wl/2)**2 ! prefactor for partial temporal coherence (focus spread)
  xtpf = (pi*sc*fs/wl/2)**2 ! prefactor in corss-terms
  ! allocate wave arrays
  allocate(wtmp(nft,nft), wtmpc(nft,nft), wtmp2(nft,nft), stat=nerr)
  allocate(cctf(nwy,nwx), dchix(nwy,nwx), dchiy(nwy,nwx), stat=nerr)
  allocate(wxs(nwx), wys(nwy), stat=nerr)
  wtmp = cmplx(0.0,0.0)
  wxs = 0.0
  wys = 0.0
  
  ! transfer wave data to temp array
  do j=1, nwy
    do i=1, nwx
      wtmp(i,j) = wave(i,j)
    end do
  end do
  
  ! forward FFT
  call FSFFT(wtmp,nwx,nwy,'for')
  wtmpc = conjg(wtmp) ! save the conjugate extra
  
  ! fourier-space wave is now in wtmp, scrambled and transposed
  
  call PostMessage("Preparing temporary array data.")
  
  ! clear up data in array wave
  wave = cmplx(0.0,0.0)
  wtmp2 = cmplx(0.0,0.0)
  
  ! setup scrambled frequencies x
  wxs = ic_twx*ic_twx
  wys = ic_twy*ic_twy
  ! setup aberration function and gradients arrays
  do i=1, nwx
    wi = ic_twx(i)
    do j=1, nwy
      wj = ic_twy(j)
      chi1 = AF_AberrationFunction(wi,wj) ! coherent aberration function
      cctf(j,i) = ic_objaper(j,i)*exp( cmplx(0.0,-chi1) ) ! coherent phase shift * aperture
      call AF_AberrationFunctionGradient(wi, wj, chi1, chi2)
      dchix(j,i) = chi1 ! gradient x of aberration function
      dchiy(j,i) = chi2 ! gradient y of aberration function
    end do
  end do
  
  !
  write(unit=stmp,fmt=*) nwx*nwy
  sdbgmsg = " Starting 4D-loop ("//trim(adjustl(stmp))// &
          & " steps) in Fourier space, this may take a while."
  call PostMessage(trim(sdbgmsg))
  !
  ! non-linear Image formation in Fourier space
  ! I(g1) = Sum_g2 TCC(g1,g2)*Psi(g1+g2)*PsiStar(g2)
  !    g1 = k + g
  !    g2 = k
  
  ! Progress bar
  ! get number of steps
  i = 1+nwx*nwy
  call CSPROG_INIT(i,60)
  
  ! calculation of image in Fourier-space is a weighted auto-correlation
  ! 2D loop over the spatial image frequencies g
  do i=1, nwx ! i -> wx
    wx2 = wxs(i) ! -> wx^2
    do j=1, nwy ! j -> wy
      w2 = wx2 + wys(j) ! -> w^2
      ! get the cross-term u, it depends on (g+k - k) = g
      u = 1. + 4.*xtpf*w2
      sqrtu = sqrt(u)
          
      ! 2D loop over the spatial frequencies
      do i1=1, nwx ! kx: i1 -> g2x = kx
        ! get conjugate frequency g1x = kx + gx
        id = ic_iwx(i1) + ic_iwx(i) ! direct frequency kx + gx
        ! WARNING: these frequencies may be double the highest wave frequency
        !          cut any frequency outsides the range -nyq:nyq-1
        if (id > nx2-1) cycle ! do not calculate frequencies higher than positive wave nyquist
        if (id <  -nx2) cycle ! do not calculate frequencies lower than negative wave nyquist
        ! calculate index corresponding to the direct frequency
        i2 = modulo(id+nwx,nwx)+1 ! index i2 -> g1x = gx + kx
                
        do j1=1, nwy ! ky: j1 -> g2y = ky
          ! get conjugate frequency g1y = ky + gy
          jd = ic_iwy(j1) + ic_iwy(j) ! direct frequency ky + gy
          ! WARNING: these frequencies may be double the highest wave frequency
          !          cut any frequency outsides the range -nyq:nyq-1
          if (jd > ny2-1) cycle ! do not calculate frequencies higher than positive wave nyquist
          if (jd <  -ny2) cycle ! do not calculate frequencies lower than negative wave nyquist
          ! change the direct frequency back to scrambled array index
          j2 = modulo(jd+nwy,nwy)+1 ! index j2 -> g1y = gy + ky
             
          ! get the coherent transfer function for g+k and k
          ctf = cctf(j2,i2)*conjg(cctf(j1,i1))
          
          ! get the cross-term phase
          ! P  = Exp( 2*ii/u)*(pi*sc*wl*fs/2)**2 *((g+k)**2 - k**2) *((g+k-k).(GradChi(g+k) - GradChi(k))) )
          xph = exp( cmplx(0., 2*xtpf/u *( wxs(i2)+wys(j2)-wxs(i1)-wys(j1) ) &
              &                *(  ic_twx(i)*(dchix(j2,i2)-dchix(j1,i1)) & 
              &                  + ic_twy(j)*(dchiy(j2,i2)-dchiy(j1,i1))   ) ) )
          
          ! get the temporal coherence envelope
          Envt = exp( -fspf/u*( ( wxs(i2)+wys(j2) -wxs(i1)-wys(j1) )**2) )
          ! get the spatial coherence envelope
          Envs = exp( -scpf/u*( (dchix(j2,i2)-dchix(j1,i1))**2 + (dchiy(j2,i2)-dchiy(j1,i1))**2 ) )
          ! get the cross-term envelope
          ! Ex = 1/Sqrt(u)*Exp( -sc**2/u *(pi*sc*wl*fs/2)**2 *((g+k-k) x (GradChi(g+k) - GradChi(k)))**2 )
          !                                                        ( g x dgc ) = (gx*dgcy - gy*dgcx)
          Envx = exp( -sc2*xtpf/u *(  ic_twx(i)*(dchiy(j2,i2)-dchiy(j1,i1)) &
               &                    - ic_twy(j)*(dchix(j2,i2)-dchix(j1,i1)) )**2 )/sqrtu
      
          ! add up to the image in Fourier space
          wtmp2(j,i) = wtmp2(j,i) + wtmp(j2,i2)*wtmpc(j1,i1)*ctf*xph*Envs*Envt*Envx
        
        end do
      end do
      
      FSTEP_CUR = FSTEP_CUR+1
      call CSPROG_UPDATE()
      
    end do
  end do
  
  ! backward FFT
  call FSFFT(wtmp2,nwx,nwy,'bac')
  
  ! transfer image data to wave array
  do j=1, nwy
    do i=1, nwx
      wave(i,j) = real(wtmp2(i,j))
    end do
  end do
  
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  call CSPROG_FINISH()
  
  deallocate(wtmp,wtmp2,wtmpc, stat=nerr)
  deallocate(cctf,dchix,dchiy, stat=nerr)
  deallocate(wxs, wys, stat=nerr)
  
  return

end subroutine EnvelopeTCC
! --------------------------------------------------------------------






! --------------------------------------------------------------------
!
!  subroutine LinearEnvelopeTCC
!
! Envelope-based image calculation by transmission cross-coefficients
!
! The calculation is very similar to the non-linear routine above,
! except that the full convolution over all nonlinear terms is skipped.
! only the 2 linear terms Psi(g)*PsiStar(0) + Psi(0)*PsiStar(-g) are
! taken into account.
! 
! The image calculation is done in Fourier-space, linear!
! The envelope for the partial spatial coherence is calculated from 
! a linearized form of the complete transfer function, where cross-terms
! between partial spatial and partial temporal coherence are not considered.
!
! !!! call PrepareIC(1) before using this routine
!
! - input: complex*8 real-space wave, 2d array of arbitrary size and sampling
! - 1: FFT -> Fourier-space wave
! - 2: Array preparation -> wx, wy, chi, dchix, dchiy
! - 3: 2D-loop over image frequencies g
!       - a: sum up on image frequency g
!            TCC(g,0)*Psi(g)*PsiStar(0) + TCC(0,-g)*Psi(0)*PsiStar(-g)
!            where TCC(g1,g2) = EnvTemporal(g1,g2)*EnvSpatial(g1,g2)*Exp(-ii*Chi(g1)+ii*Chi(g2)),
!                  EnvTemporal = Exp( -( g1**2 - g2**2 )**2 *(pi*wl*fs/2)**2 )
!                  EnvSpatial  = Exp( -( GradChi(g1) - GradChi(g2) )**2 *(sc/2)**2 )
!            -> Fourier space image
! - 4: iFFT -> real-space image
! - output: complex*8 real-space image, 2d array of same size and sampling as input
!           the image data is found in the real part of the complex data
!
! In order to speed up the computation of envelope values, the spatial coherence
! envelope is not calculated by an external routine. Instead, the required gradients
! of the aberration function are pre-calculated in extra arrays. The same
! pre-calculation is done for the aberration function.
!
subroutine LinearEnvelopeTCC(wave)

  use wavimgprm
  use AberrationFunctions
  use CSPROG
  
  implicit none
  
  complex*8, intent(inout) :: wave(nwx,nwy) ! real space wave function (input)
                                            ! real space image (output) (in real-part of wave)
  
  integer*4 :: nx2, ny2                     ! data half sizes
  integer*4 :: i, j, i1, j1, i2, j2         ! std. iterators
  real*4 :: chi1, chi2                      ! aberration value
  real*4 :: scpf                            ! prefactor in spatial coherence envelope
  real*4 :: fspf                            ! prefactor in temporal coherence envelope
  real*4 :: Envt, Envs                      ! Envelope values
  real*4, allocatable :: wxs(:), wys(:)     ! prepared squared scrambled frequencies [rad^2]
  real*4, allocatable :: dchix(:,:)         ! prepared gradient of aberration function (X)
  real*4, allocatable :: dchiy(:,:)         ! prepared gradient of aberration function (Y)
  complex*8 :: ctf                          ! transfer function
  complex*8, allocatable :: cctf(:,:)       ! pre-calculated transfer function values
  complex*8, allocatable :: wtmp(:,:)       ! temporary wave data
  complex*8, allocatable :: wtmp2(:,:)      ! temporary wave data
  complex*8, allocatable :: wtmpc(:,:)      ! conjugate of temporary wave data
  
  external :: FSFFT ! link 'wavimgsubs.f90'
  
  ! init
  call PostMessage("Linear image calculation by envelope TCC integration.")
  ! calculate parameters
  nx2 = nwx/2 ! wave function nyquist index x
  ny2 = nwy/2 ! wave function nyquist index y
  scpf = (sc/2)**2 ! prefactor for partial spatial coherence (semi convergence angle)
  fspf = (pi*fs/wl/2)**2 ! prefactor for partial temporal coherence (focus spread)
  
  ! allocate wave arrays
  allocate(wtmp(nft,nft), wtmpc(nft,nft), wtmp2(nft,nft), stat=nerr)
  ! allocate TCC helper arrays
  allocate(cctf(nwy,nwx), dchix(nwy,nwx), dchiy(nwy,nwx), stat=nerr)
  ! allocate diffraction angle arrays
  allocate(wxs(nwx), wys(nwy), stat=nerr)
  wtmp = cmplx(0.0,0.0)
  wtmpc = cmplx(0.0,0.0)
  wtmp2 = cmplx(0.0,0.0)
  cctf = cmplx(0.0,0.0)
  dchix = 0.0
  dchiy = 0.0
  wxs = ic_twx*ic_twx
  wys = ic_twy*ic_twy
  
  ! get number of steps for progress indicator
  i = 3+nwx*nwy
  call CSPROG_INIT(i,60)
  
  ! transfer wave data to temp array
  do j=1, nwy
    do i=1, nwx
      wtmp(i,j) = wave(i,j)
    end do
  end do
  
  ! forward FFT
  call FSFFT(wtmp,nwx,nwy,'for')
  wtmpc = conjg(wtmp) ! save PsiStar, the complex conjugate of Psi
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! fourier-space wave is now in wtmp, scrambled and transposed
  
  ! clear up data in array wave
  wave = cmplx(0.0,0.0)
  
  ! setup aberration function and gradients arrays
  do i=1, nwx
    do j=1, nwy
      chi1 = AF_AberrationFunction(ic_twx(i), ic_twy(j)) ! coherent aberration function
      cctf(j,i) = ic_objaper(j,i)*exp( cmplx(0.0,-chi1) ) ! coherent phase shift * aperture
      call AF_AberrationFunctionGradient(ic_twx(i), ic_twy(j), chi1, chi2)
      dchix(j,i) = chi1 ! gradient x of aberration function
      dchiy(j,i) = chi2 ! gradient y of aberration function
    end do
  end do
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! linear Image formation in Fourier space
  ! I(g) =  TCC(g,0)*Psi(g)*PsiStar(0)+TCC(0,-g)*Psi(0)*PsiStar(-g)
  
  ! calculation of image in Fourier-space is a weighted auto-correlation
  ! 2D loop over the spatial frequencies g1
  do i=1, nwx
    
    i1 = i ! gx index
    i2 = modulo(nwx-ic_iwx(i),nwx)+1  ! -gx index
    
    do j=1, nwy
    
      j1 = j ! gy index
      j2 = modulo(nwy-ic_iwy(j),nwy)+1  ! -gy index
    
      ! interference for g1=g and g2=0
      ! - ctf
      ctf = cctf(j1,i1)*conjg(cctf(1,1))
      ! - temporal envelope
      Envt = exp( -fspf*( ( wxs(i1) + wys(j1) )**2) )
      ! - spatial envelope
      Envs = exp( -scpf*( (dchix(j1,i1)-dchix(1,1))**2 + (dchiy(j1,i1)-dchiy(1,1))**2 ) )
      ! - add up to the image in Fourier space
      wtmp2(j,i) = wtmp2(j,i) + wtmp(j1,i1)*wtmpc(1,1)*ctf*Envs*Envt
        
      ! interference for g1=0 and g2=-g
      ! - ctf
      ctf = cctf(1,1)*conjg(cctf(j2,i2))
      ! - temporal envelope
      Envt = exp( -fspf*( ( wxs(i2) +wys(j2) )**2) )
      ! - spatial envelope
      Envs = exp( -scpf*( (dchix(1,1)-dchix(j2,i2))**2 + (dchiy(1,1)-dchiy(j2,i2))**2 ) )
      ! - add up to the image in Fourier space
      wtmp2(j,i) = wtmp2(j,i) + wtmp(1,1)*wtmpc(j2,i2)*ctf*Envs*Envt
      
      FSTEP_CUR = FSTEP_CUR+1
      call CSPROG_UPDATE()
      
    end do
  end do
  
  ! backward FFT
  call FSFFT(wtmp2,nwx,nwy,'bac')
  
  ! transfer image data to wave array
  do j=1, nwy
    do i=1, nwx
      wave(i,j) = real(wtmp2(i,j))
    end do
  end do
  
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  call CSPROG_FINISH()
  
  deallocate(wtmp,wtmp2,wtmpc, stat=nerr)
  deallocate(cctf,dchix,dchiy, stat=nerr)
  deallocate(wxs, wys, stat=nerr)
  
  return

end subroutine LinearEnvelopeTCC
! --------------------------------------------------------------------






! --------------------------------------------------------------------
!
!  subroutine QuasiCoherent
!
! Envelope-based image calculation (quasi-coherent approximation)
!
! Contrast dampening is applied to the beams of the wave function directly
! in form of the linear envelopes for partial temporal (Et) and
! partial spatial (Es) coherence. A dampened wave function is calculated,
! then inverse Fourier transformed and the absolute square is taken
! in real space.
!
! Warning: This approximation leads to a loss of intensity, and a
!          too strong dampening of the non-linear interference terms.
! 
! !!! call PrepareIC(1) before using this routine
!
! - input: complex*8 real-space wave, 2d array of arbitrary size and sampling
! - 1: FFT -> Fourier-space wave
! - 2: Array preparation -> wx, wy, chi, dchix, dchiy
! - 3: 2D-loop over image frequencies g
!       - multiply envelopes and aberration function
!         tau(g) = Exp[ -i*Chi(g))
!         Et(g) = Exp( -g**4 *(pi*wl*fs/2)**2 )
!         Es(g) = Exp( -(GradChi(g))**2 *(sc/2)**2 )
! - 4: iFFT -> image wave function
! - 5: Absolute square -> image intensity
! - output: complex*8 real-space image, 2d array of same size and sampling as input
!           the image data is found in the real part of the complex data
!
! In order to speed up the computation of envelope values, the spatial coherence
! envelope is not calculated by an external routine. Instead, the required gradients
! of the aberration function are pre-calculated in extra arrays. The same
! pre-calculation is done for the aberration function.
!
subroutine QuasiCoherent(wave)

  use wavimgprm
  use AberrationFunctions
  use CSPROG
  
  implicit none
  
  complex*8, intent(inout) :: wave(nwx,nwy) ! real space wave function (input)
                                            ! real space image (output) (in real-part of wave)
  
  integer*4 :: nx2, ny2                     ! data half sizes
  integer*4 :: i, j                         ! iterators
  real*4 :: chi1, chi2                      ! aberration value
  real*4 :: scpf                            ! prefactor in spatial coherence envelope
  real*4 :: fspf                            ! prefactor in temporal coherence envelope
  real*4 :: Envt, Envs                      ! Envelope values
  real*4, allocatable :: wxs(:), wys(:)     ! prepared squared scrambled frequencies [rad^2]
  real*4, allocatable :: dchix(:,:)         ! prepared gradient of aberration function (X)
  real*4, allocatable :: dchiy(:,:)         ! prepared gradient of aberration function (Y)
  complex*8 :: ctf                          ! transfer function
  complex*8 :: cval                         ! a complex number
  complex*8, allocatable :: cctf(:,:)       ! pre-calculated transfer function values
  complex*8, allocatable :: wtmp(:,:)       ! temporary wave data
  complex*8, allocatable :: wtmp2(:,:)      ! temporary wave data
  
  external :: FSFFT ! link 'wavimgsubs.f90'
  
  ! init
  call PostMessage("Quasi-coherent image calculation.")
  ! calculate parameters
  nx2 = nwx/2 ! wave function nyquist index x
  ny2 = nwy/2 ! wave function nyquist index y
  scpf = (sc/2)**2 ! prefactor for partial spatial coherence (semi convergence angle)
  fspf = (pi*fs/wl/2)**2 ! prefactor for partial temporal coherence (focus spread)
  
  ! allocate wave arrays
  allocate(wtmp(nft,nft), wtmp2(nft,nft), stat=nerr)
  ! allocate TCC helper arrays
  allocate(cctf(nwy,nwx), dchix(nwy,nwx), dchiy(nwy,nwx), stat=nerr)
  ! allocate diffraction angle arrays
  allocate(wxs(nwx), wys(nwy), stat=nerr)
  wtmp = cmplx(0.0,0.0)
  wtmp2 = cmplx(0.0,0.0)
  cctf = cmplx(0.0,0.0)
  dchix = 0.0
  dchiy = 0.0
  wxs = ic_twx*ic_twx
  wys = ic_twy*ic_twy
  
  ! get number of steps for progress indicator
  i = 3 + nwx*nwy
  call CSPROG_INIT(i,60)
  
  ! transfer wave data to temp array
  do j=1, nwy
    do i=1, nwx
      wtmp(i,j) = wave(i,j)
    end do
  end do
  
  ! 1) forward FFT
  call FSFFT(wtmp,nwx,nwy,'for')
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! fourier-space wave is now in wtmp, scrambled and transposed
  
  ! clear up data in array wave
  wave = cmplx(0.0,0.0)
  
  ! 2) setup aberration function and gradients arrays
  do i=1, nwx
    do j=1, nwy
      chi1 = AF_AberrationFunction(ic_twx(i), ic_twy(j)) ! chi(g) coherent aberration function
      cctf(j,i) = ic_objaper(j,i)*exp( cmplx(0.0,-chi1) ) ! tau(g) coherent phase shift * aperture
      call AF_AberrationFunctionGradient(ic_twx(i), ic_twy(j), chi1, chi2)
      dchix(j,i) = chi1 ! gradient x of aberration function
      dchiy(j,i) = chi2 ! gradient y of aberration function
    end do
  end do
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  
  ! 3) Calculation of a dampened image wave function in Fourier space
  !    Psi'(g) = Psi(g)*tau(g)*Es(g)*Et(g)
  !    2D loop over the spatial frequencies g
  do i=1, nwx
    
    do j=1, nwy
    
      ! - ctf
      ctf = cctf(j,i)
      ! - temporal envelope
      Envt = exp( -fspf*( ( wxs(i) + wys(j) )**2) )
      ! - spatial envelope
      Envs = exp( -scpf*( (dchix(j,i)-dchix(1,1))**2 + (dchiy(j,i)-dchiy(1,1))**2 ) )
      ! - image wave coefficient
      wtmp2(j,i) = wtmp(j,i) * ctf * Envs * Envt
      
      FSTEP_CUR = FSTEP_CUR+1
      call CSPROG_UPDATE()
      
    end do
  end do
  
  ! 4) backward FFT
  call FSFFT(wtmp2,nwx,nwy,'bac')
  !    absolute square -> image
  do j=1, nwy
    do i=1, nwx
      cval = wtmp2(i,j)
      wave(i,j) = cmplx( real( cval*conjg(cval) ), 0. )
    end do
  end do
  
  FSTEP_CUR = FSTEP_CUR+1
  call CSPROG_UPDATE()
  call CSPROG_FINISH()
  
  deallocate(wtmp,wtmp2, stat=nerr)
  deallocate(cctf,dchix,dchiy, stat=nerr)
  deallocate(wxs, wys, stat=nerr)
  
  return

end subroutine QuasiCoherent
! --------------------------------------------------------------------






