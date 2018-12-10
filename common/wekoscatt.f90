! ******************************************************************** !
!
! Implementation of Weickenmeier & Kohl electron atomic form factors
!
! Acta Cryst. A 47 (1991) 590 - 597
!
! by J. Barthel, Forschungszentrum Jülich GmbH, Jülich, Germany
!                RWTH Aachen University, Aachen, Germany
! 2018-May-29
!
! requires extended source lines ( option "/extend_source:132" )
!
! depends on "integration.f90"
!
! ******************************************************************** !
!
! Original implementation by A. Weickenmeier (wscatt.f)
! Modified: 
!		A.W.	07.08.90
!		P.S.	01.14.91
!       J.B.    23.11.2010
!       J.B.    03.03.2011 (added parameters for H, fitted from E. Kirklands table.)
!       J.B.    14.11.2014 (added function "dwfjbr" returning a dwf.)
!       J.B.    30.06.2017 (added parameters for L, z=0, vacancy, V&B-parameters play no role)
!       J.B.    29.05.2018 (added f^2(theta) integrals, wekof2)
! 
! ******************************************************************** !
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

! ******************************************************************** !
!
! Interfaces
!
! Scattering factors are returned in [nm] and after relativistic
! correction.
!
! wekosca(g,dw,z,akv,dwfflg,absflg) : (function) (complex*8) [nm]
!                     returns the atomic form factor for electron
!                     scattering for a given atom type z at spatial
!                     frequency g [1/nm], with Debye-Waller parameter
!                     dw [nm^2] = Biso and electron energy akv [keV]
!                     real-part = atomic form factor [nm]
!                     imag-part = absorptive form factor [nm]
! wekoscar(g,dw,z,akv,dwfflg) : (function) (real*4) [nm]
!                     returns the atomic form factor for electron
!                     scattering for a given atom type z at spatial
!                     frequency g [1/nm], with Debye-Waller parameter
!                     dw [nm^2] = Biso and electron energy akv [keV]
! wekoscar1(s)      : (function) (real*8) [A] (only for internal use)
!                     returns the atomic form factor for electron
!                     scattering for scattering vector magnitude s.
!                     s = g/2 [1/A] (real*8)
!                     no relativistic correction and no DWF applied.
!                     atomic form factor parameters wakia are
!                     defined in common block /waki1/
! wekoabs(g,dw,a,b,k0) : function) (real*4) [A]
!                     returns the absorptive form factor for given
!                     Weickenmeier & Kohl parameters a(2) and b(6),
!                     Debye-Waller parameter dw [A^2] = Biso, and
!                     spatial frequency g [1/A] = scattering vector
!                     g = 2*s [1/A] (real*4)
!                     wave number k0 [1/A]
! wekof2(tmax,dw,a,b,k0) : (function) (real*4) [A^2]
!                     Returns the integral of f^2 from theta=0 to
!                     theta=tmax for given
!                     Weickenmeier & Kohl parameters a(2) and b(6),
!                     Debye-Waller parameter dw [A^2] = Biso and
!                     wave number k0 [1/A]
! wekosymb(z, symb) : returns the atom symbol symb for atomic number z
! getweko(z, a, b)  : returns the parameterization coefficients
!                     a(14) for atomic number z
! real*4 :: weko(a, b, s) : (function) (real*4)
!                     returns the electron atomic form factor fe(s) [A]
!                     for given parameterization coefficients a(14)
!                     and scattering vector magnitude s [1/A]
!                     s = k*arcsin(theta/2) = g/2
!                     doesn't contain the ionic charge term for ions.
! wekoimag(g,ul,a,b) : (function) (real*4)
!                     returns the absorptive form factor using
!                     Weickenmeiers analytical formula.
!
! First use getweko to obtain the coefficients of the parameterization.
! Then use weko to calculate fe(s).
! Don't forget to add the ionic charge potential later.
!
! ******************************************************************** !

! ******************************************************************** !
!	
!	Note:
!	=====
!
!	The scattering amplitudes are multiplied by an additional 
!	factor 4pi.
!	
! ******************************************************************** !


! ******************************************************************** !	
!
! complex*8 function wekosca
!
! input
!   - real     g       spatial frequency [1/nm]
!   - real     dw      debye-waller parameter in nm^2
!   - integer  z       atomic number
!   - akv      akv     high tension in kV
!   - logical  dwflg   flag: use debye-waller factor
!   - logical  absflg  flag: calculate absorptiion potential
!
! output / return value
!   - complex*8        complex scattering factor for the given
!                      atom type at the given spatial frequency
!
function wekosca(g,dw,z,akv,dwfflg,absflg)
!
  implicit none
!
  real*4, parameter :: twopi	= 6.28318530718 ! 2*Pi
  real*4, parameter :: fourpi	= 12.56637061436 ! 4*Pi
  real*4, parameter :: r8pi2	= 0.0126651479553 ! 1/(8*Pi^2)
  real*4, parameter :: tiny     = 1.0E-8 ! small value threshold
  real*4, parameter :: elrmkv   = 510.998947 ! m0*c^2 in keV (electron rest mass)
!  
  real*4, intent(in) :: g, dw, akv
  integer*4, intent(in) :: z
  logical, intent(in) :: dwfflg, absflg ! calculation flags
!  
  complex*8 :: wekosca
!  
  real*4 :: dwa, ga, sa, ua, dewa
  real*4 :: k0, rc
  real*4 :: fr, fi
  real*4 ::	a(2), b(6) ! coefficients of the weko parameterization
  real*4, external :: weko, wekoimag, getdwf, wekoabs
  external :: getweko
!	
! relativistic correction
  rc	= (elrmkv + akv) / elrmkv
! Debye-Waller factor
  dewa  = getdwf(g, dw, dwfflg)
! [nm] to [ang]
  ga	= 0.100 * g * twopi ! translate g from 1/nm to 1/A and times 2*pi
  sa	= 0.050 * g ! sa is ga / 4pi = g / 20
  dwa   = 100.0 * dw
  ua	= sqrt( dwa * r8pi2) ! translate from DW-parameter to RMS vibr amp [A]
! get the WK-coefficients for atom z
  call	getweko(z, a, b) ! 
! form factor, relativistic correction, 4*Pi and DWF
  fr 	= rc * weko(a, b, sa) * fourpi * dewa ! real part of the scattering factor
  fi	= 0.0 ! preset imag part to zero
! calculate imag scatt      
  if ((dwfflg).and.(absflg).and.(abs(ua)>=tiny).and.(abs(fr)>=tiny)) then
    !wave number in [A^-1] 0.506774 = 2*Pi*e/(h*c) * 1E-7 [A^-1 * kV^-1]
    k0  = 0.506774 * sqrt((2*elrmkv + akv) * akv) ! 2*Pi*k ! andere k-Notation hier: Aufpassen!
    fi  = rc*rc*wekoimag (ga, ua, a, b) / k0
    !fi  = rc*rc*wekoabs(0.1*g, dwa, a, b, k0/twopi)*k0
  end if
! [ang] to [nm]
  wekosca = 0.1 * cmplx (fr, fi)
  return
!
end function wekosca
! ******************************************************************** !
	
	
! ******************************************************************** !
!
! real*4 function wekoscar
!
! input
!   - real     g       spatial frequency [1/nm]
!   - real     dw      debye-waller parameter in nm^2
!   - integer  z       atomic number
!   - akv      akv     high tension in kV
!   - logical  dwflg   flag: use debye-waller factor
!
! output / return value
!   - real*4           elastic scattering factor for the given
!                      atom type at the given spatial frequency
!
!
function wekoscar(g,dw,z,akv,dwfflg)
!
  implicit none
!
  real*4, parameter :: twopi	= 6.28318530718 ! 2*Pi
  real*4, parameter :: fourpi	= 12.56637061436 ! 4*Pi
  real*4, parameter :: r8pi2	= 0.0126651479553 ! 1/(8*Pi^2)
  real*4, parameter :: elrmkv   = 510.998947 ! m0*c^2 in keV (electron rest mass)
!
!  
  logical, intent(in) :: dwfflg ! calculation flags
  real*4, intent(in) :: g, dw, akv
  integer*4, intent(in) :: z
!  
  real*4 :: wekoscar
!  
  real*4 :: sa, dewa
  real*4 :: rc
  real*4 :: fr
  real*4 ::	a(2), b(6) ! coefficients of the weko parameterization
  real*4, external :: weko, getdwf
  external :: getweko
!
! relativistic correction
  rc	= (elrmkv + akv) / elrmkv
! debye-waller factor
  dewa  = getdwf(g, dw, dwfflg)
! scattering vector magnitude
  sa	= 0.050 * g  ! s = g/2 in [ang]
! coefficients of the parameterization
  call	getweko(z, a, b) ! get the WK-coefficients for atom z
! form factor with relativistic correction, DWF  and  x 4*Pi
  fr 	= rc * weko(a, b, sa) * fourpi * dewa 
! [ang] to [nm]
  wekoscar = 0.1 * fr
!
  return
!
end function wekoscar
! ******************************************************************** !





! ******************************************************************** !
!
function getdwf(g,dw,dwfflg)
  implicit none
  logical, intent(in) :: dwfflg
  real*4, intent(in) :: g, dw ! g [nm] and B [nm^2]
  real*4 :: getdwf
  real*4 :: dewa
!
  dewa  = 1.0
  if (dwfflg) then
    dewa = exp( -0.25 * dw * g * g )
  end if
  getdwf = dewa
  return
end
! ******************************************************************** !

	
	
! ******************************************************************** !
!
function weko(a, b, s)
!
!	electron scattering apmlitude f(s) [A]
  implicit none
  real*4, intent(inout) :: a(2)
  real*4, intent(inout) :: b(6)
  real*4, intent(in) :: s
  real*4 :: s2, argu
  real*4 :: weko
  integer*4 :: i, j
!
  weko	= 0.0
  if (s .gt. 0.0) then
    s2 = 1./(s*s)
  end if
  do i=1,6
	j = 1+(i-1)/3
	argu = b(i)*s*s
	if (argu < 0.1) then
	  weko = weko + a(j)*b(i) * (1.-.5*argu)
	else if (argu > 20.0) then
	  weko = weko + a(j)*s2
	else
	  weko = weko + a(j)*(1.-exp(-argu))*s2
	end if
  end do
!  
  return
end function weko
! ******************************************************************** !



! ******************************************************************** !
!
function wekoimag(g,ul,a,b)
!
! calculate absorptive form factors for
! "g" = g * 2Pi [ang]
! "ul" = sqrt(B/(8 Pi^2)) ! RMS vibr amp [A]
! a, b = coefficients of We&Ko
  implicit none
  real*4, parameter :: fourpi	= 12.56637061436 ! 4*Pi
  real*4, parameter :: fp2	    = fourpi*fourpi ! (4*Pi)^2
  real*4, intent(in) :: g, ul, a(2), b(6)
  real*4 :: wekoimag
  real*4 :: a1(2), b1(6), u2, fi, g2, dewa
  integer*4 :: i, j, jj, ii
  real*4, external :: ri1, ri2
  u2 = ul*ul
  do i=1,2
    a1(i) = a(i) * fp2 
  end do
  do i=1,6
    b1(i) = b(i) / fp2
  end do
  fi = 0.
  g2 = g*g
  dewa  = exp(-.5*u2*g2)
  do j=1,6
	jj = 1+(j-1)/3
	do i=1,6
	  ii = 1+(i-1)/3
	  fi = fi + a1(jj)*a1(ii)*(  dewa * ri1(b1(i),b1(j),g) - &
                              &  ri2(b1(i),b1(j),g,ul)  )
    end do
  end do
  wekoimag = fi	
  return
end function wekoimag
! ******************************************************************** !




! ******************************************************************** !
!
! wekoabs
!
! This function starts a numerical integration of the absorptive
! form factor for electron scattering form factors given in the
! Weickenmeier & Kohl parameterization.
!
! input: (real*4)
!   ga = diffraction vector magnitude [1/A]
!   dwa = Debye-Waller parameter (Biso) [A^2]
!   a(2) = parameters ai from We&Ko's table (defines the atom)
!   b(6) = parameters bi from We&Ko's table (defines the atom)
!   k0 = incident beam vacuum wave vector magnitude [1/A] = 1/Lambda
!      = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
! output: (real*4)
!   wakiabs = absorptive form factor without the required pre-factor
!             DWF * gamma^2 / (2*pi*k0)
!
function wekoabs(ga,dwa,a,b,k0)
!
  implicit none
! 
  real*4, intent(in) :: ga, dwa, a(2), b(6), k0
! 
  real*4 :: wekoabs
! 
  real*4 :: wekoa(2), wekob(6)
  real*8 :: si0, si1, fa, pi0, pi1
  real*8 :: wekog, wekodw, wekok
! common blocks ...
  common /weko1/ wekoa, wekob !  for wekoscar1
  common /weko2/ wekog, wekodw, wekok !  for integrands wakimu0 and wakimug
! 
  real*8, external :: dsgrid2d ! link integration.f90 !
  real*8, external :: wekomug ! integrand function
! 
  wekoabs = 0.0
  wekoa = a ! set the We&Ko parameters ai to be used
  wekob = b ! set the We&Ko parameters bi to be used
  wekog = dble(ga) ! store ga for the intagrand wekomug
  wekodw = dble( dwa ) ! store Debye-Waller parameter for both integrands
  wekok = k0 ! store teh incident electron wave vector [1/A]
! set integration range
  si0 = 0.d+0 ! 0
  si1 = 3.1415926535898D+0 ! Pi
  pi0 = 0.d+0 ! 0
  pi1 = 3.1415926535898D+0 ! Pi
! call the numerical integrator
! use a grid of 128 pixels along theta with a square sampling
!           and  64 pixels along phi with linear sampling
!           phi on the half side only (symmteric) therefore:  * 2
  fa = dsgrid2d(wekomug,si0,si1,pi0,pi1,2.0d+0,1.0d+0,128,64) * 2.0D+0
  !
  wekoabs = real(fa, kind=4 )
  return
end function wekoabs
!
! ******************************************************************** !



!**********************************************************************!
!
! function wekomug(theta,phi)
!
! returns the integrand for the absorptive form factor at g=wekog
! for debye-waller parameter b=wekob
! for incident electron wave vector k=wekok
! for atom defined by ai = wekoa and bi=wekob
!
! Set the values of common block parameters before using this function!
! uses common blocks /weko1/ and /weko2/
!
function wekomug(theta,phi)
  
  implicit none
  
  real*8, parameter :: twopi = 6.283185307179586476925286766559
  
  real*8, intent(in) :: theta, phi
  real*8 :: wekomug
  real*8 :: k, q, g, qg, biso, ct, st, cp, twok, sg, sq, sqg
  real*8 :: fq, fqg, dwfg, dwfq, dwfqg
  real*8 :: wekog, wekodw, wekok
  common /weko2/ wekog, wekodw, wekok
  
  real*8, external :: wekoscar1 ! form factor function
                                ! uses common block /weko1/
  ct = dcos(theta)
  st = dsin(theta)
  cp = dcos(phi)
  ! q = scattering vector length of Q in the ewald sphere
  ! Q = (qx, qy, qz) = K' - K
  ! qx = k * sin(theta) * cos(phi)
  ! qy = k * sin(theta) * sin(phi)
  ! qz = k * ( cos(theta) - 1 )
  ! k = |K| = wekok
  k = wekok
  twok = k+k
  ! q = Sqrt[ 2 k^2 (1 - Cos[q]) ]
  q = k * dsqrt( 2.D+0 - 2.D+0 * ct )
  ! g = length of some reciprocal space vector G = (gx, gy, gz)
  !     with gx = g, gy = 0, gz = 0 (obda)
  g = wekog
  ! qg = length of the difference vector Q - G
  qg = dsqrt( g*g + twok*k - twok*(k*ct + g*cp*st) )
  biso = wekodw
  sg = g * 0.5D+0
  sq = q * 0.5D+0
  sqg = qg * 0.5D+0
  fq = wekoscar1(sq) ! f(g/2)
  fqg = wekoscar1(sqg)
  dwfg = dexp( -biso*sg*sg ) ! exp(-biso*sg^2) = exp( -1/4 *biso*g^2)
  dwfq = dexp( -biso*sq*sq )
  dwfqg = dexp( -biso*sqg*sqg )
  wekomug = fq*fqg * (dwfg - dwfq*dwfqg) * st
  
  return
    
end function wekomug




! ******************************************************************** !
!
! wekof2
!
! This function starts a numerical integration of the squared form
! factor for electron scattering given in the
! We & Ko parameterization.
!
! input: (real*4)
!   tmax = max. scattering angle [rad] (upper limit is Pi)
!   dwa = Debye-Waller parameter (Biso) [A^2]
!   a(2) = a parameters from We&Ko's table (defines the atom)
!   b(6) = b parameters from We&Ko's table (defines the atom)
!   k0 = incident beam vacuum wave vector magnitude [1/A] = 1/Lambda
!      = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
! output: (real*8)
!   wekof2 = integrated squared form factor [A^2] no rel. corr.
!
function wekof2(tmax,dwa,z,k0)
!
  implicit none
!
  integer*4, parameter :: ngrid = 2048
  real*8, parameter :: twopi = 6.2831853071796D+0
!
  real*4, intent(in) :: tmax, dwa, k0
  integer*4, intent(in) :: z
!
  real*8 :: wekof2
!
  real*4 :: wekoa(2), wekob(6)
  real*8 :: si0, si1, fa
  real*8 :: wekog, wekodw, wekok
! common blocks ...
  common /weko1/ wekoa, wekob        !  for wekoscar1
  common /weko2/ wekog, wekodw, wekok !  for integrand wekof2g
!
  real*8, external :: dsgrid1d ! link integration.f90 !
  real*8, external :: wekof2g ! integrand functions

!
  wekof2 = 0.0
  call getweko(z, wekoa, wekob) ! set the Weickenmeier & Kohl parameters to be used
  wekog = 0.D+0 ! not used
  wekodw = dble( dwa ) ! store Debye-Waller parameter
  wekok = k0 ! store wave number [1/A]
! set integration range
  si0 = 0.d+0 ! 0
  si1 = min(3.1415926535898D+0, dble(abs(tmax)) ) ! min(Pi,abs(tmax))
! call the numerical integrator
! - f2 is constant with phi, therefore:  * 2*Pi
  fa = dsgrid1d(wekof2g,si0,si1,2.0d+0,ngrid) * twopi
  !
  wekof2 = fa
  return
end function wekof2
!
! ******************************************************************** !



!**********************************************************************!
!
! function wekof2g(theta)
!
! returns the integrand f(theta)^2 * sin(theta)
! for debye-waller parameter b=wekodw (set to zero if not used)
! for incident electron wave vector k=wekok
!     -> s = q/2 = k*sin(theta/2)
!                = k*sqrt( (1-cos(theta))/2 )
! for atom defined by a = wekoa and b = wekob
!
! Set the values of common block parameters before using this function!
! uses common blocks /weko1/ and /weko2/
!
function wekof2g(theta)
  
  implicit none
  
  real*8, intent(in) :: theta
  real*8 :: wekof2g
  real*8 :: k, q, biso, ct, st, sq
  real*8 :: fq
  real*8 :: wekog, wekodw, wekok
  common /weko2/ wekog, wekodw, wekok
  
  real*8, external :: wekoscar1 ! form factor function
                                ! uses common block /weko1/
  
  ct = dcos(theta)
  st = dsin(theta)
  ! q = scattering vector length of Q in the ewald sphere
  !     for scattering angle theta
  ! Q = (qx, qy, qz) = K' - K
  ! qx = k * sin(theta) * cos(phi)
  ! qy = k * sin(theta) * sin(phi)
  ! qz = k * ( cos(theta) - 1 )
  ! k = |K| = wakik
  k = wekok
  ! q = Sqrt[ 2 k^2 (1 - Cos[q]) ]
  q = k * dsqrt( 2.D+0 - 2.D+0 * ct )
  biso = wekodw
  sq = 0.5D+0 * q ! q/2
  fq = wekoscar1(sq) ! f(s) = f(q/2)
  if (dabs(biso)>0.D+0) then
    fq = fq * dexp( -biso*sq*sq ) ! DWF = EXP( -1/4 * BISO * q^2 ) = EXP( -BISO * s^2 )
  end if
  wekof2g = st*fq*fq ! Sin(theta)*(f(theta)*DWF))^2
  
  return
    
end function wekof2g






! ******************************************************************** !
!
! real*8 function wekoscar1
!
! used only internally for the integration of absorptive form factors
!
! input
!   - real*8    s       scattering vector amplitude [1/A] = g/2
!
! output / return value
!   - real*8            elastic scattering factor for the given s
!
! the coeffcients of the parameterization are taken
! from common block /weko1/ wekoa, wekob
! - neutral atomic form factors only
! - no relativistic correction
! - no factor of 4*Pi
! - no DWF
!
function wekoscar1(s)
!
  implicit none
!
  real*8, intent(in) :: s
!  
  real*8 :: wekoscar1
!  
  real*4 ::	sf, fr, wekoa(2), wekob(6) ! coefficients of the weko parameterization
  common /weko1/ wekoa, wekob
  real*4, external :: weko
!  
  sf = real( s, kind=4 )
  fr = weko(wekoa, wekob, sf) 
  wekoscar1 = dble(fr)
!
  return
!
end function wekoscar1
! ******************************************************************** !







! ******************************************************************** !
!
function ri1(bi,bj,g)
!
! erstes integral fuer die absorptionspotentiale
!
  implicit none
  real*4, parameter :: pi = 3.1415927
  real*4, parameter :: c  = 0.5772157
  real*4, intent(in) :: bi,bj, g
  real*4 :: big2, bjg2
  real*4 :: x1, x2, x3, g2
  real*4, external :: ei, rih1
  real*4 :: ri1
  if (g == 0.) then
	ri1 = bi * log( (bi+bj)/bi ) + bj * log( (bi+bj)/bj )
	ri1 = ri1 * pi
    return
  end if
  g2    = g*g
  big2  = bi*g2
  bjg2  = bj*g2
  ri1   = 2.*c + log( big2 ) + log( bjg2 ) - 2.*ei( -bi*bj*g2/(bi+bj) )
  x1    = big2
  x2    = big2*bi/(bi+bj)
  x3    = big2
  ri1   = ri1 + rih1(x1,x2,x3)
  x1    = bjg2
  x2    = bjg2*bj/(bi+bj)
  x3    = bjg2
  ri1   = ri1 + rih1(x1,x2,x3)
  ri1   = ri1 * pi / g2
  return
end function ri1
! ******************************************************************** !



! ******************************************************************** !
!
function ri2(bi,bj,g,u)
!
! zweites integral fuer die absorptionspotentiale
!
  implicit none
  real*4, parameter :: pi = 3.1415927
  real*4, intent(in) :: bi,bj, g, u
  real*4 :: u2, g2, biuh, bjuh, biu, bju
  real*4 :: x1, x2, x3
  real*4, external :: ei, rih1
  real*4 :: ri2
  u2 = u*u
  if (g==0.) then
	ri2 = (bi+u2) * log( (bi+bj+u2)/(bi+u2) )
	ri2 = ri2 + bj * log( (bi+bj+u2)/(bj+u2) )
	ri2 = ri2 + u2 * log( u2/(bj+u2) )
	ri2 = ri2 * pi
	return
  end if
  
  g2    = g*g
  biuh  = bi + .5*u2
  bjuh  = bj + .5*u2
  biu   = bi + u2
  bju   = bj + u2

  ri2   = ei( -.5*u2*g2*biuh/biu ) + ei( -.5*u2*g2*bjuh/bju )
  ri2   = ri2 - ei( -biuh*bjuh*g2/(biuh+bjuh) ) - ei( -.25*u2*g2 )
  ri2   = 2.*ri2
  x1    = .5*u2*g2
  x2    = .25*u2*g2
  x3    = .25*u2*u2*g2/biu
  ri2   = ri2 + rih1(x1,x2,x3)

  x1    = .5*u2*g2
  x2    = .25*u2*g2
  x3    = .25*u2*u2*g2/bju
  ri2   = ri2 + rih1(x1,x2,x3)

  x1    = biuh*g2
  x2    = biuh*biuh*g2/(biuh+bjuh)
  x3    = biuh*biuh*g2/biu
  ri2   = ri2 + rih1(x1,x2,x3)

  x1    = bjuh*g2
  x2    = bjuh*bjuh*g2/(biuh+bjuh)
  x3    = bjuh*bjuh*g2/bju
  ri2   = ri2 + rih1(x1,x2,x3)

  ri2   = ri2 * pi / g2

  return
end function ri2
! ******************************************************************** !




! ******************************************************************** !
!
function rih1(x1,x2,x3)
!
! wertet den ausdruck exp(-x1) * ( ei(x2)-ei(x3) ) aus
!
  implicit none
  real*4, intent(in) :: x1, x2, x3
  real*4, external :: ei, rih2
  real*4 :: rih1
  if (x2 <= 20.  .and.  x3 <= 20.) then
	rih1 = exp(-x1) * ( ei(x2)-ei(x3) )
    return
  end if

  if (x2 > 20) then
 	rih1 = exp(x2-x1)*rih2(x2)/x2
  else 
    rih1 = exp(-x1)*ei(x2)
  end if

  if (x3 > 20) then
    rih1 = rih1 - exp(x3-x1)*rih2(x3)/x3
  else 
    rih1 = rih1 - exp(-x1)*ei(x3)
  end if

  return
end
! ******************************************************************** !




! ******************************************************************** !
!
function rih2(x)
!
!	wertet x*exp(-x)*ei(x) aus fuer grosse x
!	durch interpolation der tabelle ... aus abramowitz
!
  implicit none
  real*4, intent(in) :: x
  real*4 :: rih2
  real*4 :: f(0:20)
  data f / 1.000000,1.005051,1.010206,1.015472,1.020852, &
     &     1.026355,1.031985,1.037751,1.043662,1.049726, &
     &     1.055956,1.062364,1.068965,1.075780,1.082830, &
     &     1.090140,1.097737,1.105647,1.113894,1.122497, 1.131470 /
  real*4 :: x1
  integer*4 :: i, i1
  x1 = 1./x
  i  = int( 200.*x1 )
  i1 = i+1
  rih2 = f(i) + 200.*( f(i1)-f(i) ) * ( x1-.5e-2*real(i) )
  !rih2 = f(i) + 200.*( f(i1)-f(i) ) * ( x1-.5e-3*real(i) ) ! original by Weickenmeier, corrected by JB (26.08.2017)

  return
end function rih2
! ******************************************************************** !



! ******************************************************************** !
!
function ei(x)
!
!	exponentialintegral
!     getestet -60 < x < 60
  implicit none
  real*4, parameter :: a1=8.57332, a2=18.05901, a3=8.63476 , a4=.26777
  real*4, parameter :: b1=9.57332, b2=25.63295, b3=21.09965, b4=3.95849
  real*4, intent(in) :: x
  real*4 :: ei
  real*4 :: xp, ri, ri1, si, summ
  integer*4 :: i
!
  if (x > 60.) then 
    write (6,*) '>>> ei fuer x= ',x,' nicht getestet <<<'
	stop
  end if
  if (x < -60.) then
    ei = 0.
    return
  end if
  if (x < -1.) then
!   abramowitz (5.1.56)
    xp = abs(x)
	ei = -( a4+xp*(a3+xp*(a2+xp*(a1+xp))) ) / &
     &      ( b4+xp*(b3+xp*(b2+xp*(b1+xp))) ) * exp(-xp)/xp
  return
  end if
  ei   = .577216 + log( abs(x) )
  i    = 1
  si   = x
  summ = si

  do
    ri   = real(i)
	ri1  = ri + 1.
	si   = si * x * ri/(ri1*ri1)
	summ = summ + si
	if (abs(si/x) .gt. 1.e-6) then
	   i = i+1	   
	   ! continue with loop
	else
	   ei = ei + summ
	   return
	end if
  end do

  return
end
! ******************************************************************** !



! ******************************************************************** !
!
! modified 2017-06-30 by J. Barthel, Forschungzentrum Juelich GmbH
! - added H (Z=1, hydrogen) and L (Z=0, vacancy)
!
subroutine getweko(z, a, b)
!
! returns coefficients for the weickenmeier & kohl parameterization
! of atomic form factors for electron scattering
!
! neutral atoms of z=0 ... 98
!
  implicit none
  integer*4, intent(in) :: z
  real*4, intent(out) :: a(2), b(6)
!
  integer*4 :: i, j
  real*4 ::	v(0:98)
  real*4 :: bb(1:6,0:98)
!
  data v / 0.5, &
     &     0.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, &
     &     0.5, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.3, &
     &     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, &
     &     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.3, 0.5, 0.5, &
     &     0.5, 0.5, 0.5, 0.4, 0.5, 0.5, 0.5, 0.3, 0.4, 0.6, &
     &     0.6, 0.6, 0.4, 0.4, 0.1, 0.1, 0.3, 0.3, 0.2, 0.2, &
     &     0.2, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.1, &
     &     0.1, 0.1, 0.4, 0.2, 0.5, 0.4, 0.5, 0.5, 0.4, 0.4, &
     &     0.4, 0.3, 0.4, 0.4, 0.4, 0.4, 0.1, 0.2, 0.2, 0.3, &
     &     0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2, 0.3 /
!
	data ((bb(i,j),i=1,6),j=0,0) / &
     &     48.75740,  4.96588, 18.24440, 18.24440, 18.24440, 18.24440/
!
	data ((bb(i,j),i=1,6),j=1,10) / &
     &     48.75740,  4.96588, 18.24440, 18.24440, 18.24440, 18.24440, &
     &      2.54216,  8.74302, 12.69098,  0.43711,  5.29446, 28.25045, &
     &      0.68454,  3.06497,  6.23974,126.17816,131.20160,131.76538, &
     &      0.53996,  3.38752, 55.62340, 50.78098, 67.00502, 96.36635, &
     &      0.33138,  2.97485, 34.01118, 35.98365, 36.68364, 60.80991, &
     &      0.29458,  3.93381, 24.97836, 25.27916, 25.46696, 46.70328, &
     &      0.23925,  4.93515, 18.11895, 15.69698, 15.81922, 40.24150, &
     &      6.37582,  8.03744, 27.20649,  0.11157,  0.38686, 10.89944, &
     &      0.21800,  6.76987,  7.05056,  6.67484, 12.38148, 28.08398, &
     &      0.20055,  5.49814,  6.28052,  7.19211,  7.54763, 23.26388/
!
	data ((bb(i,j),i=1,6),j=11,20) / &
     &      0.21902,  5.30022,  5.31938,  5.28281,  5.28546,128.18391, &
     &      1.97633,  2.80902, 16.39184,  0.05494,  2.06121,121.70512, &
     &      2.29692,  2.35822, 24.98576,  0.07462,  0.55953,128.50104, &
     &      1.73656,  3.04329, 30.57191,  0.05070,  0.99181, 86.18340, &
     &      0.17949,  2.63250,  2.67559, 34.57098, 36.77888, 54.06180, &
     &      1.00609,  4.90414, 31.34909,  0.03699,  0.98700, 44.94354, &
     &      0.18464,  1.47963,  5.20989, 24.79470, 32.06184, 39.09933, &
     &      0.20060,  6.53262, 22.72092,  1.20022,  1.27398, 36.25907, &
     &      0.44424,  3.36735, 19.63031,  0.01824, 23.51332,212.86819, &
     &      0.18274,  2.06638, 16.99062, 11.57795, 13.97594,186.10446/

	data ((bb(i,j),i=1,6),j=21,30) / &
     &      0.14245,  1.46588, 15.46955,  4.24287,  9.80399,121.46864, &
     &      0.12782,  1.45591, 12.09738,  4.61747, 11.96791,105.00546, &
     &      0.13126,  1.39923,  8.00762,  7.98129, 13.41408, 95.30811, &
     &      0.12311,  2.38386,  9.92149,  1.64793, 11.00035, 68.45583, &
     &      0.48173,  3.78306,  8.47337,  0.04690,  8.74544, 77.44405, &
     &      0.44704,  6.89364,  6.90335,  0.05691,  3.02647, 70.86599, &
     &      0.10705,  3.63573,  7.55825,  1.27986,  5.14045, 67.16051, &
     &      0.11069,  1.61889,  6.00325,  5.97496,  6.06049, 59.41419, &
     &      0.11293,  1.89077,  5.08503,  5.07335,  5.09928, 46.38955, &
     &      0.10209,  1.73365,  4.78298,  4.80706,  5.64485, 51.21828/
!
	data ((bb(i,j),i=1,6),j=31,40) / &
     &      0.10642,  1.53735,  5.13798,  4.74298,  4.99974, 61.42872, &
     &      0.09583,  1.67715,  4.70275,  2.91198,  7.87009, 64.93623, &
     &      0.09428,  2.21409,  3.95060,  1.52064, 15.81446, 52.41380, &
     &      0.09252,  1.60168,  3.04917,  3.18476, 18.93890, 47.62742, &
     &      0.09246,  1.77298,  3.48134,  1.88354, 22.68630, 40.69434, &
     &      0.49321,  2.08254, 11.41282,  0.03333,  2.09673, 42.38068, &
     &      0.15796,  1.71505,  9.39164,  1.67464, 23.58663,152.53635, &
     &      0.36052,  2.12757, 12.45815,  0.01526,  2.10824,133.17088, &
     &      0.09003,  1.41396,  2.05348, 10.25766, 10.74831, 90.63555, &
     &      0.10094,  1.15419,  2.34669, 10.58145, 10.94962, 82.82259/
!
	data ((bb(i,j),i=1,6),j=41,50) / &
     &      0.09243,  1.16977,  5.93969,  1.30554, 13.43475, 66.37486, &
     &      0.43543,  1.24830,  7.45369,  0.03543,  9.91366, 61.72203, &
     &      0.45943,  1.18155,  8.31728,  0.03226,  8.32296, 64.97874, &
     &      0.08603,  1.39552, 11.69728,  1.39552,  3.45200, 55.55519, &
     &      0.09214,  1.11341,  7.65767,  1.12566,  8.32517, 48.38017, &
     &      0.09005,  1.12460,  9.69801,  1.08539,  5.70912, 33.48585, &
     &      0.08938,  3.19060,  9.10000,  0.80898,  0.81439, 41.34453, &
     &      0.28851,  1.61312,  8.99691,  0.01711,  9.46666, 58.13256, &
     &      0.08948,  1.23258,  8.23129,  1.22390,  7.06201, 59.69622, &
     &      0.07124,  0.85532,  6.40081,  1.33637,  6.38240, 50.92361/
!
	data ((bb(i,j),i=1,6),j=51,60) / &
     &      0.35749,  1.32481,  6.51696,  0.03550,  6.51913, 50.80984, &
     &      0.50089,  3.95301,  7.62830,  0.03005,  0.50737, 49.62628, &
     &      0.08429,  1.12959,  8.86209,  1.12981,  9.13243, 56.01965, &
     &      0.27796,  1.62147, 11.45200,  0.02032,  3.27497, 51.44078, &
     &      0.12045,  1.53654,  9.81569, 41.21656, 42.62216,224.34816, &
     &      0.12230,  1.44909,  9.50159, 49.40860, 74.94942,217.04485, &
     &      0.08930,  1.26225,  8.09703,  1.20293, 17.65554,116.61481, &
     &      0.08504,  1.28286, 11.22123,  1.32741,  4.61040,112.19678, &
     &      0.09805,  1.52628,  8.58953,  1.23893, 22.49126,140.02856, &
     &      0.09413,  1.26616,  5.98844, 17.78775, 18.14397,132.59305/
!
	data ((bb(i,j),i=1,6),j=61,70) / &
     &      0.09447,  1.25111,  5.91205, 16.28675, 16.73089,127.90916, &
     &      0.09061,  1.59281, 10.64077,  1.78861,  2.22148,124.56328, &
     &      0.10485,  1.54396,  8.65223,  7.09290, 53.36537,183.69014, &
     &      0.09338,  1.38681,  7.35883,  1.55122, 20.81916,111.03201, &
     &      0.10190,  1.52368,  7.16923, 20.86269, 49.29465,166.09206, &
     &      0.08402,  1.40890,  7.14042,  1.34848, 11.42203,108.01204, &
     &      0.09441,  1.61807,  6.27142, 40.34946, 42.82722,130.59616, &
     &      0.08211,  1.25106,  4.81241, 10.84493, 10.90164,100.07855, &
     &      0.09662,  1.60236,  5.67480, 30.59014, 31.12732,138.69682, &
     &      0.09493,  1.60220,  5.43916, 28.31076, 29.27660,138.08665/
!
	data ((bb(i,j),i=1,6),j=71,80) / &
     &      0.09658,  1.56751,  5.32170, 34.18217, 35.25187,121.42893, &
     &      0.09294,  1.55499,  5.25121, 37.51883, 38.88302,105.16978, &
     &      0.06298,  0.81950,  2.89124,  5.54290,  5.98101, 54.42459, &
     &      0.07902,  1.37096,  8.23364,  1.38300,  1.39219, 77.11813, &
     &      0.05266,  0.90718,  4.43830,  0.94590,  4.37477, 43.97909, &
     &      0.22700,  1.56975,  6.34451,  0.01564,  1.61769, 46.15815, &
     &      0.05055,  0.86775,  5.09325,  0.88123,  3.56919, 39.77390, &
     &      0.05253,  0.83773,  3.95899,  0.81515,  6.44217, 34.21146, &
     &      0.54927,  1.72752,  6.71952,  0.02637,  0.07253, 35.45745, &
     &      0.21941,  1.41611,  6.68241,  0.01472,  1.57578, 37.15826/
!
	data ((bb(i,j),i=1,6),j=81,90) / &
     &      0.22459,  1.12822,  4.30289,  0.01485,  7.15607, 43.08737, &
     &      0.06432,  1.19406,  7.39342,  1.14160,  1.28905, 51.13401, &
     &      0.05380,  0.86719,  1.87540,  7.64796,  7.86794, 45.63897, &
     &      0.50112,  1.63784,  6.78551,  0.02187,  0.08602, 46.72951, &
     &      0.22321,  1.10827,  3.59116,  0.01011, 11.63732, 45.06839, &
     &      0.21152,  1.14015,  3.41473,  0.01188, 13.41211, 43.11389, &
     &      0.09435,  1.02649,  6.25480, 32.51444, 36.29119,149.11722, &
     &      0.07300,  1.01825,  5.89629,  1.03089, 20.37389,115.34722, &
     &      0.07515,  0.94941,  3.72527, 17.58346, 19.75388,109.12856, &
     &      0.06385,  0.90194,  4.65715,  0.90253, 15.70771, 83.69695/
!
	data ((bb(i,j),i=1,6),j=91,98) / &
     &      0.07557,  0.84920,  4.00991, 16.95003, 17.78767,100.20415, &
     &      0.07142,  1.14907,  9.21231,  0.95923,  1.20275,104.32746, &
     &      0.06918,  0.98102,  5.95437,  0.99086, 22.06437, 90.98156, &
     &      0.07136,  0.95772,  6.13183,  0.97438, 15.67499, 89.86625, &
     &      0.07301,  0.93267,  6.34836,  0.91032, 13.26179, 86.85986, &
     &      0.05778,  0.72273,  3.01146,  9.21882,  9.53410, 65.86810, &
     &      0.07088,  0.77587,  6.14295,  1.79036, 15.12379, 83.56983, &
     &      0.06164,  0.81363,  6.56165,  0.83805,  4.18914, 61.41408/
!
!     for z = 0 a(i) will be zero and v, b(i) won't matter ! 
  a(1)	= 0.023933659 * real(z) / (3.0 * (1.0 + v(z)))
  a(2)	= v(z)    * a(1)
!
  do i = 1, 6
	b(i)	= bb(i, z)
  end do
!
  return
end
! ******************************************************************** !



! ******************************************************************** !
!     added 23.11.2010, J. Barthel, Forschungszentrum Juelich GmbH
!     modified 30.06.2017, J. Barthel, Forschungszentrum Juelich GmbH
!     - added 'L ' for vacancy at z = 0
subroutine wekosymb(z, symbol)
  implicit none
  integer*4, intent(in) :: z
  character*2, intent(out) :: symbol
!
  character*2 :: sy(0:98)
!
  data sy / 'L ', &
     &  'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
     &  'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
     &  'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
     &  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
     &  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
     &  'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
     &  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
     &  'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
     &  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
     &  'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf' /
!
  symbol	= sy(z)
!
  return
end subroutine wekosymb
! ******************************************************************** !
	
	
