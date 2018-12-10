! ******************************************************************** !
!
! Implementation of Waasmaier & Kirfel X-ray atomic form factors
!
! Acta Cryst. A 51 (1995) 416 - 431
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
! wakisca (g,dw,z,akv,dwfflg) : (function) (complex*8) [nm]
!                     returns the atomic form factor for electron
!                     scattering for a given atom type z at spatial
!                     frequency g [1/nm], with Debye-Waller parameter
!                     dw [nm^2] = Biso and electron energy akv [keV]
!                     real-part = atomic form factor [nm]
!                     imag-part = absorptive form factor [nm]
! wakiscar(g,dw,z,akv,dwfflg) : (function) (real*4) [nm]
!                     returns the atomic form factor for electron
!                     scattering for a given atom type z at spatial
!                     frequency g [1/nm], with Debye-Waller parameter
!                     dw [nm^2] = Biso and electron energy akv [keV]
! wakiscar1(s)      : (function) (real*8) [A] (only for internal use)
!                     returns the atomic form factor for electron
!                     scattering for scattering vector magnitude s.
!                     s = g/2 [1/A] (real*8)
!                     no relativistic correction and no DWF applied.
!                     atomic form factor parameters wakia are
!                     defined in common block /waki1/
! wakiabs(g,dw,a)   : (function) (real*4) [A]
!                     returns the absorptive form factor for given
!                     Waasmaier & Kirfel parameters a(14), Debye-
!                     Waller parameter dw [A^2] = Biso and spatial
!                     frequency g [1/A] = scattering vector magnitude
!                     g = 2*s [1/A] (real*4)
! wakif2(tmax,dw,z,k0) : (function) (real*4) [A^2]
!                     Returns the integral of f^2 from theta=0 to
!                     theta=tmax for given
!                     atomic number z,
!                     Debye-Waller parameter dw [A^2] = Biso and
!                     wave number k0 [1/A]
! wakisymb(z, symb) : returns the atom symbol symb for atomic number z
! getwaki(z, a)     : returns the x-ray parameterization coefficients
!                     a(14) for atomic number z
! real*4 :: waki(a, s) : (function) (real*4)
!                     returns the electron atomic form factor fe(s) [A]
!                     for given parameterization coefficients a(14)
!                     and scattering vector magnitude s [1/A]
!                     s = k*arcsin(theta/2) = g/2
!                     doesn't contain the ionic charge term for ions.
!
! First use getwaki to obtain the coefficients of the parameterization.
! Then use waki to calculate fe(s).
! Don't forget to add the ionic charge potential later.
!
! ******************************************************************** !



! ******************************************************************** !
!
! wakiscar
!
! returns the atomic form factor for electron scattering including
! an absorptive form factor on demand
function wakisca(g,dw,z,akv,dwfflg,absflg)
!--- input
!-   - real     g       spatial frequency [1/nm]
!-   - real     dw      debye-waller parameter in nm^2
!-   - integer  z       atomic number
!-   - akv      akv     high tension in kV
!-   - logical  dwfflg  flag: use debye-waller factor
!-   - logical  absflg  flag: use debye-waller factor
!-
!--- output / return value
!-   - real*4           atomic form factor for the given
!-                      atom type at the given spatial frequency
!-
  implicit none
! parameters
  real*4, parameter	:: twopi    = 6.283185307180 ! 2*pi
  real*4, parameter :: fourpi   = 12.56637061436 ! 4*pi
  real*4, parameter :: e0kv     = 510.998947 ! electron rest mase energy equivalent m0*c^2 in [keV]
  real*4, parameter :: thrima   = 1.0E-8 ! threshold for imaginary part calculation
  real*4, parameter :: kpre     = 0.506773075855 ! 2*pi * e / (h*c) *10^-7 [1/(Angström * kV)]
! ínterface variable
  real*4, intent(in) :: g, dw, akv
  integer*4, intent(in) :: z
  logical, intent(in) :: dwfflg, absflg
  complex*8 :: wakisca
  real*4, external :: waki, wakiabs
! local variables
  real*4 :: dwa, sa, dewa
  real*4 :: k0, rc
  real*4 :: fr, fi
  real*4 :: a2(14)
!
  wakisca = 0.0
! relativistic correction
  rc    = (e0kv + akv) / e0kv
!
! Debye-Waller Temperature Coeff dwa = 8 pi**2 usa**2
!   usa = rms amplitude of atom displacement in [ang]
!   [nm] to [ang]
!
  dewa  = 1.0 ! preset DWF to 1
  dwa   = 100.0 * dw ! translate DW-parameter to angstrom squared
  sa    = 0.050 * g ! sa is ga / 4pi = g / 20
!
! get the Waasmaier & Kirfel table parameters for atom z
  call getwaki(z,a2)
!
! calculate the atomic form factor based on the parameterization
  fr    = rc * waki(a2, sa) * fourpi ! fr is in [A] and relativistically corrected now
!
  if (dwfflg) then ! use DWF
    dewa = exp( -dwa * sa*sa ) ! calculate the DWF as exp( -B * s^2 ) = exp( -0.25 * B * g^2 )
    fr = fr * dewa
  end if
!
! calculate imag form factor component
  if ((dwfflg).and.(absflg).and.(abs(dwa)>=thrima).and.(abs(fr)>=thrima)) then
!   wave number in [A^-1] 
    k0  = kpre * sqrt((2.0*e0kv + akv) * akv) ! 2*Pi*k ! andere k-Notation hier: Aufpassen!
	fi  = rc * rc * wakiabs(0.1*g, dwa, a2, k0/twopi) * k0 ! apply double relativistic correction and divide by (2*pi*k)
    !if (dwfflg) fi = fi * dewa ! apply the Debye Waller factor to the absorptive part
  end if
!
!-------[ang] to [nm]
!
  wakisca = 0.1 * cmplx (fr, fi)
!
  return
!
end function wakisca
!
! ******************************************************************** !


! ******************************************************************** !
!
! wakiscar
!
! returns the atomic form factor for electron scattering
function wakiscar(g,dw,z,akv,dwfflg)
!--- input
!-   - real     g       spatial frequency [1/nm]
!-   - real     dw      debye-waller parameter in nm^2
!-   - integer  z       atomic number
!-   - akv      akv     high tension in kV
!-   - logical  dwflg   flag: use debye-waller factor
!-
!--- output / return value
!-   - real*4           atomic form factor for the given
!-                      atom type at the given spatial frequency
!-
  implicit none
! parameters
  real*4, parameter	:: twopi = 6.28318530718 ! 2*pi
  real*4, parameter :: fourpi = 12.56637061436 ! 4*pi
  real*4, parameter :: r8pi2 = 0.0126651479553 ! 1/(8*pi^2)
  real*4, parameter :: e0kv = 510.998947 ! electron rest mase energy equivalent m0*c^2 in [keV]
! ínterface variable
  real*4, intent(in) :: g, dw, akv
  integer*4, intent(in) :: z
  logical, intent(in) :: dwfflg
  real*4 :: wakiscar
  real*4, external :: waki
! local variables
  real*4 :: dwa, sa, dewa
  real*4 :: rc
  real*4 :: fr
  real*4 :: a2(14)
!
  wakiscar = 0.0
! relativistic correction
  rc    = (e0kv + akv) / e0kv
!
! Debye-Waller Temperature Coeff dwa = 8 pi**2 usa**2
!   usa = rms amplitude of atom displacement in [ang]
!   [nm] to [ang]
!
  dewa  = 1.0 ! preset DWF to 1
  dwa   = 100.0 * dw ! translate DW-parameter to angstrom squared
  sa    = 0.050 * g ! sa is ga / 4pi = g / 20
!
! get the Waasmaier & Kirfel table parameters for atom z
  call getwaki(z,a2)
!
! calculate the atomic form factor based on the parameterization
  fr    = rc * waki (a2, sa) * fourpi ! fr is in [A] and relativistically corrected now
!
  if (dwfflg) then ! use DWF
    dewa = exp( -dwa * sa*sa) ! calculate the DWF as exp( -B * s^2 ) = exp( -0.25 * B * g^2 )
    fr = fr * dewa
  end if
!
!-------[ang] to [nm]
!
  wakiscar = 0.1 * fr
!
  return
!
end function wakiscar
!
! ******************************************************************** !






! ******************************************************************** !
!
! wakiabs
!
! This function starts a numerical integration of the absorptive
! form factor for electron scattering form factors given in the
! Waasmaier and Kirfel parameterization.
!
! input: (real*4)
!   ga = diffraction vector magnitude [1/A]
!   dwa = Debye-Waller parameter (Biso) [A^2]
!   a(14) = parameters from Wa&Ki's table (defines the atom)
!   k0 = incident beam vacuum wave vector magnitude [1/A] = 1/Lambda
!      = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
! output: (real*4)
!   wakiabs = absorptive form factor without the required pre-factor
!             DWF * gamma^2 / (2*pi*k0)
!
function wakiabs(ga,dwa,a,k0)
!
  implicit none
!
  real*4, intent(in) :: ga, dwa, a(14), k0
!
  real*4 :: wakiabs
!
  real*4 :: wakia(14)
  real*8 :: si0, si1, fa, pi0, pi1
  real*8 :: wakig, wakib, wakik
! common blocks ...
  common /waki1/ wakia        !  for wakiscar1
  common /waki2/ wakig, wakib, wakik !  for integrand wakimug
!
  real*8, external :: dsgrid2d ! link integration.f90 !
  real*8, external :: wakimug ! integrand functions

!
  wakiabs = 0.0
  wakia = a ! set the Waasmaier and Kirfel parameters to be used
  wakig = dble(ga) ! store sa for the integrand wakimug
  wakib = dble( dwa ) ! store Debye-Waller parameter for both integrands
  wakik = k0
! set integration range
  si0 = 0.d+0 ! 0
  si1 = 3.1415926535898D+0 ! Pi
  pi0 = 0.d+0 ! 0
  pi1 = 3.1415926535898D+0 ! Pi
! call the numerical integrator
! use a grid of 128 pixels along theta with a square sampling
!           and  64 pixels along phi with linear sampling
!           phi on the half side only (symmteric) therefore:  * 2
  fa = dsgrid2d(wakimug,si0,si1,pi0,pi1,2.0d+0,1.0d+0,128,64) * 2.0D+0
  !
  wakiabs = real(fa, kind=4 )
  return
end function wakiabs
!
! ******************************************************************** !




! ******************************************************************** !
!
! wakif2
!
! This function starts a numerical integration of the squared form
! factor for electron scattering given in the
! Waasmaier and Kirfel parameterization.
!
! input: (real*4)
!   tmax = max. scattering angle [rad] (upper limit is Pi)
!   dwa = Debye-Waller parameter (Biso) [A^2]
!   z = atomic number
!   k0 = incident beam vacuum wave vector magnitude [1/A] = 1/Lambda
!      = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
! output: (real*8)
!   wakif2 = integrated squared form factor
!
function wakif2(tmax,dwa,z,k0)
!
  implicit none
!
  integer*4, parameter :: ngrid = 2048
  real*8, parameter :: twopi = 6.2831853071796D+0
!
  real*4, intent(in) :: tmax, dwa, k0
  integer*4, intent(in) :: z
!
  real*8 :: wakif2
!
  real*4 :: wakia(14)
  real*8 :: si0, si1, fa
  real*8 :: wakig, wakib, wakik
! common blocks ...
  common /waki1/ wakia        !  for wakiscar1
  common /waki2/ wakig, wakib, wakik !  for integrand wakimug
!
  real*8, external :: dsgrid1d ! link integration.f90 !
  real*8, external :: wakif2g ! integrand functions

!
  wakif2 = 0.0
  call getwaki(z, wakia) ! set the Waasmaier and Kirfel parameters to be used
  wakig = 0.D+0 ! not used
  wakib = dble( dwa ) ! store Debye-Waller parameter
  wakik = k0 ! store wave number [1/A]
! set integration range
  si0 = 0.d+0 ! 0
  si1 = min(3.1415926535898D+0, dble(abs(tmax)) ) ! min(Pi,abs(tmax))
! call the numerical integrator
! - f2 is constant with phi, therefore:  * 2*Pi
  fa = dsgrid1d(wakif2g,si0,si1,2.0d+0,ngrid) * twopi
  !
  wakif2 = fa
  return
end function wakif2
!
! ******************************************************************** !




!**********************************************************************!
!
! function wakimug(theta,phi)
!
! returns the integrand for the absorptive form factor at g=wakig
! for debye-waller parameter b=wakib
! for incident electron wave vector k=wakik
! for atom defined by ai = wakia
!
! Set the values of common block parameters before using this function!
! uses common blocks /waki1/ and /waki2/
!
function wakimug(theta,phi)
  
  implicit none
  
  real*8, intent(in) :: theta, phi
  real*8 :: wakimug
  real*8 :: k, q, g, qg, biso, ct, st, cp, twok, sg, sq, sqg
  real*8 :: fq, fqg, dwfg, dwfq, dwfqg
  real*8 :: wakig, wakib, wakik
  common /waki2/ wakig, wakib, wakik
  
  real*8, external :: wakiscar1 ! form factor function
                                ! uses common block /waki1/
  
  ct = dcos(theta)
  st = dsin(theta)
  cp = dcos(phi)
  ! q = scattering vector length of Q in the ewald sphere
  ! Q = (qx, qy, qz) = K' - K
  ! qx = k * sin(theta) * cos(phi)
  ! qy = k * sin(theta) * sin(phi)
  ! qz = k * ( cos(theta) - 1 )
  ! k = |K| = wakik
  k = wakik
  twok = k+k
  ! q = Sqrt[ 2 k^2 (1 - Cos[q]) ]
  q = k * dsqrt( 2.D+0 - 2.D+0 * ct )
  ! g = length of some reciprocal space vector G = (gx, gy, gz)
  !     with gx = g, gy = 0, gz = 0 (obda)
  g = wakig
  ! qg = length of the difference vector Q - G
  qg = dsqrt( g*g + twok*k - twok*(k*ct + g*cp*st) )
  biso = wakib
  sg = 0.5D+0 * g ! g/2
  sq = 0.5D+0 * q ! q/2
  sqg = 0.5D+0 * qg ! qg/2
  fq = wakiscar1(sq) ! f(s) = f(g/2)
  fqg = wakiscar1(sqg) ! f(s) = f(g/2)
  dwfg = dexp( -biso*sg*sg ) ! DWF = EXP( -1/4 * BISO * g^2 ) = EXP( -BISO * s^2 )
  dwfq = dexp( -biso*sq*sq )
  dwfqg = dexp( -biso*sqg*sqg )
  wakimug = fq*fqg * (dwfg - dwfq*dwfqg) * st
  
  return
    
end function wakimug



!**********************************************************************!
!
! function wakif2g(theta)
!
! returns the integrand f(theta)^2 * sin(theta)
! for debye-waller parameter b=wakib (set to zero if not used)
! for incident electron wave vector k=wakik
!     -> s = q/2 = k*sin(theta/2)
!                = k*sqrt( (1-cos(theta))/2 )
! for atom defined by ai = wakia
!
! Set the values of common block parameters before using this function!
! uses common blocks /waki1/ and /waki2/
!
function wakif2g(theta)
  
  implicit none
  
  real*8, intent(in) :: theta
  real*8 :: wakif2g
  real*8 :: k, q, biso, ct, st, sq
  real*8 :: fq
  real*8 :: wakig, wakib, wakik
  common /waki2/ wakig, wakib, wakik
  
  real*8, external :: wakiscar1 ! form factor function
                                ! uses common block /waki1/
  
  ct = dcos(theta)
  st = dsin(theta)
  ! q = scattering vector length of Q in the ewald sphere
  !     for scattering angle theta
  ! Q = (qx, qy, qz) = K' - K
  ! qx = k * sin(theta) * cos(phi)
  ! qy = k * sin(theta) * sin(phi)
  ! qz = k * ( cos(theta) - 1 )
  ! k = |K| = wakik
  k = wakik
  ! q = Sqrt[ 2 k^2 (1 - Cos[q]) ]
  q = k * dsqrt( 2.D+0 - 2.D+0 * ct )
  biso = wakib
  sq = 0.5D+0 * q ! q/2
  fq = wakiscar1(sq) ! f(s) = f(q/2)
  if (dabs(biso)>0.D+0) then
    fq = fq * dexp( -biso*sq*sq ) ! DWF = EXP( -1/4 * BISO * q^2 ) = EXP( -BISO * s^2 )
  end if
  wakif2g = st*fq*fq ! Sin(theta)*(f(theta)*DWF))^2
  
  return
    
end function wakif2g



! ******************************************************************** !
!
! this function is used internally only for the absorptive
! potential integration, the interface is via common block /waki1/+
!
! assumes that - wakia(14) is set via getwaki
!
! differences to wakiscar:
! - no relativistic correction
! - no factor of 4*pi
! - no debye-waller factor
!
function wakiscar1(s)
  implicit none
  real*8, parameter :: alph2 = 4.0d-04 ! square of inverse yukawa range in A (only used for small s2 in ions) 
  real*8, parameter :: c1 = 0.0239336609804d+0 ! m0 * e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) [A]
  real*8, intent(in) :: s
  real*4 :: sf
  real*4, external :: waki
  real*8 :: wakiscar1
  real*8 :: fa, faio
  real*4 :: wakia(14), dz
  common /waki1/ wakia
  dz = 0.0
  dz = wakia(13) - wakia(12) - sum( wakia(2:6) )
  dz = real( nint(dz) ) ! support integer ionic charges only
  fa = 0.0d+0
  faio = 0.0d+0
  sf = real(s, kind=4 )
  fa = dble( waki(wakia,sf) ) ! core form factor
  if (abs(dz)>0.0) then ! add ionic charge potential
    faio = c1 * dble(dz) / (s*s+alph2)
  end if
  wakiscar1 = fa + faio
end function wakiscar1


! ******************************************************************** !
! 
! wakisymb
!
! returns the atom symbol for a given z
subroutine wakisymb(z, symbol)
  implicit none
  integer*4, intent(in) :: z
  character*2, intent(out) :: symbol
!
  character*2 :: sy(0:92)
  character(len=10) :: snum
!
  data sy / &
     &  'L ', &   
     &  'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
     &  'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
     &  'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
     &  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
     &  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
     &  'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
     &  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
     &  'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
     &  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
     &  'Pa', 'U '/
!
  if (z<0 .or. z>92) then
    write(unit=snum,fmt=*) z
    write(*,*) "Atomic number Z="//trim(adjustl(snum))// &
             & " is not in the table of Waasmaier and Kirfel."
    stop
  end if
!
  symbol = sy(z)
!
  return
end subroutine wakisymb
	

! ******************************************************************** !	
!
!     fx parameter table of Waasmaier and Kirfel (X-ray form factors)
!     copied from muSTEM
!     original table from: Acta Cryst. A 51 (1995) 416 - 431
!
!     fx(s) = Sum[ a_i * Exp[ -b_i * s^2 ], (i,1,5) ] + c
!             for s = 0 ... 6 1/A
!     fx(s) = d_1 / (d_2 + s^2)
!             for s > 6 1/A
!
!     by J. Barthel, Forschungszentrum Jülich GmbH, Jülich, Germany
!                    RWTH Aachen University, Aachen, Germany
!     2017-July-07
!
!     Interface:
!            getwaki(z,a)
!     Input: 
!            integer z = atomic number
!     Output:
!            real a(14) = fx parameters of Waasmaier and Kirfel 
!                         + large s extension parameters
!
!     Returns the coefficients of the extended table as a list a
!     with in total 14 entries.
!     Entry key:
!     a(1) = Z (atomic number)
!     a(2) = a1 (parameterization coefficient)
!     ...
!     a(6) = a5 (parameterization coefficient)
!     a(7) = b1 (parameterization coefficient)
!     ...
!     a(11) = b5 (parameterization coefficient)
!     a(12) = c  (parameterization coefficient)
!     a(13) = d_1 ( = Z, parameterization coefficient)
!     a(14) = d_2 (parameterization coefficient)
!
! ---------------------------------------------------------------------

subroutine getwaki(z, a)
  implicit none
  integer*4, intent(in) :: z
  real*4, intent(out) :: a(14)
!
  integer*4 :: i, j
  real*4 :: bb(1:14,0:92)
  character(len=10) :: snum
!
! 0 -> vacancy z = 0 -> b's won't matter, all a's set to zero
  data ((bb(i,j),i=1,14),j=0,0) / &
  &   0.0,      0.0,      0.0,      0.0,      0.0,      0.0,11.553918,4.59583100,1.54629900,26.4639640,0.37752300,       0.0,0.0,8.769017454422158E-3/
!
! 1 - 10
  data ((bb(i,j),i=1,14),j=1,10) / &
  &  1.0,0.73235400,0.7538960,0.2838190,0.1900030,0.0391390,11.553918,4.59583100,1.54629900,26.4639640,0.37752300,0.00048700,2.0,8.769017454422158E-3,&
  &  2.0,0.73235400,0.7538960,0.2838190,0.1900030,0.0391390,11.553918,4.59583100,1.54629900,26.4639640,0.37752300,0.00048700,2.0,8.769017454422158E-3,&
  &  3.0,0.97463700,0.1584720,0.8118550,0.2624160,0.7901080,4.3349460,0.34245100,97.1029690,201.363824,1.40923400,0.00254200,3.0,3.053829815542358E-2,&
  &  4.0,1.53371200,0.6382830,0.6010520,0.1061390,1.1184140,42.662078,0.59542000,99.1065010,0.15134000,1.84309300,0.00251100,4.0,2.673011412693727E-2,&
  &  5.0,2.08518500,1.0645800,1.0627880,0.1405150,0.6417840,23.494069,1.13789400,61.2389750,0.11488600,0.39903600,0.00382300,5.0,4.375691177284862E-2,&
  &  6.0,2.65750600,1.0780790,1.4909090,-4.241070,0.7137910,14.780758,0.77677500,42.0868430,-0.0002940,0.23953500,4.29798300,6.0,7.163539380752063E-2,&
  &  7.0,11.8937800,3.2774790,1.8580920,0.8589270,0.9129850,0.0001580,10.2327230,30.3446900,0.65606500,0.21728700,-11.804902,7.0,0.112381925459201,&
  &  8.0,2.96042700,2.5088180,0.6378530,0.7228380,1.1427560,14.182259,5.93685800,0.11272600,34.9584810,0.39024000,0.02701400,8.0,0.171987724092682,&
  &  9.0,3.51195400,2.7722440,0.6783850,0.9151590,1.0892610,10.687859,4.38046600,0.09398200,27.2552030,0.31306600,0.03255700,9.0,0.223745513843300,&
  & 10.0,4.1837490,2.9057260,0.5205130,1.1356410,1.2280650,8.1754570,3.25253600,0.06329500,21.8139090,0.22495200,0.02557600,10.0,0.287627317325155/
!
! 11 - 20
  data ((bb(i,j),i=1,14),j=11,20) / &
  & 11.0,4.9101270,3.0817830,1.2620670,1.0989380,0.5609910,3.2814340,9.11917800,0.10276300,132.013942,0.40587800,0.07971200,11.0,0.366747421044202,&
  & 12.0,4.7089710,1.1948140,1.5581570,1.1704130,3.2394030,4.8752070,108.506079,0.11151600,48.2924070,1.92817100,0.12684200,12.0,0.470985327996563,&
  & 13.0,4.7307960,2.3139510,1.5419800,1.1175640,3.1547540,3.6289310,43.0511660,0.09596000,108.932389,1.55591800,0.13950900,13.0,0.528931702148052,&
  & 14.0,5.2753290,3.1910380,1.5115140,1.3568490,2.5191140,2.6313380,33.7307280,0.08111900,86.2886400,1.17008700,0.14507300,14.0,0.592196010700308,&
  & 15.0,1.9505410,4.1469300,1.4945600,1.5220420,5.7297110,0.9081390,27.0449530,0.07128000,67.5201900,1.98117300,0.15523300,15.0,0.660050129218602,&
  & 16.0,6.3721570,5.1545680,1.4737320,1.6350730,1.2093720,1.5143470,22.0925280,0.06137300,55.4451760,0.64692500,0.15472200,16.0,0.726458197799532,&
  & 17.0,1.4460710,6.8706090,6.1518010,1.7503470,0.6341680,0.0523570,1.19316500,18.3434160,46.3983940,0.40100500,0.14677300,17.0,0.792912007492641,&
  & 18.0,7.1880040,6.6384540,0.4541800,1.9295930,1.5236540,0.9562210,15.3398770,15.3398620,39.0438240,0.06240900,0.26595400,18.0,0.874904035952563,&
  & 19.0,8.1639910,7.1469450,1.0701400,0.8773160,1.4864340,12.816323,0.80894500,210.327009,39.5976510,0.05282100,0.25361400,19.0,0.924257331263367,&
  & 20.0,8.5936550,1.4773240,1.4362540,1.1828390,7.1132580,10.460644,0.04189100,81.3903820,169.847839,0.68809800,0.19625500,20.0,0.967132503306553/
!
! 21 - 30
  data ((bb(i,j),i=1,14),j=21,30) / &
  & 21.0,1.4765660,1.4872780,1.6001870,9.1774630,7.0997500,53.131022,0.03532500,137.319495,9.09803100,0.60210200,0.15776500,21.0,1.012975070461970,&
  & 22.0,9.8185240,1.5226460,1.7031010,1.7687740,7.0825550,8.0018790,0.02976300,39.8854230,120.158000,0.53240500,0.10247300,22.0,1.050876425199770,&
  & 23.0,10.473575,1.5478810,1.9863810,1.8656160,7.0562500,7.0819400,0.02604000,31.9096720,108.022844,0.47488200,0.06774400,23.0,1.086704198113510,&
  & 24.0,11.007069,1.5554770,2.9852930,1.3478550,7.0347790,6.3662810,0.02398700,23.2448380,105.774500,0.42936900,0.06551000,24.0,1.115649280845480,&
  & 25.0,11.709542,1.7334140,2.6731410,2.0233680,7.0031800,5.5971200,0.01780000,21.7884190,89.5179150,0.38305400,-0.1472930,25.0,1.137904401534530,&
  & 26.0,12.311098,1.8766230,3.0661770,2.0704510,6.9751850,5.0094150,0.01446100,18.7430410,82.7678740,0.34650600,-0.3049310,26.0,1.157781969494130,&
  & 27.0,12.914510,2.4819080,3.4668940,2.1063510,6.9608920,4.5071380,0.00912600,16.4381300,76.9873170,0.31441800,-0.9365720,27.0,1.170782166833100,&
  & 28.0,13.521865,6.9472850,3.8660280,2.1359000,4.2847310,4.0772770,0.28676300,14.6226340,71.9660780,0.00443700,-2.7626970,28.0,1.181463309062120,&
  & 29.0,14.014192,4.7845770,5.0568060,1.4579710,6.9329960,3.7382800,0.00374400,13.0349820,72.5547930,0.26566600,-3.2544770,29.0,1.189120722007890,&
  & 30.0,14.741002,6.9077480,4.6423370,2.1917660,38.424042,3.3882320,0.24331500,11.9036890,63.3121300,0.00039700,-36.915828,30.0,1.195267016288100/
!
! 31 - 40
  data ((bb(i,j),i=1,14),j=31,40) / &
  & 31.0,15.758946,6.8411230,4.1210160,2.7146810,2.3952460,3.1217540,0.22605700,12.4821960,66.2036220,0.00723800,-0.8473950,31.0,1.200512778463350,&
  & 32.0,16.540614,1.5679000,3.7278290,3.3450980,6.7850790,2.8666180,0.01219800,13.4321630,58.8660460,0.21097400,0.01872600,32.0,1.200655847244250,&
  & 33.0,17.025643,4.5034410,3.7159040,3.9372000,6.7901750,2.5977390,0.00301200,14.2721190,50.4379970,0.19301500,-2.9841170,33.0,1.198314141346290,&
  & 34.0,17.354071,4.6532480,4.2594890,4.1364550,6.7491630,2.3497870,0.00255000,15.5794600,45.1812010,0.17743200,-3.1609820,34.0,1.198528926857700,&
  & 35.0,17.550570,5.4118820,3.9371800,3.8806450,6.7077930,2.1192260,16.5571850,0.00248100,42.1640090,0.16212100,-2.4920880,35.0,1.199161541836410,&
  & 36.0,17.655279,6.8481050,4.1710040,3.4467600,6.6852000,1.9082310,16.6062350,0.00159800,39.9174710,0.14689600,-2.8105920,36.0,1.199687791703070,&
  & 37.0,8.1231340,2.1380420,6.7617020,1.1560510,17.679547,15.142385,33.5426660,0.12937200,224.132506,1.71336800,1.13954800,37.0,1.210575215851800,&
  & 38.0,17.730219,9.7958670,6.0997630,2.6200250,0.6000530,1.5630600,14.3108680,0.12057400,135.771318,0.12057400,1.14025100,38.0,1.201749013781440,&
  & 39.0,17.792040,10.253252,5.7149490,3.1705160,0.9182510,1.4296910,13.1328160,0.11217300,108.197029,0.11217300,1.31787000,39.0,1.375021858691620,&
  & 40.0,17.859771,10.911038,5.8211150,3.5125130,0.7469650,1.3106920,12.3192850,0.10435300,91.7775440,0.10435300,1.12485900,40.0,1.188442917505480/
!
! 41 - 50
  data ((bb(i,j),i=1,14),j=41,50) / &
  & 41.0,17.958398,12.063054,5.0070150,3.2876670,1.5310190,1.2115900,12.2466870,0.09861500,75.0119440,0.09861500,1.12345200,41.0,1.189359603719730,&
  & 42.0,6.2362180,17.987711,12.973127,3.4514260,0.2108990,0.0907800,1.10831000,11.4687200,66.6841530,0.09078000,1.10877000,42.0,1.199482060769910,&
  & 43.0,17.840964,3.4282360,1.3730120,12.947364,6.3354690,1.0057290,41.9013830,119.320541,9.78154200,0.08339100,1.07478400,43.0,1.202197031328830,&
  & 44.0,6.2716240,17.906739,14.123269,3.7460080,0.9082350,0.0770400,0.92822200,9.55534500,35.8606780,123.552247,1.04399200,44.0,1.214228346945880,&
  & 45.0,6.2166480,17.919738,3.8542520,0.8403260,15.173498,0.0707890,0.85612100,33.8894840,121.686688,9.02951700,0.99545200,45.0,1.225662002307150,&
  & 46.0,6.1215110,4.7840630,16.631683,4.3182580,13.246773,0.0625490,0.78403100,8.75139100,34.4899830,0.78403100,0.88309900,46.0,1.236213451657990,&
  & 47.0,6.0738740,17.155437,4.1733440,0.8522380,17.988685,0.0553330,7.89651200,28.4437390,110.376108,0.71680900,0.75660300,47.0,1.256597538627620,&
  & 48.0,6.0809860,18.019468,4.0181970,1.3035100,17.974669,0.0489900,7.27364600,29.1192830,95.8312080,0.66123100,0.60350400,48.0,1.278260006400250,&
  & 49.0,6.1964770,18.816183,4.0504790,1.6389290,17.962912,0.0420720,6.69566500,31.0097910,103.284350,0.61071400,0.33309700,49.0,1.290470357489140,&
  & 50.0,19.325171,6.2815710,4.4988660,1.8569340,17.917318,6.1181040,0.03691500,32.5290470,95.0371820,0.56565100,0.11902400,50.0,1.330549377714000/
!
! 51 - 60
  data ((bb(i,j),i=1,14),j=51,60) / &
  & 51.0,5.3949560,6.5495700,19.650681,1.8278200,17.867833,33.326523,0.03097400,5.56492900,87.1309650,0.52399200,-0.2905060,51.0,1.360402686938670,&
  & 52.0,6.6603020,6.9407560,19.847015,1.5571750,17.802427,33.031656,0.02575000,5.06554700,84.1016130,0.48766000,-0.8066680,52.0,1.395171636062690,&
  & 53.0,19.884502,6.7365930,8.1105160,1.1709530,17.548715,4.6285910,0.02775400,31.8490960,84.4063910,0.46355000,-0.4488110,53.0,1.434934988906490,&
  & 54.0,19.978920,11.774945,9.3321820,1.2447490,17.737501,4.1433560,0.01014200,28.7961990,75.2806880,0.41361600,-6.0659020,54.0,1.461928150660960,&
  & 55.0,17.418675,8.3144440,10.323193,1.3838340,19.876252,0.3998280,0.01687200,25.6058280,233.339674,3.82691500,-2.3228020,55.0,1.504734672998740,&
  & 56.0,19.747344,17.368476,10.465718,2.5926020,11.003653,3.4818230,0.37122400,21.2266410,173.834271,0.01071900,-5.1834970,56.0,1.540051791581430,&
  & 57.0,19.966018,27.329654,11.018425,3.0866960,17.335454,3.1974080,0.00344600,19.9554920,141.381979,0.34181700,-21.745489,57.0,1.579460053005520,&
  & 58.0,17.355121,43.988498,20.546650,3.1306700,11.353665,0.3283690,0.00204700,3.08819600,134.907661,18.8329610,-38.386017,58.0,1.606421481188560,&
  & 59.0,21.551311,17.161729,11.903859,2.6791030,9.5641970,2.9956750,0.31249100,17.7167050,152.192827,0.01046800,-3.8710680,59.0,1.720026152321630,&
  & 60.0,17.331244,62.783923,12.160097,2.6634830,22.239951,0.3002690,0.00132000,17.0260010,148.748986,2.91026800,-57.189844,60.0,1.683652354338470/
!
! 61 - 70
  data ((bb(i,j),i=1,14),j=61,70) / &
  & 61.0,17.286388,51.560161,12.478557,2.6755150,22.960947,0.2866200,0.00155000,16.2237550,143.984513,2.79648000,-45.973681,61.0,1.724693404362830,&
  & 62.0,23.700364,23.072215,12.777782,2.6842170,17.204366,2.6895390,0.00349100,15.4954370,139.862475,0.27453600,-17.452166,62.0,1.764012158926800,&
  & 63.0,17.186195,37.156839,13.103387,2.7072460,24.419271,0.2616780,0.00199500,14.7873600,134.816293,2.58188300,-31.586687,63.0,1.797805735087610,&
  & 64.0,24.898118,17.104951,13.222581,3.2661520,48.995214,2.4350280,0.24696100,13.9963250,110.863093,0.00138300,-43.505684,64.0,1.840118637821660,&
  & 65.0,25.910013,32.344139,13.765117,2.7514040,17.064405,2.3739120,0.00203400,13.4819690,125.836511,0.23691600,-26.851970,65.0,1.871313100571290,&
  & 66.0,26.671785,88.687577,14.065445,2.7684970,17.067782,2.2825930,0.00066500,12.9202300,121.937188,0.22553100,-83.279831,66.0,1.903720643077100,&
  & 67.0,27.150190,16.999819,14.059334,3.3869790,46.546471,2.1696600,0.21541400,12.2131480,100.506781,0.00121100,-41.165253,67.0,1.926235156557560,&
  & 68.0,28.174886,82.493269,14.624002,2.8027560,17.018515,2.1209950,0.00064000,11.9152560,114.529936,0.20751900,-77.135221,68.0,1.946913272263270,&
  & 69.0,28.925894,76.173796,14.904704,2.8148120,16.998117,2.0462030,0.00065600,11.4653750,111.411979,0.19937600,-70.839813,69.0,1.963724205103120,&
  & 70.0,29.676760,65.624068,15.160854,2.8302880,16.997850,1.9776300,0.00072000,11.0446220,108.139150,0.19211000,-60.313812,70.0,1.979258983014140/
!
! 71 - 80
  data ((bb(i,j),i=1,14),j=71,80) / &
  & 71.0,30.122865,15.099346,56.314899,3.5409800,16.943730,1.8830900,10.3427640,0.00078000,89.5592480,0.18384900,-51.049417,71.0,1.995446352857520,&
  & 72.0,30.617033,15.145351,54.933548,4.0962530,16.896157,1.7956130,9.93446900,0.00073900,76.1897070,0.17591400,-49.719838,72.0,2.006728379573300,&
  & 73.0,31.066358,15.341823,49.278296,4.5776650,16.828321,1.7087320,9.61845500,0.00076000,66.3462020,0.16800200,-44.119025,73.0,2.014810640008000,&
  & 74.0,31.507901,15.682498,37.960127,4.8855090,16.792113,1.6294850,9.44644800,0.00089800,59.9806750,0.16079800,-32.864576,74.0,2.024233087057100,&
  & 75.0,31.888456,16.117103,42.390296,5.2116690,16.767591,1.5492380,9.23347400,0.00068900,54.5163710,0.15281500,-37.412681,75.0,2.032203945323960,&
  & 76.0,32.210298,16.678440,48.559907,5.4558390,16.735532,1.4735310,9.04969500,0.00051900,50.2102010,0.14577100,-43.677954,76.0,2.037562299152530,&
  & 77.0,32.004437,1.9754540,17.070104,15.939454,5.9900030,1.3537670,81.0141720,0.12809300,7.66119600,26.6594030,4.01889300,77.0,2.070931471354530,&
  & 78.0,31.273891,18.445441,17.063745,5.5559330,1.5752700,1.3169920,8.79715400,0.12474100,40.1779940,1.31699700,4.05039400,78.0,2.070307656542580,&
  & 79.0,16.777389,19.317156,32.979682,5.5954530,10.576854,0.1227370,8.62157000,1.25690200,38.0088210,0.00060100,-6.2790780,79.0,2.058834893082310,&
  & 80.0,16.839889,20.023823,28.428565,5.8815640,4.7147060,0.1159050,8.25692700,1.19525000,39.2472260,1.19525000,4.07647800,80.0,2.063026432193550/
!
! 81 - 90
  data ((bb(i,j),i=1,14),j=81,90) / &
  & 81.0,16.630795,19.386615,32.808570,1.7471910,6.3568620,0.1107040,7.18140100,1.11973000,90.6602620,26.0149780,4.06693900,81.0,2.055981336260000,&
  & 82.0,16.419567,32.738592,6.5302470,2.3427420,19.916475,0.1054990,1.05504900,25.0258900,80.9065960,6.66444900,4.04982400,82.0,2.050009704660440,&
  & 83.0,16.282274,32.725137,6.6783020,2.6947500,20.576559,0.1011800,1.00228700,25.7141450,77.0575500,6.29188200,4.04091400,83.0,2.047846539848800,&
  & 84.0,16.289164,32.807170,21.095164,2.5059010,7.2545890,0.0981210,0.96626500,6.04662200,76.5980710,28.0961280,4.04655600,84.0,2.048640858343320,&
  & 85.0,16.011461,32.615549,8.1138990,2.8840820,21.377867,0.0926390,0.90441600,26.5432560,68.3729610,5.49951200,3.99568400,85.0,2.043577109058100,&
  & 86.0,16.070228,32.641105,21.489659,2.2992180,9.4801840,0.0904370,0.87640900,5.23968700,69.1884770,27.6326400,4.02097700,86.0,2.053343320831720,&
  & 87.0,16.007386,32.663830,21.594351,1.5984970,11.121192,0.0870310,0.84018700,4.95446700,199.805805,26.9051060,4.00347200,87.0,2.056402858312870,&
  & 88.0,32.563691,21.396671,11.298093,2.8346880,15.914965,0.8019800,4.59066600,22.7589730,160.404392,0.08354400,3.98177300,88.0,2.062350469636460,&
  & 89.0,15.914053,32.535042,21.553976,11.433394,3.6124090,0.0805110,0.77066900,4.35220600,21.3816220,130.500748,3.93921200,89.0,2.059606329205960,&
  & 90.0,15.784024,32.454898,21.849222,4.2390770,11.736191,0.0770670,0.73513700,4.09797600,109.464113,20.5121380,3.92253300,90.0,2.076093709574790/
!
! 91 - 92
  data ((bb(i,j),i=1,14),j=91,92) / &
  & 91.0,32.740208,21.973674,12.957398,3.6838320,15.744058,0.7095450,4.05088100,19.2315420,117.255006,0.07404000,3.88606600,91.0,2.084768702621660,&
  & 92.0,15.679275,32.824305,13.660459,3.6872610,22.279435,0.0712060,0.68117700,18.2361570,112.500040,3.93032500,3.85444400,92.0,2.096286630784310/
!
  if (z<0 .or. z>92) then
    write(unit=snum,fmt=*) z
    write(*,*) "Atomic number Z="//trim(adjustl(snum))// &
             & " is not in the table of Waasmaier and Kirfel."
    stop
  end if
!
  do i = 1, 14
	a(i) = bb(i, z)
  end do
!
  return
end subroutine getwaki


! ******************************************************************** !	
!
!     fx parameter table of Waasmaier and Kirfel (X-ray form factors)
!     copied from muSTEM
!     original table from: Acta Cryst. A 51 (1995) 416 - 431
!
!     fx(s) = Sum[ a_i * Exp[ -b_i * s^2 ], (i,1,5) ] + c
!             for s = 0 ... 6 1/A
!     fx(s) = d_1 / (d_2 + s^2)
!             for s > 6 1/A
!
!     Coefficients of the extended table as a list a
!     with in total 14 entries.
!     Entry key:
!     a(1) = Z (atomic number)
!     a(2) = a1 (parameterization coefficient)
!     ...
!     a(6) = a5 (parameterization coefficient)
!     a(7) = b1 (parameterization coefficient)
!     ...
!     a(11) = b5 (parameterization coefficient)
!     a(12) = c  (parameterization coefficient)
!     a(13) = d_1 ( = Z, parameterization coefficient)
!     a(14) = d_2 (parameterization coefficient)
!
! ---------------------------------------------------------------------
!
! real*4 function waki(a,z) returns electron form factors
!     doesn't add the ionic charge term for ions
!
! ---------------------------------------------------------------------
function waki(a,s)
  ! no variable implication allowed
  implicit none
  ! parameters              
  real*4, parameter :: c1 = 0.0239336609804 ! m0 e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) *10^-10 [ -> A^-1 ]
  real*4, parameter :: al = 0.02 ! inverse yukawa range in A (only used for small s2 in ions to calculate 1/(s2+al**2)) 
  real*4, parameter :: s20 = 0.001 ! small s^2 regime limit [1/A^2]
  ! interface variables
  real*4, intent(in) :: a(14)
  real*4, intent(in) :: s
  ! function declaration
  real*4 :: waki
  ! local variables
  integer*4 :: n, m
  real*4 :: s2, total, xbs
  real*4 :: deltak, deltak1 !, deltak2 ! charge difference to atomic number
  real*4 :: fxray, s2l !, a2
  ! 
  s2 = s*s ! s^2
  ! distinguish between small s^2 and large s^2
  if ( s2 < s20 ) then ! small s range
    total = a(12) ! c
    do n=2, 6
      total = total + a(n) ! sum over all a_i
    end do
    ! total = c + Sum[a_i]_i=1,5
    !
    ! if 0 == deltak = Z - sum a_i - c then we have a neutral atom
    deltak = a(13) - total
    deltak1 = real(nint(deltak)) ! integer part = ionic charge
    !         Z - sum_ai - INT( Z - sum_ai )
    ! deltak2 = deltak - deltak1   ! fractional part = table inconsistency
    ! calculate fe
    total = 0.0
    do n=2, 6
	  m = n + 5 ! index for b_i
	  ! approximative form of a_i * Exp[ -b_i * s^2 ] / s^2 for small s^2
	  !                       a_i + a_i * b_i * ( 1 - 0.5*b_i*s^2 )
	  total = total + a(n)*a(m)*(1.0 - 0.5*a(m)*s2)
    end do
    ! add ionic charge contribution (integer charge deltak1)
    !total = total + deltak1/(s2 + al**2)  ! ionic charge contribution
    total = total + deltak1/(s2 + al**2)  ! ionic charge contribution
    ! add a correction resolving the discontinuity at s2 = s20
    ! total = total + deltak2/s20
    ! apply prefactor  m0 e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) *10^-10 [ -> A^-1 ]
    waki = c1 * total
    !
  else ! large s range
    !
    if ( s2 > 36.0 ) then ! out of parameterization range ?
	  ! Large scattering vector limit
	  s2l = s2 + a(14)
      waki = c1 * a(13) / s2l
    else ! scattering vector in parameterisation range
      fxray = a(12) ! c
      do n=2, 6 ! indices to a_i
	    m = n + 5 ! indices to b_i
	    xbs = -a(m) * s2 ! argument of exponential 
	    xbs = exp(xbs)   ! exponential
	    fxray = fxray + a(n) * xbs ! accumulation of a_i * exp(-b_i s^2)
      end do
      ! apply Mott formula to fx to get fe
      waki = c1 * (a(13) - fxray) / s2
    end if
  end if
  !
  return
end function waki