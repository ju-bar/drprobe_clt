!**********************************************************************!
!**********************************************************************!
!
! F90 Module - CellSlice (CS)
!
! Author: Dr. J. Barthel
!         Forschungszentrum Juelich GmbH
!         Jülich, Germany
!         ju.barthel@fz-juelich.de
!         first version: 13.12.2008
!         last version: 21.03.2018
!
! Purpose: Data and Methods to create Potential slices
!          for TEM, from given atomic data in form of supercells
!          usable by a multislice algorithm
!
!          - supercell atomic structure loading & interpretation
!            formats: CEL, CIF
!          - slice potential calculation according to EMS by
!            P. Stadelman (EPFL,CH)
!          - frozen lattice generation included
!          - shared code for projects ProbeInspectorDlls and simepw
!            -> Dr. Probe, CELSLC
!
! Linkage: random.f90 (external subs for random number generation)
!          wekoscatt.f90 (function to calculate scattering amplitudes
!                    by weikenmeier and kohl)
!          wakiscatt.f90 (atomic form factors according to Waasmaier
!                    and Kirfel)
!          integration.f90 (numerical integration routines)
!          fscatab.f90 (module handling external scattering table data)
!          fitfeprm.f90 (module handling external scattering data parameterization)
!          SFFTs.f  (sub for FFT)
!          cifio.f90 (module handling input and output with CIF structure data files)
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
!**********************************************************************!
MODULE CellSlicer

  implicit none
  
  SAVE


!**********************************************************************!
!
! SUBS AND FUNCTIONS
!
!  public ::
  public :: CS_INIT, CS_UNINIT
  public :: CS_LOAD_EMSCELL
  public :: CS_LOAD_EMSCELLINFO
  public :: CS_LOAD_EMSCELLINFO2
  public :: CS_SAVE_EMSCELL
  public :: CS_LOAD_CIFCELL
  public :: CS_LOAD_CIFCELLINFO
  public :: CS_LOAD_CIFCELLINFO2
  public :: CS_SAVE_CIFCELL
  public :: CS_DEALLOC_CELLMEM
  public :: CS_ALLOC_CELLMEM
  public :: CS_DEALLOC_SLICEMEM
  public :: CS_ALLOC_SLICEMEM
  public :: CS_DEALLOC_SCAMPMEM
  public :: CS_SETSLICE_AUTO
  public :: CS_SETSLICE_EQUIDIST
  public :: CS_GETSLICE_PRO
  public :: CS_GETSLICE_PRO2
  public :: CS_GETSLICE_POT
  public :: CS_GETSLICE_POT2
  public :: CS_GETSLICE_PGR
  public :: CS_GETCELL_POT
  public :: CS_GETCELL_VOL
  public :: CS_PREPARE_SCATTAMPS
  public :: CS_WRITE_FFDEC
  public :: CS_WRITE_F2DEC
  public :: CS_GET_MEANINNERPOT
  public :: CS_POLINT
  public :: CS_GETSCATTAMP
  public :: CS_GETSLICETITLE
  public :: CS_GET_SAMESITE
  public :: CS_GET_NNAT
  public :: CS_ORIENT_CELL
  public :: CS_SHIFT_ATOMS
  public :: CS_DICE_PARTIALOCC
  public :: CS_SUGGEST_NSLCEQUI
!
!  private ::
  private :: CS_ERROR
  private :: CS_MESSAGE
  private :: CS_GETFREELUN
  private :: CS_WL2HT, CS_HT2WL
  private :: CS_CCFFT2D
  private :: CS_PROG_START
  private :: CS_PROG_STOP
  private :: CS_PROG_UPDATE
  private :: CS_PROG_ABORT

!
!
!**********************************************************************!

!**********************************************************************!
!  
! PARAMETERS
!
! logical unit number limits
  integer*4, public, parameter :: CS_LUNMIN = 20
  integer*4, public, parameter :: CS_LUNMAX = 99
!
! number of cell file record end characters
  integer*4, public, parameter :: CS_CFRECNUM = 1
!
! math and physics constants used by CellSlicer routines:
! pi
  real*4, public, parameter :: CS_pi    = 3.1415926536 ! Pi
  real*4, public, parameter :: CS_tpi   = 6.2831853072 ! 2*Pi
  real*4, public, parameter :: CS_fpi   = 12.566370614 ! 4*Pi
  real*4, public, parameter :: CS_pi4   = 0.7853981634 ! Pi/4
  real*4, public, parameter :: CS_rr8p2 = 0.11253953952 ! 1/(Sqrt(8*Pi^2))
! scale degree to radian
  real*4, public, parameter :: CS_rd2r  = 0.01745329252 ! Pi/180
! electron rest mass energy in keV
  real*4, public, parameter :: CS_elm0  = 510.998946269 ! m0*c^2 [keV]
! electron wave length scale in nm/kV
  real*4, public, parameter :: CS_elwl  = 1.23984197396 ! c*h/e [nm/kV]
! inverse Yukawa range to handle divergence in ionic charge potentials. (This is a fudge.)
  real*4, public, parameter :: CS_iyr   = 0.2 ! [1/nm]
! pre-factor for form factor calculations  
  real*4, public, parameter :: CS_scaprea = 0.0239336609804 ! m0 e^2 / ( 2 h^2 ) / ( 4 Pi eps0 ) *10^-10 [ -> A^-1 ]
! prefector for potential calculations: V0
  real*4, public, parameter :: CS_v0    = 0.03809982080 ! hbar^2 / ( 2 * m0 * e ) * (10^9)^2 ) [ V nm^2 ]
! interaction constant: sigma
!                     = m0 * e / (2*Pi * hbar^2) * (10^-9)^2  [ V^-1 nm^-2 ]
!                     = (1/(4*Pi)) * (2 * m0 * e) / hbar^2 * (10^-9)^2  [ V^-1 nm^-2 ]
!                     The (1/(4*Pi)) compensates again the 4*Pi multiplied to the atomic form factors
!                     in the code fscatt.f taken from EMS. As the original factors by Weickenmeier and
!                     Kohl do not contain this factor, sigma should be referenced as
!                     sigma = 2 * m0 * e / hbar^2 [ V / m**2] = 2 * m0 / hbar^2 [ J / m**2]
  real*4, public, parameter :: CS_sig   = 2.0886573708 ! m0 * e / (2*Pi * hbar^2) * (10^-9)^2  [ V^-1 nm^-2 ]
!
! Note on the essential physical constants used in this code:
! * Atomic form factors (scattering factors):
!        CS_scaprea is the pre-factor for calculating
!        atomic form factors in units of [A] from the Mott formula.
! * Scattering potentials:
!        Calculating Uf = fsca * CS_v0 / cellvolume / rc
!        gives potentials in Volt units without relativistic correction.
!        cellvolume = slice or unit cell volume in nm^3
!        rc = relativistic correction = (E0 + U)/E0
!        E0 = CS_elm0 electron rest mass in keV
!        U  = ht electron kinetic energy in keV
! * Phase shift (phase gratings):
!        Calculating
!          pha  = fsca * CS_v0 * CS_sig / cellvolume * dZ * Lambda 
!        gives the phase shift to the electron.
!        dZ = thickness of slice or cell in nm.
!        Lambda = electron wave length (vacuum) in nm.
! * Side note:
!        (CS_v0*CS_sig = 1/(4*Pi)
!        This cancels the 4*Pi from the scattering factors as produced by
!          fsca = FST_GETSCAR(sym, g, dwf, ht, dwfflg )
!          fsca = wakiscar(g, dwf, z, ht, dwfflg)
!          fsca = wekoscar(g, dwf, z, ht, dwfflg)
!        which are all in units of [nm] and contain relativistic correction.
!
! atomic distance threshold in 3D in nm.
  real*4, parameter, public :: CS_atdst_thrnm = 0.05 ! 50 pm is below the shortest bond length (H-H 74 pm)
!
!
! std string or table length
  integer*4, public, parameter :: CS_ll = 1024
!
! atom identification string length
  integer*4, public, parameter :: CS_atsl = 4
!
!! FFT min/max size
!  integer*4, public :: CS_FFT_BOUND_MAX = 8192
!  integer*4, public :: CS_FFT_BOUND_MIN = 128
!  integer*4, public :: CS_FFT_BOUND
!  DATA CS_FFT_BOUND /2048/
  !integer*4, private, parameter :: FFT_BOUND = 4096
  !integer*4, private, parameter :: FFT_BOUND = 8192
! INTERNAL FFT CODE min and max size / see routines made available by FFTs.f
  integer*4, private, parameter :: NFT_MIN = 128
  integer*4, private, parameter :: NFT_MAX = 8192
!
! Propagator fix relative size aperture (to cut of scattering wrap arounds)
  real*4, private, parameter :: CS_PROSIZE_THR = 2.0 / 3.0
!
! Max. atomic number
  integer*4, private, parameter :: CS_maxz = 98
!
! scattering amplitude rel. threshold
  real*4, private, parameter :: CS_scamprthr = 1.0e-5
!
! switch to bandwidth limiting of potentials and propagators
! potential bandwidth limitation
! CS_do_bwl_pot = .true. : apply bandwidth limiting aperture (default)
! CS_do_bwl_pot = .false. : do not apply bandwidth limiting aperture (for test purpose only)
! phase grating bandwisth limitation
! CS_do_bwl_pgr = .true. : apply bandwidth limiting aperture (default)
! CS_do_bwl_pgr = .false. : do not apply bandwidth limiting aperture (for test purpose only)
! Remark: Bandwidth limiting takes place anyway by the choice of
!         the potential sampling. Without the extra bandwidth limiting
!         aperture, the limitation is by the rectangular array form and may
!         be asymmetric if asymmetric reciprocal space sampling is applied.
!         With sampling rates chosen fine enough, the b.w.l. should not play
!         a significant role.
!         The extra bandwidth limitation applies a circular aperture in
!         the (physical) reciprocal space. It may appear elliptical
!         in case of asymmetric reciprocal space samplings, however, it guarantees
!         that the rotational symmetry of the round atomic potentials is preserved.
  logical, private, parameter :: CS_do_bwl_pot = .false. ! .false. is default
!         The bandwidth limitation of phase gratings seems unphysical
!         at first glance and is in fact a numerical measure applied
!         solely to avoid alias beams to appear after scattering.
  logical, private, parameter :: CS_do_bwl_pgr = .true.
!
! Atomic form factor oversampling
! We precalculate the spherical symmetrical atomic form factors on a
! linear grid to save time, especially since the calculation of absorptive
! form factors requires more computation time.
! In order to provide additional smoothness for transferring this data
! to 2D and 3D potential grids, a significant oversampling needs to be
! applied.
! Keep the oversampling rate above 10!
  integer*4, private, parameter :: CS_scaf_oversamp = 12
!
! Atomic form factor interpolation mode
! An interpolation is applied when transferring the atomic form factors
! from the precalculated linear grid to a 2D slice or 3D volume grid.
! This parameter triggers, which type of interpolation is used
!   CS_scaf_iptype = 0 -> nearest neighbor interpolation
!   CS_scaf_iptype = 1 -> linear interpolation (default)
!   CS_scaf_iptype > 1 -> polynomial interpolation order
  integer*4, private, parameter :: CS_scaf_iptype = 4
!
! Atomic form factor table to be used, can be triggered as an option
  integer*4, public :: CS_scaf_table
  DATA CS_scaf_table /1/
!   CS_scaf_table = 1 ! -> Weickenmeier & Kohl (default & fallback)
!   CS_scaf_table = 2 ! -> Waasmaier & Kirfel
!
!**********************************************************************!

!**********************************************************************!
!
! VARIABLES
!
! error counter  
  integer*4, public :: CS_err_num, CS_warn_num
  DATA CS_err_num /0/
  DATA CS_warn_num /0/
!  
! status
  integer*4, public :: CS_status
  DATA CS_status /0/
!
! last error message
  character(len=CS_ll), public :: CS_errmsg, CS_warnmsg
!
! console message flag
  integer*4, public :: CS_doconsolemsg
  DATA CS_doconsolemsg /0/
!
! cell file record end characters
  character*1, dimension(1:CS_CFRECNUM), public :: CS_CFREC
  DATA CS_CFREC(1) /Z'A'/
!
! workspace
!  complex*8, allocatable, dimension(:,:), private :: cw
!
! progress output control
  integer*4, private :: CS_prog_on
  DATA CS_prog_on /0/
  integer*4, private :: CS_prog_num
  DATA CS_prog_num /100/
  integer*4, private :: CS_prog_idx
  DATA CS_prog_idx /0/
  real*4, private :: CS_prog_pproc
  DATA CS_prog_pproc /1.0/
  real*4, private :: CS_prog_last
  DATA CS_prog_last /-1000.0/  
!
! supercell data
!
!   number of atoms per super cell
  integer*4, public :: CS_numat
  DATA CS_numat /0/
!
!   super cell size
  real*4, public :: CS_scsx, CS_scsy, CS_scsz
  DATA CS_scsx, CS_scsy, CS_scsz /0.0,0.0,0.0/
!
!   supercell axis angles [deg]
  real*4, public :: CS_scalpha, CS_scbeta, CS_scgamma
  DATA CS_scalpha, CS_scbeta, CS_scgamma /90.0,90.0,90.0/
!
!   orthogonality flag
  integer*4, public :: CS_scorth_check
  DATA CS_scorth_check /0/
!
!
!   atomic data
!
!     memory allocation flag
  logical, public :: CS_cellmem_allocated
  DATA CS_cellmem_allocated /.FALSE./
!
!     atomic type
  character(len=CS_atsl), allocatable, dimension(:), public :: CS_attype
!
!     atomic numbers
  integer*4, allocatable, dimension(:), public :: CS_atnum
!
!     atomic charges
  real*4, allocatable, dimension(:), public :: CS_atcrg

!
!     atomic position in supercell
  real*4, allocatable, dimension(:,:), public :: CS_atpos
!
!     position occupancy
  real*4, allocatable, dimension(:), public :: CS_atocc
!
!     debye-waller factor
  real*4, allocatable, dimension(:), public :: CS_atdwf
!
!     Link list of atoms sharing the same site.
!     This is a list of indices, one index per atom.
!     Links are used to displace linked atoms in the same way
!     for frozen lattice calculations. The list is not used
!     at all for calculations without frozen lattices.
!     Initially and by default all indices are zero,
!     which means, that this atom of the cell is independent
!     from other atoms of the cell.
!     CS_atlnk(i)>0 means that atom #(i) is linked to the
!     atom #(CS_atlnk(i)), with CS_atlnk(i) < i.
!     The links are reset before slices are generated, or before
!     a 3D potential is calculated.
  integer*4, allocatable, dimension(:), public :: CS_atlnk
!
!
! scattering amplitudes
!
!   usage of external scattering table data
  integer*4, public :: CS_useextsca
  DATA CS_useextsca /0/ ! 0 = off, default, use Weickenmeier&Kohl data
                        ! 1 = try to use what is present in module FSCATAB
                        !     if not possible, switch back to W&K
  
!
!   preparation + memory allocation flag
  logical, public :: CS_scattamp_prepared
  DATA CS_scattamp_prepared /.FALSE./
!
!   number of identified atom types = number of different scattering amplitudes
  integer*4, public :: CS_scampnum
  DATA CS_scampnum /0/
!
!   pointer identifying the scattering table to be used for each atom in the cell
  integer*4, public, allocatable, dimension(:) :: CS_scampptr
!
!   atom type-wise atomic symbol for scattering data
  character(len=CS_atsl), public, allocatable, dimension(:) :: &
                                                & CS_scampsym
!
!   atom type-wise atomic numbers for scattering data
  integer*4, public, allocatable, dimension(:) :: CS_scampatn
!
!   atom type-wise ionic charge
  real*4, public, allocatable, dimension(:) :: CS_scampcrg
!
!   atom type-wise Debye-Waller factors
  real*4, public, allocatable, dimension(:) :: CS_scampdwf
  
!   atom type-wise temporary data used during calculations
  real*4, public, allocatable, dimension(:,:) :: CS_scamptmp

!!
!!   atom type-wise max. scattering vector in 1/nm
!  real*4, public, allocatable, dimension(:) :: CS_scampgmax
!
!   dimension of allocated scattering data
  integer*4, public :: CS_scadim
  DATA CS_scadim /0/
!
!   scattering data sampling rate
  real*4, public :: CS_scaitog
  DATA CS_scaitog /0.0/
!
!
!   absorption calculation parameters
  integer*4, public :: CS_useabsorption
  DATA CS_useabsorption /0/
  real*4, public :: CS_absorptionprm
  DATA CS_absorptionprm /0.0/
!
!
  integer*4, public :: CS_useppot ! flag: using projected potentials
  DATA CS_useppot /1/ ! on by default
!   complex scattering amplitudes (will be allocated when needed)
  complex*8, public, allocatable, dimension (:,:) :: CS_scampdat
  complex*8, public, allocatable, dimension (:,:,:) :: CS_scaff2d
  complex*8, public, allocatable, dimension (:,:,:) :: CS_scagn
  real*4, public, allocatable, dimension (:,:,:) :: CS_scadwf
!
!
! slice data
!
!   memory allocation flag
  logical, public :: CS_slicemem_allocated
  DATA CS_slicemem_allocated /.FALSE./
!
!   repetition of supercell data in x- and y-direction
!     (for construction of slices containing the supercell several times,
!      for multislice on a larger area)
  integer*4 :: CS_repx, CS_repy, CS_repz
!
!   numeric size of slices
  integer*4 :: CS_sdimx, CS_sdimy, CS_sdimz
!
!   reals-space sampling for slices, results from CS_scsx*CS_repx/CS_sdimx
  real*4, public :: CS_sampx, CS_sampy, CS_sampz
!
!   number of slices per supercell
  integer*4, public :: CS_nspsc
  DATA CS_nspsc /0/
!
!   slicing, z-limits of slices, thickness
!   - allows non-equidistant slicing (TODO)
  real*4, allocatable, dimension(:,:), public :: CS_slczlim
!
  real*4, public :: CS_slcztol
  DATA CS_slcztol /0.003/ ! 3 pm default tolerated minimum slice thickness
!
!   number of atoms in slice
  integer*4, allocatable, dimension(:), public :: CS_slcatnum
!
!   slicing, atom association to slices
!    (contains indices of CS_at*** arrays,
!     is thus a quasi pointer to these lists)
  integer*4, allocatable, dimension(:,:), public :: CS_slcatacc
!
!   accumulated maximum value of phase shifts (do reset this in the beginning
  real*4, public :: CS_maxphase
  DATA CS_maxphase /0.0/
!
!   potential backup functionality setup
!   - setup these variable from external code
! 
  ! flag the backup OFF (0) or ON (/=0)
  integer*4, public :: CS_backup_pot_flg            
  DATA CS_backup_pot_flg /0/
  ! potential backup data (last calculated potential)
  complex*8, allocatable, dimension(:,:), public :: CS_backup_pot 
!
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
  CONTAINS
!**********************************************************************!
!**********************************************************************!  



  
  
!**********************************************************************!
!
! CS_INIT
!
! subroutine to initialise data and allocate memory according to
! previously set variables
!
subroutine CS_INIT()

  implicit none
  
  integer*4 :: nalloc
  
  CS_err_num = 0
  CS_errmsg = ""
  CS_status = 0
  
  call InitRand()
  
  if (allocated(CS_backup_pot)) then ! deallocate previous potential backup memory
    deallocate(CS_backup_pot,stat=nalloc)
  end if
  
  CS_status = 1
  
  return

end subroutine CS_INIT
!**********************************************************************!

!**********************************************************************!
!
! CS_UNINIT
!
! subroutine to uninitialise data and deallocate memory
!
subroutine CS_UNINIT()

  implicit none
  
  integer*4 :: nerr, nalloc
  
  !
  ! deallocate atom data array
  !
  nerr = 0
  if (CS_cellmem_allocated) then
    call CS_DEALLOC_CELLMEM(nerr)
    if (nerr/=0) return
  end if
  if (CS_slicemem_allocated) then
    call CS_DEALLOC_SLICEMEM(nerr)
    if (nerr/=0) return
  end if
  if (CS_scattamp_prepared) then
    call CS_DEALLOC_SCAMPMEM(nerr)
    if (nerr/=0) return
  end if
!  if (allocated(cw)) then
!    deallocate(cw,stat=nalloc)
!  end if
  if (allocated(CS_backup_pot)) then ! deallocate previous potential backup memory
    deallocate(CS_backup_pot,stat=nalloc)
  end if
  
  CS_status = 0
  
  return

end subroutine CS_UNINIT
!**********************************************************************!

!**********************************************************************!
!
! CS_ERROR
!
! subroutine handling error event (simple)
!
subroutine CS_ERROR(sErrMsg)

  implicit none
  
  character(len=*), intent(in) :: sErrMsg
  
  CS_err_num = CS_err_num + 1
  CS_errmsg = "Error: "//trim(sErrMsg)
  call CS_MESSAGE(trim(CS_errmsg))
  
  return

end subroutine CS_ERROR
!**********************************************************************!

!**********************************************************************!
!
! CS_MESSAGE
!
! subroutine handling console message event (simple, flagged)
!
subroutine CS_MESSAGE(smsg)

  implicit none
  
  character(len=*), intent(in) :: smsg
  
  if (CS_doconsolemsg<=0) return
  write(unit=6,fmt='(A)') trim(smsg)
  
  return

end subroutine CS_MESSAGE
!**********************************************************************!


!**********************************************************************!
!
! CS_PROG_START
!
! subroutine starting a progress output
! make sure not to interrupt this by other output
!
subroutine CS_PROG_START(nmax,pstep)

  implicit none
  
  integer*4, intent(in) :: nmax
  real*4, intent(in) :: pstep
  
  if (CS_doconsolemsg<=0) return
  
  if (CS_PROG_on==1) call CS_PROG_ABORT()
  
  CS_PROG_on = 1
  CS_PROG_num = nmax
  CS_PROG_idx = 0
  CS_PROG_pproc = pstep
  CS_PROG_last = 0.0
  
  OPEN (UNIT=6,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')
  write(unit=6,fmt='("+",A,F5.1,A)') "Progress: ",0.0,"%   "
  
  return

end subroutine CS_PROG_START
!**********************************************************************!


!**********************************************************************!
!
! CS_PROG_UPDATE
!
! subroutine udating a running progress output
! make sure not to interrupt this by other output
!
subroutine CS_PROG_UPDATE(idx)

  implicit none
  
  integer*4, intent(in) :: idx
  real*4 :: pcur
  
  if ( CS_doconsolemsg<=0 .or. CS_PROG_on==0 ) return
  
  CS_PROG_idx = idx
  pcur = 100.*real(idx)/real(CS_PROG_num)
  
  if ( pcur - CS_PROG_last >= CS_PROG_pproc ) then
    CS_PROG_last = CS_PROG_last + CS_PROG_pproc
    write(unit=6,fmt='("+",A,F5.1,A)') "Progress: ",pcur,"%   "
  end if
  
  return

end subroutine CS_PROG_UPDATE
!**********************************************************************!


!**********************************************************************!
!
! CS_PROG_STOP
!
! subroutine stopping a progress output
! make sure not to interrupt this by other output
!
subroutine CS_PROG_STOP(idx)

  implicit none
  
  integer*4, intent(in) :: idx
  real*4 :: pcur
  
  if ( CS_doconsolemsg<=0 .or. CS_PROG_on==0) return
  pcur = 100.*real(idx)/real(CS_PROG_num)
  write(unit=6,fmt='("+",A,F5.1,A)') "Progress: ",pcur,"%   "
  OPEN (UNIT=6,FORM='FORMATTED',CARRIAGECONTROL='LIST')
  
  CS_PROG_on = 0
  CS_PROG_num = 100
  CS_PROG_idx = 0
  CS_PROG_pproc = 1.0
  CS_PROG_last = -1000.0
  
  return

end subroutine CS_PROG_STOP
!**********************************************************************!


!**********************************************************************!
!
! CS_PROG_ABORT
!
! subroutine aborts a progress output
! make sure not to interrupt this by other output
!
subroutine CS_PROG_ABORT()

  implicit none
  
  OPEN (UNIT=6,FORM='FORMATTED',CARRIAGECONTROL='LIST')
  
  CS_PROG_on = 0
  CS_PROG_num = 100
  CS_PROG_idx = 0
  CS_PROG_pproc = 1.0
  CS_PROG_last = -1000.0
  
  return

end subroutine CS_PROG_ABORT
!**********************************************************************!


!**********************************************************************!
!
! CS_GETFREELUN
!
! function returning next free logical unit number
!
! IN/OUTPUT:
!   integer*4 :: nerr         = error flag
!
integer*4 function CS_GETFREELUN(nerr)

  implicit none
  
  integer*4, intent(inout) :: nerr
  
  logical :: isopen             ! isopen request from inquire
  integer*4 :: lun, nret        ! temporary luns
  
  nret = 0
  
  do lun = CS_LUNMIN, CS_LUNMAX
    inquire(unit=lun,opened=isopen,iostat=nerr)
    if (nerr/=0) then
      CS_GETFREELUN = 0
      return
    end if
    if (.not.isopen) then
      nret = lun
      exit
    end if
  end do
  
  if (nret>0) then
    nerr = 0
    CS_GETFREELUN = nret
  else
    nerr = -1
  end if
  
  return

end function CS_GETFREELUN
!**********************************************************************!

!**********************************************************************!
!
! CS_DEALLOC_CELLMEM
!
! subroutine for deallocation of cell data memory
!
! IN/OUTPUT:
!   integer*4 :: nerr = error flag
!
subroutine CS_DEALLOC_CELLMEM(nerr)
  implicit none
  integer*4, intent(inout) :: nerr
  nerr = 0
  if (CS_cellmem_allocated) then
    deallocate(CS_attype, CS_atnum, CS_atlnk, CS_atcrg, &
     &         CS_atpos, CS_atocc, CS_atdwf, stat=nerr)
    if (nerr==0) then
      CS_cellmem_allocated = .FALSE.
    end if
  end if
  return
end subroutine CS_DEALLOC_CELLMEM
!**********************************************************************!

!**********************************************************************!
!
! CS_ALLOC_CELLMEM
!
! subroutine for allocation of cell data memory
!
! cell data memory are array containing data for each atom defined with
! a loaded cell file.
! - CS_attype = atom symbol string (read from file)
! - CS_atnum  = atomic number assigned to symbol (set automatically)
! - CS_atcrg  = ionic charge (set automatically from symbol)
! - CS_atpos  = atomic positions relative in super cell
! - CS_atocc  = atom position occupancy (usually = 1)
! - CS_atdwf  = Debye-Waller factor (nm^2)
! - CS_atlnk  = atom link list (shared sites)
!
! INPUT:
!   integer*4 :: n    = number of atoms
!
! IN/OUTPUT:
!   integer*4 :: nerr = error flag
!
subroutine CS_ALLOC_CELLMEM(n,nerr)
  implicit none
  integer*4, intent(in) :: n
  integer*4, intent(inout) :: nerr
  integer*4 :: i
  nerr = 0
  i = n+3
  if (.not.CS_cellmem_allocated) then
    allocate(CS_attype(i), stat=nerr)
    if (nerr==0) allocate(CS_atnum(i), stat=nerr)
    if (nerr==0) allocate(CS_atcrg(i), stat=nerr)
    if (nerr==0) allocate(CS_atpos(3,i), stat=nerr)
    if (nerr==0) allocate(CS_atocc(i), stat=nerr)
    if (nerr==0) allocate(CS_atdwf(i), stat=nerr)
    if (nerr==0) allocate(CS_atlnk(i), stat=nerr)
    if (nerr==0) then
      CS_attype = ""
      CS_atnum = 0
      CS_atlnk = 0
      CS_atcrg = 0.0
      CS_atpos = 0.0
      CS_atocc = 0.0
      CS_atdwf = 0.0
      CS_numat = n
      CS_cellmem_allocated = .TRUE.
    end if
  else
    nerr = -1
  end if
  return
end subroutine CS_ALLOC_CELLMEM
!**********************************************************************!

!**********************************************************************!
!
! CS_DEALLOC_SLICEMEM
!
! subroutine for deallocation of slice data memory
!
! IN/OUTPUT:
!   integer*4 :: nerr = error flag
!
subroutine CS_DEALLOC_SLICEMEM(nerr)
  implicit none
  integer*4, intent(inout) :: nerr
  nerr = 0
  if (CS_slicemem_allocated) then
    deallocate(CS_slczlim, CS_slcatacc, CS_slcatnum, stat=nerr)
    if (nerr==0) then
      CS_slicemem_allocated = .FALSE.
    end if
  end if
  return
end subroutine CS_DEALLOC_SLICEMEM
!**********************************************************************!

!**********************************************************************!
!
! CS_ALLOC_SLICEMEM
!
! subroutine for allocation of slice data memory
!
! requires:
! - cell data memory to be allocated, if not allocated, returns error
!
! slice data memory are array containing data for each slice.
! - CS_slczlim   = z coordinate limits of slice, rel. to super cell
! - CS_slcatnum  = number of atoms in slice (list len in CS_slcatacc)
! - CS_slcatacc  = atom index-pointer to supercell data arrays
!
! INPUT:
!   integer*4 :: n    = number of slices
!
! IN/OUTPUT:
!   integer*4 :: nerr = error flag
!
subroutine CS_ALLOC_SLICEMEM(n,nerr)
  implicit none
  integer*4, intent(in) :: n
  integer*4, intent(inout) :: nerr
  integer*4 :: i
  nerr = 0
  i = n+3
  ! need size of cell data, thus check if cell data is ready
  if (.not.CS_cellmem_allocated) then
    nerr = -1
    return
  end if
  ! check allocation state
  if (CS_slicemem_allocated) then
    call CS_DEALLOC_SLICEMEM(nerr)
    if (nerr/=0) return
  end if
  ! ready to allocate
  allocate(CS_slczlim(2,i), stat=nerr)
  if (nerr==0) allocate(CS_slcatnum(i), stat=nerr)
  if (nerr==0) allocate(CS_slcatacc(CS_numat+3,i), stat=nerr)
  if (nerr==0) then
    CS_slczlim = 0.0
    CS_slcatnum = 0
    CS_slcatacc = 0
    CS_nspsc = n
    CS_slicemem_allocated = .TRUE.
  end if
  return
end subroutine CS_ALLOC_SLICEMEM
!**********************************************************************!

!**********************************************************************!
!
! CS_DEALLOC_SCAMPMEM
!
! subroutine for deallocation of scattering amplitude data memory
!
! IN/OUTPUT:
!   integer*4 :: nerr = error flag
!
subroutine CS_DEALLOC_SCAMPMEM(nerr)
  implicit none
  integer*4, intent(inout) :: nerr
  nerr = 0
  if (CS_scattamp_prepared) then
    ! data is prepared, throw away data and prepare again
    deallocate(CS_scampdat, &
     &         CS_scampsym, &
     &         CS_scampptr, &
     &         CS_scampatn, &
     &         CS_scampcrg, &
     &         CS_scampdwf, &
     &         CS_scamptmp, stat=nerr)
    if (allocated(CS_scaff2d)) deallocate(CS_scaff2d, stat=nerr)
    if (allocated(CS_scagn)) deallocate(CS_scagn, stat=nerr)
    if (allocated(CS_scadwf)) deallocate(CS_scadwf, stat=nerr)
    if (nerr==0) then
      CS_scattamp_prepared = .FALSE.
      CS_scadim = 0
      CS_scaitog = 0.0
      CS_scampnum = 0
    end if
  end if
  return
end subroutine CS_DEALLOC_SCAMPMEM
!**********************************************************************!

!**********************************************************************!
!
! CS_WL2HT
!
! real*4 function, returns HT in kV for a given electron wavelength
!
! INPUT:
!   real*4 :: wl = wavelength [nm]
!
function CS_WL2HT(wl)
  implicit none
  real*4, intent(in) :: wl
  real*4 :: CS_WL2HT
  CS_WL2HT = CS_elm0*( sqrt( 1. + (CS_elwl/wl/CS_elm0)**2 ) - 1. )
  return
end function CS_WL2HT
!**********************************************************************!

!**********************************************************************!
!
! CS_HT2WL
!
! real*4 function, returns electron wavelength in nm
!                          for a given high tension
!
! INPUT:
!   real*4 :: ht = high tension [kV]
!
function CS_HT2WL(ht)
  implicit none
  real*4, intent(in) :: ht
  real*4 :: CS_HT2WL
  CS_HT2WL = CS_elwl / sqrt( ht * (2*CS_elm0+ht) )
  return
end function CS_HT2WL
!**********************************************************************!



!**********************************************************************!
!
! CS_CCFFT2D
!
! subroutine performing out of place 2D Fourier transforms of arbitrary
! size complex arrays
!
! Fourier-space array dimensions are scrambled and transposed with
! respect to real space arrays.
!
! ndir >= 0 -> forward FFT from crs to cfs
! ndir <  0 -> backwards FFT from cfs to crs
!
subroutine CS_CCFFT2D(crs,cfs,nx,ny,ndir)

  implicit none
  
  integer*4, parameter :: nmaxft = 8192 ! max. dimensions supported
  
  integer*4, intent(in) :: nx, ny, ndir
  complex*8, intent(inout) :: crs(nx,ny), cfs(ny,nx) 
  
! external routines for 2D Fourier transform (link FFTs.f)
  external :: ODDCC128S, ODDCC256S,  ODDCC512S,  ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S
  
  integer*4 :: nft, nalloc, n1, n2, i
  complex*8, allocatable, dimension(:,:) :: cft
  character*40 :: direc
  
  nalloc = 0
  direc = 'FORWARD'
  n1 = nx
  n2 = ny
  if (ndir<0) then
    direc = 'BACKWARD'
    n1 = ny
    n2 = nx
  end if
  nft = 2**CEILING(LOG(real(max(nx,ny)))/LOG(2.0)) ! next 2^N above max(nx,ny)
  nft = max(nft, 128)
  if (nft>nmaxft) goto 14
  allocate(cft(nft,nft),stat=nalloc)
  if (nalloc/=0) goto 13
  cft(1:nft,1:nft) = cmplx(0.0,0.0)
  if (ndir<0) then ! cfs is input
    do i=1, n2
      cft(1:n1,i) = cfs(1:n1,i)
    end do
  else ! crs is input
    do i=1, n2
      cft(1:n1,i) = crs(1:n1,i)
    end do
  end if
  ! FT or iFT on cft
  select case (nft)
  case  (128)
    call  ODDCC128S(cft,nx,ny,direc)
  case  (256)
    call  ODDCC256S(cft,nx,ny,direc)
  case  (512)
    call  ODDCC512S(cft,nx,ny,direc)
  case (1024)
    call ODDCC1024S(cft,nx,ny,direc)
  case (2048)
    call ODDCC2048S(cft,nx,ny,direc)
  case (4096)
    call ODDCC4096S(cft,nx,ny,direc)
  case (8192)
    call ODDCC8192S(cft,nx,ny,direc)
  end select
  if (ndir<0) then ! crs is output
    do i=1, n1
      crs(1:n2,i) = cft(1:n2,i)
    end do
  else ! cfs is output
    do i=1, n1
      cfs(1:n2,i) = cft(1:n2,i)
    end do
  end if
  deallocate(cft, stat=nalloc)
  !
  return

13 continue !
  call CS_ERROR("FFT array allocation failed.")
  return
14 continue ! 
  call CS_ERROR("Requested FFT size is out of bounds (8192).")
  return
  
end subroutine CS_CCFFT2D
!
!**********************************************************************!



!**********************************************************************!
!
! CS_WRITE_FFDEC
!
! subroutine, writes k-values (1/nm) where the atomic form factors decay
!             to a given value
!
! INPUT:
!   real*4 :: vdectrg           = relative target decay value
!
! IN/OUTPUT:
!   integer*4 :: nerr           = error flag
!
subroutine CS_WRITE_FFDEC(vdectrg,nerr)

  use FSCATAB
  
  implicit none
  
  real*4, parameter :: threq = 0.000001 ! decay equivalence threshold
  real*4, parameter :: thrqs = 0.000001 ! qstep threshold
  
  real*4, intent(in) :: vdectrg
  integer*4, intent(inout) :: nerr

  logical :: dwfflg                 ! dwf-flag  
  integer*4 :: lfu                  ! file unit number
  integer*4 :: ioerr                ! I/O error code
  integer*4 :: i, i1                ! iterators
  integer*4 :: z                    ! atomic number
  real*4 :: vdeccur                 ! current relative decay value
  real*4 :: ff0                     ! form factor at q = 0 in nm
  real*4 :: ffq                     ! form factor at current q in nm
  real*4 :: qcur                    ! current q in 1/nm
  real*4 :: qstep                   ! current q step size in 1/nm
  real*4 :: qstepdir                ! current q step direction (+1 or -1)
  character(len=CS_atsl) :: sym     ! atom symbol
  
  integer*4, external :: getfreelun ! returns free logical file unit, from binio.f90
  
  real*4, external :: wekoscar
  complex*8, external :: wekosca    !(g,dw,z,accvlt,absflg) from wekoscatt.f90,
                                    !    calculates scatt amp complex
  real*4, external :: wakiscar
  complex*8, external :: wakisca    !(g,dw,z,accvlt,absflg) from wakiscatt.f90,
                                    !    calculates scatt amp complex
  
  !
  ! preset error flag
  !
  nerr = 0
  dwfflg = .FALSE.
  
  !
  ! check cell data presence
  !
  
  !
  ! prepare parameters
  !
  if (CS_useextsca==1) then ! use external table fit
    i1 = FST_nat
  else ! use internal table by ...
    if (CS_scaf_table == 2) then ! ... Waasmaier & Kirfel
      i1 = 92
    else ! ... Weickenmeier & Kohl
      i1 = 98
    end if
  end if
  
  !
  if (i1<1) then
    nerr = 1
    call CS_ERROR("No atomic form factor data.")
    return
  end if
  
  !
  ! open the output file
  lfu = getfreelun()
  open (unit=lfu, file='ffdec.txt', status='REPLACE', action='WRITE', iostat=ioerr)
  if (ioerr/=0) then
    nerr = 2
    call CS_ERROR("Failed to open output file 'ffdec.txt'.")
    return
  end if
  !
  ! write the header line
  !
1000 FORMAT(A4,", ",A6,", ",A12,", ",A12)
1001 FORMAT(I4,", ",A6,", ",F12.4,", ",F12.6)
  write(unit=lfu, fmt=1000) "Z", "Symbol", "kmax [1/nm]", "f(kmax)/f(0)"
  !
  ! loop over all atomic numbers starting with 1 to maximum
  !
  do i=1, i1
    !
    ! get ff0, z and symbol
    sym = "                "
    z = 0
    ff0 = 1.
    if (CS_useextsca==1) then ! use external table fit
      z = nint(FST_crg(2,i))
      sym = FST_sym(i)
      ff0 = FST_GETSCAR(sym, 0., 0., 100., dwfflg )
    else ! use internal table by ...
      if (CS_scaf_table == 2) then ! ... Waasmaier & Kirfel
        z = i
        call wakisymb(z, sym(1:2))
        ff0 = wakiscar(0., 0., z, 100., dwfflg )
      else ! ... Weickenmeier & Kohl
        z = i
        call wekosymb(z, sym(1:2))
        ff0 = wekoscar(0., 0., z, 100., dwfflg )
      end if
    end if
    !
    ! loop until the decay value is met
    !
    vdeccur = 1.0
    qcur = 0.
    qstep = 1.
    qstepdir = 1.
    !
    do while( abs(vdeccur-vdectrg)>threq .and. qstep>thrqs ) ! get closer to the target ...
      !
      ! update q
      qcur = qcur + qstepdir*qstep
      !
      ! get ffq
      if (CS_useextsca==1) then ! use external table fit
        ffq = FST_GETSCAR(sym, qcur, 0., 100., dwfflg )
      else ! use internal table by ...
        if (CS_scaf_table == 2) then ! ... Waasmaier & Kirfel
          ffq = wakiscar(qcur, 0., z, 100., dwfflg )
        else ! ... Weickenmeier & Kohl
          ffq = wekoscar(qcur, 0., z, 100., dwfflg )
        end if
      end if
      !
      ! update vdeccur
      vdeccur = ffq / ff0
      if ( abs(vdeccur-vdectrg) <= threq ) exit ! done, qcur is good enough
      !
      if (qstepdir>0.) then ! went right ...
        ! ... expect lower ff
        if (vdeccur < vdectrg) then ! ff went lower
          ! reverse
          qstepdir = -1.
          ! reduce
          qstep = qstep * 0.5
        end if
        ! the case where vdeccur is not lower than vdectrg requires further going right
      else ! went left ...
        ! ... expect higher ff
        if (vdeccur > vdectrg) then ! ff went higher
          ! reverse
          qstepdir = 1.
          ! reduce
          qstep = qstep * 0.5
        end if
        ! the case where vdeccur is not higher than vdectrg requires further going left
      end if
      !
    end do ! convergence loop
    !
    ! export the result
    write(unit=lfu, fmt=1001) z, sym, qcur, vdeccur
    !
  end do ! atomic table loop
  
  close (unit=lfu, iostat=ioerr)
  
  return
  
end subroutine CS_WRITE_FFDEC
!**********************************************************************!



!**********************************************************************!
!
! CS_WRITE_F2DEC
!
! subroutine, writes k-values (1/nm) where the integrated squared
!             atomic form factors reach a given relative threshold
!             below the total power
!
! INPUT:
!   real*4 :: vdectrg           = relative target threshold
!   real*4 :: ht                = electron kinetic energy in keV
!
! IN/OUTPUT:
!   integer*4 :: nerr           = error flag
!
subroutine CS_WRITE_F2DEC(vdectrg, ht, nerr)

  use FSCATAB
  
  implicit none
  
  real*8, parameter :: threq = 0.000001 ! decay equivalence threshold
  real*8, parameter :: thrqs = 0.0001 ! qstep threshold
  
  real*4, intent(in) :: vdectrg, ht
  integer*4, intent(inout) :: nerr

  integer*4 :: lfu                  ! file unit number
  integer*4 :: ioerr                ! I/O error code
  integer*4 :: i, i1                ! iterators
  integer*4 :: z                    ! atomic number
  real*4 :: k0                      ! electron wave number in [1/A]
  real*8 :: vdcur, vdtrg            ! current relative loss and target values
  real*8 :: f20                     ! full integrated f2 in A^2
  real*8 :: f2q                     ! integrated f2 at current q in [A^2]
  real*4 :: qcur, tcur              ! current q in 1/nm and the corresponding theta in rad
  real*4 :: qstep                   ! current q step size in 1/nm
  real*4 :: qstepdir                ! current q step direction (+1 or -1)
  character(len=CS_atsl) :: sym     ! atom symbol
  
  integer*4, external :: getfreelun ! returns free logical file unit, from binio.f90
  
  real*8, external :: wekof2        ! (tmax,dwa,a,b,k0) from wekoscatt.f90,
  real*8, external :: wakif2        ! (tmax,dwa,a,k0) from wakiscatt.f90,
  
  !
  ! preset error flag
  !
  nerr = 0
  k0 = 0.1*sqrt((2.0*CS_elm0 + ht)*ht)/CS_elwl ! k[ht] = ( e/(h*c)*10^-7 [A/kV] ) * Sqrt[ (2*E0_keV + HT_kV)*HT_kV ]
  vdtrg = dble(vdectrg)
  
  !
  ! prepare parameters
  !
  if (CS_useextsca==1) then ! use external table fit
    i1 = FST_nat
  else ! use internal table by ...
    if (CS_scaf_table == 2) then ! ... Waasmaier & Kirfel
      i1 = 92
    else ! ... Weickenmeier & Kohl
      i1 = 98
    end if
  end if
  
  !
  if (i1<1) then
    nerr = 1
    call CS_ERROR("No atomic form factor data.")
    return
  end if
  
  !
  ! open the output file
  lfu = getfreelun()
  open (unit=lfu, file='f2dec.txt', status='REPLACE', action='WRITE', iostat=ioerr)
  if (ioerr/=0) then
    nerr = 2
    call CS_ERROR("Failed to open output file 'f2dec.txt'.")
    return
  end if
  !
  ! write the header line
  !
1000 FORMAT(A4,", ",A6,", ",A12,", ",A12)
1001 FORMAT(I4,", ",A6,", ",F12.4,", ",F12.6)
  write(unit=lfu, fmt=1000) "Z", "Symbol", "kmax [1/nm]", "1-f2(kmax)/f2(0)"
  !
  ! loop over all atomic numbers starting with 1 to maximum
  !
  do i=1, i1
    !
    ! get ff0, z and symbol
    sym = "                "
    z = 0
    f20 = 1.
    if (CS_useextsca==1) then ! use external table fit
      z = nint(FST_crg(2,i))
      sym = FST_sym(i)
      f20 = FST_GETSCAF2(sym, CS_pi, 0., ht)
    else ! use internal table by ...
      if (CS_scaf_table == 2) then ! ... Waasmaier & Kirfel
        z = i
        call wakisymb(z, sym(1:2))
        f20 = wakif2(CS_pi, 0., z, k0)
      else ! ... Weickenmeier & Kohl
        z = i
        call wekosymb(z, sym(1:2))
        f20 = wekof2(CS_pi, 0., z, k0)
      end if
    end if
    !
    ! loop until the decay value is met
    !
    vdcur = 1.0D+0
    qcur = 0.
    qstep = 1.
    qstepdir = 1.
    !
    do while( abs(vdcur-vdtrg)>threq .and. qstep>thrqs ) ! get closer to the target ...
      !
      ! update q [1/A] because k0 is also in [1/A]
      qcur = qcur + qstepdir*qstep
      ! update theta
      ! q/2 = k0 * sin(theta/2)
      ! q   = k0 * sqrt( 2 - 2*cos(theta) )
      ! (q/k0)^2 = 2 - 2*cos(theta)
      ! cos(theta) = 1 - 0.5*(q/k0)^2
      ! theta = acos( 1-0.5*(q/k0)^2 )
      tcur = acos( 1.0 - 0.5*(qcur/k0)**2 )
      !
      ! get f2q
      if (CS_useextsca==1) then ! use external table fit
        f2q = FST_GETSCAF2(sym, tcur, 0., ht)
      else ! use internal table by ...
        if (CS_scaf_table == 2) then ! ... Waasmaier & Kirfel
          f2q = wakif2(tcur, 0., z, k0)
        else ! ... Weickenmeier & Kohl
          f2q = wekof2(tcur, 0., z, k0)
        end if
      end if
      !
      ! update vdeccur
      vdcur = 1.0D+0 - f2q / f20
      if ( abs(vdcur-vdtrg) <= threq ) exit ! done, qcur is good enough
      !
      if (qstepdir>0.) then ! went right ...
        ! ... expect lower ff
        if (vdcur < vdtrg) then ! ff went lower
          ! reverse
          qstepdir = -1.
          ! reduce
          qstep = qstep * 0.5
        end if
        ! the case where vdeccur is not lower than vdectrg requires further going right
      else ! went left ...
        ! ... expect higher ff
        if (vdcur > vdtrg) then ! ff went higher
          ! reverse
          qstepdir = 1.
          ! reduce
          qstep = qstep * 0.5
        end if
        ! the case where vdeccur is not higher than vdectrg requires further going left
      end if
      !
    end do ! convergence loop
    !
    ! export the result
    write(unit=lfu, fmt=1001) z, sym, qcur*10., vdcur
    !
  end do ! atomic table loop
  
  close (unit=lfu, iostat=ioerr)
  
  return
  
end subroutine CS_WRITE_F2DEC
!**********************************************************************!



!**********************************************************************!
!
! CS_PREPARE_SCATTAMPS
!
! subroutine, prepares scattering amplitudes for each atom type in a
!             supercell.
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - if nabs==2, set CS_absorptionprm to the desired absorption factor
! - oversampling depends on CS_scaf_oversamp
!
! MODIFIES:
! - CS_scattamp_prepared        = preparation flag, sets TRUE
!                                 if previously true, reallocates memory
!                                    and erases old data
! - CS_scampnum                 = number of different atom types (Z,DW)
! - CS_scampptr                 = pointer connecting atoms in super cell with
!                                 scattering amplitude data
!                                 via indices in CS_scampdat, CS_scampatn, CS_scampdwf
! - CS_scampsym                 = atomic symbols for each atom pointer
! - CS_scampatn                 = atomic numbers for each atom pointer
! - CS_scampcrg                 = atomic charges for each atom pointer
! - CS_scampdwf                 = Debye-Waller factor for each atom pointer
!! - CS_scampgmax                = max. scattering vector to consider [1/nm]
! - CS_scampdat                 = scattering amplitude data for each atom pointer
! - CS_scadim                   = dimension of CS_scampdat
! - CS_scaitog                  = sampling rate of CS_scampdat
! - CS_scaff2d                  = 2d scattering amplitude data for each atom pointer
!
! INPUT:
!   integer*4 :: nx, ny, nz     = supercell axes discretization
!   integer*4 :: ndwf           = debye-waller factor usage flag
!   integer*4 :: nabs           = absorption flag
!   real*4 :: wl                = electron wavelength [nm]
!
! IN/OUTPUT:
!   integer*4 :: nerr           = error flag
subroutine CS_PREPARE_SCATTAMPS(nx,ny,nz,ndwf,nabs,wl,nerr)

  use FSCATAB
  
  implicit none
  
  integer*4, intent(in) :: nx, ny, nz, ndwf, nabs
  real*4, intent(in) :: wl
  integer*4, intent(inout) :: nerr

  logical :: dwfflg                 ! dwf-flag  
  logical :: absflg                 ! absorption flag
  integer*4 :: i, j, k, i1, j1      ! iterators
  integer*4 :: z                    ! atomic number
  integer*4 :: nyqx, nyqy, nyqz     ! nyqist numbers
  integer*4 :: scaos                ! applied oversampling rate
  real*4 :: dwf                     ! temp DWF
  real*4 :: dz                      ! temp ionic charge
  real*4 :: sx, sy, sz, itogx, itogy, itogz ! sampling parameters
  real*4 :: gx, gy, gx2, g2         ! reciprocal space coordinates
  real*4 :: g, scgmax, gmax         ! g and max. g of prepared scattering table
  real*4 :: ht                      ! high tension [kV]
  real*4 :: ftmp                    ! max. scattering amplitude
  real*4 :: rc                      ! relativistic correction
  real*4 :: apxy                    ! aperture factor for 2d form factors
  real*4 :: pfacio                  ! prefactor for ionic form factors
  complex*8 :: cval, cval0, csf     ! complex tmp vars
  character(len=400) :: smsg        ! message string
  character(len=CS_atsl) :: sym     ! atom symbol
  
  real*4, external :: wekoscar
  complex*8, external :: wekosca    !(g,dw,z,accvlt,absflg) from wekoscatt.f90,
                                    !    calculates scatt amp complex
  real*4, external :: wakiscar
  complex*8, external :: wakisca    !(g,dw,z,accvlt,absflg) from wakiscatt.f90,
                                    !    calculates scatt amp complex
  
  !
  ! preset error flag
  !
  nerr = 0
  dwfflg = .FALSE.
  if (ndwf/=0) dwfflg = .TRUE.
  absflg = .FALSE.
  CS_useabsorption = 0
  if (nabs/=0) then 
    absflg = .TRUE.
    CS_useabsorption = nabs
  end if
  scaos = min(max(10, CS_scaf_oversamp),100) ! limits the applied oversampling rate between 10 and 100
  !
  ! check cell data presence
  !
  if (.not.CS_cellmem_allocated) goto 13
  if (CS_numat<1) goto 12
  
  !
  ! release previous preparation
  !
  call CS_DEALLOC_SCAMPMEM(nerr)
  if (nerr/=0) goto 14
  
  !
  ! allocate new memory (required in case that supercell is changed)
  !
  allocate(CS_scampptr(CS_numat),stat=nerr)
  if (nerr/=0) goto 15
  allocate(CS_scampsym(CS_numat),stat=nerr)
  if (nerr/=0) goto 15
  allocate(CS_scampatn(CS_numat),stat=nerr)
  if (nerr/=0) goto 15
  allocate(CS_scampcrg(CS_numat),stat=nerr)
  if (nerr/=0) goto 15
  allocate(CS_scampdwf(CS_numat),stat=nerr)
  if (nerr/=0) goto 15
  
  !
  ! count atom types and fill pointer list
  !
  CS_scampnum = 0   ! clear atom type number
  CS_scampsym = ""  ! clear atom symbol
  CS_scampptr = 0   ! clear atom pointer
  CS_scampatn = 0   ! clear atomic numbers
  CS_scampcrg = 0.0 ! clear ionic charge
  CS_scampdwf = 0.0 ! clear atom dwf
  dwf = 0.0
  
  do i=1, CS_numat ! loop i over all atoms in supercell
    sym = CS_attype(i)
    z = CS_atnum(i)               ! get atom number / type
    dz = CS_atcrg(i)
    dwf = 0.0
    if (dwfflg) dwf = CS_atdwf(i) ! get atom dwf
    k = 0
    if (CS_scampnum>0) then
      do j=1, CS_scampnum ! loop j over known atom types
        if ((trim(CS_scampsym(j))==trim(sym)).and.(CS_scampdwf(j)==dwf)) then ! atom is known
          k = j
          exit  ! loop j
        end if
      end do ! loop j over known atom types
    end if
    if (k==0) then ! current atom is not known
      CS_scampnum = CS_scampnum + 1
      k = CS_scampnum
      CS_scampsym(k) = sym
      CS_scampatn(k) = z
      CS_scampcrg(k) = dz
      CS_scampdwf(k) = dwf
    end if
    CS_scampptr(i) = k ! set pointer
  end do ! loop i over all atoms in supercell
  if (CS_scampnum==0) goto 16
  write(unit=smsg,fmt='(A,I6,A)') "Found ",CS_scampnum," atomic individuals to prepare."
  call CS_MESSAGE(trim(smsg))
  
  !
  ! allocate memory for scattering amplitudes
  ! !!! This is Fourier space data, thus scrambled and transposed
  !
  ! - number of values used to sample the potentials
  !   using a much finer sampling than required since we keep them
  !   in 1d although they are used in 2d or 3d later
  CS_scadim = scaos*max(nz,ny,nx)
  allocate(CS_scampdat(CS_scadim,CS_scampnum),stat=nerr)
  if (nerr/=0) goto 15
  allocate(CS_scamptmp(16,CS_scampnum), stat=nerr)
  if (nerr/=0) goto 15
  CS_scamptmp = 0.0
  
  if (CS_useppot==1) then
    allocate(CS_scaff2d(ny,nx,CS_scampnum), stat=nerr)
    if (nerr/=0) goto 17
    CS_scaff2d = cmplx(0.,0.)
    allocate(CS_scagn(ny,nx,3), stat=nerr)
    if (nerr/=0) goto 17
    CS_scagn = cmplx(0.,0.)
    if (dwfflg) then
      allocate(CS_scadwf(ny,nx,CS_scampnum), stat=nerr)
      if (nerr/=0) goto 17
      CS_scadwf = 1.0
    end if
  end if
  
  !
  ! prepare parameters
  !
  cval0 = cmplx(0.0,0.0)
  cval = cval0
  ht = CS_WL2HT(wl)                                 ! high tension in kV
  rc = (CS_elm0+ht)/CS_elm0                         ! relativistic correction
  nyqx = (nx - modulo(nx,2))/2                      ! potentials nyquist number x
  nyqy = (ny - modulo(ny,2))/2                      ! potentials nyquist number y
  nyqz = (nz - modulo(nz,2))/2                      ! potentials nyquist number z
  sx = CS_scsx/real(nx)                             ! real space sampling x
  sy = CS_scsy/real(ny)                             ! real space sampling y
  sz = CS_scsz/real(nz)                             ! real space sampling z
  itogx = 1.0/CS_scsx                               ! fourier space sampling x
  itogy = 1.0/CS_scsy                               ! fourier space sampling y
  itogz = 1.0/CS_scsz                               ! fourier space sampling z
  ! Get 1.5 times highest scattering angle = 3/2 * n/2 * itog = 3/4 * itog * n = 0.75*itog*n
  scgmax = 0.75*max(real(nx)*itogx,real(ny)*itogy,real(nz)*itogz)
  CS_scaitog = scgmax/real(CS_scadim)               ! sampling rate of scattering data
  ! calculate gmax for the 2d form factors as the min. of the two gmax in the reciprocal plane
  gmax = min(itogx*real(nyqx),itogy*real(nyqy))
  
  !
  ! clear fourier space working field
  !
  CS_scampdat = cval0  
  
  !
  ! fill imaginary long calulation data now, but only if absflg is set
  !
  select case (nabs)
  
  case (2) ! alternative msa absorption model: HASHIMOTO, HOWIE, WHELAN, Proc. R. Soc. London Ser. A, 269, 80-103 (1962)   
  
    call CS_MESSAGE("Calculating real scattering amplitudes with proportional absorption.")
    
    do k=1, CS_scampnum ! loop k over atom-type lists
    
      sym = CS_scampsym(k)
      z = CS_scampatn(k)    ! get atomic number
      !dz = CS_scampcrg(k)
      dwf = 0.0
      if (dwfflg) dwf = CS_scampdwf(k)  ! get Debye-Waller factor
          
      if (CS_doconsolemsg>0) then
        write(unit=smsg,fmt='(A,I3,A,F12.8,A)') "Preparing data for Z = ",z," and DW-prm = ",dwf,"."
        call CS_MESSAGE(trim(smsg))
      end if
      
      !call CS_PROG_START(CS_scadim,1.0)
      
      do i=1, CS_scadim ! loop g-axis
        
        g = CS_scaitog*real(i-1)
        ftmp = 0.0
        j = 0
        if (CS_useextsca==1) j = FST_GETIDX(sym)
        if (j>0) then
          ftmp = FST_GETSCAR(sym, g, dwf, ht, dwfflg )
        else
          if (CS_scaf_table == 2) then
            ftmp = wakiscar(g, dwf, z, ht, dwfflg)  ! get scattering amplitude from Wa&Ki's table
          else
            ftmp = wekoscar(g, dwf, z, ht, dwfflg) ! get scattering amplitude from default table
          end if
        end if
        CS_scampdat(i,k) = cmplx(ftmp,CS_absorptionprm*ftmp) ! save scattering amplitude as real and absorption part as imag
        
        ! update progress indicator
        !call CS_PROG_UPDATE(i)
        
      end do ! loop g-axis
      
      !call CS_PROG_STOP(CS_scadim)
    
    end do ! loop k over atom-type lists
  
  case (1) ! standard msa absorption model:
           ! - WEICKENMEYER&KOHL, Acta Cryst. A 47 (1991) p. 597
           ! or
           ! - numerical integration according to HALL&HIRSCH Proc.Roy.Soc. A 286 (1965) p165 eq. (14a+b)
  
    call CS_MESSAGE("Calculating real and imaginary scattering amplitudes.")
  
    do k=1, CS_scampnum ! loop k over atom-type lists
    
      sym = CS_scampsym(k)
      z = CS_scampatn(k)    ! get atomic number
      !dz = CS_scampcrg(k)
      dwf = 0.0
      if (dwfflg) dwf = CS_scampdwf(k)  ! get Debye-Waller factor
      
      if (CS_doconsolemsg>0) then
        write(unit=smsg,fmt='(A,I3,A,F12.8,A)') "Preparing data for Z = ",z," and DW-prm = ",dwf,"."
        call CS_MESSAGE(trim(smsg))
      end if
      
      !call CS_PROG_START(CS_scadim,1.0)
    
      do i=1, CS_scadim ! loop g-axis
          
        g = CS_scaitog*real(i-1)
        j = 0
        if (CS_useextsca==1) j = FST_GETIDX(sym)
        if (j>0) then
          csf = FST_GETSCAC(sym, g, dwf, ht, dwfflg, absflg )
        else
          if (CS_scaf_table == 2) then
            csf = wakisca(g, dwf, z, ht, dwfflg, absflg)  ! get scattering amplitude from Wa&Ki's table
          else
            csf = wekosca(g, dwf, z, ht, dwfflg, absflg) ! get scattering amplitude from default table
          end if
        end if
        CS_scampdat(i,k) = csf              ! save scattering amplitude
        
        ! update progress indicator
        !call CS_PROG_UPDATE(i)
        
      end do ! loop g-axis
      
      !call CS_PROG_STOP(CS_scadim)
    
    end do ! loop k over atom-type lists
  
  case (0) ! no absorption, use faster fscatt version
  
    !call CS_MESSAGE("Calculating real scattering amplitudes only, no absorption.")
  
    do k=1, CS_scampnum ! loop k over atom-type lists

      sym = CS_scampsym(k)
      z = CS_scampatn(k)    ! get atomic number
      !dz = CS_scampcrg(k)
      dwf = 0.0
      if (dwfflg) dwf = CS_scampdwf(k)  ! get Debye-Waller factor
      
      if (CS_doconsolemsg>0) then
        write(unit=smsg,fmt='(A,I3,A,F12.8,A)') "Preparing data for Z = ",z," and DW-prm = ",dwf,"."
        call CS_MESSAGE(trim(smsg))
      end if
      
      ! call CS_PROG_START(CS_scadim,1.0)
    
      do i=1, CS_scadim ! loop g-axis
       
        g = CS_scaitog*real(i-1)
        ftmp = 0.0
        j = 0
        if (CS_useextsca==1) j = FST_GETIDX(sym)
        if (j>0) then
          ftmp = FST_GETSCAR(sym, g, dwf, ht, dwfflg )
        else
          if (CS_scaf_table == 2) then
            ftmp = wakiscar(g, dwf, z, ht, dwfflg) ! get scattering amplitude from Wa&Ki's table
          else
            ftmp = wekoscar(g, dwf, z, ht, dwfflg) ! get scattering amplitude from default table
          end if
        end if
        CS_scampdat(i,k) = cmplx(ftmp,0.0)    ! save scattering amplitude
        
        ! update progress indicator
        !call CS_PROG_UPDATE(i)
        
      end do ! loop g-axis
      
      !call CS_PROG_STOP(CS_scadim)
    
    end do ! loop k over atom-type lists
  
  end select ! case (nabs)
  
  
12 CS_scattamp_prepared = .TRUE.


  if (CS_useppot) then
    if (CS_doconsolemsg>0) call CS_MESSAGE("Preparing 2d form factor arrays.")
    pfacio = rc * CS_scaprea * CS_fpi * 10. ! prefactor for the ionic rest charge potential: rel-corr * C * 4 Pi * 10. [-> nm^-1]
    ! prepare all form factors and the shift kernels in fourier space
    ! (scrambled and transposed)
    apxy = 1.0
    do j=1, nx
      j1 = modulo(j-1+nyqx,nx) - nyqx ! frequency index
      gx = itogx*real(j1) ! 
      gx2 = gx*gx
      do i=1, ny
        i1 = modulo(i-1+nyqy,ny) - nyqy ! frequency index
        gy = itogy*real(i1) !
        g2 = gx2+gy*gy
        g = sqrt(g2)
        if (CS_do_bwl_pot) apxy = 0.5 - 0.5*tanh( (g/gmax-0.9)*30.0 )
        CS_scagn(i,j,1) = cmplx(0., -CS_tpi*gx) ! store -2*Pi*I*gx(i,j) ! for x-translations
        CS_scagn(i,j,2) = cmplx(0., -CS_tpi*gy) ! store -2*Pi*I*gy(i,j) ! for y-translations
        CS_scagn(i,j,3) = apxy*cmplx(pfacio/(0.25*g2+CS_iyr*CS_iyr), 0.) ! store rc*C1/(s^2+al^2)(i,j) ! for ionic potentials
        do k=1, CS_scampnum ! loop k over atom-type lists
          call CS_GETSCATTAMP(csf,g,k)
          CS_scaff2d(i,j,k) = csf*apxy ! store f_gz=0(i,j,atom type)
          if (dwfflg) then
            dwf = CS_scampdwf(k)
            CS_scadwf(i,j,k) = real( exp(-0.25*dwf*g2) , kind=4 ) ! Debye-Waller factor map
          end if
        end do
      end do
    end do
  end if
  
  return
  
  !
  ! handle errors
  !
13 nerr=-1
  call CS_ERROR("Cell data not ready.")
  return
14 nerr=-2
  call CS_ERROR("Memory deallocation failed.")
  return
15 nerr=-3
  call CS_ERROR("Memory allocation failed.")
  return
16 nerr=-4
  call CS_ERROR("There are no atoms?")
  return
17 nerr=-5
  call CS_ERROR("Failed to allocate memory for projected form fractors")
  return
  
end subroutine CS_PREPARE_SCATTAMPS
!**********************************************************************!

!**********************************************************************!
!
! CS_GET_MEANINNERPOT
!
! subroutine, calculates mean inner potentials of the input structure
!             from prepared scattering amplitudes (real part)
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
!
! IN/OUTPUT:
!   real*4 :: meanpot           = mean inner potential in volts
!   real*4 :: ht                = electron kinetic energy in keV
!   integer*4 :: nerr           = error flag
subroutine CS_GET_MEANINNERPOT(meanpot, ht, ierr)

!  use celslcprm
  
  implicit none
  
  real*4, intent(inout) :: meanpot
  real*4, intent(in)    :: ht
  integer*4, intent(inout) :: ierr

  real*4, parameter :: e0 = 510.9989461 ! electron rest mass in keV
  
  integer*4 :: i, j                 ! iterators
  complex*8 :: cf0                  ! DC form factor
  real*4 :: sf0                     ! sum of DC form factors (real parts)
  real*4 :: focc                    ! atom occupancy
  real*4 :: cellvol                 ! cell volume in nm^3
  real*4 :: pfac                    ! potential scale
  real*4 :: rc                      ! relativistic correction factor
  
!  character(len=400) :: smsg        ! message string
  
  !
  ! check cell data presence
  !
  if (.not.CS_cellmem_allocated) goto 13
  if (.not.(allocated(CS_scampdat).and.allocated(CS_scampptr))) goto 14
  if (CS_numat<1.or.CS_scampnum<1) goto 12
  !
  ! init
  !
  ierr = 0
  sf0 = 0.0
  focc = 1.0
  rc = (ht + e0) / e0 ! relativistic correction
  call CS_GETCELL_VOL(cellvol)
  !      cellvol is the supercell volume in nm^3
  !      
  pfac = CS_v0 / cellvol ! [eV * nm^2] / [ nm^3 ] = eV / nm
  !      CS_v0 is the conversion from fE to U
  !
  do i=1, CS_numat ! loop i over all atoms in supercell
    j = CS_scampptr(i) ! scattering factor index
    call CS_GETSCATTAMP(cf0, 0.0, j) ! get f(0)
    focc = CS_atocc(i) ! get atom occupancy
    sf0 = sf0 + real(cf0)*focc ! sum up (f0*occ)
  end do ! loop i over all atoms in supercell
  !
  ! sf0 is in nm !
  !
  ! scale to volts and remove relativistic correction
  meanpot = pfac * sf0 / rc ! > eV  
  
12 continue

  return  
  !
  ! handle errors
  !
13 ierr=-1
  call CS_ERROR("Cell data not ready.")
  return
14 ierr=-2
  call CS_ERROR("Scattering factors not prepared.")
  return
  
end subroutine CS_GET_MEANINNERPOT
!**********************************************************************!

!**********************************************************************!
!
! CS_SETSLICE_AUTO
!
! subroutine Partitioning of the super-cell into slices reproducing the
!            atomic planes of the structure. Empty space at the beginning
!            of the super cell will be reproduced by an empty slice.
!            Allocates slice arrays and assigns atoms to each slice.
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! 
! modifies:
! - CS_slczlim, CS_slcatacc, and CS_slcatnum
!
! INPUT:
!   integer*4 :: nrev = reverse slicing sequence (flag)
!
! IN/OUTPUT:
!   integer*4 :: nerr = error code
!
subroutine CS_SETSLICE_AUTO(nrev,nerr)
  
  implicit none
  
  integer*4, intent(in) :: nrev
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nalloc ! allocation status
  integer*4 :: i, j, k ! iterators
  integer*4 :: nslc ! number of slices to generate
  integer*4 :: iminl(1) ! minimum location
  real*4 :: z ! z offset of slice
  real*4 :: slcnumscal ! rescale relative to slice number
  real*4 :: znext ! z position of the next atom in the sequence
  real*4 :: ztol ! tolerated relative coordinate distance
  real*4 :: dist, dcell(3) ! distance calculation helpers
  real*4, allocatable :: slcset(:,:), atz(:) ! slice setting and atom z coordinates
  
  !
  ! preset error flags
  !
  nerr = 0
  nalloc = 0
  
  !
  ! check allocation state of super cell data memory
  ! we need the super-cell to do the slicing
  !
  if (.not.CS_cellmem_allocated) then
    nerr = -1
    return
  end if
  
  !
  ! auto-slicing ...
  !
  ztol = CS_slcztol / CS_scsz ! get tolerate z distance for atomic planes in fractional cell coordinates
  if (CS_numat>0) then ! ... if we have atoms
    allocate(slcset(2,CS_numat), stat=nalloc)
    if (nalloc/=0) goto 13
    slcset = 0.0 ! preset with starting coordinate
    nslc = 1 ! start with one slice
    allocate(atz(CS_numat), stat=nalloc)
    if (nalloc/=0) goto 13
    atz = 0.0
    do i=1, CS_numat ! copy atom z coordinates
      if (nrev==0) then ! normal slining order (0 ... 1)
        atz(i) = CS_atpos(3,i)
      else ! reversed slicing order (1 ... 0)
        atz(i) = 1.0 - CS_atpos(3,i)
      end if
    end do
    !
    ! determine the super-cell partitions (destroys atz)
    !
    do i=1, CS_numat
      iminl = minloc(atz)
      znext = atz(iminl(1)) ! this is the next atom z position
      atz(iminl(1)) = 2.0 ! invalidate to not trigger again
      if (znext>slcset(2,nslc)+ztol) then ! the next atom is well beyond the tolerated distance
        ! terminate the current slice
        slcset(2,nslc) = znext
        if (znext<1.0) then ! initiate a new slice at this atom position
          nslc = nslc + 1
          slcset(1:2,nslc) = znext
        else
          exit ! znext == 1 -> terminate
        end if
      end if
    end do
    ! terminate the final slice to max. cell coordinate
    slcset(2,nslc) = 1.0
    deallocate(atz, stat=nalloc) ! don't need atz anymore
    !
  else ! ... for an empty slice
    allocate(slcset(2,1), stat=nalloc)
    if (nalloc/=0) goto 13
    nslc = 1 ! just one empty slice
    slcset(1,1) = 0.0
    slcset(2,1) = 1.0
  end if
  
  !
  ! check allocation state of slice memory
  !
  if (CS_slicemem_allocated) then
    ! is allocated, thus deallocate to set new data
    call CS_DEALLOC_SLICEMEM(nerr)
    if (nerr/=0) goto 14
  end if
  
  !
  ! allocate slice data memory reflecting new slice number in size
  !
  call CS_ALLOC_SLICEMEM(nslc,nerr)
  if (nerr/=0) goto 13
  
  !
  ! set z-limits for each slice
  !
  do j=1, nslc ! loop for each slice
    CS_slczlim(1,j) = slcset(1,j)*CS_scsz
    CS_slczlim(2,j) = slcset(2,j)*CS_scsz
  end do
  
  !
  ! link atoms to slices (loop once over all atoms of the supercell)
  !
  slcnumscal = real(nslc)
  CS_atlnk = 0 ! reset atom link list
  dcell = (/ CS_scsz, CS_scsy, CS_scsz/)
  do j=1, CS_numat
    z = CS_atpos(3,j)
    if (nrev/=0) z = 1.0 - CS_atpos(3,j) ! reverse order solved by mirroring the SC at z=0.5
    k = nslc ! default slice is final slice, this includes z = 1.0, which is not handled by the search loop below
             ! although we initially used modulo to all coordinates, z=1 might appear here for the reverse slicing scheme.
    do i=1, nslc ! loop i over slice indices
      if (z>=slcset(1,i) .and. z<slcset(2,i)) then ! the z position of the atom is within bounds [z_i,z_i+1[ of slice i
        k = i ! set target slice index
        exit
      end if
    end do
    ! assign atom j to slice k
    CS_slcatnum(k) = CS_slcatnum(k) + 1 ! raise atom count of slice
    CS_slcatacc(CS_slcatnum(k),k) = j ! save atom index in slice list
    !
    k = 0
    if (j>1) then
      ! check for identical atomic site in the preceeding part of the atom list
      ! use CS_GET_SAMESITE( fpos,            l_fpos,              &
      !    &                 scal,  ssthr,          idxss, rdss, nerr )
      call  CS_GET_SAMESITE( CS_atpos(1:3,j), CS_atpos(1:3,1:j-1), &
           &                 dcell, CS_atdst_thrnm, k,     dist, nerr )
      if (k>0) CS_atlnk(j) = k
    end if
  end do
  
  !
  ! done
  !
10 continue
  if (allocated(slcset)) deallocate(slcset,stat=nalloc)
  if (allocated(atz)) deallocate(atz,stat=nalloc)
  return
  
13 nerr=-1
  call CS_ERROR("Failed to allocate slice data memory.")
  goto 10
14 nerr=-2
  call CS_ERROR("Failed to deallocate slice data memory.")
  goto 10
  
  return
  
end subroutine CS_SETSLICE_AUTO
!**********************************************************************!

!**********************************************************************!
!
! CS_SETSLICE_EQUIDIST
!
! subroutine distributing nslc slices equidistantly in supercell
!            fills slice data memory with numbers :o)
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! 
! modifies:
! - CS_slczlim, CS_slcatacc, and CS_slcatnum
!
! INPUT:
!   integer*4 :: nslc = number of slices in supercell
!   integer*4 :: nrev = reverse slicing sequence (flag)
!
! IN/OUTPUT:
!   integer*4 :: nerr = error code
!
subroutine CS_SETSLICE_EQUIDIST(nslc,nrev,nerr)
  
  implicit none
  
  integer*4, intent(in) :: nslc, nrev
  integer*4, intent(inout) :: nerr
  
  integer*4 :: i, j, k ! iterators
  real*4 :: z ! z offset of slice
  real*4 :: dz ! z increment per slice
  real*4 :: slcnumscal ! rescale relative to slice number
  real*4 :: dcell(3), dist
  
  !
  ! preset error flag
  !
  nerr = 0
  
  !
  ! preset slicing start coordinates
  !
  z = 0.0 ! start coordinate
  dz = CS_scsz/real(nslc) ! set the step size per slice
  
  !
  ! check allocation state of super cell data memory
  !
  if (.not.CS_cellmem_allocated) then
    nerr = -1
    return
  end if
  
  !
  ! check allocation state of slice memory
  !
  if (CS_slicemem_allocated) then
    ! is allocated, thus deallocate to set new data
    call CS_DEALLOC_SLICEMEM(nerr)
    if (nerr/=0) goto 14
  end if
  
  !
  ! allocate slice data memory reflecting new slice number in size
  !
  call CS_ALLOC_SLICEMEM(nslc,nerr)
  if (nerr/=0) goto 13
  
  !
  ! set z-limits for each slice
  !
  do j=1, nslc ! loop for each slice
    CS_slczlim(1,j) = real(j-1)*CS_scsz/real(nslc)
    CS_slczlim(2,j) = real(j)*CS_scsz/real(nslc)
  end do
  
  !
  ! link atoms to slices (loop once over all atoms of the supercell)
  !
  slcnumscal = real(nslc)
  CS_atlnk = 0 ! reset atom link list
  dcell = (/ CS_scsx, CS_scsy, CS_scsz/)
  do j=1, CS_numat
    z = CS_atpos(3,j)
    if (nrev/=0) z = 1.0 - CS_atpos(3,j) ! reverse order solved by mirroring the SC at z=0.5
    i = int(slcnumscal*z)+1 ! convert relative position to slice number
    if (i<1) i = 1 ! limit lower slice number to 1
    if (i>nslc) i = nslc ! limit upper slice number to max slice number
    CS_slcatnum(i) = CS_slcatnum(i) + 1 ! raise atom count of slice
    CS_slcatacc(CS_slcatnum(i),i) = j ! save atom index in slice list
    !
    k = 0
    if (j>1) then
      ! check for identical atomic site in the preceeding part of the atom list
      ! use CS_GET_SAMESITE( fpos,            l_fpos,              &
      !    &                 scal,  ssthr,          idxss, rdss, nerr )
      call  CS_GET_SAMESITE( CS_atpos(1:3,j), CS_atpos(1:3,1:j-1), &
           &                 dcell, CS_atdst_thrnm, k,     dist, nerr )
      if (k>0) CS_atlnk(j) = k
    end if
  end do
  
  !
  ! done
  !
  return
  
13 nerr=-1
  call CS_ERROR("Failed to allocate slice data memory.")
  return
14 nerr=-2
  call CS_ERROR("Failed to deallocate slice data memory.")
  return
  
end subroutine CS_SETSLICE_EQUIDIST
!**********************************************************************!

!**********************************************************************!
!
! CS_GETSLICETITLE
!
! browses through the slice dependent atom lists and creates
! a short content string, like O13 Sr2
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - slice data data to be allocated and set up
!   (e.g. call CS_SETSLICE_EQUIDIST)
! 
! INPUT:
!   integer*4 :: nslc                   = slice index
!
! IN/OUTPUT:
!   character(len=*) :: stitle          = title string for output
!   integer*4 :: nerr                   = error code
!
subroutine CS_GETSLICETITLE(nslc, stitle, nerr)

  implicit none
  
  integer*4, intent(in) :: nslc
  integer*4, intent(inout) :: nerr
  character(len=*), intent(inout) :: stitle
  
  integer*4 :: na                       ! number of atoms in requested slice
  integer*4 :: i, j, ia, k, l, m        ! iterators
  integer*4 :: occ_int, occ_frc         ! occupancy count integer value and fractional value
  integer*4 :: nalloc                   ! allocation status
  integer*4 :: ncomp                    ! number of components of the composit list
  integer*4, allocatable :: composit(:) ! composition atoms
  real*4 :: occ                         ! occupancy of current atom
  real*4, allocatable :: compocc(:)     ! composition atom-wise occupancy
  character(len=CS_atsl), allocatable :: compsym(:)  ! composition atom-wise symbol
  character(len=CS_atsl) :: symbol                 ! atom symbol
  character(len=1024) :: stmp           ! temp. string
  
  !  external :: getsymb                   ! link fscatt.f

  !
  ! preset error flag
  !
  nerr = 0
  
  !
  ! preset number of components to zero
  !
  ncomp = 0
  
  !
  !
  !
  if (3>len(stitle)) then ! string too short
    goto 17
  end if
  
  !
  ! check mem alloc
  !
  if (.not.(CS_slicemem_allocated.and.CS_cellmem_allocated)) goto 13
  
  !
  ! check slice index
  !
  if (nslc<1.or.nslc>CS_nspsc) goto 14
  
  !
  ! clear up the string
  !
  stitle = repeat(char(0),len(stitle))
  
  !
  ! check for slice content
  !
  na = CS_slcatnum(nslc)
  if (na<=0) then ! nothing in slice
    stitle(1:11) = "Empty slice"
    goto 12
  else
  
    ! there is a number of atoms, count the atom types
    allocate(composit(na), compocc(na), compsym(na), stat=nalloc)
    if (nalloc/=0) goto 15
    composit = 0 ! preset composition counter to zero
    compocc = 0.0 ! preset occupancy info to zero
    compsym = "    "
  
    !
    ! iterate through the atoms of a slice
    !
    do i=1, na
      ! get atom index in supercell data
      k = CS_slcatacc(i,nslc)
      ia = CS_atnum(k)
      symbol = CS_attype(k)
      occ = CS_atocc(k)
      ! count atom
      do j=1, na
        if (trim(compsym(j))==trim(symbol)) then ! atom type found already
          compocc(j) = compocc(j) + occ ! raise slice occupancy with this atom type
          exit
        end if
        if (len_trim(compsym(j))==0) then ! end of current list reached, new entry
          compsym(j) = symbol ! set new symbol
          composit(j) = ia ! add new atom type to composition list
          compocc(j) = occ ! register 1st atom of the new type with current occupancy
          ncomp = ncomp + 1 ! increase number of components
          exit
        end if
      end do
    end do
    
    !
    ! iterate through composition list and create the string
    !
    j = 1
    do i=1, ncomp
      ! get atom index in supercell data
      ia = composit(i)
      ! get total occupancy for slice
      occ = abs(compocc(i))
      ! get atom symbol
      symbol = compsym(i)
      ! create substring
      k = len_trim(symbol)
      occ_int = int(occ)
      occ_frc = int(100.0*(occ-real(occ_int)))
      l = 1
      if (occ_int>0) then ! get length of the integer part string
        l = 1 + int(log10(real(abs(occ_int))))
      end if
      if (occ_frc<1) then ! non-fractional occupancy
        write(unit=stmp,fmt=110) trim(symbol), occ_int
110     format(A<k>,"_",I<l>)
      else ! fractional occupancy
        write(unit=stmp,fmt=111) trim(symbol), occ_int, occ_frc
111     format(A<k>,"_",I<l>,".",I2.2) 
      end if
      ! add things to the string ...
      if (j>1) then ! add separator of this is not the first entry
        stitle(j:j) = ' '
        j = j + 1 ! move one character forth
        if (j>=len(stitle)-1) goto 11 ! stop, string full
      end if
      ! determine sub-string position
      ! j ist sub-string start position
      l = len_trim(stmp) ! sub-string trim length
      k = min(j + l-1, len(stitle)-1 ) ! sub-string end position, limited to the title string length
      m = 1 + k - j ! update of the sub-string length
      stitle(j:k) = stmp(1:m) ! copy
      j = k+1 ! update start position for next things to add ...
      if (j>=len(stitle)-1) goto 11 ! stop, string full
    end do
    
11  deallocate(composit, compocc, stat=nalloc)
    if (nalloc/=0) goto 16
    
    if ( (j>=len(stitle)-1) .and. (i<=ncomp) ) then
      ! the string writing was not finished. mark this by setting the last usable character to '+'
      k = len(stitle)-1
      stitle(k:k) = '+'
    end if
    
  end if
  
  
12 k = len(stitle)
   stitle(k:k) = char(0) ! for fool proofness, always set the last character to 0.

   return
  
13 nerr=1
  call CS_ERROR("Cell and slice memory not allocated.")
  return
  
14 nerr=2
  call CS_ERROR("Invalid parameter, slice index.")
  return
  
15 nerr=3
  call CS_ERROR("Memory allocation failed.")
  return
  
16 nerr=4
  call CS_ERROR("Memory deallocation failed.")
  return

17 nerr=5
  call CS_ERROR("Invalid parameter, title string length insufficient.")
  return

    
end subroutine CS_GETSLICETITLE
!**********************************************************************!

!**********************************************************************!
!
! CS_LOAD_EMSCELL
!
! subroutine loading EMS supercell data from file
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   integer*4 :: nerr         = error flag
!
subroutine CS_LOAD_EMSCELL(sfile, nerr)

  use FSCATAB
  use cifio
  
  implicit none
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lun
  integer*4 :: nlines ! number of lines in file
  integer*4 :: natoms ! number of atom definitions
  integer*4 :: ndummy, ntmp, ntmp2, nct
  integer*4 :: nread
  character(len=CS_ll) :: sline, smsg
  character(len=CS_atsl) :: atomid ! atom id string
  real*4 :: apx, apy, apz ! atom position
  real*4 :: aoc ! atom occupancy
  real*4 :: dwf ! debye waller factor
  real*8 :: crg ! charge
  
  !
  ! preset error flag
  !
  nerr = 0
  !
  ! get next free logical unit number
  !
  lun = CS_GETFREELUN(nerr)
  !
  ! check success
  !
  if (nerr/=0) then
    goto 13
  end if
  !
  ! open file for unformatted sequential reading
  !
  open( file=trim(sfile),unit=lun,iostat=nerr, &
     &  action='READ',status='OLD',err=14)
     
  !
  ! determine number of lines in file
  !
  nlines = 0
  do
    ndummy = nlines
    read(unit=lun,iostat=nerr,fmt='(A)') sline
    if (nerr==0) then
      nlines = nlines + 1
    end if
    if (sline(1:1)=='*') exit ! done reading, CEL finalizer found , EMS eof signature
    if (ndummy==nlines) exit
  end do
  rewind(unit=lun,iostat=nerr)
  natoms = nlines-3 ! substract tilte, cell definition, and EOF line
  CS_numat = natoms
  
  !
  ! allocate atom data array
  !
  if (CS_cellmem_allocated) then
    call CS_DEALLOC_CELLMEM(nerr)
  end if
  if (nerr==0) call CS_ALLOC_CELLMEM(CS_numat,nerr)
  if (nerr/=0) goto 19
  
  !
  ! read header line (unimportant blabla)
  !
  read(unit=lun,err=15,fmt='(A)',iostat=nerr) sline
  if (nerr/=0) goto 15
  !
  ! read supercell definition line
  !
  read(unit=lun,err=15,fmt='(A)',iostat=nerr) sline
  if (nerr/=0) goto 15
  !
  ! parse supercell definitions
  !
  read(unit=sline,err=17,fmt=*,iostat=nerr) ndummy, CS_scsx, CS_scsy, CS_scsz, &
     &                          CS_scalpha, CS_scbeta, CS_scgamma
  if (nerr/=0) goto 17
  !
  ! loop through atom data and read
  !
  nct = 0
  do
    !
    ! try read next line
    !
    read(unit=lun,fmt='(A)',err=15,iostat=nread) sline
    if (nread/=0) goto 15 ! done reading, eof reached or some error occurred
    if (sline(1:1)=='*') exit ! done reading, CEL finalizer found , EMS eof signature
    sline = trim(adjustl(sline)) ! remove leading and trailing spaces
    !
    ! parse atom data
    ! (loop until) find first space character which follows a non-space character
    !
    ndummy = len_trim(sline)
    ntmp = index(sline," ")
    if (ntmp==0) goto 18 ! no space character in line, this cannot be
    if (ntmp==1.and.ntmp<ndummy) then
      ! loop
      do
        ntmp2 = index(sline(ntmp+1:ndummy)," ")
        if (ntmp2==0) goto 18 ! no following space character found -> error
        if (ntmp2>1) then ! other characters follow before the next space character -> exit
           ntmp = ntmp + ntmp2 ! step there
           exit ! done
        end if
        ! still here -> thus next space char follows directly
        ntmp = ntmp + ntmp2 ! step one position right, try again
        if (ntmp>=ndummy) goto 18 ! end of string reached -> error
      end do
    end if
    ! positionen error checks
    if (ntmp==0) goto 18
    if (ntmp==1) goto 18
    if (ntmp>=ndummy) goto 18
    
    !
    ! copy string data to local variables
    !
    read(unit=sline(1:ntmp-1),err=18,iostat=nerr,fmt='(A)') atomid
    if (nerr/=0) goto 18
    read(unit=sline(ntmp:ndummy),err=18,iostat=nerr,fmt=*) apx, apy, apz, aoc, dwf
    if (nerr/=0) goto 18
    
    !
    ! transfer data to module memory
    !
    nct = nct + 1
    CS_attype(nct) = adjustl(atomid)
    CS_atpos(1,nct) = modulo(apx,1.0)
    CS_atpos(2,nct) = modulo(apy,1.0)
    CS_atpos(3,nct) = modulo(apz,1.0)
    if (apz<0.0 .or. apz>=1.0 .or. apy<0.0 .or. apy>=1.0 .or. apx<0.0 .or. apx>=1.0) then
      write(unit=smsg,fmt='(A,I4,A,A2,1x,3F10.6)') "Atom(",nct,"): ",adjustl(atomid), apx, apy, apz
      call CS_MESSAGE(trim(smsg))
      write(unit=smsg,fmt='(A,3F10.6)') "         ->    ", CS_atpos(1,nct), CS_atpos(2,nct), CS_atpos(3,nct)
      call CS_MESSAGE(trim(smsg))
    end if
    
    CS_atocc(nct) = aoc
    CS_atdwf(nct) = dwf
!    write(unit=smsg,fmt='(A,I4,A,A2,1x,5F10.6)') "Atom data (",nct,"): ",adjustl(atomid), apx, apy, apz, aoc, dwf
!    call CS_MESSAGE(trim(smsg))
    
    !
    ! identify atom, save atomic number Z
    ! ... this is for checking that we can recognize the input
    !
    ntmp = 0
    if (CS_useextsca==1) then ! try to identify atom in external list 1
      ntmp = FST_GETIDX(CS_attype(nct))
      if (ntmp>0) then
        CS_atnum(nct) = nint(FST_crg(2,ntmp))
        CS_atcrg(nct) = FST_crg(3,ntmp)
      end if
    end if
    if (ntmp<=0) then ! fall-back to default table
      crg = 0.0D+0
      call CIF_get_atnum(CS_attype(nct),ntmp,crg)
      if (ntmp<0) goto 20
      CS_atnum(nct) = ntmp
      CS_atcrg(nct) = real(crg, kind=4)
    end if
    
  end do
  
  !
  ! close
  !
  close(unit=lun,err=16,iostat=nerr)
  !
  ! return without error
  !
  return
  
  !
  ! handle errors here, call specific error messages
  !
13 continue
  nerr = -1
  call CS_ERROR("Failed to get free logical unit number.")
  return
14 continue
  nerr = -2
  call CS_ERROR("Failed to open cell file.")
  return
15 continue
  nerr = -3
  call CS_ERROR("Failed to read data from cell file.")
  close(unit=lun,err=16,iostat=nerr)
  return
16 continue
  nerr = -4
  call CS_ERROR("Failed to close cell file.")
  return
17 continue
  nerr=-5
  call CS_ERROR("Failed to parse cell definitions.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: 0 <size-x> <size-y> <size-z> <alpha> <beta> <gamma>")
  close(unit=lun,err=16,iostat=nerr)
  return
18 continue
  nerr=-5
  call CS_ERROR("Failed to parse atom definitions.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: <Symbol> <pos-x> <pos-y> <pos-z> <occup.> <DW-param.> <not used> <not used> <not used>")
  close(unit=lun,err=16,iostat=nerr)
  return
19 continue
  nerr=-6
  call CS_ERROR("Failed to allocate cell data memory.")
  return
20 continue
  nerr=-1
  call CS_ERROR("Failed to identify atom.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: <Symbol> <pos-x> <pos-y> <pos-z> <occup.> <DW-param.> <not used> <not used> <not used>")
  close(unit=lun,err=16,iostat=nerr)
  return
  
end subroutine CS_LOAD_EMSCELL
!**********************************************************************!


!**********************************************************************!
!
! CS_LOAD_EMSCELLINFO
!
! subroutine loading EMS supercell info data from file, does not store the data
! in global module parameters
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   real*4 :: sx, sy, sz      = cell dimensions in nm
!   integer*4 :: natom        = total number of atoms
!   integer*4 :: nerr         = error flag
!
subroutine CS_LOAD_EMSCELLINFO(sx, sy, sz, natom, sfile, nerr)
  
  implicit none
  
  real*4, intent(inout) :: sx, sy, sz
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: natom, nerr
  
  integer*4 :: lun
  integer*4 :: nlines ! number of lines in file
  integer*4 :: natoms ! number of atom definitions
  integer*4 :: ndummy
  real*4 :: aa, ab, ac
  character(len=CS_ll) :: sline
  
  !
  ! preset error flag
  !
  nerr = 0
  !
  ! get next free logical unit number
  !
  lun = CS_GETFREELUN(nerr)
  !
  ! check success
  !
  if (nerr/=0) then
    goto 13
  end if
  !
  ! open file for unformatted sequential reading
  !
  open( file=trim(sfile),unit=lun,iostat=nerr, &
     &  action='READ',status='OLD',err=14)
     
  !
  ! determine number of lines in file
  !
  nlines = 0
  do
    ndummy = nlines
    read(unit=lun,iostat=nerr,fmt='(A)') sline
    if (nerr==0) then
      nlines = nlines + 1
    end if
    if (sline(1:1)=='*') exit ! done reading, CEL finalizer found , EMS eof signature
    if (ndummy==nlines) exit
  end do
  rewind(unit=lun,iostat=nerr)
  natoms = nlines-3 ! substract tilte, cell definition, and EOF line
  
  
  !
  ! read header line (unimportant blabla)
  !
  read(unit=lun,err=15,fmt='(A)',iostat=nerr) sline
  if (nerr/=0) goto 15
  !
  ! read supercell definition line
  !
  read(unit=lun,err=15,fmt='(A)',iostat=nerr) sline
  if (nerr/=0) goto 15
  !
  ! parse supercell definitions
  !
  read(unit=sline,err=17,fmt=*,iostat=nerr) ndummy, sx, sy, sz, aa, ab, ac
  if (nerr/=0) goto 17
  !
  ! close
  !
  close(unit=lun,err=16,iostat=nerr)
  !
  ! return without error
  !
  
  natom = natoms
  
  return
  
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed to get free logical unit number.")
  return
14 nerr = -2
  call CS_ERROR("Failed to open cell file.")
  return
15 nerr = -3
  call CS_ERROR("Failed to read data from cell file.")
  close(unit=lun,err=16,iostat=nerr)
  return
16 nerr = -4
  call CS_ERROR("Failed to close cell file.")
  return
17 nerr=-5
  call CS_ERROR("Failed to parse cell definitions.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: 0 <size-x> <size-y> <size-z> <alpha> <beta> <gamma>")
  close(unit=lun,err=16,iostat=nerr)
  return
18 nerr=-5
  call CS_ERROR("Failed to parse atom definitions.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: <Symbol> <pos-x> <pos-y> <pos-z> <occup.> <DW-param.> <not used> <not used> <not used>")
  close(unit=lun,err=16,iostat=nerr)
  return
19 nerr=-6
  call CS_ERROR("Failed to allocate cell data memory.")
  return
20 nerr=-1
  call CS_ERROR("Failed to identify atom.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: <Symbol> <pos-x> <pos-y> <pos-z> <occup.> <DW-param.> <not used> <not used> <not used>")
  close(unit=lun,err=16,iostat=nerr)
  return
  
end subroutine CS_LOAD_EMSCELLINFO
!**********************************************************************!


!**********************************************************************!
!
! CS_LOAD_EMSCELLINFO2
!
! subroutine loading EMS supercell info data from file
! - does not store the data in global module parameters
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   real*4 :: dims(3)         = cell dimensions in nm
!   real*4 :: angs(3)         = cell angles in degrees
!   integer*4 :: natom        = total number of atoms
!   integer*4 :: nerr         = error flag
!
subroutine CS_LOAD_EMSCELLINFO2(dims, angs, natom, sfile, nerr)
  
  implicit none
  
  real*4, intent(inout) :: dims(3), angs(3)
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: natom, nerr
  
  integer*4 :: lun
  integer*4 :: nlines ! number of lines in file
  integer*4 :: natoms ! number of atom definitions
  integer*4 :: ndummy
  character(len=CS_ll) :: sline
  
  !
  ! preset error flag & output data
  !
  nerr = 0
  dims = 0.0
  angs = 0.0
  natom = 0
  !
  ! get next free logical unit number
  !
  lun = CS_GETFREELUN(nerr)
  !
  ! check success
  !
  if (nerr/=0) then
    goto 13
  end if
  !
  ! open file for unformatted sequential reading
  !
  open( file=trim(sfile),unit=lun,iostat=nerr, &
     &  action='READ',status='OLD',err=14)
     
  !
  ! determine number of lines in file
  !
  nlines = 0
  do
    ndummy = nlines
    read(unit=lun,iostat=nerr,fmt='(A)') sline
    if (nerr==0) then
      nlines = nlines + 1
    end if
    if (sline(1:1)=='*') exit ! done reading, CEL finalizer found , EMS eof signature
    if (ndummy==nlines) exit
  end do
  rewind(unit=lun,iostat=nerr)
  natoms = nlines-3 ! substract tilte, cell definition, and EOF line
  !
  ! read header line (unimportant blabla)
  !
  read(unit=lun,err=15,fmt='(A)',iostat=nerr) sline
  if (nerr/=0) goto 15
  !
  ! read supercell definition line
  !
  read(unit=lun,err=15,fmt='(A)',iostat=nerr) sline
  if (nerr/=0) goto 15
  !
  ! close file
  !
  close(unit=lun,err=16,iostat=nerr)
  !
  !
  ! parse supercell definitions (line 2 remains in sline)
  !
  ! unformatted read version (new and default form)
  read(unit=sline,fmt=*,iostat=nerr) ndummy, dims(1), dims(2), dims(3), angs(1), angs(2), angs(3)
  if (nerr/=0) then
    ! try formatted version (old EMS input)
    read(unit=sline,fmt=501,iostat=nerr) dims(1), dims(2), dims(3), angs(1), angs(2), angs(3)
    if (nerr/=0) then
      ! all reading attempts failed
      goto 17
    end if
  end if
  !
  
  natom = natoms ! set number of atoms return value
  
  return
  
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed to get free logical unit number.")
  return
14 nerr = -2
  call CS_ERROR("Failed to open cell file.")
  return
15 nerr = -3
  call CS_ERROR("Failed to read data from cell file.")
  close(unit=lun,err=16,iostat=nerr)
  return
16 nerr = -4
  call CS_ERROR("Failed to close cell file.")
  return
17 nerr=-5
  call CS_ERROR("Failed to parse cell definitions.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: 0 <size-x> <size-y> <size-z> <alpha> <beta> <gamma>")
  close(unit=lun,err=16,iostat=nerr)
  return
18 nerr=-5
  call CS_ERROR("Failed to parse atom definitions.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: <Symbol> <pos-x> <pos-y> <pos-z> <occup.> <DW-param.> <not used> <not used> <not used>")
  close(unit=lun,err=16,iostat=nerr)
  return
19 nerr=-6
  call CS_ERROR("Failed to allocate cell data memory.")
  return
20 nerr=-1
  call CS_ERROR("Failed to identify atom.")
  call CS_MESSAGE(" last line reads: "//trim(sline))
  call CS_MESSAGE(" expected format: <Symbol> <pos-x> <pos-y> <pos-z> <occup.> <DW-param.> <not used> <not used> <not used>")
  close(unit=lun,err=16,iostat=nerr)
  return
  
501 FORMAT(4x, 6f8.4) ! taken from EMS (P. Stadelmann, file sc5ems.f, line 196)
  
end subroutine CS_LOAD_EMSCELLINFO2
!**********************************************************************!


!**********************************************************************!
!
! CS_SAVE_EMSCELL
!
! subroutine saving supercell data to CEL file
! uses specific EMS file format
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   integer*4 :: nerr         = error flag
!
subroutine CS_SAVE_EMSCELL(sfile, nerr)
  
  implicit none
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lun, i
  integer*4 :: nlines ! number of lines in file
  
  external :: createfilefolder
  
  !
  ! preset error flag
  !
  nerr = 0
  !
  ! get next free logical unit number
  !
  lun = CS_GETFREELUN(nerr)
  !
  ! check success
  !
  if (nerr/=0) then
    goto 13
  end if
  !
  ! open file for unformatted sequential reading
  !
  call createfilefolder(trim(sfile),nerr)
  if (nerr/=0) then
    goto 18
  end if
  open( file=trim(sfile),unit=lun,iostat=nerr, &
     &  action='WRITE',status='REPLACE',err=14)
     
  !
  ! determine number of lines in file
  !
  nlines = CS_numat
  if (nlines<1) goto 16
  
  !
  ! write header line (unimportant blabla)
  !
  write(unit=lun,err=15,fmt='(A)') "# generated by CELSLC: "//trim(sfile)
  !
  ! write supercell definition line
  !
  !1001 FORMAT(1x,i2,a1,6f8.4)
  write(unit=lun,err=15,fmt=1001) 0," ",CS_scsx, CS_scsy, CS_scsz, &
     &                          CS_scalpha, CS_scbeta, CS_scgamma
  !
  ! loop through atom data and write
  !
  do i=1, nlines
    !
    ! try write next line
    !
    !1002 FORMAT(1x,a2,a1,8f8.4)
    write(unit=lun,fmt=1002,err=15) trim(CS_attype(i)), " ", &
     &                             CS_atpos(1,i), CS_atpos(2,i), CS_atpos(3,i), &
     &                             CS_atocc(i), CS_atdwf(i), 0.0, 0.0, 0.0
    
  end do
  !
  ! write closing line
  !
  write(unit=lun,err=15,fmt='(A)') "*"
  
  !
  ! close
  !
  close(unit=lun,err=17,iostat=nerr)
  !
  ! return without error
  !
  return
  
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed to acquire logical file unit.")
  return
  
14 nerr = -2
  call CS_ERROR("Failed to open file.")
  return

15 nerr = -3
  call CS_ERROR("Failed writing data to file.")
  close(unit=lun,err=17,iostat=nerr)
  return

16 nerr = -4
  call CS_ERROR("Sorry, no atom data.")
  return
  
17 nerr = -5
  call CS_ERROR("Failed to close cell file.")
  return
  
18 nerr = -6
  call CS_ERROR("Failed to create output path.")
  return

  
1001 FORMAT(1x,i2,a1,6f9.5)
1002 FORMAT(1x,a4,a1,8f9.6)
  
end subroutine CS_SAVE_EMSCELL
!**********************************************************************!





!**********************************************************************!
!
! CS_LOAD_CIFCELLINFO
!
! subroutine loading CIF supercell info data from file, does not store the data
! in global module parameters
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   real*4 :: sx, sy, sz      = cell dimensions in nm
!   integer*4 :: natom        = total number of atoms
!   integer*4 :: nerr         = error flag
!
subroutine CS_LOAD_CIFCELLINFO(sx, sy, sz, natom, sfile, nerr)
  
  use cifio
  
  implicit none
  
  real*4, intent(inout) :: sx, sy, sz
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: natom, nerr
  
  !
  ! preset error flag
  !
  nerr = 0
  !
  ! use cifio module for loading
  !
  CIF_verbosity = CS_doconsolemsg
  !
  call CIF_READ(trim(sfile), nerr)
  if (nerr<=0) goto 13 ! check success
  !
  !
  ! transfer the cell data to this module
  !   super cell size -> from A to nm
  sx = real(CIF_cell_length_a,kind=4)*0.1
  sy = real(CIF_cell_length_b,kind=4)*0.1
  sz = real(CIF_cell_length_c,kind=4)*0.1
  !   number of atoms per super cell
  natom = CIF_atom_site_number
  !
  return
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed reading CIF information from file ["//trim(sfile)//"].")
  return
  
end subroutine CS_LOAD_CIFCELLINFO
!**********************************************************************!



!**********************************************************************!
!
! CS_LOAD_CIFCELLINFO2
!
! subroutine loading CIF supercell info data from file
! - does not store the data in global module parameters
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   real*4 :: dims(3)         = cell dimensions in nm
!   real*4 :: angs(3)         = cell angles in degrees
!   integer*4 :: natom        = total number of atoms
!   integer*4 :: nerr         = error flag
!
subroutine CS_LOAD_CIFCELLINFO2(dims, angs, natom, sfile, nerr)
  
  use cifio
  
  implicit none
  
  real*4, intent(inout) :: dims(3), angs(3)
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: natom, nerr
  
  !
  ! preset error flag
  !
  nerr = 0
  !
  ! use cifio module for loading
  !
  CIF_verbosity = CS_doconsolemsg
  !
  call CIF_READ(trim(sfile), nerr)
  if (nerr<=0) goto 13 ! check success
  !
  !
  ! transfer the cell data to this module
  !   super cell size -> from A to nm
  dims(1) = real(CIF_cell_length_a,kind=4)*0.1
  dims(2) = real(CIF_cell_length_b,kind=4)*0.1
  dims(3) = real(CIF_cell_length_c,kind=4)*0.1
  !   supercell axis angles [deg]
  angs(1) = real(CIF_cell_angle_alpha,kind=4)
  angs(2) = real(CIF_cell_angle_beta ,kind=4)
  angs(3) = real(CIF_cell_angle_gamma,kind=4)
  !   number of atoms per super cell
  natom = CIF_atom_site_number
  !
  return
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed reading CIF information from file ["//trim(sfile)//"].")
  return
  
end subroutine CS_LOAD_CIFCELLINFO2
!**********************************************************************!





!**********************************************************************!
!
! CS_LOAD_CIFCELL
!
! subroutine loading CIF supercell data from file
! - stores the data in global module parameters
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   integer*4 :: nerr         = error flag
!
subroutine CS_LOAD_CIFCELL(sfile, nerr)
  
  use cifio
  
  implicit none
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  integer*4 :: i, ierr
  real(SELECTED_REAL_KIND(15,307)) :: rtmp
  !
  ! preset error flag
  !
  nerr = 0
  ierr = 0
  !
  ! use cifio module for loading
  !
  CIF_verbosity = CS_doconsolemsg
  !
  call CIF_READ(trim(sfile), ierr)
  if (ierr<=0) goto 13 ! check success
  !
  !
  ! transfer the cell data to this module
  !   super cell size -> from A to nm
  CS_scsx = real(CIF_cell_length_a,kind=4)*0.1
  CS_scsy = real(CIF_cell_length_b,kind=4)*0.1
  CS_scsz = real(CIF_cell_length_c,kind=4)*0.1
  !   supercell axis angles [deg]
  CS_scalpha = real(CIF_cell_angle_alpha,kind=4)
  CS_scbeta  = real(CIF_cell_angle_beta ,kind=4)
  CS_scgamma = real(CIF_cell_angle_gamma,kind=4)
  !   number of atoms per super cell
  CS_numat = CIF_atom_site_number
  !
  ! allocate atom data array
  ierr = 0
  if (CS_cellmem_allocated) then
    call CS_DEALLOC_CELLMEM(ierr)
  end if
  if (ierr==0) call CS_ALLOC_CELLMEM(CS_numat,ierr)
  if (ierr/=0) goto 19
  !
  if (CS_numat>0) then ! are there atomic sites ?
    ! transfer all data (apply unit conversions on that way)
    do i=1, CS_numat
      call CIF_get_atom_site_type_symbol(i,CS_attype(i))
      CS_atpos(1,i) = real(CIF_atom_site(1,i),kind=4)
      CS_atpos(2,i) = real(CIF_atom_site(2,i),kind=4)
      CS_atpos(3,i) = real(CIF_atom_site(3,i),kind=4)
      CS_atnum(i)   = nint(CIF_atom_site(4,i),kind=4)
      CS_atcrg(i)   = real(CIF_atom_site(5,i),kind=4)
      CS_atocc(i)   = real(CIF_atom_site(6,i),kind=4)
      call CIF_get_atom_site_biso(i,rtmp) ! get the biso value (A^2)
      CS_atdwf(i)   = real(rtmp,kind=4)*0.01 ! scale to (nm^2)
    end do
    !
  end if
  !
  ! free memory in the cifio module
  call CIF_CLEAR()
  ! free memory in the spacegroups module
  call SG_UNINIT()
  !
  ! return without error
  !
  return
  
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed read CIF structure data from ["//trim(sfile)//"].")
  return
19 nerr = -6
  call CS_ERROR("Failed to allocate cell data memory.")
  return
  
end subroutine CS_LOAD_CIFCELL
!**********************************************************************!




!**********************************************************************!
!
! CS_SAVE_CIFCELL
!
! subroutine writing super-cell data to file in CIF format
! - uses the structure data stored in CellSlice module
!
! INPUT:
!   character(len=*) :: sfile = disk file name
!
! IN/OUTPUT:
!   integer*4 :: nerr         = error flag
!
subroutine CS_SAVE_CIFCELL(sfile, nerr)
  
  use cifio
  
  implicit none
  
  character(len=*), intent(in) :: sfile
  integer*4, intent(inout) :: nerr
  integer*4 :: i, nalloc
  !
  ! preset error flag
  !
  nerr = 0
  if (.not.CS_cellmem_allocated) then
    return ! no data, just return
  end if
  !
  ! use cifio module for writing
  !
  ! transfer data to the cif module
  ! - initialize
  CIF_verbosity = CS_doconsolemsg
  ! transfer the cell data to this module
  !   super cell size -> from nm to A
  CIF_cell_length_a = dble(CS_scsx*10.0)
  CIF_cell_length_b = dble(CS_scsy*10.0)
  CIF_cell_length_c = dble(CS_scsz*10.0)
  !   supercell axis angles [deg]
  CIF_cell_angle_alpha = dble(CS_scalpha)
  CIF_cell_angle_beta  = dble(CS_scbeta )
  CIF_cell_angle_gamma = dble(CS_scgamma)
  !   symmetry group name and number
  CIF_symmetry_Int_Tables_number = 1
  CIF_symmetry_space_group_name_HM = "P 1"
  !   prepare atom site table
  CIF_atom_site_number = 0
  if (CS_numat>0) then
    allocate(CIF_atom_site(CIF_atom_site_items_num,CS_numat), stat=nalloc) ! allocate
    if (nalloc/=0) goto 19
    CIF_atom_site = 0.0d+0 ! initialize
    CIF_atom_site_number = CS_numat
    ! transfer all data (apply unit conversions on that way)
    do i=1, CS_numat
      CIF_atom_site(1,i) = dble(CS_atpos(1,i))
      CIF_atom_site(2,i) = dble(CS_atpos(2,i))
      CIF_atom_site(3,i) = dble(CS_atpos(3,i))
      CIF_atom_site(4,i) = dble(CS_atnum(i))
      CIF_atom_site(5,i) = dble(CS_atcrg(i))
      CIF_atom_site(6,i) = dble(CS_atocc(i))
      CIF_atom_site(7,i) = dble(2.0) ! adf type = biso
      CIF_atom_site(8,i) = dble(CS_atdwf(i)*100.0) ! biso value (A^2)
    end do
    !
  end if
  !
  ! - write
  !
  call CIF_WRITE(trim(sfile), nerr)
  if (nerr<=0) goto 13 ! check success
  !
  ! - free memory in the cifio module
  call CIF_CLEAR()
  !
  ! return without error
  !
  return
  
  !
  ! handle errors here, call specific error messages
  !
13 nerr = -1
  call CS_ERROR("Failed write CIF structure to ["//trim(sfile)//"].")
  return
19 nerr = -6
  call CS_ERROR("Failed to allocate cell data memory.")
  return
  
end subroutine CS_SAVE_CIFCELL
!**********************************************************************!






!**********************************************************************!
!
! CS_GETSLICE_PRO
!
! subroutine, calculates fresnel propagator for a slice n
! 090529: JB: supports now small sample tilts (tx,ty) specified in rad
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - slice data data to be allocated and set up
!   (e.g. call CS_SETSLICE_EQUIDIST)
!
! REMARKS: 
! - Fourier space is scrambled and transposed
!   thus, the propagator is an array of size (ny,nx)
!   where the supercell is an array of size (nx,ny)
! - Sample tilt is simulated by propagator shift
!   !! TODO check orientation of tilt
!
! INPUT:
!   integer*4 :: nslc               = slice number
!   integer*4 :: nx, ny             = discretisation of supercell axes
!   integer*4 :: nrx, nry           = repetition of supercell in slice
!   real*4 :: wl                    = electron wavelength [nm]
!   real*4 :: tx, ty                = sample tilt x and y [rad]
!
! IN/OUTPUT:
!   complex*8 :: pro(nx*nrx,ny*nry) = fresnel propagator
!   integer*4 :: nerr               = error code
!
subroutine CS_GETSLICE_PRO(nslc, nx, ny, nrx, nry, wl, tx, ty, pro, nerr)

  implicit none
  
  integer*4, intent(in) :: nslc, nx, ny, nrx, nry
  real*4, intent(in) :: wl, tx, ty
  integer*4, intent(inout) :: nerr
  complex*8, intent(inout) :: pro(ny*nry,nx*nrx)
  
  integer*4 :: i, j, i1, j1 ! iterators
  integer*4 :: dimx, dimy ! full size
  integer*4 :: nyqx, nyqy ! nyquist numbers
  real*4 :: sx, sy, itowx, itowy ! sampling constants
  real*4 :: dz ! slice thickness
  real*4 :: pco, pwthr ! propagator constants
  real*4 :: wx, wy, wx2, w2, wcx2, wc2 ! fourier space coordinates
  real*4 :: chi ! phase
  real*4 :: chit ! central beam phase
  complex*8 :: cval0, cval ! complex vars
  
  !
  ! preset error flag
  !
  nerr = 0
  
  !
  ! check mem alloc
  !
  if (.not.(CS_slicemem_allocated.and.CS_cellmem_allocated)) goto 13
  
  !
  ! check parameters
  !
  if (nx*nrx*nry*ny<=0) goto 14
  if (wl<=0.0) goto 14
  
  !
  ! get sampling constants
  !
  dimx = nx*nrx
  dimy = ny*nry
  nyqx = dimx/2
  nyqy = dimy/2
  sx = CS_scsx/real(nx)
  sy = CS_scsy/real(ny)
  itowx = wl / (real(nrx) * CS_scsx)
  itowy = wl / (real(nry) * CS_scsy)
  
  !
  ! set propagator constants
  !
  ! !!! this gets negative defocus, but neg. sign is compensated
  !     by leaving out the minus in the propagator calculation below (marked by ###)
  dz = CS_slczlim(1,nslc)-CS_slczlim(2,nslc)
  pco = CS_pi/wl*dz
  pwthr = min(itowx*real(nyqx),itowy*real(nyqy))*CS_PROSIZE_THR
  pwthr = pwthr*pwthr
  cval0 = cmplx(0.0,0.0)
  chit = pco*(tx*tx+ty*ty)
  
  !
  ! loop through fourierspace of slice and set propagator
  ! fourier space is transposed and scrambled !!! as with FFTs
  !
  do j=1, dimx ! x is with rows
    j1 = mod((j+nyqx-1),dimx)-nyqx
    wx = itowx*real(j1)
    wcx2 = wx*wx ! store angle wrt central beam
    wx = wx + tx ! apply tilt
    wx2 = wx*wx
    do i=1, dimy ! y is with columns
      i1 = mod((i+nyqy-1),dimy)-nyqy
      wy = itowy*real(i1)
      wc2 = wcx2 + wy*wy    ! store angle wrt central beam for aperture check
      wy = wy + ty          ! apply tilt
      w2 = wx2 + wy*wy
      if (wc2>pwthr) then ! hard aperture
        pro(i,j) = cval0
      else
        chi = pco*w2-chit ! correct for central beam phase
        cval = cmplx(cos(chi),sin(chi)) ! (###)
        pro(i,j) = cval
      end if
    end do
  end do
  
  return
  
  !
  ! error handling
  !
13 nerr=-1
  call CS_ERROR("Propagator memory not allocated.")
  return
14 nerr=-1
  call CS_ERROR("Invalid parameters for projector calculations.")
  return
  
end subroutine CS_GETSLICE_PRO
!**********************************************************************!


!**********************************************************************!
!
! CS_GETSLICE_PRO2
!
! subroutine, calculates propagator for a slice n
!             supports sample tilts (tx,ty) specified in rad
!
! This is the large angle version of (forwards) propagation
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - slice data data to be allocated and set up
!   (e.g. call CS_SETSLICE_EQUIDIST)
!
! REMARKS: 
! - Fourier space is scrambled and transposed
!   thus, the propagator is an array of size (ny,nx)
!   where the supercell is an array of size (nx,ny)
! - Sample tilt is simulated by propagator shift
!   !! TODO check orientation of tilt
!
! INPUT:
!   integer*4 :: nslc               = slice number
!   integer*4 :: nx, ny             = discretisation of supercell axes
!   integer*4 :: nrx, nry           = repetition of supercell in slice
!   real*4 :: wl                    = electron wavelength [nm]
!   real*4 :: tx, ty                = sample tilt x and y [rad]
!
! IN/OUTPUT:
!   complex*8 :: pro(nx*nrx,ny*nry) = fresnel propagator
!   integer*4 :: nerr               = error code
!
subroutine CS_GETSLICE_PRO2(nslc, nx, ny, nrx, nry, wl, tx, ty, pro, nerr)

  implicit none
  
  real*8, parameter :: dpi  = 3.141592653589D+000
  
  integer*4, intent(in) :: nslc, nx, ny, nrx, nry
  real*4, intent(in) :: wl, tx, ty
  integer*4, intent(inout) :: nerr
  complex*8, intent(inout) :: pro(ny*nry,nx*nrx)
  
  integer*4 :: i, j, i1, j1 ! iterators
  integer*4 :: dimx, dimy ! full size
  integer*4 :: nyqx, nyqy ! nyquist numbers
  real*8 :: dwl ! wavelength
  real*8 :: sx, sy, itogx, itogy ! sampling constants
  real*8 :: gmax, gthresh ! sampling and aperture ranges
  real*8 :: otx, oty, ot, cot, sot, od, dd ! object tilt parameters
  real*8 :: pfac ! propagation phase shift prefactor
  real*8 :: gx, gx2, gy, g2 ! reciprocal coordinates and their squares
  real*8 :: dt ! theta
  real*8 :: dz ! slice thickness
  real*8 :: chi ! phase
  complex*8 :: cval0, cval ! complex vars
  
  !
  ! preset error flag
  !
  nerr = 0
  
  !
  ! check mem alloc
  !
  if (.not.(CS_slicemem_allocated.and.CS_cellmem_allocated)) goto 13
  
  !
  ! check parameters
  !
  if (nx*nrx*nry*ny<=0) goto 14
  if (wl<=0.0) goto 14
  
  !
  ! get sampling constants
  !
  dwl = dble(wl)
  dimx = nx*nrx ! # samples x
  dimy = ny*nry ! # samples y
  nyqx = dimx/2 ! nyquist x
  nyqy = dimy/2 ! nyquist y
  sx = dble(CS_scsx)/dble(nx)  ! real-space sampling rate x
  sy = dble(CS_scsy)/dble(ny)  ! ... y
  itogx = 1.d+0 / (dble(nrx) * dble(CS_scsx)) ! reciprocal space sampling rate x
  itogy = 1.d+0 / (dble(nry) * dble(CS_scsy)) ! ... y
  
  !
  ! set propagator constants
  !
  ! !!! we work with the negative slice thickness here
  !     and leave out the minus in the phase shift later ... marked by (###)
  dz = dble(CS_slczlim(1,nslc)-CS_slczlim(2,nslc)) ! neg. slice thickness [nm]
  gmax = min(dble(itogx*nyqx),dble(itogy*nyqy)) ! worst gmax [1/nm]
  gthresh = dble(CS_PROSIZE_THR)*gmax ! aperture g
  gthresh = gthresh*gthresh ! aperture g ^2
  cval0 = cmplx(0.0,0.0)
  otx = dble(tx) ! internal object tilt x
  oty = dble(ty) ! internal object tilt y
  pfac = 2.0D+0 * dpi / dwl * dz ! phase shift pre-factor times slice thickness
  ot = dsqrt(otx*otx+oty*oty) ! object tilt magnitude [rad]
  cot = dcos(ot) ! object tilt cosine
  sot = dsin(ot) ! ... and sine
  od = 0.0
  if (ot>0.D+0) od = datan2(oty,otx) ! object tilt direction [rad]
  
  !
  ! loop through fourierspace of slice and set propagator
  ! fourier space is transposed and scrambled
  !
  do j=1, dimx ! L+ all lines ! x frequencies
    j1 = mod((j+nyqx-1),dimx)-nyqx
    gx = itogx*dble(j1)
    gx2 = gx * gx
    do i=1, dimy ! L+ all columns ! y frequencies
      i1 = mod((i+nyqy-1),dimy)-nyqy
      gy = itogy*real(i1)
      g2 = gx2 + gy*gy      ! g^2
      if (g2>gthresh) then ! hard aperture in wave frame (2/3)
        pro(i,j) = cval0 ! no transfer beyond aperture
      else
        dt = 2.D+0*dasin( 0.5D+0 * dwl * dsqrt(g2) ) ! |theta|
        dd = 0.D+0
        if (dt>0.D+0) dd = datan2(gy,gx) ! theta_phi
        ! calculate phase shift to non-diffracted beam
        ! actually (-chi) since dz was negative in pfac
        chi = pfac*(  &
              &   1.D+0/(cot*dcos(dt) + sot*dcos(od-dd)*dsin(dt) ) &
              & - 1.D+0/cot )
        ! calculate complex phase factor to wavefunction, exp(I*(-chi))
        cval = cmplx( dcos(chi), dsin(chi), 4) ! (###)
        ! store the phase factor
        pro(i,j) = cval
        !
      end if
    end do
  end do
  
  return
  
  !
  ! error handling
  !
13 nerr=-1
  call CS_ERROR("Propagator memory not allocated.")
  return
14 nerr=-1
  call CS_ERROR("Invalid parameters for projector calculations.")
  return
  
end subroutine CS_GETSLICE_PRO2
!**********************************************************************!


!**********************************************************************!
!
! CS_POLINT
!
! subroutine, calculates the n-1 order polynomial interpolation y
! at x of a table of values xa(n), ya(n) including error estimation.
!
! Adopted from the routine polint of Numerical Recipes §3.1
! implementing Neville's recursive algorithm.
!
! REQUIRES
! - input arrays of length n
! - no identical values in xa
!
! INPUT:
!   real*4 :: xa(n), ya(n)          = table of function values
!   integer*4 :: n                  = table size and n-1 = order of polyn.
!   real*4 :: x                     = requested interpolation position
!
! OUTPUT:
!   real*4 :: y                     = interpolated function value
!   real*4 :: dy                    = interpolation error estimate
!
subroutine CS_POLINT(xa,ya,n,x,y,dy)
  
  implicit none
  
  real*4, intent(in) :: xa(n), ya(n), x
  integer*4, intent(in) :: n
  real*4, intent(out) :: y, dy
  
  integer*4 :: i, m, ns
  integer*4 :: nalloc
  real*4 :: den, dif, dift, ho, hp, w
  real*4, allocatable :: c(:), d(:)
  
  ! init
  ns = 1
  y = ya(1)
  dy = 0.0
  nalloc = 0
  dif = abs(x-xa(1))
  
  ! allocation of helper arrays
  allocate(c(n), d(n), stat=nalloc)
  if (nalloc/=0) return ! error
  
  ! find the index ns of the closest table value
  do i=1,n
    dift = abs(x-xa(i))
	if (dift<dif) then
	  ns = i
	  dif = dift
	end if
	! .. and init tableaus of c and d with y values
	c(i) = ya(i)
	d(i) = ya(i)
  end do
  
  ! initial nearest neighbor approximation
  y = ya(ns)
  ns = ns-1
  
  ! For each column of the tableau, update c and d
  do m=1, n-1
    do i=1, n-m
	  ho = xa(i)-x
	  hp = xa(i+m)-x
	  w = c(i+1)-d(i)
	  den = ho - hp
	  ! if (den==0.) return ! Error when two xa's are equal.
	                        ! We will make sure that this doesn't
							! happen in the calling routines.
	  den = w/den
	  ! update c's and d's here
	  d(i) = hp*den
	  c(i) = ho*den
	end do
	!
	! After each column in the tableau is completed, we descide
	! which correction, c or d, we want to add to our accumulating
	! value of y, i.e. which path to take through the tableau -
	! forking up or down. We do this in such a way as to take the
    ! most "straight line" route through the tableau to its apex,
	! updating ns accordingly to keep track of where we are. This
	! route keeps the partial approximations centered (insofar as
	! possible) on the target x. The last dy added is thus the
    ! error indication.
	!
	if (2*ns < n-m) then
	  dy = c(ns+1)
	else
	  dy = d(ns)
	  ns = ns-1
	end if
	y = y + dy
  end do
  
  deallocate(c, d, stat=nalloc)
  return
  
end subroutine CS_POLINT
!**********************************************************************!



!**********************************************************************!
!
! CS_GETSCATTAMP
!
! subroutine, calculates the scattering amplitude for a given pixel
! and a given atomic index
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - scattering amplitude data to be allocated and prepared
!   (call CS_PREPARE_SCATTAMPS)
! - the interpolation scheme depends on CS_scaf_iptype
!
! INPUT:
!   real*4 :: g                     = scattering vector [1/nm]
!   integer*4 :: nat                = atom index in scattering table
!
! IN/OUTPUT:
!   complex*8 :: csf                = scattering amplitude
!
subroutine CS_GETSCATTAMP(csf,g,nat)

  implicit none
  
  integer*4, parameter :: nmax = 11
  integer*4, intent(in) :: nat
  real*4, intent(in) :: g
  complex*8, intent(inout) :: csf
  
  integer*4 :: ix                   ! closest table point (shifted by -1)
  integer*4 :: i0,i1,i2             ! target table range (shifted by -1) or neighbouring pixels
  integer*4 :: nt, nth              ! section size and half size
  integer*4 :: nipo                 ! interpolation order
  integer*4 :: i, j, j1             ! indices
  real*4 :: scx(nmax), sca1(nmax), sca2(nmax) ! table section for interpolation
  real*4 :: scf1, scf2              ! interpolated scattering amplitudes
  real*4 :: dsc1, dsc2              ! interpolation errors (not further used)
  real*4 :: pg                      ! floating point index
  real*4 :: pf                      ! fraction distance to left neighbour
  
  nipo = max(0,min(nmax-1,CS_scaf_iptype)) ! limited internal interpolation order (0 ... 10)
  !
  ix = nint(abs(g)/CS_scaitog)      ! get nearest neighbor
  csf = CS_scampdat(ix+1,nat)       ! init with nearest neighbor interpolation   <--- 0th order interpolation done
  !
  select case (nipo)                ! switch for interpolation order >0
    !
    case (1)                        !                                            <--- 1st order interpolation
      pg = g/CS_scaitog + 1.0
      i0 = int(pg)
      i1 = i0 + 1
      pf = pg - real(i0)
      csf = CS_scampdat(i0,nat)*(1.0-pf) + CS_scampdat(i1,nat)*pf
    !
    case (2:)                       !                                            <--- 2nd order interpolation
      scx = 0.0
      sca1 = 0.0
      sca2 = 0.0
      nt = nipo + 1                 ! we need one point more than the polynomial order
      nth = (nt-mod(nt,2))/2        ! half size of the section
      i1 = ix - nth                 ! section start index (may be negative)
      i2 = i1 + nt -1               ! section end index (may also be negative)
      do i=i1, i2                   ! loop through section of the data table
	    j1 = i - ix + nth + 1       ! index in the section array
        j = 1 + abs(i)              ! CS_scampdat(*,nat) is a radial array
	                                ! so mirror at the zero point
	    j = min(j, CS_scadim)       ! limit the access to the scattering data
	    scx(j1) = real(i)*CS_scaitog        ! interpolation g values
	    sca1(j1) = real(CS_scampdat(j,nat)) ! real part of the scattering amplitude
	    sca2(j1) = imag(CS_scampdat(j,nat)) ! imag part ...
      end do
      ! do the interpolation
      call CS_POLINT(scx(1:nt),sca1(1:nt),nt,g,scf1,dsc1) ! interpolation of real part
      call CS_POLINT(scx(1:nt),sca2(1:nt),nt,g,scf2,dsc2) ! interpolation of imaginary part
      csf = cmplx( scf1 , scf2 ) ! resynthesize to complex
  end select ! case (nipo)
  
  return
  
end subroutine CS_GETSCATTAMP
!**********************************************************************!




!**********************************************************************!
!
! CS_GETCELL_VOL
!
! subroutine, calculates the volume of the current super cell in nm^3
!
subroutine CS_GETCELL_VOL(vol)
  implicit none
  real*4, intent(out) :: vol
  real*4 :: ca, cb, cc, da
  vol = 0.0
  ca = cos(CS_scalpha*CS_rd2r)
  cb = cos(CS_scbeta *CS_rd2r)
  cc = cos(CS_scgamma*CS_rd2r)
  da = sqrt(abs(1.0-ca*ca-cb*cb-cc*cc+2.0*ca*cb*cc))
  vol = CS_scsx * CS_scsy * CS_scsz * da
end subroutine CS_GETCELL_VOL
!**********************************************************************!




!**********************************************************************!
!
! CS_GETCELL_POT
!
! subroutine, calculates potential of the current supercell in 3D
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - scattering amplitude data to be allocated and prepared
!   (call CS_PREPARE_SCATTAMPS)
!
! REMARKS
! - when using frozen lattice displacements, be aware that
!   damping of scattering factors by debye-waller factors
!   is senseless and wrong! So turn them off when calculating
!   the scatterung factors by CS_PREPARE_SCATTAMPS
!
! - when using frozen lattice generation, call InitRand() before
!
! - Adds the ionic charge potentials here instead of to the
!   screened potentials in the CS_PREPARE_SCATTAMPS routine.
!   This may help to avoid some numerical round-off problems.
!
! - In order to enable parallel computing calls, we allocate
!   a local complex*8 square array for the fourier transforms.
!
! INPUT:
!   integer*4 :: nx, ny, nz         = discretisation of supercell axes
!   integer*4 :: nfl                = create frozen latice displacements
!   integer*4 :: ndw                = apply Debye-Waller factors
!   real*4 :: wl                    = electron wavelength [nm]
!
! IN/OUTPUT:
!   complex*8 :: pot(nx,ny,nz)      = 3D scattering potential
!   integer*4 :: nerr               = error code
!
subroutine CS_GETCELL_POT(nx, ny, nz, nfl, ndw, wl, pot, nerr)

  implicit none
  
  integer*4, intent(in) :: nx, ny, nz, nfl, ndw
  real*4, intent(in) :: wl
  integer*4, intent(inout) :: nerr
  complex*8, intent(out) :: pot(nx,ny,nz)
  
  integer*4 :: nft1, nft2, ncalc, ncalcmax
  integer*4 :: i, j, k, i1, j1, k1, ia, jptr ! iterators
  integer*4 :: npc, lndw
  integer*4 :: na ! atom count
  integer*4 :: nyqx, nyqy, nyqz ! sampling numbers
  integer*4 :: infl ! internal frozen lattice flag
  real*4 :: sx, sy, sz, itogx, itogy, itogz ! sampling constants
  real*4 :: gx, gy, gz, gx2, gz2, g, g2, gxy, gxy2, gmax, gmaxz ! fourier space coordinates
  real*4 :: ht ! high tension in kV
  real*8 :: trpha ! translation phase
  real*4 :: x, y, z, dx, dy, dz ! atom position
  real*4 :: dwc, fldamp, dwf ! displacement calculation helpers
  real*4 :: pocc ! atom occupancy
  real*4 :: vol ! volume in nm^3
  real*4 :: pscal ! potential scaling factor
  real*4 :: apxy, apz ! aperture function, damping factors at outer Fourier-space perimeters
  real*4 :: apthr ! potential aperture cut-off threshold
  real*8 :: crgio ! ionic charge
  real*8 :: relcor ! relativistic correction factor
  real*8 :: pfacio ! ionic contribution prefactor
  real*8 :: rs2io ! ionic contribution 1/s^2 term
  complex*8 :: cval0, cval ! some complex vars
  complex*8 :: csf
  complex*16 :: ctr
  real*4, dimension(nx) :: agx
  real*4, dimension(ny) :: agy
  real*4, dimension(nz) :: agz
  real*4, dimension(:,:), allocatable :: fld ! frozen lattice displacements and other atom propierties
  integer*4, dimension(:,:), allocatable :: ild ! index table for atoms in the slice
  integer*4, dimension(:), allocatable :: uatl ! list of used atom types in the slice
  complex*16 :: cstrfe, cstrfi
  logical :: dwflg
  complex*8, allocatable :: cpotft(:,:,:) ! 3D potential working field
  complex*8, allocatable :: ctmp2d(:,:), ctmp1d(:) ! transformation helper arrays

  ! external routines for 1D Fourier transform (link FFTs.f)
  external :: ODDCC128, ODDCC256,  ODDCC512,  ODDCC1024
  external :: ODDCC2048, ODDCC4096, ODDCC8192
  ! external routines for 2D Fourier transform (link FFTs.f)
  external :: ODDCC128S, ODDCC256S,  ODDCC512S,  ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S
  ! other external routines
  real*4, external :: UniRand, GaussRand, getdwf
  
  !
  ! --- Initializations
  !
  nerr = 0
  nft1 = 2**CEILING( LOG( real( nz ) )/LOG(2.0) ) ! next 2^N above nz
  nft1 = max(nft1, NFT_MIN)
  if (nft1>NFT_MAX) goto 14
  nft2 = 2**CEILING( LOG( real( max(nx,ny) ) )/LOG(2.0) ) ! next 2^N above max(nx,ny)
  nft2 = max(nft2, NFT_MIN)
  if (nft2>NFT_MAX) goto 14
  apxy = 1.0 ! init aperture xy
  apz = 1.0  ! init aperture z
  
  !
  ! --- Checks on preferences
  !
  if (.not.(CS_cellmem_allocated.and. &
     &      CS_scattamp_prepared)) goto 13
  
  !
  ! --- Check input parameters
  !
  if (nx<=0.or.ny<=0.or.nz<=0) goto 14
  if (wl<=0.0) goto 14
  
  !
  ! handle frozen lattice flag
  !
  infl = nfl
  if (infl<0) infl = 0
  if (infl>1) infl = 1
  fldamp = 0.0
  dwc = 0.0
  lndw = ndw
  if (infl==1) then
    lndw = 0
    fldamp = CS_rr8p2
  end if
  
  !
  ! handle DWF flag (for ionic potential calculations only)
  !
  if (lndw==0) then
    dwflg = .FALSE.
  else
    dwflg = .TRUE.
  end if
  
  !
  ! get sampling constants
  !
  nyqx = nx/2
  nyqy = ny/2
  nyqz = nz/2
  sx = CS_scsx/real(nx)
  sy = CS_scsy/real(ny)
  sz = CS_scsz/real(nz)
  itogx = 1.0 / CS_scsx
  itogy = 1.0 / CS_scsy
  itogz = 1.0 / CS_scsz
  cval0 = cmplx(0.0,0.0)
  cval = cval0
  
  !
  ! save as global
  !
  CS_sdimx = nx
  CS_sdimy = ny
  CS_sdimz = nz
  CS_repx = 1
  CS_repy = 1
  CS_repz = 1
  CS_sampx = sx
  CS_sampy = sy
  CS_sampz = sz
  
  !
  ! other parameters
  !
  ht = CS_WL2HT(wl)                                 ! high tension in kV
  na = CS_numat                                     ! number of atoms in cell
  vol = CS_scsx*CS_scsy*CS_scsz                     ! supercell volume in nm^3 (assuming orthogonal)
  pscal = CS_v0 / vol                               ! scattering potential pre-factor
  relcor = dble((CS_elm0 + ht)/CS_elm0)             ! relativistic correction, gamma
  pfacio = relcor * dble(CS_scaprea*CS_fpi*10.)     ! pre-factor for ionic part: rel-corr * C * 4 Pi * 10. -> 1/nm
  
  !
  ! clear fourier space working field and prepare helper arrays
  !
!  cw = cval0
  gx2 = 0.0
  g2 = 0.0
  do j=1, nx
    j1 = mod((j+nyqx-1),nx)-nyqx
    agx(j) = itogx*real(j1)
  end do
  do i=1, ny
    i1 = mod((i+nyqy-1),ny)-nyqy
    agy(i) = itogy*real(i1)
  end do
  do k=1, nz
    k1 = mod((k+nyqz-1),nz)-nyqz
    agz(k) = itogz*real(k1)
  end do
  gmax = min(itogx*real(nyqx),itogy*real(nyqy))      ! smallest max. diffraction in x and y
  gmaxz = itogz*real(nyqz)                           ! max. diffraction in z
  apthr = min(itogx*real(nyqx-1),itogy*real(nyqy-1)) ! hard aperture cut-off for saving calculation time
  
  !
  ! clear potential
  !
  pot = cval0
  
  !
  ! return in case of empty cell
  !
  if (na<=0) then 
    return
  end if
  
  !
  ! precalculate atom displacements for frozen lattice
  !
  allocate(fld(6,na),ild(2,na),uatl(CS_scampnum+1),stat=nerr)
  if (nerr/=0) goto 15
  fld = 0.0
  ild = 0
  uatl = 0
  do ia=1, na ! loop ia over all atoms in cell
    ild(1,ia) = ia
    ild(2,ia) = CS_scampptr(ia)
    uatl(ild(2,ia)) = 1
    fld(1,ia) = CS_atpos(1,ia)*CS_scsx              ! get atom x-position
    fld(2,ia) = CS_atpos(2,ia)*CS_scsy              ! get atom y-position
    fld(3,ia) = CS_atpos(3,ia)*CS_scsz              ! get atom z-position
    dwc = sqrt(CS_atdwf(ia))                ! get debye-waller parameter from (nm**2) to (nm)
    if (infl==1) then                       ! dice frozen lattice displacements (x,y) and apply
      if (CS_atlnk(ia)==0) then ! independent site
        dx = CS_rr8p2*dwc*GaussRand()
        fld(1,ia) = fld(1,ia) + dx            ! add random x displacement
        dy = CS_rr8p2*dwc*GaussRand()
        fld(2,ia) = fld(2,ia) + dy            ! add random y displacement
        dz = CS_rr8p2*dwc*GaussRand()
        fld(3,ia) = fld(3,ia) + dz            ! add random z displacement
      else ! dependent site (copy position from linked site)
        fld(1:3,ia) = fld(1:3,CS_atlnk(ia))
      end if
    end if
    fld(4,ia) = dwc                         ! Biso
    fld(5,ia) = CS_atocc(ia)                ! Occupancy
    fld(6,ia) = CS_atcrg(ia)                ! ionic charge
  end do ! loop ia over all atoms in slice
  if (CS_useextsca==0) then
    fld(6,:) = 0.0 ! reset ionic charge to zero for standard potentials, they are not for ions
    !call CS_MESSAGE("- setting atomic charges to zero. Using neutral atom scattering tables.")
  end if

  
  !
  ! calculate super-cell potential in fourier space
  ! !!! fourierspace is scrambled and transposed
  !
  ! - allocate the fourier-space potential memory
  allocate(cpotft(ny,nx,nz), stat=nerr)
  if (nerr/=0) goto 15
  cpotft = cmplx(0.0,0.0)
  !
  ncalc = 0
  ncalcmax = nx*ny*nz
  if (CS_doconsolemsg>0) then
    call CS_PROG_START(ncalcmax,1.0)
  end if
  do k=1, nz ! loop planes (wz-axis)
    gz = agz(k)
    gz2 = gz*gz
    !
    ! --> Precalculate an aperture that smoothly limits the diffraction angles z
    !     This funcion is sigmoid and drops from 1 to 0 between 0.8*Nyq and 1*Nyq
    !     independent on sampling
    if (CS_do_bwl_pot) apz = 0.5 - 0.5*tanh( (sqrt(gz2)/gmaxz-0.9)*30.0 )
    !
    do j=1, nx ! loop rows (wx-axis, transposed)
      gx = agx(j)
      gx2 = gx*gx
      do i=1, ny ! loop columns (wy-axis, transposed)
        gy = agy(i)
        gxy2 = gx2 + gy*gy
        g2 = gz2 + gxy2
        g = sqrt(g2)
        gxy = sqrt(gxy2)
        !
!        ! hard cut-off to save calculation time (round in gx-gy, no limit in gz)
!        ! we have to make this round cut-off in order to obtain the correct
!        ! x-y relation of the resulting structure factors 
!        if (g>apthr) cycle ! scattering amplitude aperture
        !
        ! ionic part 1/s^2 term
        rs2io = dble( 1.0 / ( 0.25*g2 + CS_iyr*CS_iyr) ) ! 1/(s^2+a^2) [nm^-2]
        !
        ! --> Precalculate an aperture that smoothly limits the diffraction angles (x,y)
        !     This funcion is sigmoid and drops from 1 to 0 between 0.8*Nyq and 1*Nyq
        !     independent on sampling
        if (CS_do_bwl_pot) apxy = 0.5 - 0.5*tanh( (gxy/gmax-0.9)*30.0 )
        !
        ! --> Precalculate atom type dependent values which would be too
        !     Costly for repeating it for each individual atom
        !
        if (CS_scampnum>0) then
          do jptr=1, CS_scampnum ! loop ia over all atom types in the slice
            if (uatl(jptr)==0) cycle ! this atom type is not used in this slice
            dwc = CS_scampdwf(jptr)
            ! calculcation of scattering amplitude for the current atom type and g
            call CS_GETSCATTAMP(csf,g,jptr)
            ! get the DWF for the current atom type and the current g
            dwf = getdwf( g, dwc, dwflg )
            ! store as temporary data
            CS_scamptmp(1,jptr) = real(csf)
            CS_scamptmp(2,jptr) = imag(csf)
            CS_scamptmp(3,jptr) = dwf
          end do
        end if
        !
        ! --> Calculate the structure factor for the current diffraction vector (gx,gy,gz)
        !
        cstrfe = dcmplx( 0.0, 0.0 ) ! reset the screened potential structure factor
        cstrfi = dcmplx( 0.0, 0.0 ) ! reset the ionic potential structure factor
        do ia=1, na ! loop ia over all atoms in the cell

        ! ia = atom index in super cell (in this routine, not as in the projected potential routine)
          jptr = ild(2,ia)                        ! get atom type index in scattering data
          x = fld(1,ia)                           ! get atom x-position
          y = fld(2,ia)                           ! get atom y-position
          z = fld(3,ia)                           ! get atom z-position
          pocc = fld(5,ia)                        ! get atom occupancy
          crgio = dble( fld(6,ia) )               ! get ionic charge of the atom
        
          csf = cmplx(CS_scamptmp(1,jptr),CS_scamptmp(2,jptr))*pocc ! rescale scattering amplitude to occupancy
          ! calculation of the phase factor for the current atom and g
          trpha = dble(CS_tpi*(x*gx+y*gy+z*gz))   ! get translation phase = 2*pi*g*r
          ctr = dcmplx(dcos(trpha),-dsin(trpha))  ! get translation phase factor = exp(-2*pi*I*g*r)
          cstrfe = cstrfe + dcmplx(csf)*ctr       ! sum up to the screened potential structure factor
          if (crgio/=0.0) then                    ! sum up to the ionic structure factor
            dwf = 1.0
            if (infl==0) dwf = CS_scamptmp(3,jptr) ! get the DWF for the ionic contribution (only for non-FL calculations)
            cstrfi = cstrfi + pocc * crgio * dble(dwf) * ctr
          end if
        
        end do ! loop ia over all atoms in the cell
      
        cstrfi = cstrfi * rs2io * pfacio          ! multiply the ionic structure factor by 1/s^2 (only for s>0)
                                                  ! ... and with the scattering prefactor which includes
                                                  !     a factor of 4 Pi, a relativistic correction and other
                                                  !     constants of the Coulmob interaction.
      
        cpotft(i,j,k) = apz*apxy*(cstrfe + cstrfi)! store the combined structure factors of the screened
                                                  ! potential and the ionic charge potential
      
        ncalc = ncalc + 1 ! increase calculated pixel count
        
        
  
      end do ! loop columns (wy-axis, transposed)
    end do ! loop rows (wx-axis, transposed)
    if (CS_doconsolemsg>0) then
      ! update progress indicator
      call CS_PROG_UPDATE(ncalc)
    end if
  end do ! loop planes (wz-axis)
  
  if (CS_doconsolemsg>0) then
    call CS_PROG_STOP(ncalcmax)
  end if
  deallocate(fld,ild,uatl,stat=nerr)
  if (nerr/=0) goto 16
  
!  !
!  ! debug table of scattering factors
!  !
!  ja = 148                                ! get atom index in super cell
!  dwc = CS_atdwf(ja)                      ! get debye-waller parameter B (nm**2)
!  jptr = CS_scampptr(ja)                  ! get atom type index in scattering data
!  crgio = dble( CS_atcrg(ja) )            ! get ionic charge of the atom
!  write(*,*) "s, f[el]"
!  do i=1, 35
!    g = real(i)
!    call CS_GETSCATTAMP(csf,g,jptr)         ! get scattering amplitude from the screened potential
!    dwf = dwfjbr( g, dwc, dwflg )
!    rs2io = dble( 20.0 / g )**2
!    write(*,'(F5.2,", ",F8.3)') g*0.05, (real(csf) + crgio * dble(dwf) * rs2io * pfacio ) / relcor / 0.1 / d2pi / 2 / dble(dwf)
!    
!  end do
!  !
  
  !
  ! transform back to real space
  !
  write(unit=6,fmt='(A)') "  3D Fourier Transform ... " 
  ! - prepare transform helper arrays
  allocate(ctmp2d(nft2,nft2),ctmp1d(nft1),stat=nerr)
  if (nerr/=0) goto 15

  !
  ! - transform all gz-lines to z-lines
  ctmp1d = cmplx(0.0,0.0)
  do j=1, nx
    do i=1, ny
      do k=1, nz
        ctmp1d(k) = cpotft(i,j,k)
      end do
      select case (nft1)
      case (128)
        call ODDCC128 (ctmp1d,nz,'BACK')
      case (256)
        call ODDCC256 (ctmp1d,nz,'BACK')
      case (512)
        call ODDCC512 (ctmp1d,nz,'BACK')
      case (1024)
        call ODDCC1024(ctmp1d,nz,'BACK')
      case (2048)
        call ODDCC2048(ctmp1d,nz,'BACK')
      case (4096)
        call ODDCC4096(ctmp1d,nz,'BACK')
      case (8192)
        call ODDCC8192(ctmp1d,nz,'BACK')
      end select ! case (nft1)
      do k=1, nz
        cpotft(i,j,k) = ctmp1d(k)
      end do
    end do
  end do
  !
  ! - transform all gx-gy-planes to x-y-planes
  ctmp2d = cmplx(0.0,0.0)
  do k=1, nz
    do j=1, nx
      ctmp2d(1:ny,j) = cpotft(1:ny,j,k)
    end do
    select case (nft2)
    case (128)
      call ODDCC128S (ctmp2d,nx,ny,'BACK')
    case (256)
      call ODDCC256S (ctmp2d,nx,ny,'BACK')
    case (512)
      call ODDCC512S (ctmp2d,nx,ny,'BACK')
    case (1024)
      call ODDCC1024S(ctmp2d,nx,ny,'BACK')
    case (2048)
      call ODDCC2048S(ctmp2d,nx,ny,'BACK')
    case (4096)
      call ODDCC4096S(ctmp2d,nx,ny,'BACK')
    case (8192)
      call ODDCC8192S(ctmp2d,nx,ny,'BACK')
    end select ! case (nft2)
    do j=1, ny
      pot(1:nx,j,k) = ctmp2d(1:nx,j)*pscal
    end do
  end do
  
 
  ! - get rid of the helper arrays
  deallocate(ctmp2d,ctmp1d,stat=nerr)
  deallocate(cpotft,stat=nerr)
  return
  
  !
  ! error handling
  !
13 nerr=-1
  call CS_ERROR("Potential memory not allocated.")
  return
14 nerr=-1
  call CS_ERROR("Invalid parameters for 3D potential calculations.")
  return
15 nerr=-1
  call CS_ERROR("Memory allocation failed.")
  return
16 nerr=-1
  call CS_ERROR("Memory deallocation failed.")
  return
  
end subroutine CS_GETCELL_POT
!**********************************************************************!




!**********************************************************************!
!
! CS_GETSLICE_POT
!
! subroutine, calculates projected potential of slice n
!             applies thermal displacements in case of frozen-lattice
!             simulation -> generates 1 FL variant
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - slice data data to be allocated and set up
!   (e.g. call CS_SETSLICE_EQUIDIST)
! - scattering amplitude data to be allocated and prepared
!   (call CS_PREPARE_SCATTAMPS)
!
! REMARKS
! - when using frozen lattice displacements, be aware that
!   damping of scattering factors by debye-waller factors
!   is senseless and wrong! So turn them off when calculating
!   the scatterung factors by CS_PREPARE_SCATTAMPS
!
! - when using frozen lattice generation, call InitRand() before
!
! - Adds the ionic charge potentials here instead of to the
!   screened potentials in the CS_PREPARE_SCATTAMPS routine.
!   This may help to avoid some numerical round-off problems.
!
! - In order to enable parallel computing calls, we allocate
!   a local complex*8 square array for the fourier transforms.
!   The size of the transformation array is adjusted to the
!   needs of the call, keeping the memory consumption as low
!   as possible.
!
! INPUT:
!   integer*4 :: nslc               = slice number
!   integer*4 :: nx, ny             = discretisation of supercell axes
!   integer*4 :: nrx, nry           = repetition of supercell in slice
!   integer*4 :: nfl                = create frozen latice displacements
!   integer*4 :: ndw                = apply Debye-Waller factors
!   real*4 :: wl                    = electron wavelength [nm]
!
! IN/OUTPUT:
!   complex*8 :: pot(nx*nrx,ny*nry) = projected potential of slice
!   integer*4 :: nerr               = error code
!
subroutine CS_GETSLICE_POT(nslc, nx, ny, nrx, nry, nfl, ndw, wl, pot, nerr)

  implicit none
  
  integer*4, intent(in) :: nslc, nx, ny, nrx, nry, nfl, ndw
  real*4, intent(in) :: wl
  integer*4, intent(inout) :: nerr
  complex*8, intent(inout) :: pot(nx*nrx,ny*nry)
  
  integer*4 :: i, j, i1, j1, ia, ja, jptr ! iterators
  integer*4 :: npc, lndw
  integer*4 :: na ! atom count
  integer*4 :: nft, dimx, dimy, nyqx, nyqy ! sampling numbers
  integer*4 :: infl ! internal frozen lattice flag
  real*4 :: sx, sy, itogx, itogy ! sampling constants
  real*4 :: gx, gy, gx2, g, g2, gmax ! fourier space coordinates
  real*4 :: ht ! high tension in kV
  real*8 :: trpha ! translation phase
  real*4 :: x, y, dx, dy ! atom position
  real*4 :: dwc, fldamp, dwf ! displacement calculation helpers
  real*4 :: pocc ! atom occupancy
  real*4 :: dsz, vol ! slice thickness and volume in nm^3
  real*4 :: pscal ! potential scaling factor
  real*4 :: apdmp ! aperture function (damping factor at outer Fourier-space perimeter >2/3 gmax)
  real*4 :: apthr ! potential aperture cut-off threshold
  real*8 :: crgio ! ionic charge
  real*8 :: relcor ! relativistic correction factor
  real*8 :: pfacio ! ionic contribution prefactor
  real*8 :: rs2io ! ionic contribution 1/s^2 term
  complex*8 :: cval0, cval ! some complex vars
  complex*8 :: csf
  complex*16 :: ctr, cstrfe, cstrfi
  real*4, dimension(:,:), allocatable :: scamptmp ! temporary scattering amplitudes used here
  real*4, dimension(:), allocatable :: agx, agy ! prepared field of spatial frequencies
  real*4, dimension(:,:), allocatable :: fld ! frozen lattice displacements and other atom properties
  integer*4, dimension(:,:), allocatable :: ild ! index table for atoms in the slice
  integer*4, dimension(:), allocatable :: jld ! index table linking from atom list to the slice idx
  integer*4, dimension(:), allocatable :: uatl ! list of used atom types in the slice
  complex*8, dimension(:,:), allocatable :: lcw ! local working array for FFTs
  
  logical :: dwflg

  ! external routines for 2D Fourier transform (link FFTs.f)
  external :: ODDCC128S, ODDCC256S,  ODDCC512S,  ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S
  real*4, external :: UniRand, GaussRand, getdwf
  
  !
  ! preset error flag
  !
  nerr = 0
  apdmp = 1.0 ! init without dampening
  
  !
  ! check mem alloc
  !
  if (.not.(CS_slicemem_allocated.and. &
     &      CS_cellmem_allocated.and. &
     &      CS_scattamp_prepared)) goto 13

  
  !
  ! check parameters
  !
  if (nx<=0.or.nrx<=0.or.nry<=0.or.ny<=0) goto 14
  if (wl<=0.0) goto 14
  if (nslc<1.or.nslc>CS_nspsc) goto 14
  
  !
  ! handle frozen lattice flag
  !
  infl = nfl
  if (infl<0) infl = 0
  if (infl>1) infl = 1
  fldamp = 0.0
  dwc = 0.0
  lndw = ndw
  if (infl==1) then
    lndw = 0
    fldamp = CS_rr8p2
  end if
  !
  ! handle DWF flag (for ionic potential calculations only)
  !
  if (lndw==0) then
    dwflg = .FALSE.
  else
    dwflg = .TRUE.
  end if
  
  
  !
  ! get sampling constants
  !
  dimx = nx*nrx ! number of potential output samples (x)
  dimy = ny*nry ! number of potential output samples (y)
  nyqx = nx/2 ! Nyquist number (x) for potential generation
  nyqy = ny/2 ! Nyquist number (y) for potential generation
  sx = CS_scsx/real(nx) ! real-space sampling rate (x) for potential generation
  sy = CS_scsy/real(ny) ! real-space sampling rate (y) ...
  itogx = 1.0 / CS_scsx ! reciprocal space sampling rate (x) for potential generation
  itogy = 1.0 / CS_scsy ! reciprocal space sampling rate (y) for potential generation
  cval0 = cmplx(0.0,0.0)
  cval = cval0
  ! set size and allocate FFT working array
  nft = 2**CEILING( LOG( real( max(nx,ny) ) )/LOG(2.0) ) ! next 2^N above max(nx,ny)
  nft = max(nft, NFT_MIN)
  if (nft>NFT_MAX) goto 14
  allocate(lcw(nft,nft),stat=nerr)
  if (nerr/=0) goto 15
  lcw(1:nft,1:nft) = cmplx(0.0,0.0) ! this is the working array (square 2^N size)
  
  !
  ! save as global
  !
  CS_sdimx = dimx
  CS_sdimy = dimy
  CS_repx = nrx
  CS_repy = nry
  CS_sampx = sx
  CS_sampy = sy
  
  !
  ! other parameters
  !
  ht = CS_WL2HT(wl)                                 ! high tension in kV
  na = CS_slcatnum(nslc)                            ! number of atoms in slice
  dsz = CS_slczlim(2,nslc)-CS_slczlim(1,nslc)       ! thickness of the slice in nm
  vol = CS_scsx*CS_scsy*dsz                         ! supercell volume in nm^3
  pscal = CS_v0 / vol                               ! potential scaling factor Vamp0 / Volume
  relcor = dble((CS_elm0 + ht)/CS_elm0)             ! relativistic correction for high energy electrons
  pfacio = relcor * dble(CS_scaprea*CS_fpi*10.)     ! prefector for the ionic rest charge potential: rel-corr * C * 4 Pi * 10. [-> nm^-1]
  
  !
  ! clear fourier space working field and prepare helper arrays
  !
!  cw = cval0
  gx2 = 0.0
  g2 = 0.0
  gmax = 0.0
  ! agx and agy will hold reciprocal space coordinates
  allocate(agx(nx),agy(ny),stat=nerr)
  if (nerr/=0) goto 15
  do j=1, nx
    j1 = mod((j+nyqx-1),dimx)-nyqx
    agx(j) = itogx*real(j1)
    ! gx2 = max(gx2,agx(j))
  end do
  do i=1, ny
    i1 = mod((i+nyqy-1),dimy)-nyqy
    agy(i) = itogy*real(i1)
    ! gmax = max(gmax,agy(i))
  end do
  gmax = min(itogx*real(nyqx),itogy*real(nyqy)) ! store gmax
  apthr = min(itogx*real(nyqx-1),itogy*real(nyqy-1)) ! set g threshold
  
  !
  ! clear output potential array
  !
  pot = cval0
  
  !
  ! return in case of empty slice
  !
  if (na<=0) then 
    goto 10
  end if
  
  !
  ! setup atom positions and indices for this slice +
  ! pre-calculate atomic displacements for frozen lattice (infl==1)
  !
  allocate( fld(6,na),ild(2,na),jld(CS_numat),uatl(CS_scampnum+1), &
          & scamptmp(3,CS_scampnum),stat=nerr)
  if (nerr/=0) goto 15
  fld = 0.0 ! this keeps atom site data (pos + displacement, Biso, occ, ionic charge)
  ild = 0 ! this is a hash table to re-access atomic data and form factors
  jld = 0 ! this is a hash table that links the slice atom id with the cell atom id
  uatl = 0 ! usage flag of the atomic data
  scamptmp = 0.0 ! temporary copies of form factors
  do ia=1, na ! loop ia over all atoms in slice
    ja = CS_slcatacc(ia,nslc)               ! get atom index in super cell
    jld(ja) = ia                            ! store atom slice ID
    ild(1,ia) = ja                          ! store atom cell ID
    ild(2,ia) = CS_scampptr(ja)             ! store form factor ID
    uatl(ild(2,ia)) = 1                     ! flag the use of the form factor
    fld(1,ia) = CS_atpos(1,ja)*CS_scsx              ! get atom avg. x-position
    fld(2,ia) = CS_atpos(2,ja)*CS_scsy              ! get atom avg. y-position
 !   fld(3,ia) = CS_atpos(3,ja)*CS_scsz - CS_slczlim(1,nslc) ! get atom z-position in slice
    dwc = sqrt(CS_atdwf(ja))                ! get debye-waller parameter from (nm**2) to (nm)
    if (infl==1) then                       ! dice frozen lattice displacements (x,y) and apply
      if (CS_atlnk(ja)==0) then ! independent site
        dx = CS_rr8p2*dwc*GaussRand()         !   CS_rr8p2*dwc = sqrt( Biso / 8*pi^2) = sqrt( <u_s^2> )
        fld(1,ia) = fld(1,ia) + dx            ! add random x displacement
        dy = CS_rr8p2*dwc*GaussRand()
        fld(2,ia) = fld(2,ia) + dy            ! add random y displacement
!        dz = CS_rr8p2*dwc*GaussRand()
!        fld(3,ia) = dz !+ fld(3,ia)            ! random z displacments ! are ignored -> projected
      else ! dependent site (copy positions from linked site)
        fld(1:3,ia) = fld(1:3,jld(CS_atlnk(ja)))
      end if
    end if
    fld(4,ia) = dwc                         ! sqrt(Biso)
    fld(5,ia) = CS_atocc(ja)                ! Occupancy
    fld(6,ia) = CS_atcrg(ja)                ! ionic charge
  end do ! loop ia over all atoms in slice
  if (CS_useextsca==0) then
    fld(6,:) = 0.0 ! reset ionic charge to zero for standard potentials, they are not for ions
    !call CS_MESSAGE("- setting atomic charges to zero. Using neutral atom scattering tables.")
  end if
  
  !
  ! calculate super-cell potential in fourier space
  ! !!! fourierspace is scrambled and transposed
  !
  !if (CS_doconsolemsg>0) then
  !  npc = 0
  !  write(unit=6,fmt='(A,$)') "  > Progress: "
  !end if
  !
  do j=1, nx ! loop rows (wx-axis, transposed)
    gx = agx(j) ! spatial frequency
    gx2 = gx*gx
    do i=1, ny ! loop columns (wy-axis, transposed)
      gy = agy(i) ! spatial frequency
      g2 = gx2 + gy*gy
      g = sqrt(g2)
      !
      ! ionic part 1/s^2 term
      rs2io = dble( 1.0 / ( 0.25*g2 + CS_iyr*CS_iyr) ) ! 1/(s^2+a^2) [nm^-2]
      !
      ! --> Precalculate an aperture that smoothly limits the diffraction angles (x,y)
      !     This funcion is sigmoid and drops from 1 to 0 between 0.8*Nyq and 1*Nyq
      !     independent on sampling
      if (CS_do_bwl_pot) apdmp = 0.5 - 0.5*tanh( (g/gmax-0.9)*30.0 )
      !
      ! --> Precalculate atom type dependent values which would be too
      !     Costly for repeating it for each individual atom, however,
      !     they also depend on g
      if (CS_scampnum>0) then
        do jptr=1, CS_scampnum ! loop ia over all atom types in the slice
          if (uatl(jptr)==0) cycle ! this atom type is not used in this slice
          dwc = CS_scampdwf(jptr) ! get the DW-parameter
          ! calculcation of scattering amplitude for the current atom type and g
          call CS_GETSCATTAMP(csf,g,jptr)
          ! get the DWF for the current atom type and the current g
          dwf = getdwf( g, dwc, dwflg )
          ! store as temporary data
          scamptmp(1,jptr) = real(csf)
          scamptmp(2,jptr) = imag(csf)
          scamptmp(3,jptr) = dwf
        end do
      end if
      !
      ! --> Calculate the structure factor for the current diffraction vector (gx,gy)
      !
      cstrfe = dcmplx( 0.0, 0.0 ) ! reset the screened potential structure factor
      cstrfi = dcmplx( 0.0, 0.0 ) ! reset the ionic potential structure factor
      do ia=1, na ! loop ia over all atoms in slice
        !
        ja = ild(1,ia)                          ! get atom index in super cell
        jptr = ild(2,ia)                        ! get atom type index in scattering data
        x = fld(1,ia)                           ! get atom x-position
        y = fld(2,ia)                           ! get atom y-position
!        z = fld(3,ia)                           ! get atom z-position
!        dwc = fld(4,ia)                         ! get debye-waller parameter B (nm**2)
        pocc = fld(5,ia)                        ! get atom occupancy
        crgio = dble( fld(6,ia) )               ! get ionic charge of the atom
        !
        csf = cmplx(scamptmp(1,jptr),scamptmp(2,jptr))*pocc ! rescale scattering amplitude to occupancy
        !
        ! calculation of the phase factor for the current atom and g
        trpha = dble(CS_tpi*(x*gx+y*gy)) !+z*wl*g2*0.5)            ! get translation phase = 2*pi*g*r !+ 2*pi*Z*g^2*lambda/2
        ctr = dcmplx(dcos(trpha),-dsin(trpha))  ! get translation phase factor = exp(-2*pi*I*g*r)
        cstrfe = cstrfe + dcmplx(csf)*ctr       ! sum up to the screened potential structure factor
        if (crgio/=0.0) then                    ! sum up to the ionic structure factor
          dwf = 1.0
          if (infl==0) dwf = scamptmp(3,jptr)   ! get the DWF for the ionic contribution (only for non-fl calculatins)
          cstrfi = cstrfi + pocc*crgio*dble(dwf)*ctr ! add to the ionic structure factor occ*charge*DWF*(translation term)
        end if
        !
      end do ! loop ia over all atoms in slice
      !
      cstrfi = cstrfi * rs2io * pfacio          ! multiply the ionic structure factor by 1/s^2 (only for s>0)
                                                ! ... and with the scattering prefactor which includes
                                                !     a factor of 4 Pi, a relativistic correction and other
                                                !     constants of the Coulomb interaction.
                                                ! Note: the ionic part is added here to the real-part of the
                                                !        atomic form factor. It is already included in the
                                                !        numeric integration of absorptive form factors.
                                                !
      !
      lcw(i,j) = apdmp*(cstrfe + cstrfi)        ! store the combined structure factors of the screened
                                                ! potential and the ionic charge potential, bandwidth limited by aperture
      !
      !if (CS_doconsolemsg>0) then
      !  ! update progress indicator
      !  if (npc < int( 20.0*real(i+(j-1)*ny)/real(nx*ny) ) ) then
      !    npc = npc + 1
      !    write(unit=6,fmt='(A,$)') "."
      !  end if
      !end if
      !
    end do ! loop columns (wy-axis, transposed)
  end do ! loop rows (wx-axis, transposed)
  !
  !if (CS_doconsolemsg>0) then
  !  write(unit=6,fmt='(A)') ". "
  !end if
  !
  ! transform back to real space
  !
  select case (nft)
  case  (128)
    call  ODDCC128S(lcw,nx,ny,'BACK')
  case  (256)
    call  ODDCC256S(lcw,nx,ny,'BACK')
  case  (512)
    call  ODDCC512S(lcw,nx,ny,'BACK')
  case (1024)
    call ODDCC1024S(lcw,nx,ny,'BACK')
  case (2048)
    call ODDCC2048S(lcw,nx,ny,'BACK')
  case (4096)
    call ODDCC4096S(lcw,nx,ny,'BACK')
  case (8192)
    call ODDCC8192S(lcw,nx,ny,'BACK')
  end select
  !
  ! transfer to output array, use periodic repeat by (nrx, nry)
  !
  !CS_lastslice_potmax_im = 0.0
  do j=1, dimy
    j1 = modulo(j-1,ny)+1
    do i=1, dimx
      i1 = modulo(i-1,nx)+1
      pot(i,j) = lcw(i1,j1)*pscal
      !CS_lastslice_potmax_im = max(CS_lastslice_potmax_im,imag(pot(i,j)))
    end do
  end do
  !
  ! finish off, deallocate all arrays
  !  
10 continue
  if (allocated(agx)) deallocate(agx,stat=nerr)
  if (allocated(agy)) deallocate(agy,stat=nerr)
  if (allocated(fld)) deallocate(fld,stat=nerr)
  if (allocated(ild)) deallocate(ild,stat=nerr)
  if (allocated(jld)) deallocate(jld,stat=nerr)
  if (allocated(uatl)) deallocate(uatl,stat=nerr)
  if (allocated(lcw)) deallocate(lcw,stat=nerr)
  return
  !
  ! error handling
  !
13 nerr=-1
  call CS_ERROR("Potential memory not allocated.")
  goto 10
14 nerr=-1
  call CS_ERROR("Invalid parameters for potential calculations.")
  goto 10
15 nerr=-1
  call CS_ERROR("Memory allocation failed.")
  goto 10
16 nerr=-1
  call CS_ERROR("Memory deallocation failed.")
  goto 10
  
  return
  
end subroutine CS_GETSLICE_POT
!**********************************************************************!



!**********************************************************************!
!
! CS_GETSLICE_POT2
!
! subroutine, calculates projected potential of slice n
!             applies thermal displacements in case of frozen-lattice
!             simulation -> generates 1 FL variant
!
! REQUIRES
! - supercell data to be allocated and set up
!   (e.g. call CS_LOAD_EMSCELL)
! - slice data data to be allocated and set up
!   (e.g. call CS_SETSLICE_EQUIDIST)
! - scattering amplitude data to be allocated and prepared
!   (call CS_PREPARE_SCATTAMPS)
!
! REMARKS
! - when using frozen lattice displacements, be aware that
!   damping of scattering factors by debye-waller factors
!   is senseless and wrong! So turn them off when calculating
!   the scatterung factors by CS_PREPARE_SCATTAMPS
!
! - when using frozen lattice generation, call InitRand() before
!
! - Adds the ionic charge potentials here instead of to the
!   screened potentials in the CS_PREPARE_SCATTAMPS routine.
!   This may help to avoid some numerical round-off problems.
!
! - In order to enable parallel computing calls, we allocate
!   a local complex*8 square array for the fourier transforms.
!   The size of the transformation array is adjusted to the
!   needs of the call, keeping the memory consumption as low
!   as possible.
!
! - This is a second version, for testing the speed of different 
!   approaches. It is currently the prefered verions. Use it!
!
! INPUT:
!   integer*4 :: nslc               = slice number
!   integer*4 :: nx, ny             = discretisation of supercell axes
!   integer*4 :: nrx, nry           = repetition of supercell in slice
!   integer*4 :: nfl                = create frozen latice displacements
!   integer*4 :: ndw                = apply Debye-Waller factors
!   real*4 :: wl                    = electron wavelength [nm]
!
! IN/OUTPUT:
!   complex*8 :: pot(nx*nrx,ny*nry) = projected potential of slice
!   integer*4 :: nerr               = error code
!
subroutine CS_GETSLICE_POT2(nslc, nx, ny, nrx, nry, nfl, ndw, wl, pot, nerr)

  implicit none
  
  integer*4, intent(in) :: nslc, nx, ny, nrx, nry, nfl, ndw
  real*4, intent(in) :: wl
  integer*4, intent(inout) :: nerr
  complex*8, intent(inout) :: pot(nx*nrx,ny*nry)
  
  integer*4 :: i, j, i1, j1, ia, ia2, ja, jptr ! iterators
  integer*4 :: lndw
  integer*4 :: na ! atom count
  integer*4 :: dimx, dimy, nyqx, nyqy ! sampling numbers
  integer*4 :: infl ! internal frozen lattice flag
  integer*4, dimension(:), allocatable :: jld ! hash table to access slice IDs from cell IDs
  real*4 :: sx, sy, itogx, itogy ! sampling constants
  real*4 :: ht ! high tension in kV
  real*4 :: dx, dy ! atom position
  real*4 :: dwc, fldamp, biso ! displacement calculation helpers
  real*4 :: pocc ! atom occupancy
  real*4 :: dsz, vol ! slice thickness and volume in nm^3
  real*4 :: pscal ! potential scaling factor
  real*4 :: apdmp ! aperture function (damping factor at outer Fourier-space perimeter >2/3 gmax)
  real*4 :: crgio ! ionic charge
  real*4, dimension(:,:), allocatable :: lxy ! list of calculated coordinates
  complex*8 :: cval0, cval ! some complex vars
  complex*8, dimension(:,:), allocatable :: lcw, lcwrs ! local working array for FFTs
  
  logical :: dwflg

  ! external routines for 2D Fourier transform (link FFTs.f)
  external :: ODDCC128S, ODDCC256S,  ODDCC512S,  ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S
  real*4, external :: UniRand, GaussRand, getdwf
  
  !
  ! preset error flag
  !
  nerr = 0
  apdmp = 1.0 ! init without dampening
  
  !
  ! check mem alloc
  !
  if (.not.(CS_slicemem_allocated.and. &
     &      CS_cellmem_allocated.and. &
     &      CS_scattamp_prepared)) goto 13

  
  !
  ! check parameters
  !
  if (nx<=0.or.nrx<=0.or.nry<=0.or.ny<=0) goto 14
  if (wl<=0.0) goto 14
  if (nslc<1.or.nslc>CS_nspsc) goto 14
  if (CS_useppot==0) goto 17
  
  !
  ! handle frozen lattice flag
  !
  infl = nfl
  if (infl<0) infl = 0
  if (infl>1) infl = 1
  fldamp = 0.0
  dwc = 0.0
  lndw = ndw
  if (infl==1) then
    lndw = 0
    fldamp = CS_rr8p2
  end if
  !
  ! handle DWF flag (for ionic potential calculations only)
  !
  if (lndw==0) then
    dwflg = .FALSE.
  else
    dwflg = .TRUE.
  end if
  
  
  !
  ! init sampling
  !
  dimx = nx*nrx ! number of potential output samples (x)
  dimy = ny*nry ! number of potential output samples (y)
  nyqx = nx/2 ! Nyquist number (x) for potential generation
  nyqy = ny/2 ! Nyquist number (y) for potential generation
  sx = CS_scsx/real(nx) ! real-space sampling rate (x) for potential generation
  sy = CS_scsy/real(ny) ! real-space sampling rate (y) ...
  itogx = 1.0 / CS_scsx ! reciprocal space sampling rate (x) for potential generation
  itogy = 1.0 / CS_scsy ! reciprocal space sampling rate (y) for potential generation
  cval0 = cmplx(0.0,0.0)
  cval = cval0
  allocate(lcw(ny,nx), lcwrs(nx,ny),stat=nerr)
  if (nerr/=0) goto 15
  lcw = cmplx(0.0,0.0) ! this is the working array
  lcwrs = cmplx(0.0,0.0) ! this is the working array (real-space)
  ! CS_scagx and CS_scagy hold reciprocal space coordinates
  
  !
  ! save sampling as global
  !
  CS_sdimx = dimx
  CS_sdimy = dimy
  CS_repx = nrx
  CS_repy = nry
  CS_sampx = sx
  CS_sampy = sy
  
  !
  ! other parameters
  !
  ht = CS_WL2HT(wl)                                 ! high tension in kV
  na = CS_slcatnum(nslc)                            ! number of atoms in slice
  dsz = CS_slczlim(2,nslc)-CS_slczlim(1,nslc)       ! thickness of the slice in nm
  vol = CS_scsx*CS_scsy*dsz                         ! supercell volume in nm^3
  pscal = CS_v0 / vol                               ! potential scaling factor Vamp0 / Volume
  
  
  !
  ! clear output potential array
  !
  pot = cval0
  
  !
  ! return in case of empty slice
  !
  if (na<=0) then 
    goto 10
  end if
  
  !
  ! calculate super-cell potential in fourier space
  ! !!! fourierspace is scrambled and transposed
  allocate( lxy(2,na), jld(CS_numat) , stat=nerr)
  if (nerr/=0) goto 15
  lxy = 0.0
  jld = 0
  ia2 = 0
  !
  !call CS_PROG_START(na,1.0)
  !
  do ia=1, na ! loop ia over all atoms in slice
    !
    ja = CS_slcatacc(ia,nslc)               ! get atom index in super cell
    jld(ja) = ia                            ! store slice ID in list of cell IDs
    jptr = CS_scampptr(ja)                  ! get atom type ID in scattering data
    lxy(1,ia) = CS_atpos(1,ja)*CS_scsx      ! get atom avg. x-position
    lxy(2,ia) = CS_atpos(2,ja)*CS_scsy      ! get atom avg. y-position
    biso = CS_atdwf(ja)                     ! get debye-waller parameter (nm**2) = Biso
    dwc = sqrt(biso)                        ! get debye-waller parameter from (nm**2) to (nm) (sqrt(Biso))
    if (infl==1) then                     ! dice frozen lattice displacements (x,y) and apply
      if (CS_atlnk(ja)==0) then ! independent atom site
        dx = fldamp*dwc*GaussRand()         !   fldamp*dwc = sqrt( Biso / 8*pi^2) = sqrt( <u_s^2> )
        dy = fldamp*dwc*GaussRand()
        lxy(1,ia) = lxy(1,ia) + dx          ! add random x displacement
        lxy(2,ia) = lxy(2,ia) + dy          ! add random y displacement
      else ! dependent atom site
        ia2 = jld(CS_atlnk(ja))
        lxy(1,ia) = lxy(1,ia2)
        lxy(2,ia) = lxy(2,ia2)
      end if
    end if
    pocc = CS_atocc(ja)                     ! Occupancy
    crgio = CS_atcrg(ja)                    ! ionic charge
    !
    if ( crgio==0.0  .or. CS_useextsca==0 ) then ! neutral atoms (faster)
      ! accumulate                      occupancy * form-factor (core) (DWF is handled already in CS_scaff2d) 
      lcw(1:ny,1:nx) = lcw(1:ny,1:nx) + pocc * CS_scaff2d(1:ny,1:nx,jptr) &
                     & * cexp( lxy(1,ia)*CS_scagn(1:ny,1:nx,1) &
                     &       + lxy(2,ia)*CS_scagn(1:ny,1:nx,2) ) ! * translation phase factor
      !
    else ! ion (slower, need to multiply ionic potential)
      !
      if (dwflg) then ! apply DWF to ionic part (slower)
        ! accumulate                      occupancy * ( form-factor (core)
        lcw(1:ny,1:nx) = lcw(1:ny,1:nx) + pocc * ( CS_scaff2d(1:ny,1:nx,jptr) &
                       &   + crgio* CS_scadwf(1:ny,1:nx,jptr) * CS_scagn(1:ny,1:nx,3) ) & ! + form-factor (ionic charge) with DWF )
                       & * cexp( lxy(1,ia)*CS_scagn(1:ny,1:nx,1) &
                       &       + lxy(2,ia)*CS_scagn(1:ny,1:nx,2) ) ! * translation phase factor
      else ! no DWF in ionic part (faster)
        ! accumulate                      occupancy * ( form-factor (core)
        lcw(1:ny,1:nx) = lcw(1:ny,1:nx) + pocc * ( CS_scaff2d(1:ny,1:nx,jptr) &
                       &                         + crgio*CS_scagn(1:ny,1:nx,3) ) & ! + form-factor (ionic charge) )
                       & * cexp( lxy(1,ia)*CS_scagn(1:ny,1:nx,1) &
                       &       + lxy(2,ia)*CS_scagn(1:ny,1:nx,2) ) ! * translation phase factor
      end if
    end if
    !
    ! update progress indicator
    !call CS_PROG_UPDATE(ia)
      !
  end do ! loop ia over all atoms in slice
  !
  !call CS_PROG_STOP(na)
  !
  ! transform back to real space
  !
  call CS_CCFFT2D(lcwrs,lcw,nx,ny,-1)
  !
  ! transfer to output array, use periodic repeat by (nrx, nry)
  !
  !CS_lastslice_potmax_im = 0.0
  do j=1, dimy
    j1 = modulo(j-1,ny)+1
    do i=1, dimx
      i1 = modulo(i-1,nx)+1
      pot(i,j) = lcwrs(i1,j1)*pscal
      !CS_lastslice_potmax_im = max(CS_lastslice_potmax_im,imag(pot(i,j)))
    end do
  end do
  !
  ! finish off, deallocate all arrays
  !  
10 continue
  if (allocated(lcw)) deallocate(lcw,stat=nerr)
  if (allocated(lcwrs)) deallocate(lcwrs,stat=nerr)
  if (allocated(jld)) deallocate(jld,stat=nerr)
  if (allocated(lxy)) deallocate(lxy,stat=nerr)
  return
  !
  ! error handling
  !
13 nerr=-1
  call CS_ERROR("Potential memory not allocated.")
  goto 10
14 nerr=-1
  call CS_ERROR("Invalid parameters for potential calculations.")
  goto 10
15 nerr=-1
  call CS_ERROR("Memory allocation failed.")
  goto 10
16 nerr=-1
  call CS_ERROR("Memory deallocation failed.")
  goto 10
17 nerr=-1
  call CS_ERROR("2D projected form factors not prepared.")
  goto 10
  
  return
  
end subroutine CS_GETSLICE_POT2
!**********************************************************************!






!**********************************************************************!
!
! CS_GETSLICE_PGR
!
! subroutine, calculates phase grating of slice n
!
! INPUT:
!   integer*4 :: nslc               = slice number
!   integer*4 :: nx, ny             = discretisation of supercell axes
!   integer*4 :: nrx, nry           = repetition of supercell in slice
!   integer*4 :: nabs               = flag for usage of absorption data
!   integer*4 :: nfl                = flag for creation of frozen lattice
!
! IN/OUTPUT:
!   complex*8 :: pgr(nx*nrx,ny*nry) = phase grating
!   integer*4 :: nerr               = error code
!
subroutine CS_GETSLICE_PGR(nslc, nx, ny, nrx, nry, nabs, nfl, ndw, wl, pgr, nerr)

  implicit none
  
  integer*4, intent(in) :: nslc, nx, ny, nrx, nry, nabs, nfl, ndw
  real*4, intent(in) :: wl
  integer*4, intent(inout) :: nerr
  complex*8, intent(inout) :: pgr(nx*nrx,ny*nry)
  
  integer*4 :: nalloc
  integer*4 :: i, j, i1, j1, i2, j2 ! iterators
  integer*4 :: dimx, dimy ! final dimensions
  integer*4 :: nyqx, nyqy ! nyquist numbers
  !real*4 :: relati ! relativistic correction
  real*4 :: sigmae ! interaction constant
  real*4 :: dz ! slice thickness in nm
  real*4 :: ht ! high tension [kV]
  real*4 :: pscal ! potential scaling
  real*4 :: rV ! potential data
  real*4 :: rmaxphase ! maximum phase shift
  real*4 :: relati ! relativistic correction
  real*4 :: i2gx, i2gy, gx, gy, gx2, g2, gthr2
  complex*8 :: cpot ! potential value
  complex*8 :: cpgr ! phase factor
  complex*8 :: cabsorp ! absorption factor
  complex*8, allocatable, dimension(:,:) :: pot ! potential
  complex*8, allocatable, dimension(:,:) :: pgrft ! fourier transform of the phase grating
  
  
  !
  ! init parameters of the calculation
  !
  dimx = nx*nrx
  dimy = ny*nry
  nyqx = (dimx-modulo(dimx,2))/2
  nyqy = (dimy-modulo(dimy,2))/2
  ht = CS_WL2HT(wl)                                 ! high tension in kV
  dz = abs(CS_slczlim(1,nslc)-CS_slczlim(2,nslc))   ! slice thickness
  relati = 1.0 + ht / CS_elm0                       ! relativistic correction
  sigmae = CS_sig * wl                              ! interaction constant
  pscal = sigmae * dz                               ! product of interaction constant and projected thickness
  
  !
  ! get projected potential
  !
  allocate(pot(nx,ny),stat=nalloc)
  if (nalloc/=0) then
    call CS_ERROR("Failed to allocate memory for projected potential.")
    goto 13
  end if
  pot = cmplx(0.0,0.0)
  !                   (nslc, nx, ny, nrx, nry, nfl, wl, pot, nerr)
  !call CS_GETSLICE_POT(nslc, nx, ny, 1, 1, nfl, ndw, wl, pot, nerr)
  call CS_GETSLICE_POT2(nslc, nx, ny, 1, 1, nfl, ndw, wl, pot, nerr)
  if (nerr/=0) then
    deallocate(pot,stat=nalloc)
    goto 13
  end if
  
  !
  ! optional potential backup, memory allocations
  !
  if (CS_backup_pot_flg==1) then
    if (allocated(CS_backup_pot)) then ! deallocate previous potential backup memory
      deallocate(CS_backup_pot,stat=nalloc)
    end if
    ! allocate new potentail backup memory
    allocate(CS_backup_pot(dimx,dimy),stat=nalloc)
    if (nalloc/=0) then
      call CS_ERROR("Failed to allocate memory for projected potential backup.")
      goto 13
    end if
    CS_backup_pot = cmplx(0.0,0.0)
  end if
  
  !
  ! get parameters
  !
  rmaxphase = 0.0
  CS_useabsorption = 0
  cabsorp = cmplx(0.0,1.0) * pscal
  if (nabs/=0) then 
    CS_useabsorption = 1
    ! cabsorp = cmplx(-CS_absorptionprm,1.0) * pscal, obsolete, 100128 JB
  end if
  
  !
  ! generate the phase grating // apply periodic repeats
  !
  ! pgr = EXP( ii*pot*pscal )
  !       pscal = sigmae * dz = psig * wl * dz
  !               psig = 2*pi*m0*e / h**2 *10e-18 nm-2 = 2.088656
  !               wl = electron wavelength [nm]
  !               dz = slice thickness [nm]
  !
  do j=1, ny
    do i=1, nx
      ! calculate exp( ii*potential )
      cpot = pot(i,j)*cabsorp ! cpot = ii*potential
      cpgr = exp(cpot) ! phase grating = exp( ii*potential )
      ! store dbg info on max phase
      rV = real(pot(i,j))*pscal ! get real part of potential to check max. phase shift
      rmaxphase = max(rmaxphase,abs(rV)) ! save max. phase shift
      ! repeater loop
      do j1=0, nry-1
        j2 = j+ny*j1
        do i1=0, nrx-1
          i2 = i + nx*i1
          pgr(i2,j2) = cpgr ! set phase grating value to output array
          if (CS_backup_pot_flg==1) then ! store projected potential
            ! store potential after removing the relativistic correction
            ! relati = 1.0 + ht [keV] / 511 [keV]
            CS_backup_pot(i2,j2) = pot(i,j)/relati
          end if
        end do
      end do
    end do
  end do
  
  ! update the maximum phase accumulator
  CS_maxphase = max(CS_maxphase, rmaxphase)
  ! REMOVED 181128, v0.70, JB -> changed to accumulation of max. phase
  !if (rmaxphase>CS_pi4) then ! check max. phase
  !  CS_warn_num = CS_warn_num + 1
  !  write(unit=CS_warnmsg,fmt='(A,G13.4,A)') "Diffraction phase-shift of current slice exceeds pi/4 (",rmaxphase,")."
  !  call CS_MESSAGE("Warning: "//trim(CS_warnmsg))
  !end if
  
  deallocate(pot,stat=nalloc)
  
  if (CS_do_bwl_pgr) then ! apply bwl to phase-grating
    allocate(pgrft(dimy,dimx),stat=nalloc)
    if (nerr/=0) goto 13
    pgrft = cmplx(0.0,0.0)
    call CS_CCFFT2D(pgr, pgrft, dimx, dimy, 1) ! forward FFT from pgr to pgrft
    i2gx = 1. / (real(nrx) * CS_scsx) ! fourier space sampling rate x
    i2gy = 1. / (real(nry) * CS_scsy) ! ... y
    gthr2 = ( min(i2gx*real(nyqx),i2gy*real(nyqy))*CS_PROSIZE_THR)**2 ! g2 aperture threshold
    do j=1, dimx ! x is with rows
      j1 = mod((j+nyqx-1),dimx)-nyqx
      gx = i2gx*real(j1)
      gx2 = gx*gx
      do i=1, dimy ! y is with columns
        i1 = mod((i+nyqy-1),dimy)-nyqy
        gy = i2gy*real(i1)
        g2 = gx2 + gy*gy
        if (g2>gthr2) then ! apply hard aperture
          pgrft(i,j) = cmplx(0.0,0.0)
        end if
      end do
    end do
    call CS_CCFFT2D(pgr, pgrft, dimx, dimy, -1) ! backward FFT from pgrft to pgr
    deallocate(pgrft, stat=nalloc)
  end if
  
  
  return
  
  !
  ! handle errors
  !
13 nerr=-1
  call CS_ERROR("Failed to allocate phase-grating memory.")
  return
  
end subroutine CS_GETSLICE_PGR
!**********************************************************************!





!**********************************************************************!
!
! CS_GET_SAMESITE
!
! subroutine, Determines the first position in a list for which the
!             distance to a test position is below a given threshold.
!             Poistions are given in fractional coordinates, but the
!             coordinate scale of an orthorhombic system is used to
!             calculate a real distance.
!             Periodic wrap-around is considered in the distance check.
!
! INPUT:
!   real*4 :: fpos(3)     = test position
!   real*4 :: l_fpos(:,:) = list of fractional coordinates
!   real*4 :: scal(3)     = coordinate scale
!   real*4 :: ssthr       = distance threshold for close sites
!   integer*4 :: idxss    = index of the first close site
!   real*4 :: rdss        = distance of the first close site (scaled)
!   integer*4 :: nerr     = error code (0 = success)
!
subroutine CS_GET_SAMESITE(fpos, l_fpos, scal, ssthr, idxss, rdss, nerr)

  use symops

  implicit none
  
! interface variables
  real*4, intent(in) :: fpos(3), l_fpos(:,:), scal(3), ssthr
  integer*4, intent(out) :: idxss, nerr
  real*4, intent(out) :: rdss
  
! local parameters
  
! local vaiables
  integer*4 :: i, n
  real*4 :: rdist, v_fdif(3), v_dif(3)
  
! initialization
  nerr = 0
  idxss = 0
  rdss = -1.0
  n = SIZE(l_fpos,2)
  if (3/=SIZE(l_fpos,1)) goto 802
  if (scal(1)<=0.0 .or. scal(2)<=0.0 .or. scal(3)<=0.0) goto 803
  rdist = 0.0
  v_fdif = 0.0
  v_dif = 0.0
  
! test
  do i=1, n
    v_fdif(1:3) = modulo( l_fpos(1:3,i) - fpos(1:3) + 0.5 , 1.0 ) - 0.5
    v_dif = v_fdif*scal
    rdist = sqrt( sum( v_dif*v_dif ) )
    if (rdist<ssthr) then
      rdss = rdist
      idxss = i
      goto 800
    end if
  end do
  
800 continue
  return
  
! error handling
802 continue
  nerr = 2
  call CS_ERROR("Invalid dimension (1) of input array (l_fpos).")
  return
803 continue
  nerr = 3
  call CS_ERROR("Invalid scale of orthorhombic structure (scal).")
  return
  
end subroutine CS_GET_SAMESITE





!**********************************************************************!
!
! CS_GET_NNAT
!
! subroutine, Determines the nearest neighbor 3D position from a list
!             of fractional coordinates, considering the coordinate
!             scale of an orthorhombic system.
!             Periodic wrap-around is considered in the distance.
!
! INPUT:
!   real*4 :: fpos(3)     = test position
!   real*4 :: l_fpos(:,:) = list of fractional coordinates
!   real*4 :: scal(3)     = coordinate scale
!   integer*4 :: l_use(:) = list of postions in use (use=1, ignore=0)
!   integer*4 :: idxnn    = index of the nearest neighbor
!   real*4 :: rdnn        = distance of the nearest neihbor (scaled)
!   integer*4 :: nerr     = error code (0 = success)
!
subroutine CS_GET_NNAT(fpos, l_fpos, scal, l_use, idxnn, rdnn, nerr)

  use symops

  implicit none
  
! interface variables
  real*4, intent(in) :: fpos(3), l_fpos(:,:), scal(3)
  integer*4, intent(in) :: l_use(:)
  integer*4, intent(out) :: idxnn, nerr
  real*4, intent(out) :: rdnn
  
! local parameters
  
! local vaiables
  integer*4 :: i, n
  real*4 :: rdist, v_fdif(3), v_dif(3), rdmin
  
! initialization
  nerr = 0
  idxnn = 0
  rdnn = 0.0
  n = SIZE(l_fpos,2)
  if (n/=SIZE(l_use,1)) goto 801
  if (3/=SIZE(l_fpos,1)) goto 802
  if (scal(1)<=0.0 .or. scal(2)<=0.0 .or. scal(3)<=0.0) goto 803
  rdist = 0.0
  v_fdif = 0.0
  v_dif = 0.0
  rdmin = -1.0
  
! test
  do i=1, n
    if (l_use(i)/=1) cycle ! skip not used positions
    v_fdif(1:3) = modulo( l_fpos(1:3,i) - fpos(1:3) + 0.5 , 1.0 ) - 0.5
    v_dif = v_fdif*scal
    rdist = sqrt( sum( v_dif*v_dif ) )
    if (rdmin<0.0 .or. rdmin>rdist) then
      rdmin = rdist
      idxnn = i
    end if
  end do
  rdnn = rdmin
  
800 continue
  return
  
! error handling
801 continue
  nerr = 1
  call CS_ERROR("Inconsistent length of input arrays (l_fpos, l_use).")
  return
802 continue
  nerr = 2
  call CS_ERROR("Invalid dimension (1) of input array (l_fpos).")
  return
803 continue
  nerr = 3
  call CS_ERROR("Invalid scale of orthorhombic structure (scal).")
  return
  
end subroutine CS_GET_NNAT




!**********************************************************************!
!
! CS_ORIENT_CELL
!
! subroutine, transforms the super-cell to new projection axes
!
! INPUT:
!   real*4 :: oh,ok,ol              = new z in the basis of the current
!                                     super-cell
!   real*4 :: yu,yv,yw              = new y in the basis of the current
!                                     super-cell
!   real*4 :: sa,sb,sc              = new super-cell size in nm
!   
! IN/OUTPUT:
!   integer*4 :: nerr               = error code
!
subroutine CS_ORIENT_CELL(oh,ok,ol, yu,yv,yw, sa,sb,sc, nerr)

  use symops

  implicit none
  
! interface variables
  real*4, intent(in) :: oh,ok,ol
  real*4, intent(in) :: yu,yv,yw
  real*4, intent(in) :: sa,sb,sc
  integer*4, intent(inout) :: nerr
  
! local parameters
  real*8, parameter :: pos_prec = 1.D-5 ! position precision (rounding)
  real*8, parameter :: min_blen_nm = 1.D-4 ! minimum length of basis vectors in nm
  real*8, parameter :: min_latdev_tol = 1.D-2 ! tolerated deviation of lattices (warning trigger)
  real*4, parameter :: min_posdif = 0.02 ! 0.2 A distance between atomic positions?

! local vaiables
  integer*4 :: ierr ! local error code
  integer*4 :: nalloc ! allocation status
  integer*4 :: i, j, k, l, n, m, n3 ! some iterators and numbers
  integer*4, dimension(3) :: nitf ! number of old cells in new cell per dimension (in old basis)
  integer*4, dimension(2,3) :: itf ! space filling iterator
  real*8 :: dtmp1, dtmp2, dtmp3 ! some temporary values
  real*8, dimension(3,3) :: base_0, base_0_inv ! old SC basis and inverse
  real*8, dimension(3,3) :: base_1, base_1_inv ! new SC basis and inverse
  real*8, dimension(3,3) :: m_tmp ! some 3x3 matrix
  real*8, dimension(3,1) :: v_tmp1, v_tmp2, v_tmp3 ! some 3D vectors
  real*8, dimension(3,8) :: sc1_crn ! corners of the new super-cell
  real*8, dimension(3,8) :: sc1_crn_dev ! deviations of the new super-cell corners from the old lattice
!
! temporaray atomic data 
!     atomic type
  character(len=CS_atsl), allocatable, dimension(:) :: l_attype
  integer*4, allocatable, dimension(:) :: l_atnum
  integer*4, allocatable, dimension(:) :: l_atuse
  real*4, allocatable, dimension(:) :: l_atcrg
  real*4, allocatable, dimension(:,:) :: l_atpos
  real*4, allocatable, dimension(:) :: l_atocc
  real*4, allocatable, dimension(:) :: l_atdwf
  real*4 :: rtmp(3), rtmp1 ! some temporary values
!
  
! initialization
  ierr = 0
  nerr = 0
  nalloc = 0
  
!  open ( unit=25, file="E:\Simulationen\bin\drprobe-orient-debug.log", &
!       & action='WRITE', status='REPLACE', form='FORMATTED' )
!  write(unit=25,fmt='(A)') "Log of routine CS_ORIENT_CELL"
  
!  write(unit=25,fmt='(A)') "- init matrix"
  base_0 = 0.D+0
  base_0_inv = 0.D+0
  call CONVMAT(dble(CS_scsx), dble(CS_scsy), dble(CS_scsz), &
     & dble(CS_scalpha), dble(CS_scbeta), dble(CS_scgamma), base_0)
  base_0 = TRANSPOSE(base_0) ! need to transpose because CONVMAT is column oriented
  call CHOP_MAT(base_0, 1.D-6)
!  write(unit=25,fmt='(A)') "- matrix inversion"
  call INVMAT(base_0, base_0_inv)
  
  

! check input
  if (sa<min_blen_nm .or. sb<min_blen_nm .or. sc<min_blen_nm) goto 801
  
! get the new projection vector in cartesian coordinates
!  write(unit=25,fmt='(A)') "- init projection vector"
! this will be our new z axis of the super-cell
  v_tmp1 = 0.D+0
  v_tmp1(1:3,1) = (/ oh,ok,ol /)
  v_tmp3 = 0.D+0
  v_tmp3 = MATMUL(base_0, v_tmp1)
  dtmp3 = VECLENGTH(v_tmp3(1:3,1))
  if (dtmp3<min_blen_nm) goto 802 ! error
  v_tmp3 = v_tmp3 / dtmp3 ! normalize V3 to N3
! get the new projected y axis in cartesian coordinates
!  write(unit=25,fmt='(A)') "- init projection y-orientation vector"
  v_tmp1 = 0.D+0
  v_tmp1(1:3,1) = (/ yu,yv,yw /)
  v_tmp2 = 0.D+0
  v_tmp2 = MATMUL(base_0, v_tmp1)
  dtmp2 = VECLENGTH(v_tmp2(1:3,1))
! determine the in-plane component of the chosen y axis
  dtmp1 = DOT_PRODUCT(v_tmp3(1:3,1), v_tmp2(1:3,1))
  v_tmp2 = v_tmp2 - v_tmp3*dtmp1 ! V2 = V2_0 - N3*(V2_0.N3)
  dtmp2 = VECLENGTH(v_tmp2(1:3,1))
! check the length of the new y-axis
  if (dtmp2<min_blen_nm) then
    call CS_MESSAGE("Warning: The chosen new projected y-axis is invalid.")
    call CS_MESSAGE("         choosing a new y-axis ...")
    ! - try the cartesian y-axis
    v_tmp2(1:3,1) = (/ 0.D+0, 1.D+0, 0.D+0 /)
    dtmp1 = DOT_PRODUCT(v_tmp3(1:3,1), v_tmp2(1:3,1))
    v_tmp2 = v_tmp2 - v_tmp3*dtmp1 ! V2 = V2_0 - N3*(V2_0.N3)
    dtmp2 = VECLENGTH(v_tmp2(1:3,1))
    ! check again ...
    if (dtmp2<min_blen_nm) then
      ! - try the cartesian z-axis
      v_tmp2(1:3,1) = (/ 0.D+0, 0.D+0, 1.D+0 /)
      dtmp1 = DOT_PRODUCT(v_tmp3(1:3,1), v_tmp2(1:3,1))
      v_tmp2 = v_tmp2 - v_tmp3*dtmp1 ! V2 = V2_0 - N3*(V2_0.N3)
      dtmp2 = VECLENGTH(v_tmp2(1:3,1))
    end if
    ! one of the two should always produce some non-zero component in the projection plane
  end if
! normalize the new Y vector
  v_tmp2 = v_tmp2 / dtmp2
! calculate the new X axis from the cross product N2 x N3 = N1
!  write(unit=25,fmt='(A)') "- init projection x-orientation vector"
  v_tmp1 = 0.D+0
  call INN_PRODUCT( v_tmp2, v_tmp3, v_tmp1 )
! set the new basis matrix (row oriented)
!  write(unit=25,fmt='(A)') "- set new basis vectors"
  base_1 = 0.D+0
  base_1(1:3,1) = v_tmp1(1:3,1)*dble(sa)
  base_1(1:3,2) = v_tmp2(1:3,1)*dble(sb)
  base_1(1:3,3) = v_tmp3(1:3,1)*dble(sc)
  call CHOP_MAT(base_1, 1.D-6)
! get its inverse
!  write(unit=25,fmt='(A)') "- determine inverse of the new basis"
  base_1_inv = 0.D+0
  call INVMAT(base_1, base_1_inv)

! get the coordinates of the new super-cell corners in the old basis
! these coordinates determine the repeats of the original cell in the
! new cell.
!  write(unit=25,fmt='(A)') "- determine projection range by cell corner locations"
  sc1_crn = 0.D+0
  sc1_crn(1:3,2:2) = MATMUL(base_0_inv, base_1(1:3,1:1))
  sc1_crn(1:3,3:3) = MATMUL(base_0_inv, base_1(1:3,2:2))
  sc1_crn(1:3,4:4) = MATMUL(base_0_inv, base_1(1:3,3:3))
  sc1_crn(1:3,5:5) = MATMUL(base_0_inv, base_1(1:3,1:1)+base_1(1:3,2:2))
  sc1_crn(1:3,6:6) = MATMUL(base_0_inv, base_1(1:3,1:1)+base_1(1:3,3:3))
  sc1_crn(1:3,7:7) = MATMUL(base_0_inv, base_1(1:3,2:2)+base_1(1:3,3:3))
  sc1_crn(1:3,8:8) = MATMUL(base_0_inv, base_1(1:3,1:1)+base_1(1:3,2:2)+base_1(1:3,3:3))
! get the deviation from the old lattice
  sc1_crn_dev = 0.D+0
  sc1_crn_dev = sc1_crn - dble((NINT(sc1_crn)))
! check the total deviation and calculate space filling iterators (itf)
! due to round-off errors, the iterator may be by one cell larger or
! smaller than actually required.
  dtmp1 = 0.D+0
  itf = 0 ! puh, we can start from zero, since we keep the same origin
  do i=1, 8
    dtmp1 = dtmp1 + VECLENGTH(sc1_crn_dev(1:3,i))
!    write(unit=25,fmt='(A,I1)') "  - corner #",i
    do j=1, 3
      itf(1,j) = MIN(itf(1,j),FLOOR(sc1_crn(j,i)))
      itf(2,j) = MAX(itf(2,j),CEILING(sc1_crn(j,i)))
    end do
!    write(unit=25,fmt='(A,3G15.5)') "    - projection: ", sc1_crn(1:3,i)
!    write(unit=25,fmt='(A,3I3)')    "    - min. range: ", itf(1,1:3)
!    write(unit=25,fmt='(A,3I3)')    "    - max. range: ", itf(2,1:3)
  end do
  if (dtmp1>min_latdev_tol) then
    call CS_MESSAGE("Warning: The new super-cell lattice deviates from the old lattice.")
    call CS_MESSAGE("         The resulting super-cell may not define a periodic structure.")
  end if
! get the number of repeats of the old cell in the new cell
  nitf(1) = itf(2,1) - itf(1,1)
  nitf(2) = itf(2,2) - itf(1,2)
  nitf(3) = itf(2,3) - itf(1,3)
  n3 = nitf(1)*nitf(2)*nitf(3)
  if (n3<=0) goto 803 ! each dimension should be represented
!
! When this point is reached, the new cell is valid.
!  write(unit=25,fmt='(A)') "+ found valid cell"
!
  if (CS_cellmem_allocated .and. CS_numat>0) then
    ! There are atoms in the current (old) super-cell.
    ! Fill the new super cell with atoms from the old cell.
    !
!    write(unit=25,fmt='(A)') "- filling new cell with atoms"
    !
    ! prepare memory for the new atom data
    ! - the required number of atoms (may be less in the end)
    n = n3*CS_numat
    m = 0 ! init the used number of atoms (what really is in the new cell)
    ! - allocate the temporary memory for the creation of the new structure
!    write(unit=25,fmt='(A,I6)') "- allocating local helper arrays, n = ", n
    allocate(l_attype(n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_attype = ""
    allocate(l_atnum(n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_atnum = 0
    allocate(l_atcrg(n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_atcrg = 0.0
    allocate(l_atpos(1:3,n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_atpos = 0.0
    allocate(l_atocc(n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_atocc = 0.0
    allocate(l_atdwf(n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_atdwf = 0.0
    allocate(l_atuse(n), stat=nalloc)
    if (nalloc/=0) goto 804
    l_atuse = 0
    
    !
    ! Create the new atomic positions by periodic repeat of the old
    ! structure in the basis of the new structure.
    ! And by clipping the positions to the fractional range [0,1[.
    !
    ! get the transformation matrix
!    write(unit=25,fmt='(A)') "- calculating transformation matrix"
    m_tmp = MATMUL(base_1_inv, base_0)
    !
!    write(unit=25,fmt='(A)') "- projecting all atoms"
    do i=1, CS_numat ! loop over all atoms of the old cell
      !
      ! repeater loops (dimensions of the old cell)
      do l=itf(1,3), itf(2,3)-1
        v_tmp1(3,1) = dble(l)
        do k=itf(1,2), itf(2,2)-1
          v_tmp1(2,1) = dble(k)
          do j=itf(1,1), itf(2,1)-1
            v_tmp1(1,1) = dble(j)
            v_tmp2(1:3,1) = dble(CS_atpos(1:3,i)) + v_tmp1(1:3,1) ! periodic repeat of atom (i) position in old basis
            v_tmp3 = MATMUL(m_tmp, v_tmp2) ! ... in the new basis
            call ROUND_ARR(v_tmp3(1:3,1), pos_prec) ! take care of round-off errors
            ! check whether the position is in the target cell
            if ( v_tmp3(1,1)>=0.D+0 .and. v_tmp3(1,1)<=1.D+0 .and. &
               & v_tmp3(2,1)>=0.D+0 .and. v_tmp3(2,1)<=1.D+0 .and. &
               & v_tmp3(3,1)>=0.D+0 .and. v_tmp3(3,1)<=1.D+0 ) then
              !
              ! yes, this is in the cell
              !
              m = m + 1
              l_attype(m) = CS_attype(i)
              l_atnum(m) = CS_atnum(i)
              l_atcrg(m) = CS_atcrg(i)
              l_atpos(1:3,m) = real( v_tmp3(1:3,1), kind=4 )
              l_atocc(m) = CS_atocc(i)
              l_atdwf(m) = CS_atdwf(i)
              l_atuse(m) = 1
!              write(unit=25,fmt='(A,I4,A,I4,A,3F8.4,A)') "  - atom #",m," Z = ",l_atnum(m)," at (",l_atpos(1:3,m),")"
            end if
          end do
        end do
      end do
      !
    end do
    !
    m = SUM(l_atuse)
!    write(unit=25,fmt='(A,I6,A)') "+ added ",m," atoms to the new cell"
    
    rtmp(1:3) = (/sa,sb,sc/) ! temp. store cell dimensions
!    write(unit=25,fmt='(A)') "- removing close atoms"
    if (m>1) then ! check for close atoms
      do i=1, m-1 ! check all atoms ( i = 1 .. m-1 ) agains all other atoms ( j = i+1 .. m ) ! ALWAYS !
        if (l_atuse(i)==0) then ! this atom is not used, skip it
          cycle ! skip inactive items
        end if
        do j=i+1, m ! loop over the "other" atoms (no self test)
          if (l_atuse(j)==0) then ! this atom is not used, skip it
            cycle
          end if
          !
          ! get the distance between atom position i and j
          v_tmp1(1:3,1) = l_atpos(1:3,i) - l_atpos(1:3,j)
          ! check for the shortest possibly distance within periodic boundary condistions
          v_tmp2(1:3,1) = modulo( v_tmp1(1:3,1) + 0.5 , 1.0 ) - 0.5
          v_tmp3(1:3,1) = v_tmp2(1:3,1) * rtmp ! scale the minimum distance to physical units (nm)
          rtmp1 = sqrt( sum( v_tmp3(1:3,1)*v_tmp3(1:3,1) ) ) ! calculate the distance
          ! now comes the check agains minimum required distance and atom type difference
          if (rtmp1<min_posdif .and. l_atnum(i)==l_atnum(j)) then ! two atoms of same type on same position
!          write(unit=25,fmt='(A,I5)') "  - removing atom #",j
            l_atuse(j) = 0 ! switch off atom j, keep atom i
          end if
          !
        end do ! loop over j
      end do ! loop over i
    end if
    
    m = SUM(l_atuse)
!    write(unit=25,fmt='(A,I6,A)') "+ ",m," atoms remaining in the new cell"
    
    !
    ! copy the content of the new cell to the module arrays
    ! - prepare the module arrays
    call CS_DEALLOC_CELLMEM(ierr)
    if (ierr/=0) goto 805
    
    call CS_ALLOC_CELLMEM(m, ierr) ! sets CS_numat=m
!    write(unit=25,fmt='(A,I6,A,I6)') "- CS_ALLOC_CELLMEM, m = ", m, ", err = ", ierr
    if (ierr/=0) goto 804
    !
    if (m>0) then
      ! - copy
!      write(unit=25,fmt='(A)') "- storing final atomic structure data"
      j = 0
      do i=1, n
        if (l_atuse(i)==1) then
          j = j + 1
          CS_attype(j) = l_attype(i)
          CS_atnum(j)  = l_atnum(i)
          CS_atcrg(j)  = l_atcrg(i)
          CS_atpos(1:3,j) = l_atpos(1:3,i)
          CS_atocc(j)  = l_atocc(i)
          CS_atdwf(j)  = l_atdwf(i)
        end if
      end do
      !
    end if
    !
    ! clean up the mess
!    write(unit=25,fmt='(A)') "- deallocating local helper arrays"
    deallocate(l_attype, stat=nalloc)
    if (nalloc/=0) goto 805
    deallocate(l_atnum, stat=nalloc)
    if (nalloc/=0) goto 805
    deallocate(l_atcrg, stat=nalloc)
    if (nalloc/=0) goto 805
    deallocate(l_atpos, stat=nalloc)
    if (nalloc/=0) goto 805
    deallocate(l_atocc, stat=nalloc)
    if (nalloc/=0) goto 805
    deallocate(l_atdwf, stat=nalloc)
    if (nalloc/=0) goto 805
    !
  end if
  !
  ! Set the new super-cell parameters
!  write(unit=25,fmt='(A)') "- storing new cell parameters"
  CS_scsx = VECLENGTH(base_1(1:3,1))
  CS_scsy = VECLENGTH(base_1(1:3,2))
  CS_scsz = VECLENGTH(base_1(1:3,3))
  CS_scalpha = 9.D+1
  CS_scbeta  = 9.D+1
  CS_scgamma = 9.D+1
!
! exit
800 continue
!  write(unit=25,fmt='(A,I4.4,A)') "End of CS_ORIENT_CELL (", nerr,")"
!  close ( unit=25 )
  n3 = 0
  return
  
! error handlings
801 continue
  nerr = 1
  call CS_ERROR("Invalid size of new super-cell (zero length).")
  goto 800
802 continue
  nerr = 2
  call CS_ERROR("Invalid projection vector (zero length).")
  goto 800
803 continue
  nerr = 3
  call CS_ERROR("Invalid size of the new super-cell (dimension).")
  goto 800
804 continue
  nerr = 4
  call CS_ERROR("Memory allocation failed.")
  goto 800
805 continue
  nerr = 5
  call CS_ERROR("Memory de-allocation failed.")
  goto 800
  
  return

end subroutine CS_ORIENT_CELL



!**********************************************************************!
!
! CS_SHIFT_ATOMS
!
! subroutine, translates all atoms in the super-cell by a given
!             fractional shift vector
!
! INPUT:
!   real*4 :: tx, ty, tz            = shift vector components along
!                                     the cell axes (fractional)
!   
! IN/OUTPUT:
!   integer*4 :: nerr               = error code
!
subroutine CS_SHIFT_ATOMS(tx,ty,tz, nerr)

  implicit none
  
! interface variables
  real*4, intent(in) :: tx,ty,tz
  integer*4, intent(inout) :: nerr
  
! local vaiables
  integer*4 :: ierr ! local error code
  integer*4 :: nalloc ! allocation status
  integer*4 :: i ! some iterators and numbers
  real*4 :: rtmp(3) ! some temporary values
  real*4 :: rshf(3)
!
  
! initialization
  ierr = 0
  nerr = 0
  nalloc = 0
  rshf(1) = tx
  rshf(2) = ty
  rshf(3) = tz
!
  if (CS_cellmem_allocated .and. CS_numat>0) then
    ! There are atoms in the current super-cell.
    !
    ! Create the new atomic positions by shifting
    ! all atoms and applying periodic boundary conditions.
    ! Final positions are wrapped modulo 1.0 to the fractional range [0,1[.
    !
    do i=1, CS_numat ! loop over all atoms
      !
      rtmp(1:3) = CS_atpos(1:3,i) + rshf
      !
      CS_atpos(1:3,i) = modulo( rtmp(1:3) , 1.0 )
      !
    end do
    !
  end if
!
! exit
800 continue
  return
!  
! error handlings
!
end subroutine CS_SHIFT_ATOMS


!**********************************************************************!
!
! CS_DICE_PARTIALOCC
!
! subroutine, randomized the occupancy of partially occupied sites
!             taking into account possible pair site occupations
!             the routine can handle up to 10 partial occupations per site
!
! INPUT:
!   real*4 :: dmax      = max. distance of pair sites (nm)
!   
! IN/OUTPUT:
!   integer*4 :: nerr   = error code (0 = success)
!
subroutine CS_DICE_PARTIALOCC(dmax, nerr)

  implicit none
  
! parameters
  integer*4, parameter :: npairmax = 10
  
! interface variables
  real*4, intent(in) :: dmax
  integer*4, intent(inout) :: nerr
  
! local vaiables
  integer*4 :: ierr ! local error code
  integer*4 :: nalloc ! allocation status
  integer*4 :: i, j, k, l, m, n ! some iterators and numbers
  integer*4, allocatable, dimension(:,:) :: pairs
  real*4 :: mocc, scell(3), rloc(3), dfra(3), dvec(3), dist ! some temporary values
  real*4 :: pmax, pran, pcum ! random number things
! temporaray atomic data 
!     atomic type
  character(len=CS_atsl), allocatable, dimension(:) :: l_attype
  integer*4, allocatable, dimension(:) :: l_atnum
  integer*4, allocatable, dimension(:) :: l_atuse
  real*4, allocatable, dimension(:) :: l_atcrg
  real*4, allocatable, dimension(:,:) :: l_atpos
  real*4, allocatable, dimension(:) :: l_atocc
  real*4, allocatable, dimension(:) :: l_atdwf
! external references
  external :: InitRand
  real*4, external :: UniRand
  
! initialization
  call InitRand()
  ierr = 0
  nerr = 0
  nalloc = 0
  scell(1) = CS_scsx
  scell(2) = CS_scsy
  scell(3) = CS_scsz
  !
  if (.not.CS_cellmem_allocated .or. CS_numat<=0) goto 800 ! nothing to do
  !
  ! There are atoms in the current super-cell.
  n = CS_numat
  ! - allocate the temporary memory for the creation of the new structure
  allocate(l_attype(n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_attype = ""
  allocate(l_atnum(n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_atnum = 0
  allocate(l_atcrg(n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_atcrg = 0.0
  allocate(l_atpos(3,n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_atpos = 0.0
  allocate(l_atocc(n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_atocc = 0.0
  allocate(l_atdwf(n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_atdwf = 0.0
  allocate(l_atuse(n), stat=nalloc)
  if (nalloc/=0) goto 804
  l_atuse = 1 ! use all
  !
  ! Copy the current data
  do i=1, n
    l_attype(i) = CS_attype(i)
    l_atnum(i) = CS_atnum(i)
    l_atcrg(i) = CS_atcrg(i)
    l_atpos(1:3,i) = CS_atpos(1:3,i)
    l_atocc(i) = CS_atocc(i)
    l_atdwf(i) = CS_atdwf(i)
  end do
  !
  ! Allocate an index field for the partial occupation calculation
  ! "pairs(0:m+1,1:n)"
  ! n = number of atoms in the input cell
  ! m = npairmax = 5
  ! 0   = flag, ==0 -> not changed, ==1 -> to be randomized
  ! 1   = count of pair candidates 
  ! 2:m+1 = up to five occupation candidates
  !
  allocate(pairs(0:npairmax+1,n),stat=nalloc)
  if (nalloc/=0) goto 804
  pairs = 0
  !
  ! pick up candidates
  do i=1, n
    if (l_atocc(i)<1.0) then ! candidate for dicing
      pairs(0,i) = 1
    end if
  end do
  !
  ! pick up candidate pairs
  do i=1, n
    if (pairs(0,i) == 1) then ! candidate for dicing
    rloc(1:3) = l_atpos(1:3,i)
      do j=1, n
        if (j==i) cycle ! no self-pairing snail
        if (pairs(0,j) == 1) then ! candidate for dicing
          ! check distance
          ! - get distance vector
          !   - min distance accross periodic boundaries
          !   - d' = modulo( abs(d) + 0.5, 1.0 ) - 0.5
          !   - 0.8 -> 0.3
          !   - 0.3 -> 0.3
          dfra(1:3) = rloc(1:3)-l_atpos(1:3,j)
          dfra(1:3) = modulo( abs(dfra) + 0.5, 1.0 ) - 0.5
          dvec(1:3) = dfra*scell(1:3)
          ! - get distance
          dist = sqrt(sum(dvec*dvec))
          ! - compare to dmax
          if (dist<dmax) then
            k = pairs(1,i) ! get current number of pair candidates
            if (k<npairmax) then ! max. number of pairs not exceeded
              ! add new pair candidate
              pairs(1,i) = k+1
              pairs(k+2,i) = j
            else
              ! replace the current pair candidate with lower and lowest occupancy number
              m = 0
              mocc = 1.0
              do l=1, k
                if (l_atocc(pairs(l+1,i))<mocc) then
                  mocc = l_atocc(pairs(l+1,i))
                  m = l
                end if
              end do
              if (m>0) then ! found smaller and smallest pair candidate to be replaced
                pairs(m+1,i) = j
              end if ! m>0
            end if ! k<npairmax
          end if ! dist<dmax
        end if ! pairs(0,j)==1
      end do
    end if
  end do
  !
  ! randomize the partial occupancies
  do i=1,n
    if (pairs(0,i)==0) cycle ! skip this guy, no randomization here
    k = pairs(1,i) ! get number of other candidates
    pmax = 1.0 ! max probability of occupation
    pran = UniRand() ! get the random choice
    if (k==0) then ! randomize against vacancy
      if (pran<=l_atocc(i)) then ! occupy
        l_atocc(i) = 1.0 ! full atom
      else
        l_atocc(i) = 0.0 ! no atom
        l_atuse(i) = 0   ! clear from list later
      end if
      pairs(0,i) = 0 ! no longer indicate as partial
    else ! randomize against other candidates
      m = 0 ! final candidate chosen ! preset as vacancy
      pcum = l_atocc(i) ! preset cumulated probability of occupation
      j = i ! preset current test atom with "this" atom (i)
      l = 0 ! preset alternative loop index to 0 (none of the alternatives)
      do ! loop through candidates, starting with this atom (i), followed by the alternatives
        if (pran<=pcum) then ! found the candidate which will occupy the site
          m = j ! remember it in m
          exit ! stop the search loop
        end if
        ! proceed to next candidate
        l = l + 1 ! increase list index
        if (l>k) exit ! last candidate checked, get out of here
        j = pairs(l+1,i) ! this is the next candidate index
        pcum = pcum + l_atocc(j) ! update cumulated probability
      end do
      ! occupy now with m, clear all alternative lists then
      if (m==i) then ! occupy with this atom (i)
        l_atocc(i) = 1.0
        l_atuse(i) = 1
      else
        l_atocc(i) = 0.0
        l_atuse(i) = 0
      end if
      do l=1,k ! loop through the alternative list and set the occupancie according to the made choice
        j = pairs(l+1,i) ! get the candidate index
        if (m==j) then ! occupy with this atom
          l_atocc(j) = 1.0
          l_atuse(j) = 1
        else ! de-occupy this atom
          l_atocc(j) = 0.0
          l_atuse(j) = 0
        end if
        pairs(:,j) = 0 ! clear atom j's pair list, ensures that there is no other check done with atom j
      end do
      pairs(:,i) = 0 ! clear atom i's pair list in any case
      ! if neither P#1 nor P#2 applied, a vacancy is generated (atoms are removed from this site)
      ! more than npairsmax partial occupations per site will cause errors here
    end if
    !
  end do
  !
  ! re-set atom list in the module arrays
  ! copy the content of the new cell to the module arrays
  ! - prepare the module arrays
  call CS_DEALLOC_CELLMEM(ierr)
  if (ierr/=0) goto 805
  m = SUM(l_atuse)
  call CS_ALLOC_CELLMEM(m, ierr) ! sets CS_numat = m
  if (ierr/=0) goto 804
  !
  if (m>0) then
    ! - copy
    j = 0
    do i=1, n
      if (l_atuse(i)==1) then
        j = j + 1
        CS_attype(j) = l_attype(i)
        CS_atnum(j)  = l_atnum(i)
        CS_atcrg(j)  = l_atcrg(i)
        CS_atpos(1:3,j) = l_atpos(1:3,i)
        CS_atocc(j)  = l_atocc(i)
        CS_atdwf(j)  = l_atdwf(i)
      end if
    end do
    !
  end if
  !
  !
  ! clean up the mess
  deallocate(l_attype, stat=nalloc)
  deallocate(l_atnum, stat=nalloc)
  deallocate(l_atcrg, stat=nalloc)
  deallocate(l_atpos, stat=nalloc)
  deallocate(l_atocc, stat=nalloc)
  deallocate(l_atdwf, stat=nalloc)
  deallocate(pairs, stat=nalloc)
  !
!
! exit
800 continue
  return
!  
! error handlings
804 continue
  nerr = 4
  call CS_ERROR("Memory allocation failed.")
  return
805 continue
  nerr = 5
  call CS_ERROR("Memory de-allocation failed.")
  return
end subroutine CS_DICE_PARTIALOCC



!**********************************************************************!
!
! CS_SUGGEST_NSLCEQUI
!
! subroutine, suggests a number "nslc" of slices for partitioning the
!             current super cell equidistantly along the c axis
!             the decision is made based on a minimum slice thickness
!             "dmin" and the distribution of atoms in the cell
!
! INPUT:
!   real*4 :: dmin      = min. distance of slices (nm)
!   
! IN/OUTPUT:
!   integer*4 :: nslc   = returned number of slices
!   integer*4 :: nerr   = error code (0 = success)
!
subroutine CS_SUGGEST_NSLCEQUI(dmin, nslc, nerr)

  implicit none

! routine parameters
  integer*4, parameter :: nft_max = 8192 ! max. fft size
  real*4, parameter :: dmin_default = 0.02 ! 0.2 A (very thin slices)
  
! interface variables
  real*4, intent(in) :: dmin
  integer*4, intent(inout) :: nslc, nerr
  
! local vaiables
  integer*4 :: ierr ! local error code
  integer*4 :: nalloc ! allocation status
  integer*4 :: i, j, m, n ! some iterators and numbers
  integer*4 :: nhis ! histogram bins
  real*4 :: zpos, pcur, rfac, rpos, vcur
  real*4 :: dmini ! internall minimum size used
  real*4 :: shis ! histogram sampling
  complex*8, allocatable :: ahis(:) ! histogram array
  external :: ODDCC8192
  
  
! initialization
  ierr = 0
  nerr = 0
  nalloc = 0
  nslc = 1
  dmini = dmin
  !
  if (CS_scsz<=0.0) goto 801 ! error cell not defined
  if (.not.CS_cellmem_allocated .or. CS_numat<=0) goto 800 ! nothing to do
  if (dmini<=0.0) dmini = dmin_default ! invalid minimum slice thickness, run with default parameter
  !
  ! There are atoms in the current super-cell.
  n = CS_numat
  ! Preset the histogram sampling.
  shis = 0.5*dmini
  nhis = ceiling(CS_scsz / shis)
  if (nhis<2) goto 802
  m = ceiling(real(nhis)*0.5)
  if (nhis>nft_max) goto 803
  ! Prepare the histogram array
  allocate(ahis(nft_max),stat=nalloc)
  if (nalloc/=0) goto 804
  ahis = cmplx(0.0,0.0)
  !
  ! Reset the histogram sampling rate to match the chosen binning.
  shis = CS_scsz / real(nhis)
  rfac = -1.0/(dmini*dmini)
  ! Setup the histogram data by adding a Gaussian for each atom
  do j=1, n ! loop over all atoms
    zpos = CS_atpos(3,j)*CS_scsz
    do i=1, nhis ! loop through each of the histogram pixels
      pcur = real(i-1)*shis
      rpos = pcur-zpos
      vcur = exp( rfac*rpos*rpos )
      ahis(i) = ahis(i) + cmplx( vcur, 0.0 ) ! add to the histogram
    end do
  end do
  ! Fourier transform the histogram
  call ODDCC8192(ahis,nhis,'for')
  ! Analyse the Fourier-Transform, find the component of highest periodicity
  nslc = 1 ! init the max. finder (this is the frequency number = i-1 )
  pcur = cabs(ahis(1))
  pcur = cabs(ahis(2)) ! ... also
  do i=3, m ! start searching from k=2 (i=3), k=1 is preset as default solution
    rpos = cabs(ahis(i))
    if (rpos>pcur) then ! found better periodicity
      nslc = i - 1
      pcur = rpos
    end if
  end do
! done.
! exit
800 continue
  if (allocated(ahis)) deallocate(ahis,stat=nalloc)
  return
!  
! error handlings
801 continue
  nerr = 1
  call CS_ERROR("Invalid size of super-cell (zero length).")
  goto 800
802 continue
  nerr = 2
  call CS_ERROR("Too few number of samples (<2). Aborting.")
  goto 800
803 continue
  nerr = 3
  call CS_ERROR("Too many samples (>8192). Aborting.")
  goto 800
804 continue
  nerr = 4
  call CS_ERROR("Memory allocation failed.")
  goto 800
!
  return
end subroutine CS_SUGGEST_NSLCEQUI



END MODULE CellSlicer
!**********************************************************************!
!**********************************************************************!