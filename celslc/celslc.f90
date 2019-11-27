!**********************************************************************!
!**********************************************************************!
!
! FILE: "celslc.f90"
!
! AUTHOR: Dr. J. Barthel
!         Ernst Ruska-Centre
!         Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
!           and
!         RWTH Aachen University, 52074 Aachen, Germany
!
! PURPOSE: Implementations of program CELSLC
!
! VERSION: 1.0.1 (20191127)
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
!
! When using this program, please consider to cite the following
! article as an appreciation of the many hours of work put into
! this project:
!
! J. Barthel, "Dr. Probe - A software for high-resolution STEM
! image simulation", Ultramicroscopy 193 (2018) 1-11.
! doi: 10.1016/j.ultramic.2018.06.003
!
! Read 'celslc howto.txt' for detailed explanations about the program.
!
!----------------------------------------------------------------------



!**********************************************************************!
!
!  PROGRAM: CELSLC
!
!**********************************************************************!

program celslc

  use celslcprm
  use CellSlicer
  use EMSdata
  use FSCATAB
  use m3dpot

  implicit none
  
  
  ! ****************************************************************** !
  !
  ! *** DECLARATIONS *** !
  !
  integer*4 :: i, j, k, nat
  character(len=1024) :: smsg, sfil, stmp1, stmp2
  real*4 :: fz0, fz1
  complex*8, allocatable, dimension(:,:,:) :: slcdat
  
  external :: InitRand
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** INITIALIZATION
  !
  nerr = 0
  nverbose = 1
  nfl = 1
  ndwf = 0
  nabs = 0
  nabf = 0
  nfx = 0
  nfe = 0
  i = 0
  j = 0
  k = 0
  nat = 0
  call InitRand()
!  call CS_INIT()
  call EMS_INIT()
  call ParseCommandLine()
  CS_doconsolemsg = nverbose
  FST_domsg = nverbose
  if (nerr/=0) call CriticalError("Failed to determine input parameters.")
  call Introduce()
  call CheckCommandLine()
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** INPUT *** !
  !
  call csprm_initcl()
  !
!  ! determine max size of FFT
!  k = max(nx, ny)
!  CS_FFT_BOUND = max(CS_FFT_BOUND_MIN,2**CEILING(LOG(real(k))/LOG(2.0)))
!  if (CS_FFT_BOUND>CS_FFT_BOUND_MAX) call CriticalError("FFT plan exceeds maximum size (8192).")
  call CS_INIT()
  !
  ! determine output file name number digits
  write(unit=stmp1,fmt=*) nz
  ndigsl = max(3, len_trim(adjustl(stmp1)))
  write(unit=stmp1,fmt=*) nv
  ndigvr = max(3, len_trim(adjustl(stmp1)))
  !
  ! handle use of external scattering factors
  CS_useextsca = 0
  if (nfx>0) then
    call FST_INIT()
    call FST_LoadFX(trim(sfxfile), nerr)
    if (nerr/=0) then
      call CriticalError("Failed to read x-ray form factors from file ["//trim(sfxfile)//"].")
    else
      call FST_SaveFEP(trim(sfxfile)//".fep", nerr, 1)
      CS_useextsca = 1
    end if
  end if
  if (nfe>0) then
    call FST_INIT()
    call FST_LoadFEP(trim(sfefile), nerr)
    if (nerr/=0) then
      call CriticalError("Failed to read el. form factors from file ["//trim(sfefile)//"].")
    else
      CS_useextsca = 1
    end if
  end if
  if (CS_useextsca/=1 .and. nfin<10) then ! internal form-factor table switching
    if (nsca==2) CS_scaf_table = 2
    if (CS_scaf_table==1) call PostMessage("Using atomic form factors of Weickenmeier & Kohl.")
    if (CS_scaf_table==2) call PostMessage("Using atomic form factors of Waasmaier & Kirfel.")
  end if
  !
  ! handle external potentials
  if (nextpot>0) then
    call PrepExtPotentials(nerr)
    if (nerr/=0) then
      call CriticalError("Failed to prepare additional external potentials.")
    end if
  end if
  !
  call PostRuntime("program initialized")
  !
  if (nffdec==1) then ! output ff decay data
    call CS_WRITE_FFDEC(vffdec, nerr)
    if (nerr==0) then
      call PostMessage("Form-factor decay data written to file 'ffdec.txt'.")
    else
      call PostMessage("Failed to output form-factor decay data.")
    end if
  end if
  !
  if (nf2dec==1) then ! output ff decay data
    call CS_WRITE_F2DEC(vf2dec, ht, nerr)
    if (nerr==0) then
      call PostMessage("Analysis for loss of scattering power written to file 'f2dec.txt'.")
    else
      call PostMessage("Failed to output analysis for loss of scattering power.")
    end if
  end if
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** CALCULATION *** !
  ! 
  select case (nfin) ! distinguish structure data input forms
  
  !--------------------------------------------------------------------!
  
  case (0:1) ! **** DEFAULT VERSION :: CEL/CIF FILE **** !
  
    call PostMessage("Loading data from super-cell file ["// &
     &               trim(scellfile)//"].")
    
    if (nfin==0) then ! CEL
      call CS_LOAD_EMSCELL(trim(scellfile), nerr)
    elseif (nfin==1) then ! CIF
      call CS_LOAD_CIFCELL(trim(scellfile), nerr)
    end if
    if (nerr/=0) then
      call CriticalError("Failed to read super-cell data from file [" &
     &               //trim(scellfile)//"].")
    end if
    !
    call PostRuntime("loaded structure file")
    !
    ! Post loading operations
    !
    if (buni==1 .and. CS_numat>0) then ! universal Biso values
      write(unit=smsg,fmt=*) buniv
      call PostMessage("Replacing individual B_ISO values by "// &
     &    trim(adjustl(smsg))//" nm^2 for all atoms.")
      CS_atdwf = buniv
    end if
    !
    if (blk==1) then ! re-orientations
      !
      call PostMessage("Re-orienting the super-cell to new projection axes.")
      call CS_ORIENT_CELL(bloh,blok,blol, blyh,blyk,blyl, blsa,blsb,blsc, nerr)
      if (nerr/=0) then
        ! orientation error
      end if
      
    end if
    !
    !
    ! check the current super cell
    call CheckCell()
    
    if (CS_scorth_check<0) then ! need to orthogonalize the super-cell
      !
      call PostMessage("Try calling CELSLC again with the -prj option.")
      call CriticalError("I refuse to run with a non-orthogonal super-cell.")
      !
    end if
    !
    !
    if (ntla==1) then ! translate atoms
      !
      call PostMessage("Shifting the atoms along the cell axes.")
      call CS_SHIFT_ATOMS(tlax,tlay,tlaz, nerr)
      !
    end if
    !
    ! store / backup used super-cell
    if (nfin==0) then ! CEL
      call CS_SAVE_EMSCELL("CELSLC_used.cel",nerr)
    elseif (nfin==1) then ! CIF
      call CS_SAVE_CIFCELL("CELSLC_used.cif",nerr)
    end if
    
    ! post cell info to console (also calculates sampling rates sdx, sdy, sdz)
    call PostCellInfo()
    
    !
    call PostRuntime("prepared structure data")
    
    !
    ! --- make & save slices ---<<<
    !
    if (0==n3dp) then
      call CEL2SLC()
      call PostRuntime("generated phase gratings from projected potentials")
    else
      call CEL2POT3D2SLC()
      call PostRuntime("generated phase gratings from 3D potential")
    end if
    !
    !
  
    
  !--------------------------------------------------------------------!
    
  case (10:19) ! **** 3D POTENTIAL VERSION :: ASC FILE ******** !
               !      The last digit tells this program from
               !      which column the potential values should
               !      read in.
    
    call PostMessage("Loading 3D potential parameters from file ["// &
       & trim(scellfile)//"].")
    call M3D_readasctxt(trim(scellfile), nfin-10, nerr)
    if (nerr/=0) then
      call CriticalError("Failed to read 3D potential parameters "// &
         & "from file ["//trim(scellfile)//"].")
    end if
    call PostRuntime("loaded external 3D potential")
!    ! dump the potential back to HD
!    call savedata1( "pot3d.dat", 2*M3D_n1*M3D_n2*M3D_n3, M3D_pot, nerr)
    ! post cell info to console
    call M3D_PostCellInfo()
    
    ! orthogonalize the cell if required
    call M3D_Othogonalize(fft_dmax, nerr)
    !
    ! apply apterture and dwf
    call M3D_DiffractionFilter(1, buni*ndwf, buniv)
    ! 
    ! dump the potential back to HD
    call savedata1( "pot3d-2.dat", 2*M3D_n1*M3D_n2*M3D_n3, M3D_pot, nerr)
    ! post cell info to console
    call M3D_PostCellInfo()
    call PostRuntime("prepared potential data")
    
    !
    nx = M3D_n1
    ny = M3D_n2
    if (nz<=0) nz = M3D_n3 ! use potential z-sampling when no nz option is given
    if (nz>M3D_n3) then
      call PostWarning("Specified number of slices is larger than "// &
     &     "the z-sampling of the 3d potential.")
      call PostWarning("Number of slices is automatically reduced to"//&
     &     "match the z-sampling of the 3d potential.")
      nz = M3D_n3
    end if
    !
    call M3D_CellDimension(CS_scsx, CS_scsy, CS_scsz )
    !
    ! determine sampling rates
    sdx = CS_scsx / real(nx)
    sdy = CS_scsy / real(ny)
    sdz = CS_scsz / real(nz)
    !
    ! check single slice calculation index
    ssc = min(ssc,nz) ! limit slice index to number of slices
    !
    
    !
    ! Potentials from external sources need a relativistic correction.
    M3D_pot = M3D_pot * (1.0 + ht / 511.9989)
    
    !
    ! --- make & save slices from the 3d potential ---<<<
    !
    ! prepare slicing
    call CS_ALLOC_CELLMEM(0,nerr)
    call CS_ALLOC_SLICEMEM(nz, nerr)
    call CS_SETSLICE_EQUIDIST(nz,nrev, nerr)
    !
    call POT3D2SLC()
    call PostRuntime("generated phase gratings from 3D potential")
    !
    !
    
  !--------------------------------------------------------------------!
   
  case default ! *** UNKNOWN VERSION :: ERROR ***************** !
    
    if (nerr/=0) then
      call CriticalError("Unknown format specified for input file ["// &
     &     trim(scellfile)//"].")
    end if
  
  !--------------------------------------------------------------------!
  
  end select ! case (nfin)
  !
  ! ****************************************************************** !
  
  
  
  ! ****************************************************************** !
  !
  ! *** FINISH *** !
  !
  call PostMessage("Finished.")
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** EXIT *** !
  !
  call PostRuntime("terminating")
  call EMS_UNINIT()
  call CS_UNINIT()

end program celslc

!**********************************************************************!
!**********************************************************************!