!**********************************************************************!
!**********************************************************************!
!
! FILE: "simepw.f90"
!
! AUTHOR: Dr. J. Barthel
!         Ernst Ruska-Centre
!         Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
!           and
!         RWTH Aachen University, 52074 Aachen, Germany
!
! PURPOSE: Implementations of program CELSLC
!
! VERSION: 0.41b (20150921)
!
!**********************************************************************!
!**********************************************************************!


!**********************************************************************!
!
!  PROGRAM: CELSLC
!
!**********************************************************************!

program celslc

  use CellSlicer
  use EMSdata
  use FSCATAB
  use m3dpot

  implicit none
  
  
  ! ****************************************************************** !
  !
  ! *** DECLARATIONS *** !
  !
  include "global.fi"
  
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
  CS_doconsolemsg = 1
  call InitRand()
  call CS_INIT()
  call EMS_INIT()
  call Introduce()
  nerr = 0
  nfl = 1
  ndwf = 0
  nabs = 0
  nabf = 0
  nzd = 0
  nfx = 0
  nfe = 0
  i = 0
  j = 0
  k = 0
  nat = 0
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** INPUT *** !
  !
  call PostMessage("Parsing command line.")
  call ParseCommandLine()
  if (nerr/=0) call CriticalError("Failed to determine input parameters.")
  !
  ! determine output file name number digits
  write(unit=stmp1,fmt=*) nz
  ndigsl = max(3, len_trim(adjustl(stmp1)))
  write(unit=stmp1,fmt=*) nv
  ndigvr = max(3, len_trim(adjustl(stmp1)))
  !
  ! handle use of external scattering factors
  if (nfx>0) then
    call FST_INIT()
    call FST_LoadFX(trim(sfxfile), nerr)
    if (nerr/=0) then
      call CriticalError("Failed to read scattering factors from file ["//trim(sfxfile)//"].")
    end if
    call FST_SaveFEP(trim(sfxfile)//".fep", nerr, 1)
    CS_useextsca = 1
  end if
  if (nfe>0) then
    call FST_INIT()
    call FST_LoadFEP(trim(sfefile), nerr)
    if (nerr/=0) then
      call CriticalError("Failed to read scattering factors from file ["//trim(sfefile)//"].")
    end if
    CS_useextsca = 1
  end if
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** CALCULATION *** !
  ! 
  select case (nfin) ! distinguish structure data input forms
  
  !--------------------------------------------------------------------!
  
  case (0) ! **** DEFAULT VERSION :: CEL FILE **** !
  
    call PostMessage("Loading data from super-cell file ["// &
     &               trim(scellfile)//"].")
    call CS_LOAD_EMSCELL(trim(scellfile), nerr)
    if (nerr/=0) then
      call CriticalError("Failed to read super-cell data from file [" &
     &               //trim(scellfile)//"].")
    end if
    !
    if (ndwf==1 .and. buni==1 .and. CS_numat>0) then
      write(unit=smsg,fmt=*) buniv
      call PostMessage("Replacing individual B_ISO values by "// &
     &    trim(adjustl(smsg))//" nm^2 for all atoms.")
      CS_atdwf = buniv
    end if
    
    ! determine sampling rates
    sdx = CS_scsx / real(nx)
    sdy = CS_scsy / real(ny)
    sdz = CS_scsz / real(nz)
  
    ! post cell info to console
    call PostCellInfo()
    
    ! check single slice calculation index
    ssc = min(ssc,nz) ! limit slice index to number of slices
    !
    
    !
    ! --- make & save slices ---<<<
    !
    if (0==n3dp) then
      call CEL2SLC()
    else
      call CEL2POT3D2SLC()
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
!    ! dump the potential back to HD
!    call savedata1( "pot3d.dat", 2*M3D_n1*M3D_n2*M3D_n3, M3D_pot, nerr)
    ! post cell info to console
    call M3D_PostCellInfo()
    
    ! orthogonalize the cell if required
    call M3D_Othogonalize(fft_dmax, nerr)
    !
    ! apply apterture and dwf
    call M3D_DiffractionFilter(1, buni, buniv)
    ! 
    ! dump the potential back to HD
    call savedata1( "pot3d-2.dat", 2*M3D_n1*M3D_n2*M3D_n3, M3D_pot, nerr)
    ! post cell info to console
    call M3D_PostCellInfo()
    
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
    call POT3D2SLC()
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
  call Outroduce()
  !
  ! ****************************************************************** !
  
  
  ! ****************************************************************** !
  !
  ! *** EXIT *** !
  !
  call EMS_UNINIT()
  call CS_UNINIT()

end program celslc

!**********************************************************************!
!**********************************************************************!