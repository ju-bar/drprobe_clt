!**********************************************************************!
!**********************************************************************!
!
! FILE: "3dpot.f90"
!
! AUTHOR: Dr. J. Barthel, ju.barthel@fz-juelich.de
!         Ernst Ruska-Centre
!         Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
!         RWTH Aachen University, 52074 Aachen, Germany
!
! PURPOSE: Implementations for input of 3d potentials from files
! 
! MODULE: 3dpot
!
! VERSION: 1.1, J.B., 09.06.2015
!
! LINK:   "binio2.f90", "FFTs.f"
!         link a file which implements the routines <UPPERCASE> and 
!         <LOWERCASE>, e.g. "simsubs.f90"
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

MODULE m3dpot

  implicit none
  
  SAVE
  
  !
  ! *** MODULE ROUTINES *** !
  !
  public :: M3D_readasctxt
  !
  public :: M3D_CrossProduct
  public :: M3D_CellDimension
  public :: M3D_CellVolume
  public :: M3D_CellAngles
  public :: M3D_PostCellInfo
  !
  private :: M3D_FFT1D
  private :: M3D_FFT2D
  private :: M3D_FFT3D
  !
  public :: M3D_OrthoGamma
  public :: M3D_OrthoBeta
  public :: M3D_OrthoAlpha
  public :: M3D_Othogonalize
  !
  public :: M3D_DiffractionFilter
  public :: M3D_ApplyDiffApAuto
  public :: M3D_ApplyDWF
  !
  public :: M3D_getslice_pgr
  !
  private :: M3D_Message
  private :: M3D_ErrorMessage
  
  !
  ! *** MODULE PARAMETERS *** !
  !
  integer*4, public, parameter :: M3D_ll = 2048 ! default string/line length
  integer*4, public, parameter :: M3D_stdout = 6 ! standard output unit
  integer*4, public, parameter :: M3D_FFTBOUND = 8192 ! max. size of implemented FFT routines
  real*4, public, parameter :: M3D_orththr = 0.001 ! orthogonality threshold (rad)
  real*8, public, parameter :: M3D_diophthr = 1.0D-04 ! diophantine threshold (rel. mismatch)
  
  !
  ! *** MODULE VARIABLES *** !
  !
  integer*4, public :: M3D_intype ! input type / set on call of respective input routine
  DATA M3D_intype /0/             ! 0 = default = internal
                                  ! 1 = asc or asc.txt file
                                  
  integer*4, public :: M3D_out    ! output unit
  DATA M3D_out /M3D_stdout/       ! default = console output unit
                                  
  integer*4, public :: M3D_n1     ! 1st dimension sampling
  DATA M3D_n1 /0/
  integer*4, public :: M3D_n2     ! 2nd dimension sampling
  DATA M3D_n2 /0/
  integer*4, public :: M3D_n3     ! 3rd dimension sampling
  DATA M3D_n3 /0/
  
  real*8, public :: M3D_b1(3)     ! 1st dimension physical basis vector
  DATA M3D_b1 /0.0,0.0,0.0/
  real*8, public :: M3D_b2(3)     ! 2nd dimension physical basis vector
  DATA M3D_b2 /0.0,0.0,0.0/
  real*8, public :: M3D_b3(3)     ! 3rd dimension physical basis vector
  DATA M3D_b3 /0.0,0.0,0.0/
  
  complex*8, public, allocatable :: M3D_pot(:,:,:) ! potential, allocated on call of input routine
  
  integer*4, public :: M3D_backup_slcpot ! flag: set to 1 to switch slice potential backup saving
  DATA M3D_backup_slcpot /0/
  complex*8, public, allocatable :: M3D_slcpot(:,:) ! memory to store the current slice potential
  
  CONTAINS

!**********************************************************************!
!
subroutine M3D_readasctxt(sfile, ncol, nerr)
!
! Purpose: Reads the data from the file specified by file name "sfile".
!          Extracts the cell parameters (size and sampling) and the
!          potential data.
!
! Parameters:
!          character(len=*) :: sfile = file name
!          integer*4 :: ncol = potential input column (read potential from item 4+ncol of each line)
!          integer*4 :: nerr = error code
!
! Remarks: Expects a file in text form consisting of 3 header lines
!          followed by data lines.
!          Header line content:
!          1) cell title definition
!          2) variable and unit definition
!          3) cell sampling definition
!          Here, ...
!          - we ignore the title line.
!          - we assume variables "X"  "Y"  "Z" (all three in Angstrom)
!            and "DENSITY or POTENTIAL" (atomic units)
!          - Sampling information of the form
!            ZONE I=    36 J=    36 K=   144 F=POINT
!          - Data line 4-tuple format
!         --> 2.461500000000000E+00 6.316211944934507E-01 1.440000000000000E+00-5.471180606696066E-01<--
!          The sampling-rate is assumed to be equidistant and
!          extracted from the steps in the 3 coordinates.
!

  implicit none
  
  character(len=*), intent(in) :: sfile
  integer*4 :: ncol
  integer*4, intent(inout) :: nerr
  
  integer*4 :: lun, ioerr, nalloc
  integer*8 :: ndata
  integer*4 :: i, j, k, di, dj, dk, nsd(3)
  logical :: fex
  character(len=M3D_ll) :: sline, sform, sinfo
  real*8 :: r(13), r0(3), dr(3), vp, d1, d2
  real*8 :: sd(3,3), ssd(3,3), rserr(3), vltmp
  integer*4, external :: getfreelun
  
1001 FORMAT(A8,A) ! format of the title line
1002 FORMAT(A12,A) ! incomplete form of the variable line
1003 FORMAT(A7,I6,A3,I6,A3,I6,A3,A) ! format of the sampling info line

  nerr = 0
  ioerr = 0
  nalloc = 0
  lun = getfreelun()
  if (lun<0) goto 101
  inquire( file=trim(sfile), EXIST=fex )
  if (.not.fex) goto 102
  open( file=trim(sfile), unit=lun, iostat=ioerr, action='READ', &
     &  status='OLD' )
  if (ioerr/=0) goto 103
  
  ! read & forget the first 2 header lines
  read(unit=lun, fmt=1001, iostat=ioerr) sinfo, sline
  if (ioerr/=0) goto 104
  read(unit=lun, fmt=1002, iostat=ioerr) sinfo, sline
  if (ioerr/=0) goto 104
  ! read & interprete the 3rd header line
  read(unit=lun, fmt=1003, iostat=ioerr) sinfo(1:7), M3D_n1, &
     & sinfo(8:10), M3D_n2, sinfo(11:13), M3D_n3, sinfo(14:16), sform
  if (ioerr/=0) goto 104
  
  if ( M3D_n1<=0 .or. M3D_n2<=0 .or. M3D_n3<=0 ) goto 105
  if ( trim(sform)/="POINT" ) goto 106
  !
  ! allocate the potential memory
  if (allocated(M3D_pot)) deallocate(M3D_pot, stat=nalloc)
  if (nalloc/=0) goto 108
  allocate(M3D_pot(M3D_n1,M3D_n2,M3D_n3), stat=nalloc)
  if (nalloc/=0) goto 109
  !
  ! loop max. over all points of the data
  ndata = 0
  nsd = 0
  sd = 0.0D+00
  ssd = 0.0D+00
  !
  write(unit=sinfo,fmt='(A,I1)') &
     &  "(M3D_readasctxt): " // &
     &  "Reading potential from column #", 4+ncol
  call M3D_Message(trim(sinfo))
  !
  do k=1, M3D_n3
    dk = k-1
    do j=1, M3D_n2
      dj = j-1
      do i=1, M3D_n1
        di = i-1
        ! read the next line of data
        read(unit=lun, fmt=*, iostat=ioerr) r(1:4+ncol)
        if (ioerr/=0) goto 107
        ndata = ndata + 1 ! number of data read so far
        if (ndata==1) then ! 1st point (1,1,1) used as coordinate reference
          r0 = r(1:3)
        else
          dr = r(1:3) - r0
          if (dj==0.and.dk==0) then ! 1st basis vector accumulation
            sd(:,1) = sd(:,1) + dr/dble(di)
            ssd(:,1) = ssd(:,1) + dr*dr/dble(di)/dble(di)
            nsd(1) = nsd(1) + 1
          end if
          if (di==0.and.dk==0) then ! 2nd basis vector accumulation
            sd(:,2) = sd(:,2) + dr/dble(dj)
            ssd(:,2) = ssd(:,2) + dr*dr/dble(dj)/dble(dj)
            nsd(2) = nsd(2) + 1
          end if
          if (di==0.and.dj==0) then ! 3rd basis vector accumulation
            sd(:,3) = sd(:,3) + dr/dble(dk)
            ssd(:,3) = ssd(:,3) + dr*dr/dble(dk)/dble(dk)
            nsd(3) = nsd(3) + 1
          end if
        end if
        ! store potential value
        ! inverse sign -> scattering factor
        M3D_pot(i,j,k) = cmplx( -real( r(4+ncol), kind=4 ), 0.0)
      end do
    end do
  end do
  !
  ! analyse the sampling rate basis
  rserr = 0.0D+00
  do i=1, 3
    ! - normalize and calculate mean values
    if (nsd(i)>0) then
      sd(:,i) = sd(:,i) / dble(nsd(i))
    end if
    ! - normalize and calculate variances & relative errors
    vltmp = SUM( sd(:,i)*sd(:,i) )
    if (nsd(i)>1.and.vltmp>1.0D-100) then
      ssd(:,i) = ( ssd(:,i) / dble(nsd(i)) - sd(:,i)*sd(:,i) ) &
     &           * dble(nsd(i)) / dble( nsd(i)-1 )
      rserr(i) = SUM( ssd(:,i) ) / vltmp
      if (rserr(i)>1.0D-03) then
        write(unit=sinfo,fmt='(A,I1)') &
     &  "Warning (M3D_readasctxt): " // &
     &  "Non-equidistant sampling found in dimension ", i
        call M3D_Message(trim(sinfo))
      end if
    end if
  end do
  M3D_b1(:) = sd(:,1)*0.1 ! to module and to nm
  M3D_b2(:) = sd(:,2)*0.1 ! to module and to nm
  M3D_b3(:) = sd(:,3)*0.1 ! to module and to nm
  
99 close(unit=lun, iostat=ioerr)

  M3D_intype = 1 ! mark input format type

  return
  
! --- Error handling
!
! ------ Error #1 : No free logical file unit. Unable to open file.
101 nerr = 1
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "No free logical file unit. Unable to open file.", nerr )
  return
! ------ Error #2 : File doesn't exist.
102 nerr = 2
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "File doesn't exist.", nerr )
  return
! ------ Error #3 : Failed to open file for reading.
103 nerr = 3
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Failed to open file for reading.", nerr )
  return
! ------ Error #4 : Failed to read header data from file.
104 nerr = 4
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Failed to read header data from file.", nerr )
  return
! ------ Error #5 : Invalid cell sampling
105 nerr = 5
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Invalid cell sampling.", nerr )
  return
! ------ Error #6 : Unsupported data format
106 nerr = 6
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Unsupported data format.", nerr )
  return
! ------ Error #7 : Failed to read data from file.
107 nerr = 7
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Failed to read data from file.", nerr )
  return
! ------ Error #8 : Failed to deallocate an allocated array.
108 nerr = 8
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Failed to deallocate an allocated array.", nerr )
  return
! ------ Error #9 : Failed to allocate array memory.
109 nerr = 9
  call M3D_ErrorMessage("(M3D_readasctxt) "// &
     & "Failed to allocate array memory.", nerr )
  return
  
end subroutine M3D_readasctxt
!
!**********************************************************************!




!**********************************************************************!
!
subroutine M3D_CrossProduct(v1,v2,vcp)

  implicit none
  
  real*8, intent(in)  :: v1(3), v2(3)
  real*8, intent(out) :: vcp(3)
  
  vcp(1) = v1(2)*v2(3) - v1(3)*v2(2)
  vcp(2) = v1(3)*v2(1) - v1(1)*v2(3)
  vcp(3) = v1(1)*v2(2) - v1(2)*v2(1)
  
  return
  
end subroutine M3D_CrossProduct
!
!**********************************************************************!



!**********************************************************************!
!
subroutine M3D_CellDimension(s1,s2,s3)

  implicit none
  
  real*4, intent(out) :: s1, s2, s3
  
  s1 = real( dsqrt( sum( M3D_b1*M3D_b1 ) ) * dble(M3D_n1), kind=4 )
  s2 = real( dsqrt( sum( M3D_b2*M3D_b2 ) ) * dble(M3D_n2), kind=4 )
  s3 = real( dsqrt( sum( M3D_b3*M3D_b3 ) ) * dble(M3D_n3), kind=4 )
  
  return
  
end subroutine M3D_CellDimension
!
!**********************************************************************!



!**********************************************************************!
!
subroutine M3D_CellVolume(vol)

  implicit none
  
  real*4, intent(out) :: vol
  real*8 :: v12(3)
  
  vol = 0.0
  v12 = 0.0D+00
  call M3D_CrossProduct( M3D_b1*dble(M3D_n1), M3D_b2*dble(M3D_n2), v12)
  vol = real( sum( v12*M3D_b3 ) * dble(M3D_n3), kind=4 )
  
  return
  
end subroutine M3D_CellVolume
!
!**********************************************************************!



!**********************************************************************!
!
subroutine M3D_CellAngles(alpha,beta,gamma)

  implicit none
  
  real*4, intent(out) :: alpha, beta, gamma
  
  real*8 :: lb1, lb2, lb3, sp12, sp23, sp13
  
  lb1 = dsqrt( sum( M3D_b1*M3D_b1 ) )
  lb2 = dsqrt( sum( M3D_b2*M3D_b2 ) )
  lb3 = dsqrt( sum( M3D_b3*M3D_b3 ) )
  
  sp12 = sum( M3D_b1*M3D_b2 ) / lb1 / lb2
  sp23 = sum( M3D_b2*M3D_b3 ) / lb2 / lb3
  sp13 = sum( M3D_b1*M3D_b3 ) / lb1 / lb3
  
  gamma = real( dacos( sp12 ), kind=4 )
  alpha = real( dacos( sp23 ), kind=4 )
  beta  = real( dacos( sp13 ), kind=4 )
  
  return
  
end subroutine M3D_CellAngles
!
!**********************************************************************!




!**********************************************************************!
!
subroutine M3D_PostCellInfo()

  implicit none
  
  real*4 :: rad2deg = 57.29578
  real*4 :: sx, sy, sz, sdx, sdy, sdz, a, b, c
  character(len=M3D_ll) :: smsg
  
  call M3D_CellDimension( sx, sy, sz )
  sdx = sx / real(M3D_n1)
  sdy = sy / real(M3D_n2)
  sdz = sz / real(M3D_n3)
  call M3D_CellAngles(a,b,c)
  
  call M3D_Message("Information on the current 3D potential cell.")
  
  write(unit=smsg,fmt='(A,3F12.8)') "Cell size (x,y,z) in nm: ", &
     &  sx, sy, sz 
  call M3D_Message(trim(smsg))
  
  write(unit=smsg,fmt='(A,3F9.4)') &
     & "Cell angles (alpha,beta,gamma) in deg: ", a*rad2deg, &
     & b*rad2deg, c*rad2deg
  call M3D_Message(trim(smsg))
  
  write(unit=smsg,fmt='(A,3I5)') &
     & "Cell sampling (nx,ny,nz) in pixels: ", M3D_n1, M3D_n2, M3D_n3 
  call M3D_Message(trim(smsg))
  
  write(unit=smsg,fmt='(A,3F12.8)') &
     & "Cell sampling rate (dx,dy,dz) in nm/pix: ", sdx, sdy, sdz 
  call M3D_Message(trim(smsg))
  
  
  return

end subroutine M3D_PostCellInfo
!
!**********************************************************************!







!**********************************************************************!
!
subroutine M3D_FFT1D(dat, n, nft, dir)

  implicit none
  
  integer*4, intent(in) :: n, nft
  character(len=*), intent(in) :: dir
  complex*8, intent(inout) :: dat(nft)
    
  integer*4 :: i
  character*40 :: direction
  external :: ODDCC128, ODDCC256, ODDCC512, ODDCC1024
  external :: ODDCC2048, ODDCC4096, ODDCC8192 ! link FFTs.f
  
  if (n<1 .or. n>nft) then
    call M3D_ErrorMessage("(M3D_FFT1D) Invalid array size.", 1)
    return
  end if
  direction = dir
  
  i = 0
  select case (nft)
  case (128)
    call ODDCC128(dat,n,direction)
    i = nft
  case (256)
    call ODDCC256(dat,n,direction)
    i = nft
  case (512)
    call ODDCC512(dat,n,direction)
    i = nft
  case (1024)
    call ODDCC1024(dat,n,direction)
    i = nft
  case (2048)
    call ODDCC2048(dat,n,direction)
    i = nft
  case (4096)
    call ODDCC4096(dat,n,direction)
    i = nft
  case (8192)
    call ODDCC8192(dat,n,direction)
    i = nft
  end select
  
  if (i==0) call M3D_ErrorMessage("(M3D_FFT1D) Unsupported size.", 2)
  
  return
  
end subroutine M3D_FFT1D
!
!**********************************************************************!




!**********************************************************************!
!
subroutine M3D_FFT2D(dat, nx, ny, nft, dir)

  implicit none
  
  integer*4, intent(in) :: nx, ny, nft
  character(len=*), intent(in) :: dir
  complex*8, intent(inout) :: dat(nft,nft)
    
  integer*4 :: i
  character*40 :: direction
  external :: ODDCC128S, ODDCC256S, ODDCC512S, ODDCC1024S
  external :: ODDCC2048S, ODDCC4096S, ODDCC8192S ! link FFTs.f
  
  if (nx<1 .or. nx>nft .or. ny<1 .or. ny>nft) then
    call M3D_ErrorMessage("(M3D_FFT2D) Invalid array size.", 1)
    return
  end if
  direction = dir
  
  i = 0
  select case (nft)
  case (128)
    call ODDCC128S(dat,nx,ny,direction)
    i = nft
  case (256)
    call ODDCC256S(dat,nx,ny,direction)
    i = nft
  case (512)
    call ODDCC512S(dat,nx,ny,direction)
    i = nft
  case (1024)
    call ODDCC1024S(dat,nx,ny,direction)
    i = nft
  case (2048)
    call ODDCC2048S(dat,nx,ny,direction)
    i = nft
  case (4096)
    call ODDCC4096S(dat,nx,ny,direction)
    i = nft
  case (8192)
    call ODDCC8192S(dat,nx,ny,direction)
    i = nft
  end select
  
  if (i==0) call M3D_ErrorMessage("(M3D_FFT2D) Unsupported size.", 2)
  
  return
  
end subroutine M3D_FFT2D
!
!**********************************************************************!




!**********************************************************************!
!
subroutine M3D_FFT3D(datr, datk, nx, ny, nz, dir)
!
! complex*8 -> complex*8     3-dimesional
! Fourier-Transforms between datr and datk in direction dir.
! datk is scrambled and x-y tranposed.
! The respective input is left unchanged.
!

  implicit none
  
  integer*4, intent(in) :: nx, ny, nz
  character(len=*), intent(in) :: dir
  complex*8, intent(inout) :: datr(nx,ny,nz), datk(ny,nx,nz)
    
  integer*4 :: i, j, k, nft1, nft2, nalloc
  character*40 :: direction
  character(len=80) :: sline
  complex*8, allocatable :: ctmp2d(:,:), ctmp1d(:) ! transformation helper arrays
  external :: UPPERCASE
  
  if (nx<1 .or. ny<1 .or. nz<1) goto 101
  call UPPERCASE(dir, direction)
  
  nft1 = 2**CEILING( LOG( real( nz ) )/LOG(2.0) ) ! next 2^N above nz
  nft1 = max(nft1, 256)
  if (nft1>M3D_FFTBOUND) goto 114
  nft2 = 2**CEILING( LOG( real( max(nx,ny) ) )/LOG(2.0) ) ! next 2^N above max(nx,ny)
  nft2 = max(nft2, 256)
  if (nft2>M3D_FFTBOUND) goto 114
  
  ! - prepare transform helper arrays
  allocate(ctmp2d(nft2,nft2),ctmp1d(nft1),stat=nalloc)
  if (nalloc/=0) goto 115
  
  ! - different sequence for forward and backward transformations
  if (direction(1:3)=='FOR') then ! FORWARD FT
    !
    call M3D_Message("Running forward 3D FFT ...")
    ! - all x-y planes first
    !   to temp. fill the result array and keep the input array
    ctmp2d = cmplx(0.0,0.0)
    do k=1, nz ! LOOP planes
      do j=1, ny ! LOOP rows & transfer to helper
        ctmp2d(1:nx,j) = datr(1:nx,j,k)
      end do
      call M3D_FFT2D(ctmp2d, nx, ny, nft2, direction)
      ! - now x and y is swaped 
      do j=1, nx ! LOOP rows & transfer to result
        datk(1:ny,j,k) = ctmp2d(1:ny,j)
      end do
    end do
    ! - all z lines second
    ctmp1d = cmplx(0.0,0.0)
    do j=1, nx ! LOOP rows
      do i=1, ny ! LOOP columns & transfer
        do k=1, nz
          ctmp1d(k) = datk(i,j,k)
        end do
        call M3D_FFT1D(ctmp1d, nz, nft1, direction)
        do k=1, nz
          datk(i,j,k) = ctmp1d(k)
        end do
      end do
    end do
    ! - done
    call M3D_Message("Forward 3D FFT done.")
    !
  else if (direction(1:3)=='BAC') then ! BACKWARD FT
    !
    call M3D_Message("Running backward 3D FFT ...")
    ! - all x-y planes first
    !   to temp. fill the result array and keep the input array
    ctmp2d = cmplx(0.0,0.0)
    do k=1, nz ! LOOP planes
      do j=1, nx ! LOOP rows & transfer to helper
        ctmp2d(1:ny,j) = datk(1:ny,j,k)
      end do
      call M3D_FFT2D(ctmp2d, nx, ny, nft2, direction)
      ! - now x and y is swaped 
      do j=1, ny ! LOOP rows & transfer to result
        datr(1:nx,j,k) = ctmp2d(1:nx,j)
      end do
    end do
    ! - all z lines second
    ctmp1d = cmplx(0.0,0.0)
    do j=1, ny ! LOOP rows
      do i=1, nx ! LOOP columns & transfer
        do k=1, nz
          ctmp1d(k) = datr(i,j,k)
        end do
        call M3D_FFT1D(ctmp1d, nz, nft1, direction)
        do k=1, nz
          datr(i,j,k) = ctmp1d(k)
        end do
      end do
    end do
    ! - done
    call M3D_Message("Backward 3D FFT done.")
    !
  else ! UNKNOW DIRECTION
    goto 103
  end if
  !
  ! - get rid of the helper arrays
  deallocate(ctmp2d,ctmp1d,stat=nalloc)
  if (nalloc/=0) goto 116
    
  return
  
! error handling
101 call M3D_ErrorMessage("(M3D_FFT3D) Invalid array size.", 1)
  return
103 call M3D_ErrorMessage("(M3D_FFT3D) Invalid direction input.", 3)
  return
114 call M3D_ErrorMessage("(M3D_FFT3D) FFT size out of range (8192).", 14)
  return
115 call M3D_ErrorMessage("(M3D_FFT3D) Memory allocation failed.", 15)
  return
116 call M3D_ErrorMessage("(M3D_FFT3D) Memory deallocation failed.", 16)
  return
  
end subroutine M3D_FFT3D
!
!**********************************************************************!







!**********************************************************************!
!
subroutine M3D_OrthoGamma(maxdim, nerr)
  
  implicit none
  
  real*4, parameter :: pihalf = 1.5707963
  real*8, parameter :: twopi = 6.2831853071796
  
  integer*4, intent(in) :: maxdim
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nalloc
  integer*4 :: i, j, k, l, i1, j1, k1, l1
  integer*4 :: kmax, lmax, kmin, lmin, kuse, luse
  integer*4 :: n1, n2, n3, n1p, n2p, nft1
  real*4 :: gamma
  complex*8, allocatable :: potcopy(:,:,:)
  real*8 :: b1(3), b2(3), e1(3), e2(3), e2p(3)
  real*8 :: lb1, lb2, sp12, le2p, lb2p, vtmp, vpha
  real*8 :: dioph, dmin, rtmp, rmin
  complex*8, allocatable :: cshft(:)
  character(len=M3D_ll) :: sinfo
  
  nerr = 0
  nalloc = 0
  b1 = M3D_b1
  b2 = M3D_b2
  n1 = M3D_n1
  n2 = M3D_n2
  n3 = M3D_n3
  
  if ( .not.allocated(M3D_pot) ) then
    nerr = 1
    call M3D_ErrorMessage("(M3D_OrthoGamma) No 3D potential is "// &
     & "available.", nerr)
    return
  end if
  
  if ( maxdim < n1 .or. maxdim < n2 ) then
    nerr = 2
    call M3D_ErrorMessage("(M3D_OrthoGamma) At least one plane "// &
     & "sampling is >= the max. dimension parameter.", nerr)
    return
  end if
  
  lb1 = dsqrt( sum( b1*b1 ) )
  lb2 = dsqrt( sum( b2*b2 ) )
  
  if ( 1.0D-34 > lb1 .or. 1.0D-34 > lb2 ) then
    nerr = 3
    call M3D_ErrorMessage("(M3D_OrthoGamma) At least one basis "// &
     & "vector has zero length.", nerr)
    return
  end if
  
  e1 = b1 / lb1
  e2 = b2 / lb2
  sp12 = sum( e1*e2 )
  gamma = real( dacos( sp12 ), kind=4 )
  
  if ( M3D_orththr > abs(abs(gamma)-pihalf) ) then
    call M3D_Message("(M3D_OrthoGamma) Gamma is already close to pi/2.")
    call M3D_Message("(M3D_OrthoGamma) No further orthogonalization.")
    return
  end if
  
  if ( M3D_orththr > abs(gamma) ) then
    nerr = 4
    call M3D_ErrorMessage("(M3D_OrthoGamma) The two basis vectors "// &
     & "defining the plane are colinear.", nerr)
    return
  end if 
  
  ! Construct the new basis vector 2' in the original 1-2-plane,
  ! which is orthogonal to the basis vector 1.
  ! Do this by using the Gram-Smith rule
  e2p = e2 - e1*sp12
  le2p = dsqrt( sum( e2p*e2p ) )
  e2p = e2p / le2p
  
  ! Construction of a repeated periodic cell with
  ! basis vectors k*lb1*e1 and l*lb2p*e2p, where k and l are
  ! integer numbers.
  ! Find (k,l) such that 0 = k*lb1 + l*lb2*(e1.e2)
  !                      0 = k*lb1 + l*lb2*sp12
  !                      0 = k + l * sp12*lb2/lb1
  ! Counting l positive from 1 to lmax, is there a k in (-kmax,kmax)
  ! fullfilling the Diophantine equation above?
  !
  ! Or which pair (k,l) produces the smallest number dioph in
  !                  dioph = k/l + sp12*lb2/lb1 (l/=0 and lb1/=0)
  !
  ! get the maximum repeat until maxdim is reached.
  kmax = floor( real(maxdim)/real(n1) )
  lmax = floor( real(maxdim)/real(n2) )
  ! Calculate the missmatch to the solutions of the Diophantine
  ! equation.
  ! Search the acceptable solution with smallest (k*n1)^2+(l*n2)^2.
  vtmp = sp12*lb2/lb1
  rmin = (dble(n1)*dble(kmax))**2 + (dble(n2)*dble(lmax))**2
  k1 = 0
  l1 = 0
  kmin = 0
  lmin = 0
  dmin = 1.0D+100
  do l=1, lmax
    do k=-kmax, kmax
      dioph = dble(k) + dble(l)*vtmp
      if (dabs(dioph) < dmin) then ! Store the best solution
        dmin = dabs(dioph)
        kmin = k
        lmin = l
      end if
      if (dabs(dioph) < M3D_diophthr) then ! Store the accaptable solution ...
        rtmp = (dble(n1)*dble(k))**2 + (dble(n2)*dble(l))**2
        if (rtmp<rmin) then ! ... with minimum sampling
          rmin = rtmp
          k1 = k
          l1 = l
        end if
      end if
    end do
  end do
  
  if (k1==0 .and. l1==0) then ! No acceptable solution was found.
    call M3D_Message("(M3D_OrthoGamma) WARNING: Orthogonal problem.")
    call M3D_Message("(M3D_OrthoGamma) No acceptable solution of the " &
     &   //"Diophantine equation was found to orthogonalize gamma.")
    call M3D_Message("(M3D_OrthoGamma) Using the best solution:")
    write(unit=sinfo, fmt='(A,I4,A,I4)') &
     &   "- basis vector (1&2) repeats: ", kmin, ", ", lmin
    call M3D_Message("(M3D_OrthoGamma) "//trim(sinfo))
    write(unit=sinfo, fmt='(A,E12.4)') &
     &   "- Diophantine mismatch: ", dmin
    call M3D_Message("(M3D_OrthoGamma) "//trim(sinfo))
    call M3D_Message("(M3D_OrthoGamma) Quasi-orthogonalizing gamma " &
     &   //"with a non-optimum solution.")
    kuse = kmin
    luse = lmin
  else ! The best acceptable solution was identified.
    call M3D_Message("(M3D_OrthoGamma) An acceptable solution of the " &
     &   //"Diophantine equation was found to orthogonalize gamma.")
    call M3D_Message("(M3D_OrthoGamma) Using the best solution:")
    write(unit=sinfo, fmt='(A,I4,A,I4)') &
     &   "- basis vector (1&2) repeats: ", k1, ", ", l1
    call M3D_Message("(M3D_OrthoGamma) "//trim(sinfo))
    write(unit=sinfo, fmt='(A,E12.4)') &
     &   "- Diophantine mismatch: ", dble(k1) + dble(l1)*vtmp
    call M3D_Message("(M3D_OrthoGamma) "//trim(sinfo))
    kuse = k1
    luse = l1
  end if
  !
  ! Catch the no-change solution
  if (abs(kuse)==1.and.luse==1) then
    call M3D_Message("(M3D_OrthoGamma) The chosen solution will not " &
     &   //"change the sampling of the potential.")
    return
  end if
  !
  ! Calculate the physical length of the new basis vector 2 e2p
  lb2p = lb2*dsqrt(1.0D+00 - sp12*sp12)
  ! Calculate the sampling of the new array
  n1p = kuse*n1
  n2p = abs(luse)*n2
  !
  ! Make a copy of the current potential
  allocate(potcopy(n1,n2,n3), stat=nalloc)
  if (nalloc/=0) then
    nerr = 5
    call M3D_ErrorMessage("(M3D_OrthoGamma) Memory allocation "// &
     & "failed.", nerr)
    return
  end if
  ! Copy the current potential data
  potcopy = M3D_pot
  ! Re-allocate the module potential memory to hold the new sampling
  if (allocated(M3D_pot)) deallocate(M3D_pot, stat=nalloc)
  if (nalloc/=0) then
    nerr = 6
    call M3D_ErrorMessage("(M3D_OrthoGamma) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  nft1 = 2**CEILING( LOG( real( n1p ) )/LOG(2.0) ) ! next 2^N above the new 1st dimension
  nft1 = max(nft1, 256) ! min FFT size
  allocate(M3D_pot(n1p,n2p,n3), cshft(nft1), stat=nalloc)
  if (nalloc/=0) then
    nerr = 7
    call M3D_ErrorMessage("(M3D_OrthoGamma) Memory allocation "// &
     & "failed.", nerr)
    return
  end if
  M3D_pot = cmplx(0.0,0.0)
  cshft = cmplx(0.0,0.0)
  !
  ! Shift-copy the resampled data
  do k=1, n3
    do j=1, n2p
      j1 = 1 + modulo(j-1,n2) ! coordinate in the old array
      cshft = cmplx(0.0,0.0)
      do i=1, n1p
        i1 = 1 + modulo(i-1,n1) ! coordinate in the old array
        cshft(i) = potcopy( i1, j1, k ) ! repeated copy of a line
      end do
      ! We need to shift the data along dimension 1 such that the
      ! coordinate 0 resides on the first pixel.
      ! - Determine the shift value in pixels
      vtmp = 1.0D+00*dble(j-1)*sp12/real(n1p)
      ! 1D-shift by Fourier interpolation:
      ! - FT of the new line
      call M3D_FFT1D(cshft,n1p,nft1,'for')
      ! - Apply the shift by multiplication with a phase ramp EXP( -2*Pi*I*K*D )
      do i=1, n1p
        i1 = i-1 ! positive K
        if (i1>=n1p/2) i1 = i1 - n1p ! negative K
        vpha = twopi*real(i1)*vtmp ! Phase = 2*Pi*K*D
        cshft(i) = cshft(i)*exp( cmplx(0.0, -vpha) ) ! Phase ramp multiplication
      end do
      ! - iFT of the shifted line
      call M3D_FFT1D(cshft,n1p,nft1,'bac')
      ! - Copy the shifted data (real-part) into the new array.
      if (M3D_intype==0) then ! complex input potentals remain complex
        M3D_pot(1:n1p,j,k) = cshft(1:n1p)
      else ! real input potentials should remain real
        do i=1, n1p
          M3D_pot(i,j,k) = cmplx(real(cshft(i)),0.0)
        end do
      end if
    end do
  end do
  
  if (allocated(potcopy)) deallocate(potcopy, stat=nalloc)
  if (nalloc/=0) then
    nerr = 8
    call M3D_ErrorMessage("(M3D_OrthoGamma) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  if (allocated(cshft)) deallocate(cshft, stat=nalloc)
  if (nalloc/=0) then
    nerr = 9
    call M3D_ErrorMessage("(M3D_OrthoGamma) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  
  ! update the potential sampling
  M3D_n1 = n1p
  M3D_n2 = n2p
  ! M3D_n3 = M3D_n3
  
  M3D_b1 = dble(kuse)*M3D_b1
  M3D_b2 = lb2p*e2p
  ! M3D_b3 = M3D_b3
  
  return
  
end subroutine M3D_OrthoGamma
!
!**********************************************************************!







!**********************************************************************!
!
subroutine M3D_OrthoBeta(maxdim, nerr)
  
  implicit none
  
  real*4, parameter :: pihalf = 1.5707963
  real*8, parameter :: twopi = 6.2831853071796
  
  integer*4, intent(in) :: maxdim
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nalloc
  integer*4 :: i, j, k, l, i1, j1, k1, l1
  integer*4 :: kmax, lmax, kmin, lmin, kuse, luse
  integer*4 :: n1, n2, n3, n1p, n3p, nft1
  real*4 :: beta
  complex*8, allocatable :: potcopy(:,:,:)
  real*8 :: b1(3), b3(3), e1(3), e3(3), e3p(3)
  real*8 :: lb1, lb3, sp13, le3p, lb3p, vtmp, vpha
  real*8 :: dioph, dmin, rtmp, rmin
  complex*8, allocatable :: cshft(:)
  character(len=M3D_ll) :: sinfo
  
  nerr = 0
  nalloc = 0
  b1 = M3D_b1
  b3 = M3D_b3
  n1 = M3D_n1
  n2 = M3D_n2
  n3 = M3D_n3
  
  if ( .not.allocated(M3D_pot) ) then
    nerr = 1
    call M3D_ErrorMessage("(M3D_OrthoBeta) No 3D potential is "// &
     & "available.", nerr)
    return
  end if
  
  if ( maxdim < n1 .or. maxdim < n3 ) then
    nerr = 2
    call M3D_ErrorMessage("(M3D_OrthoBeta) At least one plane "// &
     & "sampling is >= the max. dimension parameter.", nerr)
    return
  end if
  
  lb1 = dsqrt( sum( b1*b1 ) )
  lb3 = dsqrt( sum( b3*b3 ) )
  
  if ( 1.0D-34 > lb1 .or. 1.0D-34 > lb3 ) then
    nerr = 3
    call M3D_ErrorMessage("(M3D_OrthoBeta) At least one basis "// &
     & "vector has zero length.", nerr)
    return
  end if
  
  e1 = b1 / lb1
  e3 = b3 / lb3
  sp13 = sum( e1*e3 )
  beta = real( dacos( sp13 ), kind=4 )
  
  if ( M3D_orththr > abs(abs(beta)-pihalf) ) then
    call M3D_Message("(M3D_OrthoBeta) Beta is already close to pi/2.")
    call M3D_Message("(M3D_OrthoBeta) No further orthogonalization.")
    return
  end if
  
  if ( M3D_orththr > abs(beta) ) then
    nerr = 4
    call M3D_ErrorMessage("(M3D_OrthoBeta) The two basis vectors "// &
     & "defining the plane are colinear.", nerr)
    return
  end if 
  
  ! Construct the new basis vector 3' in the original 1-3-plane,
  ! which is orthogonal to the basis vector 1.
  ! Do this by using the Gram-Smith rule
  e3p = e3 - e1*sp13
  le3p = dsqrt( sum( e3p*e3p ) )
  e3p = e3p / le3p
  
  ! Construction of a repeated periodic cell with
  ! basis vectors k*lb1*e1 and l*lb3p*e3p, where k and l are
  ! integer numbers.
  ! Find (k,l) such that 0 = k*lb1 + l*lb3*(e1.e3)
  !                      0 = k*lb1 + l*lb3*sp13
  !                      0 = k + l * sp13*lb3/lb1
  ! Counting l positive from 1 to lmax, is there a k in (-kmax,kmax)
  ! fullfilling the Diophantine equation above?
  !
  ! Or which pair (k,l) produces the smallest number dioph in
  !                  dioph = k/l + sp13*lb3/lb1 (l/=0 and lb1/=0)
  !
  ! get the maximum repeat until maxdim is reached.
  kmax = floor( real(maxdim)/real(n1) )
  lmax = floor( real(maxdim)/real(n3) )
  ! Calculate the missmatch to the solutions of the Diophantine
  ! equation.
  ! Search the acceptable solution with smallest (k*n1)^2+(l*n3)^2.
  vtmp = sp13*lb3/lb1
  rmin = (dble(n1)*dble(kmax))**2 + (dble(n3)*dble(lmax))**2
  k1 = 0
  l1 = 0
  kmin = 0
  lmin = 0
  dmin = 1.0D+100
  do l=1, lmax
    do k=-kmax, kmax
      dioph = dble(k) + dble(l)*vtmp
      if (dabs(dioph) < dmin) then ! Store the best solution
        dmin = dabs(dioph)
        kmin = k
        lmin = l
      end if
      if (dabs(dioph) < M3D_diophthr) then ! Store the accaptable solution ...
        rtmp = (dble(n1)*dble(k))**2 + (dble(n3)*dble(l))**2
        if (rtmp<rmin) then ! ... with minimum sampling
          rmin = rtmp
          k1 = k
          l1 = l
        end if
      end if
    end do
  end do
  
  if (k1==0 .and. l1==0) then ! No acceptable solution was found.
    call M3D_Message("(M3D_OrthoBeta) WARNING: Orthogonal problem.")
    call M3D_Message("(M3D_OrthoBeta) No acceptable solution of the " &
     &   //"Diophantine equation was found to orthogonalize beta.")
    call M3D_Message("(M3D_OrthoBeta) Using the best solution:")
    write(unit=sinfo, fmt='(A,I4,A,I4)') &
     &   "- basis vector (1&3) repeats: ", kmin, ", ", lmin
    call M3D_Message("(M3D_OrthoBeta) "//trim(sinfo))
    write(unit=sinfo, fmt='(A,E12.4)') &
     &   "- Diophantine mismatch: ", dmin
    call M3D_Message("(M3D_OrthoBeta) "//trim(sinfo))
    call M3D_Message("(M3D_OrthoBeta) Quasi-orthogonalizing beta " &
     &   //"with a non-optimum solution.")
    kuse = kmin
    luse = lmin
  else ! The best acceptable solution was identified.
    call M3D_Message("(M3D_OrthoBeta) An acceptable solution of the " &
     &   //"Diophantine equation was found to orthogonalize beta.")
    call M3D_Message("(M3D_OrthoBeta) Using the best solution:")
    write(unit=sinfo, fmt='(A,I4,A,I4)') &
     &   "- basis vector (1&3) repeats: ", k1, ", ", l1
    call M3D_Message("(M3D_OrthoBeta) "//trim(sinfo))
    write(unit=sinfo, fmt='(A,E12.4)') &
     &   "- Diophantine mismatch: ", dble(k1) + dble(l1)*vtmp
    call M3D_Message("(M3D_OrthoBeta) "//trim(sinfo))
    kuse = k1
    luse = l1
  end if
  !
  ! Catch the no-change solution
  if (abs(kuse)==1.and.luse==1) then
    call M3D_Message("(M3D_OrthoBeta) The chosen solution will not " &
     &   //"change the sampling of the potential.")
    return
  end if
  !
  ! Calculate the physical length of the new basis vector 3 e3p
  lb3p = lb3*dsqrt(1.0D+00 - sp13*sp13)
  ! Calculate the sampling of the new array
  n1p = kuse*n1
  n3p = abs(luse)*n3
  !
  ! Make a copy of the current potential
  allocate(potcopy(n1,n2,n3), stat=nalloc)
  if (nalloc/=0) then
    nerr = 5
    call M3D_ErrorMessage("(M3D_OrthoBeta) Memory allocation "// &
     & "failed.", nerr)
    return
  end if
  ! Copy the current potential data
  potcopy = M3D_pot
  ! Re-allocate the module potential memory to hold the new sampling
  if (allocated(M3D_pot)) deallocate(M3D_pot, stat=nalloc)
  if (nalloc/=0) then
    nerr = 6
    call M3D_ErrorMessage("(M3D_OrthoBeta) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  nft1 = 2**CEILING( LOG( real( n1p ) )/LOG(2.0) ) ! next 2^N above the new 1st dimension
  nft1 = max(nft1, 256) ! min FFT size
  allocate(M3D_pot(n1p,n2,n3p), cshft(nft1), stat=nalloc)
  if (nalloc/=0) then
    nerr = 7
    call M3D_ErrorMessage("(M3D_OrthoBeta) Memory allocation "// &
     & "failed.", nerr)
    return
  end if
  M3D_pot = cmplx(0.0,0.0)
  cshft = cmplx(0.0,0.0)
  !
  ! Shift-copy the resampled data
  do j=1, n2
    do k=1, n3p
      k1 = 1 + modulo(k-1,n3) ! coordinate in the old array
      cshft = cmplx(0.0,0.0)
      do i=1, n1p
        i1 = 1 + modulo(i-1,n1) ! coordinate in the old array
        cshft(i) = potcopy( i1, j, k1 ) ! repeated copy of a line
      end do
      ! We need to shift the data along dimension 1 such that the
      ! coordinate 0 resides on the first pixel.
      ! - Determine the shift value in pixels
      vtmp = 1.0D+00*dble(k-1)*sp13/real(n1p)
      ! 1D-shift by Fourier interpolation:
      ! - FT of the new line
      call M3D_FFT1D(cshft,n1p,nft1,'for')
      ! - Apply the shift by multiplication with a phase ramp EXP( -2*Pi*I*K*D )
      do i=1, n1p
        i1 = i-1 ! positive K
        if (i1>=n1p/2) i1 = i1 - n1p ! negative K
        vpha = twopi*real(i1)*vtmp ! Phase = 2*Pi*K*D
        cshft(i) = cshft(i)*exp( cmplx(0.0, -vpha) ) ! Phase ramp multiplication
      end do
      ! - iFT of the shifted line
      call M3D_FFT1D(cshft,n1p,nft1,'bac')
      ! - Copy the shifted data into the new array.
      if (M3D_intype==0) then ! complex input potentals remain complex
        M3D_pot(1:n1p,j,k) = cshft(1:n1p)
      else ! real input potentials should remain real
        do i=1, n1p
          M3D_pot(i,j,k) = cmplx(real(cshft(i)),0.0)
        end do
      end if
    end do
  end do
  
  if (allocated(potcopy)) deallocate(potcopy, stat=nalloc)
  if (nalloc/=0) then
    nerr = 8
    call M3D_ErrorMessage("(M3D_OrthoBeta) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  if (allocated(cshft)) deallocate(cshft, stat=nalloc)
  if (nalloc/=0) then
    nerr = 9
    call M3D_ErrorMessage("(M3D_OrthoBeta) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  
  ! update the potential sampling
  M3D_n1 = n1p
  ! M3D_n2 = M3D_n2
  M3D_n3 = n3p
  
  M3D_b1 = dble(kuse)*M3D_b1
  ! M3D_b2 = M3D_b2
  M3D_b3 = lb3p*e3p
  
  
  return
  
end subroutine M3D_OrthoBeta
!
!**********************************************************************!




!**********************************************************************!
!
subroutine M3D_OrthoAlpha(maxdim, nerr)
  
  implicit none
  
  real*4, parameter :: pihalf = 1.5707963
  real*8, parameter :: twopi = 6.2831853071796
  
  integer*4, intent(in) :: maxdim
  integer*4, intent(inout) :: nerr
  
  integer*4 :: nalloc
  integer*4 :: i, j, k, l, i1, j1, k1, l1
  integer*4 :: kmax, lmax, kmin, lmin, kuse, luse
  integer*4 :: n1, n2, n3, n2p, n3p, nft1
  real*4 :: alpha
  complex*8, allocatable :: potcopy(:,:,:)
  real*8 :: b2(3), b3(3), e2(3), e3(3), e3p(3)
  real*8 :: lb2, lb3, sp23, le3p, lb3p, vtmp, vpha
  real*8 :: dioph, dmin, rtmp, rmin
  complex*8, allocatable :: cshft(:)
  character(len=M3D_ll) :: sinfo
  
  nerr = 0
  nalloc = 0
  b2 = M3D_b2
  b3 = M3D_b3
  n1 = M3D_n1
  n2 = M3D_n2
  n3 = M3D_n3
  
  if ( .not.allocated(M3D_pot) ) then
    nerr = 1
    call M3D_ErrorMessage("(M3D_OrthoAlpha) No 3D potential is "// &
     & "available.", nerr)
    return
  end if
  
  if ( maxdim < n2 .or. maxdim < n3 ) then
    nerr = 2
    call M3D_ErrorMessage("(M3D_OrthoAlpha) At least one plane "// &
     & "sampling is >= the max. dimension parameter.", nerr)
    return
  end if
  
  lb2 = dsqrt( sum( b2*b2 ) )
  lb3 = dsqrt( sum( b3*b3 ) )
  
  if ( 1.0D-34 > lb2 .or. 1.0D-34 > lb3 ) then
    nerr = 3
    call M3D_ErrorMessage("(M3D_OrthoAlpha) At least one basis "// &
     & "vector has zero length.", nerr)
    return
  end if
  
  e2 = b2 / lb2
  e3 = b3 / lb3
  sp23 = sum( e2*e3 )
  alpha = real( dacos( sp23 ), kind=4 )
  
  if ( M3D_orththr > abs(abs(alpha)-pihalf) ) then
    call M3D_Message("(M3D_OrthoAlpha) Alpha is already close to pi/2.")
    call M3D_Message("(M3D_OrthoAlpha) No further orthogonalization.")
    return
  end if
  
  if ( M3D_orththr > abs(alpha) ) then
    nerr = 4
    call M3D_ErrorMessage("(M3D_OrthoAlpha) The two basis vectors "// &
     & "defining the plane are colinear.", nerr)
    return
  end if 
  
  ! Construct the new basis vector 3' in the original 2-3-plane,
  ! which is orthogonal to the basis vector 2.
  ! Do this by using the Gram-Smith rule.
  e3p = e3 - e2*sp23
  le3p = dsqrt( sum( e3p*e3p ) )
  e3p = e3p / le3p
  
  ! Construction of a repeated periodic cell with
  ! basis vectors k*lb2*e2 and l*lb3p*e3p, where k and l are
  ! integer numbers.
  ! Find (k,l) such that 0 = k*lb2 + l*lb3*(e2.e3)
  !                      0 = k*lb2 + l*lb3*sp23
  !                      0 = k + l * sp23*lb3/lb2
  ! Counting l positive from 1 to lmax, is there a k in (-kmax,kmax)
  ! fullfilling the Diophantine equation above?
  !
  ! Or which pair (k,l) produces the smallest number dioph in
  !                  dioph = k/l + sp23*lb3/lb2 (l/=0 and lb2/=0)
  !
  ! get the maximum repeat until maxdim is reached.
  kmax = floor( real(maxdim)/real(n2) )
  lmax = floor( real(maxdim)/real(n3) )
  ! Calculate the missmatch to the solutions of the Diophantine
  ! equation.
  ! Search the acceptable solution with smallest (k*n2)^2+(l*n3)^2.
  vtmp = sp23*lb3/lb2
  rmin = (dble(n2)*dble(kmax))**2 + (dble(n3)*dble(lmax))**2
  k1 = 0
  l1 = 0
  kmin = 0
  lmin = 0
  dmin = 1.0D+100
  do l=1, lmax
    do k=-kmax, kmax
      dioph = dble(k) + dble(l)*vtmp
      if (dabs(dioph) < dmin) then ! Store the best solution
        dmin = dabs(dioph)
        kmin = k
        lmin = l
      end if
      if (dabs(dioph) < M3D_diophthr) then ! Store the accaptable solution ...
        rtmp = (dble(n2)*dble(k))**2 + (dble(n3)*dble(l))**2
        if (rtmp<rmin) then ! ... with minimum sampling
          rmin = rtmp
          k1 = k
          l1 = l
        end if
      end if
    end do
  end do
  
  if (k1==0 .and. l1==0) then ! No acceptable solution was found.
    call M3D_Message("(M3D_OrthoAlpha) WARNING: Orthogonal problem.")
    call M3D_Message("(M3D_OrthoAlpha) No acceptable solution of the " &
     &   //"Diophantine equation was found to orthogonalize beta.")
    call M3D_Message("(M3D_OrthoAlpha) Using the best solution:")
    write(unit=sinfo, fmt='(A,I4,A,I4)') &
     &   "- basis vector (1&3) repeats: ", kmin, ", ", lmin
    call M3D_Message("(M3D_OrthoAlpha) "//trim(sinfo))
    write(unit=sinfo, fmt='(A,E12.4)') &
     &   "- Diophantine mismatch: ", dmin
    call M3D_Message("(M3D_OrthoAlpha) "//trim(sinfo))
    call M3D_Message("(M3D_OrthoAlpha) Quasi-orthogonalizing beta " &
     &   //"with a non-optimum solution.")
    kuse = kmin
    luse = lmin
  else ! The best acceptable solution was identified.
    call M3D_Message("(M3D_OrthoAlpha) An acceptable solution of the " &
     &   //"Diophantine equation was found to orthogonalize beta.")
    call M3D_Message("(M3D_OrthoAlpha) Using the best solution:")
    write(unit=sinfo, fmt='(A,I4,A,I4)') &
     &   "- basis vector (1&3) repeats: ", k1, ", ", l1
    call M3D_Message("(M3D_OrthoAlpha) "//trim(sinfo))
    write(unit=sinfo, fmt='(A,E12.4)') &
     &   "- Diophantine mismatch: ", dble(k1) + dble(l1)*vtmp
    call M3D_Message("(M3D_OrthoAlpha) "//trim(sinfo))
    kuse = k1
    luse = l1
  end if
  !
  ! Catch the no-change solution
  if (abs(kuse)==1.and.luse==1) then
    call M3D_Message("(M3D_OrthoAlpha) The chosen solution will not " &
     &   //"change the sampling of the potential.")
    return
  end if
  !
  ! Calculate the physical length of the new basis vector 3 e3p
  lb3p = lb3*dsqrt(1.0D+00 - sp23*sp23)
  ! Calculate the sampling of the new array
  n2p = kuse*n2
  n3p = abs(luse)*n3
  !
  ! Make a copy of the current potential
  allocate(potcopy(n1,n2,n3), stat=nalloc)
  if (nalloc/=0) then
    nerr = 5
    call M3D_ErrorMessage("(M3D_OrthoAlpha) Memory allocation "// &
     & "failed.", nerr)
    return
  end if
  ! Copy the current potential data
  potcopy = M3D_pot
  ! Re-allocate the module potential memory to hold the new sampling
  if (allocated(M3D_pot)) deallocate(M3D_pot, stat=nalloc)
  if (nalloc/=0) then
    nerr = 6
    call M3D_ErrorMessage("(M3D_OrthoAlpha) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  nft1 = 2**CEILING( LOG( real( n2p ) )/LOG(2.0) ) ! next 2^N above the new 1st dimension
  nft1 = max(nft1, 256) ! min FFT size
  allocate(M3D_pot(n1,n2p,n3p), cshft(nft1), stat=nalloc)
  if (nalloc/=0) then
    nerr = 7
    call M3D_ErrorMessage("(M3D_OrthoAlpha) Memory allocation "// &
     & "failed.", nerr)
    return
  end if
  M3D_pot = cmplx(0.0,0.0)
  cshft = cmplx(0.0,0.0)
  !
  ! Shift-copy the resampled data
  do i=1, n1
    do k=1, n3p
      k1 = 1 + modulo(k-1,n3) ! coordinate in the old array
      cshft = cmplx(0.0,0.0)
      do j=1, n2p
        j1 = 1 + modulo(j-1,n2) ! coordinate in the old array
        cshft(j) = potcopy( i, j1, k1 ) ! repeated copy of a line
      end do
      ! We need to shift the data along dimension 2 such that the
      ! coordinate 0 resides on the first pixel.
      ! - Determine the shift value in pixels
      vtmp = 1.0D+00*dble(k-1)*sp23/real(n2p)
      ! 1D-shift by Fourier interpolation:
      ! - FT of the new line
      call M3D_FFT1D(cshft,n2p,nft1,'for')
      ! - Apply the shift by multiplication with a phase ramp EXP( -2*Pi*I*K*D )
      do j=1, n2p
        j1 = j-1 ! positive K
        if (j1>=n2p/2) j1 = j1 - n2p ! negative K
        vpha = twopi*real(j1)*vtmp ! Phase = 2*Pi*K*D
        cshft(j) = cshft(j)*exp( cmplx(0.0, -vpha) ) ! Phase ramp multiplication
      end do
      ! - iFT of the shifted line
      call M3D_FFT1D(cshft,n2p,nft1,'bac')
      ! - Copy the shifted data (real-part) into the new array.
      if (M3D_intype==0) then ! complex input potentals remain complex
        do j=1, n2p
          M3D_pot(i,j,k) = cmplx(real(cshft(j)),0.0)
        end do
      else ! real input potentials should remain real
        do j=1, n2p
          M3D_pot(i,j,k) = cmplx(real(cshft(j)),0.0)
        end do
      end if
    end do
  end do
  
  if (allocated(potcopy)) deallocate(potcopy, stat=nalloc)
  if (nalloc/=0) then
    nerr = 8
    call M3D_ErrorMessage("(M3D_OrthoAlpha) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  if (allocated(cshft)) deallocate(cshft, stat=nalloc)
  if (nalloc/=0) then
    nerr = 9
    call M3D_ErrorMessage("(M3D_OrthoAlpha) Memory de-allocation "// &
     & "failed.", nerr)
    return
  end if
  
  ! update the potential sampling
  ! M3D_n1 = M3D_n1
  M3D_n2 = n2p
  M3D_n3 = n3p
  
  ! M3D_b1 = M3D_b1
  M3D_b2 = dble(kuse)*M3D_b2
  M3D_b3 = lb3p*e3p
  
  
  return
  
end subroutine M3D_OrthoAlpha
!
!**********************************************************************!




!**********************************************************************!
!
subroutine M3D_Othogonalize(maxdim, nerr)

  implicit none
  
  real*4, parameter :: pihalf = 1.5707963
  
  integer*4, intent(in) :: maxdim
  integer*4, intent(inout) :: nerr
  
  integer*4 :: n1, n2, n3, maxd(1), oseq(3), i
  real*4 :: alpha, beta, gamma, dpih(3)
  
  n1 = M3D_n1
  n2 = M3D_n2
  n3 = M3D_n3
  call M3D_CellAngles(alpha,beta,gamma)
  dpih(1) = abs( abs(alpha) - pihalf )
  dpih(2) = abs( abs(beta) - pihalf )
  dpih(3) = abs( abs(gamma) - pihalf )
  
  ! determine the orthogonalization sequence
  maxd = maxloc(dpih)
  oseq(1) = maxd(1)
  dpih(maxd(1)) = -1.0
  maxd = maxloc(dpih)
  oseq(2) = maxd(1)
  dpih(maxd(1)) = -1.0
  maxd = maxloc(dpih)
  oseq(3) = maxd(1)
  dpih(maxd(1)) = -1.0
  
  do i=1, 3
  
    select case (oseq(i))
    case (1)
      call M3D_OrthoAlpha(maxdim,nerr)
    case (2)
      call M3D_OrthoBeta(maxdim,nerr)
    case (3)
      call M3D_OrthoGamma(maxdim,nerr)
    end select ! case (oseq(i))
  
  end do
  
  return
  
end subroutine M3D_Othogonalize
!
!**********************************************************************!






!**********************************************************************!
!
subroutine M3D_DiffractionFilter(nautoap, nbuni, buniv)
!
! Applies diffraction filters to the current potential.
!
! - nautoap : flag switching an automatic aperture
! - nbuni   : flag switching a universal isotropic Debye-Waller Factor
!   buniv   : universal B_ISO parameter (nm^2)
!
  implicit none
  
  integer*4, intent(in) :: nautoap, nbuni
  real*4, intent(in) :: buniv
  
  integer*4 :: nact, ndata, nalloc
  real*4 :: sx, sy, sz, ssx, ssy, ssz
  complex*8, dimension(:,:,:), allocatable :: cpotft
  
  nact = nautoap + nbuni
  
  if (nact==0) return ! nothing to do
  
  ndata = M3D_n1 * M3D_n2 * M3D_n3
  
  if (ndata<=0) return ! no data
  
  ! allocate memory for the fourier transform of the potential
  allocate(cpotft(M3D_n2,M3D_n1,M3D_n3), stat=nalloc)
  if (nalloc/=0) goto 115
  cpotft = cmplx(0.0,0.0)
  
  ! calculate the fourier-space sampling rates
  call M3D_CellDimension(ssx, ssy, ssz )
  sx = 1.0 / ssx
  sy = 1.0 / ssy
  sz = 1.0 / ssz
  
  ! forward fourier transform
  call M3D_FFT3D(M3D_pot, cpotft, M3D_n1, M3D_n2, M3D_n3, 'FORWARD')
  
  if (nautoap/=0) then ! apply auto aperture
    call M3D_ApplyDiffApAuto(cpotft, M3D_n1, M3D_n2, M3D_n3, sx, sy, sz)
  end if
  
  if (nbuni/=0) then ! apply dwf
    call M3D_ApplyDWF(cpotft, M3D_n1, M3D_n2, M3D_n3, sx, sy, sz, buniv)
  end if
  
  ! backward fourier transform
  call M3D_FFT3D(M3D_pot, cpotft, M3D_n1, M3D_n2, M3D_n3, 'BACKWARD')
  
  deallocate(cpotft, stat=nalloc)
  if (nalloc/=0) goto 116
  
  return
  
! error handling
115 call M3D_ErrorMessage("(M3D_FFT3D) Memory allocation failed.", 15)
  return
116 call M3D_ErrorMessage("(M3D_FFT3D) Memory deallocation failed.", 16)
  return
  
end subroutine M3D_DiffractionFilter
!
!**********************************************************************!






!**********************************************************************!
!
subroutine M3D_ApplyDiffApAuto(potk, nx, ny, nz, sx, sy, sz)
!
! Applys an automatic diffraction aperture, which is round in (kx,ky).
!
  implicit none
  
  complex*8, intent(inout) :: potk(ny,nx,nz) ! fourier transform of pot
  integer*4, intent(in) :: nx, ny, nz ! sampling of potk
  real*4, intent(in) :: sx, sy, sz ! sampling rates of potk (kx, ky, kz per pixel)
  
  integer*4 :: i, j, k, i1, j1, k1, nyqx, nyqy, nyqz
  real*4 :: gz, gz2, gx, gx2, gy, gxy, gmax, gmaxz
  real*4 :: apz, apxy ! aperture values
  real*4, dimension(nx) :: agx
  real*4, dimension(ny) :: agy
  real*4, dimension(nz) :: agz
  
  nyqx = nx/2 ! x Nyquist number
  nyqy = ny/2 ! y Nyquist number
  nyqz = nz/2 ! z Nyquist number
  gmax = min(sx*real(nyqx),sy*real(nyqy)) ! smallest max. diffraction in x and y
  gmaxz = sz*real(nyqz) ! max. diffraction in z
  ! prepare sampled kx values
  do j=1, nx
    j1 = mod((j+nyqx-1),nx)-nyqx
    agx(j) = sx*real(j1)
  end do
  ! prepare sampled ky values
  do i=1, ny
    i1 = mod((i+nyqy-1),ny)-nyqy
    agy(i) = sy*real(i1)
  end do
  ! prepare sampled kz values
  do k=1, nz
    k1 = mod((k+nyqz-1),nz)-nyqz
    agz(k) = sz*real(k1)
  end do
  !
  do k=1, nz ! loop planes (wz-axis)
    gz = agz(k)
    gz2 = gz*gz
    !
    apz = 0.5 - 0.5*tanh( (sqrt(gz2)/gmaxz-0.9)*30.0 )
    !
    do j=1, nx ! loop rows (wx-axis, transposed)
      gx = agx(j)
      gx2 = gx*gx
      do i=1, ny ! loop columns (wy-axis, transposed)
        gy = agy(i)
        gxy = sqrt(gx2 + gy*gy)
        !
        apxy = 0.5 - 0.5*tanh( (gxy/gmax-0.9)*30.0 )
        !
        potk(i,j,k) = potk(i,j,k) * apxy * apz
        !
      end do
    end do
  end do
  
  return
  
end subroutine M3D_ApplyDiffApAuto
!
!**********************************************************************!





!**********************************************************************!
subroutine M3D_ApplyDWF(potk, nx, ny, nz, sx, sy, sz, Bprm)
!
! Applys a Debye-Waller factor in 3 dimensions.
!
! Bprm is given in nm^2 (isotropic) B_ISO
!
! DWF = EXP( - B_ISO * g^2 / 4 )
!
  implicit none
  
  complex*8, intent(inout) :: potk(ny,nx,nz) ! fourier transform of pot
  integer*4, intent(in) :: nx, ny, nz ! sampling of potk
  real*4, intent(in) :: sx, sy, sz ! sampling rates of potk (kx, ky, kz per pixel)
  real*4, intent(in) :: Bprm ! B_ISO parameter in nm^2
  
  integer*4 :: i, j, k, i1, j1, k1, nyqx, nyqy, nyqz
  real*4 :: gz, gz2, gx, gx2, gy, gxy2, g2, dwprm
  real*4 :: dwf ! value of the Debye-Waller factor
  real*4, dimension(nx) :: agx
  real*4, dimension(ny) :: agy
  real*4, dimension(nz) :: agz
  
  nyqx = nx/2 ! x Nyquist number
  nyqy = ny/2 ! y Nyquist number
  nyqz = nz/2 ! z Nyquist number
  dwprm = -0.25 * Bprm
  ! prepare sampled kx values
  do j=1, nx
    j1 = mod((j+nyqx-1),nx)-nyqx
    agx(j) = sx*real(j1)
  end do
  ! prepare sampled ky values
  do i=1, ny
    i1 = mod((i+nyqy-1),ny)-nyqy
    agy(i) = sy*real(i1)
  end do
  ! prepare sampled kz values
  do k=1, nz
    k1 = mod((k+nyqz-1),nz)-nyqz
    agz(k) = sz*real(k1)
  end do
  !
  do k=1, nz ! loop planes (wz-axis)
    gz = agz(k)
    gz2 = gz*gz
    !
    do j=1, nx ! loop rows (wx-axis, transposed)
      gx = agx(j)
      gx2 = gx*gx
      do i=1, ny ! loop columns (wy-axis, transposed)
        gy = agy(i)
        g2 = gz2 + gx2 + gy*gy
        !
        dwf = EXP( dwprm*g2 )
        !
        potk(i,j,k) = potk(i,j,k) * dwf
        !
      end do
    end do
  end do
  
  return
  
end subroutine M3D_ApplyDWF
!**********************************************************************!





!**********************************************************************!
!
subroutine M3D_getslice_pgr(islc, nslc, nrx, nry, fabs, ht, pgr, nerr)
!
! - Use "fabs" = 0 if  IMAG( "M3D_pot" ) can be /= 0 !
! - Use "fabs" > 0 ONLY IF  IMAG( "M3D_pot" ) == 0 !

  implicit none
  
  ! PARAMETER DECLARATIONS
  real*4, parameter :: erest = 510.998928   ! electron rest energy in eV
  real*4, parameter :: psig = 2.088656      ! 2*pi*m0*e / h**2 *10e-18 nm-2
  real*4, parameter :: pi4 = 0.7853981634   ! pi/4
  
  ! INTERFACE DECLARATIONS
  integer*4, intent(in) :: islc, nslc, nrx, nry
  real*4, intent(in) :: fabs, ht
  complex*8, intent(out) :: pgr(M3D_n1*nrx,M3D_n2*nry)
  integer*4, intent(inout) :: nerr
  
  ! LOCAL VARIABLE DECLARATIONS
  integer*4 :: ndimx, ndimy
  integer*4 :: nalloc
  integer*4 :: i, j, i1, j1, i2, j2, i3, k, k0, k1, k2
  real*4 :: csz, fz0, fz1, fdz, dz, wl, sigmae, pscal
  real*4 :: fk0, fk1, fk, fki, fkt
  complex*8 :: cabsorb, cimag, cpgr
  
  ! INITIALIZATIONS
  nerr = 0
  nalloc = 0
  ndimx = M3D_n1*nrx
  ndimy = M3D_n2*nry
  cabsorb = cmplx(1.0, fabs) ! absorption potential transfer to imaginary part
  cimag = cmplx(0.0, 1.0) ! factor to transform *I (imaginary constant)
  wl = 1.239842447 / sqrt( ht * ( 1022.0 + ht ) ) ! lambda, electron wavelength [nm]
  sigmae = psig * wl ! interaction constant [nm-1]
    
  ! PRE-CHECK PARAMETERS
  if (nslc>M3D_n3) goto 101
  if (islc<1.or.islc>nslc) goto 102
  if (fabs<0.0.or.fabs>1.0) goto 103
  
  ! ALLOCATIONS
  ! - projected potential
  if (allocated(M3D_slcpot)) deallocate(M3D_slcpot, stat=nalloc)
  if (nalloc/=0) goto 108
  allocate(M3D_slcpot(ndimx,ndimy), stat=nalloc)
  if (nalloc/=0) goto 109
  M3D_slcpot = cmplx(0.0,0.0)
  
  ! PROJECTION
  ! - calculate fractional z coordinates of the slice start and end
  csz = sqrt( sum( M3D_b3*M3D_b3 ) ) * real( M3D_n3 ) ! cell size along the 3rd dimension [nm]
  fdz = 1.0 / real(nslc) ! fractional slice thickness
  dz  = fdz * csz ! real slice thickness [nm]
  fz0 = real(islc-1) * fdz ! fractional slice start z-coordinate
  fz1 = real(islc) * fdz   ! fractional slice stop  z-coordinate
  pscal = sigmae * dz ! projected potential phase action pre-factor
  fk0 = fz0 * real( M3D_n3 ) ! fractional 3D potential start z-coordinate
  fk1 = fz1 * real( M3D_n3 ) ! fractional 3D potential stop  z-coordinate
  k0 = floor(fk0)   ! integer z coordinate below the start coordinate
  k1 = ceiling(fk1) ! integer z coordinate above the start coordinate
  ! -----> M3D_slcpot <-----| SUM( M3D_pot(:,:,k+1), k=fk0...fk1 )
  !
  ! Example of the projection algorithm for the case of 8 3D potential planes
  ! and 3 projected slice potentials.
  !
  ! k=  8  1  2  3  4  5  6  7  8  1  2
  !        |                    |  |  
  ! ----+--+--+--+--+--+--+--+--+--+--+--> fz
  !    -1  0  1  2  3  4  5  6  7  8  9
  !        |       |       |       |
  ! fkx=   0.000   2.667   5.333   1.000
  !
  ! sl1=  (1)*1.000
  !         +(2)*1.000
  !            +(3)*0.667
  ! sl2=        (3)*0.333
  !               +(4)*1.000
  !                  +(5)*1.000
  !                     +(6)*0.333
  ! sl3=                 (6)*0.667
  !                        +(7)*1.000
  !                           +(8)*1.000
  !
  do k=k0, k1
    fki = max(0.0, min(1.0, fk1-real(k)))
    fkt = min(1.0, max(0.0, real(k+1)-fk0))
    fk  = min(fki,fkt) ! contribution of plane k to islc
    if (fk==0.0) cycle ! skip this plane
    k2 = 1 + modulo(k, M3D_n3)
    do j=1, M3D_n2
      do i=1, M3D_n1
        M3D_slcpot(i,j) = M3D_slcpot(i,j) + fk*cabsorb*M3D_pot(i,j,k2)
      end do
    end do
  end do
  
  ! PERIODIC REPEAT (nrx,nry)
  do j1=1, nry
    do i1=1, nrx
      if (i1==1.and.j1==1) cycle ! this part contains the data to repeat
      i2 = 1 + (i1-1)*M3D_n1
      i3 = i1*M3D_n1
      do j=1, M3D_n2
        j2 = j + (j-1)*M3D_n2  
        M3D_slcpot(i2:i3,j2) = M3D_slcpot(1:M3D_n1,j)
      end do
    end do
  end do
  
  ! PHASE-GRATING CALCULATION
  !
  ! save as phase grating
  !
  do j=1, ndimy
    do i=1, ndimx
      cpgr = exp( cimag*pscal*M3D_slcpot(i,j) ) ! phase grating = exp( I*V - ABF*V )
      pgr(i,j) = cpgr ! set phase grating value
    end do
  end do
  
  return

! ------ Error #1 : Invalid input parameter (number of slices).
101 nerr = 1
  call M3D_ErrorMessage("(M3D_getslice_pgr) "// &
     & "Invalid input parameter (number of slices).", nerr )
  return
! ------ Error #2 : Invalid input parameter (slice index).
102 nerr = 2
  call M3D_ErrorMessage("(M3D_getslice_pgr) "// &
     & "Invalid input parameter (slices index).", nerr )
  return
! ------ Error #3 : Invalid input parameter (absorption constant).
103 nerr = 3
  call M3D_ErrorMessage("(M3D_getslice_pgr) "// &
     & "Invalid input parameter (absorption constant).", nerr )
  return
! ------ Error #8 : Failed to deallocate an allocated array.
108 nerr = 8
  call M3D_ErrorMessage("(M3D_getslice_pgr) "// &
     & "Failed to deallocate an allocated array.", nerr )
  return
! ------ Error #9 : Failed to allocate array memory.
109 nerr = 9
  call M3D_ErrorMessage("(M3D_getslice_pgr) "// &
     & "Failed to allocate array memory.", nerr )
  return
  
end subroutine M3D_getslice_pgr
!
!**********************************************************************!



!**********************************************************************!
!
subroutine M3D_Message(smsg)
  implicit none
  character(len=*) :: smsg
  write(unit=M3D_out,fmt='(A)') trim(smsg)
  return
end subroutine M3D_Message
!
!**********************************************************************!



!**********************************************************************!
!
subroutine M3D_ErrorMessage(smsg,code)
  implicit none
  character(len=*), intent(in) :: smsg
  integer*4, intent(in) :: code
  write(unit=M3D_out,fmt='("ERROR ",A," (",I6,")")') trim(smsg),code
  return
end subroutine M3D_ErrorMessage
!
!**********************************************************************!



END MODULE m3dpot

!**********************************************************************!
!**********************************************************************!
