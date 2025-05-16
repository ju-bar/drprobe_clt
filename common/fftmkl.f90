!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  fftmkl.f90                                            !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2021        !
!                                                                      !
!**********************************************************************!
!                                                                      !
!    MODULE fftmkl                                                     !
!    -------------                                                     !
!                                                                      !
!    Purpose  : Wrapper to the Intel Math Kernel Library FFT routines  !
!    Version  : 1.0.0, Oct 22, 2021                                    !
!                                                                      !
!   Linked Libs: Intel Math Kernel Library (sequential)                !
!   Includes   : mkl_dfti.f90 (from mkl/include into source project)   !
!                MKL doesn't provide object code for this module, so   !
!                the source must be included in the source project     !
!                to generate the object code.                          !
!                                                                      !
!**********************************************************************!
!                                                                       
!  Author:  Juri Barthel                                
!           Ernst Ruska-Centre                                          
!           Forschungszentrum Jülich GmbH, 52425 Jülich, Germany        
!                                                                       
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!  
! REMARKS:
! - The transformations never apply re-normalization of the data.
!   With each call, there is a factor of sqrt(nx*ny) applied to each
!   coefficient.
!                                                                       
!-----------------------------------------------------------------------
  
  
subroutine fft_cc(cdata, n, dir)
! Calculates a one-dimensional complex-to-complex FFT by the Intel MKL
! in place of a given array cdata (complex*8) with length n.
! The parameter dir determines forward (>=0) or backward (<0) transform.
! single precision
  use mkl_dfti
  implicit none
  integer*4, intent(in) :: n, dir
  complex*8, intent(inout) :: cdata(n)   
  integer*4 :: ftstatus
  type(DFTI_DESCRIPTOR), POINTER :: descr
  if (n <= 0) return ! skip invalid length
  ftstatus = DftiCreateDescriptor(descr, DFTI_SINGLE,&
         & DFTI_COMPLEX, 1, n)
  ftstatus = DftiCommitDescriptor(descr)
  if (dir >= 0) then ! forward transform
    ftstatus = DftiComputeForward(descr, cdata)
  else ! backward transform
    ftstatus = DftiComputeBackward(descr, cdata)
  end if
  ftstatus = DftiFreeDescriptor(descr)
  return
end 

subroutine fft2_cc(cdata, nx, ny, dir)
! Calculates a two-dimensional complex-to-complex FFT by the Intel MKL
! in place of a given array cdata (complex*8) with dimensions nx, ny.
! The parameter dir determines forward (>=0) or backward (<0) transform.
! single precision
  use mkl_dfti
  implicit none
  integer*4, intent(in) :: nx, ny, dir
  complex*8, intent(inout) :: cdata(nx,ny)   
  integer*4 :: ftstatus, ftdims(2)
  type(DFTI_DESCRIPTOR), POINTER :: descr
  if (nx <= 0 .or. ny <= 0) return ! skip invalid dimensions
  ftdims = (/ nx, ny /)
  ftstatus = DftiCreateDescriptor(descr, DFTI_SINGLE,&
         & DFTI_COMPLEX, 2, ftdims)
  ftstatus = DftiCommitDescriptor(descr)
  if (dir >= 0) then ! forward transform
    ftstatus = DftiComputeForward(descr, cdata(:,1)) ! was interfacing cdata(:,1) in original code
  else ! backward transform
    ftstatus = DftiComputeBackward(descr, cdata(:,1))
  end if
  ftstatus = DftiFreeDescriptor(descr)
  return
end 
  
subroutine fft3_cc(cdata, nx, ny, nz, dir)
! Calculates a three-dimensional complex-to-complex FFT by the Intel MKL
! in place of a given array cdata (complex*8) with dimensions nx, ny, nz.
! The parameter dir determines forward (>=0) or backward (<0) transform.
! single precision
  use mkl_dfti
  implicit none
  integer*4, intent(in) :: nx, ny, nz, dir
  complex*8, intent(inout) :: cdata(nx,ny,nz)   
  integer*4 :: ftstatus, ftdims(3)
  type(DFTI_DESCRIPTOR), POINTER :: descr
  if (nx <= 0 .or. ny <= 0 .or. nz <=0) return ! skip invalid dimensions
  ftdims = (/ nx, ny, nz /)
  ftstatus = DftiCreateDescriptor(descr, DFTI_SINGLE,&
         & DFTI_COMPLEX, 3, ftdims)
  ftstatus = DftiCommitDescriptor(descr)
  if (dir >= 0) then ! forward transform
    ftstatus = DftiComputeForward(descr, cdata(:,1,1))
  else ! backward transform
    ftstatus = DftiComputeBackward(descr, cdata(:,1,1))
  end if
  ftstatus = DftiFreeDescriptor(descr)
  return
end 