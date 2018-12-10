
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Version 19.11.01, A.Thust * modified by J.Barthel 24.01.15                 C
C                                                                              C
C   Orginal code: FFTPACK                                                      C
C     P.N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations        C
C     (G. Rodrigue, ed.), Academic Press (1982) 51–83.                         C
C     https://doi.org/10.1016/B978-0-12-592101-5.50007-5                       C
C                                                                              C
C  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   C
C  C                                                                       C   C
C  C     This file contains machine independent 2-dimensional FFTs         C   C
C  C     in SINGLE PRECISION.                                              C   C
C  C                                                                       C   C
C  C     Also the appropriate scrambling/unscrambling and                  C   C
C  C     transposition routines are stored here.                           C   C
C  C                                                                       C   C
C  C     Optimum performance on DEC ALPHA and IBM RS6000 systems.          C   C
C  C                                                                       C   C
C  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   C
C                                                                              C
C                                                                              C
C  Watch out: The pure FFTs without the use of the appropriate                 C
C             scrambling/unscrambling and transposition routine                C
C             (one routine for each FFT) yields a scrambled                    C
C             and transposed result with x and y interchanged.                 C
C                                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                       
!  Authors: Andreas Thust & Juri Barthel                                
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




C General complex-to-complex variable x/y dimension FFTs (rectangular shape!)
C ODDCC128S  : Complex*8 to Complex*8 FFT, max. dimension 128, FOR/BACK
C ODDCC256S  : Complex*8 to Complex*8 FFT, max. dimension 256, FOR/BACK
C ODDCC512S  : Complex*8 to Complex*8 FFT, max. dimension 512, FOR/BACK
C ODDCC1024S  : Complex*8 to Complex*8 FFT, max. dimension 1024, FOR/BACK
C ODDCC2048S  : Complex*8 to Complex*8 FFT, max. dimension 2048, FOR/BACK
C ODDCC4096S  : Complex*8 to Complex*8 FFT, max. dimension 4096, FOR/BACK
C ODDCC8192S  : Complex*8 to Complex*8 FFT, max. dimension 8192, FOR/BACK

C General complex-to-complex variable 1-dimensional FFTs (arbitrary sub-length!)
C ODDCC128  : Complex*8 to Complex*8 FFT, max. dimension 128, FOR/BACK
C ODDCC256  : Complex*8 to Complex*8 FFT, max. dimension 256, FOR/BACK
C ODDCC512  : Complex*8 to Complex*8 FFT, max. dimension 512, FOR/BACK
C ODDCC1024  : Complex*8 to Complex*8 FFT, max. dimension 1024, FOR/BACK
C ODDCC2048  : Complex*8 to Complex*8 FFT, max. dimension 2048, FOR/BACK
C ODDCC4096  : Complex*8 to Complex*8 FFT, max. dimension 4096, FOR/BACK
C ODDCC8192  : Complex*8 to Complex*8 FFT, max. dimension 8192, FOR/BACK

C SCRODD128S: Scram./Transp. for FFT named ODDCC128S, max. dimension 128 
C SCRODD256S: Scram./Transp. for FFT named ODDCC256S, max. dimension 256 
C SCRODD512S: Scram./Transp. for FFT named ODDCC512S, max. dimension 512 
C SCRODD1024S: Scram./Transp. for FFT named ODDCC1024S, max. dimension 1024 
C SCRODD2048S: Scram./Transp. for FFT named ODDCC2048S, max. dimension 2048 
C SCRODD4096S: Scram./Transp. for FFT named ODDCC4096S, max. dimension 4096 
C SCRODD8192S: Scram./Transp. for FFT named ODDCC8192S, max. dimension 8192



C CALCULATE LARGEST PRIME FACTOR : ROUTINE "PRIME"
C FOR USE (OPTIONAL) WITH VARIABLE DIMENSION FFTS 









C***************************************************************************
C***
C*** 2D FFTs ***
C***
C***************************************************************************







C***************************************************************************
C********************** BEGIN subroutine ODDCC128S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 128 x 128  pixels            *
C    *****************************************************************            
C    *                        ODDCC128S                              *
C    *****************************************************************

    



      subroutine ODDCC128S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=128)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC128S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC128S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 128.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"
      
      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"

      
      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC128S *********************
C***************************************************************************








C***************************************************************************
C********************** BEGIN subroutine ODDCC256S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 256 x 256  pixels            *
C    *****************************************************************            
C    *                        ODDCC256S                              *
C    *****************************************************************

    



      subroutine ODDCC256S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=256)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC256S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC256S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 256.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"
      
      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"

      
      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC256S *********************
C***************************************************************************







C***************************************************************************
C********************** BEGIN subroutine ODDCC512S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 512 x 512  pixels            *
C    *****************************************************************            
C    *                        ODDCC512S                              *
C    *****************************************************************

    



      subroutine ODDCC512S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=512)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC512S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC512S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 512.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"

      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"


      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC512S *********************
C***************************************************************************






C***************************************************************************
C********************** BEGIN subroutine ODDCC1024S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 1024 x 1024  pixels          *
C    *****************************************************************            
C    *                        ODDCC1024S                             *
C    *****************************************************************

    



      subroutine ODDCC1024S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=1024)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC1024S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC1024S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 1024.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"

      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"


      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC1024S *********************
C***************************************************************************







C***************************************************************************
C********************** BEGIN subroutine ODDCC2048S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 2048 x 2048  pixels          *
C    *****************************************************************            
C    *                        ODDCC2048S                             *
C    *****************************************************************

    



      subroutine ODDCC2048S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=2048)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC2048S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC2048S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 2048.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"

      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"


      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC2048S *********************
C***************************************************************************






C***************************************************************************
C********************** BEGIN subroutine ODDCC4096S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 4096 x 4096  pixels          *
C    *****************************************************************            
C    *                        ODDCC4096S                             *
C    *****************************************************************

    



      subroutine ODDCC4096S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=4096)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC4096S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC4096S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 4096.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"

      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"


      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC4096S *********************
C***************************************************************************






C***************************************************************************
C********************** BEGIN subroutine ODDCC8192S ************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY image dimensions *
C    *       For complex*8 arrays up to 8192 x 8192  pixels          *
C    *****************************************************************            
C    *                        ODDCC8192S                             *
C    *****************************************************************

    



      subroutine ODDCC8192S(A,nnx,nny,direction)



C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)


      PARAMETER (ndim=8192)


      character*40 direction
      character*3  dir

     
      complex*8 A(ndim,ndim)   
      real*4    WORK(4*ndim+15)
      complex*8 temp


C******************* Initialize Parameters ********************************


      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC8192S.'
      STOP
      endif



C     Exchange x- and y-dimension in case of backtransform
C     (Before a backtransform the array is scrambled and transposed)
C     The exchange of dimensions is just for a more convenient
C     handling, since then the x- and y- dimension stay the same
C     in the calling routine, regardless if forward or backward
C     transform. This is inverse to the actual physical array size 
C     in case of a backward transform of a scrambled/transposed array !!! 
C     The nomenclature addresses therefore always the unscrambled and
C     untransposed state of the array.
C     The dimensions "nnx" and "nny" of the input parameter list
C     should NEVER be changed by this routine, otherwise external
C     routines called afterwards see exchanged parameters "nnx" and
C     "nny". Use internally "nx" and "ny" instead.



      if(isign .eq. 1)then
      nx = nnx
      ny = nny
      scale =  1./(float(nx) * float(ny))
      else
      ntemp = nny
      ny = nnx
      nx = ntemp
      scale = 1.
      endif




C     Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC8192S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 8192.'
      write(unit=6,fmt=*)
      STOP
      endif




C     Find maximum of the x- and y- dimensions

      if(nx .ge. ny)then
      nxy = nx
      else
      nxy = ny
      endif


C     Initialize work array "WORK"

      call CCFFTI(nx,WORK)



C********************** Do FFT ********************************************


C--------------------------------------------------------------------------
C---> First do sucessive 1D FFTs on the ROWS of "A"  (do 10 loop)


      do 10 j=1,ny

C     Transform j-th row of 2-dim array "A"

      if (isign .eq. 1)then
      call CCFFTF(nx,A(1,j),WORK)
      else
      call CCFFTB(nx,A(1,j),WORK)
      endif

10    continue   
C--------------------------------------------------------------------------




C--------------------------------------------------------------------------
C---> Transpose 2-dim array "A" 


      do 20 j=1,nxy
      do 30 i=1,j
      temp   = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
30    continue
20    continue
C--------------------------------------------------------------------------




C---> Re-initialize work array "WORK"


      call CCFFTI(ny,WORK)



C---------------------------------------------------------------------------
C---> Repeat sequence of 1-dim row FFTs on transposed array "A" (do 40 loop)

      do 40 j=1,nx

C     Transform j-th row of 2-dim array "A"

      if(isign .eq. 1)then
      call CCFFTF(ny,A(1,j),WORK)
      else
      call CCFFTB(ny,A(1,j),WORK)
      endif

C     Scale j-th row of 2-dim array "A" in case of forward FFT

      if(isign .eq. 1)then
      do 50 i=1,ny
      A(i,j) = A(i,j)*scale
50    continue
      endif

40    continue   
C---------------------------------------------------------------------------



      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC8192S *********************
C***************************************************************************










C***************************************************************************
C***
C*** 1D FFTs ***
C***
C***************************************************************************




C***************************************************************************
C********************** BEGIN subroutine ODDCC128 **************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 128 pixels                      *
C    *****************************************************************            
C    *                        ODDCC128                               *
C    *****************************************************************


      subroutine ODDCC128(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=128)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC128.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC128.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 128.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC128 ***********************
C***************************************************************************






C***************************************************************************
C********************** BEGIN subroutine ODDCC256 **************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 256 pixels                      *
C    *****************************************************************            
C    *                        ODDCC256                               *
C    *****************************************************************


      subroutine ODDCC256(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=256)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC256.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC256.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 256.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC256 ***********************
C***************************************************************************




C***************************************************************************
C********************** BEGIN subroutine ODDCC512 **************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 512 pixels                      *
C    *****************************************************************            
C    *                        ODDCC512                               *
C    *****************************************************************


      subroutine ODDCC512(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=512)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC512.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC512.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 512.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC512 ***********************
C***************************************************************************



C***************************************************************************
C********************** BEGIN subroutine ODDCC1024 *************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 1024 pixels                     *
C    *****************************************************************            
C    *                        ODDCC1024                              *
C    *****************************************************************


      subroutine ODDCC1024(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=1024)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC1024.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC1024.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 1024.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC1024 **********************
C***************************************************************************


C***************************************************************************
C********************** BEGIN subroutine ODDCC2048 *************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 2048 pixels                     *
C    *****************************************************************            
C    *                        ODDCC2048                              *
C    *****************************************************************


      subroutine ODDCC2048(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=2048)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC2048.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC2048.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 2048.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC2048 **********************
C***************************************************************************




C***************************************************************************
C********************** BEGIN subroutine ODDCC4096 *************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 4096 pixels                     *
C    *****************************************************************            
C    *                        ODDCC4096                              *
C    *****************************************************************


      subroutine ODDCC4096(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=4096)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC4096.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC4096.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 4096.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC4096 **********************
C***************************************************************************



C***************************************************************************
C********************** BEGIN subroutine ODDCC8192 *************************
C***************************************************************************


C    *****************************************************************
C    * General complex to complex FFT for ARBITRARY 1-dim arrys      *
C    *       For complex*8 arrays up 8192 pixels                     *
C    *****************************************************************            
C    *                        ODDCC8192                              *
C    *****************************************************************


      subroutine ODDCC8192(A,nn,direction)


C********************** Declarations **************************************

      IMPLICIT REAL*4 (A-H,O-Z)

      PARAMETER (ndim=8192)

      character*40 direction
      character*3  dir

      complex*8 A(ndim)   
      real*4    WORK(4*ndim+15)

C******************* Initialize Parameters ********************************

      dir = direction(1:3)

      if( dir .eq. 'for' .or. dir .eq. 'FOR') then
      isign =  1
      else if ( dir .eq. 'bac' .or. dir .eq. 'BAC') then
      isign = -1
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of ODDCC8192.'
      STOP
      endif

      n = nn
      if(isign .eq. 1)then
      scale =  1./(float(n))
      else
      scale = 1.
      endif

C     Check if given x- and y-dimensions "n" does not exceed
C     array boundaries defined by parameter "ndim".

      if(n .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine ODDCC8192.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 8192.'
      write(unit=6,fmt=*)
      STOP
      endif

C     Initialize work array "WORK"

      call CCFFTI(n,WORK)

C********************** Do FFT ********************************************

C--------------------------------------------------------------------------
C     Transform array "A"

      if (isign .eq. 1)then
      call CCFFTF(n,A,WORK)
      A = A * scale
      else
      call CCFFTB(n,A,WORK)
      endif
C--------------------------------------------------------------------------

      return
      end 


C***************************************************************************
C************************ END of subroutine ODDCC8192 **********************
C***************************************************************************




C***************************************************************************
C********* END GENERAL COMPLEX-TO-COMPLEX ARBITRARY SIZE FFTs **************
C***************************************************************************













C******************************************************************************
C******************************************************************************
C********* 1-dim. FFT required by all FFT routines listed above ***************
C*********        Algorithm according to FFTPACK suite          ***************
C******************************************************************************
C******************************************************************************

C  "CCFFTf"   |
C  "CCFFTf1"  |
C  "PISSF"    |
C  "PISSF2"   |    Routines for 1-dim forward transform
C  "PISSF3"   |
C  "PISSF4"   |
C  "PISSF5"   |


C  "CCFFTb"    |
C  "CCFFTb1"   |
C  "PISSB"    |
C  "PISSB2"   |    Routines for 1-dim backward transform
C  "PISSB3"   |
C  "PISSB4"   |
C  "PISSB5"   |


C  "CCCFFTi"    |    Initialization routines
C  "CCFFTi1"   |

C-------------------------------------------------------------------------------


      SUBROUTINE CCFFTF (N,C,WSAVE)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       C(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CCFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

C-------------------------------------------------------------------------------


      SUBROUTINE CCFFTF1 (N,C,CH,WA,IFAC)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)
      murks2 = 2
      NF = IFAC(murks2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PISSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PISSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PISSF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PISSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PISSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PISSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PISSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PISSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PISSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PISSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSF2 (IDO,L1,CC,CH,WA1)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(1)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSF3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(1)     ,WA2(1)
      DATA TAUR,TAUI /-0.5,-.866025403784439/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------


      SUBROUTINE PISSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,
     1   -.951056516295154,-.809016994374947,
     1   -.587785252292473/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------


      SUBROUTINE CCFFTB (N,C,WSAVE)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       C(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CCFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE CCFFTB1 (N,C,CH,WA,IFAC)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)
      murks2 = 2
      NF = IFAC(murks2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PISSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PISSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PISSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PISSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PISSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PISSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PISSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PISSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PISSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PISSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END


C-------------------------------------------------------------------------------

      SUBROUTINE PISSB2 (IDO,L1,CC,CH,WA1)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(1)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C-------------------------------------------------------------------------------

      SUBROUTINE PISSB3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(1)     ,WA2(1)
      DATA TAUR,TAUI /-0.5,0.866025403784439/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE PISSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,
     1   .951056516295154, -.809016994374947,
     1   .587785252292473/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE CCFFTI (N,WSAVE)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
C     WSAVE(1:4*N) = 0
      CALL CCFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

C-------------------------------------------------------------------------------

      SUBROUTINE CCFFTI1 (N,WA,IFAC)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION       WA(1)      ,IFAC(1)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      murks3  = 3
      IFAC(murks3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      murks2  = 2
      IFAC(murks2) = NF
      TPI = 6.28318530717959
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD+L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
C-------------------------------------------------------------------------------



C*******************************************************************************
C**************************** END 1-D FFTs *************************************
C*******************************************************************************






C*******************************************************************************
C*******************************************************************************
C*******************************************************************************
C*******************************************************************************



C*******************************************************************************
C** BEGIN arbitrary-dimension complex*8 SCRAMBLING + TRANSPOSITION ROUTINES ****
C*******************************************************************************



C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD128S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC128S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  128 * 128.
C     ANY arbitrary array size with nx,ny smaller or equal
C     128 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD128S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=128)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD128S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD128S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 128.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD128S **********************************
C*******************************************************************************







C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD256S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC256S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  256 * 256.
C     ANY arbitrary array size with nx,ny smaller or equal
C     256 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD256S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=256)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD256S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD256S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 256.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD256S **********************************
C*******************************************************************************






C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD512S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC512S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  512 * 512.
C     ANY arbitrary array size with nx,ny smaller or equal
C     512 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD512S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=512)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD512S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD512S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 512.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD512S **********************************
C*******************************************************************************








C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD1024S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC1024S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  1024 * 1024.
C     ANY arbitrary array size with nx,ny smaller or equal
C     1024 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD1024S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=1024)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD1024S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD1024S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 1024.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD1024S **********************************
C*******************************************************************************







C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD2048S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC2048S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  2048 * 2048.
C     ANY arbitrary array size with nx,ny smaller or equal
C     2048 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD2048S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=2048)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD2048S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD2048S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 2048.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD2048S *********************************
C*******************************************************************************





C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD4096S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC4096S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  4096 * 4096.
C     ANY arbitrary array size with nx,ny smaller or equal
C     4096 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD4096S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=4096)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD4096S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD4096S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 4096.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD4096S *********************************
C*******************************************************************************





C*******************************************************************************
C*******************************************************************************
C
C           ***************************************************
C           *                                                 *
C           * subroutine SCRODD8192S(SC,USC,nx,ny,direction)  *
C           *                                                 *
C           ***************************************************
C
C
C
C     This routine scrambles/unscrambles and TRANSPOSES 
C     2D-FFT arbitrary rectangular arrays for/from routine 
C     "ODDCC8192S".
C     Two things have to be done:
C
C     o Transposition
C     o conventional scrambling/unscrambling
C
C     This version for arbitrary "nx" * "ny" sized arrays.
C     Maximum array size is  8192 * 8192.
C     ANY arbitrary array size with nx,ny smaller or equal
C     8192 can be handled. 
C
C     SC   :     Storage for scrambled array
C     USC  :     Storage for unscrambled array
C     DIRECTION: If direction = "scramble" scrambling will take place
C                If direction = "unscramble" unscrambling will take place
C
C*******************************************************************************
C*******************************************************************************




      subroutine SCRODD8192S(SC,USC,nx,ny,direction) 



      PARAMETER(ndim=8192)




      parameter  (ndim2=ndim/2)    
      complex*8  SC(ndim,ndim)
      complex*8  USC(-ndim2:ndim2-1,-ndim2:ndim2-1)
      character*40 direction
      character*3 dir
      logical scram

C*******************************************************************************

      dir=direction(1:3)

      if( dir .eq. 'scr' .or. dir .eq. 'SCR')then
      scram = .true.
      else if ( dir .eq. 'uns' .or. dir .eq. 'UNS') then
      scram = .false.
      else
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)' ERROR in the parameter list of SCRODD8192S.'
      stop
      endif


C*******************************************************************************




C---> Check if given x- and y-dimensions "nx" and "ny" do not exceed
C     array boundaries defined by parameter "ndim".

      if(nx .gt. ndim .or. ny .gt. ndim)then
      write(unit=6,fmt=*)
      write(unit=6,fmt=*)'   ERROR message from routine SCRODD8192S.'
      write(unit=6,fmt=*)'   Maximum allowed array dimension is 8192.'
      write(unit=6,fmt=*)
      STOP
      endif


C********************************************************************************


C---> Prepare array dimensions depending on if nx,ny are even or odd


      if(mod(nx,2) .eq. 0)then
      nxlo  = -nx/2
      nxhi  =  nx/2 - 1
      nxha  =  nx/2
      else
      nxlo  = -(nx-1)/2
      nxhi  =  (nx-1)/2
      nxha  =  (nx+1)/2  
      endif



      if(mod(ny,2) .eq. 0)then
      nylo  = -ny/2
      nyhi  =  ny/2 - 1
      nyha  =  ny/2
      else
      nylo  = -(ny-1)/2
      nyhi  =  (ny-1)/2
      nyha  =  (ny+1)/2
      endif


C*****************************************************************************


      if(scram)then


C     Scramble + Transpose


      do 300 kj= nylo,nyhi

      if (kj .ge. 0)then
      j=kj+1
      else
      j=kj+1+ny
      endif

      do 400 ki=nxlo,nxhi

      if (ki .ge. 0)then
      i=ki+1
      else
      i=ki+1+nx
      endif


      SC(j,i)=USC(ki,kj)

400   continue
300   continue



      endif


C*******************************************************************************

      if(.not. scram)then
     

C     Unscramble + Transpose

      do 30 j=1,nx

      if (j .le. nxha)then
      kj=j-1
      else
      kj=j-1-nx
      endif


      do 40 i=1,ny

      if (i .le. nyha)then
      ki=i-1
      else
      ki=i-1-ny
      endif

      USC(kj,ki)=SC(i,j)

40    continue
30    continue




      endif

C*******************************************************************************


      return
      end



C*******************************************************************************
C****************** END SUBROUTINE SCRODD8192S *********************************
C*******************************************************************************





C*******************************************************************************
C** END arbitrary-dimension complex*8 SCRAMBLING + TRANSPOSITION ROUTINES ******
C*******************************************************************************





*******************************************************************************
************************* BEGIN ROUTINE PRIME *********************************
*******************************************************************************
      subroutine prime(num,n)
      
C     PURPOSE: calculate largest prime factor "n" for a number "num".
C     For use to predict calculation time of variable dimension FFT
C     and to warn user, if array dimesnion "num" is a large prime.
     

      n = num
      
      
5     continue      
                  
      
      do 10 i=2,n
            
      frac=float(n)/float(i)
      diff = abs(frac - float(nint(frac)))
      
      if(diff .lt. 1.e-6 .and. nint(frac) .gt. 1)then
      n = nint(frac)
      if(n .ne. num .and. n .ne. 1)then
      goto 5
      endif
      endif
            
10    continue
     
      
      
      
      end
*******************************************************************************
************************* END ROUTINE PRIME ***********************************
******************************************************************************* 
