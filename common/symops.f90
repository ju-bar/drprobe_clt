MODULE symops
!
!
!**********************************************************************************
!*  SYMOPS                                                                        *
!**********************************************************************************
!* This module handles symmetry operations of linear type.                        *
!* A symmetry operation transforms a 3D position P into a new position P' with    *
!*     P' = M.P + S,                                                              *
!* where M is a 3x3 matrix and S is a 3D shift vector.                            *
!* Internally a symmetry operation is defined by a 12*REAL(dp) array, where       *
!* the items 1-3 define S and the items 4-12 define M, with 4-6 for the x',       *
!* 7-9 for the y', and 10-12 for the z' component.                                *
!**********************************************************************************
!* Linkage: - routines for messaging (PostMessage, PostWarning, CriticalError)    *
!**********************************************************************************
!* Last modification: J. Barthel - 07 Jan. 2016                                   *
!*                    - migrated from ATOMSK                                      *
!*                    - added tranformation functions from ATOMSK lib             *
!*                    - rearranged data structures by combining P and AUX         *
!*                      P contains all atomic site data:                          *
!*                        1:3  = coordinates                                      *
!*                        4    = atomic number                                    *
!*                        5    = charge                                           *
!*                        6    = occupancy                                        *
!*                        7:14 = thermal vibration parameter(s)                   *
!*                      P is transposed, compared to the ATOMSK code              *
!*                      each row contains data of one atomic site.                *
!**********************************************************************************
!* Symmetry operation handling and parsing, added by J. Barthel, July 2015        *
!*     SUBROUTINE SYMOPS_INIT initializes the array symops_trf to identity        *
!*                operartions. Since we leave the allocation of the array to      *
!*                other routines, and keep the array for public access accross    *
!*                the program, this routine is a pure initialization routine.     *
!*                Note:  The allocation state of symops_trf will be checked,      *
!*                but not altered. No error message will occur if not allocated.  *
!*     SYMOPS_TRF_UANI transforms a second order coefficient matrix like the      *
!*                anisotropic thermal displacements for a given transformation    *
!*                of real-space coordinates in the basis of the super-cell.       *
!*                See http://ww1.iucr.org/comm/cnom/adp/finrep/node8.html,        *
!*                equation (2.1.33), representing U' = A U A^T.                   *
!*     SUBROUTINE SYMOPS_APPLY applies the currently defined symmetry operations  *
!*                from the array "symops_trf" to the passed array P.              *
!*                P(:,1:3) are assumed to be fractional atom coordinates.         *
!*                Note: P may be re-allocated at the end of this routine          *
!*                since the number of atoms may have increased. Update size       *
!*                variables in the calling routine after calling SYMOPS_APPLY!    *
!*     SUBROUTINE SYMOPS_CHECK_STR checks if a string defines a symmetry          *
!*                operation.                                                      *
!*     SUBROUTINE SYMOPS_SET_STR translates a string symmetry operation to        *
!*                numbers defining linear transformations and stores these        *
!*                numbers in the module symops.                                   *
!*     SUBROUTINE SYMOPS_PARSE_STR_LINTRF translates a 1D operation string to     *
!*                linear transformation parameters by parsing the string.         *
!*     SUBROUTINE SYMOPS_SET_SGNAME sets symmetry operations for a spacegroup     *
!*                given by the space group Hermann-Mauguin symbol.                *
!*     SUBROUTINE SYMOPS_SET_SGNUM sets symmetry operations for a spacegroup      *
!*                given by the space group number (1 - 230).                      *
!* REMARK: SUBROUTINE SYMOPS_PARSE_STR_LINTRF is a recursive routine.             *
!*                Use the respective compiler option to enable recursice routines *
!*                (ifort: /recursive)                                             *
!**********************************************************************************
!* This program is free software: you can redistribute it and/or modify           *
!* it under the terms of the GNU General Public License as published by           *
!* the Free Software Foundation, either version 3 of the License, or              *
!* (at your option) any later version.                                            *
!*                                                                                *
!* This program is distributed in the hope that it will be useful,                *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
!* GNU General Public License for more details.                                   *
!*                                                                                *
!* You should have received a copy of the GNU General Public License              *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
!**********************************************************************************
!
USE spacegroups
!
IMPLICIT NONE
!
!
! MODULE DATA
!
INTEGER,PARAMETER,PRIVATE:: dp = SELECTED_REAL_KIND(15,307)  !reals with 64-bits precision
INTEGER,PUBLIC:: symops_verbosity ! verbosity level
DATA symops_verbosity /0/ ! default verbosity level is root level (almost nothing)
INTEGER,PUBLIC:: symops_nerr, symops_warn
DATA symops_nerr /0/ ! number of errors
DATA symops_warn /0/ ! number of warnings
!
CHARACTER(LEN=3),PARAMETER,PUBLIC:: symops_chanstr = 'xyz'
INTEGER,PARAMETER,PUBLIC:: symops_nltrf = 12
REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: symops_trf ! symmetry operations list
!
!
!
! MODULE ROUTINES
!
CONTAINS
!
!
SUBROUTINE SYMOPS_INIT()
!
INTEGER:: n, i, j
!
! Initialization
IF(.NOT.ALLOCATED(symops_trf)) RETURN ! nothing to do, exit
j=SIZE(symops_trf,1) ! Get number of row entries
n=SIZE(symops_trf,2) ! Get number of transformations
!
! Check number of row entries for consistency:
IF(j.NE.symops_nltrf) RETURN ! inconsistent row size, exit
                             ! Leave initialization to external routine.
IF(n<=0) RETURN ! Invalid number of rows, exit
                ! Leave initialization to external routine.
symops_trf = 0.d0 ! set all entries to zero
DO i=1, n ! loop through rows and set M to Identity
  symops_trf( 4,i) = 1.d0
  symops_trf( 8,i) = 1.d0
  symops_trf(12,i) = 1.d0
ENDDO
!
! reset error and warning counts
symops_nerr = 0
symops_warn = 0
!
END SUBROUTINE SYMOPS_INIT
!
!
!
!
!********************************************************
! SYMOPS_MSG
! This subroutine posts a message depending on the
! module verbosity level.
! Currently this is a wrapper function calling a program
! specific messaging routine.
!********************************************************
SUBROUTINE SYMOPS_MSG(smsg,nlev)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: smsg
  INTEGER,OPTIONAL :: nlev
  INTEGER :: ulev
  ulev = 0 ! default is to post base level messages (almost nothing)
  IF(PRESENT(nlev)) ulev=nlev
  IF(ulev<=symops_verbosity) THEN
  !
    call PostMessage(trim(smsg))
  !
  ENDIF
  !
END SUBROUTINE SYMOPS_MSG
!
!
!********************************************************
! SYMOPS_WRN
! This subroutine posts a warning.
! Currently this is a wrapper function calling a program
! specific messaging routine.
!********************************************************
SUBROUTINE SYMOPS_WRN(smsg)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: smsg
  !
  call PostWarning(trim(smsg))
  symops_warn = symops_warn+1
  !
END SUBROUTINE SYMOPS_WRN
!
!
!
!
!********************************************************
! SYMOPS_ERR
! This subroutine posts an error message.
!********************************************************
SUBROUTINE SYMOPS_ERR(smsg,error_code)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: smsg
  INTEGER, OPTIONAL :: error_code
  CHARACTER(LEN=20) :: scode
  !
  scode=""
  IF( PRESENT(error_code) ) THEN
    write (unit=scode,fmt=*) error_code
    write (unit=6,fmt='(A)') "ERROR: "//trim(smsg)// &
        & " ("//trim(adjustl(scode))//")"
  ELSE
    write (unit=6,fmt='(A)') "ERROR: "//trim(smsg)
  ENDIF
  symops_nerr = symops_nerr+1
  !
END SUBROUTINE SYMOPS_ERR
!
!
!
!
!
!
!!********************************************************
!! INVMAT
!! This subroutine inverts a NxN matrix M
!! and outputs the result into the matrix G.
!! requires LAPACK
!!********************************************************
!SUBROUTINE INVMAT(M,G,status)
!!
!IMPLICIT NONE
!REAL(dp),DIMENSION(:,:),INTENT(IN):: M
!REAL(dp),DIMENSION(:,:),INTENT(OUT):: G
!INTEGER,INTENT(OUT),OPTIONAL:: status
!INTEGER:: i
!INTEGER:: LWORK !for LAPACK routine DGETRI
!INTEGER,DIMENSION(SIZE(M,1)):: IPIV !for LAPACK routine DGETRI
!REAL(dp):: det
!REAL(dp),DIMENSION(SIZE(M,1)):: WORK !for LAPACK routine DGETRI
!!
!i=0
!!
!IF( SIZE(M,1).NE.SIZE(M,2) ) THEN
!  !Non-square matrix: cancel
!  i=1
!  !
!ELSE
!  IF( SIZE(M,1)==3 .AND. SIZE(M,2)==3 ) THEN
!    !3x3 matrix: simple enough, let's do it by hand
!    det =   M(1,1)*M(2,2)*M(3,3) - M(1,1)*M(3,2)*M(2,3) &
!        & - M(2,1)*M(1,2)*M(3,3) + M(2,1)*M(3,2)*M(1,3) &
!        & + M(3,1)*M(1,2)*M(2,3) - M(3,1)*M(2,2)*M(1,3)
!    !
!    G(1,1) = (M(2,2)*M(3,3) - M(2,3)*M(3,2))/det
!    G(2,1) = (M(2,3)*M(3,1) - M(2,1)*M(3,3))/det
!    G(3,1) = (M(2,1)*M(3,2) - M(2,2)*M(3,1))/det
!    !
!    G(1,2) = (M(3,2)*M(1,3) - M(3,3)*M(1,2))/det
!    G(2,2) = (M(3,3)*M(1,1) - M(3,1)*M(1,3))/det
!    G(3,2) = (M(3,1)*M(1,2) - M(3,2)*M(1,1))/det
!    !
!    G(1,3) = (M(1,2)*M(2,3) - M(1,3)*M(2,2))/det
!    G(2,3) = (M(1,3)*M(2,1) - M(1,1)*M(2,3))/det
!    G(3,3) = (M(1,1)*M(2,2) - M(1,2)*M(2,1))/det
!    !
!  ELSEIF( SIZE(G,1)==SIZE(M,1) .AND. SIZE(G,2)==SIZE(M,2) ) THEN
!    !general NxN matrix: call LAPACK routines DGETRF and DGETRI
!    G(:,:) = M(:,:)
!    CALL DGETRF( SIZE(M,1), SIZE(M,2), G, SIZE(M,1), IPIV, i)
!    IF( i==0 ) THEN
!      CALL DGETRI( SIZE(M,1), G, SIZE(M,1), IPIV, WORK, SIZE(M,1), i )
!    ENDIF
!    !
!  ELSE
!    !non-consistent array sizes: some programmer
!    CALL SYMOPS_ERR("Matrix inversion failed, incosistent array sizes.")
!    i=1
!  ENDIF
!  !
!ENDIF
!!
!IF(PRESENT(status)) status=i
!!
!END SUBROUTINE INVMAT
!**********************************************************************!
! INVMAT
! This contains the Gauss-Jordan Algorithm copied from Numerical Recipes
! The algorithm applies full pivoting, but requires NxN matrices.
! The input matrix A is inverted and the result is stored in B.
! The status of inversion is written to status.
!   0: success,   >0: error code
! - modified by to work on arbitrary size NxN
!  (2015-12-16: J. Barthel, ju.barthel@fz-juelich.de)
!**********************************************************************!
SUBROUTINE INVMAT(A,B,status)
  implicit none
  REAL(dp),DIMENSION(:,:),INTENT(IN):: A
  REAL(dp),DIMENSION(:,:),INTENT(OUT):: B
  INTEGER,INTENT(OUT),OPTIONAL:: status
  INTEGER:: istatus, nalloc
  INTEGER:: I, ICOL, IROW, J, K, L, LL, N
  INTEGER,DIMENSION(:),ALLOCATABLE:: IPIV,INDXR,INDXC
  REAL(dp):: BIG, DUM, PIVINV
  !
  istatus = 0 ! preset internal status number
  nalloc = 0
  N = SIZE(A,1) ! get input array number of columns
  if (N<1) then
    istatus = 1
    goto 99
  end if
  if (N/=SIZE(A,2).or.N/=SIZE(B,1).or.N/=SIZE(B,2)) then ! check that all matrices have the same NxN size
    istatus = 2
    goto 99
  end if
  ! inversions
  select case(N)
  case (1) ! scalar
    if (A(1,1)/=0.0d+0) then
      B(1,1) = 1.0d+0/A(1,1)
    else
      B = 0.0d+0
      call SYMOPS_ERR("Matrix inversion failed, singular matrix (DET=0).", 3)
      istatus = 3
    end if
  case (2) ! 2x2 matrix
    DUM = A(1,1)*A(2,2)-A(1,2)*A(2,1)
    if (DUM/=0.d+0) then
      PIVINV = 1.0d+0/DUM
      B(1,1) = A(2,2)*PIVINV
      B(2,1) = -A(2,1)*PIVINV
      B(1,2) = -A(1,2)*PIVINV
      B(2,2) = A(1,1)*PIVINV
    else
      B = 0.0d+0
      call SYMOPS_ERR("Matrix inversion failed, singular matrix (DET=0).", 4)
      istatus = 4
    end if
  case (3) ! 3x3 matrix
    DUM =    A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) &
        &  - A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) &
        &  + A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
    if (DUM/=0.d+0) then
      PIVINV = 1.0d+0/DUM
      B(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2))*PIVINV
      B(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))*PIVINV
      B(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))*PIVINV
      B(1,2) = (A(3,2)*A(1,3) - A(3,3)*A(1,2))*PIVINV
      B(2,2) = (A(3,3)*A(1,1) - A(3,1)*A(1,3))*PIVINV
      B(3,2) = (A(3,1)*A(1,2) - A(3,2)*A(1,1))*PIVINV
      B(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2))*PIVINV
      B(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3))*PIVINV
      B(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))*PIVINV
    else
      B = 0.0d+0
      call SYMOPS_ERR("Matrix inversion failed, singular matrix (DET=0).", 5)
      istatus = 5
    end if
  case default ! general NxN case, solved by Gauss-Jordan
    B=A ! copy A to B, do all transformations on B only, leaving A untouched
    ALLOCATE(IPIV(N),INDXR(N),INDXC(N),STAT=nalloc) ! allocate bookeeping arrays for pivoting
    do J=1,N
      IPIV(J)=0
    end do
    do I=1,N ! main loop over columns to be reduced
      BIG=0.d+0
      do J=1,N ! outer loop of the search for a pivot elements
        if (IPIV(J).NE.1) then
          do K=1,N
            if (IPIV(K).EQ.0) then
              if (DABS(B(J,K)).GE.BIG) then
                BIG=DABS(B(J,K))
                IROW=J
                ICOL=K
              end if
            else if (IPIV(K).GT.1) then
              call SYMOPS_ERR("Matrix inversion failed, singular matrix (IPIV).", 6)
              istatus = 6
              goto 99
            end if
          end do
        end if
      end do
      IPIV(ICOL)=IPIV(ICOL)+1
      ! We now have the pivot element, so we interchange rows, if needed, to put the pivot
      ! element on the diagonal. The columns are not physically interchanged, only relabeled
      ! indxc(i): the column of the ith pivot element, is the ith column that is reduced
      ! indxr(i): the row in which that pivot element was originally located.
      ! If indxr(i) /= indxc(i) there is an implied column interchange.
      ! With this form of bookkeeping we keep track of the scambling of the inverse matrix.
      if (IROW.NE.ICOL) then
        do L=1,N
          DUM=B(IROW,L)
          B(IROW,L)=B(ICOL,L)
          B(ICOL,L)=DUM
        end do
      end if
      INDXR(I)=IROW ! We are now ready to divide the pivot row by the pivot
      INDXC(I)=ICOL !    element, located at irow and icol.
      if (B(ICOL,ICOL).EQ.0.d+0) then
        call SYMOPS_ERR("Matrix inversion failed, singular matrix (DIAG).", 7)
        istatus = 7
        goto 99
      end if
      PIVINV=1.d+0/B(ICOL,ICOL)
      B(ICOL,ICOL)=1.D+0
      do L=1,N
        B(ICOL,L)=B(ICOL,L)*PIVINV
      end do
      do LL=1,N ! Next, we reduce the rows except for the pivot one, of course.
        if (LL.NE.ICOL) then
          DUM=B(LL,ICOL)
          B(LL,ICOL)=0.d+0
          do L=1,N
            B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
          end do
        end if
      end do ! outer loop of the search for a pivot elements
    end do ! main loop over columns to be reduced
    ! Done.
    ! It only remains to unscramble the solution in view of the column interchanges.
    ! We do this by interchanging pairs of columns in the reverse order that the
    ! permutation was built up.
    do L=N,1,-1
      if (INDXR(L).NE.INDXC(L)) then
        do K=1,N
          DUM=A(K,INDXR(L))
          B(K,INDXR(L))=B(K,INDXC(L))
          B(K,INDXC(L))=DUM
        end do
      end if
    end do
    DEALLOCATE(IPIV,INDXR,INDXC,STAT=nalloc)
  end select ! case(n)
!
99 if(present(status)) status=istatus ! finish by transfering optional status
!
END SUBROUTINE INVMAT
!**********************************************************************!




!
!
!********************************************************
!  CONVMAT
!  This subroutine converts conventional vectors
!  defined by a b c alpha beta gamma,
!  into a lower triangular matrix:
!        |  H(1,1)    0       0     |
!   H =  |  H(2,1)  H(2,2)    0     |
!        |  H(3,1)  H(3,2)  H(3,3)  |
!  The array H contains the basis vectors in its columns.
!********************************************************
!
SUBROUTINE CONVMAT(a,b,c,alpha,beta,gamma,H)
!
IMPLICIT NONE
REAL(dp),PARAMETER:: d2r = 0.0174532925199
REAL(dp),INTENT(IN):: a, b, c, alpha, beta, gamma
REAL(dp),DIMENSION(3,3):: H
REAL(dp):: al, be, ga
!
al = alpha*d2r
be = beta *d2r
ga = gamma*d2r
H(:,:) = 0.d0
H(1,1) = a
H(2,1) = b*DCOS(ga)
H(2,2) = b*DSIN(ga)
H(3,1) = c*DCOS(be)
H(3,2) = c*(DCOS(al)-DCOS(be)*DCOS(ga))/DSIN(ga)
H(3,3) = c*DSQRT(                                       &
       &     DSIN(ga)**2 - DCOS(be)**2 - DCOS(al)**2    &
       &     +2.d0*DCOS(al)*DCOS(be)*DCOS(ga)           &
       &        ) / DSIN(ga)
!
!
END SUBROUTINE CONVMAT
!
!
!********************************************************
!  MATCONV
!  This subroutine converts a matrix into conventional
!  vectors defined by a b c alpha beta gamma.
!********************************************************
!
SUBROUTINE MATCONV(H,a,b,c,alpha,beta,gamma)
!
IMPLICIT NONE
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
!
 a = VECLENGTH(H(1,:))
 b = VECLENGTH(H(2,:))
 c = VECLENGTH(H(3,:))
 alpha = ANGVEC( H(2,:),H(3,:) )
  beta = ANGVEC( H(3,:),H(1,:) )
 gamma = ANGVEC( H(1,:),H(2,:) )
!
!
END SUBROUTINE MATCONV
!
!
!
!
!********************************************************
! CART2FRAC
! Converts the 3 first columns of an array A, assumed
! to contain cartesian coordinates, into fractional
! coordinates.
!********************************************************
SUBROUTINE CART2FRAC(A,H)
!
IMPLICIT NONE
INTEGER:: i
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:):: A
!
IF( SIZE(A(:,1)).NE.0 .AND. SIZE(A(1,:))>=3 ) THEN
    CALL INVMAT(H,G)
    DO i=1,SIZE(A(:,1))
      P1 = A(i,1)
      P2 = A(i,2)
      P3 = A(i,3)
      A(i,1) = P1*G(1,1) + P2*G(2,1) + P3*G(3,1)
      A(i,2) = P1*G(1,2) + P2*G(2,2) + P3*G(3,2)
      A(i,3) = P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
    END DO
ELSE
    CALL SYMOPS_ERR("Could not transform to fractional, inconsistent array size.")
ENDIF
!
END SUBROUTINE CART2FRAC
!
!
!********************************************************
! FRAC2CART
! Converts the 3 first columns of an array A, assumed
! to contain fractional coordinates, into cartesian
! coordinates.
!********************************************************
SUBROUTINE FRAC2CART(A,H)
!
IMPLICIT NONE
INTEGER:: i
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(:,:):: A
!
IF( SIZE(A(:,1)).NE.0 .AND. SIZE(A(1,:))>=3 ) THEN
    DO i=1,SIZE(A(:,1))
      P1 = A(i,1)
      P2 = A(i,2)
      P3 = A(i,3)
      A(i,1) = P1*H(1,1) + P2*H(2,1) + P3*H(3,1)
      A(i,2) = P1*H(1,2) + P2*H(2,2) + P3*H(3,2)
      A(i,3) = P1*H(1,3) + P2*H(2,3) + P3*H(3,3)
    ENDDO
ELSE
    call SYMOPS_ERR("Could not transform to cartesian, inconsistent array size.")
ENDIF
!
END SUBROUTINE FRAC2CART
!
!
!
!
!********************************************************
!  VECLENGTH
!  This function calculates the length of a vector.
!********************************************************
FUNCTION VECLENGTH(V) RESULT(Vlength)
!
IMPLICIT NONE
REAL(dp), DIMENSION(3),INTENT(IN):: V
REAL(dp):: Vlength
!
Vlength = DSQRT( V(1)**2 + V(2)**2 + V(3)**2 )
!
END FUNCTION VECLENGTH
!
!
!
!********************************************************
!  INN_PRODUCT
!  This routine calculates the inner product of two 3D
!  vectors. ( C = A (x) B )
!********************************************************
SUBROUTINE INN_PRODUCT(A, B, C)
!
IMPLICIT NONE
REAL(dp), DIMENSION(3), INTENT(IN) :: A, B
REAL(DP), DIMENSION(3), INTENT(OUT) :: C
!
C(1) = A(2)*B(3)-A(3)*B(2)
C(2) = A(3)*B(1)-A(1)*B(3)
C(3) = A(1)*B(2)-A(2)*B(1)
!
END SUBROUTINE INN_PRODUCT
!
!
!
!********************************************************
!  CHOP_ARR
!  This routine chops small values from an array.
!********************************************************
SUBROUTINE CHOP_ARR(A, VCHOP)
!
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(INOUT) :: A
REAL(DP), INTENT(IN) :: VCHOP
INTEGER :: N, I
REAL(DP) :: DCHOP
!
N = SIZE(A,1)
DCHOP = DABS(VCHOP)
IF (N>0) THEN
  DO I=1, N
    IF (ABS(A(I))<DCHOP) A(I) = 0.D+0
  END DO
END IF
!
END SUBROUTINE CHOP_ARR
!
!
!
!********************************************************
!  CHOP_MAT
!  This routine chops small values from a matrix (2-rank).
!********************************************************
SUBROUTINE CHOP_MAT(A, VCHOP)
!
IMPLICIT NONE
REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: A
REAL(DP), INTENT(IN) :: VCHOP
INTEGER :: N, M, I, J
REAL(DP) :: DCHOP
!
N = SIZE(A,1)
M = SIZE(A,2)
DCHOP = DABS(VCHOP)
IF (N>0.AND.M>0) THEN
  DO J=1, M
    DO I=1, N
      IF (ABS(A(I,J))<DCHOP) A(I,J) = 0.D+0
    END DO
  END DO
END IF
!
END SUBROUTINE CHOP_MAT
!
!
!
!********************************************************
!  ROUND_VAL
!  This routine rounds a value with given precision.
!  Warning: The value of the input "A" is changed.
!********************************************************
SUBROUTINE ROUND_VAL(A, VPREC)
!
IMPLICIT NONE
REAL(DP), INTENT(INOUT) :: A
REAL(DP), INTENT(IN) :: VPREC
REAL(DP) :: DPREC
!
DPREC = DABS(VPREC)
A = DNINT( A / DPREC ) * DPREC
!
END SUBROUTINE ROUND_VAL
!
!
!
!********************************************************
!  ROUND_ARR
!  This rounds 1D array values with given precision.
!  Warning: The values of the input "A" is changed.
!********************************************************
SUBROUTINE ROUND_ARR(A, VPREC)
!
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(INOUT) :: A
REAL(DP), INTENT(IN) :: VPREC
REAL(DP) :: DPREC, DTMP
INTEGER :: I, N
!
N = SIZE(A,1)
DPREC = DABS(VPREC)
IF (N>0) THEN
  DO I=1, N
    DTMP = DNINT( A(I) / DPREC )
    A(I) = DNINT( A(I) / DPREC ) * DPREC
  END DO
END IF
!
END SUBROUTINE ROUND_ARR
!
!
!
!********************************************************
!  ROUND_MAT
!  This rounds matrix values with given precision.
!  Warning: The values of the input "A" is changed.
!********************************************************
SUBROUTINE ROUND_MAT(A, VPREC)
!
IMPLICIT NONE
REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: A
REAL(DP), INTENT(IN) :: VPREC
REAL(DP) :: DPREC
INTEGER :: I, J, N, M
!
N = SIZE(A,1)
M = SIZE(A,2)
DPREC = DABS(VPREC)
IF (N>0 .AND. M>0) THEN
  DO J=1, M
    DO I=1, N
      A(I,J) = DNINT( A(I,J) / DPREC ) * DPREC
    END DO
  END DO
END IF
!
END SUBROUTINE ROUND_MAT
!
!
!
!
!********************************************************
!  ANGVEC
!  This function calculates the angle between two vectors
!  V1 and V2. Angle theta is returned in radians.
!********************************************************
FUNCTION ANGVEC(V1,V2) RESULT(theta)
!
IMPLICIT NONE
REAL(dp),DIMENSION(3),INTENT(IN):: V1, V2
REAL(dp):: theta
!
theta = DOT_PRODUCT(V1,V2) /                                    &
      & ( DSQRT(DOT_PRODUCT(V1,V1))*DSQRT(DOT_PRODUCT(V2,V2)) )
!
IF(theta>1.d0) THEN
  theta = 1.d0
ELSEIF(theta<-1.d0) THEN
  theta = -1.d0
ENDIF
!
theta = DACOS(theta)
!
RETURN
!
END FUNCTION ANGVEC
!
!
!
SUBROUTINE SYMOPS_TRF_UANI(uani0,uani1,H,TRF)
!
IMPLICIT NONE
!
REAL(dp),DIMENSION(6),INTENT(IN)  :: uani0 ! original Uani coefficients (U11,U22,U33,U12,U13,U23)
REAL(dp),DIMENSION(6),INTENT(OUT) :: uani1 ! transformed Uani coefficients (U11',U22',U33',U12',U13',U23')
REAL(dp),DIMENSION(3,3),INTENT(IN):: H     ! Base vectors of the supercell ((a11,a12,a13),(a21,... )
REAL(dp),DIMENSION(3,3),INTENT(IN):: TRF   ! Transformation matrix in the supercell basis
!                                          ! TRF transforms fractional real-space coordinates like
!                                          ! X' = TRF.X
!
CHARACTER(LEN=128) :: temp, msg            ! strings
REAL(dp),DIMENSION(3):: lh                 ! lengths of the basis vectors
REAL(dp),DIMENSION(3,3):: Q                ! modified transformation matrix
INTEGER:: i,j
INTEGER,DIMENSION(3,3):: li                ! index look-up table
!
uani1 = 0.d+0
li = RESHAPE( (/ 1,4,5, 4,2,6, 5,6,3 /), (/ 3, 3 /) ) ! look-up hash -> Uani matrix
! get the lengts of the super-cell's basis vectors
DO i=1,3
  lh(i) = VECLENGTH(H(i,1:3))
END DO
! tranform the transformation matrix into the basis of the Uani matrix
DO j=1,3
  DO i=1,3
    Q(i,j) = lh(i)*TRF(i,j)/lh(j)
  END DO
END DO
! do the second order parameter transformation
DO j=1,3
  DO i=1,3
    uani1(1) = uani1(1) + uani0(li(i,j)) * Q(1,i) * Q(j,1)
    uani1(2) = uani1(2) + uani0(li(i,j)) * Q(2,i) * Q(j,2)
    uani1(3) = uani1(3) + uani0(li(i,j)) * Q(3,i) * Q(j,3)
    uani1(4) = uani1(4) + uani0(li(i,j)) * Q(2,i) * Q(j,1)
    uani1(5) = uani1(5) + uani0(li(i,j)) * Q(3,i) * Q(j,1)
    uani1(6) = uani1(6) + uani0(li(i,j)) * Q(3,i) * Q(j,2)
  END DO
END DO
!
RETURN
!
END SUBROUTINE SYMOPS_TRF_UANI
!
!
!
SUBROUTINE SYMOPS_APPLY(H,P,dmindist,nchk)
!
IMPLICIT NONE
!
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE::P   !Structure data
REAL(DP),INTENT(IN) :: dmindist ! minimum atom distance
INTEGER,INTENT(OUT) :: nchk ! success code: 0=failure, 1=success
!
CHARACTER(LEN=128) :: temp, msg ! strings
INTEGER::i,j,isym,occ,idup,j3
INTEGER::n1,n2
INTEGER::MP,NP,Nsym
REAL(dp)::a,b,c,alpha,beta,gamma
REAL(dp)::trf(symops_nltrf),rtmp,rfdif(3),rdif(3),rdist,p1(3),p2(3),ds
REAL(dp)::socc1,socc2,soccs
INTEGER,DIMENSION(:),ALLOCATABLE::SUSE ! intermediate array for new data
REAL(dp),DIMENSION(:,:),ALLOCATABLE::S1,S2,SAUX1,SAUX2 ! intermediate arrays for new data 
!REAL(dp),DIMENSION(3,3):: G ! inverse of H
!
msg = 'entering SYMOPS_APPLY'
CALL SYMOPS_MSG(msg,1)
!
!Initial checks.
nchk=0
IF (.NOT.ALLOCATED(P)) RETURN ! Just exit, since there are no atoms.
IF (.NOT.ALLOCATED(symops_trf)) RETURN ! Just exit, since no symmetries are defined.
!
!Initialization
!CALL MATCONV(H,a,b,c,alpha,beta,gamma) ! not needed here, because H is input
!CALL INVMAT(H,G) ! not needed here
Nsym=SIZE(symops_trf,2)
IF (Nsym<=0) RETURN ! Just exit, since no symmetries are defined.
MP=SIZE(P,1)
NP=SIZE(P,2) 
IF (NP<=0) RETURN ! Just exit, since there are no atoms
WRITE(msg,*) NP
msg = '- input number of atoms: '//TRIM(ADJUSTL(msg))
CALL SYMOPS_MSG(msg,2)
WRITE(msg,*) Nsym
msg = '- number of symmetry operations: '//TRIM(ADJUSTL(msg))
CALL SYMOPS_MSG(msg,2)
idup=0
n1=NP
n2=2*n1
!Prepare initial arrays
ALLOCATE(S1(MP,n1),S2(MP,n2),SUSE(n2))
S1=P ! copy P to S1
!CALL CART2FRAC(S1,H) ! transform to fractional coordinates ! we work on fractional coordinates only
S1(1:3,1:n1)=MODULO(S1(1:3,1:n1),1.0) ! wrap into [0,1[
S2=0.d0 ! reset S2
S2(1:MP,1:n1)=S1(1:MP,1:n1) ! fill first half of S2 with S1
!
SUSE=0 ! reset SUSE
SUSE(1:n1)=1 ! mark the copied content of P to be used
!
DO isym=1,Nsym ! apply each symmetry op. ...
  !Get current transformation parameters
  trf = symops_trf(:,isym)
  !Apply this transformation from S1 to S2
  DO i=1,n1 ! ... to each of the current positions
    ! get current atomic position in fractional coordinates
    p1(1:3) = MODULO( S1(1:3,i), 1.0 ) ! wrapped into the unit cell
    DO j=1,3
      j3=3*(j-1)
      ! get new fractional coordinate from symmetry operation
      rtmp = trf(j)+trf(4+j3)*p1(1)+trf(5+j3)*p1(2)+trf(6+j3)*p1(3)
      ! wrap it periodically into the cell
      S2(j,i+n1) = MODULO( rtmp, 1.d0 )
    ENDDO
    IF (MP>3) THEN
      !
      S2(4:MP,i+n1) = S1(4:MP,i) ! copy the rest of the data row (at 4 we have the atomic numbers)
      !
      ! The second rank tensor of anisotropic thermal displacements
      ! should be invariant upon the crystallographic symmetry operations
      ! in a proper structure definition.
      ! I think this is not entirely true, because this applies to
      ! the symmetry of the specific atomic site only and not to the
      ! symmetry of the whole structure. Think of mirror operations for
      ! example. Nevertheless, we currently need isotropic displacement
      ! parameters for the diffraction calculations, which depend only
      ! on the unit cell basis vectors. We thus leave the Uani table as
      ! is.
      ! 
    END IF
    !
    ! Now that we have a new atomic site, we should check if this atom is
    ! a duplicate, which is absolutely identical to one of the old atoms.
    ! Identical means S2(1:4,i)==S2(1:4,j), including the atomic number.
    ! We will not allow duplicates here and leave them unmarked (0) in the
    !  array "SUSE". We will mark all unique new positions (1) in "SUSE".
    !
    ! Note on partial site occupancies:
    ! The above check allows partial site occupancies by different atomic
    ! species, since their atomic number is different. It doesn't allow a
    ! site sharing by two identical species, regardless of their occupancy
    ! factor. This would be a very strange structure definition anyway.
    !
    idup=0 ! reset the duplicate indicator
    p1 = S2(1:3,i+n1) ! get the new symmetry position
    DO j=1,n1+i-1 ! look through the S2 data up to the new position
      IF (SUSE(j)==0) CYCLE ! skip any unmarked previous position.
      p2 = S2(1:3,j) ! get fractional coordinate (should already be in [0,1[ )
      rfdif = MODULO(p2-p1,1.0) ! fractional distance vector (modulo unit cell)
      ! Note: We are looking for the shortest fractional distance
      !       under periodic boundary conditions. This may require to look
      !       for shorter vectors across the cell boundary:
      IF (rfdif(1)>0.5d0) rfdif(1) = rfdif(1)-1.0d0
      IF (rfdif(2)>0.5d0) rfdif(2) = rfdif(2)-1.0d0
      IF (rfdif(3)>0.5d0) rfdif(3) = rfdif(3)-1.0d0
      ! calculate cartesian distance vector
      rdif(1) = SUM( H(1:3,1)*rfdif )
      rdif(2) = SUM( H(1:3,2)*rfdif )
      rdif(3) = SUM( H(1:3,3)*rfdif )
      rdist = VECLENGTH(rdif) ! get distance in A
      ! check if this is a "duplicate"
      ! ( distance to other atom < min. distance and atomic numbers equal )
      IF (rdist<=dmindist .AND. DABS(S2(4,j)-S2(4,i))<0.5d0) THEN
        idup = idup + 1 ! raise the duplicate count
        EXIT ! Exit the loop over other atoms, we found one duplicate, that's enough.
      ENDIF
    ENDDO
    ! handle non-duplicates
    IF (idup==0) SUSE(i+n1)=1 ! mark new position as unique
    !
  ENDDO ! (positions)
  !
  ! There is a new set of positions in S2.
  ! We transfer it now back to S1 and reinitiate everything.
  ! Update number of atoms n1
  n1 = SUM(SUSE)
  IF (ALLOCATED(S1)) DEALLOCATE(S1)
  ALLOCATE(S1(1:MP,1:n1))
  ! Transfer the usable data from S2 to the new S1
  j = 0
  DO i=1, n2
    IF (SUSE(i)==0) CYCLE ! skip unused atoms
    j = j + 1 ! raise atom list index
    S1(1:MP,j) = S2(1:MP,i) ! positions and species
  ENDDO
  ! Update the secondary array set
  n2 = 2*n1 ! possible number of new atoms
  IF (ALLOCATED(S2)) DEALLOCATE(S2)
  IF (ALLOCATED(SUSE)) DEALLOCATE(SUSE)
  ALLOCATE(S2(1:MP,1:n2),SUSE(1:n2)) ! new S2 and SUSE
  S2 = 0.d0
  SUSE = 0
  S2(1:MP,1:n1) = S1(1:MP,1:n1) ! copy S1
  SUSE(1:n1) = 1 ! mark S1 as usable
  !
ENDDO ! (symmetry operations)
!
! Finally update P
DEALLOCATE(P)
NP = n1 ! set new number of atoms
ALLOCATE(P(1:MP,1:NP))
P = S1 ! set new atomic site data
!CALL FRAC2CART(P,H) ! transform P to cartesian coordinates ! we work on fractional coordinates only
! Clean up the heap.
IF (ALLOCATED(S1)) DEALLOCATE(S2)
IF (ALLOCATED(SUSE)) DEALLOCATE(SUSE)
WRITE(msg,*) NP
msg = '- output number of atoms: '//TRIM(ADJUSTL(msg))
CALL SYMOPS_MSG(msg,2)
nchk=1
!
END SUBROUTINE SYMOPS_APPLY
!
!
!
!
SUBROUTINE SYMOPS_CHECK_STR(instr,nchk)
!
INTEGER,PARAMETER:: numtests = 2 ! NUMBER OF CHECKS: INCREASE IF YOU IMPLEMENT MORE!
CHARACTER(LEN=*), INTENT(IN) :: instr ! input string
INTEGER, INTENT(OUT) :: nchk ! output check result:
                             !   0: no symmetry operation
                             !   1: symmetry operation
INTEGER:: i, l, c1, c2, ichk ! iterators, helpers, and subcheck results
CHARACTER(LEN=128) :: temp, msg, symopchars ! strings
!
! Initialize
nchk=0
ichk=0
i=0
c1=0
c2=0
temp=ADJUSTL(instr)
l=LEN_TRIM(temp)
symopchars='xyz1234567890+-*/.'
!
msg = 'entering SYMOPS_CHECK_STR'
CALL SYMOPS_MSG(msg,2)
msg = "- checking '"//TRIM(temp)//"' for CIF symmetry operation style"
CALL SYMOPS_MSG(msg,3)
!
!Check #1: A well defined operation contains the right sides
!          of 3 equations, separated by extacly 2 "," characters
c1=INDEX(TRIM(temp(1:l)),",",BACK=.FALSE.)
c2=INDEX(TRIM(temp(1:l)),",",BACK=.TRUE.)
IF (c1>0 .AND. c2>0 .AND. c2>c1) ichk = ichk + 1
IF (ichk==0) THEN
  msg = "- failed the 2-comma test."
  CALL SYMOPS_MSG(msg,3)
  RETURN
END IF
!Remove the commas.
temp(c1:c1) = " "
temp(c2:c2) = " "
!
!Check #2: After removing the 2 commas, the string should contain
!          only the characters in symopchars.
!Replace all the occurrences of symopchars in temp by " ":
c1=0
DO i=1, l
  c1=INDEX(TRIM(symopchars),temp(i:i))
  IF (c1>0) temp(i:i) = " "
END DO
c2=LEN_TRIM(temp(1:l))
IF (c2>0) THEN
  msg = "- failed the restricted characters test."
  CALL SYMOPS_MSG(msg,3)
  RETURN
ENDIF
ichk = ichk+1 ! 2nd check passed.
!
IF (ichk==numtests) THEN
  nchk=1 ! success, this is a symmetry operation
  temp=ADJUSTL(instr)
  msg = "- '"//TRIM(temp)// &
      & "' passed all tests for a CIF symmetry op."
  CALL SYMOPS_MSG(msg,3)
ENDIF
!
END SUBROUTINE SYMOPS_CHECK_STR
!
!
!
!
SUBROUTINE SYMOPS_SET_STR(instr,irow,nchk)
! This routine transforms instr to linear transformation parameters
! of the form M.(x,y,z)+S, where S is a 3D shift vector and and M
! is a 3x3 tranformation matrix.
! The values of S and M a stored in sequence (S,M) in the array
! symops_trf(:,irow) of the module symops (link 'symops.f90'!).
! symops_trf must be preallocated and initialized.
! 'instr' should be of the form 'x-y, y+1/2, -z/2', which would result
! in S = (/ 0, 0.5, 0 /) and
! and M = (/ (/ 1, -1, 0 /), (/ 0, 1, 0 /), (/ 0, 0, -0.5 /) /).
! The input parameter 'irow' determines the row of the data storage.
! The output parameter nchk signalizes the success or failur of the
! translation in this routine.
!
CHARACTER(LEN=*),INTENT(IN):: instr ! input string
INTEGER,INTENT(IN):: irow ! input row index in the transformation table
INTEGER, INTENT(OUT):: nchk ! output success code: 0: failure, 1: success
!
CHARACTER(LEN=128) :: temp, msg, substr(3) ! strings
INTEGER:: i, j, k, l ! iterators
INTEGER:: c1, c2 ! positions of substring separation
REAL(dp):: trf(symops_nltrf) ! local transformation
!
! Initialize
trf = 0.d0
trf(4)  = 1.d0
trf(8)  = 1.d0
trf(12) = 1.d0
substr=""
c1=0
c2=0
i=0
j=0
k=0
temp=ADJUSTL(instr)
l=LEN_TRIM(temp)
nchk=0
IF (.NOT.ALLOCATED(symops_trf)) GOTO 850
IF (SIZE(symops_trf,1).NE.symops_nltrf) GOTO 850
IF (SIZE(symops_trf,2)<irow) GOTO 850
! The minimum valid string is 'x,y,z'.
IF (l<5) GOTO 801 ! warn & report invalid CIF symmetry operation string.
!
msg = 'entering SYMOPS_SET_STR'
CALL SYMOPS_MSG(msg,2)
!
! Determine the two comma positions
c1=INDEX(temp(1:l),",",BACK=.FALSE.)
c2=INDEX(temp(1:l),",",BACK=.TRUE.)
! There should be two commas in the string. The first comma should be at
! a position >1, such that the first substring fits in before it. The
! second comma should come at least 2 characters beyond the first comma,
! such that the second substring fits in between. The second comma
! shouldn't be at the end of the string, such that there is still space
! for the third substring.
IF (c1<2 .OR. c2-c1<2 .or. c2>=l) GOTO 801
! Extract the 3 substrings
substr(1) = TRIM(ADJUSTL(temp(1   :c1-1)))
substr(2) = TRIM(ADJUSTL(temp(c1+1:c2-1)))
substr(3) = TRIM(ADJUSTL(temp(c2+1:l   )))
! Parsing the 3 substrings
temp = TRIM(substr(1))
call SYMOPS_PARSE_STR_LINTRF(temp,trf(1),trf(4:6),k)
IF (k==0) GOTO 801
temp = TRIM(substr(2))
call SYMOPS_PARSE_STR_LINTRF(temp,trf(2),trf(7:9),k)
IF (k==0) GOTO 801
temp = TRIM(substr(3))
call SYMOPS_PARSE_STR_LINTRF(temp,trf(3),trf(10:12),k)
IF (k==0) GOTO 801
!
! Store the transformation in the symops module
symops_trf(:,irow) = trf(:)
nchk=1 ! success
! debug out
msg = "- found transformation for '"//TRIM(ADJUSTL(instr))//"':"
CALL SYMOPS_MSG(msg,3)
WRITE(msg,901) trf(1), trf(4), trf(5), trf(6)
msg = "  x' = "//TRIM(ADJUSTL(msg))
!CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,901) trf(2), trf(7), trf(8), trf(9)
msg = "  y' = "//TRIM(ADJUSTL(msg))
!CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,901) trf(3), trf(10), trf(11), trf(12)
msg = "  z' = "//TRIM(ADJUSTL(msg))
!CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
901 FORMAT(F8.3," + ",F8.3,"*x + ",F8.3,"*y + ",F8.3,"*z")

GOTO 1000 ! skip warnings / errors.
!
801 CONTINUE
CALL SYMOPS_WRN(TRIM(temp))

GOTO 1000
850 CONTINUE
msg = 'Unable to store transformation in symops module. skipping'
CALL SYMOPS_ERR(trim(msg),1707)
!
1000 CONTINUE
!
END SUBROUTINE SYMOPS_SET_STR
!
!
!
!
SUBROUTINE SYMOPS_PARSE_STR_LINTRF(instr,shift,slope,nchk)
!
CHARACTER(LEN=*),INTENT(IN):: instr
REAL(dp),INTENT(OUT)::shift,slope(3)
INTEGER,INTENT(OUT)::nchk
!
INTEGER:: i,j,l
CHARACTER(LEN=128):: temp,msg
INTEGER:: ichan, k1, k2
REAL(dp):: rtmp, sh1, sh2, sl1(3), sl2(3)
!
nchk=0
shift=0.d0
slope=0.d0
temp=ADJUSTL(instr)
l=LEN_TRIM(temp)
!
msg = 'entering SYMOPS_PARSE_STR_LINTRF'
CALL SYMOPS_MSG(msg,2)
msg = "- parsing '"//temp(1:l)//"'..."
CALL SYMOPS_MSG(msg,3)
!
IF (l<=2) GOTO 400 ! TRIVIAL CASE: numbers or characters, no operation
!
! Check for basic operations: the operation character should be beyond
! position 1 of the string "temp". Store the position in "j".
j=INDEX(temp(1:l),'+')
IF (j>1) GOTO 500
j=INDEX(temp(1:l),'-')
IF (j>1) GOTO 520
j=INDEX(temp(1:l),'*')
IF (j>1) GOTO 540
j=INDEX(temp(1:l),'/')
IF (j>1) GOTO 560
! No operation character found. We assume now, that temp is a trivial
! number or character, which would mean something like "0.5" or "-10"
! since temp is longer than 2.
!
400 CONTINUE ! TRIVIAL CASE: single numbers or single characters, no operation
!                          ! Solve it here !
! Determine the output channel: shift=0, x=1, y=2, z=3
ichan=0
j=0
DO i=1,3
  j=INDEX(temp(1:l),symops_chanstr(i:i))
  IF(j>0) THEN
    ichan=i ! got the character channel
    temp(j:j)='1' ! replace the character by a number
    EXIT ! done, exit this loop
  ENDIF
END DO
! "temp" should now contain just a trivial number
! Translate to a numerical value
READ(temp,*,ERR=801,END=801) rtmp
! Put numerical value into channel
IF(ichan==0) THEN
  shift = rtmp
ELSE
  slope(ichan) = rtmp
ENDIF
! Success of trivial case
GOTO 800
!
500 CONTINUE ! ADDITION CASE       : operation '+' with two substrings
!            ! Bifurcate here !
call SYMOPS_PARSE_STR_LINTRF(temp(1:j-1),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(temp(j+1:l),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine
shift = sh1+sh2
slope = sl1+sl2
GOTO 800
!
520 CONTINUE ! SUBTRACTION CASE    : operation '-' with two substrings
!            ! Bifurcate here !
call SYMOPS_PARSE_STR_LINTRF(temp(1:j-1),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(temp(j+1:l),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine
shift = sh1-sh2
slope = sl1-sl2
GOTO 800
!
540 CONTINUE ! MULTIPLICATION CASE : operation '*' with two substrings
!            ! Bifurcate here !
call SYMOPS_PARSE_STR_LINTRF(temp(1:j-1),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(temp(j+1:l),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine: Possible scenarios
! 1) S1=A,M1=0, S2=0,M2=X -> S=0, M=A*X -> S=S1*S2, M=S1*M2+S2*M1
! 2) S1=0,M1=X, S2=A,M2=0 -> S=0, M=X*A -> S=S1*S2, M=S1*M2+S2*M1
! 3) S1=A,M1=0, S2=B,M2=0 -> S=A*B, M=0 -> S=S1*S2, M=S1*M2+S2*M1
! - no scenario where M1/=0 and M2/=0, this would be a non-linear form
! - no scenario like A*(B+X), since brackets are not supported.
shift = sh1*sh2
slope = sh1*sl2+sh2*sl1
GOTO 800
!
560 CONTINUE ! DIVISION CASE       : operation '/' with two substrings
!            ! Bifurcate here !
call SYMOPS_PARSE_STR_LINTRF(temp(1:j-1),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(temp(j+1:l),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine: Possible scenarios
! 1) S1=0,M1=X, S2=A,M2=0 -> S=0, M=X/A -> S=S1/S2, M=M1/S2
! 2) S1=B,M1=0, S2=A,M2=0 -> S=B/A, M=0 -> S=S1/S2, M=M1/S2
! - requires S2/=0 and M2==0 ! Only division by pure numbers
! - no scenario with A/M2 or M1/M2, this would be a non-linear form
IF (sh2==0.d0) GOTO 802 ! error, division by zero
IF ( sl2(1)/=0.d0 .OR. sl2(2)/=0.d0 .OR. sl2(3)/=0.d0 ) GOTO 802 ! error, division by variable
shift = sh1/sh2
slope = sl1/sh2
GOTO 800
!
!
800 CONTINUE ! success
nchk = 1
! debug out
!WRITE(msg,'(F8.3," + ",F8.3,"*x + ",F8.3,"*y + ",F8.3,"*z")') &
!     &    shift, slope(1), slope(2), slope(3)
!msg = "- transformation: "//TRIM(ADJUSTL(msg))
!CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
GOTO 1000

801 CONTINUE ! error, failed to convert to number
CALL SYMOPS_ERR(temp(1:l),808)
GOTO 1000

802 CONTINUE ! error, failed to convert operation
CALL SYMOPS_ERR(temp(1:l),1813)
GOTO 1000
!
1000 CONTINUE
!
END SUBROUTINE SYMOPS_PARSE_STR_LINTRF
!
!
!
!
SUBROUTINE SYMOPS_SET_SGNAME(sgname,nchk)
!
CHARACTER(LEN=*),INTENT(IN)::sgname ! Hermann-Mauguin symbol identifying the space group
INTEGER,INTENT(OUT)::nchk ! output success code: 0: failure, 1: success
!
INTEGER::nsgnum,i
CHARACTER(LEN=128)::temp,msg
!
!Initialization
nchk=0
nsgnum=0
i=0
msg = 'entering SYMOPS_SET_SGNAME'
CALL SYMOPS_MSG(msg,2)
!
!Get the space group number. This call will initialize the spacegroup
!module if it isn't initialized yet.
temp=ADJUSTL(sgname)
CALL SG_NAMGETNUM(TRIM(temp),nsgnum)
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 801 ! invalid space group name
WRITE(msg,*) nsgnum
msg = "- identified space group '"//TRIM(temp)//"' as number "// &
    & TRIM(ADJUSTL(msg))//"."
CALL SYMOPS_MSG(msg,3)
!
!Apply set the symmetry operations based on the SG number
CALL SYMOPS_SET_SGNUM(nsgnum,i)
IF (i.NE.1) THEN
  nchk=0
  GOTO 1000 ! pipe error through, error messages were posted by SYMOPS_SET_SGNUM
ENDIF
!
800 CONTINUE ! success
nchk=1
GOTO 1000
!
!Error handling
801 CONTINUE ! error, invalid space group name
CALL SYMOPS_ERR("Invalid space group name: "//TRIM(temp),809)
nchk = 0
GOTO 1000
!
!Routine exit.
1000 CONTINUE
!
END SUBROUTINE SYMOPS_SET_SGNAME
!
!
!
!
SUBROUTINE SYMOPS_SET_SGNUM(nsgnumber,nchk)
!
INTEGER,INTENT(IN)::nsgnumber ! space group number (1 - 230)
INTEGER,INTENT(OUT)::nchk ! output success code: 0: failure, 1: success
!
!
INTEGER::nsgnum,nsymnum,ichk,i,nfail
CHARACTER(LEN=128)::temp,msg
!
!Initialization
nchk=0
ichk=0
nsgnum=nsgnumber
nsymnum=0
nfail=0
msg = 'entering SYMOPS_SET_SGNUM'
CALL SYMOPS_MSG(msg,2)
!
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 801 ! invalid space group number
!Get number of symmetry operations for the space group
CALL SG_NUMGETSYMNUM(nsgnum,nsymnum,ichk)
IF (ichk.NE.1) THEN
  IF (ichk==-1) GOTO 802
  IF (ichk==-2) GOTO 801
ENDIF
IF (nsymnum.LE.0) GOTO 802
!Prepare the symmetry data array
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
ALLOCATE(symops_trf(symops_nltrf,nsymnum))
!Initialize symmetry operations
CALL SYMOPS_INIT()
!
!Set the operations from the space group data
DO i=1, nsymnum
  !Get the symmetry operation string
  temp=''
  CALL SG_NUMGETSYMOP(nsgnum,i,temp(1:sg_soplen),ichk)
  IF (ichk.LE.0) GOTO 802
  !Set the operation parameters from the string
  CALL SYMOPS_SET_STR(TRIM(temp),i,ichk)
  IF (ichk.LE.0) THEN ! report interpretation error
                      ! do not stop here, go on
    CALL SYMOPS_ERR("Failed to set symmetry OP: "//TRIM(temp),812)
    nfail = nfail+1
  ENDIF
ENDDO
!
!
800 CONTINUE 
IF (nfail==0) nchk=1 ! success
GOTO 1000
!
!Error handling
801 CONTINUE ! error, invalid space group number
WRITE(unit=temp,fmt=*) nsgnum
CALL SYMOPS_ERR("Invalid space group number: "//TRIM(temp),810)
nchk = 0
GOTO 1000
!
802 CONTINUE ! error, failed to access space group data
WRITE(unit=temp,fmt=*) nsgnum
CALL SYMOPS_ERR("Failed to access space group: "//TRIM(temp),811)
nchk = 0
GOTO 1000
!
803 CONTINUE ! error, failed to access space group data
CALL SYMOPS_ERR("Failed to access space group: "//TRIM(temp),812)
nchk = 0
GOTO 1000
!
!Routine exit.
1000 CONTINUE
!
END SUBROUTINE SYMOPS_SET_SGNUM
!
!
!
!
END MODULE symops
