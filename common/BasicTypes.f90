!******************************************************************************!
!******************************************************************************!
!                                                                              !
!                            file   "BasicTypes".f95                           !
!                                                                              !
!    Copyright:  J.Barthel, Forschungszentrum Juelich                          !
!    Version  :  1.0.0, October 06, 2003                                       !
!                                                                              !
!                                                                              !
!******************************************************************************!


!******************************************************************************!
!                                                                              !
!   Purpose: implementing MODULE BasicTypes                                    !
!            accessable from F90/95                                            !
!                                                                              !
!   External  sources to link: none                                            !
!                                                                              !
!******************************************************************************!
!                                                                              !
!   CONTAINS:                                                                  !
!      1) MODULE Basic Types                                                   !
!         t_point            | integer*4 2d-point type                         !
!         t_rect             | integer*4 2d-rectangle type                     !
!                                                                              !
!******************************************************************************!

!******************************************************************************!
!*********************** MODULE BasicTypes ************************************!
!******************************************************************************!
!******************************************************************************!
!                                                                              !
!   Purpose: Defining often used data types and procedures                     !
!                                                                              !
!                                                                              !
!******************************************************************************!
MODULE BasicTypes

    implicit none
    
    save
    
    public :: InfltRect, SgnRectX, SgnRectY, IsRectPositive

	!****************************************************************************!
	!********************* TYPE definitions *************************************!
	!****************************************************************************!
	! interger*4 2d-point
	TYPE t_point
		INTEGER*4 :: x
		INTEGER*4 :: y
	END TYPE t_point	
	! integer*4 2d-rectangle defined by low-bottom and right-top points
	TYPE t_rect
		TYPE(t_point) :: lb
		TYPE(t_point) :: rt
	END TYPE t_rect


	!****************************************************************************!
	!********************* SUBROUTINE IMPLEMENTATION ****************************!
	!****************************************************************************!
	CONTAINS
		
	! rect clipping in rcBase coordinates
	TYPE(t_rect) FUNCTION InfltRect(rcBase,rcStamp)
		TYPE(t_rect), INTENT(IN) :: rcBase, rcStamp
		
		TYPE(t_rect) :: rcResult
		
		rcResult%lb%x=MAX(rcBase%lb%x,rcStamp%lb%x)
		rcResult%lb%y=MAX(rcBase%lb%y,rcStamp%lb%y)
		rcResult%rt%x=MIN(rcBase%rt%x,rcStamp%rt%x)
		rcResult%rt%y=MIN(rcBase%rt%y,rcStamp%rt%y)

		InfltRect = rcResult

	END FUNCTION InfltRect
	

	! returns sign of a rectangle in x-direction
	INTEGER FUNCTION SgnRectX(rcBase)
	
		TYPE(t_rect), INTENT(IN) :: rcBase
		
		INTEGER*4 :: nDiff
		
		nDiff = rcBase%rt%x-rcBase%lb%x
		
		IF (nDiff/=0) THEN
			nDiff = SIGN(1,nDiff)
		END IF
		
		SgnRectX = nDiff
		
	END FUNCTION SgnRectX


	! returns sign of a rectangle in y-direction
	INTEGER FUNCTION SgnRectY(rcBase)
	
		TYPE(t_rect), INTENT(IN) :: rcBase
		
		INTEGER*4 :: nDiff
		
		nDiff = rcBase%rt%y-rcBase%lb%y
		
		IF (nDiff/=0) THEN
			nDiff = SIGN(1,nDiff)
		END IF
		
		SgnRectY = nDiff
		
	END FUNCTION SgnRectY
	
	! returns "TRUE" if both directions of rcBase are positivly defined	
	LOGICAL FUNCTION IsRectPositive(rcBase)

		TYPE(t_rect), INTENT(IN) :: rcBase
		
		LOGICAL :: lResult
		
			lResult = .TRUE.
			
			IF ( (SgnRectX(rcBase)==-1).OR.(SgnRectY(rcBase)==-1) ) THEN
				lResult = .FALSE.
			END IF
		
		IsRectPositive = lResult

	END FUNCTION IsRectPositive

END MODULE BasicTypes
!******************************************************************************!
