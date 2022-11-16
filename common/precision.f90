!**********************************************************************!
!**********************************************************************!
!                                                                      !
!    File     :  precision.f90                                         !
!                                                                      !
!    Copyright:  (C) J. Barthel (ju.barthel@fz-juelich.de) 2022        !
!                                                                      !
!**********************************************************************!
!                                                                      !
!    MODULE precision                                                  !
!    -----------------                                                 !
!                                                                      !
!    Purpose  : parameters concerning numerical precision              !
!    Version  : 1.0.0, Nov 9, 2022                                     !
!    To Link  : none                                                   !
!                                                                      !
!**********************************************************************!
!                                                                       
!  Author:  Juri Barthel                                
!           Ernst Ruska-Centre (ER-C2)
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
  

!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
  
  
MODULE precision

  implicit none
  
  integer*4, parameter, public :: fpp_ac = 8 ! accumulator precision

END MODULE precision


!**********************************************************************!
!**********************************************************************!
!**********************************************************************!