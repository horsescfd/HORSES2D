!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      Qsort.f90
!      Created: 2016-08-25 14:30 (GMT+0)
!      By: Juli Rew (juliana@ucar.edu)
! 		Recursive Fortran 95 quicksort rOUTINe
! 			sorts INTEGER numbers INto ascENDing numerical order
!  		Based on algorithm from Cormen et al., Introduction to Algorithms,
! 			1997 prINting
!
!      ptau2.5.preprocessINg Code
!	
!      NNATAC Project
!
!////////////////////////////////////////////////////////////////////////////////////////
!

MODULE Sorting

IMPLICIT NONE
PUBLIC  :: Qsort
PRIVATE :: Partition

CONTAINS

RECURSIVE SUBROUTINE Qsort(A)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER                               :: iq

  IF(size(A) > 1) THEN
     call Partition(A, iq)
     call Qsort(A(:iq-1))
     call Qsort(A(iq:))
  ENDIF
END SUBROUTINE Qsort

SUBROUTINE Partition(A, marker)
  INTEGER, INTENT(IN OUT), DIMENSION(:) :: A
  INTEGER, INTENT(OUT)                  :: marker
  INTEGER                               :: i, j
  INTEGER                               :: temp
  INTEGER                               :: x      ! pivot poINt
  x = A(1)
  i = 0
  j = size(A) + 1

  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseIF (i == j) THEN
        marker = i+1
        return
     else
        marker = i
        return
     ENDIF
  END DO

END SUBROUTINE Partition

END MODULE Sorting
