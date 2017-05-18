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
MODULE Headers
   use SMConstants
   implicit none

#include "Defines.h"

        INTEGER,        SAVE     :: iter
        INTEGER,        SAVE     :: loop_size
        INTEGER,        SAVE     :: no_of_points
        INTEGER,        SAVE     :: prev

        PRIVATE
        PUBLIC  Main_header, Ruler_Reset, Ruler_Update, Section_Header
        PUBLIC  SubSection_Header, Ruler_Header_Reset , Ruler_Header_Update
!
!       ========
        CONTAINS
!       ========
!       
        SUBROUTINE Main_header(title)
                IMPLICIT NONE
                CHARACTER(LEN = *)      :: title
                INTEGER, PARAMETER      :: siz = 103
                CHARACTER(LEN = siz)    :: ast1,ast2
                CHARACTER(LEN = siz)    :: compilationDate
                CHARACTER(LEN = siz)    :: astCompilation
                CHARACTER(LEN = siz)    :: astBox
                INTEGER                 :: i

                ast1 = ''
                DO i = 1 , siz
                        ast1 = TRIM(ast1) // '#'
                END DO
                ast2 = ''
                ast2(1:1) = '#'
                ast2(siz:siz) = '#'

                write(compilationDate,'(A,A,A,A)') "Compiled at ",__DATE__ , ", ",  __TIME__
                do i = 2 , siz-3
                  astCompilation(i:i) = " "
                end do
                astCompilation(1:1) = '#'
                astCompilation(siz-4-len_trim(compilationDate) : siz-2) = "  " // TRIM(compilationDate)
                astCompilation(siz-1:siz) = " #"
                do i = 1 , siz
                  astBox(i:i) = ","
                end do

                do i = siz-4-len_trim(compilationDate) , siz-1
                  astBox(i:i) = "-"
                end do

                astBox(1:1) = "#"
                astBox(siz:siz) = "#"
                WRITE(*,'(A)') ast1
                WRITE(*,'(A)') ast2
                WRITE(*,'(A)') ast2
#include "./header.incf"
                WRITE(*,'(A)') ast2
                WRITE(*,'(A)') astCompilation
                WRITE(*,'(A)') ast1

        END SUBROUTINE Main_header
     
        SUBROUTINE Section_header(title)
           IMPLICIT NONE
           CHARACTER(LEN=*)     :: title
           INTEGER, PARAMETER   :: siz = 86
           CHARACTER(LEN=siz)      :: dotted_title
           INTEGER        :: i
     
           dotted_title = ''
           
           dotted_title(1:LEN_TRIM(title)) = TRIM(title)
     
           DO i = 1 , siz - LEN_TRIM(title)
              dotted_title = TRIM(dotted_title) // '.'
           END DO   
           WRITE(*,'(5X,A,A)') "\\\\ ",TRIM(dotted_title)
           FLUSH(6)
           
        END SUBROUTINE Section_header
      
        SUBROUTINE SubSection_header(title)
           IMPLICIT NONE
           CHARACTER(LEN=*)     :: title
           INTEGER        :: i
           character(len=len_trim(title))  :: auxstring
     
           do i = 1 , len_trim(title)
            auxstring(i:i) = "-"
           end do
           WRITE(STD_OUT,'(15X,A)') TRIM(title)
           write(STD_OUT,'(15X,A)') trim(auxstring)
           FLUSH(6)
           
        END SUBROUTINE SubSection_header
   
        SUBROUTINE Ruler_Header_Reset(title , loop_in)
                IMPLICIT NONE
                character(len=*)          :: title
                integer                   :: loop_in
                integer, parameter        :: siz = 86
      
                iter          = 0
                loop_size     = loop_in
                no_of_points  = siz - len_trim(title)
                prev          = 0
                
                call flush(STD_OUT)               
                write(* , '(/,5X,A,A)',advance = "no") "\\\\ ",trim(title)
                call flush(STD_OUT)               

        end subroutine Ruler_Header_Reset


        SUBROUTINE Ruler_Header_Update()
                IMPLICIT NONE
        
                iter = iter + 1     
                IF(FLOOR(1.0d0*iter*no_of_points/loop_size) .GT. prev) THEN
                        WRITE(*,'(A)',ADVANCE='NO') '.'
                        prev = FLOOR(1.0d0*iter*no_of_points/loop_size)
                        call flush(6)
                END IF
                
                
        END SUBROUTINE Ruler_Header_Update
                
 
        SUBROUTINE Ruler_Reset(title,n_in,loop_in)
                IMPLICIT NONE
                CHARACTER(LEN=*)  :: title
                INTEGER                 :: loop_in
                INTEGER                 :: n_in

                WRITE(* , '(20X,A)' , ADVANCE = 'NO' ) TRIM(title)
                FLUSH(6)
                iter = 0
                loop_size = loop_in
                no_of_points    = n_in - LEN_TRIM(title)
                prev    = 0
        END SUBROUTINE Ruler_Reset

        SUBROUTINE Ruler_Update()
                IMPLICIT NONE
        
                iter = iter + 1     
                IF(FLOOR(1.0d0*iter*no_of_points/loop_size) .GT. prev) THEN
                        WRITE(*,'(A)',ADVANCE='NO') '.'
                        prev = FLOOR(1.0d0*iter*no_of_points/loop_size)
                        FLUSH(6)
                END IF
                
                
        END SUBROUTINE Ruler_Update
                
                
END MODULE Headers


