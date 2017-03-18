module Utilities
   use SMConstants
   implicit none
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////
!
      pure LOGICAL FUNCTION AlmostEqual( a, b ) 
!
!        *************************************+
!           Function by David A. Kopriva
!        *************************************+
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), intent(in) :: a, b
         IF ( a == 0.0_RP .OR. b == 0.0_RP )     THEN
            IF ( ABS(a-b) <= 2*EPSILON(b) )     THEN
               AlmostEqual = .TRUE.
            ELSE
               AlmostEqual = .FALSE.
            END IF
         ELSE
            IF( ABS( b - a ) <= 2*EPSILON(b)*MAX(ABS(a), ABS(b)) )     THEN
               AlmostEqual = .TRUE.
            ELSE
               AlmostEqual = .FALSE.
            END IF
         END IF
   
      END FUNCTION AlmostEqual

      pure function ThirdDegreeRoots(a,b,c) result (val)
!
!           ----------------------------------------------------------
!              Solves the equation x^3 + a x^2 + b x + c = 0
!              Just one physical (real) solution is expected.
!           ----------------------------------------------------------
!
            implicit none  
            real(kind=RP), intent(in) :: a
            real(kind=RP), intent(in) :: b
            real(kind=RP), intent(in) :: c
            real(kind=RP)             :: val
!           ------------------------------------------------------
            real(kind=RP)              :: p , q
            complex(kind=CP)           :: z1


            p = b - a * a / 3.0_RP
            q = 2.0_RP * a * a * a / 27.0_RP - a * b / 3.0_RP + c
         
            z1 = (( -q + sqrt(q*q + 4.0_RP * p * p * p / 27.0_RP))/2.0_RP)**(1.0_RP / 3.0_RP)         

            val = real(z1 - p / (3.0_RP * z1) - a / 3.0_RP , kind=RP)

         end function ThirdDegreeRoots

         pure function solveTwoEquationLinearSystem( A , b ) result ( x )
            implicit none
            real(kind = RP),  intent(in)     :: A(2,2)
            real(kind = RP),  intent(in)     :: b(2)
            real(kind = RP)                  :: x(2)
!           -----------------------------------------------------
            real(kind = RP)                  :: invdetA

            invdetA = 1.0_RP  / ( A(1,1) * A(2,2) - A(1,2) * A(2,1) )
            x(1) = (b(1) * A(2,2) - b(2) * A(1,2) ) * invdetA
            x(2) = (b(2) * A(1,1) - b(1) * A(2,1) ) * invdetA

         end function solveTwoEquationLinearSystem
      
         pure subroutine solveSecondDegreeEquation(a,b,c,flag , x )
            implicit none
            real(kind=RP), intent(in)            :: a, b, c
            integer,       intent(out), optional :: flag
            real(kind=RP), intent(out)           :: x(2)
!           ---------------------------------------------------
            real(kind=RP)                    :: disc

            if ( almostEqual ( a , 0.0_RP ) ) then
               x = -c / b

               if ( present ( flag ) ) then
                  flag = 1

               end if

            else
               disc = b*b - 4.0_RP * a * c

               if ( disc .ge. 0.0_RP ) then
                  x(1) = (-b + sqrt(disc)) / (2.0_RP * a)
                  x(2) = (-b - sqrt(disc)) / (2.0_RP * a)

                  if ( present ( flag ) ) then
                     flag = 2

                  end if

               else
                  x = 0.0_RP

                  if ( present ( flag ) ) then
                     flag = 0

                  end if
               end if
            end if
            
         end subroutine solveSecondDegreeEquation
   
         pure function cross ( x , y ) result ( val )
            use Physics
            implicit none
            real(kind=RP), intent(in) :: x(NDIM)
            real(kind=RP), intent(in) :: y(NDIM)
            real(kind=RP)             :: val
!           -----------------------------------
            
            val = x(IX) * y(IY) - x(IY) * y(IX)

         end function cross

         pure function newSign ( x ) result ( val ) 
            implicit none
            real(kind=RP), intent(in) :: x
            integer                   :: val
!           ------------------------------------
            integer           :: gt0 , lt0

            lt0 = (x .lt. 0)
            gt0 = (x .gt. 0)

            val = gt0 - lt0

         end function newSign

end module Utilities


