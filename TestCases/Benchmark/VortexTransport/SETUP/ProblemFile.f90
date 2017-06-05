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
#include "Defines.h"

function getProblemFileName() result ( Name )
   implicit none
   character(len=LINE_LENGTH)    :: Name

   Name = "Inviscid taylor vortex transport"

end function getProblemFileName

function UserDefinedInitialCondition(x , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)                        :: x(NDIM)
   class(Thermodynamics_t), intent(in)  :: thermodynamics_
   class(Setup_t),          intent(in)  :: Setup_
   class(RefValues_t),      intent(in)  :: refValues_
   class(Dimensionless_t),  intent(in)  :: dimensionless_
   real(kind=RP)                        :: val(NCONS)
!
!  ---------------
!  Local variables
!  ---------------
!
   real(kind=RP), parameter :: R = 0.1_RP                  ! This is the "vortex radius" (dimensionless)
   real(kind=RP), parameter :: Beta = 0.05_RP              ! This is the "vortex strength"
   real(kind=RP), parameter :: XC = 0.0_RP                 ! Vortex X position (in dimensionless coordinates)
   real(kind=RP), parameter :: YC = 0.0_RP                 ! Vortex Y position (in dimensionless coordinates)
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP
   real(kind=RP)            :: r2 , rho , u , v , T

   associate ( gamma => Thermodynamics_ % Gamma , Mach => Dimensionless_ % Mach , cv => Dimensionless_ % cv )

   r2 = ((x(IX) - XC)*(x(IX) - XC) + (x(IY) - YC)*(x(IY) - YC)) / (R*R)

   u = refValues_ % V * (cos(AngleOfAttack) - Beta * (x(IY) - YC) / R * exp(-0.5_RP * r2))
   v = refValues_ % V * (sin(AngleOfAttack) + Beta * (x(IX) - XC) / R * exp(-0.5_RP * r2))
   T = refValues_ % T * (1.0_RP - gamma * Mach * Mach * beta * beta / (2.0_RP * Dimensionless_ % cp) * exp(-r2) )
   rho = refValues_ % rho * (T/refValues_ % T)**( Thermodynamics_ % invgm1 ) 

   val(IRHO)  = rho
   val(IRHOU) = rho * u
   val(IRHOV) = rho * v
   val(IRHOE) = thermodynamics_ % cv * (rho * T) + 0.5_RP * rho * ( u*u + v*v )

   end associate

end function UserDefinedInitialCondition

function Finalize( sem_ , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ , Monitors_) result(exit_code)
    use SMConstants
    use DGSEM_Class
    use Setup_class
    use QuadMeshClass
    use QuadElementClass
    use Headers
    use Physics
    use MonitorsClass
    implicit none
    class(DGSEM_t)                      :: sem_
    class(Thermodynamics_t), intent(in) :: thermodynamics_
    class(Setup_t),          intent(in) :: Setup_
    class(RefValues_t),      intent(in) :: refValues_
    class(Dimensionless_t),  intent(in) :: dimensionless_
    class(Monitor_t),       intent(in)  :: Monitors_
    integer                             :: exit_code
!
!   ************************************************
!   Compute the error w.r.t. the analytical solution
!   ************************************************
!
    integer       :: i , j  , eID
    real(kind=RP)    :: errors(NCONS) , localErrors(NCONS)
    real(kind=RP)    :: analytical(NCONS)
    real(kind=RP)    :: x(NDIM)
    integer,       parameter  :: expectedIter = 3237
    real(kind=RP), parameter  :: expectedErrors(NCONS) = &
            [ 1.7564464390384948E-005_RP, 1.3312909555291963E-003_RP, 3.8656265138769554E-004_RP, 9.5352754637234582E-004_RP ]
    interface
      function UserDefinedInitialCondition(x, Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (val)
         use SMConstants
         use Setup_class
         use Physics
         implicit none
         real(kind=RP),           intent(in)           :: x(NDIM)
         class(Thermodynamics_t), intent(in)           :: thermodynamics_
         class(Setup_t),          intent(in)           :: Setup_
         class(RefValues_t),      intent(in)           :: refValues_
         class(Dimensionless_t),  intent(in)           :: dimensionless_
         real(kind=RP)                                 :: val(NCONS)
      end function UserDefinedInitialCondition
    end interface

    errors = 0.0_RP

    do eID = 1 , sem_ % mesh % no_of_elements
      do i = 0 , sem_ % mesh % elements(eID) % spA % N
         do j = 0 , sem_ % mesh % elements(eID) % spA % N

            x = sem_ % mesh % elements(eID) % x(i,j,1:NDIM)
            analytical = UserDefinedInitialCondition( x , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) 

            analytical(IRHO) = analytical(IRHO) / refValues_ % rho
            analytical(IRHOU) = analytical(IRHOU) / ( refValues_ % rho * refValues_ % a )
            analytical(IRHOV) = analytical(IRHOV) / ( refValues_ % rho * refValues_ % a )
            analytical(IRHOE) = analytical(IRHOE) / refValues_ % p

            localErrors = abs(analytical - sem_ % mesh % elements(eID) % Q(i,j,1:NCONS) )

            if ( localErrors(IRHO) .gt. errors(IRHO) ) then
               errors(IRHO) = localErrors(IRHO)
            end if

            if ( localErrors(IRHOU) .gt. errors(IRHOU) ) then
               errors(IRHOU) = localErrors(IRHOU)
            end if

            if ( localErrors(IRHOV) .gt. errors(IRHOV) ) then
               errors(IRHOV) = localErrors(IRHOV)
            end if

            if ( localErrors(IRHOE) .gt. errors(IRHOE) ) then
               errors(IRHOE) = localErrors(IRHOE)
            end if
   
         end do
      end do
    end do

    write(STD_OUT , '(/ , 30X , A , A50 , ES10.3)') "-> " , "Error w.r.t. analytical solution in density: "    , errors(IRHO)
    write(STD_OUT , '(    30X , A , A50 , ES10.3)') "-> " , "Error w.r.t. analytical solution in x-momentum: " , errors(IRHOU)
    write(STD_OUT , '(    30X , A , A50 , ES10.3)') "-> " , "Error w.r.t. analytical solution in y-momentum: " , errors(IRHOV)
    write(STD_OUT , '(    30X , A , A50 , ES10.3)') "-> " , "Error w.r.t. analytical solution in energy: "     , errors(IRHOE)
    write(STD_OUT , '(    30X , A , A50 , ES10.3)') "-> " , "Difference w.r.t. expected errors: "     , maxval(abs(errors - expectedErrors))
    write(STD_OUT , '(    30X , A , A35 , I0,A,I0)') "-> " , "Error in the final iteration: "     , expectedIter , "/",sem_ % Integrator % iter

    if ( (maxval(abs(errors-expectedErrors)) .gt. 1.0e-12_RP) .or. ( expectedIter .ne. sem_ % Integrator % iter ) ) then
      write(STD_OUT , '(    30X , A , A50,A,A)') "-> " , "Inviscid vortex transport benchmark test" , " : " ,  "failed."
      exit_code = FAILED
    else
      write(STD_OUT , '(    30X , A , A50,A,A)') "-> " , "Inviscid vortex transport benchmark test" , " : " ,  "succeeded."
      exit_code = SUCCESSFUL
    end if
   
     
end function Finalize

function BoundaryConditionFunction1(x,time, Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (state)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP),           intent(in)           :: x(NDIM)
   real(kind=RP),           intent(in)           :: time
   class(Thermodynamics_t), intent(in)           :: thermodynamics_
   class(Setup_t),          intent(in)           :: Setup_
   class(RefValues_t),      intent(in)           :: refValues_
   class(Dimensionless_t),  intent(in)           :: dimensionless_
   real(kind=RP)                                 :: state(NCONS)
!
!  ---------------
!  Local variables
!  ---------------
!
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP

   state(IRHO)  = refValues_ % rho
   state(IRHOU) = refValues_ % rho * refValues_ % V * cos(AngleOfAttack)
   state(IRHOV) = refValues_ % rho * refValues_ % V * sin(AngleOfAttack)
   state(IRHOE) = thermodynamics_ % cv * refValues_ % rho * refValues_ % T + 0.5_RP * refValues_ % rho * refValues_ % V * refValues_ % V

end function BoundaryConditionFunction1

function BoundaryConditionFunction2(x,time, Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (state)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP),           intent(in)           :: x(NDIM)
   real(kind=RP),           intent(in)           :: time
   class(Thermodynamics_t), intent(in)           :: thermodynamics_
   class(Setup_t),          intent(in)           :: Setup_
   class(RefValues_t),      intent(in)           :: refValues_
   class(Dimensionless_t),  intent(in)           :: dimensionless_
   real(kind=RP)                                 :: state(NCONS)
!
!  ---------------
!  Local variables
!  ---------------
!
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP

   state(IRHO)  = refValues_ % rho
   state(IRHOU) = refValues_ % rho * refValues_ % V * cos(AngleOfAttack)
   state(IRHOV) = refValues_ % rho * refValues_ % V * sin(AngleOfAttack)
   state(IRHOE) = thermodynamics_ % cv * refValues_ % rho * refValues_ % T + 0.5_RP * refValues_ % rho * refValues_ % V * refValues_ % V

end function BoundaryConditionFunction2

function BoundaryConditionFunction3(x,time, Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (state)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP),           intent(in)           :: x(NDIM)
   real(kind=RP),           intent(in)           :: time
   class(Thermodynamics_t), intent(in)           :: thermodynamics_
   class(Setup_t),          intent(in)           :: Setup_
   class(RefValues_t),      intent(in)           :: refValues_
   class(Dimensionless_t),  intent(in)           :: dimensionless_
   real(kind=RP)                                 :: state(NCONS)
!
!  ---------------
!  Local variables
!  ---------------
!
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP

   state(IRHO)  = refValues_ % rho
   state(IRHOU) = refValues_ % rho * refValues_ % V * cos(AngleOfAttack)
   state(IRHOV) = refValues_ % rho * refValues_ % V * sin(AngleOfAttack)
   state(IRHOE) = thermodynamics_ % cv * refValues_ % rho * refValues_ % T + 0.5_RP * refValues_ % rho * refValues_ % V * refValues_ % V

end function BoundaryConditionFunction3

function BoundaryConditionFunction4(x,time, Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (state)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP),           intent(in)           :: x(NDIM)
   real(kind=RP),           intent(in)           :: time
   class(Thermodynamics_t), intent(in)           :: thermodynamics_
   class(Setup_t),          intent(in)           :: Setup_
   class(RefValues_t),      intent(in)           :: refValues_
   class(Dimensionless_t),  intent(in)           :: dimensionless_
   real(kind=RP)                                 :: state(NCONS)
!
!  ---------------
!  Local variables
!  ---------------
!
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP

   state(IRHO)  = refValues_ % rho
   state(IRHOU) = refValues_ % rho * refValues_ % V * cos(AngleOfAttack)
   state(IRHOV) = refValues_ % rho * refValues_ % V * sin(AngleOfAttack)
   state(IRHOE) = thermodynamics_ % cv * refValues_ % rho * refValues_ % T + 0.5_RP * refValues_ % rho * refValues_ % V * refValues_ % V

end function BoundaryConditionFunction4

function AnalyticalSolution(x,time, Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP),           intent(in)           :: x(NDIM)
   real(kind=RP),           intent(in)           :: time
   class(Thermodynamics_t), intent(in)           :: thermodynamics_
   class(Setup_t),          intent(in)           :: Setup_
   class(RefValues_t),      intent(in)           :: refValues_
   class(Dimensionless_t),  intent(in)           :: dimensionless_
   real(kind=RP)                                 :: val(NCONS)

   val = 0.0_RP

end function AnalyticalSolution


