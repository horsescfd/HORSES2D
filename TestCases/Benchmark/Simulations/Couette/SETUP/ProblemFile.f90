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

   Name = "Couette flow problem file"

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
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP
   real(kind=RP)            :: T

   T = refValues_ % T - (refValues_ % V*refValues_ % V) * dimensionless_ % Pr / (8.0_RP * Thermodynamics_ % cp) * (x(IY)+1.0_RP)*(x(IY)-1.0_RP)

   val(IRHO)  = Setup_ % pressure_ref / ( Thermodynamics_ % R * T )
   val(IRHOU) = val(IRHO) * refValues_ % V * 0.5_RP * ( x(IY) + 1.0_RP ) 
   val(IRHOV) = 0.0_RP
   val(IRHOE) = thermodynamics_ % cv * val(IRHO) * T + 0.5_RP * val(IRHOU) * val(IRHOU)/ val(IRHO)

   if ( (x(IX)*x(IX) + x(IY)*x(IY)) .lt. 0.005_RP ) then
      val(IRHOU) = 1.005_RP * val(IRHOU) 
   end if

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
    real(kind=RP)            :: errors(NCONS) , localErrors(NCONS)
    real(kind=RP)            :: analytical(NCONS)
    real(kind=RP)            :: x(NDIM)
    integer, parameter       :: finalIteration = 19585
    real(kind=RP), parameter :: expectedResiduals(NCONS) &
            = [2.9153825904001754E-05_RP , 1.0168055570704469E-05_RP , 2.9823607695025403E-05_RP , 9.9806687831331889E-05_RP]
    interface
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
      end function AnalyticalSolution
    end interface
      
      
    errors = 0.0_RP

    do eID = 1 , sem_ % mesh % no_of_elements
      do i = 0 , sem_ % mesh % elements(eID) % spA % N
         do j = 0 , sem_ % mesh % elements(eID) % spA % N

            x = sem_ % mesh % elements(eID) % x(i,j,1:NDIM)
            analytical = AnalyticalSolution( x , 0.0_RP , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) 

            analytical(IRHO) = analytical(IRHO) / refValues_ % rho
            analytical(IRHOU) = analytical(IRHOU) / ( refValues_ % rho * refValues_ % a )
            analytical(IRHOV) = analytical(IRHOV) / ( refValues_ % rho * refValues_ % a )
            analytical(IRHOE) = analytical(IRHOE) / refValues_ % p

            localErrors = abs(analytical - sem_ % mesh % elements(eID) % Q(1:NCONS,i,j) )

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
    write(STD_OUT , '(    30X , A , A50 , I0,A,I0)') "-> " , "Error in the final iteration: "     , finalIteration , "/",sem_ % Integrator % iter
    write(STD_OUT,'(30X,A,A50,A,ES10.3)') "-> ", "Error found in continuity residuals" , " : " , abs(expectedResiduals(1) - Monitors_ % residuals % values(IRHO ,1)) 
    write(STD_OUT,'(30X,A,A50,A,ES10.3)') "-> ", "Error found in x-momentum residuals" , " : " , abs(expectedResiduals(2) - Monitors_ % residuals % values(IRHOU,1)) 
    write(STD_OUT,'(30X,A,A50,A,ES10.3)') "-> ", "Error found in y-momentum residuals" , " : " , abs(expectedResiduals(3) - Monitors_ % residuals % values(IRHOV,1)) 
    write(STD_OUT,'(30X,A,A50,A,ES10.3)') "-> ", "Error found in energy residuals" , " : " , abs(expectedResiduals(4) - Monitors_ % residuals % values(IRHOE,1))

    if ( (maxval(abs(expectedResiduals - Monitors_ % residuals % values(:,1))) .gt. 1.0e-12_RP) .or. ( finalIteration .ne. sem_ % Integrator % iter ) ) then
      write(STD_OUT , '(    30X , A , A50,A,A)') "-> " , "Couette benchmark test" , " : " , "failed."
      exit_code = FAILED
    else
      write(STD_OUT , '(    30X , A , A50,A,A)') "-> " , "Couette benchmark test" , " : " ,  "succeeded."
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
   real(kind=RP)            :: T

   T = refValues_ % T - (refValues_ % V*refValues_ % V) * dimensionless_ % Pr / (8.0_RP * Thermodynamics_ % cp) * (x(IY)+1.0_RP)*(x(IY)-1.0_RP)

   state(IRHO)  = Setup_ % pressure_ref / ( Thermodynamics_ % R * T )
   state(IRHOU) = state(IRHO) * refValues_ % V * 0.5_RP * ( x(IY) + 1.0_RP ) 
   state(IRHOV) = 0.0_RP
   state(IRHOE) = thermodynamics_ % cv * state(IRHO) * T + 0.5_RP * state(IRHOU) * state(IRHOU)/ state(IRHO)

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
!
!  ---------------
!  Local variables
!  ---------------
!
   real(kind=RP)            :: T

   T = refValues_ % T - (refValues_ % V*refValues_ % V) * dimensionless_ % Pr / (8.0_RP * Thermodynamics_ % cp) * (x(IY)+1.0_RP)*(x(IY)-1.0_RP)

   val(IRHO)  = Setup_ % pressure_ref / ( Thermodynamics_ % R * T )
   val(IRHOU) = val(IRHO) * refValues_ % V * 0.5_RP * ( x(IY) + 1.0_RP ) 
   val(IRHOV) = 0.0_RP
   val(IRHOE) = thermodynamics_ % cv * val(IRHO) * T + 0.5_RP * val(IRHOU) * val(IRHOU)/ val(IRHO)

end function AnalyticalSolution
