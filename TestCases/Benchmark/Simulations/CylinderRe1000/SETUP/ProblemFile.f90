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

   Name = "Cylinder Re 1000 IP problem file"

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

   val(IRHO)  = refValues_ % rho
   val(IRHOU) = refValues_ % rho * refValues_ % V * cos(AngleOfAttack)
   val(IRHOV) = refValues_ % rho * refValues_ % V * sin(AngleOfAttack)
   val(IRHOE) = thermodynamics_ % cv * refValues_ % rho * refValues_ % T + 0.5_RP * refValues_ % rho * refValues_ % V * refValues_ % V

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
!   ------------------------------------------------------------------------------
!   This subroutine checks that the residuals, lift, and drag, match stored values
!   ------------------------------------------------------------------------------
!
    real(kind=RP), parameter :: residuals(4) = [1.0120032175167813E+02_RP , 1.1656904056056065E+02_RP , 5.7372439386641894E+01_RP , 3.5383131992651778E+02_RP]
    real(kind=RP), parameter :: cl = -3.2989364784973684E-08_RP
    real(kind=RP), parameter :: cd =  2.8170926604773996E+01_RP
    integer,       parameter :: finalIteration = 512

    write(STD_OUT,*)
    write(STD_OUT,'(30X,A,A35,A,ES10.3)') "-> ", "Error found in continuity residuals" , " : " , abs(residuals(1) - Monitors_ % residuals % values(IRHO ,1)) 
    write(STD_OUT,'(30X,A,A35,A,ES10.3)') "-> ", "Error found in x-momentum residuals" , " : " , abs(residuals(2) - Monitors_ % residuals % values(IRHOU,1)) 
    write(STD_OUT,'(30X,A,A35,A,ES10.3)') "-> ", "Error found in y-momentum residuals" , " : " , abs(residuals(3) - Monitors_ % residuals % values(IRHOV,1)) 
    write(STD_OUT,'(30X,A,A35,A,ES10.3)') "-> ", "Error found in energy residuals" , " : " , abs(residuals(4) - Monitors_ % residuals % values(IRHOE,1)) 
    write(STD_OUT,'(30X,A,A35,A,ES10.3)') "-> ", "Error found in lift coefficient" , " : " , abs(cl - Monitors_ % surfaceMonitors(1) % values(1))
    write(STD_OUT,'(30X,A,A35,A,ES10.3)') "-> ", "Error found in drag coefficient" , " : " , abs(cd - Monitors_ % surfaceMonitors(2) % values(1)) 
    write(STD_OUT , '(    30X , A , A35 , A,I0,A,I0)') "-> " , "Error in the final iteration" , " : "      , finalIteration , "/",sem_ % Integrator % iter

    if ( (maxval(abs(residuals - Monitors_ % residuals % values(:,1))) .gt. 1.0e-10_RP) .or. ( finalIteration .ne. sem_ % Integrator % iter ) &
                           .or. ( abs(cl-Monitors_ % surfaceMonitors(1) % values(1)) .gt. 1.0e-12 ) .or. &
                             ( abs(cd - Monitors_ % surfaceMonitors(2) % values(1) ) .gt. 1.0e-12 ) ) then
      write(STD_OUT , '(    30X , A , A35,A,A)') "-> " , "Cylinder Re 1000 IP benchmark test" , " : " ,  "failed."
      exit_code = FAILED
    else
      write(STD_OUT , '(    30X , A , A35,A,A)') "-> " , "Cylinder Re 1000 IP benchmark test" , " : " ,  "succeeded."
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
