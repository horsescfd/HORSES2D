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

   Name = "Shock wave"

end function getProblemFileName

function HasAnalyticalSolution() 
   implicit none
   logical        :: hasAnalyticalSolution

   hasAnalyticalSolution = .false.

end function HasAnalyticalSolution

function UserDefinedInitialCondition(x , thermodynamics_ , Setup_ , refValues_ , dimensionless_ )  result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)           :: x(NDIM)
   class(Thermodynamics_t), intent(in)  :: thermodynamics_
   class(Setup_t),          intent(in)  :: Setup_
   class(RefValues_t),      intent(in)  :: refValues_
   class(Dimensionless_t),  intent(in)  :: dimensionless_
   real(kind=RP)           :: val(NCONS)

    if ( x(IX) .lt. 0.0_RP ) then
         val = [0.537037037037037_RP,   0.829539386177859_RP,                   0.0_RP,   1.657627118644068_RP]
    else
         val = [1.000000000000000_RP,   0.829539386177859_RP,                   0.0_RP,   2.844067796610170_RP]

    end if

    val(IRHO) = val(IRHO) * refValues_ % rho
    val(IRHOU) = val(IRHOU) * refValues_ % rho * refValues_ % a
    val(IRHOV) = val(IRHOV) * refValues_ % rho * refValues_ % a
    val(IRHOE) = val(IRHOE) * refValues_ % p

end function UserDefinedInitialCondition

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

