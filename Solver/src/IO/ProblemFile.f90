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

function UserDefinedInitialCondition(x , argin , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)                        :: x(NDIM)
   real(kind=RP), optional              :: argin
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
   real(kind=RP)            :: pressure
   real(kind=RP), parameter :: AngleOfAttack = 0.0_RP


   pressure = Setup_ % pressure_ref / RefValues_ % p

   associate ( gamma => Thermodynamics_ % gamma , Mach => dimensionless_ % Mach ) 

#ifdef _DIMENSIONLESS_TAU
   val(IRHO)  = 1.0_RP
   val(IRHOU) = sqrt(gamma) * Mach * cos ( AngleOfAttack )
   val(IRHOV) = sqrt(gamma) * Mach * sin ( AngleOfAttack )
   val(IRHOE) = Dimensionless_ % cv * pressure + 0.5_RP * gamma * Mach * Mach
#else
   val(IRHO)  = 1.0_RP
   val(IRHOU) = cos ( AngleOfAttack )
   val(IRHOV) = sin ( AngleOfAttack )
   val(IRHOE) = Dimensionless_ % cv * pressure + 0.5_RP 
#endif

   end associate


end function UserDefinedInitialCondition


