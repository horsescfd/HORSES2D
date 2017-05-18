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

function UserDefinedInitialCondition(x , argin) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)           :: x(NDIM)
   real(kind=RP), optional :: argin
   real(kind=RP)           :: val(NCONS)


   associate ( gamma => thermodynamics % gamma , Mach => dimensionless % Mach ) 


!   val(IRHO) = 1.0_RP
!   val(IRHOV) = 0.0_RP
!   val(IRHOE) = dimensionless % cv * 1.0_RP + 0.5_RP * gamma * Mach * Mach
!
!   if ( x(IY) .lt. 0.0_RP ) then
!      val(IRHOU) = sqrt(gamma) * Mach
!   else
!      val(IRHOU) = -sqrt(gamma) * Mach
!   
!   end if
!
!   val(IRHOU) = val(IRHOU) + 1.0e-3_RP * sqrt(gamma) * Mach * exp(-(x(IY) / (28.0_RP))**2.0_RP) * ( cos(2.0_RP * PI * x(IX)) + cos(20.0_RP * PI * x(IX) ) ) 
!

    if ( x(IX) .lt. 0.0_RP ) then
         val = [0.537037037037037_RP,   0.829539386177859_RP,                   0.0_RP,   1.657627118644068_RP]

    else
         val = [1.000000000000000_RP,   0.829539386177859_RP,                   0.0_RP,   2.844067796610170_RP]


    end if
   end associate

end function UserDefinedInitialCondition


