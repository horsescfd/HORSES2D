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


