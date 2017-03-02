function UserDefinedInitialCondition(x , argin) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)           :: x(NDIM)
   real(kind=RP), optional :: argin
   real(kind=RP)           :: val(NCONS)


   associate ( gamma => thermodynamics % gamma , Mach => dimensionless % Mach ) 

   val(IRHOU ) = sqrt(gamma) * Mach

   if ( x(IY) .lt. 0.0_RP ) then
      val(IRHO) = 1.0_RP
      val(IRHOE) = dimensionless % cv * 1.0_RP + 0.5_RP * val(IRHO) * val(IRHOU) * val(IRHOU) 
   else
      val(IRHO) = 0.9_RP
      val(IRHOE) = dimensionless % cv * 1.0_RP + 0.5_RP * val(IRHO) * val(IRHOU) * val(IRHOU) 
   
   end if

   end associate
      

end function UserDefinedInitialCondition


