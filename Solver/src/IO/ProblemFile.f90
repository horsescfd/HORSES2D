function UserDefinedInitialCondition(x , argin) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)           :: x(NDIM)
   real(kind=RP), optional :: argin
   real(kind=RP)           :: val(NCONS)

   val = 0.0_RP * x(1) + 0.0_RP * x(2)


end function UserDefinedInitialCondition


