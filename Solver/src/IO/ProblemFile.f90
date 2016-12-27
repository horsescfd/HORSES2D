function UserDefinedInitialCondition(x) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)     :: x
   real(kind=RP)     :: val(NEC)

   val = 1.0_RP

end function UserDefinedInitialCondition


