function UserDefinedInitialCondition(x) result (val)
   use SMConstants
   use Setup_class
   use Physics
   implicit none
   real(kind=RP)     :: x
   real(kind=RP)     :: val(NEC)

   val = sin(pi * x / Setup % T)

end function UserDefinedInitialCondition


