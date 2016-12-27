module InitialConditions
   use SMConstants
   use Physics
   implicit none
!
!  *******
   private
   public ICFcn , InitialCondition
!  *******
!
   character(len = *), parameter :: ConstantIC    = "Uniform"
   character(len = *), parameter :: WaveIC        = "Wave"
   character(len = *), parameter :: UserDefinedIC = "UserDefined"
! 
!  ******************
   abstract interface
!  ******************
!
      function ICFcn (x) result (val)
         use SMConstants
         use Physics
         implicit none
         real(kind=RP)     :: x
         real(kind=RP)     :: val(NEC)
      end function ICFcn
!
!  *************
   end interface
!  *************
!
!
!  ========   
   contains
!  ========   
!
      function UniformInitialCondition(x) result(val)
!        ******************************************
!           Loads an uniform initial condition
!        ******************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x
         real(kind=RP)        :: val(NEC)

         val = 1.0_RP

      end function UniformInitialCondition

      subroutine InitialCondition( fcn )
         use SMConstants
         use Setup_class
         implicit none
         procedure(ICFcn), pointer     :: fcn
         interface
            function UserDefinedInitialCondition(x) result (val)
               use SMConstants
               use Setup_class
               use Physics
               implicit none
               real(kind=RP)     :: x
               real(kind=RP)     :: val(NEC)
            end function UserDefinedInitialCondition
         end interface
!
!        ********************************
         select case ( trim(Setup % IC) )
!        ********************************
!
!           =========================
            case ( trim(ConstantIC) ) 
!           =========================
!
               fcn => UniformInitialCondition
!
!           =========================
            case ( trim(waveIC) ) 
!           =========================
!
!
!
!           ============================
            case ( trim(UserDefinedIC) )
!           ============================
!
               fcn => UserDefinedInitialCondition
!
!           ============               
            case default
!           ============               
!
               STOP "Unknown Initial condition" 
!
!        **********
         end select
!        **********
!


      end subroutine InitialCondition
            

         

end module InitialConditions
