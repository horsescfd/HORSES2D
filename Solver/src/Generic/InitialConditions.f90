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
   type Parameters_t
      real(kind=RP), pointer     :: a0          ! IC Uniform value
      real(kind=RP)              :: thetaIC     ! IC phase [rad]
      real(kind=RP)              :: knIC        ! IC wavenumber
      real(kind=RP), pointer     :: ampIC       ! IC amplitude
   end type Parameters_t
   
   

   type(Parameters_t)            :: params
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

         val = params % a0

      end function UniformInitialCondition

      function WaveInitialCondition(x) result(val)
!        ************************************************
!           Loads a sinusoidal initial condition
!              u = a0 + sin(k_n x + theta)
!                 · k_n is a wavenumber k_n = pi n / T
!                 · theta is the initial phase
!        ************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x
         real(kind=RP)        :: val(NEC)

         val = params % a0 + params % ampIC * sin( params % knIC * x + params % thetaIC )

      end function WaveInitialCondition
   
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
               params % a0 => Setup % a0
               fcn => UniformInitialCondition
!
!           =========================
            case ( trim(waveIC) ) 
!           =========================
!
               params % a0      => Setup % a0
               params % thetaIC =  Setup % thetaIC * pi / 180.0_RP
               params % knIC    =  Setup % nIC * pi / Setup % T
               params % ampIC   => Setup % ampIC
               fcn => WaveInitialCondition
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
