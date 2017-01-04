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
   character(len = *), parameter :: SteadyIC      = "Steady"
   character(len = *), parameter :: VortexIC      = "Vortex transport"
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
         real(kind=RP)     :: x(NDIM)
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
               real(kind=RP)     :: x(NDIM)
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
!           ============================
            case ( trim(SteadyIC) )
!           ============================
!
               fcn => SteadyInitialCondition
!
!           ============================
            case ( trim(VortexIC) )
!           ============================
!
               fcn => InviscidVortexTransportInitialCondition
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
            
      function UniformInitialCondition(x) result(val)
!        ***************************************************************
!           Loads an uniform initial condition from reference values
!        ***************************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NEC)
         real(kind=RP), parameter   :: AngleOfAttack = 0.0_RP

         val(IRHO)  = 1.0_RP
         val(IRHOU) = Dimensionless % Mach * cos ( AngleOfAttack)
         val(IRHOV) = Dimensionless % Mach * sin ( AngleOfAttack )
         val(IRHOE) = Dimensionless % cv + 0.5_RP * Thermodynamics % Gamma * Dimensionless % Mach * Dimensionless % Mach

      end function UniformInitialCondition

      function SteadyInitialCondition(x) result(val)
!        ***************************************************************
!           Loads a steady flow from the reference values
!        ***************************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NEC)

         val(IRHO) = 1.0_RP
         val(IRHOU) = 0.0_RP
         val(IRHOV) = 0.0_RP
         val(IRHOE) = Dimensionless % cv

      end function SteadyInitialCondition
         
      function InviscidVortexTransportInitialCondition(x) result(val)
!        ****************************************************************
!           Loads an initial condition consistent with the Taylor Green 
!          vortex problem
!        ****************************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NEC)
         real(kind=RP), parameter      :: R = 0.25_RP                         ! This is the "vortex radius" (dimensionless)
         real(kind=RP), parameter      :: Beta = 0.25_RP                      ! This is the "vortex strength"
         real(kind=RP), parameter      :: XC = 0.5_RP                         ! Vortex X position (in dimensionless coordinates)
         real(kind=RP), parameter      :: YC = 0.5_RP                         ! Vortex Y position (in dimensionless coordinates)
         real(kind=RP), parameter      :: AngleOfAttack = 0.0_RP
         real(kind=RP)                 :: r2 , rho , u , v , T

         r2 = ((x(iX) - XC)*(x(iX) - XC) + (x(iY) - YC)*(x(iY) - YC)) / (R*R)
      
         u = Dimensionless % Mach * (cos(AngleOfAttack) - Beta * (x(iY) - YC) / R )* exp(-0.5_RP * r2)
         v = Dimensionless % Mach * (sin(AngleOfAttack) + Beta * (x(iX) - XC) / R )* exp(-0.5_RP * r2)
         T = 1.0_RP - Dimensionless % Mach * Dimensionless % Mach * beta * beta / (2.0_RP * Dimensionless % cv) * exp(-r2)
         rho = T**( Thermodynamics % invgm1 ) 

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = Dimensionless % cv * rho * T + 0.5_RP * Thermodynamics % Gamma * rho * ( u*u + v*v )
            

         
      end function InviscidVortexTransportInitialCondition

end module InitialConditions
