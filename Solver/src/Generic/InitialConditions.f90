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

   integer, parameter            :: STR_LEN_IC  = 128
   character(len = *), parameter :: ConstantIC    = "Uniform"
   character(len = *), parameter :: SteadyIC      = "Steady"
   character(len = *), parameter :: VortexIC      = "Vortex transport"
   character(len = *), parameter :: UserDefinedIC = "UserDefined"
   character(len = *), parameter :: ChecksIC      = "Checks"

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

      subroutine InitialCondition( fcn , which )
         use SMConstants
         use Setup_class
         implicit none
         procedure(ICFcn), pointer     :: fcn
         character(len=*), optional    :: which
         character(len=STR_LEN_IC)     :: ICName
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
         
         if (present(which)) then
            ICName = which
         else
            ICName = Setup % IC
         end if

!
!        ********************************
         select case ( trim(ICName) )
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
!           ============================
            case ( trim(ChecksIC) )
!           ============================
!
               fcn => ChecksInitialCondition
!
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

         associate ( gamma => Thermodynamics % gamma ) 

         val(IRHO)  = 1.0_RP
         val(IRHOU) = sqrt( gamma ) * Dimensionless % Mach * cos ( AngleOfAttack )
         val(IRHOV) = sqrt( gamma ) * Dimensionless % Mach * sin ( AngleOfAttack )
         val(IRHOE) = 0.5_RP * gamma * Dimensionless % Mach * Dimensionless % Mach

         end associate

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
         val(IRHOE) = 0.0_RP

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

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         r2 = ((x(iX) - XC)*(x(iX) - XC) + (x(iY) - YC)*(x(iY) - YC)) / (R*R)
      
         u = sqrt(gamma) * Mach * (cos(AngleOfAttack) - Beta * (x(iY) - YC) / R * exp(-0.5_RP * r2))
         v = sqrt(gamma) * Mach * (sin(AngleOfAttack) + Beta * (x(iX) - XC) / R * exp(-0.5_RP * r2))
         T = 1.0_RP - gamma * Mach * Mach * beta * beta / (2.0_RP * Dimensionless % cp) * exp(-r2)
         rho = T**( Thermodynamics % invgm1 ) 

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = Dimensionless % cv * (rho * T - 1.0_RP ) + 0.5_RP * rho * ( u*u + v*v )

         end associate
            
      end function InviscidVortexTransportInitialCondition

      function ChecksInitialCondition (x) result (val)
         use SMConstants
         implicit none
         real(kind=RP)              :: x(NDIM)
         real(kind=RP)              :: val(NEC)

!         val(IRHO)   = x(iX) ** 4.0_RP
!         val(IRHOU)  = x(iY) ** 4.0_RP
!         val(IRHOV)  = (x(iX)** 2.0_RP + x(iY)**2.0_RP)
!         val(IRHOE)  = x(iX) ** 4.0_RP * x(iY) ** 4.0_RP
          val(IRHO) = x(iX)
          val(IRHOU) = x(iY)
          val(IRHOV) = x(iX) + x(iY)
          val(IRHOE) = x(iX) - x(iY)
         
      end function ChecksInitialCondition

end module InitialConditions
