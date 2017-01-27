module InitialConditions
   use SMConstants
   use Physics
   implicit none
!
!  *******
   private
   public ICFcn , InitialCondition , InitialCondition_Describe
!  *******
!

   integer,            parameter :: STR_LEN_IC     = 128
   character(len = *), parameter :: ConstantIC     = "Uniform"
   character(len = *), parameter :: PerturbationIC = "Perturbation"
   character(len = *), parameter :: SteadyIC       = "Steady"
   character(len = *), parameter :: VortexIC       = "Vortex transport"
   character(len = *), parameter :: TaylorIC       = "Taylor vortex"
   character(len = *), parameter :: TurbulenceIC   = "Turbulence2D"
   character(len = *), parameter :: UserDefinedIC  = "UserDefined"
   character(len = *), parameter :: RestartIC      = "Restart"
   character(len = *), parameter :: ChecksIC       = "Checks"
   character(len = *), parameter :: WaveIC         = "Wave"

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
!           =========================
            case ( trim(PerturbationIC) ) 
!           =========================
!
               fcn => PerturbedUniformInitialCondition

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
!           ============================
            case ( trim(TaylorIC) )
!           ============================
!
               fcn => TaylorVortexInitialCondition
!
!           ============================
            case ( trim(TurbulenceIC) )
!           ============================
!
               fcn => Turbulence2DInitialCondition


!
!           ============================
            case ( trim(WaveIC) )
!           ============================
!
               fcn => WaveInitialCondition

!
!           ============================
            case ( trim(RestartIC) )
!           ============================
!
               fcn => NULL()
!
!
!
!           ============               
            case default
!           ============               
!
               print*, 'Error: Unknown initial condition "',trim(ICName),'".'
               print*, "   Options available:"
               print*, "      * ",trim(ConstantIC)
               print*, "      * ",trim(SteadyIC)
               print*, "      * ",trim(VortexIC)
               print*, "      * ",trim(ChecksIC)
               print*, "      * ",trim(RestartIC) 
               print*, "      * ",trim(PerturbationIC) 
               STOP "Stopped." 
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
         val(IRHOE) = Dimensionless % cv + 0.5_RP * gamma * Dimensionless % Mach * Dimensionless % Mach

         end associate

      end function UniformInitialCondition

      function PerturbedUniformInitialCondition(x) result(val)
!        ***************************************************************
!           Loads an uniform initial condition from reference values
!        ***************************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NEC)
         real(kind=RP), parameter   :: AngleOfAttack = 0.0_RP
         real(kind=RP), parameter   :: perturbation = 0.1_RP
         real(kind=RP)              :: per

         associate ( gamma => Thermodynamics % gamma ) 
         per = 0.0_RP
         if ( (x(IX) .gt. 0.4_RP) .and. (x(IX) .lt. 0.6_RP) ) then
            if ( (x(IY) .gt. 0.4_RP) .and. (x(IY) .lt. 0.6_RP) ) then
               per = perturbation
            end if
         end if

         val(IRHO)  = 1.0_RP
         val(IRHOU) = sqrt( gamma ) * Dimensionless % Mach * (cos ( AngleOfAttack ) + per  )
         val(IRHOV) = sqrt( gamma ) * Dimensionless % Mach * sin ( AngleOfAttack )
         val(IRHOE) = Dimensionless % cv + 0.5_RP * gamma * Dimensionless % Mach * Dimensionless % Mach
         end associate

      end function PerturbedUniformInitialCondition



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
         real(kind=RP), parameter      :: R = 0.1_RP                         ! This is the "vortex radius" (dimensionless)
         real(kind=RP), parameter      :: Beta = 0.01_RP                      ! This is the "vortex strength"
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
         val(IRHOE) = Dimensionless % cv * (rho * T) + 0.5_RP * rho * ( u*u + v*v )

         end associate
            
      end function InviscidVortexTransportInitialCondition

      function TaylorVortexInitialCondition(x) result(val)
!        ****************************************************************
!           Loads an initial condition consistent with the Taylor Green 
!          vortex problem
!        ****************************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NEC)
         real(kind=RP)        :: rho , u , v , p

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         rho = 1.0_RP
         u = sqrt(gamma) * Mach * sin(PI*x(1)) * cos(PI*x(2))
         v = -sqrt(gamma) * Mach *  cos(PI*x(1)) * sin(PI*x(2))
         p = 1.0_RP + 0.125_RP * gamma * Mach * Mach * (cos(2.0_RP * PI * x(1)) + cos(2.0_RP * PI * x(2)) )

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
            
      end function TaylorVortexInitialCondition

      function Turbulence2DInitialCondition(x) result(val)
!        ****************************************************************
!           Loads an initial condition to generate 2D Turbulence
!        ****************************************************************
         use SMConstants
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NEC)
         real(kind=RP)        :: rho , u , v , p
         real(kind=RP), parameter      :: k = 10.0_RP

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         rho = 1.0_RP
         u = -sqrt(gamma) * Mach * cos(2.0_RP * k * PI*x(1)) * sin(2.0_RP * k * PI*x(2))
         v = sqrt(gamma) * Mach *  sin(2.0_RP * k * PI*x(1)) * cos(2.0_RP * k * PI*x(2))
         p = 1.0_RP 

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
            
      end function Turbulence2DInitialCondition

      function ChecksInitialCondition (x) result (val)
         use SMConstants
         use Physics
         implicit none
         real(kind=RP)              :: x(NDIM)
         real(kind=RP)              :: val(NEC)
         real(kind=RP)        :: rho , u , v , p

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         rho = 1.0_RP
         u = sqrt(gamma) * Mach * sin(PI*x(1)) * cos(PI*x(2))
         v = -sqrt(gamma) * Mach *  cos(PI*x(1)) * sin(PI*x(2))
         p = 1.0_RP + 0.125_RP * gamma * Mach * Mach * (cos(2.0_RP * PI * x(1)) + cos(2.0_RP * PI * x(2)) )

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
 
         
      end function ChecksInitialCondition

      function WaveInitialCondition ( x ) result( val )
         use SMConstants
         implicit none
         real(kind=RP)              :: x(NDIM)
         real(kind=RP)              :: val(NEC)

         associate ( gamma => Thermodynamics % Gamma  , Mach => Dimensionless % Mach ) 
         val(IRHO) = 1.0_RP
         val(IRHOU) = sqrt(gamma) * Mach * ( 1.0_RP + 0.05_RP * sin(2.0_RP * pi * x(1) ) ) 
         val(IRHOV) = 0.0_RP
         val(IRHOE) = Dimensionless % cv + 0.5_RP * 1.0_RP * ( val(IRHOU)**2.0_RP )
         end associate
      end function WaveInitialCondition
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        DESCRIBE
!        --------
!/////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine InitialCondition_Describe()
         use Headers
         use Setup_class
         implicit none
   
         write(STD_OUT ,'(/)')
         call Section_header("Initial condition description")
         write(STD_OUT ,*)
         
         write(STD_OUT , '(30X,A,A25,A,A)') "-> ", "Initial condition type: " , trim(Setup % IC) , "."
   
         if ( trim(Setup % IC) .eq. "Restart" ) then
            write(STD_OUT , '(30X,A,A25,A,A)') "-> ", "Data restart from file: " , trim(Setup % restart_file) , "."
         end if

         write(STD_OUT , '(30X,A,A25,ES10.4,A)' ) "-> " , "Initial time: " , Setup % initialTime

      end subroutine InitialCondition_Describe

end module InitialConditions
