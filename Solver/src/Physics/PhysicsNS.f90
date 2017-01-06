module PhysicsNS
    use SMConstants
    use Setup_class

    private
    public :: NEC , NDIM , IX , IY , IRHO , IRHOU , IRHOV , IRHOE , solver
    public :: RefValues , Dimensionless , Thermodynamics
    public :: RiemannSolverFunction , InitializePhysics
    public :: InviscidFlux
!
!   *****************************************
!        Definitions
!   *****************************************
!
    integer, parameter              :: NEC = 4
    integer, parameter              :: NDIM = 2
    integer, parameter              :: STR_LEN_PHYSICS = 128
!
!   *****************************************
!        Parameter to control dimensions
!   *****************************************
!
    integer, parameter :: IX = 1
    integer, parameter :: IY = 2
!
!   ******************++**********************
!        Parameters to select variables
!   ******************++**********************
!
    integer, parameter :: IRHO  = 1
    integer, parameter :: IRHOU = 2
    integer, parameter :: IRHOV = 3
    integer, parameter :: IRHOE = 4
!
!   ********************************************
!        Current solver
!   ********************************************
!
    character(LEN=*), PARAMETER     :: solver = "Navier-Stokes"

    type Thermodynamics_t  
      character(len=STR_LEN_PHYSICS)       :: fluid
      real(kind=RP)        :: R
      real(kind=RP)        :: gamma       ! Specific heat ratio
      real(kind=RP)        :: gm1         ! (gamma - 1)
      real(kind=RP)        :: gogm1       ! gamma / ( gamma - 1 )
      real(kind=RP)        :: invgm1      ! 1.0_RP / (gamma - 1 )
      real(kind=RP)        :: cp
      real(kind=RP)        :: cv
    end type Thermodynamics_t

    type RefValues_t
      real(kind=RP)        :: L
      real(kind=RP)        :: T
      real(kind=RP)        :: p
      real(kind=RP)        :: rho
      real(kind=RP)        :: V
      real(kind=RP)        :: a
      real(kind=RP)        :: Mach
      real(kind=RP)        :: mu
      real(kind=RP)        :: kappa
      real(kind=RP)        :: tc
    end type RefValues_t

    type Dimensionless_t
      real(kind=RP)        :: mu
      real(kind=RP)        :: cp
      real(kind=RP)        :: cv
      real(kind=RP)        :: Re
      real(kind=RP)        :: Pr
      real(kind=RP)        :: kappa
      real(kind=RP)        :: Mach
    end type Dimensionless_t

    abstract interface
      function RiemannSolverFunction( QL , QR , n ) result ( val )
         use SMConstants
         import NEC , NDIM
         real(kind=RP), dimension(NEC)       :: QL
         real(kind=RP), dimension(NEC)       :: QR
         real(kind=RP), dimension(NDIM)      :: n
         real(kind=RP), dimension(NEC)       :: val
      end function RiemannSolverFunction
    end interface

    type(Thermodynamics_t), target  :: thermodynamicsAir = Thermodynamics_t("Air",287.15_RP , 1.4_RP , 0.4_RP , 3.5_RP , 2.5_RP , 287.0_RP*3.5_RP , 287.0_RP*2.5_RP)
    type(Thermodynamics_t), pointer, protected            :: thermodynamics
    type(RefValues_t), protected       :: refValues      
    type(Dimensionless_t), protected   :: dimensionless

!
    interface inviscidFlux
      module procedure inviscidFlux0D , inviscidFlux1D , inviscidFlux2D
    end interface inviscidFlux
!
!    interface viscousFlux
!      module procedure viscousFlux0D , viscousFlux1D , viscousFlux2D
!    end interface viscousFlux
!
    contains

       subroutine InitializePhysics()
         implicit none
!
!        *******************************
!        Select fluid
!        *******************************
!
         if (trim(Setup % Gas) == "Air") then
            thermodynamics => thermodynamicsAir
         else
            print*, "Unknown fluid ", trim(Setup % Gas),"."
            stop "Stopped."
         end if
!
!        *******************************
!        Initialize the reference values
!        *******************************
!
         refValues % L     = Setup % reynolds_length
         refValues % T     = Setup % temperature_ref
         refValues % p     = Setup % pressure_ref
         refValues % rho   = refValues % p / (thermodynamics % R * refValues % T)
         refValues % a     = sqrt( refValues % p / refValues % rho )
         refValues % Mach  = Setup % Mach_number
         refValues % V     = sqrt( Thermodynamics % gamma ) * refValues % Mach * refValues % a
         refValues % mu    = refValues % rho * refValues % V * refValues % L / Setup % reynolds_number
         refValues % kappa = refValues % mu * thermodynamics % cp / Setup % prandtl_number
         refValues % tc    = refValues % L / (refValues % V )
!
!        ***********************************
!        Initialize the dimensionless values
!        ***********************************
!
         dimensionless % mu         = sqrt( thermodynamics % gamma ) * refValues % Mach / Setup % reynolds_number
         dimensionless % kappa      = dimensionless % mu * thermodynamics % gogm1 / Setup % prandtl_number
         dimensionless % cp         = thermodynamics % gogm1
         dimensionless % cv         = thermodynamics % invgm1
         dimensionless % Re         = Setup % reynolds_number
         dimensionless % Pr         = Setup % prandtl_number
         dimensionless % Mach       = Setup % Mach_number


         call Describe

       end subroutine InitializePhysics
!
      function inviscidFlux0D(u) result(val)
         implicit none
         real(kind=RP)          :: u(NEC)
         real(kind=RP), target  :: val(NEC,NDIM)
         real(kind=RP), pointer :: F(:) , G(:)

         F(1:NEC)    => val(1:NEC,iX)
         G(1:NEC)    => val(1:NEC,iY)

         associate ( Gamma => Thermodynamics % Gamma , gm1 => Thermodynamics % gm1 , Mach => Dimensionless % Mach ) 
    
         F(IRHO)  = u(IRHOU)
         F(IRHOU) = 0.5_RP * (3.0_RP - Gamma) * u(IRHOU)*u(IRHOU) / u(IRHO) - 0.5_RP * gm1 * u(IRHOV)*u(IRHOV)/u(IRHO) + gm1 * u(IRHOE)
         F(IRHOV) = u(IRHOU)*u(IRHOV) / u(IRHO)
         F(IRHOE) = (Gamma * u(IRHOE) - 0.5_RP*gm1*( u(IRHOU)*u(IRHOU) + u(IRHOV)*u(IRHOV) ) / u(IRHO) + Dimensionless % cp) * u(IRHOU) / u(IRHO)

         G(IRHO)  = u(IRHOV)
         G(IRHOU) = u(IRHOU)*u(IRHOV) / u(IRHO)
         G(IRHOV) = 0.5_RP * (3.0_RP - Gamma) * u(IRHOV)*u(IRHOV) / u(IRHO) - 0.5_RP * gm1 * u(IRHOU)*u(IRHOU)/u(IRHO) + gm1 * u(IRHOE)
         G(IRHOE) = (Gamma * u(IRHOE) - 0.5_RP*gm1*( u(IRHOU)*u(IRHOU) + u(IRHOV)*u(IRHOV) ) / u(IRHO) + Dimensionless % cp) * u(IRHOV) / u(IRHO)

         F = F / (sqrt(gamma) * Mach)
         G = G / (sqrt(gamma) * Mach)

         end associate

      end function inviscidFlux0D

      function inviscidFlux1D(u) result(val)
         implicit none
         real(kind=RP)                      :: u(0:,:)
         real(kind=RP), allocatable, target :: val(:,:,:)
         real(kind=RP), pointer             :: F(:,:) , G(:,:)

         allocate( val(0:size(u,1)-1 , 1:NEC , 1:NDIM ) )

         F(0:,1:)    => val(0:,1:,iX)
         G(0:,1:)    => val(0:,1:,iY)

         associate ( Gamma => Thermodynamics % Gamma , gm1 => Thermodynamics % gm1 , Mach => Dimensionless % Mach ) 
    
         F(:,IRHO)  = u(:,IRHOU)
         F(:,IRHOU) = 0.5_RP * (3.0_RP - Gamma) * u(:,IRHOU)*u(:,IRHOU) / u(:,IRHO) - 0.5_RP * gm1 * u(:,IRHOV)*u(:,IRHOV)/u(:,IRHO) + gm1 * u(:,IRHOE)
         F(:,IRHOV) = u(:,IRHOU)*u(:,IRHOV) / u(:,IRHO)
         F(:,IRHOE) = (Gamma * u(:,IRHOE) - 0.5_RP*gm1*( u(:,IRHOU)*u(:,IRHOU) + u(:,IRHOV)*u(:,IRHOV) )/ u(:,IRHO) + Dimensionless % cp) * u(:,IRHOU) / u(:,IRHO)

         G(:,IRHO)  = u(:,IRHOV)
         G(:,IRHOU) = u(:,IRHOU)*u(:,IRHOV) / u(:,IRHO)
         G(:,IRHOV) = 0.5_RP * (3.0_RP - Gamma) * u(:,IRHOV)*u(:,IRHOV) / u(:,IRHO) - 0.5_RP * gm1 * u(:,IRHOU)*u(:,IRHOU)/u(:,IRHO) + gm1 * u(:,IRHOE)
         G(:,IRHOE) = (Gamma * u(:,IRHOE) - 0.5_RP*gm1*( u(:,IRHOU)*u(:,IRHOU) + u(:,IRHOV)*u(:,IRHOV) ) / u(:,IRHO) + Dimensionless % cp) * u(:,IRHOV) / u(:,IRHO)

         F = F / (sqrt(gamma) * Mach)
         G = G / (sqrt(gamma) * Mach)

         end associate


      end function inviscidFlux1D

      function inviscidFlux2D(u) result(val)
         implicit none
         real(kind=RP)                      :: u(0:,0:,:)
         real(kind=RP), allocatable, target :: val(:,:,:,:)
         real(kind=RP), pointer             :: F(:,:,:) , G(:,:,:)

         allocate( val(0:size(u,1)-1 , 0:size(u,2)-1 , 1:NEC , 1:NDIM ) )

         F(0:,0:,1:)    => val(0:,0:,1:,iX)
         G(0:,0:,1:)    => val(0:,0:,1:,iY)


         associate ( Gamma => Thermodynamics % Gamma , gm1 => Thermodynamics % gm1 , Mach => Dimensionless % Mach ) 
    
         F(:,:,IRHO)  = u(:,:,IRHOU)
         F(:,:,IRHOU) = 0.5_RP * (3.0_RP - Gamma) * u(:,:,IRHOU)*u(:,:,IRHOU) / u(:,:,IRHO) - 0.5_RP * gm1 * u(:,:,IRHOV)*u(:,:,IRHOV)/u(:,:,IRHO) + gm1 * u(:,:,IRHOE)
         F(:,:,IRHOV) = u(:,:,IRHOU)*u(:,:,IRHOV) / u(:,:,IRHO)
         F(:,:,IRHOE) = (Gamma * u(:,:,IRHOE) - 0.5_RP*gm1*( u(:,:,IRHOU)*u(:,:,IRHOU) + u(:,:,IRHOV)*u(:,:,IRHOV) )/ u(:,:,IRHO) + Dimensionless % cp)  & 
                              * u(:,:,IRHOU) / u(:,:,IRHO)

         G(:,:,IRHO)  = u(:,:,IRHOV)
         G(:,:,IRHOU) = u(:,:,IRHOU)*u(:,:,IRHOV) / u(:,:,IRHO)
         G(:,:,IRHOV) = 0.5_RP * (3.0_RP - Gamma) * u(:,:,IRHOV)*u(:,:,IRHOV) / u(:,:,IRHO) - 0.5_RP * gm1 * u(:,:,IRHOU)*u(:,:,IRHOU)/u(:,:,IRHO) + gm1 * u(:,:,IRHOE)
         G(:,:,IRHOE) = (Gamma * u(:,:,IRHOE) - 0.5_RP*gm1*( u(:,:,IRHOU)*u(:,:,IRHOU) + u(:,:,IRHOV)*u(:,:,IRHOV) ) / u(:,:,IRHO) + Dimensionless % cp) & 
                              * u(:,:,IRHOV) / u(:,:,IRHO)

         F = F / (sqrt(gamma) * Mach)
         G = G / (sqrt(gamma) * Mach)

         end associate

      end function inviscidFlux2D
!
!      function inviscidFlux1D(u) result(val)
!         implicit none
!         real(kind=RP)     :: u(:)
!         real(kind=RP),allocatable     :: val(:)
!
!         allocate ( val(size(u,1) ) )
!         val = 0.5_RP*u*u
!   
!      end function inviscidFlux1D
!
!      function inviscidFlux2D(u) result(val)
!         implicit none
!         real(kind=RP)     :: u(:,:)
!         real(kind=RP),allocatable     :: val(:,:)
!
!         allocate ( val(size(u,1) , size(u,2) ) )
!         val = 0.5_RP*u*u
!   
!      end function inviscidFlux2D
!
!      function viscousFlux0D(g,u) result(val)
!         implicit none
!         real(kind=RP), optional     :: u
!         real(kind=RP)               :: g
!         real(kind=RP)               :: val
!
!
!         val = Setup % nu * g
!
!      end function viscousFlux0D
!
!      function viscousFlux1D(g,u) result(val)
!         implicit none
!         real(kind=RP), optional     :: u(:)
!         real(kind=RP)               :: g(:)
!         real(kind=RP), allocatable     :: val(:)
!
!         allocate(val( size(g,1)) )
!
!         val = Setup % nu * g
!
!      end function viscousFlux1D
!
!      function viscousFlux2D(g,u) result(val)
!         implicit none
!         real(kind=RP), optional     :: u(:,:)
!         real(kind=RP)               :: g(:,:)
!         real(kind=RP), allocatable     :: val(:,:)
!
!         allocate(val( size(g,1) , size(g,2) ) )
!
!         val = Setup % nu * g
!
!      end function viscousFlux2D
!!
!!     ****************************************************
!!        Riemann solvers
!!     ****************************************************
!!
!      function RoeFlux(uL , uR , n) result(val)
!         implicit none
!         real(kind=RP), dimension(NEC)       :: uL
!         real(kind=RP), dimension(NEC)       :: uR
!         real(kind=RP), dimension(NEC)       :: val
!         real(kind=RP)                       :: n
!         real(kind=RP), dimension(NEC)       :: ustar
!         
!         ustar = (uL + uR)*n
!
!         if (ustar(1) .gt. 0.0_RP) then
!            val = 0.5_RP * uL*uL
!         elseif( ustar(1) .lt. 0.0_RP) then
!            val = 0.5_RP * uR*uR
!         else
!            val = 0.0_RP
!         end if
!
!      end function RoeFlux
!
!      function ECONFlux(uL , uR , n) result(val)
!         implicit none
!         real(kind=RP), dimension(NEC)    :: uL
!         real(kind=RP), dimensioN(NEC)    :: uR
!         real(kind=RP), dimension(NEC)    :: val
!         real(kind=RP)                    :: n
!         
!         val = 0.25_RP*( uL*uL + uR*uR ) - 1.0_RP / 12.0_RP * (uL - uR)*(uL-uR)
!
!      end function ECONFlux
!
!      function LocalLaxFriedrichsFlux(uL , uR , n) result(val)
!         implicit none
!         real(kind=RP), dimension(NEC)    :: uL
!         real(kind=RP), dimensioN(NEC)    :: uR
!         real(kind=RP), dimension(NEC)    :: val
!         real(kind=RP)                    :: n
!
!         val = 0.25_RP * (uL*uL + uR*uR) - 0.5_RP*max(abs(uL),abs(uR))*(uR-uL)*n
!
!      end function LocalLaxFriedrichsFlux

!
!      ***************************************************************************
!           Subroutine for the module description
!      ***************************************************************************
!
       subroutine Describe()
         use Headers
         implicit none
!
!        *************************
!        Show reference quantities
!        *************************
!
         write(STD_OUT,'(/,/)')
         call Section_header("Loading " // solver // " physics module")
         write(STD_OUT,'(/)')
         call Subsection_header("Fluid data")
         write(STD_OUT,'(30X,A,A25,A15,A)') "-> ","Gas: ",trim(thermodynamics % fluid)
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Gas constant: ",thermodynamics % R," S.I."
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Specific heat ratio: ",thermodynamics % gamma,"."
         write(STD_OUT,'(/)')
         call Subsection_header("Reference quantities")
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Reference Temperature: ",refValues % T," K."
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Reference pressure: ",refValues % p," Pa."
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Reference density: ",refValues % rho," kg/m^3."
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Reference velocity: ",refValues % V," m/s."
         write(STD_OUT,'(30X,A,A25,E15.3,A)') "-> ","Reynolds length: ",refValues % L," m."
         write(STD_OUT,'(30X,A,A25,E15.3,A)') "-> ","Reference viscosity: ",refValues % mu," Pa·s."
         write(STD_OUT,'(30X,A,A25,E15.3,A)') "-> ","Reference conductivity: ",refValues % kappa," W/(m·K)."
         write(STD_OUT,'(30X,A,A25,E15.3,A)') "-> ","Reference time: ",refValues % tc," s."
         write(STD_OUT,'(/)')
         call Subsection_header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A25,F15.3)') "-> ","Reynolds number: ",dimensionless % Re
         write(STD_OUT,'(30X,A,A25,F15.3)') "-> ","Prandtl number: ",dimensionless % Pr
         write(STD_OUT,'(30X,A,A25,F15.3)') "-> ","Mach number: " , dimensionless % Mach
 
       end subroutine Describe
end module PhysicsNS
