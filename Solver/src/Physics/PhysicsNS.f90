module PhysicsNS
    use SMConstants
    use Setup_class

#include "Defines.h"

    private
    public :: solver
    public :: RefValues , Dimensionless , Thermodynamics
    public :: RiemannSolverFunction , InitializePhysics
    public :: InviscidFlux , ViscousFlux , AdiabaticViscousFlux
    public :: ComputeViscousTensor
    public :: HLLFlux, RoeFlux , AUSMFlux , ExactRiemannSolver
    public :: getPressure , getSoundSpeed , getStrainTensor , getTemperatureGradient

!   ********************************************
!        Current solver
!   ********************************************
!
    character(LEN=*), PARAMETER     :: solver = "Navier-Stokes"
    integer, parameter              :: STR_LEN_PHYSICS = 128

    type Thermodynamics_t  
      character(len=STR_LEN_PHYSICS)       :: fluid
      real(kind=RP)        :: R
      real(kind=RP)        :: gamma       ! Specific heat ratio
      real(kind=RP)        :: gm1         ! (gamma - 1)
      real(kind=RP)        :: gogm1       ! gamma / ( gamma - 1 )
      real(kind=RP)        :: invgm1      ! 1.0_RP / (gamma - 1 )
      real(kind=RP)        :: cp
      real(kind=RP)        :: cv
      real(kind=RP)        :: lambda
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
      real(kind=RP)        :: sqrtGammaMach
      real(kind=RP)        :: gammaMach2
      real(kind=RP)        :: invsqrtGammaMach
    end type Dimensionless_t

    abstract interface
      pure function RiemannSolverFunction( QL , QR , n ) result ( val )
         use SMConstants
         real(kind=RP), dimension(NCONS), intent(in) :: QL
         real(kind=RP), dimension(NCONS), intent(in) :: QR
         real(kind=RP), dimension(NDIM) , intent(in) :: n
         real(kind=RP), dimension(NCONS)             :: val
      end function RiemannSolverFunction
    end interface



    type(Thermodynamics_t), target  :: thermodynamicsAir = Thermodynamics_t("Air",287.15_RP , 1.4_RP , &
                              0.4_RP , 3.5_RP , 2.5_RP , 287.0_RP*3.5_RP , 287.0_RP*2.5_RP , -2.0_RP / 3.0_RP )
    type(Thermodynamics_t), pointer, protected            :: thermodynamics
    type(RefValues_t), protected       :: refValues      
    type(Dimensionless_t), protected   :: dimensionless

!
    interface inviscidFlux
      module procedure inviscidFlux0D , inviscidFlux1D , inviscidFlux2D
    end interface inviscidFlux

    interface viscousFlux
      module procedure viscousFlux0D , viscousFlux1D , viscousFlux2D , viscousFluxBC0D , viscousFluxBC1D
    end interface viscousFlux

    interface adiabaticViscousFlux
      module procedure adiabaticViscousFlux0D , adiabaticViscousFlux1D , adiabaticViscousFlux2D
    end interface adiabaticViscousFlux

    interface getPressure
      module procedure getPressure0D , getPressure1D , getPressure2D
    end interface getPressure

    interface getTemperature
      module procedure getTemperature0D , getTemperature1D , getTemperature2D
    end interface getTemperature

    interface getSoundSpeed
      module procedure getSoundSpeed0D , getSoundSpeed1D , getSoundSpeed2D 
    end interface getSoundSpeed
#ifdef NAVIER_STOKES
    interface getStrainTensor
      module procedure getStrainTensor0D , getStrainTensor1D , getStrainTensor2D   
    end interface getStrainTensor

    interface getTemperatureGradient
      module procedure getTemperatureGradient0D , getTemperatureGradient1D , getTemperatureGradient2D
    end interface getTemperatureGradient
#endif
!
!   //////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           Define module procedures
!
!   //////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
    interface
      module pure function getPressure0D(Q) result (p)
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)                       :: p
      end function getPressure0D

      module pure function getPressure1D(N,Q) result (p)
         implicit none
         integer,       intent (in)  :: N
         real(kind=RP), intent (in)  :: Q(0:N,1:NCONS)
         real(kind=RP)               :: p(0:N)
      end function getPressure1D

      module pure function getPressure2D(N,Q) result (p)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)               :: p(0:N,0:N)
      end function getPressure2D

      module pure function getTemperature0D(Q) result (T)
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)                       :: T
      end function getTemperature0D

      module pure function getTemperature1D(N,Q) result (T)
         implicit none
         integer,       intent (in)  :: N
         real(kind=RP), intent (in)  :: Q(0:N,1:NCONS)
         real(kind=RP)               :: T(0:N)
      end function getTemperature1D

      module pure function getTemperature2D(N,Q) result (T)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)               :: T(0:N,0:N)
      end function getTemperature2D

      module pure function getSoundSpeed0D(Q) result (a)
         implicit none
         real(kind=RP), intent(in)        :: Q(1:NCONS)
         real(kind=RP)                    :: a
      end function getSoundSpeed0D

      module pure function getSoundSpeed1D(N,Q) result (a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(0:N,1:NCONS)
         real(kind=RP)                    :: a(0:N)
      end function getSoundSpeed1D

      module pure function getSoundSpeed2D(N,Q) result(a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)                    :: a(0:N,0:N)
      end function getSoundSpeed2D
#ifdef NAVIER_STOKES
      module pure function getStrainTensor0D ( Q , dQ ) result ( du )
         implicit none
         real(kind=RP), intent(in)     :: Q(1:NCONS)
         real(kind=RP), intent(in)     :: dQ(1:NDIM,1:NCONS)
         real(kind=RP)                 :: du(1:NDIM , 1:NDIM)     
      end function getStrainTensor0D
 
      module pure function getStrainTensor1D ( N ,  Q , dQ ) result ( du )
         implicit none
         integer,       intent(in)     :: N
         real(kind=RP), intent(in)     :: Q(0:N,1:NCONS)
         real(kind=RP), intent(in)     :: dQ(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                 :: du(0:N,1:NDIM , 1:NDIM)     
      end function getStrainTensor1D
 
      module pure function getStrainTensor2D ( N ,  Q , dQ ) result ( du )
         implicit none
         integer,       intent(in)     :: N
         real(kind=RP), intent(in)     :: Q(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)     :: dQ(0:N,0:N,1:NDIM,1:NCONS)
         real(kind=RP)                 :: du(0:N,0:N,1:NDIM , 1:NDIM)     
      end function getStrainTensor2D

      module pure function getTemperatureGradient0D( Q , dQ ) result ( gradT )
         implicit none
         real(kind=RP),    intent(in)        ::  Q(1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(1:NDIM , 1:NCONS)
         real(kind=RP)                       :: gradT(1:NDIM) 
      end function getTemperatureGradient0D

      module pure function getTemperatureGradient1D( N , Q , dQ ) result ( gradT )
         implicit none
         integer      ,    intent(in)        :: N 
         real(kind=RP),    intent(in)        ::  Q(0:N,1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(0:N,1:NDIM , 1:NCONS)
         real(kind=RP)                       :: gradT(0:N,1:NDIM) 
      end function getTemperatureGradient1D

      module pure function getTemperatureGradient2D( N , Q , dQ ) result ( gradT )
         implicit none
         integer      ,    intent(in)        :: N 
         real(kind=RP),    intent(in)        ::  Q(0:N,0:N,1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(0:N,0:N,1:NDIM , 1:NCONS)
         real(kind=RP)                       :: gradT(0:N,0:N,1:NDIM) 
      end function getTemperatureGradient2D
#endif 
    end interface

    interface
      module pure function ExactRiemannSolver(qL , qR , n) result (Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS), intent(in) :: qL
         real(kind=RP), dimension(NCONS), intent(in) :: qR
         real(kind=RP), dimension(NDIM) , intent(in) :: n
         real(kind=RP), dimension(NCONS) :: Fstar
      end function ExactRiemannSolver
         
      module pure function RoeFlux(qL, qR , n) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS), intent(in)     :: qL
         real(kind=RP), dimension(NCONS), intent(in)     :: qR
         real(kind=RP), dimension(NDIM) , intent(in)     :: n
         real(kind=RP), dimension(NCONS)     :: Fstar
      end function RoeFlux

      module pure function HLLFlux(qL , qR , n) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS), intent(in)     :: qL
         real(kind=RP), dimension(NCONS), intent(in)     :: qR
         real(kind=RP), dimension(NDIM) , intent(in)     :: n
         real(kind=RP), dimension(NCONS)     :: Fstar
      end function HLLFlux
    end interface

    interface
      module pure function inviscidFlux0D(u) result(val)
         implicit none
         real(kind=RP), intent(in) :: u(NCONS)
         real(kind=RP)             :: val(NCONS,NDIM)
      end function inviscidFlux0D

      module pure function inviscidFlux1D(N,u) result(val)
         implicit none
         integer, intent(in)       :: N
         real(kind=RP), intent(in) :: u(0:N,1:NCONS)
         real(kind=RP)             :: val(0:N,1:NCONS,1:NDIM)
      end function inviscidFlux1D

      module pure function InviscidFlux2D(N,u) result(val)
         implicit none
         integer, intent(in)      :: N
         real(kind=RP),intent(in) :: u(0:N,0:N,1:NCONS)
         real(kind=RP)            :: val(0:N,0:N,1:NCONS,1:NDIM)
      end function inviscidFlux2D
     
      module pure function F_inviscidFlux(rho,u,v,p,H) result(F)
         implicit none
         real(kind=RP), intent(in) :: rho
         real(kind=RP), intent(in) :: u
         real(kind=RP), intent(in) :: v
         real(kind=RP), intent(in) :: p
         real(kind=RP), intent(in) :: H
         real(kind=RP)             :: F(NCONS)
      end function F_inviscidFlux
    end interface

    interface
      module pure function viscousFlux0D( q , dq) result(F)
         implicit none
         real(kind=RP), intent(in)    :: q(NCONS)
         real(kind=RP), intent(in)    :: dq(NDIM , NCONS)
         real(kind=RP)                :: F(1:NCONS,NDIM)
      end function viscousFlux0D

      module pure function viscousFlux1D( N , q , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,1:NCONS,1:NDIM)
      end function viscousFlux1D

      module pure function viscousFlux2D( N , q , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: Q(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dQ(0:N,0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,0:N,1:NCONS,1:NDIM)
      end function viscousFlux2D

      module pure function viscousFluxBC0D( q , qB , dq ) result ( F )
         implicit none
         real(kind=RP), intent(in)          :: q(1:NCONS)
         real(kind=RP), intent(in)          :: qB(1:NCONS)
         real(kind=RP), intent(in)          :: dq(1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(1:NCONS,1:NDIM)
      end function viscousFluxBC0D

      module pure function viscousFluxBC1D( N , q , qB , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: qB(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,1:NCONS,1:NDIM)
      end function viscousFluxBC1D

      module pure function adiabaticViscousFlux0D( q , qb , dq) result(F)
         implicit none
         real(kind=RP), intent(in)    :: q(NCONS)
         real(kind=RP), intent(in)    :: qB(NCONS)
         real(kind=RP), intent(in)    :: dq(NDIM , NCONS)
         real(kind=RP)                :: F(1:NCONS,NDIM)
      end function adiabaticViscousFlux0D

      module pure function adiabaticViscousFlux1D( N , q , qb , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: qb(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,1:NCONS,1:NDIM)
      end function adiabaticViscousFlux1D

      module pure function adiabaticViscousFlux2D( N , q , qb , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: Q(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)          :: Qb(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dQ(0:N,0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,0:N,1:NCONS,1:NDIM)
      end function adiabaticViscousFlux2D

      module pure function ComputeViscousTensor ( N , Q , dQ ) result ( tau )
         implicit none
         integer,          intent(in)     :: N
         real(kind=RP),    intent(in)     ::  Q(0:N,1:NCONS)
         real(kind=RP),    intent(in)     :: dQ(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                    :: tau(0:N,1:NDIM,1:NDIM)
      end function ComputeViscousTensor
   end interface

!   ========
    contains
!   ========
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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
#ifdef _DIMENSIONLESS_TAU
         refValues % L     = Setup % reynolds_length
         refValues % T     = Setup % temperature_ref
         refValues % p     = Setup % pressure_ref
         refValues % rho   = refValues % p / ( thermodynamics % R * refValues % T )
         refValues % a     = sqrt( refValues % p / refValues % rho )
         refValues % Mach  = Setup % Mach_number
         refValues % V     = sqrt( thermodynamics % gamma ) * refValues % Mach * refValues % a
         if ( setup % reynolds_number .eq. 0.0_RP ) then
            refValues % mu = 0.0_RP
         else
            refValues % mu    = refValues % rho * refValues % V * refValues % L / Setup % reynolds_number
         end if
         refValues % kappa = refValues % mu * thermodynamics % cp / Setup % prandtl_number
         refValues % tc    = refValues % L / (refValues % V )
#else
         refValues % L     = Setup % reynolds_length
         refValues % T     = Setup % temperature_ref
         refValues % rho   = Setup % pressure_ref / (thermodynamics % R * refValues % T)
         refValues % Mach  = Setup % Mach_number
         refValues % V     = refValues % Mach * sqrt( thermodynamics % gamma * Setup % pressure_ref / refValues % rho )
         refValues % a     = refValues % V
         refValues % p     = refValues % Mach * refValues % Mach * thermodynamics % gamma * Setup % pressure_ref
         refValues % mu    = refValues % rho * refValues % V * refValues % L / Setup % reynolds_number
         refValues % kappa = refValues % mu * thermodynamics % cp / Setup % prandtl_number
         refValues % tc    = refValues % L / (refValues % V )
#endif
!
!        ***********************************
!        Initialize the dimensionless values
!        ***********************************
!
         dimensionless % Mach             = refValues % Mach
         dimensionless % sqrtGammaMach    = sqrt(thermodynamics % gamma) * dimensionless % Mach
         dimensionless % gammaMach2       = thermodynamics % gamma * dimensionless % Mach * dimensionless % Mach
         dimensionless % invSqrtGammaMach = 1.0_RP / (sqrt(thermodynamics % gamma) * dimensionless % Mach)

         if ( setup % reynolds_number .eq. 0.0_RP ) then
            dimensionless % mu = 0.0_RP
            dimensionless % kappa = 0.0_RP

         else
            dimensionless % mu               = 1.0_RP / Setup % reynolds_number
            dimensionless % kappa            = thermodynamics % gogm1 / (Setup % prandtl_number * Setup % reynolds_number)

         end if

         dimensionless % cp               = thermodynamics % gogm1
         dimensionless % cv               = thermodynamics % invgm1
         dimensionless % Re               = Setup % reynolds_number
         dimensionless % Pr               = Setup % prandtl_number

         call Describe

       end subroutine InitializePhysics

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
         write(STD_OUT,'(30X,A,A25,F15.3,A)') "-> ","Reference sound speed: ",refValues % a," m/s."
         write(STD_OUT,'(30X,A,A25,ES15.3,A)') "-> ","Reynolds length: ",refValues % L," m."
         write(STD_OUT,'(30X,A,A25,ES15.3,A)') "-> ","Reference viscosity: ",refValues % mu," Pa·s."
         write(STD_OUT,'(30X,A,A25,ES15.3,A)') "-> ","Reference conductivity: ",refValues % kappa," W/(m·K)."
         write(STD_OUT,'(30X,A,A25,ES15.3,A)') "-> ","Reference time: ",refValues % tc," s."
         write(STD_OUT,'(/)')
         call Subsection_header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A25,F15.3)') "-> ","Reynolds number: ",dimensionless % Re
         write(STD_OUT,'(30X,A,A25,F15.3)') "-> ","Prandtl number: ",dimensionless % Pr
         write(STD_OUT,'(30X,A,A25,F15.3)') "-> ","Mach number: " , dimensionless % Mach
 
       end subroutine Describe
end module PhysicsNS
