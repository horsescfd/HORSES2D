module PhysicsNS
    use SMConstants
    use Setup_class

    private
    public :: NEC , NDIM , IX , IY , IRHO , IRHOU , IRHOV , IRHOE , solver
!
!   *****************************************
!        Definitions
!   *****************************************
!
    integer, parameter              :: NEC = 4
    integer, parameter              :: NDIM = 2
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
      real(kind=RP)        :: R
      real(kind=RP)        :: gamma       ! Specific heat ratio
      real(kind=RP)        :: gm1         ! (gamma - 1)
      real(kind=RP)        :: gogm1       ! gamma / ( gamma - 1 )
      real(kind=RP)        :: invgm1      ! 1.0_RP / (gamma - 1 )
    end type Thermodynamics_t

    type RefValues_t
      real(kind=RP)        :: L
      real(kind=RP)        :: T
      real(kind=RP)        :: p
      real(kind=RP)        :: rho
      real(kind=RP)        :: V
      real(kind=RP)        :: Mach
      real(kind=RP)        :: mu
      real(kind=RP)        :: kappa
      real(kind=RP)        :: t
    end type RefValues_t

    type Dimensionless_t
      real(kind=RP)        :: mu
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

    type(Thermodynamics_t), parameter  :: thermodynamics = Thermodynamics_t(287.0_RP , 1.4_RP , 0.4_RP , 3.5_RP , 2.5_RP)
    type(RefValues_t), protected       :: refValues      

!
!    interface inviscidFlux
!      module procedure inviscidFlux0D , inviscidFlux1D , inviscidFlux2D
!    end interface inviscidFlux
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
!        Initialize the reference values
!        *******************************
!
         refValues % p = Setup % pressure_ref
         refValues % T = Setup % temperature_ref
         refValues % rho = refValues % p / (thermodynamics % R * refValues % T)
         refValues % L = Setup % reynolds_length
         refValues %  = Setup % reynolds_number

       end subroutine InitializePhysics
!
!      function inviscidFlux0D(u) result(val)
!         implicit none
!         real(kind=RP)     :: u
!         real(kind=RP)    :: val
!
!         val = 0.5_RP * u * u
!
!      end function inviscidFlux0D
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
end module PhysicsNS
