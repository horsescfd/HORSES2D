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
!
!    interface inviscidFlux
!      module procedure inviscidFlux0D , inviscidFlux1D , inviscidFlux2D
!    end interface inviscidFlux
!
!    interface viscousFlux
!      module procedure viscousFlux0D , viscousFlux1D , viscousFlux2D
!    end interface viscousFlux
!
!    contains
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
