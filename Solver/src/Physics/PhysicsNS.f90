module PhysicsNS
    use SMConstants
    use Setup_class

    private
    public :: NEC , NDIM , IX , IY , IRHO , IRHOU , IRHOV , IRHOE , solver
    public :: RefValues , Dimensionless , Thermodynamics
    public :: RiemannSolverFunction , InitializePhysics
    public :: InviscidFlux
    public :: RoeFlux
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
      function RiemannSolverFunction( QL , QR , T , Tinv ) result ( val )
         use SMConstants
         import NEC , NDIM
         real(kind=RP), dimension(NEC)     :: QL
         real(kind=RP), dimension(NEC)     :: QR
         real(kind=RP), dimension(NEC,NEC) :: T
         real(kind=RP), dimension(NEC,NEC) :: Tinv
         real(kind=RP), dimension(NEC)     :: val
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
     
      function F_inviscidFlux(u) result(F)
         implicit none
         real(kind=RP)        :: u(NEC)
         real(kind=RP)        :: F(NEC)
   
         associate ( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )    

         F(IRHO)  = u(IRHOU)
         F(IRHOU) = 0.5_RP * (3.0_RP - Gamma) * u(IRHOU)*u(IRHOU) / u(IRHO) - 0.5_RP * gm1 * u(IRHOV)*u(IRHOV)/u(IRHO) + gm1 * u(IRHOE)
         F(IRHOV) = u(IRHOU)*u(IRHOV) / u(IRHO)
         F(IRHOE) = (Gamma * u(IRHOE) - 0.5_RP*gm1*( u(IRHOU)*u(IRHOU) + u(IRHOV)*u(IRHOV) ) / u(IRHO) + Dimensionless % cp) * u(IRHOU) / u(IRHO)

         end associate
      end function F_inviscidFlux

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

!
!     ****************************************************
!        Riemann solvers
!     ****************************************************
!
      function RoeFlux(qL3D , qR3D , T , Tinv) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NEC)     :: qL3D
         real(kind=RP), dimension(NEC)     :: qR3D
         real(kind=RP), dimension(NEC,NEC) :: T
         real(kind=RP), dimension(NEC,NEC) :: Tinv
         real(kind=RP), dimension(NEC)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP), dimension(NEC) :: qL , qR
         real(kind=RP)                 :: rhoL , uL , vL , HL
         real(kind=RP)                 :: rhoR , uR , vR , HR
         real(kind=RP)                 :: rho , u , v , H , a
         real(kind=RP)                 :: dq(NEC)
         real(kind=RP)                 :: lambda(NEC)
         real(kind=RP)                 :: K(NEC,NEC)
         real(kind=RP)                 :: alpha(NEC)
         integer                       :: eq

!        0/ Gather variables
!           ----------------

            qL = MatrixTimesVector_F( A=T , X=qL3D )
            qR = MatrixTimesVector_F( A=T , X=qR3D )

            associate( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )
            rhoL = sqrt(qL(IRHO))
            uL   = qL(IRHOU) / qL(IRHO)
            vL   = qL(IRHOV) / qL(IRHO)
            HL   = gamma * qL(IRHOE) / qL(IRHO) - 0.5_RP * gm1 * ( uL*uL + vL*vL )

            rhoR = sqrt(qR(IRHO))
            uR   = qR(IRHOU) / qR(IRHO)
            vR   = qR(IRHOV) / qR(IRHO)
            HR   = gamma * qR(IRHOE) / qR(IRHO) - 0.5_RP * gm1 * ( uR*uR + vR*vR )
            end associate

!        1/ Compute Roe averages
!           --------------------
            rho = rhoL + rhoR
            u   = (rhoL*uL + rhoR*uR)/rho
            v   = (rhoL*vL + rhoR*vR)/rho
            H   = (rhoL*HL + rhoR*HR)/rho
            associate( gm1 => Thermodynamics % gm1 ) 
            a   = sqrt(gm1*(H - 0.5_RP*(u*u + v*v) ) )
            end associate

!
!        2/ Compute Roe matrix eigenvalues
!           ------------------------------
            lambda(1) = u - a
            lambda(2) = u
            lambda(3) = u
            lambda(4) = u + a 

!
!        3/ Compute the averaged right eigenvectors
!           ---------------------------------------
            K(1:NEC , 1)  = reshape( (/ 1.0_RP , u-a    , v      , H-u*a                /)  ,  (/ NEC /) )
            K(1:NEC , 2)  = reshape( (/ 1.0_RP , u      , v      , 0.5_RP * (u*u + v*v) /)  ,  (/ NEC /) )
            K(1:NEC , 3)  = reshape( (/ 0.0_RP , 0.0_RP , 1.0_RP , v                    /)  ,  (/ NEC /) )
            K(1:NEC , 4)  = reshape( (/ 1.0_RP , u + a  , v      , H + u*a              /)  ,  (/ NEC /) )
!
!        4/ Compute the wave strengths
!           --------------------------
            associate( gm1 => Thermodynamics % gm1 ) 
            dq = qR - qL
            alpha(3) = dq(IRHOV) - v*dq(IRHO)
            alpha(2) = (dq(IRHO) * (H - u*u) + u*dq(IRHOU) - dq(IRHOE) + (dq(IRHOV) - v*dq(IRHO))*v) / ( a*a )
            alpha(1) = (dq(IRHO) * (u + a) - dq(IRHOU) - a * alpha(2)) / (2.0_RP*a)
            alpha(4) = dq(IRHO) - (alpha(1) + alpha(2)) 
            end associate
!
!        5/ Compute the flux
!           ----------------
            Fstar = 0.5_RP * ( F_inviscidFlux(qL) + F_inviscidFLux(qR) )
               
            do eq = 1 , NEC
               Fstar = Fstar - 0.5_RP * alpha(eq) * abs(lambda(eq)) * K(1:NEC , eq)
            end do

!
!        6/ Scale it with the Mach number
!           -----------------------------
            associate( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
            Fstar = Fstar / ( sqrt(gamma) * Mach)
            end associate
            
!  
!        7/ Return to the 3D Space
!           ----------------------
            Fstar = MatrixTimesVector_F( A = Tinv , X = Fstar )

      end function RoeFlux
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
