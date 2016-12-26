module Setup_class
    use SMConstants
    implicit none

    private
    public  :: setup
  
    integer, parameter          :: STR_LEN_SETUP = 128
    type Setup_t
!
!       -------------------------------------------------------------------------------------
!              Spatial discretization parameters
!       -------------------------------------------------------------------------------------
!
        integer                      :: nodes                     =  LG       ! Interpolation / Integration nodes strategy
        integer                      :: N                         =  5         ! Polynomial order (generic)
        real(kind=RP)                :: nu                        =  0.1_RP ! Viscous coefficient
!
!       --------------------------------------------------------------------------------------
!              Initialization
!       --------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: IC                        =  "Wave"    ! Initial condition type
!
!       -----------------------------------------------------------------------------------------
!              Advective flux discretization
!       -----------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: inviscid_discretization   = "Standard"
        character(len=STR_LEN_SETUP) :: inviscid_flux             = "Roe"
        integer                      :: integration_points        =  ceiling((3.0_RP *5.0_RP + 1.0_RP)/2.0_RP)
!     
!       -------------------------------------------------------------------------------------------
!              Viscous discretization
!       -------------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: viscous_discretization    =  "IP"
        character(len=STR_LEN_SETUP) :: IPMethod                  =  "SIPG"
        real(kind=RP)                :: sigma0IP                  =  10.0_RP
        real(kind=RP)                :: sigma1IP                  =  00.0_RP
!
!       ---------------------------------------------------------------------------------------
!              Boundary conditions parameters
!       ---------------------------------------------------------------------------------------
!
        integer, dimension(4)        :: markers                   =  [1,2,3,4]
        integer, dimension(2)        :: BCTypes                   =  [PERIODIC_BC , PERIODIC_BC]
        integer, dimension(4)        :: periodicBCFaces           =  [2,1,3,4]
        real(kind=RP), dimension(2)  :: dirichletBC               =  [0.0_RP , 0.0_RP]
!
!       ------------------------------------------------------------------------------
!              Integration parameters
!       ------------------------------------------------------------------------------
!
        integer                      :: integrationMode           =  STEADY
        real(kind=RP)                :: dt                        =  1.0e-5_RP
        real(kind=RP)                :: simulationTime            =  1.0_RP
        integer                      :: no_of_iterations          =  100000
        real(kind=RP)                :: initialTime               = 0.0_RP
        character(len=STR_LEN_SETUP) :: integrationMethod         = "Explicit-Euler"
!
!       ------------------------------------------------------------------------------
!             Output parameters
!       ------------------------------------------------------------------------------
!
        integer                      :: autosaveInterval          = 100
        character(len=STR_LEN_SETUP) :: saveVariables             = "Q_QDot_dQ"
    end type Setup_t

    type(Setup_t), target       :: setup


end module Setup_class
