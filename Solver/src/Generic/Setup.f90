module Setup_class
    use SMConstants
    use ParamfileIO
    implicit none

    private
    public  :: setup
  
    integer, parameter          :: STR_LEN_SETUP = 128
    type Setup_t
!
!       -------------------------------------------------------------------------------------
!              Reference quantities
!       -------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: mesh_file           
        character(len=STR_LEN_SETUP) :: bdry_file            
!
!       -------------------------------------------------------------------------------------
!              Reference quantities
!       -------------------------------------------------------------------------------------
!
        real(kind=RP), allocatable                :: pressure_ref
        real(kind=RP), allocatable                :: temperature_ref           
        real(kind=RP), allocatable                :: reynolds_length          
        real(kind=RP), allocatable                :: reynolds_number          
        real(kind=RP), allocatable                :: prandtl_number           
        real(kind=RP), allocatable                :: Mach_number              
        character(len=STR_LEN_SETUP) :: Gas                    
!
!       -------------------------------------------------------------------------------------
!              Spatial discretization parameters
!       -------------------------------------------------------------------------------------
!
        integer, allocatable                      :: nodes                          ! Interpolation / Integration nodes strategy
        integer, allocatable                      :: N                              ! Polynomial order (generic)
!
!       --------------------------------------------------------------------------------------
!              Initialization
!       --------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: IC                                          ! Initial condition type
!
!       -----------------------------------------------------------------------------------------
!              Advective flux discretization
!       -----------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: inviscid_discretization   
        integer                      :: inviscid_formulation      
        character(len=STR_LEN_SETUP) :: inviscid_flux             
        integer, allocatable         :: integration_points        
!     
!       -------------------------------------------------------------------------------------------
!              Viscous discretization
!       -------------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: viscous_discretization    
        character(len=STR_LEN_SETUP) :: IPMethod                  
        real(kind=RP), allocatable   :: sigma0IP                  
        real(kind=RP), allocatable   :: sigma1IP                  
!
!       ------------------------------------------------------------------------------
!              Integration parameters
!       ------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: integrationMode           
        real(kind=RP), allocatable   :: dt                        
        real(kind=RP), allocatable   :: simulationTime            
        integer, allocatable         :: no_of_iterations          
        real(kind=RP), allocatable   :: initialTime               
        character(len=STR_LEN_SETUP) :: integrationMethod         
!
!       ------------------------------------------------------------------------------
!             Output parameters
!       ------------------------------------------------------------------------------
!
        integer, allocatable         :: autosaveInterval         
        integer, allocatable         :: output_interval          
        character(len=STR_LEN_SETUP) :: saveVariables            

        contains
            procedure :: Initialization => Setup_Initialization
            procedure :: Default        => Setup_DefaultValues
    end type Setup_t

    type(Setup_t), protected, target       :: setup

    contains
    
      subroutine Setup_Initialization( self )
         implicit none
         class(Setup_t)               :: self
         integer                      :: nArgs
         character(len=STR_LEN_SETUP) :: arg
         integer                      :: iArg
         integer                      :: pos
         character(len=STR_LEN_SETUP) :: case_name
         character(len=STR_LEN_SETUP) :: interp_nodes
         character(len=STR_LEN_SETUP) :: inviscid_form

!
!         Get case file from command line
!         -------------------------------
          nArgs = command_argument_count()
          do iArg = 1 , nArgs
            call get_command_argument(iArg , arg)
     
            pos = index(trim(arg) , '.HiOCase')
          
            if ( pos .gt. 0 ) then
               case_name = trim(arg)
               exit
            elseif ( iArg .eq. nArgs ) then
               write(STD_OUT,'(/)')
               print*, "No input file specified."
               stop "Stopped."
            end if
           
          end do


!
!         Read from case file
!         -------------------
          call readValue(trim(case_name) , "Mesh file" , Setup % mesh_file )
          call readValue(trim(case_name) , "Boundary file" , Setup % bdry_file )

          if ( trim(Setup % bdry_file) .eq. "_this_" ) Setup % bdry_file = case_name

          call readValue ( trim ( case_name )  , "Gas"                   , Setup % Gas             ) 
          call readValue ( trim ( case_name )  , "Reference pressure"    , Setup % pressure_ref    ) 
          call readValue ( trim ( case_name )  , "Reference Temperature" , Setup % Temperature_ref ) 
          call readValue ( trim ( case_name )  , "Reynolds length"       , Setup % reynolds_length ) 
          call readValue ( trim ( case_name )  , "Reynolds number"       , Setup % reynolds_number ) 
          call readValue ( trim ( case_name )  , "Prandtl number"        , Setup % prandtl_number  ) 
          call readValue ( trim ( case_name )  , "Mach number"           , Setup % Mach_number     ) 

          call readValue ( trim ( case_name )  , "Interpolation nodes"   , interp_nodes ) 
          
          if ( trim(interp_nodes) .eq. "Legendre-Gauss" ) then
            Setup % nodes = LG
          elseif ( trim(interp_nodes) .eq. "Legendre-Gauss-Lobatto" ) then
            Setup % nodes = LGL
          end if

          call readValue ( trim ( case_name )  , "Polynomial order" , Setup % N ) 
          call readValue ( trim ( case_name )  , "Initial condition" , Setup % IC ) 
          call readValue ( trim ( case_name )  , "Inviscid strategy" , Setup % inviscid_discretization ) 
          call readValue ( trim ( case_name )  , "Inviscid formulation" , inviscid_form ) 
      
          if ( trim(inviscid_form) .eq. "Form I" ) then
            Setup % inviscid_formulation = 1
          elseif ( trim(inviscid_form) .eq. "Form II" ) then
            Setup % inviscid_formulation = 2
          end if

          call readValue ( trim ( case_name )  , "Inviscid Riemann Flux"        , Setup % inviscid_flux          ) 
          call readValue ( trim ( case_name )  , "Number of integration points" , Setup % integration_points     ) 
          call readValue ( trim ( case_name )  , "Viscous strategy"             , Setup % viscous_discretization ) 
         
          if ( Setup % viscous_discretization .eq. "Interior-penalty" ) then
            Setup % viscous_discretization = "IP"
          end if

          call readValue ( trim ( case_name )  , "Interior penalty method"          , Setup % IPMethod          ) 
          call readValue ( trim ( case_name )  , "Jumps penalty parameter"          , Setup % sigma0IP          ) 
          call readValue ( trim ( case_name )  , "Gradient jumps penalty parameter" , Setup % sigma1IP          ) 
          call readValue ( trim ( case_name )  , "Integration mode"                 , Setup % integrationMode   ) 
          call readValue ( trim ( case_name )  , "Integration scheme"               , Setup % integrationMethod ) 
          call readValue ( trim ( case_name )  , "Time step"                        , Setup % dt                ) 
          call readValue ( trim ( case_name )  , "Simulation time"                  , Setup % simulationTime    ) 
          call readValue ( trim ( case_name )  , "Number of iterations"             , Setup % no_of_iterations  ) 
          call readValue ( trim ( case_name )  , "Initial time"                     , Setup % initialTime       ) 
          call readValue ( trim ( case_name )  , "Autosave interval"                , Setup % AutosaveInterval  ) 
          call readValue ( trim ( case_name )  , "Output interval"                  , Setup % Output_Interval   ) 
          call readValue ( trim ( case_name )  , "Save variables"                   , Setup % saveVariables     ) 
         
          call Setup % Default

      end subroutine Setup_Initialization

      subroutine Setup_DefaultValues( self )
         implicit none
         class(Setup_t)             :: self


          if ( .not. allocated ( Setup % integration_points ) ) then
            allocate( Setup % integration_points ) 
            Setup % integration_points = Setup % N 
          end if

      end subroutine Setup_DefaultValues


end module Setup_class
