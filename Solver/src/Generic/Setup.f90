!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
module Setup_class
    use SMConstants
    use ParamfileIO
    implicit none

#include "Defines.h"

    private
    public  :: setup , Setup_t
  
    integer, parameter          :: STR_LEN_SETUP = 128
    type Setup_t
!
!       -------------------------------------------------------------------------------------
!              Reference quantities
!       -------------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: case_file
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
!              Artificial dissipation
!       ------------------------------------------------------------------------------
!
        logical         :: artificialDissipation
        real(kind=RP), allocatable   :: artificialDissipationIntensity
        character(len=STR_LEN_SETUP)   :: artificialDissipationIndicator
        character(len=STR_LEN_SETUP)   :: artificialDissipationType
!
!       ------------------------------------------------------------------------------
!              Integration parameters
!       ------------------------------------------------------------------------------
!
        character(len=STR_LEN_SETUP) :: integrationMode           
        real(kind=RP), allocatable   :: Ccfl
        real(kind=RP), allocatable   :: dt
        real(kind=RP), allocatable   :: simulationTime            
        real(kind=RP), allocatable   :: residualTarget 
        integer, allocatable         :: no_of_iterations          
        real(kind=RP)                :: initialTime = 0.0_RP          
        integer                      :: initialIteration = 0
        character(len=STR_LEN_SETUP) :: integrationMethod         
!
!       ------------------------------------------------------------------------------
!             Output parameters
!       ------------------------------------------------------------------------------
!
        integer, allocatable         :: autosaveInterval         
        integer, allocatable         :: output_interval          
        integer, allocatable         :: no_of_plotPoints
        character(len=STR_LEN_SETUP) :: saveVariables            
        character(len=STR_LEN_SETUP) :: solution_file
        character(len=STR_LEN_SETUP) :: restart_file
        character(len=STR_LEN_SETUP) :: outputType
        character(len=STR_LEN_SETUP) :: exportFormat

        contains
            procedure, nopass  :: Initialization => Setup_Initialization
            procedure :: SetInitialTime => Setup_SetInitialTime
    end type Setup_t

    type(Setup_t), protected, target       :: setup

    interface Setup_CheckWithError
      module procedure Setup_CheckWithError_Character , Setup_CheckWithError_Integer , Setup_CheckWithError_Real
    end interface Setup_CheckWithError

    interface Setup_CheckWithDefault
      module procedure Setup_CheckWithDefault_Character , Setup_CheckWithDefault_Integer , Setup_CheckWithDefault_Real
    end interface Setup_CheckWithDefault

    contains
    
      subroutine Setup_Initialization
         implicit none
         integer                      :: nArgs
         character(len=STR_LEN_SETUP) :: arg
         integer                      :: iArg
         integer                      :: pos
         character(len=STR_LEN_SETUP) :: case_name
         character(len=STR_LEN_SETUP) :: interp_nodes
         character(len=STR_LEN_SETUP) :: inviscid_form
         integer, allocatable         :: artificialDissipation

!
!         Get case file from command line
!         -------------------------------
          nArgs = command_argument_count()
   
          if (nArgs .eq. 0) then
            print*,""
            print*,""
            print*, "No case file(s) selected"
            stop "Stopped"
         end if
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
          Setup % case_file = trim(case_name)
!
!         Request mesh file
!         -----------------
          call readValue(trim(case_name) , "Mesh file" , Setup % mesh_file )
          call Setup_CheckWithError( Setup % mesh_file , "Mesh file" )
!
!         Request boundary file
!         ---------------------
          call readValue(trim(case_name) , "Boundary file" , Setup % bdry_file )
          call Setup_CheckWithError( Setup % mesh_file , "Boundary file" )
          if ( trim(Setup % bdry_file) .eq. "_this_" ) Setup % bdry_file = case_name
!
!         Request gas
!         -----------
          call readValue ( trim ( case_name )  , "Gas"                   , Setup % Gas             ) 
          call Setup_CheckWithDefault( Setup % Gas , "Air" , "Gas" )
!
!         Request reference pressure          
!         --------------------------
          call readValue ( trim ( case_name )  , "Reference pressure"    , Setup % pressure_ref    ) 
          call Setup_CheckWithDefault( Setup % pressure_ref , 101325.0_RP , "Reference pressure" )
!
!         Request reference temperature
!         -----------------------------
          call readValue ( trim ( case_name )  , "Reference Temperature" , Setup % Temperature_ref ) 
          call Setup_CheckWithDefault( Setup % Temperature_ref , 273.15_RP , "Reference Temperature" )
!
!         Request Reynolds length
!         -----------------------
          call readValue ( trim ( case_name )  , "Reynolds length"       , Setup % reynolds_length ) 
          call Setup_CheckWithDefault( Setup % reynolds_length , 1.0_RP , "Reynolds length" )
!
!         Request Reynolds number
!         -----------------------
          call readValue ( trim ( case_name )  , "Reynolds number"       , Setup % reynolds_number ) 
          call Setup_CheckWithError( Setup % reynolds_number , "Reynolds number" )
!
!         Request Prandtl number
!         ----------------------
          call readValue ( trim ( case_name )  , "Prandtl number"        , Setup % prandtl_number  ) 
          call Setup_CheckWithDefault( Setup % prandtl_number , 0.72_RP , "Prandtl number" )
!
!         Request Mach number
!         -------------------
          call readValue ( trim ( case_name )  , "Mach number"           , Setup % Mach_number     ) 
          call Setup_CheckWithError( Setup % mach_number , "Mach number" )
!
!         Request Interpolation nodes
!         --------------------------- 
          call readValue ( trim ( case_name )  , "Interpolation nodes"   , interp_nodes ) 
          call Setup_CheckWithDefault( interp_nodes , "Legendre-Gauss" , "Interpolation nodes" )
          
          if ( trim(interp_nodes) .eq. "Legendre-Gauss" ) then
            Setup % nodes = LG
          elseif ( trim(interp_nodes) .eq. "Legendre-Gauss-Lobatto" ) then
            Setup % nodes = LGL
          else
            print*, "Unknown option for the interpolation nodes."
            print*, "Options available are:"
            print*, "   * Legendre-Gauss"
            print*, "   * Legendre-Gauss-Lobatto"
            errorMessage(STD_OUT)
            stop 
          end if
!
!         Request polynomial order
!         ------------------------
          call readValue ( trim ( case_name )  , "Default polynomial order" , Setup % N ) 
          call Setup_CheckWithError( Setup % N , "Default polynomial order" )
!
!         Request initial condition
!         -------------------------
          call readValue ( trim ( case_name )  , "Initial condition" , Setup % IC ) 
          call Setup_CheckWithError( Setup % IC , "Initial condition" )
!
!         Request inviscid formulation
!         ----------------------------
          call readValue ( trim ( case_name )  , "Inviscid discretization" , Setup % inviscid_discretization ) 
          call Setup_CheckWithDefault( Setup % inviscid_discretization , "Standard" , "Inviscid discretization" )
!
!         Request inviscid formulation
!         ----------------------------
          call readValue ( trim ( case_name )  , "Inviscid formulation" , inviscid_form ) 
          call Setup_CheckWithDefault( inviscid_form , "Green form" , "Inviscid formulation" )
      
          if ( trim(inviscid_form) .eq. "Green form" ) then
            Setup % inviscid_formulation = FORMI
          elseif ( trim(inviscid_form) .eq. "Divergence form" ) then
            Setup % inviscid_formulation = FORMII
          else
            print*, "Unknown option for the inviscid formulation."
            print*, "Options available are:"
            print*, "   * Green form"
            print*, "   * Divergence form"
            errorMessage(STD_OUT)
            stop 
          end if
!
!         Request Riemann flux
!         --------------------
          call readValue ( trim ( case_name )  , "Inviscid Riemann solver"        , Setup % inviscid_flux          ) 
          call Setup_CheckWithDefault( Setup % inviscid_flux , "Roe" , "Inviscid Riemann solver" )
!
!         Number of integration points
!         ----------------------------
          call readValue ( trim ( case_name )  , "Number of integration points" , Setup % integration_points     ) 
          
          if ( trim(Setup % inviscid_discretization) .eq. "Over-integration" ) then
            call Setup_CheckWithError( Setup % integration_points , "Number of integration points" ) 
          else
            call Setup_CheckWithDefault( Setup % integration_points , Setup % N , "Number of integration points" ,  .false. ) 
          end if
!
!         Viscous discretization
!         ----------------------
          call readValue ( trim ( case_name )  , "Viscous discretization"             , Setup % viscous_discretization ) 
          call Setup_CheckWithDefault( Setup % viscous_discretization , "BR1" , "Viscous discretization" )
         
          if ( Setup % viscous_discretization .eq. "Interior-penalty" ) then
            Setup % viscous_discretization = "IP"
          end if
!
!         Select the interior penalty method type
!         ---------------------------------------
          call readValue ( trim ( case_name )  , "Interior penalty method"          , Setup % IPMethod          ) 

          if ( Setup % viscous_discretization .eq. "IP" ) then
            call Setup_CheckWithDefault( Setup % IPMethod , "SIPG" , "Interior penalty method" )
          end if
!
!         Request the Jumps penalty parameter
!         -----------------------------------
          call readValue ( trim ( case_name )  , "Jumps penalty parameter"          , Setup % sigma0IP          ) 

          if ( Setup % viscous_discretization .eq. "IP" ) then
            call Setup_CheckWithDefault( Setup % sigma0IP , 1.0_RP , "Jumps penalty parameter" )
          end if
!
!         Request the Gradients jumps penalty parameter
!         ---------------------------------------------
          call readValue ( trim ( case_name )  , "Gradient jumps penalty parameter" , Setup % sigma1IP          ) 

          if ( Setup % viscous_discretization .eq. "IP" ) then
            call Setup_CheckWithDefault( Setup % sigma1IP , 0.0_RP , "Gradient jumps penalty parameter" )
          end if
!
!         Request the integration mode
!         ----------------------------
          call readValue ( trim ( case_name )  , "Integration mode"                 , Setup % integrationMode   ) 
          call Setup_CheckWithDefault( Setup % integrationMode , "Steady" , "Integration mode" )
!
!         Request the artificial dissipation
!         ----------------------------------
          call readValue ( trim ( case_name ) , "Artificial dissipation (0/1)" , artificialDissipation )
          call Setup_CheckWithDefault ( artificialDissipation , 0 , "Artificial dissipation (0/1)" )
          if ( artificialDissipation .eq. 1 ) then
            Setup % artificialDissipation = .true.

          else
            Setup % artificialDissipation = .false.

          end if 

          call readValue ( trim ( case_name ) , "Artificial dissipation intensity" , Setup % artificialDissipationIntensity )
          call Setup_CheckWithDefault ( Setup % artificialDissipationIntensity , 1.0_RP , "Artificial dissipation intensity" )

          call readValue ( trim ( case_name ) , "Artificial dissipation indicator" , Setup % artificialDissipationIndicator )
          call Setup_CheckWithDefault ( Setup % artificialDissipationIndicator , "Jumps-based" , "Artificial dissipation indicator" )

          call readValue ( trim ( case_name ) , "Artificial dissipation type" , Setup % artificialDissipationType )
          call Setup_CheckWithDefault ( Setup % artificialDissipationType , "Physical" , "Artificial dissipation type" )
!
!         Request the integration scheme          
!         ------------------------------
          call readValue ( trim ( case_name )  , "Integration scheme"               , Setup % integrationMethod ) 
          call Setup_CheckWithDefault( Setup % integrationMethod , "Williamson RK3" , "Integration scheme" )
!
!         Request the CFL number
!         ----------------------
          call readValue ( trim ( case_name )  , "CFL Number"                       , Setup % Ccfl              ) 
          call Setup_CheckWithDefault( Setup % Ccfl , 0.1_RP , "CFL Number" )
!
!         Request the time step
!         ---------------------
          call readValue ( trim ( case_name )  , "Time step"                        , Setup % dt                ) 
          call Setup_CheckWithDefault( Setup % dt , 0.01_RP , "Time step" ) 
!
!         Request the simulation time
!         ---------------------------
          call readValue ( trim ( case_name )  , "Simulation time"                  , Setup % simulationTime    ) 
          call Setup_CheckWithDefault( Setup % simulationTime , 1.0_RP , "Simulation time" ) 
!
!         Request the residual convergence target
!         ---------------------------------------
          call readValue ( trim ( case_name )  , "Residual target"                  , Setup % residualTarget    ) 
          call Setup_CheckWithDefault( Setup % residualTarget , 1.0e-12_RP , "Residual target" ) 
!
!         Request the number of iterations
!         --------------------------------
          call readValue ( trim ( case_name )  , "Number of iterations"             , Setup % no_of_iterations  ) 
          call Setup_CheckWithDefault( Setup % no_of_iterations , 1000 , "Number of iterations" ) 
!
!         Autosave interval
!         -----------------
          call readValue ( trim ( case_name )  , "Autosave interval"                , Setup % AutosaveInterval  ) 
          call Setup_CheckWithDefault( Setup % AutosaveInterval , 0 , "Autosave interval" ) 
!
!         Output interval
!         ---------------         
          call readValue ( trim ( case_name )  , "Output interval"                  , Setup % Output_Interval   ) 
          call Setup_CheckWithDefault( Setup % Output_Interval , 1 , "Output interval" ) 
!
!         Save variables
!         --------------
          call readValue ( trim ( case_name )  , "Save variables"                   , Setup % saveVariables     ) 
          call Setup_CheckWithDefault( Setup % saveVariables , "rho_u_v_p" , "Save variables" ) 
!
!         Restart file
!         ------------
          call readValue ( trim ( case_name )  , "Restart file"                     , Setup % restart_file      ) 

          if ( Setup % IC .eq. "Restart" ) then
            call Setup_CheckWithError( Setup % restart_file , "Restart file" ) 
          end if
!
!         Solution file
!         -------------
          call readValue ( trim ( case_name )  , "Solution file"                    , Setup % solution_file     ) 
          call Setup_CheckWithDefault( Setup % solution_file , "Solution.HiORst" , "Solution file" ) 

          pos = index( trim(Setup % solution_file) , ".HiORst" )  

          if (pos .eq. 0) then
            Setup % solution_file = trim(Setup % solution_file) // ".HiORst"
          end if
!
!         Output file type
!         ----------------
          call readValue ( trim ( case_name )  , "Output file type"                 , Setup % outputType        ) 
          call Setup_CheckWithDefault( Setup % outputType , "Interpolated" , "Output file type" ) 
!
!         Number of representation points
!         -------------------------------
          call readValue ( trim ( case_name )  , "Number of representation points"  , Setup % no_of_plotPoints  ) 
          call Setup_CheckWithDefault( Setup % no_of_plotPoints , 2 * Setup % N , "Number of representation points" ) 
!
!         Export format
!         -------------
          call readValue ( trim ( case_name )  , "Export format"  , Setup % exportFormat  ) 
          call Setup_CheckWithDefault( Setup % exportFormat , "Tecplot" , "Export format" ) 

      end subroutine Setup_Initialization

      subroutine Setup_SetInitialTime( self , t , iter ) 
         implicit none
         class(Setup_t)          :: self
         real(kind=RP)           :: t
         integer                 :: iter

         self % initialTime = t
         self % initialIteration = iter

      end subroutine Setup_SetInitialTime
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!                 CHECK SUBROUTINES
!                 -----------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Setup_CheckWithError_Character( variable , label )
         implicit none
         character(len=*),    intent(in)     :: variable
         character(len=*),    intent(in)     :: label

         if ( len_trim(variable) .eq. 0 ) then
            print*, 'Variable "',trim(adjustl(label)),'" was not found in the case file.'
            stop "Fatal error"
         end if

       end subroutine Setup_CheckWithError_Character

      subroutine Setup_CheckWithError_Integer( variable , label )
         implicit none
         integer,          allocatable,  intent(in)     :: variable
         character(len=*),               intent(in)     :: label

         if ( .not. allocated(variable) ) then
            print*, 'Variable "',trim(adjustl(label)),'" was not found in the case file.'
            stop "Fatal error"
         end if

       end subroutine Setup_CheckWithError_Integer

      subroutine Setup_CheckWithError_Real( variable , label )
         implicit none
         real(kind=RP),    allocatable,  intent(in)     :: variable
         character(len=*),               intent(in)     :: label

         if ( .not. allocated(variable) ) then
            print*, 'Variable "',trim(adjustl(label)),'" was not found in the case file.'
            stop "Fatal error"
         end if

       end subroutine Setup_CheckWithError_Real

      subroutine Setup_CheckWithDefault_Character( variable , default_value , label )
         implicit none
         character(len=*),    intent(inout)  :: variable
         character(len=*),    intent(in)     :: default_value
         character(len=*),    intent(in)     :: label

         if ( len_trim(variable) .eq. 0 ) then
            print*, 'Variable "',trim(adjustl(label)),'" was not found in the case file.'
            print*, '      >> Assigned default value: ' , trim(adjustl(default_value))
            variable = trim(adjustl(default_value))
         end if

       end subroutine Setup_CheckWithDefault_Character

      subroutine Setup_CheckWithDefault_Integer( variable , default_value , label  , verbose)
         implicit none
         integer, allocatable,   intent(inout)  :: variable
         integer,                intent(in)     :: default_value
         character(len=*),       intent(in)     :: label
         logical, optional    ,  intent(in)     :: verbose
         logical                                :: verb

         if ( present(verbose) ) then
            verb = verbose
         else
            verb = .true.
         end if

         if ( .not. allocated(variable) ) then
            if ( verb ) then
               print*, 'Variable "',trim(adjustl(label)),'" was not found in the case file.'
               print*, '      >> Assigned default value: ' , default_value
            end if
            allocate( variable )
            variable = default_value
         end if

       end subroutine Setup_CheckWithDefault_Integer

      subroutine Setup_CheckWithDefault_Real( variable , default_value , label )
         implicit none
         real(kind=RP), allocatable,   intent(inout)  :: variable
         real(kind=RP),                intent(in)     :: default_value
         character(len=*),       intent(in)     :: label

         if ( .not. allocated(variable) ) then
            print*, 'Variable "',trim(adjustl(label)),'" was not found in the case file.'
            print*, '      >> Assigned default value: ' , default_value
            allocate( variable )
            variable = default_value
         end if

       end subroutine Setup_CheckWithDefault_Real
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module Setup_class
