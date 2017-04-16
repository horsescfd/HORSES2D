!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      DGSEM_class.f90
!      Created: 2015-09-26 12:05:17
!      REV:     2016-02-24 21:16:23
!      By: Juan MANZANERO
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!   *****************
    module DGSEM_class
!   *****************

!   ------------
!   Modules used
!   ------------
    use SMConstants
    use Physics
    use NodesAndWeights_class
    use QuadMeshClass
    use MeshFileClass
    use DGSpatialDiscretizationMethods
    use DGTimeIntegrator
    use Storage_module
    use DGBoundaryConditions
    use Plotter
    implicit none
!
#include "Defines.h"
!
    private
    public DGSEM_t , DGSEM_Initialize
!
!                                   *******************
    integer, parameter           :: STR_LEN_DGSEM = 128
!                                   *******************
!
!   ---------------------
!   DGSEM type DEFINITION
!   ---------------------
!
    type DGSEM_t
        type(QuadMesh_t)                    :: mesh
        type(NodalStorage)                  :: spA             ! Interpolation nodes and weights structure
        class(NodesAndWeights_t),   pointer :: spI => NULL()   ! Integration nodes and weights structure
        type(Storage_t)                     :: Storage
        class(BoundaryCondition_t), pointer :: BoundaryConditions(:)
        type(TimeIntegrator_t)              :: Integrator
        class(Plotter_t), allocatable       :: Plotter
        contains
            procedure       :: Construct => DGSEM_construct
            procedure       :: SetInitialCondition => DGSEM_SetInitialCondition
            procedure       :: Integrate => DGSEM_Integrate
            procedure       :: LoadRestartFile => DGSEM_LoadRestartFile
    end type DGSEM_t

!
!   ========
    contains
!   ========
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        function DGSEM_Initialize() result(sem)
            type(DGSEM_t)         :: sem
!
!           Initialize mesh
!           ---------------
            sem % mesh = InitializeMesh()
!
!           Initialize Spectral Approximation object
!           ----------------------------------------
            sem % spA  = newNodalStorage()
!
!           Initialize Storage
!           ------------------ 
            sem % Storage = newStorage()
             
        end function DGSEM_Initialize
         
        subroutine DGSEM_construct( self )
            use Setup_class
            use QuadElementClass
            implicit none
            class(DGSEM_t)   :: self
!
!           ---------------
!           Local variables
!           ---------------
!
            type(MeshFile_t) :: meshFile
            integer          :: totalPolynomialOrder
!
!           Read the mesh file
!           ------------------
            call meshFile % Read
!
!           Allocate memory for the solution, time derivative, and gradients
!           ----------------------------------------------------------------
            totalPolynomialOrder = meshFile % cumulativePolynomialOrder( meshFile % no_of_elements )
            call self % Storage % AllocateMemory( totalPolynomialOrder )
!
!           Construct the spectral Integration class if Over-Integration is selected
!           ------------------------------------------------------------------------
            if (Setup % inviscid_discretization .eq. "Over-Integration") then
               allocate( self % spI )
               call self % spI % init( Setup % integration_points , Setup % nodes )
            else
               self % spI => NULL()
            end if
!
!           Construct the spectral element mesh object
!           ------------------------------------------
            call self % mesh % constructFromFile(meshFile , self % spA , self % Storage , self % spI)
!
!           Prepare the spectral Approximation structures generated for Over-Integration
!           ----------------------------------------------------------------------------
            if (Setup % inviscid_discretization .eq. "Over-Integration") then
!          
!              Set the interpolation matrices
!              ------------------------------
               call self % spA % computeInterpolationMatrices( self % spI )

            end if
!
!           Construct the domain zones
!           --------------------------              
            call self % mesh % ConstructZones( meshFile )
!
!           Initialize Inviscid and Viscous discretization methods
!           ------------------------------------------------------
            call DGSpatial_Initialization() 
!
!           Construct plotter and Export the mesh            
!           -------------------------------------
            call ConstructPlotter( self % Plotter )
            call self % Plotter % ExportMesh( self % mesh , Setup % mesh_file )   
!
!           Destruct the mesh file object
!           -----------------------------
            call meshFile % Destruct
!
!           Set the initial condition to all flow variables
!           -----------------------------------------------
            call self % SetInitialCondition()
            call self % Plotter % Export( self % mesh , './RESULTS/InitialCondition')     
            
        end subroutine DGSEM_construct
            
        subroutine DGSEM_SetInitialCondition( self , verbose )
!
!           *********************************************************
!                 This routine sets the initial condition from one 
!              previously selected from the library, or from the
!              Restart file *.HiORst
!           *********************************************************
!              
            use InitialConditions
            use Setup_class
            implicit none
            class(DGSEM_t)                   :: self
            logical, optional                :: verbose

            if ( Setup % IC .eq. "Restart" ) then
               call self % loadRestartFile( trim ( Setup % Restart_file ) ) 
            else
               call self % mesh % SetInitialCondition ()
            end if
!
!           Describe the initial condition
!           ------------------------------
            if ( present ( verbose ) ) then
               if ( verbose ) then
                  call InitialCondition_Describe
   
               end if
            else
               call InitialCondition_Describe

            end if

            call DGSpatial_ComputeTimeDerivative( self % mesh )

        end subroutine DGSEM_SetInitialCondition
      
        subroutine DGSEM_Integrate( self )
!
!           *********************************************************
!                 This routine constructs the Time Integrator,
!              and performs the time integration of the problem.
!           *********************************************************
!
            use Setup_Class
            implicit none
            class(DGSEM_t)                   :: self
            character(len=STR_LEN_DGSEM)     :: solutionpltName
            integer                          :: pos
!
!           Construct Time Integrator
!           -------------------------
            self % Integrator = NewTimeIntegrator( self % mesh )
            call self % Integrator % Describe()
!
!           Integrate
!           ---------
            call self % Integrator % Integrate( self % mesh , self % Storage)
!
!           Save the solution file
!           ----------------------
            pos = index( Setup % solution_file , '.HiORst')

            if ( pos .gt. 0 ) then
               solutionpltName = Setup % solution_file(1:pos-1) 
            end if
   
            call self % Plotter % Export ( self % mesh , trim(solutionpltname))  

        end subroutine DGSEM_Integrate
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!            
!              AUXILIAR FUNCTIONS
!              ------------------ 
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        subroutine DGSEM_loadRestartFile( self , fileName ) 
!     
!           ************************************************************************
!                 Loads a previous state from a Restart *.HiORst file
!           ************************************************************************
!
            use NetCDFInterface
            use Setup_class
            implicit none
            class(DGSEM_t)                :: self
            character(len=*)              :: fileName
            real(kind=RP), allocatable    :: t(:)
            integer, allocatable          :: iter(:)
            real(kind=RP), allocatable    :: Q(:)
!
!           Get the current time and iteration
!           ----------------------------------
            call NetCDF_getVariable ( trim ( fileName )  , "t"    , t    ) 
            call NetCDF_getVariable ( trim ( fileName )  , "iter" , iter ) 
            call Setup % SetInitialTime ( t(1) , iter(1)) 
!
!           Load the state vector
!           ---------------------
            call NetCDF_getVariable ( trim ( fileName )  , "Q"    , Q    ) 
            self % Storage % Q = Q
      
            deallocate( t , Q ) 
   
        end subroutine DGSEM_loadRestartFile

   end module DGSEM_class
