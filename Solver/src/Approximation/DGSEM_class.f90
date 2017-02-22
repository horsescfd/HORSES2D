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
    implicit none
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
         
        subroutine DGSEM_construct( self ,  meshFile )
            use Setup_class
            use QuadElementClass
            implicit none
            class(DGSEM_t)                                               :: self
            class(MeshFile_t)                                            :: meshFile
!           ----------------------------------------------------------------------------------------------
            integer                                                      :: total_polyorder
            integer                                                      :: eID
            integer                                                      :: iBC
            integer                                                      :: fID
            class(Edge_t), pointer                                       :: face
!
!           Allocate memory for Q , W , QDot , dQ , and F
!              The sizes are the following:
!                 Q    -> NCONS * (N+1) * (N+1) * no_of_elements
!                 W    -> NPRIM * (N+1) * (N+1) * no_of_elements
!                 QDot -> NCONS * (N+1) * (N+1) * no_of_elements
!                 dQ   -> NGRAD * NDIM * (N+1) * (N+1) * no_of_elements
!                 F    -> NCONS * NDIM * (N+1) * (N*1) * no_of_elements (N = integration_points for over integration)
!           ---------------------------------------------------------------------------------------------------------
            allocate ( self % Storage % Q    ( NCONS *         meshFile % cumulativePolynomialOrder ( meshFile % no_of_elements )  )  ) 
            allocate ( self % Storage % W    ( NPRIM *         meshFile % cumulativePolynomialOrder ( meshFile % no_of_elements )  )  ) 
            allocate ( self % Storage % QDot ( NCONS *         meshFile % cumulativePolynomialOrder ( meshFile % no_of_elements )  )  ) 
            allocate ( self % Storage % dQ   ( NDIM  * NGRAD * meshFile % cumulativePolynomialOrder ( meshFile % no_of_elements )  )  ) 
                
            if (Setup % inviscid_discretization .eq. "Over-Integration") then
               allocate ( self % Storage % F   ( NDIM * NCONS * meshFile % no_of_elements * ( setup % integration_points + 1)**2    ) )

            else
               allocate ( self % Storage % F   ( NDIM * NCONS * meshFile % cumulativePolynomialOrder ( meshFile % no_of_elements )  )  ) 
!
            end if
!
!           Construct the spectral Integration class if Over-Integration is selected
!           ------------------------------------------------------------------------
            if (Setup % inviscid_discretization .eq. "Over-Integration") then
               allocate( self % spI )
               call self % spI % init( Setup % integration_points , Setup % nodes )

            end if
!
!           Construct the mesh
!           ------------------
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
            use Tecplot
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
               solutionpltName = Setup % solution_file(1:pos-1) // ".plt"
            end if
   
            call ExportToTecplot( self % mesh , trim(solutionpltname))  

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
