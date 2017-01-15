module DGTimeIntegrator
   use SMConstants
   use DGSpatialDiscretizationMethods
   use QuadMeshClass
   use FileWriting

   private
   public TimeIntegrator_t , NewTimeIntegrator

   integer, parameter         :: STR_LEN_TIMEINTEGRATOR = 128

   type TimeIntegrator_t
      integer                               :: iter
      real(kind=RP)                         :: t
      real(kind=RP)                         :: dt
      real(kind=RP)                         :: t_end
      integer                               :: no_of_iterations
      integer                               :: initial_iteration 
      integer                               :: plot_interval
      integer                               :: autosave_interval
      integer                               :: output_interval
      character(len=STR_LEN_TIMEINTEGRATOR) :: method
      integer                               :: mode
      procedure(TimeScheme), pointer, nopass :: TimeStep => NULL()
      contains
         procedure ::  Describe => TimeIntegrator_Describe
         procedure ::  Integrate => TimeIntegrator_Integrate
         procedure ::  Display  => TimeIntegrator_Display
         procedure ::  Autosave => TimeIntegrator_Autosave
   end type TimeIntegrator_t

   interface NewTimeIntegrator
      module procedure TimeIntegrator_NewTimeIntegrator 
   end interface NewTimeIntegrator

   abstract interface
      subroutine TimeScheme( mesh , dt , Storage)
         use SMConstants
         use Storage_module
         use QuadMeshClass
         implicit none
         class(QuadMesh_t)      :: mesh
         real(kind=RP)        :: dt
         class(Storage_t)     :: Storage
      end subroutine TimeScheme
   end interface

   
!  ========
   contains
!  ========
      function TimeIntegrator_NewTimeIntegrator() result (Integrator)
         use Setup_class
         implicit none
         type(TimeIntegrator_t)           :: Integrator

         
         if ( trim(Setup % IntegrationMode) .eq. "Steady" ) then
            Integrator % mode = STEADY
            Integrator % dt = Setup % dt
            Integrator % no_of_iterations = Setup % no_of_iterations 
            Integrator % output_interval = Setup % output_interval
            Integrator % autosave_interval = Setup % autosaveInterval
            Integrator % t_end = (Setup % dt) * (Setup % no_of_iterations)
            Integrator % initial_iteration = Setup % initialIteration
   
         elseif ( trim(Setup % IntegrationMode) .eq. "Transient") then
            Integrator % mode = TRANSIENT
            Integrator % dt = Setup % dt
            Integrator % no_of_iterations = ceiling( Integrator % t_end / Setup % dt )
            Integrator % output_interval = Setup % output_interval
            Integrator % autosave_interval = Setup % autosaveInterval
            Integrator % t_end = Setup % no_of_iterations * Setup % dt
            Integrator % initial_iteration = Setup % initialIteration
            
         else
            write(STD_OUT , '(/,/)') 
            write(STD_OUT , *) "Mode " , trim (Setup % integrationMode) , " not implemented yet."
            write(STD_OUT , *) "Options available are: "
            write(STD_OUT , '(5X,A)') "* Steady"
            write(STD_OUT , '(5X,A)') "* Transient"
            Stop "Stopped."

         end if

         Integrator % t = Setup % initialTime

         if ( trim(Setup % integrationMethod) .eq. "Explicit-Euler" ) then
            Integrator % TimeStep  => TimeIntegrator_ExplicitEuler

         elseif ( trim(Setup % integrationMethod) .eq. "Williamson RK3" ) then
            Integrator % TimeStep => TimeIntegrator_WilliamsonRK3

         elseif ( trim(Setup % integrationMethod) .eq. "Williamson RK5" ) then
            Integrator % TimeStep => TimeIntegrator_WilliamsonRK5

         else
            write(STD_OUT , '(/,/)') 
            write(STD_OUT , *) "Method " , trim (Setup % integrationMethod) , " not implemented yet."
            write(STD_OUT , *) "Options available are: "
            write(STD_OUT , '(5X,A)') "* Explicit-Euler"
            write(STD_OUT , '(5X,A)') "* Williamson RK3"
            write(STD_OUT , '(5X,A)') "* Williamson RK5"
            stop "Stopped."

         end if
           

      end function TimeIntegrator_NewTimeIntegrator

      subroutine TimeIntegrator_Integrate( self , mesh , Storage )
         use Storage_module
         use Setup_class
         implicit none
         class(TimeIntegrator_t)          :: self
         class(QuadMesh_t)                  :: mesh
         class(Storage_t)                 :: Storage
         integer                          :: iter

         do iter = self % initial_iteration + 1 , self % initial_iteration + self % no_of_iterations

            self % iter = iter

            call self % TimeStep( mesh , self % dt , Storage)
            self % t    = self % t + self % dt

            if ( iter .eq. self % initial_iteration + 1 ) then
               call self % Display( mesh ) 
            elseif ( mod(iter , self % output_interval) .eq. 0 ) then
               call self % Display( mesh ) 
            end if

            if ( mod(iter , self % autosave_interval) .eq. 0 ) then
               call self % Autosave( Storage , mesh )
            end if


         end do
!
!        Save solution file
!        ------------------
         call self % Autosave( Storage , mesh , trim(Setup % solution_file) ) 

      end subroutine TimeIntegrator_Integrate
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!              TIME INTEGRATION METHODS LIBRARY
!              --------------------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeIntegrator_ExplicitEuler( mesh , dt , Storage)
         use Storage_module
         implicit none
         class(QuadMesh_t)         :: mesh
         real(kind=RP)           :: dt
         class(Storage_t)        :: Storage
!
!        Compute the time derivative
!  
         call DGSpatial_computeTimeDerivative( mesh )
!
!        Perform a step in the explicit Euler method
!
         Storage % Q = Storage % Q + dt * Storage % QDot

      end subroutine TimeIntegrator_ExplicitEuler

      subroutine TimeIntegrator_WilliamsonRK3( mesh , dt , Storage )
         use Storage_module
         implicit none
         class(QuadMesh_t)          :: mesh
         real(kind=RP)              :: dt
         class(Storage_t)           :: Storage
!        -----------------------------------------
         real(kind=RP), allocatable :: G(:)
         integer                    :: m 
         integer, parameter         :: N_STAGES = 3
         real(kind=RP), parameter   :: am(3) = [0.0_RP , -5.0_RP / 9.0_RP , -153.0_RP / 128.0_RP]
         real(kind=RP), parameter   :: bm(3) = [0.0_RP , 1.0_RP / 3.0_RP  , 3.0_RP / 4.0_RP ]
         real(kind=RP), parameter   :: gm(3) = [1.0_RP / 3.0_RP , 15.0_RP / 16.0_RP , 8.0_RP / 15.0_RP ]
         
         allocate ( G , source  = Storage % QDot )

         do m = 1 , N_STAGES
!
!           Compute time derivative
!           -----------------------
            call DGSpatial_ComputeTimeDerivative( mesh )

            if (m .eq. 1) then
               G = Storage % QDot
            else
               G = am(m) * G + Storage % QDot
               Storage % Q = Storage % Q + gm(m) * dt * G
            end if

         end do 
         

      end subroutine TimeIntegrator_WilliamsonRK3

      subroutine TimeIntegrator_WilliamsonRK5( mesh , dt , Storage )
!  
!        *****************************************************************************************
!           These coefficients have been extracted from the paper: "Fourth-Order 2N-Storage
!          Runge-Kutta Schemes", written by Mark H. Carpented and Christopher A. Kennedy
!        *****************************************************************************************
!
         use Storage_module
         implicit none
         class(QuadMesh_t)          :: mesh
         real(kind=RP)              :: dt
         class(Storage_t)           :: Storage
!        -----------------------------------------
         real(kind=RP), allocatable :: G(:)
         integer                    :: m 
         integer, parameter         :: N_STAGES = 5
         real(kind=RP), parameter  :: am(N_STAGES) = [0.0_RP , -0.4178904745_RP, -1.192151694643_RP , -1.697784692471_RP , -1.514183444257_RP ]
         real(kind=RP), parameter  :: gm(N_STAGES) = [0.1496590219993_RP , 0.3792103129999_RP , 0.8229550293869_RP , 0.6994504559488_RP , 0.1530572479681_RP]
         
         allocate ( G , source  = Storage % QDot )

         do m = 1 , N_STAGES
!
!           Compute time derivative
!           -----------------------
            call DGSpatial_ComputeTimeDerivative( mesh )

            if (m .eq. 1) then
               G = dt * Storage % QDot
            else
               G = am(m) * G + dt * Storage % QDot
               Storage % Q = Storage % Q + gm(m) * G
            end if

         end do 
         

      end subroutine TimeIntegrator_WilliamsonRK5

!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!              EXTRA ROUTINES
!              --------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeIntegrator_Display( self, mesh )
         use Physics
         implicit none  
         class(TimeIntegrator_t)          :: self
         class(QuadMesh_t)                  :: mesh
         integer, parameter               :: ShowLabels = 50
         integer, save                    :: shown = 0
         real(kind=RP)                    :: residuals(NEC)

         if ( mod( shown , ShowLabels) .eq. 0 ) then     ! Show labels
            write(STD_OUT , '(/)')
            write(STD_OUT , '(/)')
            write(STD_OUT , '(10X,A20,5X,A10,5X,A10,5X,A10,5X,A10,5X,A10)') "Iteration" , "time" , "continuity" , "x-momentum" , "y-momentum", "energy"
            write(STD_OUT , '(10X,A20,5X,A10,5X,A10,5X,A10,5X,A10,5X,A10)') "---------" , "--------" , "----------" , "----------" , "----------", "--------"
         end if
         shown = shown + 1

         residuals = mesh % computeResiduals()

         write(STD_OUT , '(10X,I20,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3)') self % iter ,"|", self % t ,"|", residuals(IRHO) , "|" , residuals(IRHOU) , &
                                          "|", residuals(IRHOV) , "|" , residuals(IRHOE)
      end subroutine TimeIntegrator_Display

      subroutine TimeIntegrator_Describe( self )
         use Headers
         implicit none
         class(TimeIntegrator_t)          :: self

         write(STD_OUT , *)
         call Section_header("Time integrator description")
         write(STD_OUT , *)
         if (self % mode .eq. STEADY) write(STD_OUT , '(30X,A,A30,A)') "-> ","Mode: " , "steady"
         if (self % mode .eq. TRANSIENT) write(STD_OUT , '(30X,A,A30,A)') "-> ","Mode: " , "transient"
         write(STD_OUT ,'(30X,A,A30,ES10.3)') "-> ","Time step dt: " , self % dt
         write(STD_OUT , '(30X,A,A30,I0)') "-> ","Number of iterations: " , self % no_of_iterations
         write(STD_OUT , '(30X,A,A30,ES10.3)') "-> ", "Final simulation time: " , self % t_end
         
      end subroutine TimeIntegrator_Describe
   
      subroutine TimeIntegrator_Autosave( self , Storage , mesh , fileName_in)
         use Setup_class
         use Storage_module
         use NetCDFInterface
         use Physics
         class ( TimeIntegrator_t ) , intent  ( in )  :: self
         class ( Storage_t        ) , intent  ( in )  :: Storage
         class ( QuadMesh_t       ) , intent  ( in )  :: mesh
         character(len=*) , intent(in) , optional     :: fileName_in
!        ----------------------------------------------------------------
         character(len=STR_LEN_TIMEINTEGRATOR)            :: fileName
         integer                                          :: pos

         if ( present(fileName_in) ) then
            fileName = fileName_in
            write(STD_OUT,'(/,20X,A,A,A)',advance="no") "** Saving solution file as ", trim(fileName) , "................."

         else
            fileName = trim(Setup % solution_file)
   
            pos = index(trim(fileName) , '.HiORst' )
            
            write(fileName, '(A,A,I0,A)') fileName(1:pos-1) , "_" , self % iter , ".HiORst" 
            
            write(STD_OUT,'(/,20X,A,A,A)',advance="no") "** Saving restart file as ", trim(fileName) , "................."
   
         end if

         call NetCDF_CreateFile ( trim(fileName) ) 
         call NetCDF_putDimension( trim(fileName) , "one" , 1 )
         call NetCDF_putDimension( trim(fileName) , "NDOF" , size(Storage % Q) )
         call NetCDF_putDimension( trim(fileName) , "N" , Setup % N )
         call NetCDF_putDimension( trim(fileName) , "NEC" , NEC )
         call NetCDF_putDimension( trim(fileName) , "no_of_elements" , mesh % no_of_elements )

         call NetCDF_putVariable( trim(fileName) , "t" , ["one"] , [self % t] )
         call NetCDF_putVariable( trim(fileName) , "Q" , ["NDOF"] , Storage % Q )
         call NetCDF_putVariable( trim(fileName) , "iter" , ["one"] , [self % iter] )
         
         write(STD_OUT , '(A,/)' ) "..  Saved"


      end subroutine TimeIntegrator_Autosave
!






end module DGTimeIntegrator
