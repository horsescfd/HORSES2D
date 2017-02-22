module DGTimeIntegrator
   use SMConstants
   use DGSpatialDiscretizationMethods
   use QuadMeshClass
   use MonitorsClass

   private
   public TimeIntegrator_t , NewTimeIntegrator
!
!                                *******************************
   integer, parameter         :: STR_LEN_TIMEINTEGRATOR    = 128
   integer, parameter         :: estimateTimeStep_interval = 100
!                                *******************************
!
!
!                                ********
   type(Monitor_t)            :: Monitors
!                                ********
!

   type TimeIntegrator_t
      integer                               :: iter
      real(kind=RP)                         :: t
      real(kind=RP)                         :: dt
      real(kind=RP)                         :: t_end
      real(kind=RP)                         :: Ccfl
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
         procedure ::  EstimateTimeStep => TimeIntegrator_EstimateTimeStep
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
      function TimeIntegrator_NewTimeIntegrator(mesh) result (Integrator)
         use Setup_class
         implicit none
         type(TimeIntegrator_t)           :: Integrator
         type(QuadMesh_t)                 :: mesh

         
         if ( trim(Setup % IntegrationMode) .eq. "Steady" ) then
            Integrator % mode = STEADY
            Integrator % no_of_iterations  = Setup % no_of_iterations
            Integrator % output_interval   = Setup % output_interval
            Integrator % autosave_interval = Setup % autosaveInterval
            Integrator % initial_iteration = Setup % initialIteration
            Integrator % Ccfl              = Setup % Ccfl
            Integrator % dt                = Setup % dt
   
         elseif ( trim(Setup % IntegrationMode) .eq. "Transient") then
            print*, "Do not use this mode."
            stop "Stopped"
            Integrator % mode = TRANSIENT
            Integrator % no_of_iterations  = 0
            Integrator % output_interval   = Setup % output_interval
            Integrator % autosave_interval = Setup % autosaveInterval
            Integrator % Ccfl              = Setup % Ccfl
            Integrator % t_end             = Setup % simulationTime
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

         call Integrator % EstimateTimeStep( mesh )

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
!
!        Create monitors
!        ---------------
         Monitors = ConstructMonitors( mesh )
           

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

            call Monitors % UpdateValues ( mesh , self % t , self % iter)

            if ( iter .eq. self % initial_iteration + 1 ) then
               call self % Display( mesh ) 
            elseif ( mod(iter , self % output_interval) .eq. 0 ) then
               call self % Display( mesh ) 
            end if

            if ( mod(iter , self % autosave_interval) .eq. 0 ) then
               call self % Autosave( Storage , mesh )
            end if

            if ( mod(iter , estimateTimeStep_interval) .eq. 0 ) then
               call self % EstimateTimeStep( mesh )
            end if

            call Monitors % WriteToFile()

         end do
!
!        Save solution file
!        ------------------
         call sleep(2)
         self % iter = self % initial_iteration + self % no_of_iterations
         call self % Autosave( Storage , mesh , trim(Setup % solution_file) ) 

         call Monitors % WriteToFile ( force = .true. )

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
         real(kind=RP), allocatable, save :: G(:)
         integer                    :: m 
         integer, parameter         :: N_STAGES = 3
         real(kind=RP), parameter   :: am(3) = [0.0_RP , -5.0_RP / 9.0_RP , -153.0_RP / 128.0_RP]
         real(kind=RP), parameter   :: bm(3) = [0.0_RP , 1.0_RP / 3.0_RP  , 3.0_RP / 4.0_RP ]
         real(kind=RP), parameter   :: gm(3) = [1.0_RP / 3.0_RP , 15.0_RP / 16.0_RP , 8.0_RP / 15.0_RP ]
         
         if (.not. allocated(G) ) allocate ( G , source  = Storage % QDot )

         do m = 1 , N_STAGES
!
!           Compute time derivative
!           -----------------------
            call DGSpatial_ComputeTimeDerivative( mesh )

            if (m .eq. 1) then
               G = Storage % QDot

            else
               G = am(m) * G + Storage % QDot

            end if

            Storage % Q = Storage % Q + gm(m) * dt * G

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
         real(kind=RP), save, allocatable :: G(:)
         integer                    :: m 
         integer, parameter         :: N_STAGES = 5
         real(kind=RP), parameter  :: am(N_STAGES) = [0.0_RP , -0.4178904745_RP, -1.192151694643_RP , -1.697784692471_RP , -1.514183444257_RP ]
         real(kind=RP), parameter  :: gm(N_STAGES) = [0.1496590219993_RP , 0.3792103129999_RP , 0.8229550293869_RP , 0.6994504559488_RP , 0.1530572479681_RP]
         
         if (.not. allocated(G) ) allocate ( G , source  = Storage % QDot )

         do m = 1 , N_STAGES
!
!           Compute time derivative
!           -----------------------
            call DGSpatial_ComputeTimeDerivative( mesh )

            if (m .eq. 1) then
               G = dt * Storage % QDot

            else
               G = am(m) * G + dt * Storage % QDot

            end if

            Storage % Q = Storage % Q + gm(m) * G

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
         real(kind=RP)                    :: residuals(NCONS)

         if ( mod( shown , ShowLabels) .eq. 0 ) then     ! Show labels
            write(STD_OUT , '(/)')
            write(STD_OUT , '(/)')

            write(STD_OUT , '(A10,3X,A10,3X,A10,3X,A10,3X,A10,3X,A10)', advance = "no" ) "Iteration" , "time" , "continuity" , "x-momentum" , "y-momentum", "energy"
            call Monitors % WriteLabel
            write(STD_OUT , * ) 

            write(STD_OUT , '(A10,3X,A10,3X,A10,3X,A10,3X,A10,3X,A10)' , advance = "no" ) "---------" , "--------" , "----------" , "----------" , "----------", "--------"
            call Monitors % WriteUnderlines
            write(STD_OUT , * ) 
         end if
         shown = shown + 1

         residuals = mesh % computeResiduals()

         write(STD_OUT , '(I10,1X,A,1X,ES10.3,1X,A,1X,ES10.3,1X,A,1X,ES10.3,1X,A,1X,ES10.3,1X,A,1X,ES10.3)' , advance = "no") self % iter ,"|", self % t ,"|", residuals(IRHO) , "|" , residuals(IRHOU) , &
                                          "|", residuals(IRHOV) , "|" , residuals(IRHOE)
         call Monitors % WriteValues
         write(STD_OUT,*)
         
      end subroutine TimeIntegrator_Display

      subroutine TimeIntegrator_Describe( self )
         use Headers
         use Setup_class
         implicit none
         class(TimeIntegrator_t)          :: self

         write(STD_OUT , *)
         call Section_header("Time integrator description")
         write(STD_OUT , *)
         write(STD_OUT , '(30X,A,A30,A)') "-> ","Method: " , trim(Setup % integrationMethod)
         if (self % mode .eq. STEADY) write(STD_OUT , '(30X,A,A30,A)') "-> ","Mode: " , "Steady"
         if (self % mode .eq. TRANSIENT) write(STD_OUT , '(30X,A,A30,A)') "-> ","Mode: " , "Transient"
         write(STD_OUT ,'(30X,A,A30,ES10.3)') "-> ","Initial time step dt: " , self % dt
         write(STD_OUT , '(30X,A,A30,I0)') "-> ","Number of iterations: " , self % no_of_iterations
         write(STD_OUT , '(30X,A,A30,ES10.3)') "-> ", "Estimated simulation time: " , self % dt * self % no_of_iterations 
         
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
         call NetCDF_putDimension( trim(fileName) , "NCONS" , NCONS )
         call NetCDF_putDimension( trim(fileName) , "no_of_elements" , mesh % no_of_elements )

         call NetCDF_putVariable( trim(fileName) , "t" , ["one"] , [self % t] )
         call NetCDF_putVariable( trim(fileName) , "Q" , ["NDOF"] , Storage % Q )
         call NetCDF_putVariable( trim(fileName) , "iter" , ["one"] , [self % iter] )
         
         write(STD_OUT , '(A,/)' ) "..  Saved"

      end subroutine TimeIntegrator_Autosave
!
!////////////////////////////////////////////////////////////////////////////
!
!           TIME STEP ESTIMATOR
!           -------------------
!///////////////////////////////////////////////////////////////////////////
!
      subroutine TimeIntegrator_EstimateTimeStep( self , mesh )
         use Setup_Class
         use Physics
         use NodeClass
         use QuadElementClass
         implicit none
         class(TimeIntegrator_t)                   :: self
         class(QuadMesh_t)                         :: mesh
         integer                                   :: eID
         real(kind=RP)                             :: dt = 0.0_RP
         real(kind=RP)                             :: dx
         real(kind=RP)                             :: umax = 0.0_RP , amax = 0.0_RP
         integer                                   :: iXi , iEta
         class(QuadElement_t), pointer             :: e
         class(Node_p),        pointer             :: nodes(:)

         dt = Setup % dt

         do eID = 1 , mesh % no_of_elements

            associate( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )
            e => mesh % elements(eID)
            nodes => mesh % elements(eID) % nodes

            do iXi = 0 , mesh % elements(eID) % spA % N
               do iEta = 0 , mesh % elements(eID) % spA % N

                  umax = max(umax , norm2( e % Q(iXi,iEta,IRHOU:IRHOV)) / e % Q(iXi,iEta,IRHO))
                  amax = max(amax , sqrt( gamma * gm1 * (e % Q(iXi,iEta,IRHOE) /  e % Q(iXi,iEta,IRHO) - 0.5_RP * ( e % Q(iXi,iEta,IRHOU)/ e % Q(iXi,iEta,IRHO))**2.0_RP &
                                                                                                       - 0.5_RP * ( e % Q(iXi,iEta,IRHOV)/ e % Q(iXi,iEta,IRHO))**2.0_RP )))
             
               end do
            end do

   
            dx = min( norm2(nodes(1) % n % X-nodes(2) % n % X) , norm2(nodes(2) % n % X - nodes(3) % n % X) , norm2(nodes(3) % n % X -nodes(4) % n % X) , norm2(nodes(1) % n % X - nodes(4) % n % X))
            dt = min( dt , Dimensionless % Mach * sqrt(Thermodynamics % gamma) * self % Ccfl * dx / (umax + amax) / (mesh % elements(eID) % spA % N+1) )

            end associate

         end do

         self % dt = dt

      end subroutine TimeIntegrator_EstimateTimeStep

end module DGTimeIntegrator
