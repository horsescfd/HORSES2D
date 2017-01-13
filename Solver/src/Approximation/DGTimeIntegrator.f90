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

         
         if (Setup % IntegrationMode .eq. STEADY) then
            Integrator % mode = STEADY
            Integrator % dt = Setup % dt
            Integrator % no_of_iterations = Setup % no_of_iterations 
            Integrator % output_interval = Setup % output_interval
            Integrator % t_end = (Setup % dt) * (Setup % no_of_iterations)
   
         elseif (Setup % IntegrationMode .eq. TRANSIENT) then
            Integrator % mode = TRANSIENT
            Integrator % dt = Setup % dt
            Integrator % no_of_iterations = ceiling( Integrator % t_end / Setup % dt )
            Integrator % output_interval = Setup % output_interval
            Integrator % t_end = Setup % no_of_iterations * Setup % dt
            
         else
            Stop "Stopped."

         end if

         Integrator % t = Setup % initialTime

         if ( trim(Setup % integrationMethod) .eq. "Explicit-Euler" ) then
            Integrator % TimeStep  => TimeIntegrator_ExplicitEuler

         elseif ( trim(Setup % integrationMethod) .eq. "RK3" ) then
            Integrator % TimeStep => TimeIntegrator_LowStorageRK3

         else
            write(STD_OUT , *) "Method " , trim (Setup % integrationMethod) , " not implemented yet."

         end if
           

      end function TimeIntegrator_NewTimeIntegrator

      subroutine TimeIntegrator_Integrate( self , mesh , Storage )
         use Storage_module
         implicit none
         class(TimeIntegrator_t)          :: self
         class(QuadMesh_t)                  :: mesh
         class(Storage_t)                 :: Storage
         integer                          :: iter

         do iter = 1 , self % no_of_iterations

            self % iter = iter

            call self % TimeStep( mesh , self % dt , Storage)
            self % t    = self % t + self % dt

            if ( iter .eq. 1 ) then
               call self % Display( mesh ) 
            elseif ( mod(iter , self % output_interval) .eq. 0 ) then
               call self % Display( mesh ) 
            end if

            call self % Autosave( mesh )


         end do


      end subroutine TimeIntegrator_Integrate
      
      subroutine TimeIntegrator_Display( self, mesh )
         use Physics
         implicit none  
         class(TimeIntegrator_t)          :: self
         class(QuadMesh_t)                  :: mesh
         integer, parameter               :: ShowLabels = 10
         integer, save                    :: shown = 0
         real(kind=RP)                    :: residuals(NEC)

         if ( mod( shown , ShowLabels) .eq. 0 ) then     ! Show labels
            write(STD_OUT , '(/)')
            write(STD_OUT , '(/)')
            write(STD_OUT , '(20X,A20,5X,A20,5X,A20,5X,A20,5X,A20,5X,A20)') "Iteration" , "time" , "continuity" , "x-momentum" , "y-momentum", "energy"
            write(STD_OUT , '(20X,A20,5X,A20,5X,A20,5X,A20,5X,A20,5X,A20)') "---------" , "----" , "----------" , "----------" , "----------", "------"
         end if
         shown = shown + 1

         residuals = mesh % computeResiduals()

         write(STD_OUT , '(20X,I20,2X,A,2X,ES20.8,2X,A,2X,ES20.8,2X,A,2X,ES20.8,2X,A,2X,ES20.8,2X,A,2X,ES20.8)') self % iter ,"|", self % t ,"|", residuals(IRHO) , "|" , residuals(IRHOU) , &
                                          "|", residuals(IRHOV) , "|" , residuals(IRHOE)
      end subroutine TimeIntegrator_Display

      subroutine TimeIntegrator_Describe( self )
         implicit none
         class(TimeIntegrator_t)          :: self

         write(STD_OUT , *) "Time integrator description: "
         if (self % mode .eq. STEADY) write(STD_OUT , '(20X,A,A)') "Mode: " , "steady"
         if (self % mode .eq. TRANSIENT) write(STD_OUT , '(20X,A,A)') "Mode: " , "transient"
         write(STD_OUT ,'(20X,A,E10.3)') "Time step dt: " , self % dt
         write(STD_OUT , '(20X,A,I0)') "Number of iterations: " , self % no_of_iterations
         write(STD_OUT , '(20X,A,E10.3)') "Final simulation time: " , self % t_end
         
      end subroutine TimeIntegrator_Describe
   
      subroutine TimeIntegrator_Autosave( self , mesh )
         class(TimeIntegrator_t)          :: self 
         class(QuadMesh_t)                  :: mesh

      end subroutine TimeIntegrator_Autosave
!
!     *******************************************************************************************************
!           Integration methods library
!     *******************************************************************************************************
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

      subroutine TimeIntegrator_LowStorageRK3( mesh , dt , Storage )
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
         

      end subroutine TimeIntegrator_LowStorageRK3

end module DGTimeIntegrator
