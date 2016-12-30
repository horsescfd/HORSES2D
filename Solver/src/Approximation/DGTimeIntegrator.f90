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
            Integrator % t_end = (Setup % dt) * (Setup % no_of_iterations)
   
         elseif (Setup % IntegrationMode .eq. TRANSIENT) then
            Integrator % mode = TRANSIENT
            Integrator % dt = Setup % dt
            Integrator % no_of_iterations = ceiling( Integrator % t_end / Setup % dt )
            Integrator % t_end = Setup % no_of_iterations * Setup % dt
            
         else
            Stop "Stopped."

         end if

         Integrator % t = Setup % initialTime

         if ( trim(Setup % integrationMethod) .eq. "Explicit-Euler" ) then
            Integrator % TimeStep  => TimeIntegrator_ExplicitEuler

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

            call self % Display( mesh , Storage ) 
            call self % Autosave( mesh )


         end do


      end subroutine TimeIntegrator_Integrate
      
      subroutine TimeIntegrator_Display( self, mesh, Storage )
         use Storage_module
         implicit none  
         class(TimeIntegrator_t)          :: self
         class(QuadMesh_t)                  :: mesh
         class(Storage_t)                 :: Storage
         integer, parameter               :: ShowLabels = 10

         if ( mod( self % iter-1 , ShowLabels) .eq. 0 ) then     ! Show labels
            write(STD_OUT , '(/)')
            write(STD_OUT , '(/)')
            write(STD_OUT , '(20X,A20,5X,A20,5X,A20)') "Iteration" , "Time" , "Residual"
            write(STD_OUT , '(20X,A20,5X,A20,5X,A20)') "---------" , "----" , "--------"
         end if

         write(STD_OUT , '(20X,I20,2X,A,2X,F20.8,2X,A,2X,F20.8)') self % iter ,"|", self % t ,"|", maxval(abs(Storage % QDot) )
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
end module DGTimeIntegrator
