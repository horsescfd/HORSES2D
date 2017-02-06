!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
module DGBoundaryConditions
   use SMConstants
   use Physics
   use QuadElementClass
   use ParamfileIO
   implicit none

   private
   public BoundaryCondition_t , Construct
   public PeriodicBC_t , DirichletBC_t , FarfieldBC_t , EulerWall_t , PressureOutletBC_t , PressureInletBC_t
   public RiemannBC_t

   integer, parameter         :: STR_LEN_BC             = 128
   integer, parameter         :: SPECIFY_SPEED          = 1
   integer, parameter         :: SPECIFY_TOTAL_PRESSURE = 2
!
!  **********************************
!  Base class for boundary conditions
!  **********************************
!
   type BoundaryCondition_t
      character(len=STR_LEN_BC)                         :: Name
      integer                                           :: marker
      integer                                           :: BCType
      character(len=STR_LEN_BC)                         :: RiemannSolverName
      procedure(RiemannSolverFunction), pointer, nopass :: RiemannSolver => NULL()
      contains
         procedure :: Construct        => BaseClass_Construct
         procedure :: SetRiemannSolver => BoundaryConditions_SetRiemannSolver
         procedure :: Associate        => BaseClass_Associate
         procedure :: Update           => BaseClass_Update
         procedure :: Describe         => BaseClass_Describe
   end type BoundaryCondition_t
!
!  *********************************
!  Periodic boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: PeriodicBC_t
      logical                          :: associated
      integer                          :: direction
      integer                          :: connected_marker
      contains
         procedure ::     Construct => PeriodicBC_Construct
         procedure ::     Describe  => PeriodicBC_Describe
   end type PeriodicBC_t
!
!  **********************************
!  Dirichlet boundary condition class
!  **********************************
!
   type, extends(BoundaryCondition_t)           :: DirichletBC_t
      real(kind=RP), dimension(NCONS)       :: q
      contains
         procedure ::      Construct => DirichletBC_Construct
         procedure ::      Associate => DirichletBC_Associate
         procedure ::      Describe  => DirichletBC_Describe
   end type DirichletBC_t
!
!  *********************************
!  Farfield boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: FarfieldBC_t
      real(kind=RP), dimension(NCONS)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      contains
         procedure ::      Construct => FarfieldBC_Construct
         procedure ::      Associate => FarfieldBC_Associate
         procedure ::      Update    => FarfieldBC_Update
         procedure ::      Describe  => FarfieldBC_Describe
   end type FarFieldBC_t
!
!  ********************************
!  Outflow boundary condition class
!  ********************************
!
   type, extends(BoundaryCondition_t)           :: PressureOutletBC_t
      real(kind=RP), dimension(NCONS)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      contains
         procedure ::      Construct => PressureOutletBC_Construct
         procedure ::      Associate => PressureOutletBC_Associate
         procedure ::      Update    => PressureOutletBC_Update
         procedure ::      Describe  => PressureOutletBC_Describe
   end type PressureOutletBC_t
!
!  ***************************************
!  Inflow/Outflow boundary condition class
!  ***************************************
!
   type, extends(BoundaryCondition_t)           :: PressureInletBC_t
      real(kind=RP), dimension(NCONS)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      real(kind=RP)                       :: rhot
      real(kind=RP)                       :: st
      real(kind=RP)                       :: Ht
      contains
         procedure ::      Construct => PressureInletBC_Construct
         procedure ::      Associate => PressureInletBC_Associate
         procedure ::      Update    => PressureInletBC_Update
         procedure ::      Describe  => PressureInletBC_Describe
   end type PressureInletBC_t
!
!  ********************************
!  Riemann boundary condition class
!  ********************************
!
   type, extends(BoundaryCondition_t)           :: RiemannBC_t
      real(kind=RP), dimension(NCONS)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      real(kind=RP)                       :: Ht
      real(kind=RP)                       :: Rminus
      integer                             :: mode
      contains
         procedure ::      Construct => RiemannBC_Construct
         procedure ::      Associate => RiemannBC_Associate
         procedure ::      Update    => RiemannBC_Update
         procedure ::      Describe  => RiemannBC_Describe
   end type RiemannBC_t

!
!  ***********************************
!  Euler wall boundary condition class
!  ***********************************
!
   type, extends(BoundaryCondition_t)           :: EulerWall_t
      contains
         procedure   ::    Associate => EulerWall_Associate
         procedure   ::    Update    => EulerWall_Update
         procedure   ::    Describe  => EulerWall_Describe
   end type EulerWall_t
!
!  *************************************
!  Viscous wall boundary condition class
!  *************************************
!
   type, extends(BoundaryCondition_t)           :: ViscousWall_t
      contains
         procedure   ::    Associate => ViscousWall_Associate
         procedure   ::    Update    => ViscousWall_Update
         procedure   ::    Describe  => ViscousWall_Describe
   end type ViscousWall_t


!
!  *******************
!  Construct procedure
!  *******************
!
   interface Construct
      module procedure BoundaryConditions_construct
   end interface Construct
!
!///////////////////////////////////////////////////////////////////////////////////
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////
!
      include './BoundaryConditions/DirichletBC.incf'
      include './BoundaryConditions/EulerWallBC.incf'
      include './BoundaryConditions/FarfieldBC.incf'
      include './BoundaryConditions/PeriodicBC.incf'
      include './BoundaryConditions/PressureInletBC.incf'
      include './BoundaryConditions/PressureOutletBC.incf'
      include './BoundaryConditions/RiemannBC.incf'
      include './BoundaryConditions/ViscousWall.incf'


      

      subroutine BoundaryConditions_construct( self , marker)
         use Setup_class
         implicit none
         class(BoundaryCondition_t), pointer :: self
         integer                             :: marker
         character(len=STR_LEN_BC)           :: BCType
         character(len=STR_LEN_BC)           :: in_label


         write(in_label,'(A,I0)') "# define zone ",marker


         call ReadValueInRegion( trim(Setup % bdry_file) , "Type" , BCType , in_label , "# end" )

         if (BCType .eq. "Dirichlet") then
            allocate( DirichletBC_t    :: self )
            self % BCType = DIRICHLET_BC
   
         elseif (BCType .eq. "Farfield") then
            allocate( FarfieldBC_t     :: self )
            self % BCType = FARFIELD_BC

         elseif (BCType .eq. "Pressure outlet") then
            allocate( PressureOutletBC_t     :: self )
            self % BCType = PRESSUREOUTLET_BC
      
         elseif (BCType .eq. "Pressure inlet") then
            allocate( PressureInletBC_t     :: self )
            self % BCType = PRESSUREINLET_BC

         elseif ( BCType .eq. "Riemann") then
            allocate( RiemannBC_t     :: self )
            self % BCType = RIEMANN_BC

         elseif ( BCType .eq. "Periodic") then
            allocate( PeriodicBC_t     :: self )
            self % BCType = PERIODIC_BC

         elseif ( BCType .eq. "Euler wall") then
            allocate( EulerWall_t      :: self )
            self % BCType = EULERWALL_BC

         elseif ( BCType .eq. "Viscous wall") then
            allocate( ViscousWall_t    :: self ) 
            self % BCType = VISCOUSWALL_BC

         else
            print*, 'Boundary condition "',trim(BCType),'" not implemented yet.'
            print*, "Options available are:"
            print*, "   * Dirichlet"
            print*, "   * Farfield"
            print*, "   * Pressure outlet"
            print*, "   * Pressure inlet"
            print*, "   * Riemann"
            print*, "   * Periodic"
            print*, "   * Euler wall"
            print*, "   * Viscous wall"
            stop "Stopped."

         end if
   
         self % marker = marker

         call readValueInRegion( trim(Setup % bdry_file) , "Name" , self % Name , in_label , "# end" )

         call self % SetRiemannSolver()

         call self % Construct( marker , in_label)

      end subroutine BoundaryConditions_construct

      subroutine BoundaryConditions_SetRiemannSolver( self , RiemannSolver )
         use Setup_class
         implicit none
         class(BoundaryCondition_t)                :: self
         character(len=*), optional                :: RiemannSolver
         character(len=STR_LEN_BC)                 :: in_label
         
         if ( present ( RiemannSolver ) ) then

            self % RiemannSolverName = trim(RiemannSolver)

         else

            write(in_label,'(A,I0)') "# define zone ", self % marker
            call readValueInRegion( trim(Setup % bdry_file) , "Riemann solver" , self % RiemannSolverName ,  in_label , "# end" )

         end if

         if (trim(self % RiemannSolverName ) .eq. "Exact") then
            self % RiemannSolver => ExactRiemannSolver

         elseif ( trim(self % RiemannSolverName) .eq. "Roe" ) then
            self % RiemannSolver => RoeFlux

         elseif ( trim(self % RiemannSolverName) .eq. "HLL" ) then
            self % RiemannSolver => HLLFlux

         elseif ( trim(self % RiemannSolverName) .eq. "Interior" ) then
            self % RiemannSolver => NULL()

         else     ! Default: The interior Riemann solver
            self % RiemannSolver => NULL()
            self % RiemannSolverName = "Interior"

         end if

      end subroutine BoundaryConditions_SetRiemannSolver
   
      subroutine BaseClass_Construct( self , marker , in_label)
         implicit none
         class(BoundaryCondition_t)                :: self
         integer                                   :: marker
         character(len=*)                          :: in_label
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Construct

      subroutine BaseClass_Associate( self , edge )
         implicit none
         class(BoundaryCondition_t)          :: self
         class(Edge_t)                       :: edge
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Associate

      subroutine BaseClass_Update( self , edge )
         implicit none
         class(BoundaryCondition_t)          :: self
         class(Edge_t)                       :: edge
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Update
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              Initialization subroutines
!              --------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_WeakRiemann(self , edge)
         implicit none
         class(BoundaryCondition_t)       :: self
         class(Edge_t)                    :: edge
         integer                          :: i

         associate ( N => edge % spA % N )
         select type ( edge )
         
            type is (Edge_t)
               print*, "Only boundary edges are expected."
               stop "Stopped"
      
            type is (StraightBdryEdge_t)
               allocate( edge % uB(0:N,NCONS) )
               allocate( edge % gB(0:N,NCONS,NDIM) )   ! Normal gradients
               
               do i = 0 , N
                  edge % uB(i , 1:NCONS) = 0.0_RP      ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NCONS , 1:NDIM ) = 0.0_RP
               
   
            type is (CurvedBdryEdge_t)
               allocate( edge % uB(0:N,NCONS) )
               allocate( edge % gB(0:N,NCONS,NDIM) )

               do i = 0 , N
                  edge % uB(i , 1:NCONS) = 0.0_RP    ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NCONS , 1:NDIM ) = 0.0_RP

         end select

         end associate
      end subroutine Initialize_WeakRiemann

      subroutine Initialize_WeakPrescribed(self , edge)
         implicit none
         class(BoundaryCondition_t)       :: self
         class(Edge_t)                    :: edge
         integer                          :: i

         associate ( N => edge % spA % N )
         select type ( edge )
         
            type is (Edge_t)
               print*, "Only boundary edges are expected."
               stop "Stopped"
      
            type is (StraightBdryEdge_t)
               allocate( edge % FB(0:N,NCONS) )
               
               do i = 0 , N
                  edge % FB(i , 1:NCONS) = 0.0_RP      ! Its value is not given until the update routine is invoked
               end do
   
            type is (CurvedBdryEdge_t)
               allocate( edge % FB(0:N,NCONS) )

               do i = 0 , N
                  edge % FB(i , 1:NCONS) = 0.0_RP    ! Its value is not given until the update routine is invoked
               end do

         end select

         end associate
      end subroutine Initialize_WeakPrescribed

!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              Describe subroutines
!              --------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_Describe( self ) 
         implicit none
         class(BoundaryCondition_t)          :: self

      end subroutine BaseClass_Describe

      subroutine DirichletBC_Describe( self )
         use Physics
         implicit none
         class(DirichletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Dirichlet."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         end associate

      end subroutine DirichletBC_Describe

      subroutine FarfieldBC_Describe( self )
         use Physics
         implicit none
         class(FarfieldBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Farfield."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine FarfieldBC_Describe

      subroutine PressureOutletBC_Describe( self )
         use Physics
         implicit none
         class(PressureOutletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Pressure outlet."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine PressureOutletBC_Describe

      subroutine PressureInletBC_Describe( self )
         use Physics
         implicit none
         class(PressureInletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Pressure inlet."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt * refValues % p 
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total temperature: " , self % Tt * refValues % T
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total enthalpy: " , self % Ht * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine PressureInletBC_Describe

      subroutine RiemannBC_Describe( self )
         use Physics
         implicit none
         class(RiemannBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Pressure inlet."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt * refValues % p 
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total temperature: " , self % Tt * refValues % T
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total enthalpy: " , self % Ht * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine RiemannBC_Describe

      subroutine PeriodicBC_Describe( self ) 
         implicit none
         class(PeriodicBC_t)          :: self

         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Periodic."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X,A,A25,I0,A)') "-> ", "Connecting zone: ",self % connected_marker
         if ( self % direction .eq. IX) then
            write(STD_OUT , '(30X,A,A25,A)') "-> " , "Direction: " ,"x"
         elseif ( self % direction .eq. IY) then
            write(STD_OUT , '(30X,A,A25,A)') "-> " , "Direction: " ,"y"
         end if
         
      end subroutine PeriodicBC_Describe

      subroutine EulerWall_Describe( self )
         implicit none
         class(EulerWall_t)                :: self

         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Euler wall."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)

      end subroutine EulerWall_Describe

      subroutine ViscousWall_Describe( self )
         implicit none
         class(ViscousWall_t)                :: self

         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Viscous wall."
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Riemann solver: " , trim(self % RiemannSolverName)

      end subroutine ViscousWall_Describe

end module DGBoundaryConditions
