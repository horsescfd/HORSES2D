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
   public BoundaryCondition_t , Construct!, DGBoundaryConditions_setFace
   public PeriodicBC_t , DirichletBC_t , Farfield_t , EulerWall_t

   integer, parameter         :: STR_LEN_BC = 128

!   type BoundaryCondition_t
!      integer        :: marker
!      integer        :: type
!      class(BoundaryCondition_t), pointer    :: periodicPair => NULL()
!      class(Edge_t), pointer                 :: edge => NULL()
!      real(kind=RP), pointer                 :: uBC => NULL()
!      real(kind=RP), pointer                 :: gBC => NULL()
!      contains
!         procedure   :: setFace => DGBoundaryConditions_setFace
!         procedure   :: construct => DGBoundaryConditions_construct
!   end type BoundaryCondition_t
!
!
!  **********************************
!  Base class for boundary conditions
!  **********************************
!
   type BoundaryCondition_t
      character(len=STR_LEN_BC)        :: Name
      integer                          :: BCType
      contains
         procedure ::     Construct => BaseClass_Construct
   end type BoundaryCondition_t
!
!  *********************************
!  Periodic boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: PeriodicBC_t
   end type PeriodicBC_t
!
!  **********************************
!  Dirichlet boundary condition class
!  **********************************
!
   type, extends(BoundaryCondition_t)           :: DirichletBC_t
      real(kind=RP), dimension(NEC)       :: q
      contains
         procedure ::      Construct => DirichletBC_Construct
   end type DirichletBC_t
!
!  *********************************
!  Farfield boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: Farfield_t
      real(kind=RP)           :: u
      real(kind=RP)           :: v
      real(kind=RP)           :: p
      real(kind=RP)           :: T
   end type FarField_t
!
!  ***********************************
!  Euler wall boundary condition class
!  ***********************************
!
   type, extends(BoundaryCondition_t)           :: EulerWall_t
   end type EulerWall_t

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
      
         elseif ( BCType .eq. "EulerWall") then
            allocate( EulerWall_t      :: self )
            self % BCType = EULERWALL_BC

         end if

         call self % Construct( marker , in_label)

      end subroutine BoundaryConditions_construct
   
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

      subroutine DirichletBC_Construct( self , marker , in_label)
         use Setup_class
         implicit none
         class(DirichletBC_t)      :: self
         integer                   :: marker
         character(len=*)          :: in_label
         real(kind=RP), allocatable             :: pressure
         real(kind=RP), allocatable         :: Temperature
         real(kind=RP), allocatable         :: Mach
         real(kind=RP), allocatable         :: AngleOfAttack
         real(kind=RP)                          :: rho

         call readValueInRegion( trim(Setup % bdry_file) , "Name" , self % Name , in_label , "# end" )
         call readValueInRegion( trim(Setup % bdry_file) , "pressure" , pressure , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Temperature", Temperature , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Mach" , Mach , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Angle of attack" , AngleOfAttack , in_label , "# end")
         
         if ( allocated(pressure) ) then
            pressure = pressure / refValues % p
         else
            allocate(pressure)
            pressure = 1.0_RP
         end if

         if ( allocated(Temperature) ) then
            Temperature = Temperature / refValues % T
         else
            allocate(Temperature)
            Temperature = 1.0_RP
         end if

         if ( .not. allocated(Mach) ) then
            allocate(Mach)
            Mach = Dimensionless % Mach
         end if

         if ( allocated(AngleOfAttack) ) then
            AngleOfAttack = AngleOfAttack * PI / 180.0_RP
         else
            allocate(AngleOfAttack)
            AngleOfAttack = 0.0_RP
         end if
!
!        Construct the state vector
!        --------------------------
         associate ( gamma => Thermodynamics % Gamma , cv => Dimensionless % cv)
         rho = pressure / Temperature
         self % q(IRHO) = rho
         self % q(IRHOU) = rho * sqrt(gamma) * Mach * cos(AngleOfAttack)
         self % q(IRHOV) = rho * sqrt(gamma) * Mach * sin(AngleOfAttack)
         self % q(IRHOE) = cv * pressure + 0.5_RP * rho * gamma * Mach * Mach
         end associate

      end subroutine DirichletBC_Construct
!
!      subroutine DGBoundaryConditions_setFace( self  )
!         use Physics
!         use Setup_class
!         implicit none
!         class(BoundaryCondition_t)       :: self
!
!         select type (f1=>self % edge)
!            type is (StraightBdryEdge_t)
!               if (self % type  .eq. PERIODIC_BC) then
!!                 -------------------------------
!!                    Look for its pairing
!!                 -------------------------------
!                  select type (f2=>self % periodicPair % edge)
!                     type is (StraightBdryEdge_t)
!                        f1 % uB => NULL()    ! TODO f2 % quads(1) % e % Qb( : , f2 % BCLocation )  
!                        f1 % gB => NULL()    ! TODO f2 % quads(1) % e % dQb( : , f2 % BCLocation )  
!                  end select
!               elseif ( self % type .eq. DIRICHLET_BC) then
!!                 -------------------------------
!!                    Allocate data and set values
!!                 -------------------------------
!                        allocate( f1 % uB (NEC , 0 : f1 % spA % N) )
!                        f1 % uB = Setup % dirichletBC( f1 % edgeType )
!                        f1 % gB => NULL()    ! TODO f1 % quads(1) % e % dQb( : , f1 % BCLocation )
!
!               end if
!         end select
!
!         
!
!
!      end subroutine DGBoundaryConditions_setFace   
!
!
end module DGBoundaryConditions
