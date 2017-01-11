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
   implicit none

   private
   public BoundaryCondition_t , Construct!, DGBoundaryConditions_setFace

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
         class(BoundaryCondition_t)                :: self
         integer                                   :: marker


      end subroutine BoundaryConditions_construct
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
