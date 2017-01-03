module DGBoundaryConditions
   use SMConstants
   use QuadElementClass
   use QuadMeshClass
   implicit none

   private
   public Zone_t , BoundaryCondition_t , DGBoundaryConditions_setFace

   integer, parameter         :: STR_LEN_BC = 128

   type Zone_t
      integer                   :: marker
      character(len=STR_LEN_BC) :: Name
      integer                   :: no_of_edges
      class(Edge_p), pointer    :: edges(:)
      contains
         procedure      :: Construct => Zone_Construct
   end type Zone_t

   type BoundaryCondition_t
      integer        :: marker
      integer        :: type
      class(BoundaryCondition_t), pointer    :: periodicPair => NULL()
      class(Edge_t), pointer                 :: edge => NULL()
      real(kind=RP), pointer                 :: uBC => NULL()
      real(kind=RP), pointer                 :: gBC => NULL()
      contains
         procedure   :: setFace => DGBoundaryConditions_setFace
         procedure   :: construct => DGBoundaryConditions_construct
   end type BoundaryCondition_t

   contains
      subroutine Zone_construct( self , mesh , marker , name)
         implicit none
         class(Zone_t)           :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: marker
         character(len=*)        :: name
         integer                 :: edID
         integer                 :: current

         self % marker = marker
         self % Name = trim(Name)

         self % no_of_edges = 0
!
!        ***************************************
!        Gather the number of edges for a marker
!        ***************************************
!
         do edID = 1 , mesh % no_of_edges
            if ( mesh % edges(edID) % f % edgeType .eq. marker) then
               self % no_of_edges = self % no_of_edges + 1
            end if
         end do
!
!        Allocate the structure
         allocate( self % edges( self % no_of_edges ) )

!
!        Point to all edges in the zone
         current = 0
         do edID = 1 , mesh % no_of_edges
            if ( mesh % edges(edID) % f % edgeType .eq. marker) then
               current = current + 1
               self % edges( current ) % f => mesh % edges(edID) % f
            end if
         end do

      end subroutine Zone_construct

      subroutine DGBoundaryConditions_construct( self , ID , edge , BCset)
         use Setup_class
         implicit none
         class(BoundaryCondition_t)          :: self
         integer                             :: ID
         class(Edge_t), pointer              :: edge
         class(BoundaryCondition_t), pointer :: BCset(:)

         self % marker = ID
         self % type   = Setup % BCTypes(ID)
         self % edge   => edge

         if ( self % type .eq. PERIODIC_BC) then

            self % periodicPair => BCset( Setup % periodicBCFaces(ID) )

         elseif (self % type .eq. DIRICHLET_BC) then
         
            allocate( self % uBC )

            self % uBC = Setup % DirichletBC(ID)
         
         end if
                  

      end subroutine DGBoundaryConditions_construct

      subroutine DGBoundaryConditions_setFace( self  )
         use Physics
         use Setup_class
         implicit none
         class(BoundaryCondition_t)       :: self

         select type (f1=>self % edge)
            type is (StraightBdryEdge_t)
               if (self % type  .eq. PERIODIC_BC) then
!                 -------------------------------
!                    Look for its pairing
!                 -------------------------------
                  select type (f2=>self % periodicPair % edge)
                     type is (StraightBdryEdge_t)
                        f1 % uB => NULL()    ! TODO f2 % quads(1) % e % Qb( : , f2 % BCLocation )  
                        f1 % gB => NULL()    ! TODO f2 % quads(1) % e % dQb( : , f2 % BCLocation )  
                  end select
               elseif ( self % type .eq. DIRICHLET_BC) then
!                 -------------------------------
!                    Allocate data and set values
!                 -------------------------------
                        allocate( f1 % uB (NEC , 0 : f1 % spA % N) )
                        f1 % uB = Setup % dirichletBC( f1 % edgeType )
                        f1 % gB => NULL()    ! TODO f1 % quads(1) % e % dQb( : , f1 % BCLocation )

               end if
         end select

         


      end subroutine DGBoundaryConditions_setFace   


end module DGBoundaryConditions
