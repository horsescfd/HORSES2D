module DGBoundaryConditions
   use SMConstants
   use QuadElementClass
   implicit none

   private
   public BoundaryCondition_t , DGBoundaryConditions_setFace

   integer, parameter         :: STR_LEN_BC = 128

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
