module QuadElementClass
    use Physics
    use SMConstants
    use NodeClass
    use NodesAndWeights_class
    use Storage_module
    use QuadMeshDefinitions
    implicit none

    private
    public  QuadElement_t , QuadElement_p , Edge_t , StraightBdryEdge_t , CurvedBdryEdge_t , Edge_p
!
!   ********************************************************************************
!           Quad Element derived type definition
!   ********************************************************************************
!
    type QuadElement_t
        integer                           :: ID
        integer                           :: address
!
!       -----------------------------------------------------------------
!              Connectivities
!       -----------------------------------------------------------------
!
        type(Node_p)                      :: nodes(POINTS_PER_QUAD)
        class(Edge_p), pointer            :: edges(:)
!
!       -----------------------------------------------------------------
!              Storage
!       -----------------------------------------------------------------
!
        class(NodesAndWeights_t), pointer :: spA
        class(NodesAndWeights_t), pointer :: spI
        real(kind=RP), allocatable        :: x(:,:,:)
        real(kind=RP), allocatable        :: jac(:,:)
        real(kind=RP), pointer            :: Q(:,:,:) , QDot(:,:,:) , F(:,:,:,:) , dQ(:,:,:,:)
        contains
            procedure      :: Construct  => QuadElement_Construct
            procedure      :: SetStorage => QuadElement_SetStorage
    end type QuadElement_t

    type QuadElement_p
        type(QuadElement_t),  pointer     :: e
    end type QuadElement_p
!
!   *********************************************************************************
!           Edge derived type definition:
!                 Edge_t:  For interior edges
!                 StraightBdryEdge_t: For straight boundary edges
!                 CurvedBdryEdge_t: For curved boundary edges
!   *********************************************************************************
!
    type Edge_t
        integer                           :: ID
        integer                           :: N
        real(kind=RP)                     :: nb(NDIM)          !   n always points from LEFT towards RIGHT and outside the domain for bdryedges
        type(Node_p)                      :: nodes(POINTS_PER_EDGE)
        class(QuadElement_p), pointer     :: quads(:)
        real(kind=RP), allocatable        :: Q(:,:,:) , dQ(:,:,:,:)  ! To store the interpolation to boundaries from elements
        integer,       allocatable        :: FaceLocation
        integer                           :: FaceType
        class(NodesAndWeights_t), pointer :: spA
        class(NodesAndWeights_t), pointer :: spI
    end type Edge_t

    type, extends(Edge_t)  :: StraightBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)           ! Solution at the boundary
        real(kind=RP), pointer            :: gB(:,:,:)         ! Solution gradient at the boundary
    end type StraightBdryEdge_t 

    type, extends(Edge_t)  :: CurvedBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)           ! Solution at the boundary
        real(kind=RP), pointer            :: gB(:,:,:)         ! Solution gradient at the boundary
    end type CurvedBdryEdge_t

    type Edge_p
        class(Edge_t),   pointer           :: f
        contains
            procedure      :: Construct => Edge_ConstructEdge
    end type Edge_p
!
!  -------------------------------------------------------------------------------------------------------------------------
!
!   ========
    contains
!   ========
!
        subroutine QuadElement_Construct(self , ID , nodes , N  , spA , address , storage , spI)
!            -----------------------------------------------------------------------------------
!                 This function performs the following operations to "self" element:
!                       -> Assign its ID and storage address position
!                       -> Assign its nodal storage spA, and spI if proceeds
!                       -> Assign its nodes
!                       -> Assign its storage
!            -----------------------------------------------------------------------------------
             use Setup_class
             use Physics
             class(QuadElement_t)              :: self
             integer                           :: ID
             class(Node_p)                     :: nodes(:)
             integer                           :: N
             class(NodalStorage)               :: spA
             integer                           :: address
             class(Storage_t)                  :: storage
             class(NodesAndWeights_t), pointer :: spI
!
             self % ID = ID
!
!            ************************
!            Point to nodes
!            ************************
!
             self % nodes = nodes
!
!            *********************
!            Get NodalStorage data       
!            *********************
!
             call spA % add( N , Setup % nodes , self % spA )
             self % spI => spI                              ! Points to NULL if Standard DG is chosen
!
!            ************
!            Linking data
!            ************
!            ---------------------------------------
!               These two arrays were before 
!               allocated inside each element.
!               Now they are just linked.
!                   allocate ( self % Q    ( 0:N , NEC )  ) 
!                   allocate ( self % QDot ( 0:N , NEC )  ) 
!                   allocate ( self % F    ( 0:N , NEC )  ) 
!                   allocate ( self % dQ   ( 0:N , NEC )  ) 
!            ---------------------------------------
!
             self % address = address
             call self % SetStorage( storage )
!
!            *************
!            Allocate data
!            *************
!
             allocate ( self % x   ( NDIM  , 0:N , 0:N  )  ) 
             allocate ( self % jac (         0:N , 0:N  )  ) 
             
        end subroutine QuadElement_Construct

        subroutine QuadElement_SetStorage( self , storage )
            use SMConstants
            use Storage_module
            use Setup_class
            use Physics
            implicit none
            class(QuadElement_t)      :: self
            class(Storage_t)        :: storage
        

            associate ( N => self % spA % N )
             self % Q    ( 0:N , 0:N , 1:NEC          ) => storage % Q    ( self % address: ) 
             self % QDot ( 0:N , 0:N , 1:NEC          ) => storage % QDot ( self % address: ) 
             self % dQ   ( 0:N , 0:N , 1:NEC , 1:NDIM ) => storage % dQ   ( self % address: ) 

             if ( trim(Setup % inviscid_discretization) .eq. "Over-Integration" ) then
               self % F (0: self % spI % N , 0: self % spI % N , 1:NEC , 1:NDIM) => storage % F ( (self % ID -1)*(self % spI % N+1)**2*NEC*NDIM + 1: self % ID * ( self % spI % N + 1)**2*NEC*NDIM )

             else

               self % F    ( 0:N , 0:N , 1:NEC , 1:NDIM )  => storage % F ( self % address: ) 

             end if

            end associate

        end subroutine QuadElement_SetStorage
            
        subroutine Edge_ConstructEdge( self , ID , curvilinear , faceType , spA , spI)
            use Setup_Class
            implicit none
            class(Edge_p)                     :: self
            integer                           :: ID
            logical                           :: curvilinear
            integer                           :: faceType
            class(NodalStorage)               :: spA
            class(NodesAndWeights_t), pointer :: spI
!
!           *************************************************
!              Allocate the edge depending on its type
!           *************************************************
!
!           Interior edges
!           --------------
            if (faceType .EQ. FACE_INTERIOR) then

                allocate(Edge_t :: self % f)

!               Allocate its elements 
                allocate( self % f % quads(QUADS_PER_EDGE) )

                self % f % faceType = FACE_INTERIOR

!
!           Boundary edges
!           --------------
            elseif (faceType .NE. FACE_INTERIOR) then
!
!              * Straight boundary edges
!                -----------------------
               if ( .NOT. curvilinear ) then
                  allocate(StraightBdryEdge_t   :: self % f)

!
!              * Curvilinear boundary edges
!                --------------------------
               else
                  allocate(CurvedBdryEdge_t     :: self % f)

               end if

               allocate( self % f % quads(1) )

               self % f % faceType = faceType

            end if

            self % f % ID = ID
            self % f % N  = Setup % N

            call spA % add( self % f % N , Setup % nodes , self % f % spA )
            self % f % spI    => spI


        end subroutine Edge_ConstructEdge 


end module QuadElementClass
