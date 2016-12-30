module Element1DClass
    use SMConstants
    use NodeClass
    use NodesAndWeights_class
    use Storage_module
    use QuadMeshDefinitions
    implicit none

    private
    public  QuadElement_t , QuadElement_p , Edge_t , StraightBdryEdge_t , CurvedBdryEdge_t
    public  ConstructQuadsAndEdges
!
!   ********************************************************************************
!           Quad Element derived type definition
!   ********************************************************************************
!
    type QuadElement_t
        type(Node_p)                      :: nodes(POINTS_PER_QUAD)
        integer                           :: ID
        integer                           :: edgesID(EDGES_PER_QUAD)
        integer                           :: address
        !       Storage
        class(NodesAndWeights_t), pointer :: Interp
        class(NodesAndWeights_t), pointer :: spI
        real(kind=RP), allocatable        :: x(:)
        real(kind=RP), pointer            :: Q(:,:,:) , QDot(:,:,:) , F(:,:,:,:) , dQ(:,:,:,:)
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
        class(QuadElement_p), pointer     :: elements(:)
        real(kind=RP)                     :: n(NDIM)          !   n always points from LEFT towards RIGHT and outside the domain for bdryedges
        type(Node_t), pointer             :: node             
        real(kind=RP), allocatable        :: Q(:,:,:) , dQ(:,:,:,:)  ! To store the interpolation to boundaries from elements
        integer,       allocatable        :: FaceLocation
        integer                           :: FaceType
    end type Edge_t

    type, extends(Edge_t)  :: StraightBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)
        real(kind=RP), pointer            :: gB(:,:,:)
    end type StraightBdryEdge_t 

    type, extends(Edge_t)  :: CurvedBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)
        real(kind=RP), pointer            :: gB(:,:,:)
    end type CurvedBdryEdge_t

    type Edge_p
        class(Edge_t),   pointer           :: f
    end type Edge_p
!
!  -------------------------------------------------------------------------------------------------------------------------


    contains
        subroutine ConstructQuadsAndEdges


        end subroutine ConstructQuadsAndEdges

        subroutine ConstructElement(self , ID , nodes , faceLeftID , faceRightID , N , nodes , spA , address , storage , spI)
             use Setup_class
             use Physics
             class(QuadElement_t)              :: self
             integer                           :: ID
             class(Node_p)                     :: nodes
             integer                           :: faceLeftID
             integer                           :: faceRightID
             integer                           :: nodes
             integer                           :: N
             class(NodalStorage)               :: spA
             integer                           :: address
             class(Storage_t)                  :: storage
             class(NodesAndWeights_t), pointer :: spI
!
             self % ID = ID
!
!            ************************
!            Point to neighbour nodes
!            ************************
!
             self % nodes = nodes
!
!            **************
!            Get faces data
!            **************
!
             self % edgesID(LEFT) = faceLeftID
             self % edgesID(RIGHT) = faceRightID
             self % edgesID(TOP)   = fa
!
!            *********************
!            Get NodalStorage data       
!            *********************
!
             call spA % add( N , nodes , self % Interp )
             self % spI => spI
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
             allocate ( self % x   ( NDIM  , 0:N                   )  ) 
             allocate ( self % Qb  ( NEC   , EDGES_PER_QUAD        )  ) 
             allocate ( self % dQb ( NEC   , EDGES_PER_QUAD , NDIM )  ) 
!
!            **********
!            Set values
!            **********
!
             self % x = 0.5_RP * (self % nodes(LEFT) % n % x + self % nodes(RIGHT) % n % x + self % Interp % xi * (self % nodes(RIGHT) % n % x - self % nodes(LEFT) % n % x ) )
            ! TODO
             self % hdiv2 = 0.0_RP !0.5_RP*abs(self % nodes(LEFT) % n % x - self % nodes(RIGHT) % n % x)
             
        end subroutine ConstructElement

        subroutine QuadElement_SetStorage( self , storage )
            use SMConstants
            use Storage_module
            use Setup_class
            use Physics
            implicit none
            class(QuadElement_t)      :: self
            class(Storage_t)        :: storage
        

            associate ( N => self % Interp % N )
             self % Q    ( 0:N , 1:NEC )  => storage % Q    ( self % address: ) 
             self % QDot ( 0:N , 1:NEC )  => storage % QDot ( self % address: ) 

             if ( trim(Setup % inviscid_discretization) .eq. "Over-Integration" ) then
               self % F (0: self % spI % N , 1:NEC) => storage % F ( (self % ID -1)*(self % spI % N+1) + 1: self % ID * ( self % spI % N + 1) )

             else

               self % F    ( 0:N , 1:NEC )  => storage % F    ( self % address: ) 

             end if

             self % dQ   ( 0:N , 1:NEC )  => storage % dQ   ( self % address: ) 

            end associate

        end subroutine QuadElement_SetStorage
            



end module Element1DClass
