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

!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!   ********************************************************************************
!           Quad Element derived type definition
!   ********************************************************************************
!
    type QuadElement_t
        integer                            :: ID                             ! ID of the element
        integer                            :: address                        ! Memory address of the first position in the mesh
        integer                            :: edgesDirection(EDGES_PER_QUAD) ! Direction (FORWARD/REVERSE) of the edges
        integer                            :: quadPosition(EDGES_PER_QUAD)   ! Position of the quad for the edge (LEFT/RIGHT)
        real(kind=RP), allocatable         :: x(:,:,:)                       ! Coordinates of the nodes ( X/Y , xi , eta )
        real(kind=RP), allocatable         :: dx(:,:,:,:)                    ! Mapping derivatives (X/Y , xi , eta , dxi / deta)
        real(kind=RP), allocatable         :: jac(:,:)                       ! Mapping jacobian ( xi , eta )
        real(kind=RP), pointer             :: Q(:,:,:)                       ! Pointers to the main storage:
        real(kind=RP), pointer             :: QDot(:,:,:)                    ! *   Q, QDot ( xi , eta , eq ): solution and time derivative
        real(kind=RP), pointer             :: F(:,:,:,:)                     ! *   F ( xi , eta , eq , X/Y) : contravariant fluxes
        real(kind=RP), pointer             :: dQ(:,:,:,:)                    ! *   dQ( xi ,eta , eq , X/Y):   solution gradient
        type(Node_p)                       :: nodes(POINTS_PER_QUAD)         ! Pointer to neighbour nodes
        class(Edge_p), pointer             :: edges(:)                       ! Pointer to neighbour eges
        class(NodesAndWeights_t), pointer  :: spA                            ! Pointer to the Nodal Storage
        class(NodesAndWeights_t), pointer  :: spI                            ! Pointer to the Over-Integration Nodal Storage (If proceeds)
!       ========
        contains
!       ========
            procedure      :: Construct   => QuadElement_Construct                                 ! Constructs/allocates data and points to the address in Storage
            procedure      :: SetStorage  => QuadElement_SetStorage                                ! Function to set the storage distribution
            procedure      :: SetMappings => QuadElement_SetMappings                               ! Function to compute the mapping data (x, dx, jac)
            procedure      :: Ja          => QuadElement_MetricMatrix
    end type QuadElement_t

!
!   ********************************************************************************
!           Pointer to Quad Element derived type
!   ********************************************************************************
!
    type QuadElement_p
        type(QuadElement_t),  pointer     :: e                                                     ! Pointer to a quad element
    end type QuadElement_p

!//////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!   *********************************************************************************
!           Edge derived type definition:
!                 Edge_t:  For interior edges
!                 StraightBdryEdge_t: For straight boundary edges
!                 CurvedBdryEdge_t: For curved boundary edges
!   *********************************************************************************
!
    type Edge_t
        integer                             :: ID                         ! Edge ID
        integer                             :: edgeType                   ! Edge Type: FACE_INTERIOR , or the marker value if boundary face
        integer,                    pointer :: edgeLocation(:)            ! Edge location for the two (or one) sharing elements (ETOP,EBOTTOM,ELEFT,ERIGHT)
        real(kind=RP),              pointer :: n(:,:)                     ! Unitary normal: points from LEFT towards RIGHT, and outside the domain for bdryedges
        real(kind=RP),              pointer :: X(:,:)                     ! Coordinates: (X/Y, xi)
        real(kind=RP),              pointer :: dX(:,:)                    ! Tangent vector: (X/Y, xi)
        real(kind=RP),              pointer :: dS(:,:)                    ! Surface differential vector (X/Y, xi)
        real(kind=RP),              pointer :: Q(:,:,:)                   ! Solution interpolation to boundaries ( xi , eq , LEFT/RIGHT )
        real(kind=RP),              pointer :: dQ(:,:,:,:)                ! Solution gradient interpolation to boundary ( xi , eq ,  X/Y , LEFT/RIGHT)
        type(Node_p)                        :: nodes(POINTS_PER_EDGE)     ! Pointer to the two nodes
        class(QuadElement_p),       pointer :: quads(:)                   ! Pointers to the two (or one) shared quads
        class(NodesAndWeights_t),   pointer :: spA                        ! Pointer to the approximation nodal storage
        class(NodesAndWeights_t),   pointer :: spI                        ! Pointer to the integration nodal storage (if over-integration is active)
        contains
            procedure      :: SetCurve    => Edge_SetCurve                    ! Procedure that computes the coordinates, the tangent, and the normal.
            procedure      :: Invert      => Edge_Invert                      ! Function to invert the edge orientation 
            procedure      :: evaluateX   => Edge_AnalyticalX                 ! Function to compute an edge point in a local coordinate "xi"
            procedure      :: evaluatedX  => Edge_AnalyticaldX                 ! Function to compute an edge point in a local coordinate "xi"
            procedure      :: evaluatedS  => Edge_AnalyticaldS
            procedure      :: getX        => Edge_getX                        ! 
            procedure      :: getdX       => Edge_getdX
            procedure      :: getdS       => Edge_getdS
    end type Edge_t

    type, extends(Edge_t)  :: StraightBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)           ! Solution at the boundary
        real(kind=RP), pointer            :: gB(:,:,:)         ! Solution gradient at the boundary
    end type StraightBdryEdge_t 

    type, extends(Edge_t)  :: CurvedBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)           ! Solution at the boundary
        real(kind=RP), pointer            :: gB(:,:,:)         ! Solution gradient at the boundary
        contains
            procedure      :: SetCurve   => CurvilinearEdge_SetCurve       ! Procedure that computes the coordinates, the tangent, and the normal
            procedure      :: evaluateX  => Curvilinear_InterpolantX
            procedure      :: evaluatedX => Curvilinear_InterpolantdX
            procedure      :: evaluatedS => Curvilinear_InterpolantdS
    end type CurvedBdryEdge_t

!
!   ********************************************************************************
!           Pointer to Edge derived type
!   ********************************************************************************
!
    type Edge_p
        class(Edge_t),   pointer           :: f
        contains
            procedure      :: Construct => Edge_ConstructEdge
            procedure      :: LinkWithElements => Edge_linkWithElements
    end type Edge_p
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------------------------------------------
!
!   ========
    contains
!   ========
!

        include 'QuadMappings.incf'
        include 'EdgeMappings.incf'
        include 'CurvedEdgeMappings.incf'
        include 'QuadAuxiliar.incf'


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
!            --------------------------------------------------------
             integer                           :: edge
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
!                   allocate ( self % Q    ( 0:N , 0:N , NEC )  ) 
!                   allocate ( self % QDot ( 0:N , 0:N , NEC )  ) 
!                   allocate ( self % F    ( 0:N , 0:N , NEC )  ) 
!                   allocate ( self % dQ   ( 0:N , 0:N , NEC )  ) 
!            ---------------------------------------
!
             self % address = address
             call self % SetStorage( storage )
!
!            *************
!            Allocate data
!            *************
!
             allocate ( self % x     ( NDIM  , 0:N , 0:N        )  ) 
             allocate ( self % dx    ( NDIM  , 0:N , 0:N , NDIM )  ) 
             allocate ( self % jac   ( 0:N , 0:N                )  ) 
             allocate ( self % edges ( EDGES_PER_QUAD           )  ) 
      
             do edge = 1 , EDGES_PER_QUAD
               self % edges(edge) % f => NULL()
             end do
             
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
             self % dQ   ( 0:N , 0:N , 1:NEC , 1:NDIM ) => storage % dQ   ( (self % address-1)*NDIM + 1: ) 

             if ( trim(Setup % inviscid_discretization) .eq. "Over-Integration" ) then
               self % F (0: self % spI % N , 0: self % spI % N , 1:NEC , 1:NDIM) => &
                           storage % F ( (self % ID -1)*(self % spI % N+1)**2*NEC*NDIM + 1: self % ID * ( self % spI % N + 1)**2*NEC*NDIM )

             else

               self % F    ( 0:N , 0:N , 1:NEC , 1:NDIM )  => storage % F ( (self % address-1)*NDIM+1: ) 

             end if

            end associate

        end subroutine QuadElement_SetStorage
            
        subroutine Edge_ConstructEdge( self , ID , curvilinear , nodes , edgeType , spA , spI)
            use Setup_Class
            implicit none
            class(Edge_p)                     :: self
            integer                           :: ID
            class(Node_p)                     :: nodes(:)
            logical                           :: curvilinear
            integer                           :: edgeType
            integer                           :: node
            class(NodalStorage)               :: spA
            class(NodesAndWeights_t), pointer :: spI
!           --------------------------------------------------
            integer                           :: quad
            integer                           :: p
!
!           *************************************************
!              Allocate the edge depending on its type
!           *************************************************
!
!           Interior edges
!           --------------
            if (edgeType .EQ. FACE_INTERIOR) then

                allocate(Edge_t :: self % f)

!               Allocate its elements 
                allocate ( self % f % quads        ( QUADS_PER_EDGE )  ) 
                allocate ( self % f % edgeLocation ( QUADS_PER_EDGE )  ) 

                do quad = 1 , QUADS_PER_EDGE
                  self % f % quads(quad) % e => NULL()
                end do

                self % f % edgeType = FACE_INTERIOR

!
!           Boundary edges
!           --------------
            elseif (edgeType .NE. FACE_INTERIOR) then
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

               allocate ( self % f % quads         ( 1 )  ) 
               allocate ( self % f % edgeLocation ( 1 )  ) 

               self % f % quads(1) % e => NULL()
               self % f % edgeType = edgeType


            end if

            self % f % ID = ID

            do node = 1 , POINTS_PER_EDGE
               self % f % nodes(node) % n => nodes(node) % n
            end do

            call spA % add( Setup % N , Setup % nodes , self % f % spA )
            self % f % spI    => spI

!           Geometry
!           --------
            allocate ( self % f % X  ( NDIM , 0 : self % f % spA % N )  ) 
            allocate ( self % f % dX ( NDIM , 0 : self % f % spA % N )  ) 
            allocate ( self % f % dS ( NDIM , 0 : self % f % spA % N )  ) 
            allocate ( self % f % n  ( NDIM , 0 : self % f % spA % N )  ) 

            if (edgeType .eq. FACE_INTERIOR) then
         
               allocate ( self % f % Q  ( 0 : self % f % spA % N , NEC , POINTS_PER_QUAD        )  ) 
               allocate ( self % f % dQ ( 0 : self % f % spA % N , NEC , NDIM , POINTS_PER_QUAD )  ) 

            else
   
               allocate ( self % f % Q  ( 0 : self % f % spA % N , NEC , 1        )  ) 
               allocate ( self % f % dQ ( 0 : self % f % spA % N , NEC , NDIM , 1 )  ) 

            end if 

        end subroutine Edge_ConstructEdge 

        subroutine Edge_LinkWithElements( self , el1 , el2 , elb)
            implicit none
            class(Edge_p)                          :: self
            class(QuadElement_t), target, optional :: el1
            class(QuadElement_t), target, optional :: el2
            class(QuadElement_t), target, optional :: elb
!           ----------------------------------------------------------------------
            integer                        :: nodesID(POINTS_PER_EDGE)
            integer                        :: nodesEl1(POINTS_PER_QUAD)
            integer                        :: nodesEl2(POINTS_PER_QUAD)
            integer                        :: nodesElb(POINTS_PER_QUAD)
            integer                        :: node
            integer                        :: edgeID
            integer                        :: edgePosition
            integer                        :: quadPosition
            integer                        :: edgeDirection

            do node = 1 , POINTS_PER_EDGE
               nodesID(node) = self % f % nodes(node) % n %ID
            end do

            if (present(el1) .and. present(el2) .and. (.not. present(elb)) ) then             ! Interior edge            

!              Gather all four nodes in both elements
!              --------------------------------------
               do node = 1 , POINTS_PER_QUAD 
                  nodesEl1(node)  = el1 % nodes(node) % n % ID 
                  nodesEl2(node)  = el2 % nodes(node) % n % ID
               end do

!              Search for the edge in element1
!              -------------------------------
               call searchEdge( nodesEl = nodesEl1 , nodesEdge = nodesID , edgePosition = edgePosition , quadPosition = quadPosition , edgeDirection = edgeDirection)
!
!              Link the items
!              --------------
               el1  % edges            ( edgePosition ) % f  => self % f
               el1  % edgesDirection   ( edgePosition ) =  edgeDirection
               el1  % quadPosition     ( edgePosition ) = quadPosition
               self % f % quads        ( quadPosition ) % e  => el1
               self % f % edgeLocation ( quadPosition ) = edgePosition
!
!              Search for the edge in element2
!              -------------------------------
               call searchEdge( nodesEl = nodesEl2 , nodesEdge = nodesID , edgePosition = edgePosition , quadPosition = quadPosition , edgeDirection = edgeDirection)
!
!              Link the items
!              --------------
               el2  % edges            ( edgePosition ) % f  => self % f
               self % f % quads        ( quadPosition ) % e  => el2
               el2  % quadPosition     ( edgePosition ) = quadPosition
               el2  % edgesDirection   ( edgePosition ) =  edgeDirection
               self % f % edgeLocation ( quadPosition ) = edgePosition
               

            elseif (present(elb) .and. (.not. present(el1)) .and. (.not. present(el2))) then ! Boundary edge

               do node = 1 , POINTS_PER_QUAD 
                  nodesElb(node)  = elb % nodes(node) % n % ID 
               end do

!              Search for the edge in element1
!              -------------------------------
               call searchEdge( nodesEl = nodesElb , nodesEdge = nodesID , edgePosition = edgePosition , quadPosition = quadPosition , edgeDirection = edgeDirection)
!
!              Link the items: Now always the direction is FORWARD
!              --------------
               elb  % edges            ( edgePosition ) % f  => self % f
               self % f % quads        ( 1            ) % e  => elb
               elb % quadPosition      ( edgePosition ) = 1
               elb  % edgesDirection   ( edgePosition ) =  FORWARD
               self % f % edgeLocation ( 1            ) =  edgePosition

               if (edgeDirection .eq. BACKWARD) then ! The edge must be inverted
                  call self % f % Invert()
               end if


            end if

         end subroutine Edge_LinkWithElements

         subroutine Edge_Invert ( self )
            implicit none
            class(Edge_t)           :: self
!           --------------------------------------
            class(Node_t), pointer  :: auxnode

!           Invert nodes
            auxnode => self % nodes(1) % n
            self % nodes(1) % n => self % nodes(2) % n
            self % nodes(2) % n => auxnode

!           Invert coordinates
            self % X  ( 1:NDIM , 0:self % spA % N ) = self % X  ( 1:NDIM , self % spA % N : 0 : -1 ) 
            self % dX ( 1:NDIM , 0:self % spA % N ) = self % dX ( 1:NDIM , self % spA % N : 0 : -1 ) 
            self % dS ( 1:NDIM , 0:self % spA % N ) = self % dS ( 1:NDIM , self % spA % N : 0 : -1 ) 

         end subroutine Edge_Invert

end module QuadElementClass
