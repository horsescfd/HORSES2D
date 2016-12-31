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
        integer                           :: edgesDirection(EDGES_PER_QUAD)
        real(kind=RP), allocatable        :: x(:,:,:)
        real(kind=RP), allocatable        :: dx(:,:,:,:)
        real(kind=RP), allocatable        :: jac(:,:)
        real(kind=RP), pointer            :: Q(:,:,:) , QDot(:,:,:) , F(:,:,:,:) , dQ(:,:,:,:)
        type(Node_p)                      :: nodes(POINTS_PER_QUAD)
        class(Edge_p), pointer            :: edges(:)
        class(NodesAndWeights_t), pointer :: spA
        class(NodesAndWeights_t), pointer :: spI
!       ========
        contains
!       ========
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
        real(kind=RP), allocatable        :: nb(:,:)          !   n always points from LEFT towards RIGHT and outside the domain for bdryedges
        type(Node_p)                      :: nodes(POINTS_PER_EDGE)
        class(QuadElement_p), pointer     :: quads(:)
        real(kind=RP), allocatable        :: X(:,:)
        real(kind=RP), allocatable        :: dS(:,:)
        real(kind=RP), allocatable        :: Q(:,:,:) , dQ(:,:,:,:)  ! To store the interpolation to boundaries from elements
        integer,       allocatable        :: edgeLocation(:)
        integer                           :: edgeType
        class(NodesAndWeights_t), pointer :: spA
        class(NodesAndWeights_t), pointer :: spI
        contains
            procedure      :: SetCurve => Edge_SetCurve
            procedure      :: Invert => Edge_Invert
            procedure      :: XF       => Edge_AnalyticalX
            procedure      :: dSF      => Edge_AnalyticaldS
    end type Edge_t

    type, extends(Edge_t)  :: StraightBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)           ! Solution at the boundary
        real(kind=RP), pointer            :: gB(:,:,:)         ! Solution gradient at the boundary
    end type StraightBdryEdge_t 

    type, extends(Edge_t)  :: CurvedBdryEdge_t
        real(kind=RP), pointer            :: uB(:,:)           ! Solution at the boundary
        real(kind=RP), pointer            :: gB(:,:,:)         ! Solution gradient at the boundary
        contains
            procedure      :: SetCurve => CurvilinearEdge_SetCurve
            procedure      :: XF       => Curvilinear_InterpolantX
 !           procedure      :: dSF      => Curvilinear_InterpolantdS
    end type CurvedBdryEdge_t

    type Edge_p
        class(Edge_t),   pointer           :: f
        contains
            procedure      :: Construct => Edge_ConstructEdge
            procedure      :: LinkWithElements => Edge_linkWithElements
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
             allocate ( self % edges ( EDGES_PER_QUAD ) )
      
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
            self % f % N  = Setup % N

            do node = 1 , POINTS_PER_EDGE
               self % f % nodes(node) % n => nodes(node) % n
            end do

            call spA % add( self % f % N , Setup % nodes , self % f % spA )
            self % f % spI    => spI

            allocate ( self % f % X  ( NDIM , 0 : self % f % spA % N )  ) 
            allocate ( self % f % dS ( NDIM , 0 : self % f % spA % N )  ) 
            allocate ( self % f % nb ( NDIM , 0 : self % f % spA % N )  ) 

!
!           If rectilinear edge, compute the points
            select type (f => self % f)
               type is (Edge_t)
                   
                  self % f % x = reshape((/( self % f % nodes(1) % n % X * (1.0_RP - self % f % spA % xi(p)) + self % f % nodes(2) % n % X * self % f % spA % xi(p) , &
                                                p = 0 , self % f % spA % N)/),(/ NDIM , self % f % spA % N + 1 /) )
               type is (StraightBdryEdge_t)

                  self % f % x = reshape((/( self % f % nodes(1) % n % X * (1.0_RP - self % f % spA % xi(p)) + self % f % nodes(2) % n % X * self % f % spA % xi(p) , &
                                                p = 0 , self % f % spA % N)/),(/ NDIM , self % f % spA % N + 1 /) )

               class default
            end select 

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
            self % X(1:NDIM , 0:self % spA % N) = self % X(1:NDIM , self % spA % N : 0 : -1)

         end subroutine Edge_Invert

         subroutine Edge_SetCurve( self , points , order )
            implicit none
            class(Edge_t)                :: self
            real(kind=RP), intent(in)    :: points(:,:)
            integer      , intent(in)    :: order
!
!           **************************************
!              The base class does nothing
!           **************************************
!
         end subroutine Edge_SetCurve

         subroutine CurvilinearEdge_SetCurve( self , points , order )
            use InterpolationAndDerivatives
            use MatrixOperations
            implicit none
            class(CurvedBdryEdge_t)                            :: self
            real(kind=RP)               , intent(in)           :: points(:,:)
            integer                     , intent(in)           :: order
!           -----------------------------------------------------------------------------
            real(kind=RP), allocatable                         :: CGLnodes(:)
            real(kind=RP), allocatable                         :: T(:,:)
            real(kind=RP), allocatable                         :: wb(:)
            integer                                            :: node

            allocate( CGLnodes(0 : order ) )
            allocate( wb(0 : order ) )
            allocate( T(0: self % spA % N , 0: order ) )

            CGLnodes = reshape ( (/(0.5_RP + 0.5_RP*cos(PI*(order - node)/(1.0_RP*order)),node = 0,order)/),(/order+1/) )

            call BarycentricWeights( N = order , x = CGLnodes , w = wb )
            call PolynomialInterpolationMatrix( N = order , M = self % spA % N, oldNodes = CGLnodes, weights = wb, newNodes = self % spA % xi , T = T)

            self % X = NormalMat_x_TransposeMat_F( points , T )
            !self % X = matmul( points , transpose(T) )

         end subroutine CurvilinearEdge_SetCurve
      
         function Edge_AnalyticalX( self , xi , direction ) result( p )
            implicit none
            class(Edge_t), intent(in)           :: self
            real(kind=RP), intent(in)           :: xi
            integer      , intent(in)           :: direction
            real(kind=RP)                       :: p(2)
!           ------------------------------------------------------------
            real(kind=RP)                       :: correctedXi
            
            if (direction .eq. BACKWARD) then
              correctedXi = 1.0_RP - xi
            elseif (direction .eq. FORWARD) then
              correctedXi = xi
            end if 

            p = self % nodes(1) % n % X * (1.0_RP - correctedXi) + self % nodes(2) % n % X * correctedXi

         end function Edge_AnalyticalX

         function Curvilinear_InterpolantX( self , xi , direction ) result( p )
            use MatrixOperations
            implicit none
            class(CurvedBdryEdge_t), intent (in) :: self
            real(kind=RP),           intent (in) :: xi
            integer      ,           intent (in) :: direction
            real(kind=RP)                        :: p(2)
!           ------------------------------------------------------------
            real(kind=RP), allocatable           :: auxp(:,:)
            real(kind=RP)                        :: correctedXi
            real(kind=RP), allocatable           :: lj(:,:)
            
            if (direction .eq. BACKWARD) then
              correctedXi = 1.0_RP - xi
            elseif (direction .eq. FORWARD) then
              correctedXi = xi
            end if 

            allocate(lj( 0 : self % spA % N , 1 ) )
            allocate(auxp(NDIM , 1) )
            lj(:,1) = self % spA % lj(correctedXi)

            auxp(1:NDIM,1:) = MatrixMultiply_F( self % X , lj )

            p = auxp(1:NDIM,1)
            
            deallocate(lj , auxp)

         end function Curvilinear_InterpolantX

         function Edge_AnalyticaldS( self , xi , direction ) result( dS )
            implicit none
            class(Edge_t), intent(in)        :: self
            real(kind=RP), intent(in)        :: xi
            integer,       intent(in)        :: direction
            real(kind=RP)                    :: dS(2)
!           --------------------------------------------------------------

            associate( n1 => self % nodes(1) % n % X , &
                       n2 => self % nodes(2) % n % X )

               dS(1) = n2(2) - n1(2)
               dS(2) = n2(1) - n1(1)

            end associate

         end function Edge_AnalyticaldS

         function Curvilinear_InterpolantdS( self , xi , direction ) result( dS )
            use MatrixOperations
            implicit none
            class(Edge_t), intent(in)        :: self
            real(kind=RP), intent(in)        :: xi
            integer,       intent(in)        :: direction
            real(kind=RP)                    :: dS(2)
!           --------------------------------------------------------------
            real(kind=RP), allocatable           :: dP(:,:)
            real(kind=RP), allocatable           :: auxdS(:,:)
            real(kind=RP), allocatable           :: lj(:,:)
            real(kind=RP)                        :: correctedXi
            
            if (direction .eq. BACKWARD) then
              correctedXi = 1.0_RP - xi
            elseif (direction .eq. FORWARD) then
              correctedXi = xi
            end if 

            allocate(dP(NDIM , 0 : self % spA % N) )
            allocate(lj(0 : self % spA % N , 1 ) )
            allocate(auxdS(NDIM , 1) )

            lj(:,1) = self % spA % lj(correctedXi)

            associate ( D => self % spA % D )
            
               dP = NormalMat_x_TransposeMat_F( self % X , D )
               auxdS(1:NDIM,1:) = MatrixMultiply_F( dP , lj ) 
            
            end associate

            dS = auxdS(1:NDIM , 1)
         
            deallocate( lj , dP , auxdS )

         end function Curvilinear_InterpolantdS




!
!        **********************************************************************************
!                 Auxiliar subroutines
!        **********************************************************************************
!
         subroutine searchEdge( nodesEl , nodesEdge , edgePosition , quadPosition , edgeDirection )
            implicit none
            integer, intent(in)        :: nodesEl(:)
            integer, intent(in)        :: nodesEdge(:)
            integer, intent(out)       :: edgePosition
            integer, intent(out)       :: quadPosition
            integer, intent(out)       :: edgeDirection
!           -----------------------------------------------------------
            integer                    :: currentEdge(POINTS_PER_EDGE)
            integer                    :: edge
            logical                    :: edges_are_equal

            do edge = 1 , EDGES_PER_QUAD
!
!              Obtain current edge of the element
!              ----------------------------------
               currentEdge(1) = nodesEl(edge)
               if (edge .eq. EDGES_PER_QUAD) then
                  currentEdge(2) = nodesEl(1)
               else
                  currentEdge(2) = nodesEl(edge+1)
               end if
!
!              Compare the element edge with the original edge
!              -----------------------------------------------
               call compareEdges( edge1 = currentEdge , edge2 = nodesEdge , edges_are_equal = edges_are_equal , edgeDirection = edgeDirection )

               if (edges_are_equal) then
                  edgePosition = edge
                  if ( edgeDirection .eq. FORWARD ) then
                     quadPosition = LEFT
                  else
                     quadPosition = RIGHT
                  end if

                  return

               end if

            end do   

         end subroutine searchEdge

         subroutine compareEdges ( edge1 , edge2 , edges_are_equal , edgeDirection )
            implicit none
            integer, intent(in)        :: edge1(:)
            integer, intent(in)        :: edge2(:)
            logical, intent(out)       :: edges_are_equal
            integer, intent(out)       :: edgeDirection

            if ((edge1(1) .eq. edge2(1)) .and. (edge1(2) .eq. edge2(2))) then
               edges_are_equal = .true.
               edgeDirection = FORWARD
            
            elseif ((edge1(1) .eq. edge2(2)) .and. (edge1(2) .eq. edge2(1))) then
               edges_are_equal = .true.
               edgeDirection = BACKWARD

            else
               edges_are_equal = .false.
               edgeDirection = 0
            end if

         end subroutine compareEdges

end module QuadElementClass
