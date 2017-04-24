module QuadElementClass
    use Physics
    use SMConstants
    use NodeClass
    use NodesAndWeights_class
    use Storage_module
    implicit none

#include "Defines.h"

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
        logical                            :: boundaryElement = .false.      ! Whether the element belongs to (at least) a boundary
        integer                            :: address                        ! Memory address of the first position in the mesh
        integer                            :: edgesDirection(EDGES_PER_QUAD) ! Direction (FORWARD/REVERSE) of the edges
        integer                            :: quadPosition(EDGES_PER_QUAD)   ! Position of the quad for the edge (LEFT/RIGHT)
        real(kind=RP), allocatable         :: Volume                         ! Volume of the element
        real(kind=RP), allocatable         :: x(:,:,:)                       ! Coordinates of the nodes ( xi , eta , X/Y )
        real(kind=RP), allocatable         :: Ja(:,:,:,:)                    ! Contravariant system metric matrix ( xi , eta , IROW , ICOL)
        real(kind=RP), allocatable         :: Jac(:,:)                       ! Mapping jacobian ( xi , eta )
        real(kind=RP), pointer             :: Q(:,:,:)                       ! Pointers to the main storage:
        real(kind=RP), pointer             :: QDot(:,:,:)                    ! *   Q, QDot ( xi , eta , eq ): solution and time derivative
#ifdef NAVIER_STOKES
        real(kind=RP), pointer             :: dQ(:,:,:,:)                    ! *   dQ( xi ,eta , X/Y , eq):   solution gradient
#endif
        type(Node_p)                       :: nodes(POINTS_PER_QUAD)         ! Pointer to neighbour nodes
        class(Edge_p), pointer             :: edges(:)                       ! Pointer to neighbour eges
        class(NodesAndWeights_t), pointer  :: spA                            ! Pointer to the Nodal Storage
        class(NodesAndWeights_t), pointer  :: spI                            ! Pointer to the Over-Integration Nodal Storage (If proceeds)
!       ========
        contains
!       ========
            procedure      :: Construct                 => QuadElement_Construct                                 ! Constructs/allocates data and points to the address in Storage
            procedure      :: SetStorage                => QuadElement_SetStorage                                ! Function to set the storage distribution
            procedure      :: SetMappings               => QuadElement_SetMappings                               ! Function to compute the mapping data (x, dx, jac)
#ifdef NAVIER_STOKES
            procedure      :: ComputeInteriorGradient   => QuadElement_ComputeInteriorGradient
#endif
            procedure      :: Compute_X                 => QuadElement_Compute_X
            procedure      :: FindPointWithCoords       => QuadElement_FindPointWithCoords
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
!   *******************************************************************************
!           Derived type to contain both non-conforming information
!   *******************************************************************************
!
   type BoundaryData_t
        class(NodesAndWeights_t) , pointer     :: spA
        real(kind=RP)            , pointer     :: Q   (:,:)        ! Solution interpolation to boundaries ( xi , eq )
#ifdef NAVIER_STOKES
        real(kind=RP)            , pointer :: dQ  (:,:,:)      ! Solution gradient interpolation to boundary ( xi , eq , X/Y )
#endif
        contains
            procedure      :: Initialize  => BoundaryData_Initialize
   end type BoundaryData_t

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
        integer(kind=1)                     :: edgeType                   ! Edge Type: FACE_INTERIOR , or the marker value if boundary face
        logical                             :: inverse                    ! Whether both edge projections are in the same or different direction
        logical                             :: transform(QUADS_PER_EDGE)  ! Whether the element contribution needs transformation to a higher degree 
        integer                             :: Nlow                       ! Lower polynomial degree
        integer(kind=1),            pointer :: edgeLocation(:)            ! Edge location for the two (or one) sharing elements (ETOP,EBOTTOM,ELEFT,ERIGHT)
        real(kind=RP)                       :: Area                       ! Area of the edge
        real(kind=RP)                       :: invh                       ! Normal minimum distance approximation (inversed)
        real(kind=RP),              pointer :: n(:,:)                     ! Unitary normal: points from LEFT towards RIGHT, and outside the domain for bdryedges
        real(kind=RP),              pointer :: X(:,:)                     ! Coordinates: (X/Y, xi)
        real(kind=RP),              pointer :: dX(:,:)                    ! Tangent vector: (X/Y, xi)
        real(kind=RP),              pointer :: dS(:)                      ! Surface differential vector (X/Y, xi)
        real(kind=RP),          allocatable :: T_forward(:,:)             ! Interpolation matrix from the low order to high order
        real(kind=RP),          allocatable :: T_backward(:,:)            ! Interpolation matrix from the high order to low order
        type(BoundaryData_t),  allocatable  :: storage(:)                 ! Solution interpolation to boundaries ( LEFT/RIGHT )
        type(Node_p)                        :: nodes(POINTS_PER_EDGE)     ! Pointer to the two nodes
        class(QuadElement_p),       pointer :: quads(:)                   ! Pointers to the two (or one) shared quads
        class(NodesAndWeights_t),   pointer :: spA                        ! Pointer to the approximation nodal storage. In this case, is the largest from the two elements
        class(NodesAndWeights_t),   pointer :: spI                        ! Pointer to the integration nodal storage (if over-integration is active)
        contains
            procedure      :: SetCurve     => Edge_SetCurve                    ! Procedure that computes the coordinates, the tangent, and the normal.
            procedure      :: Invert       => Edge_Invert                      ! Function to invert the edge orientation
            procedure      :: evaluateX    => Edge_AnalyticalX                 ! Function to compute an edge point in a local coordinate "xi"
            procedure      :: evaluatedX   => Edge_AnalyticaldX                ! Function to compute an edge point in a local coordinate "xi"
            procedure      :: evaluatedS   => Edge_AnalyticaldS
            procedure      :: getX         => Edge_getX
            procedure      :: getdX        => Edge_getdX
            procedure      :: getdS        => Edge_getdS
            procedure      :: ComputeJumps => Edge_ComputeJumps
    end type Edge_t

    type, extends(Edge_t)  :: StraightBdryEdge_t
        integer(kind=1)                   :: BCWeakType
        logical                           :: associated = .false.
        procedure(RiemannSolverFunction), nopass, pointer :: RiemannSolver => NULL()
        real(kind=RP), pointer            :: uB(:,:)  => NULL()      ! Solution at the boundary (used by the inviscid Riemann solver)
        real(kind=RP), pointer            :: FB(:,:)  => NULL()      ! Fluxes at the boundary (used for weak-prescribed type boundary conditions)
#ifdef NAVIER_STOKES
        integer(kind=1), allocatable      :: viscousBCType(:)
        real(kind=RP), pointer            :: uSB(:,:)  => NULL()     ! Solution at the boundary (used by the solution Riemann solver)
        real(kind=RP), pointer            :: gB(:,:,:) => NULL()     ! Solution gradient at the boundary
#endif
         contains
            procedure      :: ComputeJumps => StraightBdryEdge_ComputeJumps
    end type StraightBdryEdge_t 

    type, extends(Edge_t)  :: CurvedBdryEdge_t
        integer(kind=1)                   :: BCWeakType
        logical                           :: associated = .false.
        procedure(RiemannSolverFunction), nopass, pointer :: RiemannSolver => NULL()
        real(kind=RP), pointer            :: uB(:,:)  => NULL()     ! Solution at the boundary (used by the Riemann solver)
        real(kind=RP), pointer            :: FB(:,:)  => NULL()     ! Fluxes at the boundary (used for weak-prescribed type boundary conditions)
#ifdef NAVIER_STOKES
        integer(kind=1), allocatable      :: viscousBCType(:)
        real(kind=RP), pointer            :: uSB(:,:)  => NULL()    ! Solution at the boundary (used by the viscous solver)
        real(kind=RP), pointer            :: gB(:,:,:) => NULL()    ! Solution gradient at the boundary
#endif
        contains
            procedure      :: SetCurve     => CurvilinearEdge_SetCurve       ! Procedure that computes the coordinates, the tangent, and the normal
            procedure      :: getX         => Curvilinear_getX
            procedure      :: getdX        => Curvilinear_getdX
            procedure      :: getdS        => Curvilinear_getdS
            procedure      :: evaluateX    => Curvilinear_InterpolantX
            procedure      :: evaluatedX   => Curvilinear_InterpolantdX
            procedure      :: evaluatedS   => Curvilinear_InterpolantdS
            procedure      :: ComputeJumps => CurvedBdryEdge_ComputeJumps
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
!  ========
   contains
!  ========
!
#include "./QuadElement_Auxiliars.incf"
#include "./QuadElement_EdgeProcedures.incf"
#include "./QuadElement_EdgeMappings.incf"
#include "./QuadElement_QuadProcedures.incf"
#include "./QuadElement_QuadMappings.incf"

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
!                   allocate ( self % Q    ( 0:N , 0:N , NCONS        )  ) 
!                   allocate ( self % QDot ( 0:N , 0:N , NCONS        )  ) 
!                   allocate ( self % F    ( 0:N , 0:N , NCONS        )  ) 
!                   allocate ( self % dQ   ( 0:N , 0:N , NDIM , NGRAD )  ) 
!            ---------------------------------------
!
             self % address = address
             call self % SetStorage( storage )
!
!            *************
!            Allocate data
!            *************
!
             allocate ( self % x         ( 0:N   , 0:N , NDIM        )  ) 
             allocate ( self % Ja        ( 0:N   , 0:N , NDIM , NDIM )  ) 
             allocate ( self % jac       ( 0:N   , 0:N               )  ) 
             allocate ( self % edges     ( EDGES_PER_QUAD            )  ) 
      
             do edge = 1 , EDGES_PER_QUAD
               self % edges(edge) % f => NULL()
             end do
             
        end subroutine QuadElement_Construct

        subroutine Edge_ConstructEdge( self , ID , curvilinear , N , nodes , edgeType , spA , spI)
            use Setup_Class
            implicit none
            class(Edge_p)                     :: self
            integer                           :: ID
            class(Node_p)                     :: nodes(:)
            logical                           :: curvilinear
            integer                           :: N 
            integer(kind=1)                   :: edgeType
            integer                           :: node
            class(NodalStorage)               :: spA
            class(NodesAndWeights_t), pointer :: spI
!           --------------------------------------------------
            integer                           :: quad
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

            call spA % add( N , Setup % nodes , self % f % spA )
            self % f % spI    => spI
            self % f % Nlow = N

!           Geometry
!           --------
            allocate ( self % f % X    ( NDIM , 0 : self % f % spA % N       )  ) 

            select type ( f => self % f )

               type is (Edge_t)
                  allocate ( self % f % dX ( NDIM ,0:0 )  ) 
                  allocate ( self % f % dS (       0:0 )  ) 
                  allocate ( self % f % n  ( NDIM ,0:0 )  ) 

               type is (StraightBdryEdge_t)
                  allocate ( self % f % dX ( NDIM , 0:0 )  ) 
                  allocate ( self % f % dS (        0:0 )  ) 
                  allocate ( self % f % n  ( NDIM , 0:0 )  ) 
#ifdef NAVIER_STOKES
                  allocate ( f % viscousBCType( 0 : self % f % spA % N ) )
#endif

               type is (CurvedBdryEdge_t)
                  allocate ( self % f % dX ( NDIM , 0 : self % f % spA % N )  )
                  allocate ( self % f % dS (        0 : self % f % spA % N )  )
                  allocate ( self % f % n  ( NDIM , 0 : self % f % spA % N )  )
#ifdef NAVIER_STOKES
                  allocate ( f % viscousBCType( 0 : self % f % spA % N ) )
#endif

            end select

        end subroutine Edge_ConstructEdge 

end module QuadElementClass
