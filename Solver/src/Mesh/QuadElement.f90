!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        File: QuadElement.f90
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
module QuadElementClass
    use Physics
    use SMConstants
    use NodeClass
    use NodesAndWeights_class
    use Storage_module
    implicit none
!
#include "Defines.h"
!
    private
    public  QuadElement_t , QuadElement_p , Edge_t , CurvedEdge_t , SubdividedEdge_t , StraightBdryEdge_t , CurvedBdryEdge_t , Edge_p
    public  CurvedSubdividedEdge_t
    public  Edge_ProjectSolutionType1
    public  Edge_LinkWithElements , BdryEdge_LinkWithElements , SubdividedEdge_LinkWithElements
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           Quad Element derived type definition
!           ------------------------------------
!
    type QuadElement_t
        integer                            :: ID                             ! ID of the element
        logical                            :: boundaryElement = .false.      ! Whether the element belongs to (at least) a boundary
        integer                            :: address                        ! Memory address of the first position in the mesh
        integer                            :: edgesDirection(EDGES_PER_QUAD) ! Direction (FORWARD/REVERSE) of the edges
        integer                            :: quadPosition(EDGES_PER_QUAD)   ! Position of the quad for the edge (LEFT/RIGHT)
        real(kind=RP), allocatable         :: Volume                         ! Volume of the element
        real(kind=RP)                      :: mu_a
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
!           Pointer to Quad Element derived type
!           ------------------------------------
!
    type QuadElement_p
        type(QuadElement_t),  pointer     :: e  ! Pointer to a quad element
    end type QuadElement_p
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           Derived type to contain both non-conforming information
!           -------------------------------------------------------
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
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!   *********************************************************************************
!           Edge derived type definition:
!                 Edge_t                 : For straight interior edges
!                 CurvedEdge_t           : For curved interior edges
!                 SubdividedEdge_t       : For straight subdivided (interior) edges
!                 SubdividedCurvedEdge_t : For curved subdivided (interior) edges
!                 StraightBdryEdge_t     : For straight boundary edges
!                 CurvedBdryEdge_t       : For curved boundary edges
!   *********************************************************************************
!
!   ============
!   Edge_t class : Main Edge structure containing the basics.
!   ============                                                           ! 
    type Edge_t                                                            ! ------------------------------------------------------------------------------------------
        integer                               :: ID                        ! Edge ID
        integer(kind=1)                       :: edgeType                  ! Edge Type: FACE_INTERIOR , or the marker value if boundary face
        logical                               :: inverse                   ! Whether both edge projections are in the same or different direction
        logical                               :: transform(QUADS_PER_EDGE) ! Whether the element contribution needs transformation to a higher degree
        integer                               :: Nlow                      ! Lower polynomial degree
        integer(kind=1),            pointer   :: edgeLocation(:)           ! Edge location for the two (or one) sharing elements (ETOP,EBOTTOM,ELEFT,ERIGHT)
        real(kind=RP)                         :: Area                      ! Area of the edge
        real(kind=RP)                         :: mu_a                      ! Physical artificial viscosity
        real(kind=RP)                         :: invh                      ! Normal minimum distance approximation (inversed)
        real(kind=RP),              pointer   :: n(:,:)                    ! Unitary normal: points from LEFT towards RIGHT, and outside the domain for bdryedges
        real(kind=RP),              pointer   :: X(:,:)                    ! Coordinates: (X/Y, xi)
        real(kind=RP),              pointer   :: dX(:,:)                   ! Tangent vector: (X/Y, xi)
        real(kind=RP),              pointer   :: dS(:)                     ! Surface differential vector (X/Y, xi)
        real(kind=RP),          allocatable   :: T_forward(:,:)            ! Interpolation matrix from the low order to high order
        real(kind=RP),          allocatable   :: T_backward(:,:)           ! Interpolation matrix from the high order to low order
        type(BoundaryData_t),   allocatable   :: storage(:)                ! Solution interpolation to boundaries ( LEFT/RIGHT )
        type(Node_p),           allocatable   :: nodes(:)                  ! Pointer to the two nodes
        class(QuadElement_p),     pointer     :: quads(:)                  ! Pointers to the two (or one) shared quads
        class(NodesAndWeights_t), pointer     :: spA                       ! Pointer to the approximation nodal storage. In this case, is the largest from the two elements
        class(NodesAndWeights_t), pointer     :: spI                       ! Pointer to the integration nodal storage (if over-integration is active)
        procedure(ProjectSolutionFCN           ), nopass, pointer  :: ProjectSolution            => NULL()
        procedure(ProjectSolutionAndGradientFCN), nopass, pointer  :: ProjectSolutionAndGradient => NULL()
        procedure(ProjectFluxesFCN             ), nopass, pointer  :: ProjectFluxes              => NULL()
        procedure(ProjectGradientFluxesFCN     ), nopass, pointer  :: ProjectGradientFluxes      => NULL()
        contains
            procedure      :: SetCurve     => Edge_SetCurve                ! Procedure that computes the coordinates, the tangent, and the normal.
            procedure      :: Invert       => Edge_Invert                  ! Function to invert the edge orientation
            procedure      :: evaluateX    => Edge_AnalyticalX             ! Function to compute an edge point in a local coordinate "xi"
            procedure      :: evaluatedX   => Edge_AnalyticaldX            ! Function to compute the tangent vector in a local coordinate "xi"
            procedure      :: evaluatedS   => Edge_AnalyticaldS            ! Function to compute the dS vector in a local coordinate "xi"
            procedure      :: getX         => Edge_getX                    ! Function to get one nodal point X
            procedure      :: getdX        => Edge_getdX                   ! Function to get one nodal tangent vector dX
            procedure      :: getdS        => Edge_getdS                   ! Function to get one nodal dS vector 
            procedure      :: ComputeJumps => Edge_ComputeJumps            ! Function to compute the solution jumps across the edge
    end type Edge_t                                                        ! ----------------------------------------------------------------------------------------------
!                                                                          !
!   ==================
!   CurvedEdge_t class : Considers a curved interior edge
!   ==================
    type, extends(Edge_t)  :: CurvedEdge_t
      contains                                                             ! 
         procedure   :: SetCurve     => CurvedEdge_SetCurve                !  ---------------------------------------
         procedure   :: getX         => CurvedEdge_getX                    !        Curves edges just have
         procedure   :: GetdX        => CurvedEdge_GetdX                   !     different procedures to compute
         procedure   :: GetdS        => CurvedEdge_GetdS                   !     the edge coordinates
         procedure   :: EvaluateX    => CurvedEdge_InterpolantX            !
         procedure   :: EvaluatedX   => CurvedEdge_InterpolantdX           !
         procedure   :: EvaluatedS   => CurvedEdge_InterpolantdS           !
         procedure   :: ComputeJumps => CurvedEdge_ComputeJumps            !  ---------------------------------------
    end type CurvedEdge_t                                                  !
!
!   ======================
!   SubdividedEdge_t class : Considers a subdivided straight interior edge
!   ======================                                            !
    type, extends(Edge_t)  :: SubdividedEdge_t                        ! --------------------------------------------------
       integer                            :: N_N                      ! Polynomial order of the NORTH mortar
       integer                            :: N_S                      ! Polynomial order of the SOUTH mortar
       class(NodesAndWeights_t) , pointer :: spA_N                    ! Nodes and weights of the NORTH mortar
       class(NodesAndWeights_t) , pointer :: spA_S                    ! Nodes and weights of the SOUTH mortar
       real(kind=RP), allocatable         :: x_N(:,:)                 ! Coordinates of the NORTH mortar
       real(kind=RP), allocatable         :: x_S(:,:)                 ! Coordinates of the SOUTH mortar
       real(kind=RP), allocatable         :: normal_N(:,:)            ! NORTH MORTAR normal vector
       real(kind=RP), allocatable         :: normal_S(:,:)            ! SOUTH MORTAR normal vector
       real(kind=RP), allocatable         :: dS_N(:)                  ! NORTH MORTAR dS
       real(kind=RP), allocatable         :: dS_S(:)                  ! SOUTH MORTAR dS
       real(kind=RP), allocatable         :: T_LN_FWD(:,:)            ! LEFT to NORTH interpolation matrix
       real(kind=RP), allocatable         :: T_LS_FWD(:,:)            ! LEFT to SOUTH interpolation matrix
       real(kind=RP), allocatable         :: T_LN_BKW(:,:)            ! NORTH to LEFT interpolation matrix
       real(kind=RP), allocatable         :: T_LS_BKW(:,:)            ! SOUTH to LEFT interpolation matrix
       real(kind=RP), allocatable         :: T_RN_FWD(:,:)            ! RIGHT-NORTH to MORTAR interpolation matrix
       real(kind=RP), allocatable         :: T_RS_FWD(:,:)            ! RIGHT-SOUTH to MORTAR interpolation matrix
       real(kind=RP), allocatable         :: T_RN_BKW(:,:)            ! MORTAR to RIGHT-NORTH interpolation matrix
       real(kind=RP), allocatable         :: T_RS_BKW(:,:)            ! MORTAR to RIGHT-SOUTH interpolation matrix
      contains                                                        ! ---------------------------------------------------
         procedure   :: ConstructMortars                     => SubdividedEdge_ConstructMortars
         procedure   :: ComputeMortarsTransformationMatrices => SubdividedEdge_ComputeMortarsTransformationMatrices
         procedure   :: ComputeJumps                         => SubdividedEdge_ComputeJumps
    end type SubdividedEdge_t                                            
!
!   ============================
!   CurvedSubdividedEdge_t class : Considers a subdivided curved interior edge
!   ============================
!
    type, extends(CurvedEdge_t)  :: CurvedSubdividedEdge_t
       integer                            :: N_N                      ! Polynomial order of the NORTH mortar
       integer                            :: N_S                      ! Polynomial order of the SOUTH mortar
       class(NodesAndWeights_t) , pointer :: spA_N                    ! Nodes and weights of the NORTH mortar
       class(NodesAndWeights_t) , pointer :: spA_S                    ! Nodes and weights of the SOUTH mortar
       real(kind=RP), allocatable         :: x_N(:,:)                 ! Coordinates of the NORTH mortar
       real(kind=RP), allocatable         :: x_S(:,:)                 ! Coordinates of the SOUTH mortar
       real(kind=RP), allocatable         :: normal_N(:,:)            ! NORTH MORTAR normal vector
       real(kind=RP), allocatable         :: normal_S(:,:)            ! SOUTH MORTAR normal vector
       real(kind=RP), allocatable         :: dS_N(:)                  ! NORTH MORTAR dS
       real(kind=RP), allocatable         :: dS_S(:)                  ! SOUTH MORTAR dS
       real(kind=RP), allocatable         :: T_LN_FWD(:,:)            ! LEFT to NORTH interpolation matrix
       real(kind=RP), allocatable         :: T_LS_FWD(:,:)            ! LEFT to SOUTH interpolation matrix
       real(kind=RP), allocatable         :: T_LN_BKW(:,:)            ! NORTH to LEFT interpolation matrix
       real(kind=RP), allocatable         :: T_LS_BKW(:,:)            ! SOUTH to LEFT interpolation matrix
       real(kind=RP), allocatable         :: T_RN_FWD(:,:)            ! RIGHT-NORTH to MORTAR interpolation matrix
       real(kind=RP), allocatable         :: T_RS_FWD(:,:)            ! RIGHT-SOUTH to MORTAR interpolation matrix
       real(kind=RP), allocatable         :: T_RN_BKW(:,:)            ! MORTAR to RIGHT-NORTH interpolation matrix
       real(kind=RP), allocatable         :: T_RS_BKW(:,:)            ! MORTAR to RIGHT-SOUTH interpolation matrix 
       contains
         procedure   :: ConstructMortars                     => CurvedSubdividedEdge_ConstructMortars
         procedure   :: ComputeMortarsTransformationMatrices => CurvedSubdividedEdge_ComputeMortarsTransformationMatrices
         procedure   :: ComputeJumps                         => CurvedSubdividedEdge_ComputeJumps
      end type CurvedSubdividedEdge_t
!
!  ========================
!  StraightBdryEdge_t class : Considers a straight boundary edge
!  ========================                                                            !
    type, extends(Edge_t)  :: StraightBdryEdge_t                                       !  -------------------------------------------------------------------------
        integer(kind=1)                   :: inviscidBCType                            !  Inviscid BC type: WEAK_RIEMANN, WEAK_PRESCRIBED
        logical                           :: associated = .false.                      !  Just for periodic boundary conditions
        procedure(RiemannSolverFunction), nopass, pointer :: RiemannSolver => NULL()   !  Edge particular Riemann solver
        real(kind=RP), pointer            :: uB(:,:)  => NULL()                        !  Solution at the boundary (used by the inviscid Riemann solver)
        real(kind=RP), pointer            :: FB(:,:)  => NULL()                        !  Fluxes at the boundary (used for weak-prescribed type boundary conditions)
#ifdef NAVIER_STOKES                                                                   
        integer(kind=1), allocatable      :: viscousBCType(:)                          !  Viscous BC type: DIRICHLET, NEUMANN, PERIODIC, ADIABATIC
        real(kind=RP), pointer            :: uSB(:,:)  => NULL()                       !  Solution at the boundary (used by the solution Riemann solver)
        real(kind=RP), pointer            :: gB(:,:,:) => NULL()                       !  Solution gradient at the boundary
#endif
         contains                                                                      !
            procedure      :: ComputeJumps => StraightBdryEdge_ComputeJumps            !  --------------------------------------------------------------------------
    end type StraightBdryEdge_t                                                        ! 

    type, extends(CurvedEdge_t)  :: CurvedBdryEdge_t
        integer(kind=1)                   :: inviscidBCType
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
    end type Edge_p
!
!   
    abstract interface
         pure subroutine ProjectSolutionFCN( ed , QL , QR )
            use SMConstants
            import Edge_t
            implicit none
            class(Edge_t), intent(in)     :: ed
            real(kind=RP), intent(out)    :: QL( 0 : ed % spA % N , 1 : NCONS )
            real(kind=RP), intent(out)    :: QR( 0 : ed % spA % N , 1 : NCONS ) 
         end subroutine ProjectSolutionFCN
    end interface

    abstract interface
         pure subroutine ProjectSolutionAndGradientFCN( ed , QL , QR , dQL , dQR )
            use SMConstants
            import Edge_t
            implicit none
            class(Edge_t), intent(in)     :: ed
            real(kind=RP), intent(out)    :: QL( 0 : ed % spA % N , 1 : NCONS )
            real(kind=RP), intent(out)    :: QR( 0 : ed % spA % N , 1 : NCONS ) 
            real(kind=RP), intent(out)    :: dQL( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS )
            real(kind=RP), intent(out)    :: dQR( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
         end subroutine ProjectSolutionAndGradientFCN
    end interface

    abstract interface
         pure subroutine ProjectFluxesFCN( ed , F , FL , FR )
            use SMConstants
            import Edge_t
            implicit none
            class(Edge_t), intent(in)     :: ed
            real(kind=RP), intent(in)     :: F ( 0 : ed % spA % N , 1 : NCONS )
            real(kind=RP), intent(out)    :: FL( 0 : ed % storage(LEFT) % spA % N , 1 : NCONS )
            real(kind=RP), intent(out)    :: FR( 0 : ed % storage(RIGHT) % spA % N , 1 : NCONS )
         end subroutine ProjectFluxesFCN
    end interface

    abstract interface
         pure subroutine ProjectGradientFluxesFCN( ed , GL_edge , GR_edge , GL , GR)
            use SMConstants
            import Edge_t
            implicit none
            class(Edge_t), intent(in)     :: ed
            real(kind=RP), intent(in)     :: GL_edge ( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
            real(kind=RP), intent(in)     :: GR_edge ( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
            real(kind=RP), intent(out)    :: GL( 0 : ed % storage(LEFT) % spA % N , 1 : NCONS , 1 : NDIM )
            real(kind=RP), intent(out)    :: GR( 0 : ed % storage(RIGHT) % spA % N , 1 : NCONS , 1 : NDIM )
         end subroutine ProjectGradientFluxesFCN
    end interface
!
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
#include "./QuadElement_Auxiliars.incf"
#include "./QuadElement_EdgeProcedures.incf"
#include "./QuadElement_EdgeMappings.incf"
#include "./QuadElement_QuadProcedures.incf"
#include "./QuadElement_QuadMappings.incf"
!
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!                 Construct element procedure
!                 ---------------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        subroutine QuadElement_Construct(self , ID , nodes , N  , spA , address , storage , spI)
!
!            ***********************************************************************************
!
!                 This function performs the following operations to "self" element:
!                       -> Assign its ID and storage address position
!                       -> Assign its nodal storage spA, and spI if proceeds
!                       -> Assign its nodes
!                       -> Assign its storage
!
!            ***********************************************************************************
!
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
!            ---------------
!            Local variables
!            ---------------
!
             integer                           :: edge
!
!            Assign its ID
!            -------------
             self % ID = ID
!
!            Point to nodes
!            --------------
             self % nodes = nodes
!
!            Get NodalStorage data       
!            ---------------------
             call spA % add( N , Setup % nodes , self % spA )
             self % spI => spI                                              ! Points to NULL if Standard DG is chosen
!
!            ************
!            Linking data
!            ************
!
!            -------------------------------------------------------------
!               These two arrays were before 
!               allocated inside each element.
!               Now they are just linked.
!                   allocate ( self % Q    ( 0:N , 0:N , NCONS        )  ) 
!                   allocate ( self % QDot ( 0:N , 0:N , NCONS        )  ) 
!                   allocate ( self % dQ   ( 0:N , 0:N , NDIM , NGRAD )  ) 
!            -------------------------------------------------------------
!
             self % address = address
             call self % SetStorage( storage )
!
!            Allocate data
!            -------------
             allocate ( self % x         ( 0:N   , 0:N , NDIM        )  ) 
             allocate ( self % Ja        ( 0:N   , 0:N , NDIM , NDIM )  ) 
             allocate ( self % jac       ( 0:N   , 0:N               )  ) 
             allocate ( self % edges     ( EDGES_PER_QUAD            )  ) 
      
             do edge = 1 , EDGES_PER_QUAD
               self % edges(edge) % f => NULL()
             end do
             
        end subroutine QuadElement_Construct
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!                 Construct edge procedure
!                 ------------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        subroutine Edge_ConstructEdge( self , ID , curvilinear , N , NLow , nodes , edgeType , spA , spI)
            use Setup_Class
            implicit none
            class(Edge_p)                     :: self
            integer                           :: ID
            class(Node_p)                     :: nodes(:)
            logical                           :: curvilinear
            integer                           :: N 
            integer                           :: NLow
            integer(kind=1)                   :: edgeType
            integer                           :: node
            class(NodalStorage)               :: spA
            class(NodesAndWeights_t), pointer :: spI
!
!           ---------------
!           Local variables
!           ---------------
!
            integer                           :: quad
!
!           *************************************************
!              Allocate the edge depending on its type
!           *************************************************
!
            if ( (edgeType .EQ. FACE_INTERIOR) .and. (size(nodes) .eq. TWO) ) then
!
!> Simple edge
!
!             * Straight simple interior edge
!               -----------------------------
                if ( .not. curvilinear ) then
                  allocate(Edge_t :: self % f)
!
!             * Curved simple interior edge
!               ---------------------------
                else
                  allocate(CurvedEdge_t   :: self % f )

                end if

!               Allocate its elements 
!               ---------------------
                allocate ( self % f % nodes        ( POINTS_PER_EDGE )  ) 
                allocate ( self % f % quads        ( QUADS_PER_EDGE  )  ) 
                allocate ( self % f % edgeLocation ( QUADS_PER_EDGE  )  ) 
!
!              Link the nodes    TODO: check if possible to self % f % nodes = nodes (similar to elements!)
!              --------------
               do node = 1 , POINTS_PER_EDGE
                  self % f % nodes(node) % n => nodes(node) % n
               end do

                do quad = 1 , QUADS_PER_EDGE
                  self % f % quads(quad) % e => NULL()
                end do

                self % f % edgeType = FACE_INTERIOR

            elseif ( (edgeType .eq. FACE_INTERIOR) .and. (size(nodes) .eq. THREE) ) then
!
!> Subdivided edge
!
!            * Straight divided interior edge
!              ------------------------------
               if ( .not. curvilinear ) then
                  allocate( SubdividedEdge_t    :: self % f )
!
!            * Curved divided interior edge
!              ----------------------------
               else
                  allocate( CurvedSubdividedEdge_t  :: self % f )

               end if
!
!              Allocate its elements
!              ---------------------
               allocate ( self % f % nodes        ( POINTS_PER_SUBDIVIDED_EDGE )  ) 
               allocate ( self % f % quads        ( QUADS_PER_SUBDIVIDED_EDGE  )  ) 
               allocate ( self % f % edgeLocation ( QUADS_PER_SUBDIVIDED_EDGE  )  ) 
!
!              Link the nodes    TODO: check if possible to self % f % nodes = nodes (similar to elements!)
!              --------------
               do node = 1 , POINTS_PER_SUBDIVIDED_EDGE
                  self % f % nodes(node) % n => nodes(node) % n
               end do

               do quad = 1 , QUADS_PER_SUBDIVIDED_EDGE
                  self % f % quads(quad) % e => NULL()
               end do

               self % f % edgeType = FACE_INTERIOR

            elseif (edgeType .NE. FACE_INTERIOR) then
!
!> Boundary edges
!
!            * Straight boundary edges
!              -----------------------
               if ( .NOT. curvilinear ) then
                  allocate(StraightBdryEdge_t   :: self % f)

!
!            * Curvilinear boundary edges
!              --------------------------
               else
                  allocate(CurvedBdryEdge_t     :: self % f)

               end if

               allocate ( self % f % nodes        ( POINTS_PER_EDGE )  ) 
               allocate ( self % f % quads        ( ONE             )  ) 
               allocate ( self % f % edgeLocation ( ONE             )  ) 
!
!              Link the nodes    TODO: check if possible to self % f % nodes = nodes (similar to elements!)
!              --------------
               do node = 1 , POINTS_PER_EDGE
                  self % f % nodes(node) % n => nodes(node) % n
               end do

               self % f % quads(ONE) % e => NULL()
               self % f % edgeType       =  edgeType

            end if
!
!           Assign its ID
!           -------------
            self % f % ID = ID
!
!           Compute the nodes and weights
!           -----------------------------
            call spA % add( N , Setup % nodes , self % f % spA )
            self % f % spI    => spI
            self % f % Nlow = NLow
!
!           Allocate the geometry
!           ---------------------
            allocate ( self % f % X    ( NDIM , 0 : self % f % spA % N       )  ) 

            select type ( f => self % f )

               type is (Edge_t)
                  allocate ( f % dX ( NDIM ,0:0 )  ) 
                  allocate ( f % dS (       0:0 )  ) 
                  allocate ( f % n  ( NDIM ,0:0 )  ) 

               type is (CurvedEdge_t) 
                  allocate ( f % dX ( NDIM , 0 : f % spA % N )  )
                  allocate ( f % dS (        0 : f % spA % N )  )
                  allocate ( f % n  ( NDIM , 0 : f % spA % N )  )

               type is (SubdividedEdge_t)
                  allocate ( f % dX ( NDIM ,0:0 )  ) 
                  allocate ( f % dS (       0:0 )  ) 
                  allocate ( f % n  ( NDIM ,0:0 )  ) 

               type is (CurvedSubdividedEdge_t)
                  allocate ( f % dX ( NDIM , 0 : f % spA % N )  )
                  allocate ( f % dS (        0 : f % spA % N )  )
                  allocate ( f % n  ( NDIM , 0 : f % spA % N )  )

               type is (StraightBdryEdge_t)
                  allocate ( f % dX ( NDIM , 0:0 )  ) 
                  allocate ( f % dS (        0:0 )  ) 
                  allocate ( f % n  ( NDIM , 0:0 )  ) 
#ifdef NAVIER_STOKES
                  allocate ( f % viscousBCType( 0 : f % spA % N ) )
#endif

               type is (CurvedBdryEdge_t)
                  allocate ( f % dX ( NDIM , 0 : f % spA % N )  )
                  allocate ( f % dS (        0 : f % spA % N )  )
                  allocate ( f % n  ( NDIM , 0 : f % spA % N )  )
#ifdef NAVIER_STOKES
                  allocate ( f % viscousBCType( 0 : f % spA % N ) )
#endif
            end select

        end subroutine Edge_ConstructEdge 
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
end module QuadElementClass
