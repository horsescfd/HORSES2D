module DGSpatialDiscretizationMethods
   use SMConstants
   use Physics
   use QuadMeshClass
   use DGInviscidMethods
   use DGWeakIntegrals
#ifdef NAVIER_STOKES
   use DGArtificialDissipation
   use DGViscousMethods
#endif
   implicit none
!
#include "Defines.h"

   private
   public DGSpatial_Initialization  , DGSpatial_computeTimeDerivative , DGSpatial_interpolateSolutionToBoundaries
   public DGSpatial_newTimeStep
!
!  ************************************
!  Inviscid and Viscous methods objects
!  ************************************
!
   class(InviscidMethod_t), pointer :: InviscidMethod
#ifdef NAVIER_STOKES
   class(ViscousMethod_t),         pointer  :: ViscousMethod
   class(ArtificialDissipation_t), pointer  :: ArtificialDissipation
#endif

   interface DGSpatial_QDotFaceLoop
      module procedure DGSpatial_QDotFaceLoop_StraightBdry, DGSpatial_QDotFaceLoop_CurvedBdry , DGSpatial_QDotFaceLoop_Interior
   end interface DGSpatial_QDotFaceLoop
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              INITIALIZATION
!              --------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_Initialization()
         use Setup_class
         use Headers
         implicit none

         write(STD_OUT , '(/)')
         call Section_header("Spatial discretization overview")
!
!        Initialize Inviscid method
!        --------------------------
         InviscidMethod => InviscidMethod_Initialization()        
!
!        Initialize Viscous method
!        -------------------------
#ifdef NAVIER_STOKES
         ViscousMethod         => ViscousMethod_Initialization()
         ArtificialDissipation => ArtificialDissipation_Initialization()
#endif
  
      end subroutine DGSpatial_Initialization
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!     
!              COMPUTE TIME DERIVATIVE
!              -----------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_computeTimeDerivative( mesh ) 
!
!        ***************************************************
!           Subroutine that performs the spatial
!         discretization and computes the time
!         derivative QDot.
!        ***************************************************
!
         implicit none
         class(QuadMesh_t)         :: mesh
!
!        Prepare the mesh for a new iteration
!        ------------------------------------
         call DGSpatial_newTimeStep( mesh )
!
!        Compute QDot
!        ------------
         call DGSpatial_computeQDot( mesh )

      end subroutine DGSpatial_computeTimeDerivative
      
         subroutine DGSpatial_newTimeStep( mesh )
!
!        *************************************************************
!           This subroutine prepares the mesh struct
!          for a new time-step. 
!              1) Set QDot to zero
!              2) Interpolate the solution to boundaries
!              3) Update the boundary zones
!              4) Compute the solution gradient
!              5) Interpolate the gradient to boundaries
!        *************************************************************
!
         implicit none
         class(QuadMesh_t)         :: mesh
         integer                   :: zoneID
!
!        Interpolate solution to boundaries
!        ----------------------------------
         call DGSpatial_interpolateSolutionToBoundaries( mesh )
!
!        Update the zones solution
!        -------------------------
         do zoneID = 1 , size(mesh % zones) - 1
            call mesh % zones(zoneID) % UpdateSolution
         end do 

#ifdef NAVIER_STOKES
!
!        Compute the solution Q gradient dQ
!        ----------------------------------
         call ViscousMethod % ComputeGradient( mesh ) 
!
!        Interpolate gradient to boundaries
!        ----------------------------------
         call DGSpatial_interpolateGradientsToBoundaries( mesh )
#endif

      end subroutine DGSpatial_newTimeStep

      subroutine DGSpatial_computeQDot( mesh )
!
!        *************************************************************
!              Once the mesh has been prepared for the Time Derivative
!           calculations, it is performed in this routine.
!           This is performed in three stages:
!              1) Volume loop ( Inviscid & Viscous )
!              2) Face loop ( Inviscid & Viscous )
!              3) Scaling
!        *************************************************************
!  
         use QuadElementClass
         use Setup_class
         implicit none
         class(QuadMesh_t)         :: mesh
!        -------------------------------
         integer                 :: eID
         integer                 :: edID
         integer                 :: eq
!
!        ************
!        Volume loops
!        ************
!
         do eID = 1 , mesh % no_of_elements
            call DGSpatial_QDotVolumeLoop( mesh % elements(eID) ) 
         end do
!
!        **********
!        Face loops
!        **********
!
         do edID = 1 , mesh % no_of_edges
            select type ( f => mesh % edges(edID) % f ) 
               type is (Edge_t)
                  call DGSpatial_QDotFaceLoop_Interior( f ) 

               type is (StraightBdryEdge_t)
                  call DGSpatial_QDotFaceLoop_StraightBdry( f )
  
               type is (CurvedBdryEdge_t)
                  call DGSpatial_QDotFaceLoop_CurvedBdry( f )

            end select
!            call DGSpatial_QDotFaceLoop( f ) 
         end do
!
!        ***********************
!        Scale with the jacobian
!        ***********************
!
         do eID = 1 , mesh % no_of_elements 
            do eq = 1 , NCONS
               mesh % elements(eID) % QDot(:,:,eq) = mesh % elements(eID) % QDot(:,:,eq) / mesh % elements(eID) % Jac
            end do
         end do

      end subroutine DGSpatial_computeQDot
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              TIME DERIVATIVE VOLUME AND FACE LOOPS
!              -------------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_QDotVolumeLoop( e ) 
         use QuadElementClass
         implicit none
         class(QuadElement_t)          :: e 
         real ( kind=RP )       :: Fi ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real ( kind=RP )       :: Fv ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real ( kind=RP )       :: Fa ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real ( kind=RP )       :: F  ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real(kind=RP), pointer :: QDot(:,:,:)

         Fi = InviscidMethod % ComputeInnerFluxes( e )
#ifdef NAVIER_STOKES
         e % mu_a = ArtificialDissipation % ComputeElementViscosity ( e )

         Fv = ViscousMethod  % ComputeInnerFluxes( e )
         Fa = ArtificialDissipation % ComputeVolumeFluxes( e )

         F = Fi - Fv - Fa
#else
         F = Fi
#endif

         QDot(0:,0:,1:) => e % QDot
         QDot = ScalarWeakIntegrals % StdVolumeGreen(e , F)
  
      end subroutine DGSpatial_QDotVolumeLoop

      subroutine DGSpatial_QDotFaceLoop_Interior( ed )
         use QuadElementClass
         implicit none
         type(Edge_t)               :: ed
         real ( kind=RP )           :: FiL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FiR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FvL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FvR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: GL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS , 1 : NDIM)
         real ( kind=RP )           :: GR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS , 1 : NDIM)
         real ( kind=RP )           :: FL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real ( kind=RP ), pointer  :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FR) and the Gradient Riemann Solver (GL,GR)
!        --------------------------------------------------------------------------
         call ComputeRiemannSolver_InteriorEdge( ed , FL , FR , GL , GR)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(0:,0:,1:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT element
!        -----------------------------------------
         QDot(0:,0:,1:) => ed % quads(RIGHT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT , FR ) 

#ifdef NAVIER_STOKES
!
!        Add the gradient fluxes term
!        ----------------------------
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Add the contribution to the LEFT element
!           ----------------------------------------
            QDot(0:,0:,1:) => ed % quads(LEFT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , LEFT  , GL ) 
!
!           Add the contribution to the RIGHT element
!           -----------------------------------------
            QDot(0:,0:,1:) => ed % quads(RIGHT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT , GR ) 
!
         end if
#endif
           

      end subroutine DGSpatial_QDotFaceLoop_Interior

      subroutine DGSpatial_QDotFaceLoop_StraightBdry( ed )
         use QuadElementClass
         implicit none
         type(StraightBdryEdge_t)  :: ed
         real(kind=RP)             :: G ( 0 : ed % spA % N , 1:NCONS , 1:NDIM )  
         real(kind=RP)             :: F ( 0 : ed % spA % N , 1:NCONS)
         real ( kind=RP ), pointer :: QDot(:,:,:)

         call ComputeRiemannSolver_StraightBdryEdge( ed , F , G )

         QDot(0:,0:,1:) => ed % quads(1) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace( ed , 1 , F )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , 1  , G ) 

         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_StraightBdry

      subroutine DGSpatial_QDotFaceLoop_CurvedBdry( ed )
         use QuadElementClass
         implicit none
         type(CurvedBdryEdge_t)    :: ed
         real(kind=RP)             :: G ( 0 : ed % spA % N , 1:NCONS , 1:NDIM)
         real(kind=RP)             :: F ( 0 : ed % spA % N , 1:NCONS)
         real ( kind=RP ), pointer :: QDot(:,:,:)

         call ComputeRiemannSolver_CurvedBdryEdge( ed , F , G )

         QDot(0:,0:,1:) => ed % quads(1) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace( ed , 1 , F )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , 1  , G ) 

         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_CurvedBdry
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              COMPUTE INTERIOR FLUXES SUBROUTINES
!              -----------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeRiemannSolver_InteriorEdge( ed , FL , FR , GL , GR)
!
!        *****************************************************************************
!           This routine computes the Viscous Riemann problem in the "ed" edge.
!              -> The result is FStarL and FStarR

!        >> It considers several possibilities:
!
!              1/ LEFT edge needs a p-Transformation and RIGHT edge is reversed.
!              2/ LEFT edge needs a p-Transformation
!              3/ RIGHT edge needs a p-Transformation and RIGHT edge is reversed.
!              4/ RIGHT edge needs a p-Transformation
!              5/ RIGHT edge is reversed.
!              6/ No transformations are needed.
!
!           The Riemann solver is computed with the 
!>                self % RiemannSolver( N , QL , QR , dQL , dQR , normal )
!           procedure.
!
!        *****************************************************************************
!
         use QuadElementClass
         use MatrixOperations
         implicit none
         type(Edge_t)            :: ed
         real(kind=RP)           :: FL( 0 : ed % storage(LEFT ) % spA % N , NCONS )
         real(kind=RP)           :: FR( 0 : ed % storage(RIGHT) % spA % N , NCONS )
         real(kind=RP), intent (out) :: GL( 0 : ed % storage(LEFT ) % spA % N , 1:NCONS , 1:NDIM)
         real(kind=RP), intent (out) :: GR( 0 : ed % storage(RIGHT) % spA % N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real ( kind=RP)            :: Fi( 0 : ed % spA % N , 1 : NCONS )
         real ( kind=RP ), target   :: QL ( 0 : ed % spA % N , 1 : NCONS )
         real ( kind=RP ), target   :: QR ( 0 : ed % spA % N , 1 : NCONS )
#ifdef NAVIER_STOKES
         real ( kind=RP)            :: Fv( 0 : ed % spA % N , 1 : NCONS )
         real ( kind=RP)            :: Fa( 0 : ed % spA % N , 1 : NCONS )
         real ( kind=RP ), target   :: dQL ( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
         real ( kind=RP ), target   :: dQR ( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
         real ( kind=RP )           :: GauxL( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM)
         real ( kind=RP )           :: GauxR( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM)
#endif
         real ( kind=RP )           :: normal(NDIM , 0 : ed % spA % N )
         integer                    :: eq , iDim
!
!        Compute the normal
!        ------------------
         normal = spread( ed % n(IX:IY,0) , ncopies = ed % spA % N + 1 , dim = 2 )
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif

         if ( ed % transform(LEFT) .and. ed % inverse ) then
! 
!        -------------------------------------------------------
!>       First case: LEFT needs p-Adaption and RIGHT is reversed
!        -------------------------------------------------------
!
!           Transform the LEFT edge
!           -----------------------            
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
#ifdef NAVIER_STOKES
            do eq = 1 , NCONS
               call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % dQ(:,:,eq) , dQL(:,:,eq) )
            end do
#endif
!
!           Invert the RIGHT edge 
!           ---------------------
            QR  = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS )
#ifdef NAVIER_STOKES
            dQR = ed % storage(RIGHT) % dQ(ed % spA % N : 0 : -1 , 1:NDIM , 1:NCONS )
#endif
!
!           Compute the Riemann solver
!           --------------------------
            Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
            Fv = ViscousMethod % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )

            FR = ( Fi - Fv - Fa) * ed % dS(0)
#else
            FR = Fi * ed % dS(0)
#endif
!
!           Transform the LEFT edge
!           -----------------------
            call Mat_x_Mat( ed % T_backward , FR , FL )
!
!           Invert the RIGHT edge
!           ---------------------
            FR(0:ed % spA % N , 1:NCONS) = FR(ed % spA % N : 0 : -1 , 1:NCONS)

#ifdef NAVIER_STOKES
!
!           Compute the Gradient Riemann solver if proceeds
!           -----------------------------------------------
            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GauxL , GauxR ) 

               do iDim = 1 , NDIM
                  call Mat_x_Mat( ed % T_backward , GauxL(:,:,iDim) , GL(:,:,iDim) )
               end do

               GL = GL * ed % dS(0)
               GR = GauxR(ed % spA % N : 0 : -1 , 1:NCONS , 1:NDIM)  * ed % dS(0)

            end if
#endif

         elseif ( ed % transform(LEFT) ) then
! 
!        ----------------------------------
!>       Second case: LEFT needs p-Adaption.
!        ----------------------------------
!
!           Transform the LEFT edge
!           -----------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
#ifdef NAVIER_STOKES
            do eq = 1 , NCONS
               call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % dQ(:,:,eq) , dQL(:,:,eq) )
            end do
#endif
!
!           Get the RIGHT edge state
!           ------------------------
            QR = ed % storage(RIGHT) % Q
#ifdef NAVIER_STOKES
            dQR = ed % storage(RIGHT) % dQ
#endif
!
!           Compute the Riemann solver
!           --------------------------
            Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
            Fv = ViscousMethod % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )
            FR = ( Fi - Fv - Fa) * ed % dS(0)
#else
            FR = Fi * ed % dS(0)
#endif

!
!           Transform the LEFT edge
!           -----------------------
            call Mat_x_Mat( ed % T_backward , FR , FL )

#ifdef NAVIER_STOKES
!
!           Compute the Gradient Riemann solver if proceeds
!           -----------------------------------------------
            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GauxL , GR )

               do iDim = 1 , NDIM
                  call Mat_x_Mat( ed % T_backward , GauxL(:,:,iDim) , GL(:,:,iDim) )
               end do

               GL = GL * ed % dS(0)
               GR = GR * ed % dS(0)
            end if

#endif

         elseif ( ed % transform(RIGHT) .and.  ed % inverse ) then
! 
!        ----------------------------------------------------
!>       Third case: RIGHT needs p-Adaption, and its reversed.
!        ----------------------------------------------------
!
!           Get the LEFT edge 
!           -----------------
            QL = ed % storage(LEFT) % Q
#ifdef NAVIER_STOKES
            dQL = ed % storage(LEFT) % dQ
#endif
!
!           Transform the RIGHT edge 
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
#ifdef NAVIER_STOKES
            do eq = 1 , NCONS
               call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % dQ(:,:,eq) , dQR(:,:,eq) )
            end do
#endif
!
!           Invert the RIGHT edge 
!           ---------------------
            QR(0:ed % spA % N,1:NCONS) = QR(ed % spA % N : 0 : -1 , 1:NCONS)
#ifdef NAVIER_STOKES
            dQR = dQR(ed % spA % N : 0 : -1 , 1:NDIM , 1:NCONS )
#endif
!
!           Compute the Riemann solver
!           --------------------------
            Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
            Fv = ViscousMethod  % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )
            FL = ( Fi - Fv - Fa) * ed % dS(0)
#else
            FL = Fi * ed % dS(0)
#endif
!
!           Undo the transformation for the RIGHT edge
!           ------------------------------------------ 
            call Mat_x_Mat( ed % T_backward , FL , FR )
!
!           Invert the edge 
!           ---------------
            FR(0:ed % NLow,1:NCONS) = FR(ed % NLow:0:-1 , 1:NCONS)

#ifdef NAVIER_STOKES
!
!           Compute the Gradient Riemann solver if proceeds
!           -----------------------------------------------
            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GL , GauxR ) 

               do iDim = 1 , NDIM
                  call Mat_x_Mat( ed % T_backward , GauxR(:,:,iDim) , GR(:,:,iDim) )
               end do

               GL   = GL * ed % dS(0)
               GR   = GR(ed % spA % N : 0 : -1 , 1:NCONS, 1:NDIM) * ed % dS(0)
            end if
#endif

         elseif ( ed % transform(RIGHT) ) then
! 
!        -----------------------------------
!>       Fourth case: RIGHT needs p-Adaption.
!        -----------------------------------
!
!           Get the LEFT edge
!           -----------------
            QL = ed % storage(LEFT) % Q
#ifdef NAVIER_STOKES
            dQL = ed % storage(LEFT) % dQ
#endif
!
!           Transform the RIGHT edge
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
#ifdef NAVIER_STOKES
            do eq = 1 , NCONS
               call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % dQ(:,:,eq) , dQR(:,:,eq) )
            end do
#endif
!
!           Compute the Riemann solver
!           --------------------------
            Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
            Fv = ViscousMethod % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )
            FL = ( Fi - Fv - Fa) * ed % dS(0)
#else
            FL = Fi * ed % dS(0)
#endif
!
!           Undo the transformation for the RIGHT 
!           ------------------------------------------
            call Mat_x_Mat( ed % T_backward , FL , FR )

#ifdef NAVIER_STOKES
!
!           Compute the Gradient Riemann solver if proceeds
!           -----------------------------------------------
            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GL , GauxR ) 
      
               do iDim = 1 , NDIM
                  call Mat_x_Mat( ed % T_backward , GauxR(:,:,iDim) , GR(:,:,iDim) )
               end do

               GL = GL * ed % dS(0)
               GR = GR * ed % dS(0)
            end if
#endif
         
         elseif ( ed % inverse ) then
! 
!        -----------------------------
!>       Fifth case: RIGHT is reversed.
!        -----------------------------
!
!           Get the LEFT edge
!           -----------------
            QL = ed % storage(LEFT) % Q
#ifdef NAVIER_STOKES
            dQL = ed % storage(LEFT) % dQ
#endif
!
!           Invert the RIGHT edge
!           ---------------------
            QR(0:ed % spA % N , 1:NCONS ) = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS)
#ifdef NAVIER_STOKES
            dQR = ed % storage(RIGHT) % dQ(ed % spA % N : 0 : -1 , 1:NDIM , 1:NCONS )
#endif
!
!           Compute the Riemann solver
!           --------------------------
            Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
            Fv = ViscousMethod  % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )

            FL = ( Fi - Fv - Fa) * ed % dS(0)
#else
            FL = Fi * ed % dS(0)
#endif
!
!           Invert the RIGHT edge
!           ---------------------
            FR( 0 : ed % spA % N , 1:NCONS ) = FL ( ed % spA % N : 0 : -1 , 1:NCONS )

#ifdef NAVIER_STOKES
!
!           Compute the Gradient Riemann solver if proceeds
!           -----------------------------------------------
            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GL , GR ) 

               GL = GL * ed % dS(0)
               GR = GR(ed % spA % N : 0 : -1 , 1:NCONS , 1:NDIM) * ed % dS(0)
            end if
#endif

         else
! 
!        -----------------------------------------------------------------------
!>       Sixth case: Default case: neither p-Adaption nor inversion are required.
!        -----------------------------------------------------------------------
!
!           Get the LEFT edge
!           -----------------   
            QL = ed % storage(LEFT) % Q
#ifdef NAVIER_STOKES
            dQL = ed % storage(LEFT) % dQ
#endif 
!
!           Get the RIGHT edge
!           ------------------
            QR = ed % storage(RIGHT) % Q
#ifdef NAVIER_STOKES
            dQR = ed % storage(RIGHT) % dQ
#endif
!
!           Compute the Riemann solver
!           --------------------------
            Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
            Fv = ViscousMethod  % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )

            FL = ( Fi - Fv - Fa ) * ed % dS(0)
            FR = FL
#else
            FL = Fi * ed % dS(0)
            FR = FL
#endif

#ifdef NAVIER_STOKES
!
!           Compute the Gradient Riemann solver if proceeds
!           -----------------------------------------------
            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GL , GR ) 

               GL = GL * ed % dS(0)
               GR = GR * ed % dS(0)
            end if
#endif
         
         end if
!
!        Change the sign of FR so that it points towards the element outside
!        -------------------------------------------------------------------
         FR = -1.0_RP * FR

      end subroutine ComputeRiemannSolver_InteriorEdge

      subroutine ComputeRiemannSolver_StraightBdryEdge( ed , F , G )
         use QuadElementClass
         implicit none
         type(StraightBdryEdge_t)      :: ed
         real(kind=RP)                 :: F( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: G( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                 :: Fi( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: Q (0 : ed % spA % N , 1:NCONS) 
         real(kind=RP)                 :: Qb(0 : ed % spA % N , 1:NCONS)
#ifdef NAVIER_STOKES
         real(kind=RP)                 :: Fv( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: Fa( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: Gv( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                 :: Gaux( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                 :: dQ (0 : ed % spA % N , 1 : NDIM , 1 : NCONS)
         real(kind=RP)                 :: dQb(0 : ed % spA % N , 1 : NDIM , 1 : NCONS)
         real(kind=RP)                 :: Q1D  ( 1 : NCONS )
         real(kind=RP)                 :: Qb1D ( 1 : NCONS )
         real(kind=RP)                 :: dQ1D ( 1 : NDIM , 1 : NCONS ) 
         real(kind=RP)                 :: dQb1D( 1 : NDIM , 1 : NCONS ) 
#endif
         real(kind=RP)                 :: normal( NDIM , 0 : ed % spA % N )
         integer, pointer              :: N
         integer                       :: iXi
         procedure(RiemannSolverFunction), pointer    :: RiemannSolver

         N => ed % spA % N

         normal = spread( ed % n(IX:IY,0) , ncopies = ed % spA % N + 1 , dim = 2 )
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif

!
!        ===============
!>       INVISCID FLUXES
!        ===============
!
!        Select boundary condition type
!        ------------------------------
         select case ( ed % inviscidBCType ) 

            case ( WEAK_PRESCRIBED )
!
!              Prescribed boundary conditions: just grab the value
!              ---------------------------------------------------
               Fi = ed % FB

            case ( WEAK_RIEMANN )
!
!              Weak boundary conditions
!              ------------------------
               if ( associated ( ed % RiemannSolver ) ) then
                  RiemannSolver => ed % RiemannSolver

               else
                  RiemannSolver => InviscidMethod % RiemannSolver

               end if
!
!              Select LEFT and RIGHT states
!              ----------------------------
               if ( .not. ed % inverse ) then       ! This just occurs in periodic BCs
                  Q  = ed % storage(1) % Q
                  Qb = ed % uB
                  Fi = RiemannSolver ( ed % spA % N , Q , Qb , normal ) * ed % dS(0)

               else
                  Q  = ed % storage(1) % Q
                  Qb = ed % uB( ed % spA % N : 0 : -1 , 1:NCONS )
                  Fi = RiemannSolver ( ed % spA % N , Q , Qb , normal ) * ed % dS(0)

               end if
         end select

#ifdef NAVIER_STOKES
!
!        ==============
!>       VISCOUS FLUXES
!        ==============
!
!>       Select the appropriate boundary condition
!        -----------------------------------------
!

         if ( ed % viscousBCType(0) .eq. PERIODIC ) then
!
!           Periodic boundary conditions
!           ----------------------------
            if ( ed % inverse ) then
               Q  = ed % storage(1) % Q
               dQ = ed % storage(1) % dQ

               Qb  = ed % uB(N:0:-1,1:NCONS)
               dQb = ed % gB(N:0:-1,1:NDIM,1:NCONS)

               normal = spread( ed % n(IX:IY,0) , ncopies = N+1 , dim = 2 ) 
               Fv = ViscousMethod % RiemannSolver( ed , N , ed % invh , Q , Qb , dQ , dQb , normal ) * ed % dS(0)
               Fa = ArtificialDissipation % ComputeFaceFluxes( ed , Q , Qb , dQ , dQb , normal ) * ed % dS(0)

               if ( ViscousMethod % computeRiemannGradientFluxes ) then
                  call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , Q , Qb , normal , Gv , Gaux ) 
                  Gv = Gv * ed % dS(0)

               end if
            else

               Q  = ed % storage(1) % Q
               dQ = ed % storage(1) % dQ

               Qb  = ed % uB
               dQb = ed % gB

               normal = spread( ed % n(IX:IY,0) , ncopies = N+1 , dim = 2 ) 
               Fv = ViscousMethod % RiemannSolver( ed , N , ed % invh , Q , Qb , dQ , dQb , normal ) * ed % dS(0)
               Fa = ArtificialDissipation % ComputeFaceFluxes( ed , Q , Qb , dQ , dQb , normal ) * ed % dS(0)

               if ( ViscousMethod % computeRiemannGradientFluxes ) then
                   call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , Q , Qb , normal , Gv , Gaux )
                   Gv = Gv * ed % dS(0)

               end if

            end if

         elseif ( ed % viscousBCType(0) .eq. ADIABATIC ) then
!
!           Adiabatic Dirichlet boundary conditions
!           ---------------------------------------
            Q  = ed % storage(1) % Q
            dQ = ed % storage(1) % dQ
            Qb = ed % uSB

            normal = spread( ed % n(IX:IY,0) , ncopies = N+1 , dim = 2 )
            
            Fv = ViscousMethod % RiemannSolver_Adiabatic( ed , N , ed % invh , Q , dQ , Qb , normal ) * ed % dS(0) 
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , Q , Qb , dQ , dQ , normal ) * ed % dS(0)

            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               Gv = ViscousMethod % GradientRiemannSolver_Adiabatic( ed , N , Q , Qb , normal ) * ed % dS(0)

            end if

         else
!
!           Dirichlet/Neumann boundary conditions
!           -------------------------------------
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % storage(1) % Q , ed % uSB , ed % storage(1) % dQ , ed % storage(1) % dQ , normal ) * ed % dS(0)
            do iXi = 0 , ed % spA % N
               select case ( ed % viscousBCType(iXi) )
   
                  case ( DIRICHLET )
   
                     Q1D  = ed % storage(1) % Q(iXi,:)
                     dQ1D = ed % storage(1) % dQ(iXi,:,:)
                     Qb1D  = ed % uSB(iXi,:)
                     Fv(iXi,:) = ViscousMethod % RiemannSolver_Dirichlet ( ed , ed % spA % N , ed % invh , Q1D , dQ1D , Qb1D , ed % n(IX:IY,0) ) * ed % dS(0)

                     if ( ViscousMethod % computeRiemannGradientFluxes ) then
                        Gv(iXi,:,:) = ViscousMethod % GradientRiemannSolver_BoundaryCondition( ed , Q1D , Qb1D , ed % n(IX:IY,0) ) * ed % dS(0)
                     end if

                  case ( NEUMANN )

                     Fv(iXi,:)   = 0.0_RP
                     Fa(iXi,:)   = 0.0_RP
                     Gv(iXi,:,:) = 0.0_RP
         
                end select
             end do


         end if
#endif

#ifdef NAVIER_STOKES
         !Fa = 0.0_RP
         F = Fi - Fv - Fa
         G = Gv
#else
         F = Fi
         G = 0.0_RP
#endif

      end subroutine ComputeRiemannSolver_StraightBdryEdge

      subroutine ComputeRiemannSolver_CurvedBdryEdge( ed , F , G )
         use QuadElementClass
         implicit none
         type(CurvedBdryEdge_t)      :: ed
         real(kind=RP)                 :: F( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: G( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                 :: Fi( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: Q(0 : ed % spA % N , 1:NCONS) 
         real(kind=RP)                 :: Qb(0 : ed % spA % N , 1:NCONS)
#ifdef NAVIER_STOKES
         real(kind=RP)                 :: Fv( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: Fa( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: Gv( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                 :: Gaux( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                 :: dQ (0 : ed % spA % N , 1 : NDIM , 1 : NCONS)
         real(kind=RP)                 :: dQb(0 : ed % spA % N , 1 : NDIM , 1 : NCONS)
         real(kind=RP)                 :: Q1D  ( 1 : NCONS )
         real(kind=RP)                 :: Qb1D ( 1 : NCONS )
         real(kind=RP)                 :: dQ1D ( 1 : NDIM , 1 : NCONS ) 
         real(kind=RP)                 :: dQb1D( 1 : NDIM , 1 : NCONS ) 
#endif
         integer, pointer              :: N
         integer                       :: iXi
         procedure(RiemannSolverFunction), pointer    :: RiemannSolver

         N => ed % spA % N
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif
!
!        ===============
!>       INVISCID FLUXES
!        ===============
!
!        Select boundary condition type
!        ------------------------------
         select case ( ed % inviscidBCType ) 

            case ( WEAK_PRESCRIBED )
!
!              Prescribed boundary conditions: just grab the value
!              ---------------------------------------------------
               Fi = ed % FB

            case ( WEAK_RIEMANN )
!
!              Weak boundary conditions
!              ------------------------
               if ( associated ( ed % RiemannSolver ) ) then
                  RiemannSolver => ed % RiemannSolver

               else
                  RiemannSolver => InviscidMethod % RiemannSolver

               end if
!
!              Select LEFT and RIGHT states
!              ----------------------------
               if ( .not. ed % inverse ) then       ! This just occurs in periodic BCs
                  Q = ed % storage(1) % Q
                  Qb = ed % uB
                  Fi = RiemannSolver ( ed % spA % N , Q , Qb , ed % n ) * spread( ed % dS , ncopies = NCONS , dim = 2 )

               else
                  Q = ed % storage(1) % Q
                  Qb = ed % uB( ed % spA % N : 0 : -1 , 1:NCONS )
                  Fi = RiemannSolver ( ed % spA % N , Q , Qb , ed % n ) * spread( ed % dS , ncopies = NCONS , dim = 2 )

               end if
         end select

#ifdef NAVIER_STOKES
!
!        ==============
!>       VISCOUS FLUXES
!        ==============
!
!>       Select the appropriate boundary condition
!        -----------------------------------------
!
         if ( ed % viscousBCType(0) .eq. ADIABATIC ) then
!
!           Adiabatic Dirichlet boundary conditions
!           ---------------------------------------
            Q  = ed % storage(1) % Q
            dQ = ed % storage(1) % dQ
            Qb = ed % uSB

            Fv = ViscousMethod % RiemannSolver_Adiabatic( ed , N , ed % invh , Q , dQ , Qb , ed % n ) * spread ( ed % dS , ncopies = NCONS , dim = 2 )  
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , Q , Qb , dQ , dQ , ed % n ) * spread( ed % dS , ncopies = NCONS , dim = 2 ) 

            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               Gv = ViscousMethod % GradientRiemannSolver_Adiabatic( ed , N , Q , Qb , ed % n ) * spread( spread ( ed % dS , ncopies = NCONS , dim = 2 ) , ncopies = NDIM , dim = 3)

            end if

         else
!
!           Dirichlet/Neumann boundary conditions
!           -------------------------------------
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % storage(1) % Q , ed % uSB , ed % storage(1) % dQ , ed % storage(1) % dQ , ed % n ) * spread( ed % dS , ncopies = NCONS , dim = 2 ) 
            do iXi = 0 , ed % spA % N
               select case ( ed % viscousBCType(iXi) )
   
                  case ( DIRICHLET )
   
                     Q1D  = ed % storage(1) % Q(iXi,:)
                     dQ1D = ed % storage(1) % dQ(iXi,:,:)
                     Qb1D  = ed % uSB(iXi,:)
                     Fv(iXi,:) = ViscousMethod % RiemannSolver_Dirichlet ( ed , ed % spA % N , ed % invh , Q1D , dQ1D , Qb1D , ed % n(IX:IY,iXi) ) * ed % dS(iXi)

                     if ( ViscousMethod % computeRiemannGradientFluxes ) then
                        Gv(iXi,:,:) = ViscousMethod % GradientRiemannSolver_BoundaryCondition( ed , Q1D , Qb1D , ed % n(IX:IY,iXi) ) * ed % dS(iXi)
                     end if

                  case ( NEUMANN )

                     Fv(iXi,:)   = 0.0_RP
                     Fa(iXi,:)   = 0.0_RP
                     Gv(iXi,:,:) = 0.0_RP
         
                end select


             end do

         end if
#endif

#ifdef NAVIER_STOKES
         !Fa = 0.0_RP
         F = Fi - Fv - Fa
         G = Gv
#else
         F = Fi
         G = 0.0_RP
#endif

      end subroutine ComputeRiemannSolver_CurvedBdryEdge

!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              INTERPOLATE TO BOUNDARIES SUBROUTINES
!              -------------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_interpolateSolutionToBoundaries( mesh )
         use Physics
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(QuadMesh_t) :: mesh
!        --------------------------------------------------------------------
         integer                       :: eID , eq
         class(QuadElement_t), pointer :: e
         class(Edge_t), pointer        :: ed
         integer, pointer              :: N

         do eID = 1 , mesh % no_of_elements

            e => mesh % elements(eID) 
            N => e % spA % N
!
!           Prolong the BOTTOM edge
!           -----------------------   
            ed => e % edges(EBOTTOM) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(EBOTTOM)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,LEFT) , N+1 )
            end do
!
!           Prolong the RIGHT edge. TODOO: implement MatrixTimesVectorInIndex to avoid transposes.
!           -----------------------   
            ed => e % edges(ERIGHT) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(ERIGHT)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,RIGHT) , N+1 , trA = .true.)
            end do
!
!           Prolong the TOP edge
!           -----------------------   
            ed => e % edges(ETOP) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(ETOP)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,RIGHT) , N+1 )
            end do
!
!           Prolong the LEFT edge. TODOO: implement MatrixTimesVectorInIndex to avoid transposes.
!           -----------------------   
            ed => e % edges(ELEFT) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(ELEFT)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,LEFT) , N+1 , trA = .true.)
            end do
        
         end do
            
      end subroutine DGSpatial_interpolateSolutionToBoundaries

#ifdef NAVIER_STOKES
      subroutine DGSpatial_interpolateGradientsToBoundaries( mesh )
         use Physics
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(QuadMesh_t) :: mesh
!        --------------------------------------------------------------------
         integer                       :: eID , eq , iDim
         class(QuadElement_t), pointer :: e
         class(Edge_t), pointer        :: ed
         integer, pointer              :: N

         do eID = 1 , mesh % no_of_elements

            e => mesh % elements(eID) 
            N => e % spA % N
!
!           Prolong the BOTTOM edge
!           -----------------------   
            ed => e % edges(EBOTTOM) % f
            do eq = 1 , NCONS ; do iDim = 1 , NDIM
               ed % storage(e % quadPosition(EBOTTOM)) % dQ(0:N,iDim,eq) = MatrixTimesVector_F( e % dQ(0:N,0:N,iDim,eq) , e % spA % lb(0:N,LEFT) , N+1 )
            end do            ; end do
!
!           Prolong the RIGHT edge. TODOO: implement MatrixTimesVectorInIndex to avoid transposes.
!           -----------------------   
            ed => e % edges(ERIGHT) % f
            do eq = 1 , NCONS ; do iDim = 1 , NDIM
               ed % storage(e % quadPosition(ERIGHT)) % dQ(0:N,iDim,eq) = MatrixTimesVector_F( e % dQ(0:N,0:N,iDim,eq) , e % spA % lb(0:N,RIGHT) , N+1 , trA = .true.)
            end do            ; end do
!
!           Prolong the TOP edge
!           -----------------------   
            ed => e % edges(ETOP) % f
            do eq = 1 , NCONS ; do iDim = 1 , NDIM
               ed % storage(e % quadPosition(ETOP)) % dQ(0:N,iDim,eq) = MatrixTimesVector_F( e % dQ(0:N,0:N,iDim,eq) , e % spA % lb(0:N,RIGHT) , N+1 )
            end do            ; end do
!
!           Prolong the LEFT edge. TODOO: implement MatrixTimesVectorInIndex to avoid transposes.
!           -----------------------   
            ed => e % edges(ELEFT) % f
            do eq = 1 , NCONS ; do iDim = 1 , NDIM
               ed % storage(e % quadPosition(ELEFT)) % dQ(0:N,iDim,eq) = MatrixTimesVector_F( e % dQ(0:N,0:N,iDim,eq) , e % spA % lb(0:N,LEFT) , N+1 , trA = .true.)
            end do            ; end do
        
         end do
 
      end subroutine DGSpatial_interpolateGradientsToBoundaries
#endif

end module DGSpatialDiscretizationMethods
