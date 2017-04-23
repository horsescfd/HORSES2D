module DGSpatialDiscretizationMethods
   use SMConstants
   use Physics
   use QuadMeshClass
   use DGInviscidMethods
   use DGWeakIntegrals
#ifdef NAVIER_STOKES
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
   class(ViscousMethod_t),  pointer :: ViscousMethod
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
         ViscousMethod  => ViscousMethod_Initialization()
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
!        Reset QDot
!        ----------
         call DGSpatial_resetQDot( mesh )
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
!
      subroutine DGSpatial_resetQDot( mesh )
         implicit none
         class(QuadMesh_t)         :: mesh
!        --------------------------------------
         integer                   :: eID

         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % QDot = 0.0_RP
         end do

      end subroutine DGSpatial_resetQDot

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
         real ( kind=RP )       :: F  ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real(kind=RP), pointer :: QDot(:,:,:)

         Fi = InviscidMethod % ComputeInnerFluxes( e )
#ifdef NAVIER_STOKES
         Fv = ViscousMethod  % ComputeInnerFluxes( e )

         F = Fi - Fv
#else
         F = Fi
#endif

         QDot(0:,0:,1:) => e % QDot
         QDot = QDot + ScalarWeakIntegrals % StdVolumeGreen(e , F)
  
      end subroutine DGSpatial_QDotVolumeLoop

      subroutine DGSpatial_QDotFaceLoop_Interior( ed )
         use QuadElementClass
         implicit none
         type(Edge_t)               :: ed
         real ( kind=RP )           :: FiL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FiR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FvL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FvR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: GvL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS , 1 : NDIM)
         real ( kind=RP )           :: GvR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS , 1 : NDIM)
         real ( kind=RP )           :: FL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real ( kind=RP ), pointer  :: QDot(:,:,:)

         call InviscidMethod % ComputeRiemannFluxes( ed , FiL , FiR )
#ifdef NAVIER_STOKES
         call ViscousMethod % ComputeRiemannFluxes( ed , FvL , FvR , GvL , GvR )
         FL = FiL - FvL
         FR = FiR - FvR
#else
         FL = FiL
         FR = FiR
#endif

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
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Add the contribution to the LEFT element
!           ----------------------------------------
            QDot(0:,0:,1:) => ed % quads(LEFT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , LEFT  , GvL ) 
!
!           Add the contribution to the RIGHT element
!           -----------------------------------------
            QDot(0:,0:,1:) => ed % quads(RIGHT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT , GvR ) 
!
         end if
#endif
           

      end subroutine DGSpatial_QDotFaceLoop_Interior

      subroutine DGSpatial_QDotFaceLoop_StraightBdry( ed )
         use QuadElementClass
         implicit none
         type(StraightBdryEdge_t)  :: ed
         real(kind=RP)             :: Fi ( 0 : ed % spA % N , 1:NCONS)
         real(kind=RP)             :: Fv ( 0 : ed % spA % N , 1:NCONS)
         real(kind=RP)             :: Gv ( 0 : ed % spA % N , 1:NCONS , 1:NDIM )  
         real(kind=RP)             :: F ( 0 : ed % spA % N , 1:NCONS)
         real ( kind=RP ), pointer :: QDot(:,:,:)

         call InviscidMethod % ComputeRiemannFluxes( ed , Fi )
#ifdef NAVIER_STOKES
         call ViscousMethod % ComputeRiemannFluxes( ed , Fv , Gv)
         F = Fi - Fv
#else
         F = Fi
#endif

         QDot(0:,0:,1:) => ed % quads(1) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace( ed , 1 , F )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , 1  , Gv ) 

         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_StraightBdry

      subroutine DGSpatial_QDotFaceLoop_CurvedBdry( ed )
         use QuadElementClass
         implicit none
         type(CurvedBdryEdge_t)    :: ed
         real(kind=RP)             :: Fi ( 0 : ed % spA % N , 1:NCONS)
         real(kind=RP)             :: Fv ( 0 : ed % spA % N , 1:NCONS)
         real(kind=RP)             :: Gv ( 0 : ed % spA % N , 1:NCONS , 1:NDIM)
         real(kind=RP)             :: F ( 0 : ed % spA % N , 1:NCONS)
         real ( kind=RP ), pointer :: QDot(:,:,:)

         call InviscidMethod % ComputeRiemannFluxes( ed , Fi )
#ifdef NAVIER_STOKES
         call ViscousMethod % ComputeRiemannFluxes( ed , Fv , Gv )
         F = Fi - Fv
#else
         F = Fi
#endif

         QDot(0:,0:,1:) => ed % quads(1) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace( ed , 1 , F )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , 1  , Gv ) 

         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_CurvedBdry
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
