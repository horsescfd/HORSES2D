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
      subroutine DGSpatial_computeTimeDerivative( mesh , time ) 
!
!        ***************************************************
!           Subroutine that performs the spatial
!         discretization and computes the time
!         derivative QDot.
!        ***************************************************
!
         implicit none
         class(QuadMesh_t)          :: mesh
         real(kind=RP), intent(in)  :: time
!
!        Prepare the mesh for a new iteration
!        ------------------------------------
         call DGSpatial_newTimeStep( mesh , time )
!
!        Compute QDot
!        ------------
         call DGSpatial_computeQDot( mesh )

      end subroutine DGSpatial_computeTimeDerivative
      
         subroutine DGSpatial_newTimeStep( mesh , time)
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
         real(kind=RP)             :: time
         integer                   :: zoneID
!
!        Interpolate solution to boundaries
!        ----------------------------------
         call DGSpatial_interpolateSolutionToBoundaries( mesh )
!
!        Update the zones solution
!        -------------------------
         do zoneID = 1 , size(mesh % zones) - 1
            call mesh % zones(zoneID) % UpdateSolution(time)
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

               type is (CurvedEdge_t)
                  call DGSpatial_QDotFaceLoop_CurvedInterior( f ) 

               type is (SubdividedEdge_t)
                  call DGSpatial_QDotFaceLoop_SubdividedEdge( f ) 

               type is (CurvedSubdividedEdge_t)
                  call DGSpatial_QDotFaceLoop_CurvedSubdividedEdge( f ) 

               type is (StraightBdryEdge_t)
                  call DGSpatial_QDotFaceLoop_StraightBdry( f )
  
               type is (CurvedBdryEdge_t)
                  call DGSpatial_QDotFaceLoop_CurvedBdry( f )

            end select
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
         real(kind=RP)       :: Fi ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real(kind=RP)       :: Fv ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real(kind=RP)       :: Fa ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
         real(kind=RP)       :: F  ( 0:e % spA % N , 0:e % spA % N , 1:NCONS   , 1:NDIM )
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
         type(Edge_t)           :: ed
         real(kind=RP)          :: FiL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real(kind=RP)          :: FiR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real(kind=RP)          :: FvL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real(kind=RP)          :: FvR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real(kind=RP)          :: GL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS , 1 : NDIM)
         real(kind=RP)          :: GR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS , 1 : NDIM)
         real(kind=RP)          :: FL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real(kind=RP)          :: FR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real(kind=RP), pointer :: QDot(:,:,:)
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

      subroutine DGSpatial_QDotFaceLoop_CurvedInterior( ed )
         use QuadElementClass
         implicit none
         type(CurvedEdge_t)         :: ed
         real(kind=RP)           :: FiL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real(kind=RP)           :: FiR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real(kind=RP)           :: FvL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real(kind=RP)           :: FvR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real(kind=RP)           :: GL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS , 1 : NDIM)
         real(kind=RP)           :: GR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS , 1 : NDIM)
         real(kind=RP)           :: FL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real(kind=RP)           :: FR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
         real(kind=RP), pointer  :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FR) and the Gradient Riemann Solver (GL,GR)
!        --------------------------------------------------------------------------
         call ComputeRiemannSolver_InteriorCurvedEdge( ed , FL , FR , GL , GR)
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
           

      end subroutine DGSpatial_QDotFaceLoop_CurvedInterior

      subroutine DGSpatial_QDotFaceLoop_SubdividedEdge( ed )
         use QuadElementClass
         implicit none
         type(SubdividedEdge_t) :: ed
         real(kind=RP)          :: FiL  ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FiRN ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FiRS ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FvL  ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FvRN ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FvRS ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: GL   ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)          :: GRN  ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)          :: GRS  ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)          :: FL   ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FRN  ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)          :: FRS  ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS            )
         real(kind=RP), pointer :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FRN,FRS) and the Gradient Riemann Solver (GL,GRN,GRS)
!        ------------------------------------------------------------------------------------
         call ComputeRiemannSolver_SubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(0:,0:,1:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT-NORTH element
!        -----------------------------------------------
         QDot(0:,0:,1:) => ed % quads(RIGHT_NORTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_NORTH , FRN ) 
!
!        Add the contribution to the RIGHT-SOUTH element
!        -----------------------------------------------
         QDot(0:,0:,1:) => ed % quads(RIGHT_SOUTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_SOUTH , FRS ) 

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
!           Add the contribution to the RIGHT-NORTH element
!           -----------------------------------------------
            QDot(0:,0:,1:) => ed % quads(RIGHT_NORTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_NORTH , GRN ) 
!
!           Add the contribution to the RIGHT-SOUTH element
!           -----------------------------------------------
            QDot(0:,0:,1:) => ed % quads(RIGHT_SOUTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_SOUTH , GRS ) 
!
         end if
#endif
           

      end subroutine DGSpatial_QDotFaceLoop_SubdividedEdge

      subroutine DGSpatial_QDotFaceLoop_CurvedSubdividedEdge( ed )
         use QuadElementClass
         implicit none
         type(CurvedSubdividedEdge_t) :: ed
         real(kind=RP)                :: FiL  ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FiRN ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FiRS ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FvL  ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FvRN ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FvRS ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: GL   ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                :: GRN  ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                :: GRS  ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS , 1 : NDIM )
         real(kind=RP)                :: FL   ( 0 : ed % storage ( LEFT        ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FRN  ( 0 : ed % storage ( RIGHT_NORTH ) % spA % N , 1 : NCONS            )
         real(kind=RP)                :: FRS  ( 0 : ed % storage ( RIGHT_SOUTH ) % spA % N , 1 : NCONS            )
         real(kind=RP), pointer       :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FRN,FRS) and the Gradient Riemann Solver (GL,GRN,GRS)
!        ------------------------------------------------------------------------------------
         call ComputeRiemannSolver_CurvedSubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(0:,0:,1:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT-NORTH element
!        -----------------------------------------------
         QDot(0:,0:,1:) => ed % quads(RIGHT_NORTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_NORTH , FRN ) 
!
!        Add the contribution to the RIGHT-SOUTH element
!        -----------------------------------------------
         QDot(0:,0:,1:) => ed % quads(RIGHT_SOUTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_SOUTH , FRS ) 

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
!           Add the contribution to the RIGHT-NORTH element
!           -----------------------------------------------
            QDot(0:,0:,1:) => ed % quads(RIGHT_NORTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_NORTH , GRN ) 
!
!           Add the contribution to the RIGHT-SOUTH element
!           -----------------------------------------------
            QDot(0:,0:,1:) => ed % quads(RIGHT_SOUTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_SOUTH , GRS ) 
!
         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_CurvedSubdividedEdge

      subroutine DGSpatial_QDotFaceLoop_StraightBdry( ed )
         use QuadElementClass
         implicit none
         type(StraightBdryEdge_t)  :: ed
         real(kind=RP)             :: G ( 0 : ed % spA % N , 1:NCONS , 1:NDIM )  
         real(kind=RP)             :: F ( 0 : ed % spA % N , 1:NCONS)
         real(kind=RP), pointer :: QDot(:,:,:)

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
         real(kind=RP), pointer :: QDot(:,:,:)

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
         type(Edge_t)                :: ed
         real(kind=RP)               :: FL( 0 : ed % storage(LEFT ) % spA % N , NCONS )
         real(kind=RP)               :: FR( 0 : ed % storage(RIGHT) % spA % N , NCONS )
         real(kind=RP), intent (out) :: GL( 0 : ed % storage(LEFT ) % spA % N , 1:NCONS , 1:NDIM)
         real(kind=RP), intent (out) :: GR( 0 : ed % storage(RIGHT) % spA % N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: Fi    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP)          :: Fstar ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) , target :: QL    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) , target :: QR    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) , target :: dQL   ( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
         real(kind=RP) , target :: dQR   ( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
#ifdef NAVIER_STOKES
         real(kind=RP) :: Fv    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) :: Fa    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) :: GauxL ( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP) :: GauxR ( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM ) 
#endif
         real(kind=RP)           :: normal(NDIM , 0 : ed % spA % N )
         integer                    :: eq , iDim , i 
!
!        Compute the normal
!        ------------------
         do i = 0 , ed % spA % N
            normal(IX:IY,i) = ed % n(IX:IY,0) 
         end do
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif
!
!        Get the solution projection onto the edge
!        -----------------------------------------
#ifdef NAVIER_STOKES
         call ed % ProjectSolutionAndGradient(ed , QL , QR , dQL , dQR )
#else
         call ed % ProjectSolution( ed , QL , QR )
#endif
!
!        Compute the inviscid Riemann solver
!        -----------------------------------
         Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , normal ) 
#ifdef NAVIER_STOKES
!
!        Compute the viscous Riemann solver
!        ----------------------------------
         Fv = ViscousMethod % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , normal )
!
!        Compute the artificial dissipation Riemann solver
!        -------------------------------------------------
         Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , normal )
!
!        The resulting flux is: FStar = ( Inviscid - Viscous - ArtificialDissipation ) dS
!        --------------------------------------------------------------------------------
         FStar = ( Fi - Fv - Fa) * ed % dS(0)
#else
         FStar = Fi * ed % dS(0)
#endif
!
!        Return the resulting Riemann flux to each element frame
!        -------------------------------------------------------
         call ed % ProjectFluxes( ed , FStar , FL , FR )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , normal , GauxL , GauxR ) 

               call ed % ProjectGradientFluxes( ed , GauxL , GauxR , GL , GR )

               GL = GL * ed % dS(0)
               GR = GR * ed % dS(0)

         end if
#endif

      end subroutine ComputeRiemannSolver_InteriorEdge

      subroutine ComputeRiemannSolver_InteriorCurvedEdge( ed , FL , FR , GL , GR)
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
         type(CurvedEdge_t)          :: ed
         real(kind=RP)               :: FL( 0 : ed % storage(LEFT ) % spA % N , NCONS )
         real(kind=RP)               :: FR( 0 : ed % storage(RIGHT) % spA % N , NCONS )
         real(kind=RP), intent (out) :: GL( 0 : ed % storage(LEFT ) % spA % N , 1:NCONS , 1:NDIM)
         real(kind=RP), intent (out) :: GR( 0 : ed % storage(RIGHT) % spA % N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: Fi    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP)          :: Fstar ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) , target :: QL    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) , target :: QR    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP) , target :: dQL   ( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
         real(kind=RP) , target :: dQR   ( 0 : ed % spA % N , 1 : NDIM , 1 : NCONS ) 
#ifdef NAVIER_STOKES
         real(kind=RP)  :: Fv    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP)  :: Fa    ( 0 : ed % spA % N , 1 : NCONS            ) 
         real(kind=RP)  :: GauxL ( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: GauxR ( 0 : ed % spA % N , 1 : NCONS , 1 : NDIM ) 
#endif
         integer :: eq , dimID , i 
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif
!
!        Get the solution projection onto the edge
!        -----------------------------------------
#ifdef NAVIER_STOKES
         call ed % ProjectSolutionAndGradient(ed , QL , QR , dQL , dQR )
#else
         call ed % ProjectSolution(ed , QL , QR )
#endif
!
!        Compute the inviscid Riemann solver
!        -----------------------------------
         Fi = InviscidMethod % RiemannSolver( ed % spA % N , QL , QR , ed % n ) 
#ifdef NAVIER_STOKES
!
!        Compute the viscous Riemann solver
!        ----------------------------------
         Fv = ViscousMethod % RiemannSolver( ed , ed % spA % N , ed % invh , QL , QR , dQL , dQR , ed % n )
!
!        Compute the artificial dissipation Riemann solver
!        -------------------------------------------------
         Fa = ArtificialDissipation % ComputeFaceFluxes( ed , QL , QR , dQL , dQR , ed % n )
#endif
!
!        The resulting flux is: FStar = ( Inviscid - Viscous - ArtificialDissipation ) dS
!        --------------------------------------------------------------------------------
         do eq = 1 , NCONS 
#ifdef NAVIER_STOKES
            FStar(:,eq) = ( Fi(:,eq) - Fv(:,eq) - Fa(:,eq) ) * ed % dS
#else
            FStar(:,eq) = Fi(:,eq) * ed % dS 
#endif
         end do
!
!        Return the resulting Riemann flux to each element frame
!        -------------------------------------------------------
         call ed % ProjectFluxes( ed , FStar , FL , FR )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , ed % n , GauxL , GauxR ) 

               call ed % ProjectGradientFluxes( ed , GauxL , GauxR , GL , GR )

               do dimID = 1 , NDIM    ; do eq = 1 , NCONS
                  GL(:,eq,dimID) = GL(:,eq,dimID) * ed % dS
                  GR(:,eq,dimID) = GR(:,eq,dimID) * ed % dS
               end do                 ; end do

         end if
#endif

      end subroutine ComputeRiemannSolver_InteriorCurvedEdge

      subroutine ComputeRiemannSolver_SubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS )
         use QuadElementClass
         use MatrixOperations
         implicit none
         type(SubdividedEdge_t)              :: ed
         real(kind=RP) ,         intent(out) :: FL  ( 0 : ed % storage ( LEFT        )  % spA % N , NCONS            )
         real(kind=RP) ,         intent(out) :: FRN ( 0 : ed % storage ( RIGHT_NORTH )  % spA % N , NCONS            )
         real(kind=RP) ,         intent(out) :: FRS ( 0 : ed % storage ( RIGHT_SOUTH )  % spA % N , NCONS            )
         real(kind=RP) ,         intent(out) :: GL  ( 0 : ed % storage ( LEFT        )  % spA % N , 1:NCONS , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRN ( 0 : ed % storage ( RIGHT_NORTH )  % spA % N , 1:NCONS , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRS ( 0 : ed % storage ( RIGHT_SOUTH )  % spA % N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: FiN    ( 0 : ed % spA_N % N , 1 : NCONS            )      ! -----------------------------
         real(kind=RP)          :: FstarN ( 0 : ed % spA_N % N , 1 : NCONS            )      !     Variables in the NORTH
         real(kind=RP) , target :: QLN    ( 0 : ed % spA_N % N , 1 : NCONS            )      !  mortar.
         real(kind=RP) , target :: QRN    ( 0 : ed % spA_N % N , 1 : NCONS            )      !  
         real(kind=RP) , target :: dQLN   ( 0 : ed % spA_N % N , 1 : NDIM , 1 : NCONS )      ! 
         real(kind=RP) , target :: dQRN   ( 0 : ed % spA_N % N , 1 : NDIM , 1 : NCONS )      ! -----------------------------
         real(kind=RP)          :: FiS    ( 0 : ed % spA_S % N , 1 : NCONS            )      ! -----------------------------
         real(kind=RP)          :: FstarS ( 0 : ed % spA_S % N , 1 : NCONS            )      !     Variables in the SOUTH
         real(kind=RP) , target :: QLS    ( 0 : ed % spA_S % N , 1 : NCONS            )      !  mortar.
         real(kind=RP) , target :: QRS    ( 0 : ed % spA_S % N , 1 : NCONS            )      !  
         real(kind=RP) , target :: dQLS   ( 0 : ed % spA_S % N , 1 : NDIM , 1 : NCONS )      ! 
         real(kind=RP) , target :: dQRS   ( 0 : ed % spA_S % N , 1 : NDIM , 1 : NCONS )      ! -----------------------------
#ifdef NAVIER_STOKES
         real(kind=RP)  :: FvN    ( 0 : ed % spA_N % N , 1 : NCONS            ) 
         real(kind=RP)  :: FaN    ( 0 : ed % spA_N % N , 1 : NCONS            ) 
         real(kind=RP)  :: GauxLN ( 0 : ed % spA_N % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: GauxRN ( 0 : ed % spA_N % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: FvS    ( 0 : ed % spA_S % N , 1 : NCONS            ) 
         real(kind=RP)  :: FaS    ( 0 : ed % spA_S % N , 1 : NCONS            ) 
         real(kind=RP)  :: GauxLS ( 0 : ed % spA_S % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: GauxRS ( 0 : ed % spA_S % N , 1 : NCONS , 1 : NDIM ) 
#endif
         real(kind=RP) :: normal_N(NDIM , 0 : ed % spA_N % N )
         real(kind=RP) :: normal_S(NDIM , 0 : ed % spA_S % N )
         integer       :: eq , iDim , i 
!
!        Compute the normal
!        ------------------
         do i = 0 , ed % spA_N % N
            normal_N(IX:IY,i) = ed % normal_N(IX:IY,0)
         end do

         do i = 0 , ed % spA_S % N
            normal_S(IX:IY,i) = ed % normal_S(IX:IY,0) 
         end do
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif
!
!        Get the solution projection onto the edge
!        -----------------------------------------
         call Mat_x_Mat ( ed % T_LN_FWD , ed % storage ( LEFT        )  % Q , QLN ) 
         call Mat_x_Mat ( ed % T_LS_FWD , ed % storage ( LEFT        )  % Q , QLS ) 
         call Mat_x_Mat ( ed % T_RN_FWD , ed % storage ( RIGHT_NORTH )  % Q , QRN ) 
         call Mat_x_Mat ( ed % T_RS_FWD , ed % storage ( RIGHT_SOUTH )  % Q , QRS ) 
!
!        Compute the inviscid Riemann solver
!        -----------------------------------
         FiN = InviscidMethod % RiemannSolver( ed % spA_N % N , QLN , QRN , normal_N ) 
         FiS = InviscidMethod % RiemannSolver( ed % spA_S % N , QLS , QRS , normal_S ) 
#ifdef NAVIER_STOKES
!
!        Get the gradient projection onto the edge
!        -----------------------------------------
         do eq = 1 , NCONS
            call Mat_x_Mat ( ed % T_LN_FWD , ed % storage ( LEFT        )  % dQ(:,:,eq) , dQLN(:,:,eq) ) 
            call Mat_x_Mat ( ed % T_LS_FWD , ed % storage ( LEFT        )  % dQ(:,:,eq) , dQLS(:,:,eq) ) 
            call Mat_x_Mat ( ed % T_RN_FWD , ed % storage ( RIGHT_NORTH )  % dQ(:,:,eq) , dQRN(:,:,eq) ) 
            call Mat_x_Mat ( ed % T_RS_FWD , ed % storage ( RIGHT_SOUTH )  % dQ(:,:,eq) , dQRS(:,:,eq) ) 
         end do
!
!        Compute the viscous Riemann solver
!        ----------------------------------
         FvN = ViscousMethod % RiemannSolver( ed , ed % spA_N % N , ed % invh , QLN , QRN , dQLN , dQRN , normal_N )
         FvS = ViscousMethod % RiemannSolver( ed , ed % spA_S % N , ed % invh , QLS , QRS , dQLS , dQRS , normal_S )
!
!        Compute the artificial dissipation Riemann solver
!        -------------------------------------------------
         FaN = ArtificialDissipation % ComputeFaceFluxes( ed , QLN , QRN , dQLN , dQRN , normal_N )
         FaS = ArtificialDissipation % ComputeFaceFluxes( ed , QLS , QRS , dQLS , dQRS , normal_S )
!
!        The resulting flux is: FStar = ( Inviscid - Viscous - ArtificialDissipation ) dS
!        --------------------------------------------------------------------------------
         FStarN = ( FiN - FvN - FaN ) * ed % dS_N(0)
         FStarS = ( FiS - FvS - FaS ) * ed % dS_S(0)
#else
         FStarN = FiN * ed % dS_N(0)
         FStarS = FiS * ed % dS_S(0)
#endif
!
!        Return the resulting Riemann flux to each element frame
!        -------------------------------------------------------
         FL =  Mat_x_Mat_F( ed % T_LN_BKW , FStarN , ed % spA % N + 1, NCONS ) + Mat_x_Mat_F( ed % T_LS_BKW , FStarS , ed % spA % N + 1, NCONS )
         call Mat_x_Mat( ed % T_RN_BKW , FStarN , FRN )
         call Mat_x_Mat( ed % T_RS_BKW , FStarS , FRS )
      
#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Compute the gradient Riemann solver
!           -----------------------------------
            call ViscousMethod % GradientRiemannSolver ( ed , ed % spA_N % N , QLN , QRN , normal_N , GauxLN , GauxRN ) 
            call ViscousMethod % GradientRiemannSolver ( ed , ed % spA_S % N , QLS , QRS , normal_S , GauxLS , GauxRS ) 
!
!           Scale it with the surface Jacobian
!           ----------------------------------
            GauxLN = GauxLN * ed % dS_N(0)
            GauxLS = GauxLS * ed % dS_S(0)
            GauxRN = GauxRN * ed % dS_N(0)
            GauxRS = GauxRS * ed % dS_S(0)
!
!           Return the resulting Riemann flux to each element frame
!           -------------------------------------------------------
            GL(:,:,IX) = Mat_x_Mat_F( ed % T_LN_BKW , GauxLN(:,:,IX) , ed % spA % N + 1 , NCONS ) + Mat_x_Mat_F( ed % T_LS_BKW , GauxLS(:,:,IX) , ed % spA % N + 1 , NCONS )
            GL(:,:,IY) = Mat_x_Mat_F( ed % T_LN_BKW , GauxLN(:,:,IY) , ed % spA % N + 1 , NCONS ) + Mat_x_Mat_F( ed % T_LS_BKW , GauxLS(:,:,IY) , ed % spA % N + 1 , NCONS )

            call Mat_x_Mat( ed % T_RN_BKW , GauxRN(:,:,IX) , GRN(:,:,IX) )
            call Mat_x_Mat( ed % T_RN_BKW , GauxRN(:,:,IY) , GRN(:,:,IY) )

            call Mat_x_Mat( ed % T_RS_BKW , GauxRS(:,:,IX) , GRS(:,:,IX) )
            call Mat_x_Mat( ed % T_RS_BKW , GauxRS(:,:,IY) , GRS(:,:,IY) )

         end if
#endif
!
!        Change RIGHT fluxes signs to account for the normal direction
!        ------------------------------------------------------------- 
         FRN = -FRN
         FRS = -FRS

      end subroutine ComputeRiemannSolver_SubdividedEdge

      subroutine ComputeRiemannSolver_CurvedSubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS )
         use QuadElementClass
         use MatrixOperations
         implicit none
         type(CurvedSubdividedEdge_t)        :: ed
         real(kind=RP) ,         intent(out) :: FL  ( 0 : ed % storage ( LEFT        )  % spA % N , NCONS            )
         real(kind=RP) ,         intent(out) :: FRN ( 0 : ed % storage ( RIGHT_NORTH )  % spA % N , NCONS            )
         real(kind=RP) ,         intent(out) :: FRS ( 0 : ed % storage ( RIGHT_SOUTH )  % spA % N , NCONS            )
         real(kind=RP) ,         intent(out) :: GL  ( 0 : ed % storage ( LEFT        )  % spA % N , 1:NCONS , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRN ( 0 : ed % storage ( RIGHT_NORTH )  % spA % N , 1:NCONS , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRS ( 0 : ed % storage ( RIGHT_SOUTH )  % spA % N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: FiN    ( 0 : ed % spA_N % N , 1 : NCONS            )      ! -----------------------------
         real(kind=RP)          :: FstarN ( 0 : ed % spA_N % N , 1 : NCONS            )      !     Variables in the NORTH
         real(kind=RP) , target :: QLN    ( 0 : ed % spA_N % N , 1 : NCONS            )      !  mortar.
         real(kind=RP) , target :: QRN    ( 0 : ed % spA_N % N , 1 : NCONS            )      !  
         real(kind=RP) , target :: dQLN   ( 0 : ed % spA_N % N , 1 : NDIM , 1 : NCONS )      ! 
         real(kind=RP) , target :: dQRN   ( 0 : ed % spA_N % N , 1 : NDIM , 1 : NCONS )      ! -----------------------------
         real(kind=RP)          :: FiS    ( 0 : ed % spA_S % N , 1 : NCONS            )      ! -----------------------------
         real(kind=RP)          :: FstarS ( 0 : ed % spA_S % N , 1 : NCONS            )      !     Variables in the SOUTH
         real(kind=RP) , target :: QLS    ( 0 : ed % spA_S % N , 1 : NCONS            )      !  mortar.
         real(kind=RP) , target :: QRS    ( 0 : ed % spA_S % N , 1 : NCONS            )      !  
         real(kind=RP) , target :: dQLS   ( 0 : ed % spA_S % N , 1 : NDIM , 1 : NCONS )      ! 
         real(kind=RP) , target :: dQRS   ( 0 : ed % spA_S % N , 1 : NDIM , 1 : NCONS )      ! -----------------------------
#ifdef NAVIER_STOKES
         real(kind=RP)  :: FvN    ( 0 : ed % spA_N % N , 1 : NCONS            ) 
         real(kind=RP)  :: FaN    ( 0 : ed % spA_N % N , 1 : NCONS            ) 
         real(kind=RP)  :: GauxLN ( 0 : ed % spA_N % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: GauxRN ( 0 : ed % spA_N % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: FvS    ( 0 : ed % spA_S % N , 1 : NCONS            ) 
         real(kind=RP)  :: FaS    ( 0 : ed % spA_S % N , 1 : NCONS            ) 
         real(kind=RP)  :: GauxLS ( 0 : ed % spA_S % N , 1 : NCONS , 1 : NDIM ) 
         real(kind=RP)  :: GauxRS ( 0 : ed % spA_S % N , 1 : NCONS , 1 : NDIM ) 
#endif
         integer       :: eq , dimID , i 
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif
!
!        Get the solution projection onto the edge
!        -----------------------------------------
         call Mat_x_Mat ( ed % T_LN_FWD , ed % storage ( LEFT        )  % Q , QLN ) 
         call Mat_x_Mat ( ed % T_LS_FWD , ed % storage ( LEFT        )  % Q , QLS ) 
         call Mat_x_Mat ( ed % T_RN_FWD , ed % storage ( RIGHT_NORTH )  % Q , QRN ) 
         call Mat_x_Mat ( ed % T_RS_FWD , ed % storage ( RIGHT_SOUTH )  % Q , QRS ) 
!
!        Compute the inviscid Riemann solver
!        -----------------------------------
         FiN = InviscidMethod % RiemannSolver( ed % spA_N % N , QLN , QRN , ed % normal_N ) 
         FiS = InviscidMethod % RiemannSolver( ed % spA_S % N , QLS , QRS , ed % normal_S ) 
#ifdef NAVIER_STOKES
!
!        Get the gradient projection onto the edge
!        -----------------------------------------
         do eq = 1 , NCONS
            call Mat_x_Mat ( ed % T_LN_FWD , ed % storage ( LEFT        )  % dQ(:,:,eq) , dQLN(:,:,eq) ) 
            call Mat_x_Mat ( ed % T_LS_FWD , ed % storage ( LEFT        )  % dQ(:,:,eq) , dQLS(:,:,eq) ) 
            call Mat_x_Mat ( ed % T_RN_FWD , ed % storage ( RIGHT_NORTH )  % dQ(:,:,eq) , dQRN(:,:,eq) ) 
            call Mat_x_Mat ( ed % T_RS_FWD , ed % storage ( RIGHT_SOUTH )  % dQ(:,:,eq) , dQRS(:,:,eq) ) 
         end do
!
!        Compute the viscous Riemann solver
!        ----------------------------------
         FvN = ViscousMethod % RiemannSolver( ed , ed % spA_N % N , ed % invh , QLN , QRN , dQLN , dQRN , ed % normal_N )
         FvS = ViscousMethod % RiemannSolver( ed , ed % spA_S % N , ed % invh , QLS , QRS , dQLS , dQRS , ed % normal_S )
!
!        Compute the artificial dissipation Riemann solver
!        -------------------------------------------------
         FaN = ArtificialDissipation % ComputeFaceFluxes( ed , QLN , QRN , dQLN , dQRN , ed % normal_N )
         FaS = ArtificialDissipation % ComputeFaceFluxes( ed , QLS , QRS , dQLS , dQRS , ed % normal_S )
#endif
!
!        The resulting flux is: FStar = ( Inviscid - Viscous - ArtificialDissipation ) dS
!        --------------------------------------------------------------------------------
         do eq = 1 , NCONS
#ifdef NAVIER_STOKES
            FStarN(:,eq) = ( FiN(:,eq) - FvN(:,eq) - FaN(:,eq) ) * ed % dS_N
            FStarS(:,eq) = ( FiS(:,eq) - FvS(:,eq) - FaS(:,eq) ) * ed % dS_S
#else
            FStarN(:,eq) = FiN(:,eq) * ed % dS_N
            FStarS(:,eq) = FiS(:,eq) * ed % dS_S
#endif
         end do
!
!        Return the resulting Riemann flux to each element frame
!        -------------------------------------------------------
         FL =  Mat_x_Mat_F( ed % T_LN_BKW , FStarN , ed % spA % N + 1, NCONS ) + Mat_x_Mat_F( ed % T_LS_BKW , FStarS , ed % spA % N + 1, NCONS )
         call Mat_x_Mat( ed % T_RN_BKW , FStarN , FRN )
         call Mat_x_Mat( ed % T_RS_BKW , FStarS , FRS )
      
#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Compute the gradient Riemann solver
!           -----------------------------------
            call ViscousMethod % GradientRiemannSolver ( ed , ed % spA_N % N , QLN , QRN , ed % normal_N , GauxLN , GauxRN ) 
            call ViscousMethod % GradientRiemannSolver ( ed , ed % spA_S % N , QLS , QRS , ed % normal_S , GauxLS , GauxRS ) 
!
!           Scale it with the surface Jacobian
!           ----------------------------------
            do dimID = 1 , NDIM     ; do eq = 1 , NCONS
               GauxLN(:,eq,dimID) = GauxLN(:,eq,dimID) * ed % dS_N
               GauxLS(:,eq,dimID) = GauxLS(:,eq,dimID) * ed % dS_S
               GauxRN(:,eq,dimID) = GauxRN(:,eq,dimID) * ed % dS_N
               GauxRS(:,eq,dimID) = GauxRS(:,eq,dimID) * ed % dS_S
            end do                  ; end do
!
!           Return the resulting Riemann flux to each element frame
!           -------------------------------------------------------
            GL(:,:,IX) = Mat_x_Mat_F( ed % T_LN_BKW , GauxLN(:,:,IX) , ed % spA % N + 1 , NCONS ) + Mat_x_Mat_F( ed % T_LS_BKW , GauxLS(:,:,IX) , ed % spA % N + 1 , NCONS )
            GL(:,:,IY) = Mat_x_Mat_F( ed % T_LN_BKW , GauxLN(:,:,IY) , ed % spA % N + 1 , NCONS ) + Mat_x_Mat_F( ed % T_LS_BKW , GauxLS(:,:,IY) , ed % spA % N + 1 , NCONS )

            call Mat_x_Mat( ed % T_RN_BKW , GauxRN(:,:,IX) , GRN(:,:,IX) )
            call Mat_x_Mat( ed % T_RN_BKW , GauxRN(:,:,IY) , GRN(:,:,IY) )

            call Mat_x_Mat( ed % T_RS_BKW , GauxRS(:,:,IX) , GRS(:,:,IX) )
            call Mat_x_Mat( ed % T_RS_BKW , GauxRS(:,:,IY) , GRS(:,:,IY) )

         end if
#endif
!
!        Change RIGHT fluxes signs to account for the normal direction
!        ------------------------------------------------------------- 
         FRN = -FRN
         FRS = -FRS

      end subroutine ComputeRiemannSolver_CurvedSubdividedEdge

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
      
         do iXi = 0 , N
            normal(IX:IY,iXi) = ed % n(IX:IY,0) 
         end do
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
         integer                       :: iXi , dimID , eq , i
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

            Fv = ViscousMethod % RiemannSolver_Adiabatic( ed , N , ed % invh , Q , dQ , Qb , ed % n ) 
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , Q , Qb , dQ , dQ , ed % n )  

            do eq = 1 , NCONS
               Fv(:,eq) = Fv(:,eq) * ed % dS
               Fa(:,eq) = Fa(:,eq) * ed % dS
            end do

            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               Gv = ViscousMethod % GradientRiemannSolver_Adiabatic( ed , N , Q , Qb , ed % n ) 
               
               do dimID = 1 , NDIM     ; do eq = 1 , NCONS
                  Gv(:,eq,dimID) = Gv(:,eq,dimID) * ed % dS
               end do                  ; end do

            end if

         else
!
!           Dirichlet/Neumann boundary conditions
!           -------------------------------------
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % storage(1) % Q , ed % uSB , ed % storage(1) % dQ , ed % storage(1) % dQ , ed % n )
            
            do eq = 1 , NCONS
               Fa(:,eq) = Fa(:,eq) * ed % dS
            end do

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
