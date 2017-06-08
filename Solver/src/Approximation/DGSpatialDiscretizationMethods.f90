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
!
   private
   public DGSpatial_Initialization  , DGSpatial_computeTimeDerivative , DGSpatial_interpolateSolutionToBoundaries
   public DGSpatial_newTimeStep
#ifdef NAVIER_STOKES
   public ArtificialDissipation
#endif
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
!$omp parallel
!
!        Prepare the mesh for a new iteration
!        ------------------------------------
         call DGSpatial_newTimeStep( mesh , time )
!
!        Compute QDot
!        ------------
         call DGSpatial_computeQDot( mesh )
!$omp end parallel

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
!$omp do
         do zoneID = 1 , size(mesh % zones) - 1
            call mesh % zones(zoneID) % UpdateSolution(time)
         end do 
!$omp end do

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
         integer                 :: i , j 
         class(Edge_t), pointer  :: f
!
!        ************
!        Volume loops
!        ************
!
!$omp do schedule(runtime)
         do eID = 1 , mesh % no_of_elements
            call DGSpatial_QDotVolumeLoop( mesh % elements(eID) ) 
         end do
!$omp end do
!
!        **********
!        Face loops
!        **********
!
!$omp barrier
!$omp master
         do edID = 1 , mesh % no_of_edges
            f => mesh % edges(edID) % f
            select type ( f ) 
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
!$omp end master
!
!        ***********************
!        Scale with the jacobian
!        ***********************
!
!$omp barrier
!$omp do private(i,j) schedule(runtime)
         do eID = 1 , mesh % no_of_elements 
            do j = 0 , mesh % elements(eID) % spA % N ; do i = 0 , mesh % elements(eID) % spA % N
               mesh % elements(eID) % QDot(:,i,j) = mesh % elements(eID) % QDot(:,i,j) / mesh % elements(eID) % Jac(i,j)
            end do                                    ; end do
         end do
!$omp end do

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
         real(kind=RP)       :: Fi ( 1:NCONS , 0:e % spA % N , 0:e % spA % N , 1:NDIM )
         real(kind=RP)       :: Fv ( 1:NCONS , 0:e % spA % N , 0:e % spA % N , 1:NDIM )
         real(kind=RP)       :: Fa ( 1:NCONS , 0:e % spA % N , 0:e % spA % N , 1:NDIM )
         real(kind=RP)       :: F  ( 1:NCONS , 0:e % spA % N , 0:e % spA % N , 1:NDIM )
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

         QDot(1:,0:,0:) => e % QDot
         QDot = ScalarWeakIntegrals % StdVolumeGreen(e , F)
  
      end subroutine DGSpatial_QDotVolumeLoop

      subroutine DGSpatial_QDotFaceLoop_Interior( ed )
         use QuadElementClass
         implicit none
         type(Edge_t)           :: ed
         real(kind=RP)          :: FiL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  )
         real(kind=RP)          :: FiR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  )
         real(kind=RP)          :: FvL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  )
         real(kind=RP)          :: FvR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  )
         real(kind=RP)          :: GL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  , 1 : NDIM)
         real(kind=RP)          :: GR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  , 1 : NDIM)
         real(kind=RP)          :: FL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  )
         real(kind=RP)          :: FR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  )
         real(kind=RP), pointer :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FR) and the Gradient Riemann Solver (GL,GR)
!        --------------------------------------------------------------------------
         call ComputeRiemannSolver_InteriorEdge( ed , FL , FR , GL , GR)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT element
!        -----------------------------------------
         QDot(1:,0:,0:) => ed % quads(RIGHT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT , FR ) 

#ifdef NAVIER_STOKES
!
!        Add the gradient fluxes term
!        ----------------------------
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Add the contribution to the LEFT element
!           ----------------------------------------
            QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , LEFT  , GL ) 
!
!           Add the contribution to the RIGHT element
!           -----------------------------------------
            QDot(1:,0:,0:) => ed % quads(RIGHT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT , GR ) 
!
         end if
#endif
           

      end subroutine DGSpatial_QDotFaceLoop_Interior

      subroutine DGSpatial_QDotFaceLoop_CurvedInterior( ed )
         use QuadElementClass
         implicit none
         type(CurvedEdge_t)         :: ed
         real(kind=RP)           :: FiL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  )
         real(kind=RP)           :: FiR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  )
         real(kind=RP)           :: FvL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  )
         real(kind=RP)           :: FvR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  )
         real(kind=RP)           :: GL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  , 1 : NDIM)
         real(kind=RP)           :: GR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  , 1 : NDIM)
         real(kind=RP)           :: FL ( 1:NCONS , 0 : ed % storage ( LEFT  ) % spA % N  )
         real(kind=RP)           :: FR ( 1:NCONS , 0 : ed % storage ( RIGHT ) % spA % N  )
         real(kind=RP), pointer  :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FR) and the Gradient Riemann Solver (GL,GR)
!        --------------------------------------------------------------------------
         call ComputeRiemannSolver_InteriorCurvedEdge( ed , FL , FR , GL , GR)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT element
!        -----------------------------------------
         QDot(1:,0:,0:) => ed % quads(RIGHT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT , FR ) 

#ifdef NAVIER_STOKES
!
!        Add the gradient fluxes term
!        ----------------------------
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Add the contribution to the LEFT element
!           ----------------------------------------
            QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , LEFT  , GL ) 
!
!           Add the contribution to the RIGHT element
!           -----------------------------------------
            QDot(1:,0:,0:) => ed % quads(RIGHT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT , GR ) 
!
         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_CurvedInterior

      subroutine DGSpatial_QDotFaceLoop_SubdividedEdge( ed )
         use QuadElementClass
         implicit none
         type(SubdividedEdge_t) :: ed
         real(kind=RP)          :: FiL  ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N             )
         real(kind=RP)          :: FiRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N             )
         real(kind=RP)          :: FiRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N             )
         real(kind=RP)          :: FvL  ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N             )
         real(kind=RP)          :: FvRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N             )
         real(kind=RP)          :: FvRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N             )
         real(kind=RP)          :: GL   ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N  , 1 : NDIM )
         real(kind=RP)          :: GRN  ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N  , 1 : NDIM )
         real(kind=RP)          :: GRS  ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N  , 1 : NDIM )
         real(kind=RP)          :: FL   ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N             )
         real(kind=RP)          :: FRN  ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N             )
         real(kind=RP)          :: FRS  ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N             )
         real(kind=RP), pointer :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FRN,FRS) and the Gradient Riemann Solver (GL,GRN,GRS)
!        ------------------------------------------------------------------------------------
         call ComputeRiemannSolver_SubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT-NORTH element
!        -----------------------------------------------
         QDot(1:,0:,0:) => ed % quads(RIGHT_NORTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_NORTH , FRN ) 
!
!        Add the contribution to the RIGHT-SOUTH element
!        -----------------------------------------------
         QDot(1:,0:,0:) => ed % quads(RIGHT_SOUTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_SOUTH , FRS ) 

#ifdef NAVIER_STOKES
!
!        Add the gradient fluxes term
!        ----------------------------
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Add the contribution to the LEFT element
!           ----------------------------------------
            QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , LEFT  , GL ) 
!
!           Add the contribution to the RIGHT-NORTH element
!           -----------------------------------------------
            QDot(1:,0:,0:) => ed % quads(RIGHT_NORTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_NORTH , GRN ) 
!
!           Add the contribution to the RIGHT-SOUTH element
!           -----------------------------------------------
            QDot(1:,0:,0:) => ed % quads(RIGHT_SOUTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_SOUTH , GRS ) 
!
         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_SubdividedEdge

      subroutine DGSpatial_QDotFaceLoop_CurvedSubdividedEdge( ed )
         use QuadElementClass
         implicit none
         type(CurvedSubdividedEdge_t) :: ed
         real(kind=RP)                :: FiL  ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N             )
         real(kind=RP)                :: FiRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N             )
         real(kind=RP)                :: FiRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N             )
         real(kind=RP)                :: FvL  ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N             )
         real(kind=RP)                :: FvRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N             )
         real(kind=RP)                :: FvRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N             )
         real(kind=RP)                :: GL   ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N  , 1 : NDIM )
         real(kind=RP)                :: GRN  ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N  , 1 : NDIM )
         real(kind=RP)                :: GRS  ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N  , 1 : NDIM )
         real(kind=RP)                :: FL   ( 1:NCONS , 0 : ed % storage ( LEFT        ) % spA % N             )
         real(kind=RP)                :: FRN  ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH ) % spA % N             )
         real(kind=RP)                :: FRS  ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH ) % spA % N             )
         real(kind=RP), pointer       :: QDot(:,:,:)
!
!        Compute the Riemann Solver (FL,FRN,FRS) and the Gradient Riemann Solver (GL,GRN,GRS)
!        ------------------------------------------------------------------------------------
         call ComputeRiemannSolver_CurvedSubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS)
!
!        Add the contribution to the LEFT element
!        ----------------------------------------
         QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , LEFT  , FL ) 
!
!        Add the contribution to the RIGHT-NORTH element
!        -----------------------------------------------
         QDot(1:,0:,0:) => ed % quads(RIGHT_NORTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_NORTH , FRN ) 
!
!        Add the contribution to the RIGHT-SOUTH element
!        -----------------------------------------------
         QDot(1:,0:,0:) => ed % quads(RIGHT_SOUTH) % e % QDot
         QDot = QDot - ScalarWeakIntegrals % StdFace ( ed , RIGHT_SOUTH , FRS ) 

#ifdef NAVIER_STOKES
!
!        Add the gradient fluxes term
!        ----------------------------
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
!
!           Add the contribution to the LEFT element
!           ----------------------------------------
            QDot(1:,0:,0:) => ed % quads(LEFT) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , LEFT  , GL ) 
!
!           Add the contribution to the RIGHT-NORTH element
!           -----------------------------------------------
            QDot(1:,0:,0:) => ed % quads(RIGHT_NORTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_NORTH , GRN ) 
!
!           Add the contribution to the RIGHT-SOUTH element
!           -----------------------------------------------
            QDot(1:,0:,0:) => ed % quads(RIGHT_SOUTH) % e % QDot
            QDot = QDot - ScalarWeakIntegrals % StdGradientFace ( ed , RIGHT_SOUTH , GRS ) 
!
         end if
#endif

      end subroutine DGSpatial_QDotFaceLoop_CurvedSubdividedEdge

      subroutine DGSpatial_QDotFaceLoop_StraightBdry( ed )
         use QuadElementClass
         implicit none
         type(StraightBdryEdge_t)  :: ed
         real(kind=RP)             :: G ( 1:NCONS , 0 : ed % spA % N  , 1:NDIM )  
         real(kind=RP)             :: F ( 1:NCONS , 0 : ed % spA % N )
         real(kind=RP), pointer :: QDot(:,:,:)

         call ComputeRiemannSolver_StraightBdryEdge( ed , F , G )

         QDot(1:,0:,0:) => ed % quads(1) % e % QDot
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
         real(kind=RP)             :: G ( 1:NCONS , 0 : ed % spA % N  , 1:NDIM)
         real(kind=RP)             :: F ( 1:NCONS , 0 : ed % spA % N )
         real(kind=RP), pointer :: QDot(:,:,:)

         call ComputeRiemannSolver_CurvedBdryEdge( ed , F , G )

         QDot(1:,0:,0:) => ed % quads(1) % e % QDot
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
         real(kind=RP)               :: FL( 1:NCONS , 0 : ed % storage(LEFT ) % spA % N )
         real(kind=RP)               :: FR( 1:NCONS , 0 : ed % storage(RIGHT) % spA % N )
         real(kind=RP), intent (out) :: GL( 1:NCONS , 0 : ed % storage(LEFT ) % spA % N  , 1:NDIM)
         real(kind=RP), intent (out) :: GR( 1:NCONS , 0 : ed % storage(RIGHT) % spA % N  , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: Fi    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP)          :: Fstar ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) , target :: QL    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) , target :: QR    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) , target :: dQL   ( 1:NCONS , 0 : ed % spA % N , 1 : NDIM  ) 
         real(kind=RP) , target :: dQR   ( 1:NCONS , 0 : ed % spA % N , 1 : NDIM  ) 
#ifdef NAVIER_STOKES
         real(kind=RP) :: Fv    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) :: Fa    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) :: GauxL ( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM ) 
         real(kind=RP) :: GauxR ( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM ) 
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
         Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , QL , QR , dQL , dQR , normal )
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
!
!        Compute the Gradient Riemann solver
!        -----------------------------------
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
         real(kind=RP)               :: FL( 1:NCONS , 0 : ed % storage(LEFT ) % spA % N )
         real(kind=RP)               :: FR( 1:NCONS , 0 : ed % storage(RIGHT) % spA % N )
         real(kind=RP), intent (out) :: GL( 1:NCONS , 0 : ed % storage(LEFT ) % spA % N  , 1:NDIM)
         real(kind=RP), intent (out) :: GR( 1:NCONS , 0 : ed % storage(RIGHT) % spA % N  , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: Fi    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP)          :: Fstar ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) , target :: QL    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) , target :: QR    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP) , target :: dQL   ( 1:NCONS , 0 : ed % spA % N , 1 : NDIM  ) 
         real(kind=RP) , target :: dQR   ( 1:NCONS , 0 : ed % spA % N , 1 : NDIM  ) 
#ifdef NAVIER_STOKES
         real(kind=RP)  :: Fv    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP)  :: Fa    ( 1:NCONS , 0 : ed % spA % N             ) 
         real(kind=RP)  :: GauxL ( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM ) 
         real(kind=RP)  :: GauxR ( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM ) 
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
         Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , QL , QR , dQL , dQR , ed % n )
#endif
!
!        The resulting flux is: FStar = ( Inviscid - Viscous - ArtificialDissipation ) dS
!        --------------------------------------------------------------------------------
         do i = 0 , ed % spA % N 
#ifdef NAVIER_STOKES
            FStar(:,i) = ( Fi(:,i) - Fv(:,i) - Fa(:,i) ) * ed % dS(i)
#else
            FStar(:,i) = Fi(:,i) * ed % dS(i) 
#endif
         end do
!
!        Return the resulting Riemann flux to each element frame
!        -------------------------------------------------------
         call ed % ProjectFluxes( ed , FStar , FL , FR )

#ifdef NAVIER_STOKES
         if ( ViscousMethod % computeRiemannGradientFluxes ) then
               call ViscousMethod % GradientRiemannSolver ( ed , ed % spA % N , QL , QR , ed % n , GauxL , GauxR ) 

               do dimID = 1 , NDIM    ; do i = 0 , ed % spA % N 
                  GauxL(:,i,dimID) = GauxL(:,i,dimID) * ed % dS(i)
                  GauxR(:,i,dimID) = GauxR(:,i,dimID) * ed % dS(i)
               end do                 ; end do

               call ed % ProjectGradientFluxes( ed , GauxL , GauxR , GL , GR )

         end if
#endif

      end subroutine ComputeRiemannSolver_InteriorCurvedEdge

      subroutine ComputeRiemannSolver_SubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS )
         use QuadElementClass
         use MatrixOperations
         implicit none
         type(SubdividedEdge_t)              :: ed
         real(kind=RP) ,         intent(out) :: FL  ( 1:NCONS , 0 : ed % storage ( LEFT        )  % spA % N           )
         real(kind=RP) ,         intent(out) :: FRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH )  % spA % N           )
         real(kind=RP) ,         intent(out) :: FRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH )  % spA % N           )
         real(kind=RP) ,         intent(out) :: GL  ( 1:NCONS , 0 : ed % storage ( LEFT        )  % spA % N  , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH )  % spA % N  , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH )  % spA % N  , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: FiN    ( 1:NCONS , 0 : ed % spA_N % N             )      ! -----------------------------
         real(kind=RP)          :: FstarN ( 1:NCONS , 0 : ed % spA_N % N             )      !     Variables in the NORTH
         real(kind=RP) , target :: QLN    ( 1:NCONS , 0 : ed % spA_N % N             )      !  mortar.
         real(kind=RP) , target :: QRN    ( 1:NCONS , 0 : ed % spA_N % N             )      !  
         real(kind=RP) , target :: dQLN   ( 1:NCONS , 0 : ed % spA_N % N , 1 : NDIM  )      ! 
         real(kind=RP) , target :: dQRN   ( 1:NCONS , 0 : ed % spA_N % N , 1 : NDIM  )      ! -----------------------------
         real(kind=RP)          :: FiS    ( 1:NCONS , 0 : ed % spA_S % N             )      ! -----------------------------
         real(kind=RP)          :: FstarS ( 1:NCONS , 0 : ed % spA_S % N             )      !     Variables in the SOUTH
         real(kind=RP) , target :: QLS    ( 1:NCONS , 0 : ed % spA_S % N             )      !  mortar.
         real(kind=RP) , target :: QRS    ( 1:NCONS , 0 : ed % spA_S % N             )      !  
         real(kind=RP) , target :: dQLS   ( 1:NCONS , 0 : ed % spA_S % N , 1 : NDIM  )      ! 
         real(kind=RP) , target :: dQRS   ( 1:NCONS , 0 : ed % spA_S % N , 1 : NDIM  )      ! -----------------------------
#ifdef NAVIER_STOKES
         real(kind=RP)  :: FvN    ( 1:NCONS , 0 : ed % spA_N % N             ) 
         real(kind=RP)  :: FaN    ( 1:NCONS , 0 : ed % spA_N % N             ) 
         real(kind=RP)  :: GauxLN ( 1:NCONS , 0 : ed % spA_N % N  , 1 : NDIM ) 
         real(kind=RP)  :: GauxRN ( 1:NCONS , 0 : ed % spA_N % N  , 1 : NDIM ) 
         real(kind=RP)  :: FvS    ( 1:NCONS , 0 : ed % spA_S % N             ) 
         real(kind=RP)  :: FaS    ( 1:NCONS , 0 : ed % spA_S % N             ) 
         real(kind=RP)  :: GauxLS ( 1:NCONS , 0 : ed % spA_S % N  , 1 : NDIM ) 
         real(kind=RP)  :: GauxRS ( 1:NCONS , 0 : ed % spA_S % N  , 1 : NDIM ) 
#endif
         real(kind=RP) :: normal_N(NDIM , 0 : ed % spA_N % N )
         real(kind=RP) :: normal_S(NDIM , 0 : ed % spA_S % N )
         integer       :: eq , iDim , i , l , dimID
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
         QLN = 0.0_RP
         do l = 0 , ed % storage(LEFT) % spA % N    ; do i = 0 , ed % spA_N % N 
            QLN(:,i)    = QLN(:,i) + ed % T_LN_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                                     ; end do

         QLS = 0.0_RP
         do l = 0 , ed % storage(LEFT) % spA % N    ; do i = 0 , ed % spA_S % N 
            QLS(:,i)    = QLS(:,i) + ed % T_LS_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                                     ; end do

         QRN = 0.0_RP
         do l = 0 , ed % storage(RIGHT_NORTH) % spA % N    ; do i = 0 , ed % spA_N % N 
            QRN(:,i)    = QRN(:,i) + ed % T_RN_FWD(i,l) * ed % storage(RIGHT_NORTH) % Q(:,l)
         end do                                            ; end do

         QRS = 0.0_RP
         do l = 0 , ed % storage(RIGHT_SOUTH) % spA % N    ; do i = 0 , ed % spA_S % N 
            QRS(:,i)    = QRS(:,i) + ed % T_RS_FWD(i,l) * ed % storage(RIGHT_SOUTH) % Q(:,l)
         end do                                            ; end do
!
!        Compute the inviscid Riemann solver
!        -----------------------------------
         FiN = InviscidMethod % RiemannSolver( ed % spA_N % N , QLN , QRN , normal_N ) 
         FiS = InviscidMethod % RiemannSolver( ed % spA_S % N , QLS , QRS , normal_S ) 
#ifdef NAVIER_STOKES
!
!        Get the gradient projection onto the edge
!        -----------------------------------------
         dQLN = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(LEFT) % spA % N ; do i = 0 , ed % spA_N % N 
            dQLN(:,i,dimID)   = dQLN(:,i,dimID)    + ed % T_LN_FWD(i,l) * ed % storage(LEFT) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 

         dQLS = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(LEFT) % spA % N ; do i = 0 , ed % spA_S % N 
            dQLS(:,i,dimID)   = dQLS(:,i,dimID)    + ed % T_LS_FWD(i,l) * ed % storage(LEFT) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 

         dQRN = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(RIGHT_NORTH) % spA % N ; do i = 0 , ed % spA_N % N 
            dQRN(:,i,dimID)   = dQRN(:,i,dimID)    + ed % T_RN_FWD(i,l) * ed % storage(RIGHT_NORTH) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 

         dQRS = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(RIGHT_SOUTH) % spA % N ; do i = 0 , ed % spA_S % N 
            dQRS(:,i,dimID)   = dQRS(:,i,dimID)    + ed % T_RS_FWD(i,l) * ed % storage(RIGHT_SOUTH) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 
!
!        Compute the viscous Riemann solver
!        ----------------------------------
         FvN = ViscousMethod % RiemannSolver( ed , ed % spA_N % N , ed % invh , QLN , QRN , dQLN , dQRN , normal_N )
         FvS = ViscousMethod % RiemannSolver( ed , ed % spA_S % N , ed % invh , QLS , QRS , dQLS , dQRS , normal_S )
!
!        Compute the artificial dissipation Riemann solver
!        -------------------------------------------------
         FaN = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA_N % N , QLN , QRN , dQLN , dQRN , normal_N )
         FaS = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA_S % N , QLS , QRS , dQLS , dQRS , normal_S )
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
         FL = 0.0_RP
         do l = 0 , ed % spA_N % N  ; do i = 0 , ed % storage(LEFT) % spA % N  
            FL(:,i) = FL(:,i) + ed % T_LN_BKW(i,l) * FStarN(:,l)
         end do                     ; end do

         do l = 0 , ed % spA_S % N  ; do i = 0 , ed % storage(LEFT) % spA % N  
            FL(:,i) = FL(:,i) + ed % T_LS_BKW(i,l) * FStarS(:,l) 
         end do                     ; end do

         FRN = 0.0_RP
         do l = 0 , ed % spA_N % N    ; do i = 0 , ed % storage(RIGHT_NORTH) % spA % N  
            FRN(:,i) = FRN(:,i) - ed % T_RN_BKW(i,l) * FStarN(:,l)  
         end do                     ; end do

         FRS = 0.0_RP
         do l = 0 , ed % spA_S % N    ; do i = 0 , ed % storage(RIGHT_SOUTH) % spA % N  
            FRS(:,i) = FRS(:,i) - ed % T_RS_BKW(i,l) * FStarS(:,l)  
         end do                     ; end do
      
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
            GL = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_N % N    ; do i = 0 , ed % storage(LEFT) % spA % N  
               GL(:,i,dimID) = GL(:,i,dimID) + ed % T_LN_BKW(i,l) * GauxLN(:,l,dimID)
            end do              ; end do                     ; end do

            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_S % N    ; do i = 0 , ed % storage(LEFT) % spA % N  
               GL(:,i,dimID) = GL(:,i,dimID) + ed % T_LS_BKW(i,l) * GauxLS(:,l,dimID) 
            end do              ; end do                     ; end do
 
            GRN = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_N % N    ; do i = 0 , ed % storage(RIGHT_NORTH) % spA % N  
               GRN(:,i,dimID) = GRN(:,i,dimID) + ed % T_RN_BKW(i,l) * GauxRN(:,l,dimID)  
            end do               ; end do                     ; end do

            GRS = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_S % N    ; do i = 0 , ed % storage(RIGHT_SOUTH) % spA % N  
               GRS(:,i,dimID) = GRS(:,i,dimID) + ed % T_RS_BKW(i,l) * GauxRS(:,l,dimID)  
            end do               ; end do                     ; end do
 
         end if
#endif

      end subroutine ComputeRiemannSolver_SubdividedEdge

      subroutine ComputeRiemannSolver_CurvedSubdividedEdge( ed , FL , FRN , FRS , GL , GRN , GRS )
         use QuadElementClass
         use MatrixOperations
         implicit none
         type(CurvedSubdividedEdge_t)        :: ed
         real(kind=RP) ,         intent(out) :: FL  ( 1:NCONS , 0 : ed % storage ( LEFT        )  % spA % N            )
         real(kind=RP) ,         intent(out) :: FRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH )  % spA % N            )
         real(kind=RP) ,         intent(out) :: FRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH )  % spA % N            )
         real(kind=RP) ,         intent(out) :: GL  ( 1:NCONS , 0 : ed % storage ( LEFT        )  % spA % N  , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRN ( 1:NCONS , 0 : ed % storage ( RIGHT_NORTH )  % spA % N  , 1:NDIM )
         real(kind=RP) ,         intent(out) :: GRS ( 1:NCONS , 0 : ed % storage ( RIGHT_SOUTH )  % spA % N  , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: FiN    ( 1:NCONS , 0 : ed % spA_N % N             )      ! -----------------------------
         real(kind=RP)          :: FstarN ( 1:NCONS , 0 : ed % spA_N % N             )      !     Variables in the NORTH
         real(kind=RP) , target :: QLN    ( 1:NCONS , 0 : ed % spA_N % N             )      !  mortar.
         real(kind=RP) , target :: QRN    ( 1:NCONS , 0 : ed % spA_N % N             )      !  
         real(kind=RP) , target :: dQLN   ( 1:NCONS , 0 : ed % spA_N % N , 1 : NDIM  )      ! 
         real(kind=RP) , target :: dQRN   ( 1:NCONS , 0 : ed % spA_N % N , 1 : NDIM  )      ! -----------------------------
         real(kind=RP)          :: FiS    ( 1:NCONS , 0 : ed % spA_S % N             )      ! -----------------------------
         real(kind=RP)          :: FstarS ( 1:NCONS , 0 : ed % spA_S % N             )      !     Variables in the SOUTH
         real(kind=RP) , target :: QLS    ( 1:NCONS , 0 : ed % spA_S % N             )      !  mortar.
         real(kind=RP) , target :: QRS    ( 1:NCONS , 0 : ed % spA_S % N             )      !  
         real(kind=RP) , target :: dQLS   ( 1:NCONS , 0 : ed % spA_S % N , 1 : NDIM  )      ! 
         real(kind=RP) , target :: dQRS   ( 1:NCONS , 0 : ed % spA_S % N , 1 : NDIM  )      ! -----------------------------
#ifdef NAVIER_STOKES
         real(kind=RP)  :: FvN    ( 1:NCONS , 0 : ed % spA_N % N             ) 
         real(kind=RP)  :: FaN    ( 1:NCONS , 0 : ed % spA_N % N             ) 
         real(kind=RP)  :: GauxLN ( 1:NCONS , 0 : ed % spA_N % N  , 1 : NDIM ) 
         real(kind=RP)  :: GauxRN ( 1:NCONS , 0 : ed % spA_N % N  , 1 : NDIM ) 
         real(kind=RP)  :: FvS    ( 1:NCONS , 0 : ed % spA_S % N             ) 
         real(kind=RP)  :: FaS    ( 1:NCONS , 0 : ed % spA_S % N             ) 
         real(kind=RP)  :: GauxLS ( 1:NCONS , 0 : ed % spA_S % N  , 1 : NDIM ) 
         real(kind=RP)  :: GauxRS ( 1:NCONS , 0 : ed % spA_S % N  , 1 : NDIM ) 
#endif
         integer       :: eq , dimID , i , l 
!
!        Compute the edge artificial dissipation
!        ---------------------------------------
#ifdef NAVIER_STOKES
         ed % mu_a = ArtificialDissipation % ComputeEdgeViscosity( ed )
#endif
!
!        Get the solution projection onto the edge
!        -----------------------------------------
         QLN = 0.0_RP
         do l = 0 , ed % storage(LEFT) % spA % N    ; do i = 0 , ed % spA_N % N 
            QLN(:,i)    = QLN(:,i) + ed % T_LN_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                                     ; end do

         QLS = 0.0_RP
         do l = 0 , ed % storage(LEFT) % spA % N    ; do i = 0 , ed % spA_S % N 
            QLS(:,i)    = QLS(:,i) + ed % T_LS_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                                     ; end do

         QRN = 0.0_RP
         do l = 0 , ed % storage(RIGHT_NORTH) % spA % N    ; do i = 0 , ed % spA_N % N 
            QRN(:,i)    = QRN(:,i) + ed % T_RN_FWD(i,l) * ed % storage(RIGHT_NORTH) % Q(:,l)
         end do                                            ; end do

         QRS = 0.0_RP
         do l = 0 , ed % storage(RIGHT_SOUTH) % spA % N    ; do i = 0 , ed % spA_S % N 
            QRS(:,i)    = QRS(:,i) + ed % T_RS_FWD(i,l) * ed % storage(RIGHT_SOUTH) % Q(:,l)
         end do                                            ; end do
!
!        Compute the inviscid Riemann solver
!        -----------------------------------
         FiN = InviscidMethod % RiemannSolver( ed % spA_N % N , QLN , QRN , ed % normal_N ) 
         FiS = InviscidMethod % RiemannSolver( ed % spA_S % N , QLS , QRS , ed % normal_S ) 
#ifdef NAVIER_STOKES
!
!        Get the gradient projection onto the edge
!        -----------------------------------------
         dQLN = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(LEFT) % spA % N ; do i = 0 , ed % spA_N % N 
            dQLN(:,i,dimID)   = dQLN(:,i,dimID)    + ed % T_LN_FWD(i,l) * ed % storage(LEFT) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 

         dQLS = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(LEFT) % spA % N ; do i = 0 , ed % spA_S % N 
            dQLS(:,i,dimID)   = dQLS(:,i,dimID)    + ed % T_LS_FWD(i,l) * ed % storage(LEFT) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 

         dQRN = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(RIGHT_NORTH) % spA % N ; do i = 0 , ed % spA_N % N 
            dQRN(:,i,dimID)   = dQRN(:,i,dimID)    + ed % T_RN_FWD(i,l) * ed % storage(RIGHT_NORTH) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 

         dQRS = 0.0_RP 
         do dimID = 1 , NDIM  ; do l = 0 , ed % storage(RIGHT_SOUTH) % spA % N ; do i = 0 , ed % spA_S % N 
            dQRS(:,i,dimID)   = dQRS(:,i,dimID)    + ed % T_RS_FWD(i,l) * ed % storage(RIGHT_SOUTH) % dQ(:,l,dimID)
         end do               ; end do                  ; end do 
!
!        Compute the viscous Riemann solver
!        ----------------------------------
         FvN = ViscousMethod % RiemannSolver( ed , ed % spA_N % N , ed % invh , QLN , QRN , dQLN , dQRN , ed % normal_N )
         FvS = ViscousMethod % RiemannSolver( ed , ed % spA_S % N , ed % invh , QLS , QRS , dQLS , dQRS , ed % normal_S )
!
!        Compute the artificial dissipation Riemann solver
!        -------------------------------------------------
         FaN = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA_N % N , QLN , QRN , dQLN , dQRN , ed % normal_N )
         FaS = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA_S % N , QLS , QRS , dQLS , dQRS , ed % normal_S )
#endif
!
!        The resulting flux is: FStar = ( Inviscid - Viscous - ArtificialDissipation ) dS
!        --------------------------------------------------------------------------------
         do i = 0 , ed % spA_N % N 
#ifdef NAVIER_STOKES
            FStarN(:,i) = ( FiN(:,i) - FvN(:,i) - FaN(:,i) ) * ed % dS_N(i)
#else
            FStarN(:,i) = FiN(:,i) * ed % dS_N(i)
#endif
         end do

         do i = 0 , ed % spA_S % N 
#ifdef NAVIER_STOKES
            FStarS(:,i) = ( FiS(:,i) - FvS(:,i) - FaS(:,i) ) * ed % dS_S(i)
#else
            FStarS(:,i) = FiS(:,i) * ed % dS_S(i)
#endif
         end do
!
!        Return the resulting Riemann flux to each element frame
!        -------------------------------------------------------
         FL = 0.0_RP
         do l = 0 , ed % spA_N % N  ; do i = 0 , ed % storage(LEFT) % spA % N  
            FL(:,i) = FL(:,i) + ed % T_LN_BKW(i,l) * FStarN(:,l)
         end do                     ; end do

         do l = 0 , ed % spA_S % N  ; do i = 0 , ed % storage(LEFT) % spA % N  
            FL(:,i) = FL(:,i) + ed % T_LS_BKW(i,l) * FStarS(:,l) 
         end do                     ; end do

         FRN = 0.0_RP
         do l = 0 , ed % spA_N % N    ; do i = 0 , ed % storage(RIGHT_NORTH) % spA % N  
            FRN(:,i) = FRN(:,i) - ed % T_RN_BKW(i,l) * FStarN(:,l)  
         end do                     ; end do

         FRS = 0.0_RP
         do l = 0 , ed % spA_S % N    ; do i = 0 , ed % storage(RIGHT_SOUTH) % spA % N  
            FRS(:,i) = FRS(:,i) - ed % T_RS_BKW(i,l) * FStarS(:,l)  
         end do                     ; end do
 
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
            do dimID = 1 , NDIM     ; do i = 0 , ed % spA_N % N
               GauxLN(:,i,dimID) = GauxLN(:,i,dimID) * ed % dS_N(i)
               GauxRN(:,i,dimID) = GauxRN(:,i,dimID) * ed % dS_N(i)
            end do                  ; end do

            do dimID = 1 , NDIM     ; do i = 0 , ed % spA_S % N
               GauxLS(:,i,dimID) = GauxLS(:,i,dimID) * ed % dS_S(i)
               GauxRS(:,i,dimID) = GauxRS(:,i,dimID) * ed % dS_S(i)
            end do                  ; end do
!
!           Return the resulting Riemann flux to each element frame
!           -------------------------------------------------------
            GL = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_N % N    ; do i = 0 , ed % storage(LEFT) % spA % N  
               GL(:,i,dimID) = GL(:,i,dimID) + ed % T_LN_BKW(i,l) * GauxLN(:,l,dimID)
            end do              ; end do                     ; end do

            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_S % N    ; do i = 0 , ed % storage(LEFT) % spA % N  
               GL(:,i,dimID) = GL(:,i,dimID) + ed % T_LS_BKW(i,l) * GauxLS(:,l,dimID) 
            end do              ; end do                     ; end do
 
            GRN = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_N % N    ; do i = 0 , ed % storage(RIGHT_NORTH) % spA % N  
               GRN(:,i,dimID) = GRN(:,i,dimID) + ed % T_RN_BKW(i,l) * GauxRN(:,l,dimID)  
            end do               ; end do                     ; end do

            GRS = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , ed % spA_S % N    ; do i = 0 , ed % storage(RIGHT_SOUTH) % spA % N  
               GRS(:,i,dimID) = GRS(:,i,dimID) + ed % T_RS_BKW(i,l) * GauxRS(:,l,dimID)  
            end do               ; end do                     ; end do
 
         end if
#endif

      end subroutine ComputeRiemannSolver_CurvedSubdividedEdge

      subroutine ComputeRiemannSolver_StraightBdryEdge( ed , F , G )
         use QuadElementClass
         implicit none
         type(StraightBdryEdge_t)      :: ed
         real(kind=RP)                 :: F( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: G( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                 :: Fi( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: Q (1:NCONS , 0 : ed % spA % N ) 
         real(kind=RP)                 :: Qb(1:NCONS , 0 : ed % spA % N )
#ifdef NAVIER_STOKES
         real(kind=RP)                 :: Fv( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: Fa( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: Gv( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM )
         real(kind=RP)                 :: Gaux( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM )
         real(kind=RP)                 :: dQ (1:NCONS , 0 : ed % spA % N , 1 : NDIM )
         real(kind=RP)                 :: dQb(1:NCONS , 0 : ed % spA % N , 1 : NDIM )
         real(kind=RP)                 :: Q1D  ( 1 : NCONS )
         real(kind=RP)                 :: Qb1D ( 1 : NCONS )
         real(kind=RP)                 :: dQ1D ( 1 : NCONS , 1 : NDIM  ) 
         real(kind=RP)                 :: dQb1D( 1 : NCONS , 1 : NDIM  ) 
#endif
         real(kind=RP)                 :: normal( NDIM , 0 : ed % spA % N )
         integer, pointer              :: N
         integer                       :: iXi , i , dimID
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
                  do i = 0 , ed % spA % N 
                     Qb(:,i) = ed % uB(:,ed % spA % N - i )
                  end do
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

               do i = 0 , ed % spA % N 
                  Qb(:,i)  = ed % uB(:,N-i)
               end do
         
               do dimID = 1 , NDIM  ; do i = 0 , ed % spA % N 
                  dQb(:,i,dimID) = ed % gB(:,ed % spA % N - i , dimID)
               end do               ; end do

               Fv = ViscousMethod % RiemannSolver( ed , N , ed % invh , Q , Qb , dQ , dQb , normal ) * ed % dS(0)
               Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , Q , Qb , dQ , dQb , normal ) * ed % dS(0)

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
               Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , Q , Qb , dQ , dQb , normal ) * ed % dS(0)

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
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , Q , Qb , dQ , dQ , normal ) * ed % dS(0)

            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               Gv = ViscousMethod % GradientRiemannSolver_Adiabatic( ed , N , Q , Qb , normal ) * ed % dS(0)

            end if

         else
!
!           Dirichlet/Neumann boundary conditions
!           -------------------------------------
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , ed % storage(1) % Q , ed % uSB , ed % storage(1) % dQ , ed % storage(1) % dQ , normal ) * ed % dS(0)
            do iXi = 0 , ed % spA % N
               select case ( ed % viscousBCType(iXi) )
   
                  case ( DIRICHLET )
   
                     Q1D  = ed % storage(1) % Q(:,iXi)
                     dQ1D = ed % storage(1) % dQ(:,iXi,:)
                     Qb1D  = ed % uSB(:,iXi)
                     Fv(:,iXi) = ViscousMethod % RiemannSolver_Dirichlet ( ed , ed % spA % N , ed % invh , Q1D , dQ1D , Qb1D , ed % n(IX:IY,0) ) * ed % dS(0)

                     if ( ViscousMethod % computeRiemannGradientFluxes ) then
                        Gv(:,iXi,:) = ViscousMethod % GradientRiemannSolver_BoundaryCondition( ed , Q1D , Qb1D , ed % n(IX:IY,0) ) * ed % dS(0)
                     end if

                  case ( NEUMANN )

                     Fv(:,iXi)   = 0.0_RP
                     Fa(:,iXi)   = 0.0_RP
                     Gv(:,iXi,:) = 0.0_RP
         
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
         real(kind=RP)                 :: F( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: G( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                 :: Fi( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: Q(1:NCONS , 0 : ed % spA % N ) 
         real(kind=RP)                 :: Qb(1:NCONS , 0 : ed % spA % N )
#ifdef NAVIER_STOKES
         real(kind=RP)                 :: Fv( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: Fa( 1:NCONS , 0 : ed % spA % N  )
         real(kind=RP)                 :: Gv( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM )
         real(kind=RP)                 :: Gaux( 1:NCONS , 0 : ed % spA % N  , 1 : NDIM )
         real(kind=RP)                 :: dQ (1:NCONS , 0 : ed % spA % N , 1 : NDIM )
         real(kind=RP)                 :: dQb(1:NCONS , 0 : ed % spA % N , 1 : NDIM )
         real(kind=RP)                 :: Q1D  ( 1 : NCONS )
         real(kind=RP)                 :: Qb1D ( 1 : NCONS )
         real(kind=RP)                 :: dQ1D ( 1 : NCONS , 1 : NDIM  ) 
         real(kind=RP)                 :: dQb1D( 1 : NCONS , 1 : NDIM  ) 
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
                  Fi = RiemannSolver ( ed % spA % N , Q , Qb , ed % n ) 
                  do i = 0 , ed % spA % N
                     Fi(:,i) = Fi(:,i) * ed % dS(i)
                  end do

               else
                  Q = ed % storage(1) % Q

                  do i = 0 , ed % spA % N 
                     Qb(:,i) = ed % uB( : , ed % spA % N - i )
                  end do

                  Fi = RiemannSolver ( ed % spA % N , Q , Qb , ed % n ) 

                  do i = 0 , ed % spA % N
                     Fi(:,i) = Fi(:,i) * ed % dS(i)
                  end do

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
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , Q , Qb , dQ , dQ , ed % n )  

            do i = 0 , ed % spA % N 
               Fv(:,i) = Fv(:,i) * ed % dS(i)
               Fa(:,i) = Fa(:,i) * ed % dS(i)
            end do

            if ( ViscousMethod % computeRiemannGradientFluxes ) then
               Gv = ViscousMethod % GradientRiemannSolver_Adiabatic( ed , N , Q , Qb , ed % n ) 
               
               do dimID = 1 , NDIM     ; do i = 0 , ed % spA % N 
                  Gv(:,i,dimID) = Gv(:,i,dimID) * ed % dS(i)
               end do                  ; end do

            end if

         else
!
!           Dirichlet/Neumann boundary conditions
!           -------------------------------------
            Fa = ArtificialDissipation % ComputeFaceFluxes( ed , ed % spA % N , ed % storage(1) % Q , ed % uSB , ed % storage(1) % dQ , ed % storage(1) % dQ , ed % n )
            
            do i = 0 , ed % spA % N 
               Fa(:,i) = Fa(:,i) * ed % dS(i)
            end do

            do iXi = 0 , ed % spA % N
               select case ( ed % viscousBCType(iXi) )
   
                  case ( DIRICHLET )
   
                     Q1D  = ed % storage(1) % Q(:,iXi)
                     dQ1D = ed % storage(1) % dQ(:,iXi,:)
                     Qb1D  = ed % uSB(:,iXi)
                     Fv(:,iXi) = ViscousMethod % RiemannSolver_Dirichlet ( ed , ed % spA % N , ed % invh , Q1D , dQ1D , Qb1D , ed % n(IX:IY,iXi) ) * ed % dS(iXi)

                     if ( ViscousMethod % computeRiemannGradientFluxes ) then
                        Gv(:,iXi,:) = ViscousMethod % GradientRiemannSolver_BoundaryCondition( ed , Q1D , Qb1D , ed % n(IX:IY,iXi) ) * ed % dS(iXi)
                     end if

                  case ( NEUMANN )

                     Fv(:,iXi)   = 0.0_RP
                     Fa(:,iXi)   = 0.0_RP
                     Gv(:,iXi,:) = 0.0_RP
         
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
         integer                       :: i , l 

!$omp do private(l,i,ed,eq,N,e) schedule(runtime)
         do eID = 1 , mesh % no_of_elements

            e => mesh % elements(eID) 
            N => e % spA % N
!
!           Prolong the BOTTOM edge
!           -----------------------   
            ed => e % edges(EBOTTOM) % f
            ed % storage(e % quadPosition(EBOTTOM) ) % Q = 0.0_RP
            do l = 0 , e % spA % N ; do i = 0 , e % spA % N
               ed % storage(e % quadPosition(EBOTTOM)) % Q(:,i) = ed % storage(e % quadPosition(EBOTTOM)) % Q(:,i) + e % Q(:,i,l) * e % spA % lb(l,LEFT)
            end do                 ; end do
!
!           Prolong the RIGHT edge
!           ----------------------   
            ed => e % edges(ERIGHT) % f
            ed % storage(e % quadPosition(ERIGHT) ) % Q = 0.0_RP
            do i = 0 , e % spA % N ; do l = 0 , e % spA % N
               ed % storage(e % quadPosition(ERIGHT)) % Q(:,i) = ed % storage(e % quadPosition(ERIGHT)) % Q(:,i) + e % Q(:,l,i) * e % spA % lb(l,RIGHT) 
            end do                 ; end do
!
!           Prolong the TOP edge
!           --------------------   
            ed => e % edges(ETOP) % f
            ed % storage(e % quadPosition(ETOP) ) % Q = 0.0_RP
            do l = 0 , e % spA % N ; do i = 0 , e % spA % N
               ed % storage(e % quadPosition(ETOP)) % Q(:,i) = ed % storage(e % quadPosition(ETOP)) % Q(:,i) + e % Q(:,i,l) * e % spA % lb(l,RIGHT) 
            end do                 ; end do
!
!           Prolong the LEFT edge
!           ---------------------   
            ed => e % edges(ELEFT) % f
            ed % storage(e % quadPosition(ELEFT) ) % Q = 0.0_RP
            do i = 0 , e % spA % N ; do l = 0 , e % spA % N
               ed % storage(e % quadPosition(ELEFT)) % Q(:,i) = ed % storage(e % quadPosition(ELEFT)) % Q(:,i) + e % Q(:,l,i) * e % spA % lb(l,LEFT)
            end do                 ; end do
        
         end do
!$omp end do
            
      end subroutine DGSpatial_interpolateSolutionToBoundaries

#ifdef NAVIER_STOKES
      subroutine DGSpatial_interpolateGradientsToBoundaries( mesh )
         use Physics
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(QuadMesh_t) :: mesh
!        --------------------------------------------------------------------
         integer                       :: eID , eq , dimID
         class(QuadElement_t), pointer :: e
         class(Edge_t), pointer        :: ed
         integer, pointer              :: N
         integer                       :: i , l

!$omp do private(l,i,ed,eq,N,e,dimID) schedule(runtime)
         do eID = 1 , mesh % no_of_elements

            e => mesh % elements(eID) 
            N => e % spA % N
!
!           Prolong the BOTTOM edge
!           -----------------------   
            ed => e % edges(EBOTTOM) % f
            ed % storage(e % quadPosition(EBOTTOM) ) % dQ = 0.0_RP
            do dimID = 1 , NDIM  ; do l = 0 , e % spA % N ; do i = 0 , e % spA % N
               ed % storage(e % quadPosition(EBOTTOM)) % dQ(:,i,dimID) = ed % storage(e % quadPosition(EBOTTOM)) % dQ(:,i,dimID) + e % dQ(:,i,l,dimID) * e % spA % lb(l,LEFT)
            end do               ; end do                 ; end do
!
!           Prolong the RIGHT edge
!           ----------------------   
            ed => e % edges(ERIGHT) % f
            ed % storage(e % quadPosition(ERIGHT) ) % dQ = 0.0_RP
            do dimID = 1 , NDIM  ; do i = 0 , e % spA % N ; do l = 0 , e % spA % N
               ed % storage(e % quadPosition(ERIGHT)) % dQ(:,i,dimID) = ed % storage(e % quadPosition(ERIGHT)) % dQ(:,i,dimID) + e % dQ(:,l,i,dimID) * e % spA % lb(l,RIGHT) 
            end do               ; end do                 ; end do 
!
!           Prolong the TOP edge
!           --------------------   
            ed => e % edges(ETOP) % f
            ed % storage(e % quadPosition(ETOP) ) % dQ = 0.0_RP
            do dimID = 1 , NDIM ; do l = 0 , e % spA % N ; do i = 0 , e % spA % N
               ed % storage(e % quadPosition(ETOP)) % dQ(:,i,dimID) = ed % storage(e % quadPosition(ETOP)) % dQ(:,i,dimID) + e % dQ(:,i,l,dimID) * e % spA % lb(l,RIGHT) 
            end do              ; end do                 ; end do
!
!           Prolong the LEFT edge
!           ---------------------   
            ed => e % edges(ELEFT) % f
            ed % storage(e % quadPosition(ELEFT) ) % dQ = 0.0_RP
            do dimID = 1 , NDIM ; do i = 0 , e % spA % N ; do l = 0 , e % spA % N
               ed % storage(e % quadPosition(ELEFT)) % dQ(:,i,dimID) = ed % storage(e % quadPosition(ELEFT)) % dQ(:,i,dimID) + e % dQ(:,l,i,dimID) * e % spA % lb(l,LEFT)
            end do              ; end do                 ; end do
        
         end do
!$omp end do
            
      end subroutine DGSpatial_interpolateGradientsToBoundaries
#endif

end module DGSpatialDiscretizationMethods
