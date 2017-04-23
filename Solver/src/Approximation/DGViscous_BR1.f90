!
#ifdef NAVIER_STOKES
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!        This submodule contains the BASSI-REBAY I scheme procedures
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
submodule (DGViscousMethods)  DGViscous_BR1
   use SMConstants
   use QuadMeshClass
   use QuadElementClass
   use Physics
   use Setup_class
   implicit none
!
#include "Defines.h"
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////
!
      module subroutine BR1_ComputeGradient( self , mesh ) 
         use DGWeakIntegrals
         implicit none
         class(BR1Method_t)   , intent(in)      :: self
         class(QuadMesh_t)    , intent(inout)   :: mesh
!
!        ---------------
!        Local variables   
!        ---------------
!
         integer     :: eID , edID , iDim , eq
         real(kind=RP), pointer  :: dQ(:,:,:,:)
!
!        Perform volume loops
!        --------------------
         do eID = 1 , mesh % no_of_elements
            dQ(0:,0:,1:,1:)   => mesh % elements(eID) % dQ
            dQ = - VectorWeakIntegrals % StdVolumeGreen( mesh % elements(eID) , mesh % elements(eID) % Q )
            
         end do
!
!        Perform face loops
!        ------------------
         do edID = 1 , mesh % no_of_edges
            select type ( f => mesh % edges(edID) % f ) 
               type is (Edge_t)
                  call BR1_dQFaceLoop_Interior(self , f)

               type is (StraightBdryEdge_t)
                  call BR1_dQFaceLoop_StraightBdry(self , f)

               type is (CurvedBdryEdge_t)
                  call BR1_dQFaceLoop_CurvedBdry(self , f)

            end select
         end do
!
!        Perform the scaling with the jacobian
!        -------------------------------------
         do eID = 1 , mesh % no_of_elements
            do iDim = 1 , NDIM   ; do eq = 1 , NCONS
               mesh % elements(eID) % dQ(:,:,iDim,eq) = mesh % elements(eID) % dQ(:,:,iDim,eq) / mesh % elements(eID) % jac
            end do               ; end do
         end do

      end subroutine BR1_ComputeGradient 

      subroutine BR1_dQFaceLoop_Interior( ViscousMethod , ed )
         use QuadElementClass
         use DGWeakIntegrals
         implicit none
         class(ViscousMethod_t)        :: ViscousMethod
         type(Edge_t)                  :: ed
         real(kind=RP)                 :: UstarL( 0 : ed % storage(LEFT ) % spA % N , 1:NCONS )
         real(kind=RP)                 :: UstarR( 0 : ed % storage(RIGHT) % spA % N , 1:NCONS )
         real(kind=RP), pointer        :: dQ(:,:,:,:)

         call ViscousMethod % ComputeSolutionRiemann( ed , UstarL , UstarR )
!
!>       Add the contribution to the LEFT element
!        ----------------------------------------
         dQ(0:,0:,1:,1:)   => ed % quads(LEFT) % e % dQ
         dQ = dQ + VectorWeakIntegrals % StdFace( ed , LEFT , UstarL )
!
!>       Add the contribution to the RIGHT element
!        -----------------------------------------
         dQ(0:,0:,1:,1:)   => ed % quads(RIGHT) % e % dQ
         dQ = dQ + VectorWeakIntegrals % StdFace( ed , RIGHT , UstarR ) 

      end subroutine BR1_dQFaceLoop_Interior

      subroutine BR1_dQFaceLoop_StraightBdry( ViscousMethod , ed )
         use QuadElementClass
         use DGWeakIntegrals
         implicit none
         class(ViscousMethod_t)           :: ViscousMethod
         type(StraightBdryEdge_t)         :: ed
         real(kind=RP)                    :: Ustar( 0 : ed % spA % N , 1:NCONS)
         real(kind=RP), pointer           :: dQ(:,:,:,:)
      
         call ViscousMethod % ComputeSolutionRiemann( ed , Ustar ) 
!
!>       Add the contribution to the element
!        -----------------------------------
         dQ(0:,0:,1:,1:)   => ed % quads(1) % e % dQ
         dQ = dQ + VectorWeakIntegrals % StdFace( ed , 1 , Ustar )
!
      end subroutine BR1_dQFaceLoop_StraightBdry

      subroutine BR1_dQFaceLoop_CurvedBdry( ViscousMethod , ed )
         use QuadElementClass
         use DGWeakIntegrals
         implicit none
         class(ViscousMethod_t)      :: ViscousMethod
         type(CurvedBdryEdge_t)      :: ed
         real(kind=RP), pointer        :: dQ(:,:,:,:)
!
!>       Add the boundary contribution to the element
!        --------------------------------------------
         dQ(0:,0:,1:,1:)   => ed % quads(1) % e % dQ
         dQ = dQ + VectorWeakIntegrals % StdFace( ed , 1 , ed % uSB )
!
      end subroutine BR1_dQFaceLoop_CurvedBdry         

      module pure function BR1_ComputeInnerFluxes( self , e ) result (Fv)
!
!        ***************************************************************************************
!              The BR1 method computes the following fluxes:
!
!                 Fv = viscousFluxes( Q , dQ )
!
!           on the set of interpolation nodes. Q is the solution, and dQ is the solution gradient
!           computed with the extended system weak formulation.
!        ****************************************************************************************
!           
         use QuadElementClass
         implicit none
         class(BR1Method_t),   intent(in)   :: self
         class(QuadElement_t), intent(in)   :: e
         real(kind=RP)                      :: Fv(0 : e % spA % N , 0 : e % spA % N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)              :: F_cartesian(0:e % spA % N,0:e % spA % N,1:NCONS,1:NDIM)
         integer                    :: eq
         integer                    :: N 

         N = e % spA % N         
!
!        Compute the cartesian flux
!        --------------------------
         F_cartesian = ViscousFlux( e % spA % N , e % Q , e % dQ)

         do eq = 1 , NCONS
!           
!           F flux (contravariant)
!           ----------------------
            Fv(0:N,0:N,eq,IX) = F_cartesian(0:N,0:N,eq,IX) * e % Ja(0:N,0:N,1,1) + F_cartesian(0:N,0:N,eq,IY) * e % Ja(0:N,0:N,2,1)
!           
!           G flux (contravariant)
!           ----------------------
            Fv(0:N,0:N,eq,IY) = F_cartesian(0:N,0:N,eq,IX) * e % Ja(0:N,0:N,1,2) + F_cartesian(0:N,0:N,eq,IY) * e % Ja(0:N,0:N,2,2)
         end do


      end function BR1_ComputeInnerFluxes

      module pure function BR1_SolutionRiemannSolver( self , N , UL , UR ) result ( uStar )
!
!        *************************************************************************
!              The BR1 Solution Riemann solver simply averages the two states
!        *************************************************************************
!
         implicit none
         class(BR1Method_t), intent(in) :: self
         integer, intent(in)            :: N
         real(kind=RP), intent(in)      :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)      :: uR(0:N , 1:NCONS)
         real(kind=RP)                  :: uStar(0:N,1:NCONS)
!
!        Perform the average
!        -------------------
         uStar = 0.5_RP * ( uL + uR )

      end function BR1_SolutionRiemannSolver

      module pure function BR1_RiemannSolver( self , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
!
!        *****************************************************************************************
!              The BR1 Viscous Riemann Solver averages the fluxes obtained from both sides.
!           UL and UR are the LEFT and RIGHT solutions respectively, whilst dUL and dUR are 
!           the LEFT and RIGHT solution gradients computed from the extended system weak 
!           formulation. Therefore, it results in a second order stencil.
!        *****************************************************************************************
!
         implicit none
         class(BR1Method_t), intent(in)   :: self
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: invh_edge
         real(kind=RP), intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: dUL(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: dUR(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(0:N , 1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer           :: eq
         real(kind=RP)     :: FL(0:N , 1:NCONS , 1:NDIM)
         real(kind=RP)     :: FR(0:N , 1:NCONS , 1:NDIM)
!
!        Compute the LEFT and RIGHT viscous fluxes
!        -----------------------------------------
         FL = ViscousFlux( N , uL , duL )
         FR = ViscousFlux( N , uR , duR )
!
!        Perform the average and the projection along the edge normal
!        ------------------------------------------------------------
         do eq = 1 , NCONS
            Fstar(:,eq) = 0.5_RP * ( ( FL(:,eq,IX) + FR(:,eq,IX) ) * normal(IX,:) + ( FL(:,eq,IY) + FR(:,eq,IY) ) * normal(IY,:) )  
         end do

      end function BR1_RiemannSolver

      module pure function BR1_RiemannSolver_Dirichlet( self , N , invh_edge , u , g , uB , normal ) result ( Fstar )
!
!        *****************************************************************************************
!              For the Dirichlet boundary conditions, the BR1 scheme uses the interior values
!        *****************************************************************************************
!
         implicit none
         class(BR1Method_t), intent(in)     :: self
         integer      ,          intent(in)     :: N 
         real(kind=RP),          intent(in)     :: invh_edge
         real(kind=RP),          intent(in)     :: u(NCONS)
         real(kind=RP),          intent(in)     :: g(NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(NCONS)
         real(kind=RP),          intent(in)     :: normal(NDIM)
         real(kind=RP)                          :: Fstar(1:NCONS)
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: F(1:NCONS,1:NDIM)
!
!        Compute the two dimensional flux
!        --------------------------------
         F = ViscousFlux( u , g ) 
!
!        Projection along the boundary normal
!        ------------------------------------
         FStar =  F(:,IX) * normal(IX) + F(:,IY) * normal(IY) 

      end function BR1_RiemannSolver_Dirichlet

      module pure function BR1_RiemannSolver_Adiabatic( self , N , invh_edge , u , g , uB , normal ) result ( Fstar )
!
!        ************************************************************************************************
!              For the Adiabatic dirichlet boundary conditions, the BR1 scheme uses the interior values
!        ************************************************************************************************
!
         implicit none
         class(BR1Method_t), intent (in) :: self
         integer      ,      intent (in) :: N
         real(kind=RP),      intent (in) :: invh_edge
         real(kind=RP),      intent (in) :: u      ( 0:N  , NCONS      )
         real(kind=RP),      intent (in) :: g      ( 0:N  , NDIM,NCONS )
         real(kind=RP),      intent (in) :: uB     ( 0:N  , NCONS      )
         real(kind=RP),      intent (in) :: normal ( NDIM , 0:N        )
         real(kind=RP)                   :: Fstar  ( 0:N  , 1:NCONS    )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: F(0:N , 1:NCONS , 1:NDIM)
         integer           :: eq
!
!        Compute the Adiabatic viscous flux based on the interior points
!        ---------------------------------------------------------------
         F = AdiabaticViscousFlux( N , u , u , g)
!
!        Perform the projection along the boundary normal
!        ------------------------------------------------
         do eq = 1 , NCONS
            FStar(:,eq) = F(:,eq,IX) * normal(IX,:) + F(:,eq,IY) * normal(IY,:)
         end do
 
      end function BR1_RiemannSolver_Adiabatic

end submodule
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#endif
!
