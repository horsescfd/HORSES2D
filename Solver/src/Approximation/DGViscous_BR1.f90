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
         real(kind=RP)                      :: Fv(0 : e % spA % N , 0 : e % spA % N , 2:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)              :: F_cartesian(0:e % spA % N,0:e % spA % N,2:NCONS,1:NDIM)
         integer                    :: eq
         integer                    :: N 

         N = e % spA % N         
!
!        Compute the cartesian flux
!        --------------------------
         F_cartesian = ViscousFlux( e % spA % N , e % Q , e % dQ)

         do eq = IRHO+1 , NCONS
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

      module pure function BR1_RiemannSolver( self , N , UL , UR , dUL , dUR , normal ) result ( FStar )
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
         real(kind=RP), intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: dUL(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: dUR(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(0:N , 2:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer           :: eq
         real(kind=RP)     :: FL(0:N , 2:NCONS , 1:NDIM)
         real(kind=RP)     :: FR(0:N , 2:NCONS , 1:NDIM)
!
!        Compute the LEFT and RIGHT viscous fluxes
!        -----------------------------------------
         FL = ViscousFlux( N , uL , duL )
         FR = ViscousFlux( N , uR , duR )
!
!        Perform the average and the projection along the edge normal
!        ------------------------------------------------------------
         do eq = 2 , NCONS
            Fstar(:,eq) = 0.5_RP * ( ( FL(:,eq,IX) + FR(:,eq,IX) ) * normal(IX,:) + ( FL(:,eq,IY) + FR(:,eq,IY) ) * normal(IY,:) )  
         end do

      end function BR1_RiemannSolver

      module pure function BR1_RiemannSolver_Dirichlet( self , u , g , uB , n ) result ( Fstar )
!
!        *****************************************************************************************
!              For the Dirichlet boundary conditions, the BR1 scheme uses the interior values
!        *****************************************************************************************
!
         implicit none
         class(BR1Method_t), intent(in)     :: self
         real(kind=RP),          intent(in)     :: u(NCONS)
         real(kind=RP),          intent(in)     :: g(NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(NCONS)
         real(kind=RP),          intent(in)     :: n(NDIM)
         real(kind=RP)                          :: Fstar(2:NCONS)
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: F(2:NCONS,1:NDIM)
!
!        Compute the two dimensional flux
!        --------------------------------
         F = ViscousFlux( u , g ) 
!
!        Projection along the boundary normal
!        ------------------------------------
         FStar =  F(:,IX) * n(IX) + F(:,IY) * n(IY) 

      end function BR1_RiemannSolver_Dirichlet

      module pure function BR1_RiemannSolver_Adiabatic( self , N , u , g , uB , normal ) result ( Fstar )
!
!        ************************************************************************************************
!              For the Adiabatic dirichlet boundary conditions, the BR1 scheme uses the interior values
!        ************************************************************************************************
!
         implicit none
         class(BR1Method_t), intent (in) :: self
         integer      ,      intent (in) :: N
         real(kind=RP),      intent (in) :: u      ( 0:N  , NCONS      )
         real(kind=RP),      intent (in) :: g      ( 0:N  , NDIM,NCONS )
         real(kind=RP),      intent (in) :: uB     ( 0:N  , NCONS      )
         real(kind=RP),      intent (in) :: normal ( NDIM , 0:N        )
         real(kind=RP)                   :: Fstar  ( 0:N  , 2:NCONS    )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: F(0:N , 2:NCONS , 1:NDIM)
         integer           :: eq
!
!        Compute the Adiabatic viscous flux based on the interior points
!        ---------------------------------------------------------------
         F = AdiabaticViscousFlux( N , u , u , g)
!
!        Perform the projection along the boundary normal
!        ------------------------------------------------
         do eq = 2 , NCONS
            FStar(:,eq) = F(:,eq,IX) * normal(IX,:) + F(:,eq,IY) * normal(IY,:)
         end do
 
      end function BR1_RiemannSolver_Adiabatic

end submodule
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#endif
!
