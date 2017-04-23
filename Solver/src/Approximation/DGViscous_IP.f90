!
#ifdef NAVIER_STOKES
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!        This submodule contains the INTERIOR PENALTY scheme procedures
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
submodule (DGViscousMethods)  DGViscous_IP
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
      module subroutine IP_ComputeGradient( self , mesh ) 
         use DGWeakIntegrals
         implicit none
         class(IPMethod_t)   ,  intent (in)    :: self
         class(QuadMesh_t)    , intent (inout) :: mesh
!
!        ---------------
!        Local variables   
!        ---------------
!
         integer     :: eID , edID , iDim , eq
         real(kind=RP), pointer  :: dQ(:,:,:,:)
!
!        Just compute the interior gradient
!        ----------------------------------
         do eID = 1 , mesh % no_of_elements
            call mesh % elements(eID) % ComputeInteriorGradient
         end do

      end subroutine IP_ComputeGradient

      module pure function IP_ComputeInnerFluxes( self , e ) result (Fv)
!
!        ***************************************************************************************
!              The IP method computes the following fluxes:
!
!                 Fv = viscousFluxes( Q , dQ )
!
!           on the set of interpolation nodes. Q is the solution, and dQ is the solution gradient
!           computed with the extended system weak formulation.
!        ****************************************************************************************
!           
         use QuadElementClass
         implicit none
         class(IPMethod_t),   intent(in)   :: self
         class(QuadElement_t), intent(in)   :: e
         real(kind=RP)                      :: Fv(0 : e % spA % N , 0 : e % spA % N , 1:NCONS , 1:NDIM)
!
!        Compute the cartesian flux
!        --------------------------
         Fv = ViscousFlux( e % spA % N , e % Q , e % dQ)

      end function IP_ComputeInnerFluxes

      module pure function IP_RiemannSolver( self , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
!
!        *****************************************************************************************
!              The IP Viscous Riemann Solver averages the fluxes obtained from both sides.
!           UL and UR are the LEFT and RIGHT solutions respectively, whilst dUL and dUR are 
!           the LEFT and RIGHT solution gradients computed from the extended system weak 
!           formulation. Therefore, it results in a second order stencil.
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
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
         real(kind=RP)     :: penalty
!
!        Compute the penalty parameter
!        -----------------------------
         penalty = self % sigma0 * dimensionless % mu * N * N * invh_edge
         
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

         Fstar = Fstar - penalty * ( uL - uR )


      end function IP_RiemannSolver

      module pure function IP_RiemannSolver_Dirichlet( self , N , invh_edge , u , g , uB , normal ) result ( Fstar )
!
!        *****************************************************************************************
!              For the Dirichlet boundary conditions, the IP scheme uses the interior values
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)     :: self
         integer,                intent(in)     :: N
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
         real(kind=RP)  :: penalty
!
!        Compute the two dimensional flux
!        --------------------------------
         F = ViscousFlux( u , g ) 
!
!        Projection along the boundary normal
!        ------------------------------------
         FStar =  F(:,IX) * normal(IX) + F(:,IY) * normal(IY) 

         penalty = self % sigma0 * dimensionless % mu * N * N * invh_edge
         Fstar = Fstar - penalty * ( u - uB )

      end function IP_RiemannSolver_Dirichlet

      module pure function IP_RiemannSolver_Adiabatic( self , N , invh_edge , u , g , uB , normal ) result ( Fstar )
!
!        ************************************************************************************************
!              For the Adiabatic dirichlet boundary conditions, the IP scheme uses the interior values
!        ************************************************************************************************
!
         implicit none
         class(IPMethod_t), intent (in) :: self
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
         real(kind=RP)     :: penalty
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

         penalty = self % sigma0 * dimensionless % mu * N * N * invh_edge
         Fstar = Fstar - penalty * ( u - uB )
 
      end function IP_RiemannSolver_Adiabatic

      module pure subroutine IP_GradientRiemannSolver( self , N , UL , UR , normal , GstarL , GstarR ) 
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP), intent(out)       :: GstarL(0:N , 1:NCONS , 1:NDIM)
         real(kind=RP), intent(out)       :: GstarR(0:N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: falseGradient (0:N , 1:NDIM , 1:NCONS)
         integer              :: eq
!
!        The false gradients are the interfaces jumps (pointing from L -> R)
!        -------------------------------------------------------------------
         do eq = 1 , NCONS
            falseGradient(:,IX,eq) = (UL(:,eq) - UR(:,eq)) * normal(IX,:)
            falseGradient(:,IY,eq) = (UL(:,eq) - UR(:,eq)) * normal(IY,:)
         end do
!
!        Compute the fluxes
!        ------------------
         GStarL =  0.5_RP * self % epsilon * ViscousFlux( N , uL , falseGradient )
         GStarR = -0.5_RP * self % epsilon * ViscousFlux( N , uR , falseGradient )

      end subroutine IP_GradientRiemannSolver

      module pure function IP_GradientRiemannSolver_BoundaryCondition( self , u , uB , normal ) result ( Gstar ) 
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
         real(kind=RP), intent(in)        :: u(1:NCONS)
         real(kind=RP), intent(in)        :: uB(1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY)
         real(kind=RP)                    :: Gstar(1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: falseGradient (1:NDIM , 1:NCONS)
         integer              :: eq
!
!        The false gradients are the interfaces jumps (pointing from L -> R)
!        -------------------------------------------------------------------
         falseGradient(IX,:) = (u - uB) * normal(IX)
         falseGradient(IY,:) = (u - uB) * normal(IY)
!
!        Compute the fluxes
!        ------------------
         GStar = 1.0_RP * self % epsilon * ViscousFlux( u , falseGradient )

      end function IP_GradientRiemannSolver_BoundaryCondition

      module pure function IP_GradientRiemannSolver_Adiabatic( self , N , u , uB , normal )  result ( Gstar )
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: u(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uB(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Gstar(0:N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: falseGradient (0:N , 1:NDIM , 1:NCONS)
         integer              :: eq
!
!        The false gradients are the interfaces jumps (pointing from L -> R)
!        -------------------------------------------------------------------
         do eq = 1 , NCONS
            falseGradient(:,IX,eq) = (u(:,eq) - uB(:,eq)) * normal(IX,:)
            falseGradient(:,IY,eq) = (u(:,eq) - uB(:,eq)) * normal(IY,:)
         end do
!
!        Compute the fluxes
!        ------------------
!         GStar =  self % epsilon * AdiabaticViscousFlux( N , u , u , falseGradient )
         !GStar =  self % epsilon * ViscousFlux( N , u , falseGradient )
         Gstar = 0.0_RP

      end function IP_GradientRiemannSolver_Adiabatic
end submodule DGViscous_IP

#endif
