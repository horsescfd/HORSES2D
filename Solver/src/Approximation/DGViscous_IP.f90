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
!
!        Just compute the interior gradient
!        ----------------------------------
         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % dQ = mesh % elements(eID) % ComputeInteriorGradient()
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
         real(kind=RP)                      :: Fv(1 : NCONS , 0 : e % spA % N , 0 : e % spA % N , 1:NDIM)
!
!        Compute the cartesian flux
!        --------------------------
         Fv = (dimensionless % mu + e % mu_a) * ViscousFlux( e % spA % N , e % Q , e % dQ)

      end function IP_ComputeInnerFluxes

      module pure function IP_RiemannSolver( self , edge , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
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
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: invh_edge
         real(kind=RP), intent(in)        :: uL(1:NCONS , 0:N )
         real(kind=RP), intent(in)        :: uR(1:NCONS , 0:N )
         real(kind=RP), intent(in)        :: dUL(1:NCONS , 0:N, 1:NDIM )
         real(kind=RP), intent(in)        :: dUR(1:NCONS , 0:N, 1:NDIM )
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(1:NCONS , 0:N )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer           :: dimID , i
         real(kind=RP)     :: FL(1:NCONS , 0:N , 1:NDIM)
         real(kind=RP)     :: FR(1:NCONS , 0:N , 1:NDIM)
         real(kind=RP)     :: penalty
!
!        Compute the penalty parameter
!        -----------------------------
         penalty = self % sigma0 * (dimensionless % mu + edge % mu_a) * N * N * invh_edge
!
!        Compute the LEFT and RIGHT viscous fluxes
!        -----------------------------------------
         FL = (dimensionless % mu + edge % mu_a) * ViscousFlux( N , uL , duL )
         FR = (dimensionless % mu + edge % mu_a) * ViscousFlux( N , uR , duR )
!
!        Perform the average and the projection along the edge normal
!        ------------------------------------------------------------
         Fstar = 0.0_RP
         do dimID = 1 , NDIM     ; do i = 0 , N
            Fstar(:,i) = Fstar(:,i) + 0.5_RP * ( FL(:,i,dimID) + FR(:,i,dimID) ) * normal(dimID,i) 
         end do                  ; end do

         Fstar = Fstar - penalty * ( uL - uR )

      end function IP_RiemannSolver

      module pure function IP_RiemannSolver_Dirichlet( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
!
!        *****************************************************************************************
!              For the Dirichlet boundary conditions, the IP scheme uses the interior values
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)     :: self
         class(Edge_t)    , intent(in)          :: edge
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
         F = (dimensionless % mu + edge % mu_a) * ViscousFlux( u , g ) 
!
!        Projection along the boundary normal
!        ------------------------------------
         FStar =  F(:,IX) * normal(IX) + F(:,IY) * normal(IY) 

         penalty = self % sigma0 * (dimensionless % mu + edge % mu_a) * N * N * invh_edge
         Fstar = Fstar - penalty * ( u - uB )

      end function IP_RiemannSolver_Dirichlet

      module pure function IP_RiemannSolver_Adiabatic( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
!
!        ************************************************************************************************
!              For the Adiabatic dirichlet boundary conditions, the IP scheme uses the interior values
!        ************************************************************************************************
!
         implicit none
         class(IPMethod_t), intent (in)  :: self
         class(Edge_t)    ,  intent(in)  :: edge
         integer      ,      intent (in) :: N
         real(kind=RP),      intent (in) :: invh_edge
         real(kind=RP) , intent (in)     :: u      ( 1 : NCONS , 0 : N                )
         real(kind=RP) , intent (in)     :: g      ( 1 : NCONS , 0 : N     , 1 : NDIM )
         real(kind=RP) , intent (in)     :: uB     ( 1 : NCONS , 0 : N                )
         real(kind=RP) , intent (in)     :: normal ( 1 : NDIM  , 0 : N                )
         real(kind=RP)                   :: Fstar  ( 1 : NCONS , 0 : N                )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: F(1 : NCONS , 0:N , 1:NDIM)
         real(kind=RP)     :: penalty 
         integer           :: dimID , i
!
!        Compute the Adiabatic viscous flux based on the interior points
!        ---------------------------------------------------------------
         F = (dimensionless % mu + edge % mu_a) * AdiabaticViscousFlux( N , u , u , g)
!
!        Perform the projection along the boundary normal
!        ------------------------------------------------
         Fstar = 0.0_RP
         do dimID = 1 , NDIM  ; do i = 0 , N
            FStar(:,i) = Fstar(:,i) + F(:,i,dimID) * normal(dimID,i)
         end do               ; end do

         penalty = self % sigma0 * (dimensionless % mu + edge % mu_a) * N * N * invh_edge
         Fstar = Fstar - penalty * ( u - uB )
 
      end function IP_RiemannSolver_Adiabatic

      module pure subroutine IP_GradientRiemannSolver( self , edge , N , UL , UR , normal , GstarL , GstarR ) 
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: uL(1:NCONS , 0:N )
         real(kind=RP), intent(in)        :: uR(1:NCONS , 0:N )
         real(kind=RP), intent(in)        :: normal(IX:IY , 0:N)
         real(kind=RP), intent(out)       :: GstarL(1:NCONS , 0:N  , 1:NDIM)
         real(kind=RP), intent(out)       :: GstarR(1:NCONS , 0:N  , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: falseGradient (1:NCONS , 0:N , 1:NDIM)
         integer              :: i , dimID
!
!        The false gradients are the interfaces jumps (pointing from L -> R)
!        -------------------------------------------------------------------
         do dimID = 1 , NDIM  ; do i = 0 , N
            falseGradient(:,i,dimID) = (UL(:,i) - UR(:,i)) * normal(dimID,i)
         end do               ; end do
!
!        Compute the fluxes
!        ------------------
         GStarL =  0.5_RP * self % epsilon * (dimensionless % mu + edge % mu_a) * ViscousFlux( N , uL , falseGradient )
         GStarR = -0.5_RP * self % epsilon * (dimensionless % mu + edge % mu_a) * ViscousFlux( N , uR , falseGradient )

      end subroutine IP_GradientRiemannSolver

      module pure function IP_GradientRiemannSolver_BoundaryCondition( self , edge , u , uB , normal ) result ( Gstar ) 
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         real(kind=RP), intent(in)        :: u(1:NCONS)
         real(kind=RP), intent(in)        :: uB(1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY)
         real(kind=RP)                    :: Gstar(1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: falseGradient (1:NCONS , 1:NDIM)
         integer              :: eq
!
!        The false gradients are the interfaces jumps (pointing from L -> R)
!        -------------------------------------------------------------------
         falseGradient(:,IX) = (u - uB) * normal(IX)
         falseGradient(:,IY) = (u - uB) * normal(IY)
!
!        Compute the fluxes
!        ------------------
         GStar = 1.0_RP * self % epsilon * (dimensionless % mu + edge % mu_a) * ViscousFlux( u , falseGradient )

      end function IP_GradientRiemannSolver_BoundaryCondition

      module pure function IP_GradientRiemannSolver_Adiabatic( self , edge , N , u , uB , normal )  result ( Gstar )
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: u(1:NCONS , 0:N )
         real(kind=RP), intent(in)        :: uB(1:NCONS , 0:N )
         real(kind=RP), intent(in)        :: normal(IX:IY, 0:N)
         real(kind=RP)                    :: Gstar(1:NCONS , 0:N  , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: falseGradient (1:NCONS,0:N , 1:NDIM)
         integer              :: i , dimID
!
!        The false gradients are the interfaces jumps (pointing from L -> R)
!        -------------------------------------------------------------------
         do dimID = 1 , NDIM  ; do i = 0 , N
            falseGradient(:,i,dimID) = (u(:,i) - uB(:,i)) * normal(dimID,i)
         end do               ; end do
!
!        Compute the fluxes
!        ------------------
         GStar =  self % epsilon * (dimensionless % mu + edge % mu_a) * AdiabaticViscousFlux( N , u , u , falseGradient )

      end function IP_GradientRiemannSolver_Adiabatic
end submodule DGViscous_IP

#endif
