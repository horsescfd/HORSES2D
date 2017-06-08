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
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           DGViscous procedures
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#ifdef NAVIER_STOKES
#include "Defines.h"
!
module DGViscousMethods
   use SMConstants
   use QuadMeshClass
   use QuadElementClass
   use Physics
   use Setup_class
   implicit none
!
!  *******
   private
   public  ViscousMethod_t , IPMethod_t , BR1Method_t , LDGMethod_t
   public  ViscousMethod_Initialization
!  *******
!
!                                *************************
   integer, parameter         :: STR_LEN_VISCOUS = 128
!                                *************************
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!                --------------------
!                | TYPE DEFINITIONS |
!                --------------------
!
!  ---------------------------------------------------------------
!        VISCOUS METHOD GENERIC TYPE 
!  ---------------------------------------------------------------
!
   type ViscousMethod_t
      logical                                 :: computeRiemannGradientFluxes
      character(len=STR_LEN_VISCOUS)          :: method
      contains
!
!                                ***************************************
!                                   GRADIENTS PROCEDURE
!                                ***************************************
!
         procedure          ::   ComputeGradient                        => BaseClass_ComputeGradient
!
!                                ***************************************
!                                   INNER FLUXES PROCEDURE
!                                ***************************************
!
         procedure          ::   ComputeInnerFluxes                     => BaseClass_ComputeInnerFluxes
!
!                                ***************************************
!                                   SOLUTION RIEMANN PROCEDURE
!                                ***************************************
!
         generic, public    ::   ComputeSolutionRiemann                  => ComputeSolutionRiemann_Interior       , &
                                                                            ComputeSolutionRiemann_CurvedInterior , &
                                                                            ComputeSolutionRiemann_StraightBdry   , &
                                                                            ComputeSolutionRiemann_CurvedBdry   
         generic, public    ::   ComputeSubdividedSolutionRiemann        => ComputeSolutionRiemann_Subdivided     , &
                                                                            ComputeSolutionRiemann_CurvedSubdivided
         procedure, private ::   ComputeSolutionRiemann_Interior         => BaseClass_ComputeSolutionRiemann_Interior
         procedure, private ::   ComputeSolutionRiemann_CurvedInterior   => BaseClass_ComputeSolutionRiemann_CurvedInterior
         procedure, private ::   ComputeSolutionRiemann_Subdivided       => BaseClass_ComputeSolutionRiemann_Subdivided
         procedure, private ::   ComputeSolutionRiemann_CurvedSubdivided => BaseClass_ComputeSolutionRiemann_CurvedSubdivided
         procedure, private ::   ComputeSolutionRiemann_StraightBdry     => BaseClass_ComputeSolutionRiemann_StraightBdry
         procedure, private ::   ComputeSolutionRiemann_CurvedBdry       => BaseClass_ComputeSolutionRiemann_CurvedBdry
         procedure, private ::   SolutionRiemannSolver                   => BaseClass_SolutionRiemannSolver
!
!                                ***************************************
!                                    RIEMANN SOLVER PROCEDURE
!                                ***************************************
!
         procedure          ::   RiemannSolver                           => BaseClass_RiemannSolver
         procedure          ::   RiemannSolver_Dirichlet                 => BaseClass_RiemannSolver_Dirichlet
         procedure          ::   RiemannSolver_Adiabatic                 => BaseClass_RiemannSolver_Adiabatic
!
!                                ***************************************
!                                    GRADIENT RIEMANN SOLVER PROCEDURE
!                                ***************************************
!
         procedure          ::   GradientRiemannSolver                   => BaseClass_GradientRiemannSolver
         procedure          ::   GradientRiemannSolver_BoundaryCondition => BaseClass_GradientRiemannSolver_BoundaryCondition
         procedure          ::   GradientRiemannSolver_Adiabatic         => BaseClass_GradientRiemannSolver_Adiabatic
         procedure          ::   Describe                                => ViscousMethod_describe
   end type ViscousMethod_t
!
!  ---------------------------------------------------------------
!        INTERIOR PENALTY METHOD
!  ---------------------------------------------------------------
!
   type, extends(ViscousMethod_t) :: IPMethod_t
      character(len=STR_LEN_VISCOUS) :: subType
      real(kind=RP)                  :: sigma0
      real(kind=RP)                  :: sigma1
      real(kind=RP)                  :: epsilon
      contains
         procedure   :: ComputeGradient                         => IP_ComputeGradient
         procedure   :: ComputeInnerFluxes                      => IP_ComputeInnerFluxes
         procedure   :: RiemannSolver                           => IP_RiemannSolver
         procedure   :: RiemannSolver_Dirichlet                 => IP_RiemannSolver_Dirichlet
         procedure   :: RiemannSolver_Adiabatic                 => IP_RiemannSolver_Adiabatic
         procedure   :: GradientRiemannSolver                   => IP_GradientRiemannSolver
         procedure   :: GradientRiemannSolver_BoundaryCondition => IP_GradientRiemannSolver_BoundaryCondition
         procedure   :: GradientRiemannSolver_Adiabatic         => IP_GradientRiemannSolver_Adiabatic
   end type IPMethod_t
!
!  ---------------------------------------------------------------
!        BASSI-REBAY I METHOD
!  ---------------------------------------------------------------
!
   type, extends(ViscousMethod_t) ::  BR1Method_t
      contains
         procedure          :: ComputeGradient         => BR1_ComputeGradient
         procedure          :: ComputeInnerFluxes      => BR1_ComputeInnerFluxes
         procedure, private :: SolutionRiemannSolver   => BR1_SolutionRiemannSolver
         procedure          :: RiemannSolver           => BR1_RiemannSolver
         procedure          :: RiemannSolver_Dirichlet => BR1_RiemannSolver_Dirichlet
         procedure          :: RiemannSolver_Adiabatic => BR1_RiemannSolver_Adiabatic
   end type BR1Method_t
!
!  ---------------------------------------------------------------
!        LOCAL DISCONTINUOUS GALERKIN (LDG) METHOD
!  ---------------------------------------------------------------
!
   type, extends(ViscousMethod_t) ::  LDGMethod_t

   end type LDGMethod_t
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!        
!
!  --------------
!  BR1 Interfaces
!  --------------
!
   interface
      module subroutine BR1_ComputeGradient( self , mesh ) 
         implicit none
         class(BR1Method_t),    intent(in)     :: self
         class(QuadMesh_t) ,    intent(inout)  :: mesh
      end subroutine BR1_ComputeGradient

      module pure function BR1_ComputeInnerFluxes( self , e ) result (Fv)
         use QuadElementClass
         implicit none
         class(BR1Method_t),   intent(in)   :: self
         class(QuadElement_t), intent(in)   :: e
         real(kind=RP)                      :: Fv(1:NCONS , 0 : e % spA % N , 0 : e % spA % N , 1:NDIM)
      end function BR1_ComputeInnerFluxes

      module pure function BR1_SolutionRiemannSolver( self , N , UL , UR ) result ( uStar )
         implicit none
         class(BR1Method_t), intent(in) :: self
         integer, intent(in)            :: N
         real(kind=RP), intent(in)      :: uL(1:NCONS , 0:N)
         real(kind=RP), intent(in)      :: uR(1:NCONS , 0:N)
         real(kind=RP)                  :: uStar(1:NCONS , 0:N)
      end function BR1_SolutionRiemannSolver

      module pure function BR1_RiemannSolver( self , edge , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
         implicit none
         class(BR1Method_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: invh_edge
         real(kind=RP), intent(in)        :: uL(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: uR(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: dUL(1:NCONS , 0:N, 1:NDIM)
         real(kind=RP), intent(in)        :: dUR(1:NCONS , 0:N, 1:NDIM)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(1:NCONS , 0:N)
      end function BR1_RiemannSolver

      module pure function BR1_RiemannSolver_Dirichlet( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
         implicit none
         class(BR1Method_t), intent(in)     :: self
         class(Edge_t)    , intent(in)          :: edge
         integer      ,          intent(in)     :: N
         real(kind=RP),          intent(in)     :: invh_edge
         real(kind=RP),          intent(in)     :: u(NCONS)
         real(kind=RP),          intent(in)     :: g(NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(NCONS)
         real(kind=RP),          intent(in)     :: normal(NDIM)
         real(kind=RP)                          :: Fstar(1:NCONS)
      end function BR1_RiemannSolver_Dirichlet

      module pure function BR1_RiemannSolver_Adiabatic( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
         implicit none
         class(BR1Method_t), intent(in)     :: self
         integer      ,          intent(in)     :: N 
         class(Edge_t)         , intent(in)     :: edge
         real(kind=RP),          intent(in)     :: invh_edge
         real(kind=RP),          intent(in)     :: u(1:NCONS , 0:N)
         real(kind=RP),          intent(in)     :: g(1:NCONS , 0:N , 1:NDIM)
         real(kind=RP),          intent(in)     :: uB(1:NCONS , 0:N)
         real(kind=RP),          intent(in)     :: normal(NDIM , 0:N)
         real(kind=RP)                          :: Fstar(1:NCONS , 0:N)
      end function BR1_RiemannSolver_Adiabatic
   end interface
!
!  -------------
!  IP Interfaces
!  -------------
!
   interface
      module subroutine IP_ComputeGradient( self , mesh ) 
         use DGWeakIntegrals
         implicit none
         class(IPMethod_t)   ,  intent (in)    :: self
         class(QuadMesh_t)    , intent (inout) :: mesh
      end subroutine IP_ComputeGradient

      module pure function IP_ComputeInnerFluxes( self , e ) result (Fv)
         use QuadElementClass
         implicit none
         class(IPMethod_t),   intent(in)   :: self
         class(QuadElement_t), intent(in)   :: e
         real(kind=RP)                      :: Fv(1:NCONS , 0 : e % spA % N , 0 : e % spA % N , 1:NDIM)
      end function IP_ComputeInnerFluxes

      module pure function IP_RiemannSolver( self , edge , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: invh_edge
         real(kind=RP), intent(in)        :: uL(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: uR(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: dUL(1:NCONS , 0:N, 1:NDIM)
         real(kind=RP), intent(in)        :: dUR(1:NCONS , 0:N, 1:NDIM)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(1:NCONS , 0:N)
      end function IP_RiemannSolver

      module pure function IP_RiemannSolver_Dirichlet( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
         implicit none
         class(IPMethod_t), intent(in)     :: self
         class(Edge_t)    , intent(in)          :: edge
         integer      ,          intent(in)     :: N 
         real(kind=RP),          intent(in)     :: invh_edge
         real(kind=RP),          intent(in)     :: u(NCONS)
         real(kind=RP),          intent(in)     :: g(NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(NCONS)
         real(kind=RP),          intent(in)     :: normal(NDIM)
         real(kind=RP)                          :: Fstar(1:NCONS)
      end function IP_RiemannSolver_Dirichlet

      module pure function IP_RiemannSolver_Adiabatic( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
         implicit none
         class(IPMethod_t), intent (in) :: self
         class(Edge_t)    ,  intent(in)  :: edge
         integer      ,      intent (in) :: N
         real(kind=RP),      intent (in) :: invh_edge
         real(kind=RP),      intent (in) :: u      ( 1:NCONS , 0:N         )
         real(kind=RP),      intent (in) :: g      ( 1:NCONS , 0:N  , NDIM )
         real(kind=RP),      intent (in) :: uB     ( 1:NCONS , 0:N         )
         real(kind=RP),      intent (in) :: normal ( NDIM    , 0:N         )
         real(kind=RP)                   :: Fstar  ( 1:NCONS , 0:N         )
      end function IP_RiemannSolver_Adiabatic

      module pure subroutine IP_GradientRiemannSolver( self , edge , N , UL , UR , normal , GstarL , GstarR ) 
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: uL(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: uR(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP), intent(out)       :: GstarL(1:NCONS , 0:N ,1:NDIM)
         real(kind=RP), intent(out)       :: GstarR(1:NCONS , 0:N ,1:NDIM)
      end subroutine IP_GradientRiemannSolver

      module pure function IP_GradientRiemannSolver_BoundaryCondition( self , edge , u , uB , normal ) result ( Gstar ) 
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         real(kind=RP), intent(in)        :: u(1:NCONS)
         real(kind=RP), intent(in)        :: uB(1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY)
         real(kind=RP)                    :: Gstar(1:NCONS , 1:NDIM)
      end function IP_GradientRiemannSolver_BoundaryCondition

      module pure function IP_GradientRiemannSolver_Adiabatic( self , edge , N , u , uB , normal )  result ( Gstar )
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: u (1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: uB(1:NCONS , 0:N)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Gstar(1:NCONS , 0:N , 1:NDIM)
      end function IP_GradientRiemannSolver_Adiabatic
   end interface
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function ViscousMethod_Initialization() result( ViscousMethod )
         implicit none
         class(ViscousMethod_t), pointer       :: ViscousMethod
!
!        --------------------------------------
!           Prepare the second order method
!        --------------------------------------
!
         if ( trim( Setup % viscous_discretization ) .eq. "IP" ) then

            allocate(IPMethod_t  :: ViscousMethod)

         elseif ( trim( Setup % viscous_discretization ) .eq. "BR1" ) then

            allocate(BR1Method_t :: ViscousMethod)
            ViscousMethod % computeRiemannGradientFluxes = .false.      ! No gradient fluxes are needed

         else

            write(STD_OUT , * ) "Method ",trim(Setup % viscous_discretization)," not implemented yet."
            errorMessage(STD_OUT)
            STOP "Stopped."

         end if

  
         ViscousMethod % method = trim( Setup % viscous_discretization )
            
         select type (ViscousMethod)

            type is (IPMethod_t) 

               ViscousMethod % sigma0 = Setup % sigma0IP
               ViscousMethod % sigma1 = 0.0_RP

               if (trim ( Setup % IPMethod ) .eq. "SIPG" ) then

                  ViscousMethod % subType = "SIPG"
                  ViscousMethod % epsilon = -1.0_RP
                  ViscousMethod % computeRiemannGradientFluxes = .true.

              else if ( trim ( Setup % IPMethod ) .eq. "NIPG") then
      
                  ViscousMethod % subType = "NIPG"
                  ViscousMethod % epsilon = 1.0_RP
                  ViscousMethod % computeRiemannGradientFluxes = .true.

              else if ( trim ( Setup % IPMethod ) .eq. "IIPG") then

                  ViscousMethod % subType = "IIPG"
                  ViscousMethod % epsilon = 0.0_RP
                  ViscousMethod % computeRiemannGradientFluxes = .false.

              else

                  write(STD_OUT , *) "Method ",trim(Setup % IPMethod)," not implemented yet."
                  errorMessage(STD_OUT)
                  stop "Stopped."

              end if
      
            type is (BR1Method_t)

            class default

                write(STD_OUT , *) "Second order method allocation went wrong."
                errorMessage(STD_OUT)
                stop "Stopped."

         end select

         call ViscousMethod % describe()
      end function ViscousMethod_Initialization
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           BASE CLASS METHODS
!           ------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeGradient( self , mesh ) 
         implicit none
         class(ViscousMethod_t),    intent(in)     :: self
         class(QuadMesh_t)     ,    intent(inout)  :: mesh
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end subroutine BaseClass_ComputeGradient

      pure function BaseClass_ComputeInnerFluxes( self , e  ) result (Fv)
         use QuadElementClass
         implicit none
         class(ViscousMethod_t), intent(in)   :: self
         class(QuadElement_t),   intent(in)   :: e
         real(kind=RP)                        :: Fv(1 : NCONS , 0 : e % spA % N , 0 : e % spA % N , 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end function BaseClass_ComputeInnerFluxes

      pure function BaseClass_SolutionRiemannSolver( self , N , uL , uR ) result ( uStar )
         implicit none
         class(ViscousMethod_t), intent(in)        :: self
         integer,                intent(in)        :: N 
         real(kind=RP),          intent(in)        :: uL(1:NCONS , 0:N)
         real(kind=RP),          intent(in)        :: uR(1:NCONS , 0:N)
         real(kind=RP)                             :: uStar(1:NCONS , 0:N)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end function BaseClass_SolutionRiemannSolver

      pure function BaseClass_RiemannSolver( self , edge , N , invh_edge , uL , uR , duL , duR , normal ) result ( FStar )
         implicit none
         class(ViscousMethod_t), intent(in)        :: self
         class(Edge_t)    ,      intent(in)        :: edge
         integer,                intent(in)        :: N 
         real(kind=RP),          intent(in)        :: invh_edge
         real(kind=RP),          intent(in)        :: uL(1:NCONS , 0:N)
         real(kind=RP),          intent(in)        :: uR(1:NCONS , 0:N)
         real(kind=RP),          intent(in)        :: duL(1:NCONS , 0:N , 1:NDIM)
         real(kind=RP),          intent(in)        :: duR(1:NCONS , 0:N , 1:NDIM)
         real(kind=RP),          intent(in)        :: normal(1:NDIM,0:N)
         real(kind=RP)                             :: FStar(1:NCONS , 0:N)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end function BaseClass_RiemannSolver

      pure function BaseClass_RiemannSolver_Dirichlet( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
         implicit none
         class(ViscousMethod_t), intent(in)     :: self
         class(Edge_t)    , intent(in)          :: edge
         integer      ,          intent(in)     :: N 
         real(kind=RP),          intent(in)     :: invh_edge
         real(kind=RP),          intent(in)     :: u(NCONS)
         real(kind=RP),          intent(in)     :: g(NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(NCONS)
         real(kind=RP),          intent(in)     :: normal(NDIM)
         real(kind=RP)                          :: Fstar(1:NCONS)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end function BaseClass_RiemannSolver_Dirichlet

      pure function BaseClass_RiemannSolver_Adiabatic( self , edge , N , invh_edge , u , g , uB , normal ) result ( Fstar )
         implicit none
         class(ViscousMethod_t), intent(in)     :: self
         class(Edge_t)    , intent(in)          :: edge
         integer      ,          intent(in)     :: N 
         real(kind=RP),          intent(in)     :: invh_edge
         real(kind=RP),          intent(in)     :: u(1:NCONS , 0:N)
         real(kind=RP),          intent(in)     :: g(1:NCONS , 0:N , NDIM)
         real(kind=RP),          intent(in)     :: uB(1:NCONS , 0:N)
         real(kind=RP),          intent(in)     :: normal(NDIM , 0:N)
         real(kind=RP)                          :: Fstar(1:NCONS , 0:N)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end function BaseClass_RiemannSolver_Adiabatic

      pure subroutine BaseClass_GradientRiemannSolver( self , edge , N , UL , UR , normal , GstarL , GstarR ) 
!
!        *****************************************************************************************
!        *****************************************************************************************
!
         implicit none
         class(ViscousMethod_t), intent(in) :: self
         class(Edge_t)         , intent(in) :: edge
         integer, intent(in)                :: N
         real(kind=RP), intent(in)          :: uL(1:NCONS , 0:N)
         real(kind=RP), intent(in)          :: uR(1:NCONS , 0:N)
         real(kind=RP), intent(in)          :: normal(IX:IY,0:N)
         real(kind=RP), intent(out)         :: GstarL(1:NCONS , 0:N , 1:NDIM)
         real(kind=RP), intent(out)         :: GstarR(1:NCONS , 0:N , 1:NDIM)
      end subroutine BaseClass_GradientRiemannSolver

      pure function BaseClass_GradientRiemannSolver_BoundaryCondition( self , edge , u , uB , normal ) result ( Gstar )
         implicit none
         class(ViscousMethod_t), intent(in)  :: self
         class(Edge_t), intent(in)           :: edge
         real(kind=RP), intent(in)           :: u(1:NCONS)
         real(kind=RP), intent(in)           :: uB(1:NCONS)
         real(kind=RP), intent(in)           :: normal(1:NDIM)
         real(kind=RP)                       :: Gstar(1:NCONS , 1:NDIM)
      end function BaseClass_GradientRiemannSolver_BoundaryCondition

      pure function BaseClass_GradientRiemannSolver_Adiabatic( self , edge , N , u , uB , normal ) result ( Gstar )
         implicit none
         class(ViscousMethod_t), intent(in)  :: self
         class(Edge_t), intent(in)           :: edge
         integer      , intent(in)           :: N 
         real(kind=RP), intent(in)           :: u(1:NCONS , 0:N)
         real(kind=RP), intent(in)           :: uB(1:NCONS , 0:N)
         real(kind=RP), intent(in)           :: normal(1:NDIM,0:N)
         real(kind=RP)                       :: Gstar(1:NCONS , 0:N , 1:NDIM)
      end function BaseClass_GradientRiemannSolver_Adiabatic
         
      subroutine BaseClass_ComputeSolutionRiemann_Interior( self ,  ed , FuStarL , FuStarR )
         use QuadElementClass
         use MatrixOperations
         implicit none
         class(ViscousMethod_t)              :: self
         type(Edge_t)            :: ed
         real(kind=RP)           :: FuStarL( 1:NCONS , 0 : ed % storage(LEFT ) % spA % N , 1:NDIM )
         real(kind=RP)           :: FuStarR( 1:NCONS , 0 : ed % storage(RIGHT) % spA % N , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), target :: QL ( 1 : NCONS , 0 : ed % spA % N )
         real(kind=RP), target :: QR ( 1 : NCONS , 0 : ed % spA % N )
         real(kind=RP)         :: uStar(1 : NCONS , 0 : ed % spA % N ) 
         real(kind=RP)         :: uStarL(1 : NCONS , 0 : ed % storage(LEFT) % spA % N ) 
         real(kind=RP)         :: uStarR(1 : NCONS , 0 : ed % storage(RIGHT) % spA % N ) 
!
!        Project the solution onto the edge
!        ----------------------------------
         call ed % ProjectSolution( ed , QL , QR ) 
!
!        Compute the solution Riemann solver
!        ----------------------------------- 
         uStar = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!        Return to each element space
!        ----------------------------
         call ed % ProjectFluxes( ed , uStar , uStarL , uStarR ) 
!
!        Compute the solution flux
!        -------------------------
         FuStarL(:,:,IX) = uStarL * ed % dS(0) * ed % n(IX,0)
         FuStarL(:,:,IY) = uStarL * ed % dS(0) * ed % n(IY,0)
         FuStarR(:,:,IX) = uStarR * ed % dS(0) * ed % n(IX,0)
         FuStarR(:,:,IY) = uStarR * ed % dS(0) * ed % n(IY,0)

      end subroutine BaseClass_ComputeSolutionRiemann_Interior

      subroutine BaseClass_ComputeSolutionRiemann_CurvedInterior( self ,  ed , FuStarL , FuStarR )
         use QuadElementClass
         use MatrixOperations
         implicit none
         class(ViscousMethod_t) :: self
         type(CurvedEdge_t)     :: ed
         real(kind=RP)          :: FuStarL( 1 : NCONS , 0 : ed % storage(LEFT ) % spA % N , 1:NDIM)
         real(kind=RP)          :: FuStarR( 1 : NCONS , 0 : ed % storage(RIGHT) % spA % N , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), target :: QL ( 1 : NCONS , 0 : ed % spA % N )
         real(kind=RP), target :: QR ( 1 : NCONS , 0 : ed % spA % N )
         real(kind=RP)         :: uStar( 1 : NCONS , 0 : ed % spA % N ) 
         real(kind=RP)         :: FuStar( 1 : NCONS , 0 : ed % spA % N , 1 : NDIM )
         integer               :: i , eq , dimID
!
!        Project the solution onto the edge
!        ----------------------------------
         call ed % ProjectSolution( ed , QL , QR ) 
!
!        Compute the solution Riemann solver
!        ----------------------------------- 
         uStar = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!        Compute the solution flux
!        -------------------------
         do dimID = 1 , NDIM ;   do i = 0 , ed % spA % N
            FuStar(:,i,dimID) = uStar(:,i) * ed % dS(i) * ed % n(dimID,i)
         end do              ;   end do
!
!        Return its value to each element frame
!        --------------------------------------
         call ed % ProjectFluxes( ed , FuStar(:,:,IX) , FuStarL(:,:,IX) , FuStarR(:,:,IX) )
         call ed % ProjectFluxes( ed , FuStar(:,:,IY) , FuStarL(:,:,IY) , FuStarR(:,:,IY) )

      end subroutine BaseClass_ComputeSolutionRiemann_CurvedInterior

      subroutine BaseClass_ComputeSolutionRiemann_Subdivided( self ,  ed , FuStarL , FuStarRN , FuStarRS )
         use QuadElementClass
         use MatrixOperations
         implicit none
         class(ViscousMethod_t) :: self
         type(SubdividedEdge_t) :: ed
         real(kind=RP)          :: FuStarL (1 : NCONS , 0 : ed % storage(LEFT       ) % spA % N , 1:NDIM)
         real(kind=RP)          :: FuStarRN(1 : NCONS , 0 : ed % storage(RIGHT_NORTH) % spA % N , 1:NDIM)
         real(kind=RP)          :: FuStarRS(1 : NCONS , 0 : ed % storage(RIGHT_SOUTH) % spA % N , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), target :: QLN    ( 1 : NCONS , 0 : ed % spA_N % N ) 
         real(kind=RP), target :: QLS    ( 1 : NCONS , 0 : ed % spA_S % N ) 
         real(kind=RP), target :: QRN    ( 1 : NCONS , 0 : ed % spA_N % N ) 
         real(kind=RP), target :: QRS    ( 1 : NCONS , 0 : ed % spA_S % N ) 
         real(kind=RP)         :: uStarN ( 1 : NCONS , 0 : ed % spA_N % N ) 
         real(kind=RP)         :: uStarS ( 1 : NCONS , 0 : ed % spA_S % N ) 
         real(kind=RP)         :: uStarL (1 : NCONS , 0 : ed % storage(LEFT       ) % spA % N )
         real(kind=RP)         :: uStarRN(1 : NCONS , 0 : ed % storage(RIGHT_NORTH) % spA % N )
         real(kind=RP)         :: uStarRS(1 : NCONS , 0 : ed % storage(RIGHT_SOUTH) % spA % N )
         integer               :: i , l 
!
!        Get the solution projection onto the edge
!        -----------------------------------------
         QLN = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QLN(:,i)    = QLN(:,i) + ed % T_LN_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                     ; end do

         QLS = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QLS(:,i)    = QLS(:,i) + ed % T_LS_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                     ; end do

         QRN = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QRN(:,i)    = QRN(:,i) + ed % T_RN_FWD(i,l) * ed % storage(RIGHT_NORTH) % Q(:,l)
         end do                     ; end do

         QRS = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QRS(:,i)    = QRS(:,i) + ed % T_RS_FWD(i,l) * ed % storage(RIGHT_SOUTH) % Q(:,l)
         end do                     ; end do
!
!        Compute the solution Riemann solver
!        -----------------------------------
         uStarN = self % SolutionRiemannSolver( ed % spA_N % N , QLN , QRN )
         uStarS = self % SolutionRiemannSolver( ed % spA_S % N , QLS , QRS )
!
!        Return the flux to each element space
!        -------------------------------------
         uStarL = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N  
            uStarL(:,i) = uStarL(:,i) + ed % T_LN_BKW(i,l) * uStarN(:,l) + ed % T_LS_BKW(i,l) * uStarS(:,l) 
         end do                     ; end do

         uStarRN = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N  
            uStarRN(:,i) = uStarRN(:,i) - ed % T_RN_BKW(i,l) * uStarN(:,l)  
         end do                     ; end do

         uStarRS = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N  
            uStarRS(:,i) = uStarRS(:,i) - ed % T_RS_BKW(i,l) * uStarS(:,l)  
         end do                     ; end do
!
!        Compute the solution fluxes
!        ---------------------------
         FuStarL(:,:,IX)  = uStarL  * ed % dS(0)   * ed % n(IX,0)
         FuStarL(:,:,IY)  = uStarL  * ed % dS(0)   * ed % n(IY,0)
         FuStarRN(:,:,IX) = -uStarRN * ed % dS_N(0) * ed % normal_N(IX,0)
         FuStarRN(:,:,IY) = -uStarRN * ed % dS_N(0) * ed % normal_N(IY,0)
         FuStarRS(:,:,IX) = -uStarRS * ed % dS_S(0) * ed % normal_S(IX,0)
         FuStarRS(:,:,IY) = -uStarRS * ed % dS_S(0) * ed % normal_S(IY,0)

      end subroutine BaseClass_ComputeSolutionRiemann_Subdivided

      subroutine BaseClass_ComputeSolutionRiemann_CurvedSubdivided( self ,  ed , FuStarL , FuStarRN , FuStarRS )
         use QuadElementClass
         use MatrixOperations
         implicit none
         class(ViscousMethod_t)       :: self
         type(CurvedSubdividedEdge_t) :: ed
         real(kind=RP)                :: FuStarL (1 : NCONS , 0 : ed % storage(LEFT       ) % spA % N , 1:NDIM)
         real(kind=RP)                :: FuStarRN(1 : NCONS , 0 : ed % storage(RIGHT_NORTH) % spA % N , 1:NDIM)
         real(kind=RP)                :: FuStarRS(1 : NCONS , 0 : ed % storage(RIGHT_SOUTH) % spA % N , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), target :: QLN    ( 1 : NCONS , 0 : ed % spA_N % N ) 
         real(kind=RP), target :: QLS    ( 1 : NCONS , 0 : ed % spA_S % N ) 
         real(kind=RP), target :: QRN    ( 1 : NCONS , 0 : ed % spA_N % N ) 
         real(kind=RP), target :: QRS    ( 1 : NCONS , 0 : ed % spA_S % N ) 
         real(kind=RP)         :: uStarN ( 1 : NCONS , 0 : ed % spA_N % N ) 
         real(kind=RP)         :: uStarS ( 1 : NCONS , 0 : ed % spA_S % N ) 
         real(kind=RP)         :: FuStarN ( 1 : NCONS , 0 : ed % spA_N % N , 1:NDIM) 
         real(kind=RP)         :: FuStarS ( 1 : NCONS , 0 : ed % spA_S % N , 1:NDIM) 
         integer               :: dimID , eq , i , l 
!
!        Get the solution projection onto the edge
!        -----------------------------------------
         QLN = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QLN(:,i)    = QLN(:,i) + ed % T_LN_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                     ; end do

         QLS = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QLS(:,i)    = QLS(:,i) + ed % T_LS_FWD(i,l) * ed % storage(LEFT) % Q(:,l)
         end do                     ; end do

         QRN = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QRN(:,i)    = QRN(:,i) + ed % T_RN_FWD(i,l) * ed % storage(RIGHT_NORTH) % Q(:,l)
         end do                     ; end do

         QRS = 0.0_RP
         do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N 
            QRS(:,i)    = QRS(:,i) + ed % T_RS_FWD(i,l) * ed % storage(RIGHT_SOUTH) % Q(:,l)
         end do                     ; end do
!
!        Compute the solution Riemann solver
!        -----------------------------------
         uStarN = self % SolutionRiemannSolver( ed % spA_N % N , QLN , QRN )
         uStarS = self % SolutionRiemannSolver( ed % spA_S % N , QLS , QRS )
!
!        Compute the solution flux
!        -------------------------
         do dimID = 1 , NDIM ; do i = 0 , ed % spA_N % N
            FuStarN(:,i,dimID) = uStarN(:,i) * ed % dS_N(i) * ed % normal_N(dimID,i)
         end do              ; end do 

         do dimID = 1 , NDIM ; do i = 0 , ed % spA_S % N
            FuStarS(:,i,dimID) = uStarS(:,i) * ed % dS_S(i) * ed % normal_S(dimID,i)
         end do              ; end do
!
!        Return the flux to each element space
!        -------------------------------------
         FuStarL = 0.0_RP
         do dimID = 1 , NDIM  ; do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N  
            FuStarL(:,i,dimID) = FuStarL(:,i,dimID) + ed % T_LN_BKW(i,l) * FuStarN(:,l,dimID) + ed % T_LS_BKW(i,l) * FuStarS(:,l,dimID) 
         end do                     ; end do               ; end do

         FuStarRN = 0.0_RP
         do dimID = 1 , NDIM  ; do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N  
            FuStarRN(:,i,dimID) = FuStarRN(:,i,dimID) - ed % T_RN_BKW(i,l) * FuStarN(:,l,dimID)  
         end do               ; end do                     ; end do

         FuStarRS = 0.0_RP
         do dimID = 1 , NDIM  ; do l = 0 , ed % spA % N    ; do i = 0 , ed % spA % N  
            FuStarRS(:,i,dimID) = FuStarRS(:,i,dimID) - ed % T_RS_BKW(i,l) * FuStarS(:,l,dimID)
         end do               ; end do                     ; end do

      end subroutine BaseClass_ComputeSolutionRiemann_CurvedSubdivided

      subroutine BaseClass_ComputeSolutionRiemann_StraightBdry( self ,  ed , FuStar)
!
!        ****************************************************************************
!              This subroutine computes the Solution Riemann problem for a straight
!           edge "ed" -> uStar.
!
!              The solution at the boundaries is prescribed by the value uSB. 
!           However, if periodic boundary conditions are considered, an interior
!           boundary treatment is considered within the edge and the prescribed
!           value.
!        ****************************************************************************
!
         use QuadElementClass
         implicit none
         class(ViscousMethod_t)   :: self
         type(StraightBdryEdge_t) :: ed
         real(kind=RP)            :: FuStar( 1 : NCONS , 0 : ed % spA % N , 1:NDIM)
!     
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: QL(1 : NCONS , 0 : ed % spA % N) , QR(1 : NCONS , 0 : ed % spA % N)
         real(kind=RP) :: uStar( 1 : NCONS , 0 : ed % spA % N )
         integer       :: N , i 

         N = ed % spA % N

         if ( ed % associated ) then     
!
!           ----------------------------
!>          Periodic boundary conditions
!           ----------------------------
!
            if ( .not. ed % inverse ) then
               QL = ed % storage(1) % Q
               QR = ed % uB
               uStar = self % SolutionRiemannSolver( ed % spA % N , QL , QR )

            else
               QL = ed % storage(1) % Q
               do i = 0 , ed % spA % N 
                  QR(:,i) = ed % uB(:,N-i)
               end do
               uStar = self % SolutionRiemannSolver( ed % spA % N , QL , QR ) 

            end if

         else
!
!           -------------------------------------
!>          Dirichlet/Neumann boundary conditions
!           -------------------------------------
!
            uStar = ed % uSB
   
         end if
!
!        Compute the solution flux
!        -------------------------
         FuStar(:,:,IX) = uStar * ed % dS(0) * ed % n(IX,0)
         FuStar(:,:,IY) = uStar * ed % dS(0) * ed % n(IY,0)

      end subroutine BaseClass_ComputeSolutionRiemann_StraightBdry

      subroutine BaseClass_ComputeSolutionRiemann_CurvedBdry( self ,  ed , FuStar)
!
!        ****************************************************************************
!              This subroutine computes the Solution Riemann problem for a straight
!           edge "ed" -> uStar.
!
!              The solution at the boundaries is prescribed by the value uSB. 
!           However, if periodic boundary conditions are considered, an interior
!           boundary treatment is considered within the edge and the prescribed
!           value.
!        ****************************************************************************
!

         use QuadElementClass
         implicit none
         class(ViscousMethod_t) :: self
         type(CurvedBdryEdge_t) :: ed
         real(kind=RP)          :: FuStar( 1:NCONS , 0 : ed % spA % N , 1:NDIM)
!     
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: uStar( 1 : NCONS , 0 : ed % spA % N )
         real(kind=RP) :: QL(1 : NCONS , 0 : ed % spA % N) , QR(1 : NCONS , 0 : ed % spA % N )
         integer       :: N
         integer       :: i , eq , dimID

         N = ed % spA % N

         if ( ed % associated ) then     
!
!           ----------------------------
!>          Periodic boundary conditions
!           ----------------------------
!
            if ( .not. ed % inverse ) then
               QL = ed % storage(1) % Q
               QR = ed % uB
               uStar = self % SolutionRiemannSolver( ed % spA % N , QL , QR )

            else
               QL = ed % storage(1) % Q
               do i = 0 , N
                  QR(:,i) = ed % uB(:,N-i)
               end do
               uStar = self % SolutionRiemannSolver( ed % spA % N , QL , QR ) 

            end if

         else
!
!           -------------------------------------
!>          Dirichlet/Neumann boundary conditions
!           -------------------------------------
!
            uStar = ed % uSB
   
         end if
!
!        Compute the solution flux
!        -------------------------
         do dimID = 1 , NDIM ;   do i = 0 , ed % spA % N  
            FuStar(:,i,dimID) = uStar(:,i) * ed % dS(i) * ed % n(dimID,i)
         end do              ;   end do     

      end subroutine BaseClass_ComputeSolutionRiemann_CurvedBdry
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           AUXILIAR PROCEDURES
!           -------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ViscousMethod_describe( self )
         use Headers
         implicit none
         class(ViscousMethod_t)          :: self

         write(STD_OUT,'(/)') 
         call SubSection_Header("Viscous discretization")
         write(STD_OUT,'(30X,A,A15,A)') "-> ","Method: " , trim( self % method ) 
        
         select type (self)
            type is (IPMethod_t)
               write(STD_OUT,'(30X,A,A15,A)') "-> ","Sub method: " , trim( self % SubType) 
               write(STD_OUT,'(30X,A,A15,F10.4)') "-> ","Sigma0: " , self % sigma0
               write(STD_OUT,'(30X,A,A15,F10.4)') "-> ","Sigma1: " , self % sigma1
               write(STD_OUT,'(30X,A,A15,F10.4)') "-> ","Epsilon: " , self % epsilon
            
            type is (BR1Method_t)
!           ---------------------
!              Nothing to add 
!           ---------------------
            class default

         end select

      end subroutine ViscousMethod_describe

end module DGViscousMethods
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
#endif
!
