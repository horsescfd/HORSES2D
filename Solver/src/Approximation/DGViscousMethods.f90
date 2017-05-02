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
         procedure, private ::   ComputeSolutionRiemann_Interior         => BaseClass_ComputeSolutionRiemann_Interior
         procedure, private ::   ComputeSolutionRiemann_CurvedInterior   => BaseClass_ComputeSolutionRiemann_CurvedInterior
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
         class(QuadMesh_t)     ,    intent(inout)  :: mesh
      end subroutine BR1_ComputeGradient

      module pure function BR1_ComputeInnerFluxes( self , e ) result (Fv)
         use QuadElementClass
         implicit none
         class(BR1Method_t),   intent(in)   :: self
         class(QuadElement_t), intent(in)   :: e
         real(kind=RP)                      :: Fv(0 : e % spA % N , 0 : e % spA % N , 1:NCONS , 1:NDIM)
      end function BR1_ComputeInnerFluxes

      module pure function BR1_SolutionRiemannSolver( self , N , UL , UR ) result ( uStar )
         implicit none
         class(BR1Method_t), intent(in) :: self
         integer, intent(in)            :: N
         real(kind=RP), intent(in)      :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)      :: uR(0:N , 1:NCONS)
         real(kind=RP)                  :: uStar(0:N,1:NCONS)
      end function BR1_SolutionRiemannSolver

      module pure function BR1_RiemannSolver( self , edge , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
         implicit none
         class(BR1Method_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: invh_edge
         real(kind=RP), intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: dUL(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: dUR(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(0:N , 1:NCONS)
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
         real(kind=RP),          intent(in)     :: u(0:N , NCONS)
         real(kind=RP),          intent(in)     :: g(0:N , NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(0:N , NCONS)
         real(kind=RP),          intent(in)     :: normal(NDIM , 0:N)
         real(kind=RP)                          :: Fstar(0:N , 1:NCONS)
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
         real(kind=RP)                      :: Fv(0 : e % spA % N , 0 : e % spA % N , 1:NCONS , 1:NDIM)
      end function IP_ComputeInnerFluxes

      module pure function IP_RiemannSolver( self , edge , N , invh_edge , UL , UR , dUL , dUR , normal ) result ( FStar )
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: invh_edge
         real(kind=RP), intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: dUL(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: dUR(0:N, 1:NDIM , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Fstar(0:N , 1:NCONS)
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
         real(kind=RP),      intent (in) :: u      ( 0:N  , NCONS      )
         real(kind=RP),      intent (in) :: g      ( 0:N  , NDIM,NCONS )
         real(kind=RP),      intent (in) :: uB     ( 0:N  , NCONS      )
         real(kind=RP),      intent (in) :: normal ( NDIM , 0:N        )
         real(kind=RP)                   :: Fstar  ( 0:N  , 1:NCONS    )
      end function IP_RiemannSolver_Adiabatic

      module pure subroutine IP_GradientRiemannSolver( self , edge , N , UL , UR , normal , GstarL , GstarR ) 
         implicit none
         class(IPMethod_t), intent(in)   :: self
         class(Edge_t)    , intent(in)    :: edge
         integer, intent(in)              :: N
         real(kind=RP), intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP), intent(out)       :: GstarL(0:N , 1:NCONS , 1:NDIM)
         real(kind=RP), intent(out)       :: GstarR(0:N , 1:NCONS , 1:NDIM)
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
         real(kind=RP), intent(in)        :: u (0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: uB(0:N , 1:NCONS)
         real(kind=RP), intent(in)        :: normal(IX:IY,0:N)
         real(kind=RP)                    :: Gstar(0:N , 1:NCONS , 1:NDIM)
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
         real(kind=RP)                        :: Fv(0 : e % spA % N , 0 : e % spA % N , 1:NCONS , 1:NDIM)
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
         real(kind=RP),          intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP),          intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP)                             :: uStar(0:N , 1:NCONS)
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
         real(kind=RP),          intent(in)        :: uL(0:N , 1:NCONS)
         real(kind=RP),          intent(in)        :: uR(0:N , 1:NCONS)
         real(kind=RP),          intent(in)        :: duL(0:N , 1:NDIM , 1:NCONS)
         real(kind=RP),          intent(in)        :: duR(0:N , 1:NDIM , 1:NCONS)
         real(kind=RP),          intent(in)        :: normal(1:NDIM,0:N)
         real(kind=RP)                             :: FStar(0:N , 1:NCONS)
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
         real(kind=RP),          intent(in)     :: u(0:N , NCONS)
         real(kind=RP),          intent(in)     :: g(0:N , NDIM,NCONS)
         real(kind=RP),          intent(in)     :: uB(0:N , NCONS)
         real(kind=RP),          intent(in)     :: normal(NDIM , 0:N)
         real(kind=RP)                          :: Fstar(0:N , 1:NCONS)
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
         real(kind=RP), intent(in)          :: uL(0:N , 1:NCONS)
         real(kind=RP), intent(in)          :: uR(0:N , 1:NCONS)
         real(kind=RP), intent(in)          :: normal(IX:IY,0:N)
         real(kind=RP), intent(out)         :: GstarL(0:N , 1:NCONS , 1:NDIM)
         real(kind=RP), intent(out)         :: GstarR(0:N , 1:NCONS , 1:NDIM)
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
         real(kind=RP), intent(in)           :: u(0:N,1:NCONS)
         real(kind=RP), intent(in)           :: uB(0:N,1:NCONS)
         real(kind=RP), intent(in)           :: normal(1:NDIM,0:N)
         real(kind=RP)                       :: Gstar(0:N , 1:NCONS , 1:NDIM)
      end function BaseClass_GradientRiemannSolver_Adiabatic
         
      subroutine BaseClass_ComputeSolutionRiemann_Interior( self ,  ed , uStarL , uStarR )
!
!        *****************************************************************************
!           This routine computes the Solution Riemann problem in the "ed" edge.
!              -> The result is uStarL and uStarR
!           This procedure is common to all ViscousMethods, since it just provides
!           a framework to deal with the data sharing across the interface. Then,
!           each method considers its own SolutionRiemann solver procedure
!                                     -------------------
!        >> Therefore, it considers several possibilities:
!
!              1/ LEFT edge needs a p-Transformation and RIGHT edge is reversed.
!              2/ LEFT edge needs a p-Transformation
!              3/ RIGHT edge needs a p-Transformation and RIGHT edge is reversed.
!              4/ RIGHT edge needs a p-Transformation
!              5/ RIGHT edge is reversed.
!              6/ No transformations are needed.
!
!           The solution Riemann solver is computed with the 
!>                self % SolutionRiemannSolver( N , QL , QR )
!           procedure.
!
!        *****************************************************************************
!
         use QuadElementClass
         use MatrixOperations
         implicit none
         class(ViscousMethod_t)              :: self
         type(Edge_t)            :: ed
         real(kind=RP)           :: uStarL( 0 : ed % storage(LEFT ) % spA % N , 1:NCONS )
         real(kind=RP)           :: uStarR( 0 : ed % storage(RIGHT) % spA % N , 1:NCONS )
!
!        ---------------
!        Local variables
!        ---------------
!
         real ( kind=RP ), target :: QL ( 0 : ed % spA % N , 1 : NCONS )
         real ( kind=RP ), target :: QR ( 0 : ed % spA % N , 1 : NCONS )

         if ( ed % transform(LEFT) .and. ed % inverse ) then
! 
!        ---------------------------------------------------------
!>       First case: LEFT needs p-Adaption, and RIGHT is reversed.
!        ---------------------------------------------------------
!
!           Transform the LEFT edge
!           -----------------------            
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
!
!           Invert the RIGHT edge
!           ---------------------
            QR = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS )
!
!           Compute the Riemann solver
!           --------------------------
            uStarR = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the LEFT edge 
!           -----------------------
            call Mat_x_Mat( ed % T_backward , uStarR , uStarL )
!
!           Invert the RIGHT edge
!           ---------------------
            uStarR(0:ed % spA % N , 1:NCONS) = uStarR(ed % spA % N : 0 : -1 , 1:NCONS) 

         elseif ( ed % transform(LEFT) ) then
! 
!        ----------------------------------
!>       Second case: LEFT needs p-Adaption.
!        ----------------------------------
!
!           Transform the LEFT 
!           -----------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
!
!           Get the RIGHT edge state
!           ------------------------
            QR = ed % storage(RIGHT) % Q
!
!           Compute the Riemann solver
!           --------------------------
            uStarR = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the LEFT edge 
!           -----------------------
            call Mat_x_Mat( ed % T_backward , uStarR , uStarL )


         elseif ( ed % transform(RIGHT) .and.  ed % inverse ) then
! 
!        ---------------------------------------------------------
!>       Third case: RIGHT needs p-Adaption, and RIGHT is reversed.
!        ---------------------------------------------------------
!
!           Get the LEFT edge 
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Transform the RIGHT edge 
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
!
!           Invert the RIGHT edge 
!           ---------------------
            QR(0:ed % spA % N,1:NCONS) = QR(ed % spA % N : 0 : -1 , 1:NCONS)
!
!           Compute the Riemann solver
!           --------------------------
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the RIGHT edge
!           ------------------------ 
            call Mat_x_Mat( ed % T_backward , uStarL , uStarR )
!
!           Invert the RIGHT edge
!           ---------------------
            uStarR(0:ed % Nlow,1:NCONS) = uStarR(ed % Nlow:0:-1 , 1:NCONS)

         elseif ( ed % transform(RIGHT) ) then
! 
!        -----------------------------------
!>       Fourth case: RIGHT needs p-Adaption.
!        -----------------------------------
!
!           Get the LEFT edge
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Transform the RIGHT edge
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
!
!           Compute the Riemann solver
!           --------------------------
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the RIGHT edge
!           ------------------------
            call Mat_x_Mat( ed % T_backward , uStarL , uStarR )
         
         elseif ( ed % inverse ) then
! 
!        -----------------------------
!>       Fifth case: RIGHT is reversed.
!        -----------------------------
!
!           Get the LEFT edge
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Invert the RIGHT edge
!           ---------------------
            QR(0:ed % spA % N , 1:NCONS ) = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS)
!
!           Compute the Riemann solver
!           --------------------------         
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Invert the RIGHT edge
!           ---------------------
            uStarR( 0 : ed % spA % N , 1:NCONS ) = uStarL ( ed % spA % N : 0 : -1 , 1:NCONS )

         else
! 
!        -----------------------------------------------------------------------
!>       Sixth case: Default case: neither p-Adaption nor inversion are required.
!        -----------------------------------------------------------------------
!
!           Get the LEFT edge 
!           -----------------   
            QL = ed % storage(LEFT) % Q
!
!           Get the RIGHT edge
!           ------------------
            QR = ed % storage(RIGHT) % Q
!
!           Compute the Riemann solver
!           --------------------------
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Assign the value to the RIGHT edge
!           ----------------------------------
            uStarR = uStarL
         
         end if

      end subroutine BaseClass_ComputeSolutionRiemann_Interior

      subroutine BaseClass_ComputeSolutionRiemann_CurvedInterior( self ,  ed , uStarL , uStarR )
!
!        *****************************************************************************
!           This routine computes the Solution Riemann problem in the "ed" edge.
!              -> The result is uStarL and uStarR
!           This procedure is common to all ViscousMethods, since it just provides
!           a framework to deal with the data sharing across the interface. Then,
!           each method considers its own SolutionRiemann solver procedure
!                                     -------------------
!        >> Therefore, it considers several possibilities:
!
!              1/ LEFT edge needs a p-Transformation and RIGHT edge is reversed.
!              2/ LEFT edge needs a p-Transformation
!              3/ RIGHT edge needs a p-Transformation and RIGHT edge is reversed.
!              4/ RIGHT edge needs a p-Transformation
!              5/ RIGHT edge is reversed.
!              6/ No transformations are needed.
!
!           The solution Riemann solver is computed with the 
!>                self % SolutionRiemannSolver( N , QL , QR )
!           procedure.
!
!        *****************************************************************************
!
         use QuadElementClass
         use MatrixOperations
         implicit none
         class(ViscousMethod_t) :: self
         type(CurvedEdge_t)     :: ed
         real(kind=RP)          :: uStarL( 0 : ed % storage(LEFT ) % spA % N , 1:NCONS )
         real(kind=RP)          :: uStarR( 0 : ed % storage(RIGHT) % spA % N , 1:NCONS )
!
!        ---------------
!        Local variables
!        ---------------
!
         real ( kind=RP ), target :: QL ( 0 : ed % spA % N , 1 : NCONS )
         real ( kind=RP ), target :: QR ( 0 : ed % spA % N , 1 : NCONS )

         if ( ed % transform(LEFT) .and. ed % inverse ) then
! 
!        ---------------------------------------------------------
!>       First case: LEFT needs p-Adaption, and RIGHT is reversed.
!        ---------------------------------------------------------
!
!           Transform the LEFT edge
!           -----------------------            
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
!
!           Invert the RIGHT edge
!           ---------------------
            QR = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS )
!
!           Compute the Riemann solver
!           --------------------------
            uStarR = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the LEFT edge 
!           -----------------------
            call Mat_x_Mat( ed % T_backward , uStarR , uStarL )
!
!           Invert the RIGHT edge
!           ---------------------
            uStarR(0:ed % spA % N , 1:NCONS) = uStarR(ed % spA % N : 0 : -1 , 1:NCONS) 

         elseif ( ed % transform(LEFT) ) then
! 
!        ----------------------------------
!>       Second case: LEFT needs p-Adaption.
!        ----------------------------------
!
!           Transform the LEFT 
!           -----------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
!
!           Get the RIGHT edge state
!           ------------------------
            QR = ed % storage(RIGHT) % Q
!
!           Compute the Riemann solver
!           --------------------------
            uStarR = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the LEFT edge 
!           -----------------------
            call Mat_x_Mat( ed % T_backward , uStarR , uStarL )


         elseif ( ed % transform(RIGHT) .and.  ed % inverse ) then
! 
!        ---------------------------------------------------------
!>       Third case: RIGHT needs p-Adaption, and RIGHT is reversed.
!        ---------------------------------------------------------
!
!           Get the LEFT edge 
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Transform the RIGHT edge 
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
!
!           Invert the RIGHT edge 
!           ---------------------
            QR(0:ed % spA % N,1:NCONS) = QR(ed % spA % N : 0 : -1 , 1:NCONS)
!
!           Compute the Riemann solver
!           --------------------------
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the RIGHT edge
!           ------------------------ 
            call Mat_x_Mat( ed % T_backward , uStarL , uStarR )
!
!           Invert the RIGHT edge
!           ---------------------
            uStarR(0:ed % Nlow,1:NCONS) = uStarR(ed % Nlow:0:-1 , 1:NCONS)

         elseif ( ed % transform(RIGHT) ) then
! 
!        -----------------------------------
!>       Fourth case: RIGHT needs p-Adaption.
!        -----------------------------------
!
!           Get the LEFT edge
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Transform the RIGHT edge
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
!
!           Compute the Riemann solver
!           --------------------------
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Transform the RIGHT edge
!           ------------------------
            call Mat_x_Mat( ed % T_backward , uStarL , uStarR )
         
         elseif ( ed % inverse ) then
! 
!        -----------------------------
!>       Fifth case: RIGHT is reversed.
!        -----------------------------
!
!           Get the LEFT edge
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Invert the RIGHT edge
!           ---------------------
            QR(0:ed % spA % N , 1:NCONS ) = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS)
!
!           Compute the Riemann solver
!           --------------------------         
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Invert the RIGHT edge
!           ---------------------
            uStarR( 0 : ed % spA % N , 1:NCONS ) = uStarL ( ed % spA % N : 0 : -1 , 1:NCONS )

         else
! 
!        -----------------------------------------------------------------------
!>       Sixth case: Default case: neither p-Adaption nor inversion are required.
!        -----------------------------------------------------------------------
!
!           Get the LEFT edge 
!           -----------------   
            QL = ed % storage(LEFT) % Q
!
!           Get the RIGHT edge
!           ------------------
            QR = ed % storage(RIGHT) % Q
!
!           Compute the Riemann solver
!           --------------------------
            uStarL = self % SolutionRiemannSolver( ed % spA % N , QL , QR )
!
!           Assign the value to the RIGHT edge
!           ----------------------------------
            uStarR = uStarL
         
         end if

      end subroutine BaseClass_ComputeSolutionRiemann_CurvedInterior

      subroutine BaseClass_ComputeSolutionRiemann_StraightBdry( self ,  ed , uStar)
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
         real(kind=RP)            :: uStar( 0 : ed % spA % N , 1:NCONS )
!     
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                 :: QL(0 : ed % spA % N , 1:NCONS) , QR(0 : ed % spA % N , 1:NCONS)
         integer                       :: N

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
               QR = ed % uB(N : 0 : -1 , 1:NCONS)
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

      end subroutine BaseClass_ComputeSolutionRiemann_StraightBdry

      subroutine BaseClass_ComputeSolutionRiemann_CurvedBdry( self ,  ed , uStar)
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
         real(kind=RP)          :: uStar( 0 : ed % spA % N , 1:NCONS )
!     
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                 :: QL(0 : ed % spA % N , 1:NCONS) , QR(0 : ed % spA % N , 1:NCONS)
         integer                       :: N

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
               QR = ed % uB(N : 0 : -1 , 1:NCONS)
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
