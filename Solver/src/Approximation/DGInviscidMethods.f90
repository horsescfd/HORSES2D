module DGInviscidMethods
   use SMConstants
   use QuadElementClass
   use Physics
   use NodesAndWeights_class
   implicit none

#include "Defines.h"
!
!  *******************************************************************
   private
   public InviscidMethod_t , StandardDG_t , OverIntegrationDG_t , SplitDG_t
   public InviscidMethod_Initialization
!  *******************************************************************
!
!                                *************************
   integer, parameter         :: STR_LEN_INVISCID = 128
!                                *************************
!
!  *******************************************************************
   type InviscidMethod_t
      character(len=STR_LEN_INVISCID) :: method
      integer                         :: formulation
      procedure(RiemannSolverFunction), pointer, nopass :: RiemannSolver => NULL()
      contains
         procedure                  :: ComputeInnerFluxes                => StdDG_ComputeInnerFluxes
         procedure, non_overridable :: Describe                          => InviscidMethod_describe
   end type InviscidMethod_t
!  *******************************************************
!
!  *******************************************************
   type, extends(InviscidMethod_t) :: StandardDG_t
   end type StandardDG_t
!  *******************************************************
!
!  *******************************************************
   type, extends(InviscidMethod_t) ::  OverIntegrationDG_t
      integer           :: no_of_integrationPoints
   end type OverIntegrationDG_t
!  *******************************************************
!
!  *******************************************************
   type, extends(InviscidMethod_t) ::  SplitDG_t
      real(kind=RP)         :: alpha
   end type SplitDG_t
!
!
!  ========  
   contains
!  ========  
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      function InviscidMethod_Initialization() result( InviscidMethod )
         use Setup_class
         implicit none
         class(InviscidMethod_t), pointer        :: InviscidMethod
!
!        --------------------------------------
!           Prepare the first order method
!        --------------------------------------
!
!        ********************************************************
         select case ( trim ( Setup % inviscid_discretization ) ) 
!        ********************************************************
!
!           --------------------------------------------------------  
            case ( "Standard" )
               allocate(StandardDG_t   :: InviscidMethod)
!
!           --------------------------------------------------------  
            case ( "Over-Integration" )
               allocate(OverIntegrationDG_t  :: InviscidMethod)
!
!           --------------------------------------------------------  
            case ( "Split" )      
               allocate(SplitDG_t      :: InviscidMethod)
!      
!           --------------------------------------------------------  
            case default
               if ( len_trim( Setup % inviscid_discretization ) .eq. 0 ) then
                  write(STD_OUT , * ) "Inviscid discretization method not specified."
               else
                  write(STD_OUT , *) "Method ", trim(Setup % inviscid_discretization), " not implemented yet."
               end if

               write(STD_OUT , '(10X,A)') "Options available are:"
               write(STD_OUT , '(20X,A)') "* Standard"
               write(STD_OUT , '(20X,A)') "* Over-Integration"
               write(STD_OUT , '(20X,A)') "* Split"
               STOP "Stopped."
!
!        **********
         end select
!        **********
!
!        Set method
!        ----------
         InviscidMethod % method = trim( Setup % inviscid_discretization )

!        Set method parameters
!        ---------------------
         select type (InviscidMethod)

            type is (StandardDG_t)
               InviscidMethod % formulation = Setup % inviscid_formulation

            type is (OverIntegrationDG_t)
   
            type is (SplitDG_t)

            class default
         end select

!        Set the Riemann flux
!        --------------------
!
!        **********************************************
         select case ( trim ( Setup % inviscid_flux ) ) 
!        **********************************************
!
!           --------------------------------------------------------           
            case ( "Roe" )
               InviscidMethod % RiemannSolver => RoeFlux
!         
!           --------------------------------------------------------           
            case ( "HLL" )
               InviscidMethod % RiemannSolver => HLLFlux
!         
!           --------------------------------------------------------           
            case ( "Exact" ) 
               InviscidMethod % RiemannSOlver => ExactRiemannSolver
!
!           --------------------------------------------------------           
            case default
               if ( len_trim ( Setup % inviscid_flux ) .eq. 0 ) then
                  write(STD_OUT , *) "Riemann solver not specified."
               
               else
                  write(STD_OUT , *) "Riemann solver ", trim ( Setup % inviscid_flux) ," not implemented yet."
               
               end if

               write(STD_OUT , '(10X,A)') "Options available are:"
               write(STD_OUT , '(20X,A)') "* Roe"
               write(STD_OUT , '(20X,A)') "* HLL"
               write(STD_OUT , '(20X,A)') "* AUSM"
               write(STD_OUT , '(20X,A)') "* Exact"
               STOP "Stopped."
!
!        **********
         end select
!        **********
!
!        Describe the method
!        -------------------
         call InviscidMethod % describe
 

      end function InviscidMethod_Initialization

      pure function StdDG_ComputeInnerFluxes( self , e ) result (F) 
!
!        **********************************************************************
!              This subroutine computes the cartesian fluxes of the element
!        **********************************************************************
!
         use Physics
         implicit none  
         class(InviscidMethod_t), intent (in)    :: self
         class(QuadElement_t),    intent (in)    :: e
         real(kind=RP)                           :: F(0:e % spA % N , 0:e % spA % N , 1:NCONS , 1:NDIM)

         F = InviscidFlux( e % spA % N , e % Q )

      end function StdDG_ComputeInnerFluxes
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           AUXILIAR SUBROUTINES
!           --------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine InviscidMethod_describe( self )
         use Headers
         use Setup_class
         implicit none
         class(InviscidMethod_t)        :: self

         write(STD_OUT,'(/)') 
         call SubSection_Header("Inviscid discretization")
         write(STD_OUT,'(30X,A,A18,A)') "-> ","Method: " , trim( self % method ) 
!
!        ********************
         select type ( self ) 
!        ********************
!
!           -------------------------------------------------------------------------------------
            type is ( StandardDG_t )
               if ( self % formulation .eq. FORMI ) then
                  write(STD_OUT , '(30X,A,A18,A)') "-> ","Formulation: ","Green form"
               elseif ( self % formulation .eq. FORMII ) then
                  write(STD_OUT , '(30X,A,A18,A)') "-> ","Formulation: ","Divergence form"
               end if
!           
!           -------------------------------------------------------------------------------------
            type is ( OverIntegrationDG_t )
         
!           -------------------------------------------------------------------------------------
            type is ( SplitDG_t )
               write(STD_OUT , '(30X,A,F10.4)') "Split op. coefficient: " , self % alpha

!           -------------------------------------------------------------------------------------
            class default
!
!        **********
         end select
!        **********
!
         write(STD_OUT,'(30X,A,A18,A)') "-> " , "Riemann solver: " , trim ( Setup % inviscid_flux )
      end subroutine InviscidMethod_describe

end module DGInviscidMethods
