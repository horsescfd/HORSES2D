module DGInviscidMethods
   use SMConstants
   use QuadElementClass
   use Physics
   use NodesAndWeights_class
   use QuadMeshDefinitions
   implicit none
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
         generic, public            :: ComputeRiemannFluxes              => ComputeRiemannFluxes_Interior , ComputeRiemannFluxes_StraightBdry , ComputeRiemannFluxes_CurvedBdry
         procedure, private         :: ComputeRiemannFluxes_Interior     => StdDG_ComputeRiemannFluxes_Interior
         procedure, private         :: ComputeRiemannFluxes_StraightBdry => StdDG_ComputeRiemannFluxes_StraightBdry
         procedure, private         :: ComputeRiemannFluxes_CurvedBdry   => StdDG_ComputeRiemannFluxes_CurvedBdry
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
      contains
         procedure ::  QDotVolumeLoop => OIDG_QDotVolumeLoop
   end type OverIntegrationDG_t
!  *******************************************************
!
!  *******************************************************
   type, extends(InviscidMethod_t) ::  SplitDG_t
      real(kind=RP)         :: alpha
      contains
         procedure ::  QDotVolumeLoop => SplitDG_QDotVolumeLoop
   end type SplitDG_t

   interface
      module subroutine StdDG_ComputeRiemannFluxes_Interior( self , ed , FL , FR )
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(InviscidMethod_t)    :: self
         type(Edge_t)               :: ed
         real ( kind=RP )           :: FL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
      end subroutine StdDG_ComputeRiemannFluxes_Interior

      module subroutine StdDG_ComputeRiemannFluxes_StraightBdry( self , ed , F )
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(InviscidMethod_t)  :: self
         type(StraightBdryEdge_t) :: ed
         real ( kind=RP )         :: F ( 0 : ed % spA % N , 1 : NCONS )
      end subroutine StdDG_ComputeRiemannFluxes_StraightBdry

      module subroutine StdDG_ComputeRiemannFluxes_CurvedBdry( self , ed , F )
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(InviscidMethod_t) :: self
         type(CurvedBdryEdge_t)  :: ed
         real ( kind=RP )        :: F ( 0 : ed % spA % N , 1 : NCONS )
      end subroutine StdDG_ComputeRiemannFluxes_CurvedBdry
   end interface

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

      subroutine StdDG_ComputeInnerFluxes( self , e , F) 
!
!        **********************************************************************
!              This subroutine computes the contravariant fluxes of the element
!           The fluxes read:
!                 F <- F * Ja(1,1) + G * Ja(2,1)
!                 G <- F * Ja(1,2) + G * Ja(2,2)
!        **********************************************************************
!
         use Physics
         implicit none  
         class(InviscidMethod_t), intent (in)    :: self
         class(QuadElement_t),    intent (in)    :: e
         real(kind=RP),           intent (inout) :: F(0:e % spA % N , 0:e % spA % N , 1:NCONS , 1:NDIM)
!        -------------------------------------------------------------
         real(kind=RP)              :: F_cartesian(0:e % spA % N,0:e % spA % N,1:NCONS,1:NDIM)
         integer                    :: eq
         integer, pointer           :: N 

         N => e % spA % N         

         F_cartesian = InviscidFlux( e % spA % N , e % Q )

         do eq = 1 , NCONS
!           
!           F flux (contravariant)
!           ----------------------
            F(0:N,0:N,eq,IX) = F_cartesian(0:N,0:N,eq,IX) * e % Ja(0:N,0:N,1,1) + F_cartesian(0:N,0:N,eq,IY) * e % Ja(0:N,0:N,2,1)
!           
!           G flux (contravariant)
!           ----------------------
            F(0:N,0:N,eq,IY) = F_cartesian(0:N,0:N,eq,IX) * e % Ja(0:N,0:N,1,2) + F_cartesian(0:N,0:N,eq,IY) * e % Ja(0:N,0:N,2,2)
         end do

      end subroutine StdDG_ComputeInnerFluxes

      subroutine OIDG_QDotVolumeLoop( self , element , reset , compute)
         use MatrixOperations
         implicit none
         class(OverIntegrationDG_t)          :: self
         class(QuadElement_t)                  :: element
         logical, optional                   :: reset
         logical, optional                   :: compute
!
!        **********************************************
!           Still under development...
!        **********************************************
!
      end subroutine OIDG_QDotVolumeLoop

      subroutine SplitDG_QDotVolumeLoop( self , element , reset , compute)
         implicit none
         class(SplitDG_t)        :: self
         class(QuadElement_t)      :: element
         logical, optional                   :: reset
         logical, optional                   :: compute
!
!        **********************************************
!           Still under development...
!        **********************************************
!
      end subroutine SplitDG_QDotVolumeLoop
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           AUXILIAR SUBROUTINES
!           --------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure function ComputePrimitiveVariables(Q) result (W)
         implicit none
         real(kind=RP), intent(in) :: Q(NCONS)
         real(Kind=RP)             :: W(NPRIM)
         real(kind=RP)             :: invRho

         W(IRHO) = Q(IRHO)
         invRho  = 1.0_RP / Q(IRHO)
         W(IU)   = Q(IRHOU) * invRho
         W(IV)   = Q(IRHOV) * invRho
         W(IP)   = Thermodynamics % gm1 * ( Q(IRHOE) - 0.5_RP *( Q(IRHOU) * W(IU) + Q(IRHOV) * W(IV))) 
         W(IT)   = W(IP) * invRho
         W(IA)   = sqrt(Thermodynamics % gamma * W(IT))

      end function ComputePrimitiveVariables

      subroutine InviscidMethod_describe( self )
         use Headers
         implicit none
         class(InviscidMethod_t)        :: self

         write(STD_OUT,'(/)') 
         call SubSection_Header("Inviscid discretization")
         write(STD_OUT,'(30X,A,A15,A)') "-> ","Method: " , trim( self % method ) 
!
!        ********************
         select type ( self ) 
!        ********************
!
!           -------------------------------------------------------------------------------------
            type is ( StandardDG_t )
               if ( self % formulation .eq. FORMI ) then
                  write(STD_OUT , '(30X,A,A15,A)') "-> ","Formulation: ","Green form"
               elseif ( self % formulation .eq. FORMII ) then
                  write(STD_OUT , '(30X,A,A15,A)') "-> ","Formulation: ","Divergence form"
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
      end subroutine InviscidMethod_describe

end module DGInviscidMethods
