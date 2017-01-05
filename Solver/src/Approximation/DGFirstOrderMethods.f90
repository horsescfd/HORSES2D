module DGFirstOrderMethods
   use SMConstants
   use QuadElementClass
   use Physics
   use NodesAndWeights_class
   implicit none
!
!  *******************************************************************
   private
   public FirstOrderMethod_t , StandardDG_t , OverIntegrationDG_t , SplitDG_t
   public FirstOrderMethod_Initialization
!  *******************************************************************
!
!                                *************************
   integer, parameter         :: STR_LEN_FIRSTORDER = 128
!                                *************************
!
!  *******************************************************************
   type FirstOrderMethod_t
      character(len=STR_LEN_FIRSTORDER)         :: method
      procedure(RiemannSolverFunction), pointer, nopass :: RiemannSolver => NULL()
      contains
         procedure, non_overridable :: QDotFaceLoop => BaseClass_QDotFaceLoop
         procedure, non_overridable :: RiemannFlux => FirstOrderMethod_RiemannFlux
         procedure                  :: QDotVolumeLoop => BaseClass_QDotVolumeLoop
         procedure, non_overridable :: Describe => FirstOrderMethod_describe
   end type FirstOrderMethod_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(FirstOrderMethod_t) :: StandardDG_t
      contains
         procedure ::  QDotVolumeLoop => StdDG_QDotVolumeLoop
   end type StandardDG_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(FirstOrderMethod_t) ::  OverIntegrationDG_t
      contains
         procedure ::  QDotVolumeLoop => OIDG_QDotVolumeLoop
   end type OverIntegrationDG_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(FirstOrderMethod_t) ::  SplitDG_t
      real(kind=RP)              :: alpha
      contains
         procedure ::  QDotVolumeLoop => SplitDG_QDotVolumeLoop
   end type SplitDG_t
!
!  ========  
   contains
!  ========  
!
      function FirstOrderMethod_Initialization() result( FirstOrderMethod )
         use Setup_class
         implicit none
         class(FirstOrderMethod_t), pointer        :: FirstOrderMethod
!
!        --------------------------------------
!           Prepare the first order method
!        --------------------------------------
!
         if ( trim( Setup % inviscid_discretization ) .eq. "Standard" ) then
            
            allocate(StandardDG_t   :: FirstOrderMethod)

         elseif ( trim( Setup % inviscid_discretization ) .eq. "Over-Integration" ) then 

            allocate(OverIntegrationDG_t  :: FirstOrderMethod)

         elseif ( trim( Setup % inviscid_discretization ) .eq. "Split") then
      
            allocate(SplitDG_t      :: FirstOrderMethod)
      
         else
            write(STD_OUT , *) "Method ", trim(Setup % inviscid_discretization), " not implemented yet."
            STOP "Stopped."

         end if

         FirstOrderMethod % method = trim( Setup % inviscid_discretization )


         select type (FirstOrderMethod)
            type is (StandardDG_t)

            type is (OverIntegrationDG_t)
   
            type is (SplitDG_t)

            class default
         end select

!        Set the Riemann flux
         if (trim( Setup % inviscid_flux ) .eq. "Roe") then
            FirstOrderMethod % RiemannSolver => NULL()

         elseif ( trim ( Setup % inviscid_flux) .eq. "ECON") then
            FirstOrderMethod % RiemannSolver => NULL()

         elseif ( trim ( Setup % inviscid_flux ) .eq. "LLF") then
            FirstOrderMethod % RiemannSolver => NULL()

         else
            write(STD_OUT , *) "Solver ", trim ( Setup % inviscid_flux) ," not implemented yet."
         end if

         call FirstOrderMethod % describe
 

      end function FirstOrderMethod_Initialization

      subroutine BaseClass_QDotFaceLoop( self , edge )
         use MatrixOperations
         implicit none
         class(FirstOrderMethod_t)          :: self
         class(Edge_t), pointer             :: edge
         real(kind=RP), dimension(NEC)      :: Fstar
!!
!!        -------------------------------------------
!!           This is a standard fluxes term
!!        -------------------------------------------
!!
!!        Compute the averaged flux
!         Fstar = self % RiemannFlux( edge )
!!
!!        Perform the loop in both elements
!         select type ( edge )
!            type is (Edge_t)
!               associate(QDot => edge % quads(LEFT) % e % QDot, &
!                   N=> edge % quads(LEFT) % e % Interp % N, &
!                   e=> edge % quads(LEFT) % e)
!!     
!!                  This (-) sign comes from the equation ut = -fx !
!                   QDot = QDot - vectorOuterProduct(e % Interp % lb(:,RIGHT) , Fstar) * edge % nb
!
!               end associate
!
!               associate(QDot => edge % quads(RIGHT) % e % QDot, &
!                   N=> edge % quads(RIGHT) % e % Interp % N, &
!                   e=> edge % quads(RIGHT) % e)
!
!                   QDot = QDot - vectorOuterProduct(e % Interp % lb(:,LEFT) , Fstar) * (-1.0_RP * edge % n)
!
!               end associate
!
!
!            type is (BdryEdge_t)
!               associate(QDot => edge % quads(1) % e % QDot , &
!                  N => edge % quads(1) % e % Interp % N , &
!                  e => edge % quads(1) % e)
!
!                  QDot = QDot - vectorOuterProduct(e % Interp % lb(:,edge % BCLocation) , Fstar)* edge % n
!
!               end associate
!            class default
!         end select
!
 
      end subroutine BaseClass_QDotFaceLoop
!
      subroutine BaseClass_QDotVolumeLoop( self , element )
         implicit none
         class(FirstOrderMethod_t)          :: self
         class(QuadElement_t)                       :: element
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end subroutine BaseClass_QDotVolumeLoop

      subroutine StdDG_QDotVolumeLoop( self , element )
         use MatrixOperations
         implicit none
         class(StandardDG_t)     :: self
         class(QuadElement_t)      :: element
!
!        -------------------------------------------
!           The standard DG computes the volume
!        terms as:
!              tr(D)*M*F(Q) = tr(MD)*F(Q)
!        -------------------------------------------
!
!        Compute fluxes
         element % F = inviscidFlux( element % Q )
!
!        Perform the matrix multiplication
!         associate( QDot => element % QDot , &
!                    MD => element % Interp % MD)

!         QDot = QDot + TransposeMat_x_NormalMat_F( MD , element % F )

!         end associate

      end subroutine StdDG_QDotVolumeLoop

      subroutine OIDG_QDotVolumeLoop( self , element )
         use MatrixOperations
         implicit none
         class(OverIntegrationDG_t)          :: self
         class(QuadElement_t)                  :: element
!!
!!        ---------------------------------------------------
!!           The Over-Integration DG computes the volume
!!        terms as:
!!           tr(tildeM T D) * F(tildeQ)
!!        ---------------------------------------------------
!!
!!        Compute fluxes
!!        TODO: beware of this
!         element % F = 0.0_RP !inviscidFlux( matmul(element % Interp % T , element % Q) )
!!
!!        Perform the matrix multiplication
!         associate( QDot => element % QDot , &
!                  tildeMTD => element % Interp % tildeMTD)
!
!         QDot = QDot + TransposeMat_x_NormalMat_F( tildeMTD , element % F )
!
!         end associate
!
      end subroutine OIDG_QDotVolumeLoop

      subroutine SplitDG_QDotVolumeLoop( self , element )
         implicit none
         class(SplitDG_t)        :: self
         class(QuadElement_t)      :: element
!
      end subroutine SplitDG_QDotVolumeLoop


      subroutine FirstOrderMethod_describe( self )
         use Headers
         implicit none
         class(FirstOrderMethod_t)        :: self

         write(STD_OUT,'(/)') 
         call SubSection_Header("Inviscid discretization")
         write(STD_OUT,'(30X,A,A)') "Method: " , trim( self % method ) 

         select type ( self ) 
            type is ( StandardDG_t )
   
            type is ( OverIntegrationDG_t )
         
            type is ( SplitDG_t )
               write(STD_OUT , '(30X,A,F10.4)') "Split op. coefficient: " , self % alpha

            class default
      
         end select

      end subroutine FirstOrderMethod_describe
   
      function FirstOrderMethod_RiemannFlux( self , edge ) result( val )
         implicit none
         class(FirstOrderMethod_t)        :: self
         class(Edge_t), pointer           :: edge
         real(kind=RP), dimension(NEC)    :: val
         real(kind=RP), pointer  :: QL(:) , QR(:) , QBdry(:)
         real(kind=RP), pointer  :: n(:)
!
!         select type ( edge )
!            type is (StraightBdryEdge_t)
!      
!               QBdry => edge % quads(1) % e % Qb( : , edge % BCLocation  )
!!              TODO: Beware of this
!               n     => NULL()   ! face % n
!
!               val = self % RiemannSolver(QBdry , edge % uB , n)
!               
!            type is (Edge_t)
!      
!               QL => edge % quads(LEFT) % e % Qb( : , RIGHT )
!               QR => edge % quads(RIGHT) % e % Qb( : , LEFT )
!!              TODO: Beware of this
!               n  => NULL()      ! face % n
!
!               val = self % RiemannSolver(QL , QR , n)
!
!            class default
!         end select
      end function FirstOrderMethod_RiemannFlux

end module DGFirstOrderMethods
