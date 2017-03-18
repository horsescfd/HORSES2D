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
         procedure, non_overridable :: QDotFaceLoopFormI    => StdDG_QDotFaceLoopFormI
         procedure, non_overridable :: QDotFaceLoopFormII   => StdDG_QDotFaceLoopFormII
         procedure, non_overridable :: RiemannFlux          => InviscidMethod_RiemannFlux
         procedure                  :: QDotVolumeLoop       => StdDG_QDotVolumeLoop
         procedure                  :: ComputeInnerFluxes   => StdDG_ComputeInnerFluxes
         procedure, non_overridable :: Describe             => InviscidMethod_describe
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
            case ( "AUSM" )
               InviscidMethod % RiemannSOlver => AUSMFlux
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

      subroutine StdDG_QDotFaceLoopFormI( self , edge )
!
!        **************************************************************
!              This routine computes the edge loop according to the
!           "Form I" formulation:
!                 QDot +-= -\int_e F^* l_j l_i ds
!           in where F^* represents the F·n product, that is, the 
!           result is substracted to the LEFT element, and added to 
!           the RIGHT element. 
!        **************************************************************
!
         use MatrixOperations
         implicit none
         class(InviscidMethod_t)    :: self
         class(Edge_t), pointer     :: edge
!        -------------------------------------------------------
         real(kind=RP)              :: Fstar(0:edge % spA % N,1:NCONS)
         real(kind=RP)              :: FstarLowered(0:edge % Nlow, 1:NCONS)

         associate ( N => edge % spA % N )
!
!        Compute the edge Riemann Flux is proceeds, or uses the prescribed boundary flux
!        -------------------------------------------------------------------------------
         Fstar = self % RiemannFlux( edge )
!
!        Perform the loop in both elements
!        ---------------------------------
!
!        ********************
         select type ( edge )
!        ********************
!
!           --------------------------------------------------------------------------
            type is (Edge_t)
!
!              The obtained term is substracted to the LEFT element
!              ----------------------------------------------------
               if ( .not. edge % transform(LEFT) ) then
                  associate ( QDot => edge % quads(LEFT) % e % QDot )
                  QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , Fstar )
                  end associate
            
               else
                  associate ( QDot => edge % quads(LEFT) % e % QDot )
                  call Mat_x_Mat( A = edge % T_backward , B = Fstar , C = FstarLowered )
                  QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , FstarLowered )
                  end associate

               end if
!
!              The obtained term is added to the RIGHT element
!              -----------------------------------------------
               if ( .not. edge % transform(RIGHT) ) then
                  associate ( QDot => edge % quads(RIGHT) % e % QDot )
                  QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , Fstar )
                  end associate
            
               else
                  associate ( QDot => edge % quads(RIGHT) % e % QDot )
                  call Mat_x_Mat( A = edge % T_backward , B = Fstar , C = FstarLowered )
                  QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , FstarLowered )
                  end associate

               end if                 
!
!           --------------------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )

                  if ( .not. edge % inverted ) then
!
!                    If the normal points towards the domain exterior
!                    ------------------------------------------------
                     QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar )
                  else
!
!                    If the normal points towards the domain interior
!                    ------------------------------------------------
                     QDot = QDot + StdDG_QDotFaceContribution( edge , 1 , Fstar )
                  end if
               
               end associate
!
!           --------------------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )
!
!                 The normal for curved elements always points towards the domain exterior
!                 ------------------------------------------------------------------------
                  QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar ) 
               end associate
!
!           --------------------------------------------------------------------------
            class default
               STOP "Stopped."
!
!        **********
         end select
!        **********
!
        end associate
 
      end subroutine StdDG_QDotFaceLoopFormI

      subroutine StdDG_QDotFaceLoopFormII( self , edge )
!
!        **************************************************************
!              This routine computes the edge loop according to the
!           "Form II" formulation:
!                 QDot +-= -\int_e (F^*-Fe) l_j l_i ds
!           in where F^* and Fe represents the F·n product, that is, 
!           the result is substracted to the LEFT element, and added to 
!           the RIGHT element. 
!        **************************************************************
!
         use MatrixOperations
         implicit none
         class(InviscidMethod_t)    :: self
         class(Edge_t), pointer     :: edge
!        -------------------------------------------------
         real(kind=RP)              :: Fstar(0:edge % spA % N,1:NCONS)
         real(kind=RP)              :: FstarLowered(0:edge % Nlow, 1:NCONS)
!
         associate ( N => edge % spA % N )
!
!        Compute the edge Riemann Flux is proceeds, or uses the prescribed boundary flux
!        -------------------------------------------------------------------------------
         Fstar = self % RiemannFlux( edge )
!
!        Perform the loop in both elements
!        ---------------------------------
!
!        ********************
         select type ( edge )
!        ********************
!
!           --------------------------------------------------------------------------
            type is (Edge_t)
!
!              The obtained term is substracted to the LEFT element
!              ----------------------------------------------------
               if ( .not. edge % transform(LEFT) ) then
                  associate ( QDot => edge % quads(LEFT) % e % QDot )
                  QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , Fstar - edge % storage(LEFT) % F(0:N , 1:NCONS))
                  end associate
            
               else
                  associate ( QDot => edge % quads(LEFT) % e % QDot )
                  call Mat_x_Mat( A = edge % T_backward , B = Fstar , C = FstarLowered )
                  QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , FstarLowered - edge % storage(LEFT) % F(0:edge % Nlow , 1:NCONS) )
                  end associate

               end if
!
!              The obtained term is added to the RIGHT element
!              -----------------------------------------------
               if ( .not. edge % transform(RIGHT) ) then
                  associate ( QDot => edge % quads(RIGHT) % e % QDot )
                  QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , Fstar - edge % storage(RIGHT) % F(0:N , 1:NCONS))
                  end associate
            
               else
                  associate ( QDot => edge % quads(RIGHT) % e % QDot )
                  call Mat_x_Mat( A = edge % T_backward , B = Fstar , C = FstarLowered )
                  QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , FstarLowered - edge % storage(RIGHT) % F(0:edge % Nlow , 1:NCONS) )
                  end associate

               end if

!
!           --------------------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )

                  if ( .not. edge % inverted ) then
!
!                    If the normal points towards the domain exterior
!                    ------------------------------------------------
                     QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar - edge % storage(1) % F (0:N , 1:NCONS ) )

                  else
!
!                    If the normal points towards the domain interior
!                    ------------------------------------------------
                     QDot = QDot + StdDG_QDotFaceContribution( edge , 1 , Fstar - edge % storage(1) % F (0:N , 1:NCONS ) )

                  end if
               
               end associate
!
!           --------------------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )
!
!                 The normal for curved elements always points towards the domain exterior
!                 ------------------------------------------------------------------------
                  QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar - edge % storage(1) % F (0:N , 1:NCONS ) ) 
               end associate
!
!           --------------------------------------------------------------------------
            class default
               STOP "Stopped."
!
!        **********
         end select
!        **********
!
         end associate
 
 
      end subroutine StdDG_QDotFaceLoopFormII

      function StdDG_QDotFaceContribution( edge , loc , Fstar ) result ( dFJ )
!
!        *************************************************************************************
!           This subroutine computes the following integral
!
!                 ---------------------------------
!                 | dFJ = \int_e Fstar l_j l_i ds |
!                 ---------------------------------
!
!           For a given edge, and for its neighbouring element "loc" (LEFT/RIGHT)
!           The following quantities are used:
!              * direction: the edge direction in the element in its local system
!
!                    \\ Direction is FORWARD if the edge is oriented such that:
!
!                                  --------->
!                                --------------
!                             ^  |            | ^
!                             |  |            | |
!                             |  |            | |
!                             |  |            | |
!                                --------------
!                                  --------->
!
!              * pos: where the lagrange polynomials are interpolated.
!                    \\ pos is RIGHT for RIGHT and TOP edges
!                    \\ pos is LEFT  for LEFT and BOTTOM edges 
!
!              * index: where the boundary is located referred to the element local system 
!                    \\ IX for  xi = const boundaries
!                    \\ IY for eta = const boundaries
!
!        *************************************************************************************
!
         use MatrixOperations
         implicit none
         class(Edge_t)              :: edge
         integer                    :: loc
         real(kind=RP)              :: Fstar(0:edge % storage(loc) % spA % N,1:NCONS)
         real(kind=RP)              :: dFJ(0:edge % storage(loc) % spA % N,0:edge % storage(loc) % spA % N,1:NCONS)
!        ---------------------------------------------------------------------
         real(kind=RP), pointer             :: Fstar2D(:,:,:)
         real(kind=RP), target              :: Fstar2D_x(1,0:edge % storage(loc) % spA % N,1:NCONS)
         real(kind=RP), target              :: Fstar2D_y(0:edge % storage(loc) % spA % N,1,1:NCONS)
         real(kind=RP), pointer             :: lj2D(:,:)
         integer                            :: direction
         integer                            :: pos
         integer                            :: index
!

         associate( N => edge % quads(loc) % e % spA % N, &
                    e => edge % quads(loc) % e)
!
!        **************************************
         select case (edge % edgeLocation(loc))
!        **************************************
!
!           ------------------------------------------------------------------------------------------------
            case (ERIGHT)
!
!              Set the three parameters
!              ------------------------ 
               direction = e % edgesDirection( edge % edgeLocation(loc) )
               pos       = RIGHT                  
               index     = IX                     

               
               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D_x(1 , 0:N    , 1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) , rowC = N+1 , colC = NCONS )

               elseif ( direction .eq. BACKWARD ) then
!
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D_x(1 , N:0:-1 , 1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS ) , rowC = N+1 , colC = NCONS)

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if

               Fstar2D(1:,0:,1:) => Fstar2D_x
!
!           ------------------------------------------------------------------------------------------------
            case (ETOP)
!
!              Set the three parameters
!              ------------------------      
               direction = - e % edgesDirection ( edge % edgeLocation(loc) )
               pos       = RIGHT
               index     = IY

               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D_y(0:N,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) , rowC = N+1 , colC = NCONS) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D_y(N:0:-1,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) , rowC = N+1 , colC = NCONS) 

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if

               Fstar2D(0:,1:,1:) => Fstar2D_y
!
!           ------------------------------------------------------------------------------------------------
            case (ELEFT)
!
!              Set the three parameters
!              ------------------------   
               direction = - e % edgesDirection ( edge % edgeLocation(loc) )
               pos       = LEFT
               index     = IX

               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D_x(1,0:N,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D_x(1,N:0:-1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if

               Fstar2D(1:,0:,1:) => Fstar2D_x
!
!           ------------------------------------------------------------------------------------------------
            case (EBOTTOM)

               direction = e % edgesDirection ( edge % edgeLocation(loc) )
               pos       = LEFT
               index     = iY

               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D_y(0:N,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D_y(N:0:-1,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if

               Fstar2D(0:,1:,1:) => Fstar2D_y
!
!        **********
         end select
!        **********
!
!        Get the interpolated lagrange polynomials
!        ----------------------------------------- 
         lj2D(1:1,0:N) => e % spA % lb(0:N,pos)
!
!        Obtain the result             
!        -----------------
         dFJ =  MatrixMultiplyInIndex_F( Fstar2D , lj2D , N+1 , N+1 , NCONS , index ) 
!
!        Free the variables
!        ------------------
         lj2D=>NULL()
         Fstar2D => NULL()

         end associate

      end function StdDG_QDotFaceContribution

      subroutine StdDG_QDotVolumeLoop( self , element , reset , compute)
!
!        *******************************************************************************
!           This routine computes the standard DG volumetric terms according to:
!                 QDot += tr(D) M F M + M G M D          ( for Form I  )
!                 QDot -= M D F M + M G tr(D) M          ( for Form II ) 
!           The details about this matricial form is shown in the doc HiODG2DTech.pdf
!        *******************************************************************************
!
         use MatrixOperations
         implicit none
         class(InviscidMethod_t) :: self
         class(QuadElement_t)    :: element
         logical, optional       :: reset
         logical, optional       :: compute
         integer                 :: eq
         logical                 :: rst
         logical                 :: cmpt

         if ( present(reset) ) then
            rst = reset

         else
            rst = .true.

         end if

         if ( present(compute) ) then
            cmpt = compute

         else
            cmpt = .true.

         end if
!
!        Compute inner Euler fluxes
!        --------------------------
         call self % computeInnerFluxes ( element , rst)
!
!        Perform the matrix multiplication
!        ---------------------------------
         if ( cmpt ) then
            associate( QDot => element % QDot     , &
                       MD   => element % spA % MD , &
                       trMD => element % spA % trMD, &
                       M    => element % spA % M  , &
                       w    => element % spA % w  , &
                       N    => element % spA % N      )
   
            do eq = 1 , NCONS
!   
!              FORM I
!              ------
               if ( self % formulation .eq. FORMI ) then
!   
!                 F Loop
!                 ------
                  call Mat_x_Mat(A = trMD ,B = MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IX) , w , N+1 , N+1 , 2 ) , C=QDot(0:N,0:N,eq) , &
                                reset = .false. )
   
!   
!                 G Loop
!                 ------
                  call Mat_x_Mat(A = MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IY) , w , N+1 , N+1 , 1) , B = MD , C=QDot(0:N,0:N,eq) , &
                               reset = .false. )
!   
!              FORM II
!              -------
               elseif ( self % formulation .eq. FORMII ) then
!   
!                 F Loop
!                 ------
                  call Mat_x_Mat(A = -MD ,B = MatrixByVectorInIndex_F(element % F(0:N,0:N,eq,IX) , w , N+1 , N+1 ,  2 ) , C=QDot(0:N,0:N,eq) , &
                               reset = .false. )
   
!   
!                 G Loop
!                 ------
                  call Mat_x_Mat(A = -MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IY) , w , N+1 , N+1 , 1) , B = trMD , C=QDot(0:N,0:N,eq) , &
                               reset = .false. )
   
               end if
   
            end do
   
            end associate
         end if

      end subroutine StdDG_QDotVolumeLoop

      subroutine StdDG_ComputeInnerFluxes( self , element , reset) 
!
!        **********************************************************************
!              This subroutine computes the contravariant fluxes of the element
!           The fluxes read:
!                 F <- F * Ja(1,1) + G * Ja(2,1)
!                 G <- F * Ja(1,2) + G * Ja(2,2)
!        **********************************************************************
!
         implicit none  
         class(InviscidMethod_t)    :: self
         class(QuadElement_t)       :: element
         logical                    :: reset
!        -------------------------------------------------------------
         real(kind=RP)              :: F(0:element % spA % N,0:element % spA % N,1:NCONS,1:NDIM)
         integer                    :: eq

         associate( N => element % spA % N )
         
         F = InviscidFlux( element % spA % N , element % Q )

         if ( reset ) then
            do eq = 1 , NCONS
!           
!              F flux (contravariant)
!              ----------------------
               element % F(0:N,0:N,eq,IX) = F(0:N,0:N,eq,IX) * element % Ja(0:N,0:N,1,1) + F(0:N,0:N,eq,IY) * element % Ja(0:N,0:N,2,1)
!           
!              G flux (contravariant)
!              ----------------------
               element % F(0:N,0:N,eq,IY) = F(0:N,0:N,eq,IX) * element % Ja(0:N,0:N,1,2) + F(0:N,0:N,eq,IY) * element % Ja(0:N,0:N,2,2)
            end do
   
         else
            do eq = 1 , NCONS
!           
!              F flux (contravariant)
!              ----------------------
               element % F(0:N,0:N,eq,IX) = element % F(0:N,0:N,eq,IX) + F(0:N,0:N,eq,IX) * element % Ja(0:N,0:N,1,1) + F(0:N,0:N,eq,IY) * element % Ja(0:N,0:N,2,1)
!           
!              G flux (contravariant)
!              ----------------------
               element % F(0:N,0:N,eq,IY) = element % F(0:N,0:N,eq,IX) + F(0:N,0:N,eq,IX) * element % Ja(0:N,0:N,1,2) + F(0:N,0:N,eq,IY) * element % Ja(0:N,0:N,2,2)
            end do
         end if

         end associate
         
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

   
      function InviscidMethod_RiemannFlux( self , edge ) result( Fstar )
!
!        *************************************************************************
!              This function computes the Riemann flux accross one edge.
!           The result is the normal Riemann Flux F^*·n = Fstar = F·dSx + F·dSy
!
!        *************************************************************************
!
         use MatrixOperations
         use Physics
         implicit none
         class(InviscidMethod_t)    :: self
         class(Edge_t), pointer     :: edge
!        -------------------------------------------------------
         real(kind=RP)              :: Fstar(0:edge % spA % N,1:NCONS)
         real(kind=RP)              :: QL(1:NCONS) , QR(1:NCONS)
         real(kind=RP)              :: WL(1:NPRIM) , WR(1:NPRIM)
         real(kind=RP), pointer     :: T(:,:) , Tinv(:,:)
         real(kind=RP)              :: lj_forward(0 : edge % Nlow)
         integer                    :: iXi
         procedure(RiemannSolverFunction), pointer    :: RiemannSolver
!
!        ********************
         select type ( edge )
!        ********************
!
!           -----------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate( N => edge % spA % N )

               select case ( edge % BCWeakType ) 

                  case ( WEAK_PRESCRIBED )
!
!                    Prescribed boundary conditions: just pick its value
!                    ---------------------------------------------------
                     Fstar = edge % FB

                  case ( WEAK_RIEMANN )
!   
!                    Weak boundary conditions
!                    -----------------------
                     if ( associated ( edge % RiemannSolver ) ) then
                        RiemannSolver => edge % RiemannSolver

                     else
                        RiemannSolver => self % RiemannSolver

                     end if

                     if ( edge % inverted ) then

                        do iXi = 0 , N

!     
!                          Select LEFT and RIGHT states depending on the edge orientation
!                          -------------------------------------------------------------- 
                           QR = edge % storage(1) % Q(iXi , 1:NCONS)
                           QL = edge % uB(iXi, 1:NCONS)
                           WR = edge % storage(1) % W(iXi , 1:NPRIM)
                           WL = ComputePrimitiveVariables(QL)


!                          Gather edge orientation matrices
!                          --------------------------------
                           T    => edge % T    ( 1:NCONS , 1:NCONS , iXi )
                           Tinv => edge % Tinv ( 1:NCONS , 1:NCONS , iXi )
!   
!                          Compute the Riemann Flux
!                          ------------------------
                           Fstar(iXi , :) = RiemannSolver(QL , QR , WL , WR , T , Tinv)
 
                        end do

                     else

                        do iXi = 0 , N

!     
!                          Select LEFT and RIGHT states depending on the edge orientation
!                          -------------------------------------------------------------- 
                           QL = edge % storage(1) % Q(iXi , 1:NCONS)
                           QR = edge % uB(iXi, 1:NCONS)
                           WL = edge % storage(1) % W(iXi , 1:NPRIM)
                           WR = ComputePrimitiveVariables(QR)

!                          Gather edge orientation matrices
!                          --------------------------------
                           T    => edge % T    ( 1:NCONS , 1:NCONS , iXi )
                           Tinv => edge % Tinv ( 1:NCONS , 1:NCONS , iXi )
!   
!                          Compute the Riemann Flux
!                          ------------------------
                           Fstar(iXi , :) = RiemannSolver(QL , QR , WL , WR , T , Tinv)
 
                        end do

                     end if
   
                  case default
      
                     print*, "Boundary condition has undefined Type"
                     stop "Stopped."

               end select

               end associate
!               
!           -----------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate( N => edge % spA % N )

               select case ( edge % BCWeakType ) 

                  case ( WEAK_PRESCRIBED )
!
!                    Prescribed boundary conditions: just pick its value
!                    ---------------------------------------------------
                     Fstar = edge % FB

                  case ( WEAK_RIEMANN )
!
!                    Weak boundary conditions
!                    -----------------------
                     if ( associated ( edge % RiemannSolver ) ) then
                        RiemannSolver => edge % RiemannSolver

                     else
                        RiemannSolver => self % RiemannSolver

                     end if

                     if ( edge % inverted ) then

                        do iXi = 0 , N
!     
!                          Select LEFT and RIGHT states depending on the edge orientation
!                          -------------------------------------------------------------- 
                           QR = edge % storage(1) % Q(iXi , 1:NCONS)
                           QL = edge % uB(iXi, 1:NCONS)
                           WR = edge % storage(1) % W(iXi , 1:NPRIM)
                           WL = ComputePrimitiveVariables(QL)
!
!                          Gather edge orientation matrices
!                          --------------------------------
                           T    => edge % T    ( 1:NCONS , 1:NCONS , iXi )
                           Tinv => edge % Tinv ( 1:NCONS , 1:NCONS , iXi )
!   
!                          Compute the Riemann Flux
!                          ------------------------
                           Fstar(iXi , :) = RiemannSolver(QL , QR , WL , WR , T , Tinv)
 
                        end do

                     else

                        do iXi = 0 , N
!     
!                          Select LEFT and RIGHT states depending on the edge orientation
!                          -------------------------------------------------------------- 
                           QL = edge % storage(1) % Q(iXi , 1:NCONS)
                           QR = edge % uB(iXi, 1:NCONS)
                           WL = edge % storage(1) % W(iXi , 1:NPRIM)
                           WR = ComputePrimitiveVariables(QR)
!
!                          Gather edge orientation matrices
!                          --------------------------------
                           T    => edge % T    ( 1:NCONS , 1:NCONS , iXi )
                           Tinv => edge % Tinv ( 1:NCONS , 1:NCONS , iXi )
!   
!                          Compute the Riemann Flux
!                          ------------------------
                           Fstar(iXi , :) = RiemannSolver(QL , QR , WL , WR , T , Tinv)
 
                        end do

                     end if
 
                  case default
      
                     print*, "Boundary condition has undefined Type"
                     stop "Stopped."

               end select

               end associate
!
!           -----------------------------------------------------------------
            type is (Edge_t)
      
               associate( N => edge % spA % N )

               if ( edge % transform(LEFT) ) then

                  do iXi = 0 , N
!
!                    Transform left boundary values
!                    ------------------------------
                     associate ( Nlow => edge % storage(LEFT) % spA % N )
                     lj_forward = edge % T_forward(iXi,0 : NLow)
                     call MatrixTimesVector( A = edge % storage(LEFT) % Q(0 : Nlow , 1:NCONS)  , X = lj_forward , Y = QL , trA = .true. , reset = .true. )
                     call MatrixTimesVector( A = edge % storage(LEFT) % W(0 : Nlow , 1:NPRIM)  , X = lj_forward , Y = WL , trA = .true. , reset = .true. )
                     end associate
!
!                    Get right boundary values
!                    -------------------------
                     QR    = edge % storage(RIGHT) % Q(iXi , 1:NCONS)
                     WR    = edge % storage(RIGHT) % W(iXi , 1:NPRIM)

!
!                    Gather edge orientation matrices
!                    -------------------------------- 
                     T     => edge % T    ( 1 : NCONS , 1 : NCONS , iXi ) 
                     Tinv  => edge % Tinv ( 1 : NCONS , 1 : NCONS , iXi ) 
!
!                    Compute the Riemann flux
!                    ------------------------
                     Fstar(iXi , : ) = self % RiemannSolver(QL , QR , WL , WR , T , Tinv)

                  end do

               elseif ( edge % transform(RIGHT) ) then

                  do iXi = 0 , N
!
!                    Get left boundary values
!                    ------------------------
                     QL    = edge % storage(LEFT) % Q(iXi , 1:NCONS)
                     WL    = edge % storage(LEFT) % W(iXi , 1:NPRIM)
!
!                    Transform right boundary values
!                    -------------------------------
                     associate ( Nlow => edge % storage(RIGHT) % spA % N )
                     lj_forward = edge % T_forward(iXi,0 : NLow)
                     call MatrixTimesVector( A = edge % storage(RIGHT) % Q(0 : Nlow , 1:NCONS)  , X = lj_forward , Y = QR , trA = .true. , reset = .true. )
                     call MatrixTimesVector( A = edge % storage(RIGHT) % W(0 : Nlow , 1:NPRIM)  , X = lj_forward , Y = WR , trA = .true. , reset = .true. )
                     end associate
!
!                    Gather edge orientation matrices
!                    -------------------------------- 
                     T     => edge % T    ( 1 : NCONS , 1 : NCONS , iXi ) 
                     Tinv  => edge % Tinv ( 1 : NCONS , 1 : NCONS , iXi ) 
!
!                    Compute the Riemann flux
!                    ------------------------
                     Fstar(iXi , : ) = self % RiemannSolver(QL , QR , WL , WR , T , Tinv)

                  end do
            
               else

                  do iXi = 0 , N
  
                     QL    = edge % storage(LEFT) % Q(iXi , 1:NCONS)
                     WL    = edge % storage(LEFT) % W(iXi , 1:NPRIM)
                     QR    = edge % storage(RIGHT) % Q(iXi , 1:NCONS)
                     WR    = edge % storage(RIGHT) % W(iXi , 1:NPRIM)
!
!                   Gather edge orientation matrices
!                    -------------------------------- 
                     T     => edge % T    ( 1 : NCONS , 1 : NCONS , iXi ) 
                     Tinv  => edge % Tinv ( 1 : NCONS , 1 : NCONS , iXi ) 
!  
!                    Compute the Riemann flux
!                    ------------------------
                     Fstar(iXi , : ) = self % RiemannSolver(QL , QR , WL , WR , T , Tinv)

                  end do
         
               end if

               end associate
!
!           ------------------------------------------------------------------
            class default
!
!        **********
         end select
!        **********
!
      end function InviscidMethod_RiemannFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           AUXILIAR SUBROUTINES
!           --------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function ComputePrimitiveVariables(Q) result (W)
         implicit none
         real(kind=RP)           :: Q(NCONS)
         real(Kind=RP)           :: W(NPRIM)
         real(kind=RP)           :: invRho

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
