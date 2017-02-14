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
         real(kind=RP), allocatable :: Fstar(:,:)
         real(kind=RP), allocatable :: Fstar2D(:,:,:)
         real(kind=RP), pointer     :: lj2D(:,:)
         integer                    :: direction
         integer                    :: pos
         integer                    :: index

         associate ( N => edge % spA % N )
!
!        Allocate the interface flux F·n
!        -------------------------------
         allocate( Fstar( 0 : N , NCONS ) )

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
               associate ( QDot => edge % quads(LEFT) % e % QDot )
                  QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , Fstar )
               end associate
!
!              The obtained term is added to the RIGHT element
!              -----------------------------------------------
               associate ( QDot => edge % quads(RIGHT) % e % QDot ) 
                  QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , Fstar )
               end associate
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
         real(kind=RP), allocatable :: Fstar(:,:)
         real(kind=RP), allocatable :: Fstar2D(:,:,:)
         real(kind=RP), pointer     :: lj2D(:,:)
         integer                    :: direction
         integer                    :: pos
         integer                    :: index
!
         associate ( N => edge % spA % N )

         allocate( Fstar( 0 : N , NCONS ) )

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
               associate ( QDot => edge % quads(LEFT) % e % QDot )
                  QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , Fstar - edge % F(0:N , 1:NCONS , LEFT) )
               end associate
!
!              The obtained term is added to the RIGHT element
!              -----------------------------------------------
               associate ( QDot => edge % quads(RIGHT) % e % QDot ) 
                  QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , Fstar - edge % F(0:N , 1:NCONS , RIGHT ) )
               end associate
!
!           --------------------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )

                  if ( .not. edge % inverted ) then
!
!                    If the normal points towards the domain exterior
!                    ------------------------------------------------
                     QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar - edge % F (0:N , 1:NCONS , 1 ) )

                  else
!
!                    If the normal points towards the domain interior
!                    ------------------------------------------------
                     QDot = QDot + StdDG_QDotFaceContribution( edge , 1 , Fstar - edge % F (0:N , 1:NCONS , 1 ) )

                  end if
               
               end associate
!
!           --------------------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )
!
!                 The normal for curved elements always points towards the domain exterior
!                 ------------------------------------------------------------------------
                  QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar - edge % F (0:N , 1:NCONS , 1 ) ) 
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
         real(kind=RP)              :: Fstar(0:,:)
         real(kind=RP), allocatable :: dFJ(:,:,:)
!        ---------------------------------------------------------------------
         real(kind=RP), allocatable         :: Fstar2D(:,:,:)
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

               allocate(Fstar2D(1,0:N,NCONS))
               
               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D(1 , 0:N    , 1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) )

               elseif ( direction .eq. BACKWARD ) then
!
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(1 , N:0:-1 , 1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS ) )

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if
!
!           ------------------------------------------------------------------------------------------------
            case (ETOP)
!
!              Set the three parameters
!              ------------------------      
               direction = - e % edgesDirection ( edge % edgeLocation(loc) )
               pos       = RIGHT
               index     = IY

               allocate(Fstar2D(0:N,1,NCONS))
   
               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D(0:N,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(N:0:-1,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) ) 

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if
!
!           ------------------------------------------------------------------------------------------------
            case (ELEFT)
!
!              Set the three parameters
!              ------------------------   
               direction = - e % edgesDirection ( edge % edgeLocation(loc) )
               pos       = LEFT
               index     = IX

               allocate(Fstar2D(1,0:N,NCONS))
   
               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D(1,0:N,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(1,N:0:-1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) ) 

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if
!
!           ------------------------------------------------------------------------------------------------
            case (EBOTTOM)

               direction = e % edgesDirection ( edge % edgeLocation(loc) )
               pos       = LEFT
               index     = iY

               allocate(Fstar2D(0:N,1,NCONS))
   
               if ( direction .eq. FORWARD ) then
!        
!                 Introduce the result in the same order
!                 --------------------------------------
                  Fstar2D(0:N,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(N:0:-1,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS) ) 

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."

               end if
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
         allocate ( dFJ(0:N,0:N,NCONS) )
         dFJ =  MatrixMultiplyInIndex_F( Fstar2D , lj2D , index ) 
!
!        Free the variables
!        ------------------
         lj2D=>NULL()
         deallocate(Fstar2D)

         end associate

      end function StdDG_QDotFaceContribution

      subroutine StdDG_QDotVolumeLoop( self , element )
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
         integer                 :: eq
!
!        Compute inner Euler fluxes
!        --------------------------
         call self % computeInnerFluxes ( element )
!
!        Perform the matrix multiplication
!        ---------------------------------
         associate( QDot => element % QDot     , &
                    MD   => element % spA % MD , &
                    M    => element % spA % M  , &
                    w    => element % spA % w  , &
                    N    => element % spA % N      )

         do eq = 1 , NCONS
!
!           FORM I
!           ------
            if ( self % formulation .eq. FORMI ) then
!
!              F Loop
!              ------
               call Mat_x_Mat(A = MD ,B = MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IX) , w , 2 ) , C=QDot(0:N,0:N,eq) , &
                           trA = .true. , reset = .false. )

!
!              G Loop
!              ------
               call Mat_x_Mat(A = MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IY) , w , 1) , B = MD , C=QDot(0:N,0:N,eq) , &
                            reset = .false. )
!
!           FORM II
!           -------
            elseif ( self % formulation .eq. FORMII ) then
!
!              F Loop
!              ------
               call Mat_x_Mat(A = -MD ,B = MatrixByVectorInIndex_F(element % F(0:N,0:N,eq,IX) , w , 2 ) , C=QDot(0:N,0:N,eq) , &
                            reset = .false. )

!
!              G Loop
!              ------
               call Mat_x_Mat(A = -MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IY) , w , 1) , B = MD , C=QDot(0:N,0:N,eq) , &
                            trB = .true. , reset = .false. )

            end if

         end do

         end associate

      end subroutine StdDG_QDotVolumeLoop

      subroutine StdDG_ComputeInnerFluxes( self , element ) 
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
!        -------------------------------------------------------------
         real(kind=RP), allocatable :: F(:,:,:,:)
         integer                    :: eq
         integer                    :: FJa(NDIM) , GJa(NDIM)

         associate( N => element % spA % N )

         allocate( F(0:N , 0:N , NCONS , NDIM) ) 

         F = InviscidFlux( element % Q , element % W)

         do eq = 1 , NCONS
!        
!           F flux (contravariant)
!           ----------------------
            FJa = [1,1]
            GJa = [2,1]
            element % F(0:N,0:N,eq,IX) = F(0:N,0:N,eq,IX) * element % Ja(FJa) + F(0:N,0:N,eq,IY) * element % Ja(GJa)
!        
!           G flux (contravariant)
!           ----------------------
            FJa = [1,2]
            GJa = [2,2]
            element % F(0:N,0:N,eq,IY) = F(0:N,0:N,eq,IX) * element % Ja(FJa) + F(0:N,0:N,eq,IY) * element % Ja(GJa)
         end do

         deallocate( F ) 

         end associate
         
      end subroutine StdDG_ComputeInnerFluxes

      subroutine OIDG_QDotVolumeLoop( self , element )
         use MatrixOperations
         implicit none
         class(OverIntegrationDG_t)          :: self
         class(QuadElement_t)                  :: element
!
!        **********************************************
!           Still under development...
!        **********************************************
!
      end subroutine OIDG_QDotVolumeLoop

      subroutine SplitDG_QDotVolumeLoop( self , element )
         implicit none
         class(SplitDG_t)        :: self
         class(QuadElement_t)      :: element
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
         implicit none
         class(InviscidMethod_t)    :: self
         class(Edge_t), pointer     :: edge
         real(kind=RP), allocatable :: Fstar(:,:)
!        -------------------------------------------------------
         real(kind=RP), allocatable :: QL(:) , QR(:)
         real(kind=RP), pointer     :: T(:,:) , Tinv(:,:)
         integer                    :: iXi
!
!        ********************
         select type ( edge )
!        ********************
!
!           -----------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate( N => edge % spA % N )

               allocate( Fstar ( 0 : N , NCONS ) )

               if ( associated( edge % FB ) ) then
!
!                 Prescribed boundary conditions
!                 ------------------------------
                  Fstar = edge % FB

               elseif ( associated ( edge % uB ) ) then
!
!                 Weak boundary conditions
!                 -----------------------
                  do iXi = 0 , N
                     allocate ( QL(NCONS) , QR(NCONS) )
!
!                    Select LEFT and RIGHT states depending on the edge orientation
!                    -------------------------------------------------------------- 
                     if ( edge % inverted ) then
                        QR = edge % Q(iXi , 1:NCONS , 1)
                        QL = edge % uB(iXi, 1:NCONS)
                     else
                        QL = edge % Q(iXi , 1:NCONS , 1)
                        QR = edge % uB(iXi, 1:NCONS)
                     end if
!
!                    Gather edge orientation matrices
!                    --------------------------------
                     T    => edge % T    ( 1:NCONS , 1:NCONS , iXi )
                     Tinv => edge % Tinv ( 1:NCONS , 1:NCONS , iXi )
!
!                    Compute the Riemann Flux
!                    ------------------------
                     if ( associated ( edge % RiemannSolver ) ) then
!        
!                       Using an edge special Riemann Solver
!                       ------------------------------------
                        Fstar(iXi , :) = edge % RiemannSolver(QL , QR , T , Tinv)
                        
                     else
!
!                       Using the same Riemann Solver than the interior edges
!                       -----------------------------------------------------
                        Fstar(iXi , :) = self % RiemannSolver(QL , QR , T , Tinv)
   
                     end if
     
                     deallocate( QL , QR )
   
                  end do

               end if

               end associate
!               
!           -----------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate( N => edge % spA % N )

               allocate( Fstar ( 0 : N , NCONS ) )

               if ( associated( edge % FB ) ) then
!
!                 Prescribed boundary conditions
!                 ------------------------------
                  Fstar = edge % FB

               elseif ( associated ( edge % uB ) ) then
!
!                 Weak boundary conditions
!                 -----------------------
                  do iXi = 0 , N
                     allocate ( QL(NCONS) , QR(NCONS) )
!
!                    Select LEFT and RIGHT states depending on the edge orientation
!                    -------------------------------------------------------------- 
                     if ( edge % inverted ) then
                        QR = edge % Q(iXi , 1:NCONS , 1)
                        QL = edge % uB(iXi, 1:NCONS)
                     else
                        QL = edge % Q(iXi , 1:NCONS , 1)
                        QR = edge % uB(iXi, 1:NCONS)
                     end if
!
!                    Gather edge orientation matrices
!                    --------------------------------
                     T    => edge % T    ( 1:NCONS , 1:NCONS , iXi )
                     Tinv => edge % Tinv ( 1:NCONS , 1:NCONS , iXi )
!
!                    Compute the Riemann Flux
!                    ------------------------
                     if ( associated ( edge % RiemannSolver ) ) then
!        
!                       Using an edge special Riemann Solver
!                       ------------------------------------
                        Fstar(iXi , :) = edge % RiemannSolver(QL , QR , T , Tinv)
                        
                     else
!
!                       Using the same Riemann Solver than the interior edges
!                       -----------------------------------------------------
                        Fstar(iXi , :) = self % RiemannSolver(QL , QR , T , Tinv)
   
                     end if
     
                     deallocate( QL , QR )
   
                  end do

               end if

               end associate
!
!           -----------------------------------------------------------------
            type is (Edge_t)
      
               associate( N => edge % spA % N )

               allocate( Fstar ( 0 : N , NCONS ) )

               do iXi = 0 , N
!
!                 Select LEFT and RIGHT states
!                 ----------------------------
                  allocate( QL(NCONS) , QR(NCONS) )
                  QL    = edge % Q(iXi , 1:NCONS , LEFT )
                  QR    = edge % Q(iXi , 1:NCONS , RIGHT)
!
!                 Gather edge orientation matrices
!                 -------------------------------- 
                  T     => edge % T(1:NCONS , 1:NCONS , iXi)
                  Tinv  => edge % Tinv(1:NCONS , 1:NCONS , iXi)
!
!                 Compute the Riemann flux
!                 ------------------------
                  Fstar(iXi , : ) = self % RiemannSolver(QL , QR , T , Tinv)

                  deallocate( QL , QR )
  
               end do

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
