   module DGFirstOrderMethods
   use SMConstants
   use QuadElementClass
   use Physics
   use NodesAndWeights_class
   use QuadMeshDefinitions
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
      integer                                   :: formulation
      procedure(RiemannSolverFunction), pointer, nopass :: RiemannSolver => NULL()
      contains
         procedure, non_overridable :: QDotFaceLoop         => StdDG_QDotFaceLoop
         procedure, non_overridable :: RiemannFlux          => FirstOrderMethod_RiemannFlux
         procedure                  :: QDotVolumeLoop       => StdDG_QDotVolumeLoop
         procedure                  :: ComputeInnerFluxes   => StdDG_ComputeInnerFluxes
         procedure, non_overridable :: Describe             => FirstOrderMethod_describe
   end type FirstOrderMethod_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(FirstOrderMethod_t) :: StandardDG_t
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
            write(STD_OUT , '(10X,A)') "Options available are:"
            write(STD_OUT , '(20X,A)') "* Standard"
            write(STD_OUT , '(20X,A)') "* Over-Integration"
            write(STD_OUT , '(20X,A)') "* Split"
            STOP "Stopped."

         end if

         FirstOrderMethod % method = trim( Setup % inviscid_discretization )


         select type (FirstOrderMethod)
            type is (StandardDG_t)
               FirstOrderMethod % formulation = Setup % inviscid_formulation

            type is (OverIntegrationDG_t)
   
            type is (SplitDG_t)

            class default
         end select

!        Set the Riemann flux
         if (trim( Setup % inviscid_flux ) .eq. "Roe") then
            FirstOrderMethod % RiemannSolver => RoeFlux
         else
            write(STD_OUT , *) "Solver ", trim ( Setup % inviscid_flux) ," not implemented yet."
         end if

         call FirstOrderMethod % describe
 

      end function FirstOrderMethod_Initialization

      subroutine StdDG_QDotFaceLoop( self , edge )
         use MatrixOperations
         implicit none
         class(FirstOrderMethod_t)          :: self
         class(Edge_t), pointer             :: edge
         real(kind=RP), allocatable         :: Fstar(:,:)
         real(kind=RP), allocatable         :: Fstar2D(:,:,:)
         real(kind=RP), pointer             :: lj2D(:,:)
         integer                            :: direction
         integer                            :: pos
         integer                            :: index
!
!        -------------------------------------------
!           This is a standard fluxes term
!        -------------------------------------------
!
         associate ( N => edge % spA % N )

         allocate( Fstar( 0 : N , NEC ) )
!        Compute the averaged flux
         Fstar = self % RiemannFlux( edge )
!
!        Perform the loop in both elements
         select type ( edge )
            type is (Edge_t)

            associate ( QDot => edge % quads(LEFT) % e % QDot )
               QDot = QDot - StdDG_QDotFaceContribution( edge , LEFT , Fstar )
            end associate
            associate ( QDot => edge % quads(RIGHT) % e % QDot ) 
               QDot = QDot + StdDG_QDotFaceContribution( edge , RIGHT , Fstar )
            end associate

            type is (StraightBdryEdge_t)
            associate ( QDot => edge % quads(1) % e % QDot )
               QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar )
            end associate

            type is (CurvedBdryEdge_t)
            associate ( QDot => edge % quads(1) % e % QDot )
               QDot = QDot - StdDG_QDotFaceContribution( edge , 1 , Fstar ) 
            end associate

            class default
         end select

        end associate
 
      end subroutine StdDG_QDotFaceLoop

      function StdDG_QDotFaceContribution( edge , loc , Fstar ) result ( dFJ )
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

         associate( N=> edge % quads(loc) % e % spA % N, &
             e=> edge % quads(loc) % e)
 
            if ( edge % edgeLocation(loc) .eq. ERIGHT) then
   
               direction = e % edgesDirection( edge % edgeLocation(loc) )
               pos = RIGHT                    ! Where to interpolate the test function
               index = iX                     ! The coordinate (xi/eta) in which the boundary is located

               allocate(Fstar2D(1,0:N,NEC))
               
               if ( direction .eq. FORWARD ) then
                  Fstar2D(1,0:N,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) )
               elseif ( direction .eq. BACKWARD ) then
                  Fstar2D(1,N:0:-1,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC ) )
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            elseif ( edge % edgeLocation(loc) .eq. ETOP ) then
      
               direction = - e % edgesDirection ( edge % edgeLocation(loc) ) 
               pos = RIGHT
               index = iY
               allocate(Fstar2D(0:N,1,NEC))
   
               if ( direction .eq. FORWARD ) then
                  Fstar2D(0:N,1,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) ) 
               elseif ( direction .eq. BACKWARD ) then
                  Fstar2D(N:0:-1,1,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) ) 
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            elseif ( edge % edgeLocation(loc) .eq. ELEFT ) then
   
               direction = - e % edgesDirection ( edge % edgeLocation(loc) )
               pos = LEFT
               index = iX
               allocate(Fstar2D(1,0:N,NEC))
   
               if ( direction .eq. FORWARD ) then
                  Fstar2D(1,0:N,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) ) 
               elseif ( direction .eq. BACKWARD ) then
                  Fstar2D(1,N:0:-1,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) ) 
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            elseif ( edge % edgeLocation(loc) .eq. EBOTTOM ) then

               direction = e % edgesDirection ( edge % edgeLocation(loc) ) 
               pos = LEFT
               index = iY
               allocate(Fstar2D(0:N,1,NEC))
   
               if ( direction .eq. FORWARD ) then
                  Fstar2D(0:N,1,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) ) 
               elseif ( direction .eq. BACKWARD ) then
                  Fstar2D(N:0:-1,1,1:NEC) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NEC) ) 
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            end if

             lj2D(1:1,0:N) => e % spA % lb(0:N,pos)
             
             allocate ( dFJ(0:N,0:N,NEC) )
             dFJ =  MatrixMultiplyInIndex_F( Fstar2D , lj2D , index ) 

             lj2D=>NULL()
             deallocate(Fstar2D)
         end associate

      end function StdDG_QDotFaceContribution




      subroutine StdDG_QDotVolumeLoop( self , element )
         use MatrixOperations
         implicit none
         class(FirstOrderMethod_t) :: self
         class(QuadElement_t)      :: element
         integer                   :: eq
!
!        -------------------------------------------
!           The standard DG computes the volume
!        terms as:
!              tr(D)*M*F(Q) = tr(MD)*F(Q)
!        -------------------------------------------
!
!        Compute fluxes
         call self % computeInnerFluxes ( element )
!
!        Perform the matrix multiplication
         associate( QDot => element % QDot     , &
                    MD   => element % spA % MD , &
                    M    => element % spA % M  , &
                    N    => element % spA % N      )

         do eq = 1 , NEC

            if ( self % formulation .eq. FORMI ) then
!
!              F Loop
!              ------
               call Mat_x_Mat(A = MD ,B = Mat_x_Mat_F(element % F(0:N,0:N,eq,IX) , M ) , C=QDot(0:N,0:N,eq) , &
                           trA = .true. , reset = .false. )

!
!              G Loop
!              ------
               call Mat_x_Mat(A = Mat_x_Mat_F( M , element % F(0:N,0:N,eq,IY) ) , B = MD , C=QDot(0:N,0:N,eq) , &
                            reset = .false. )

            elseif ( self % formulation .eq. FORMII ) then
!
!              F Loop
!              ------
               call Mat_x_Mat(A = -MD ,B = Mat_x_Mat_F(element % F(0:N,0:N,eq,IX) , M ) , C=QDot(0:N,0:N,eq) , &
                            reset = .false. )

!
!              G Loop
!              ------
               call Mat_x_Mat(A = -Mat_x_Mat_F( M , element % F(0:N,0:N,eq,IY) ) , B = MD , C=QDot(0:N,0:N,eq) , &
                            trB = .true. , reset = .false. )

            end if

         end do
         end associate

      end subroutine StdDG_QDotVolumeLoop

      subroutine StdDG_ComputeInnerFluxes( self , element ) 
         implicit none  
         class(FirstOrderMethod_t)                :: self
         class(QuadElement_t)                     :: element
         real(kind=RP), allocatable               :: F(:,:,:,:)
         integer                                  :: eq
         integer                                  :: FJa(NDIM) , GJa(NDIM)

         associate( N => element % spA % N )

         allocate( F(0:N , 0:N , NEC , NDIM) ) 

         F = InviscidFlux( element % Q )

            do eq = 1 , NEC
               FJa = [1,1]
               GJa = [2,1]
               element % F(0:N,0:N,eq,IX) = F(0:N,0:N,eq,IX) * element % Ja(FJa) + F(0:N,0:N,eq,IY) * element % Ja(GJa)

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
               if ( self % formulation .eq. FORMI ) then
                  write(STD_OUT , '(30X,A,A)') "Formulation: Green form"
               elseif ( self % formulation .eq. FORMII ) then
                  write(STD_OUT , '(30X,A,A)') "Formulation: Divergence form"
               end if
   
            type is ( OverIntegrationDG_t )
         
            type is ( SplitDG_t )
               write(STD_OUT , '(30X,A,F10.4)') "Split op. coefficient: " , self % alpha

            class default
      
         end select

      end subroutine FirstOrderMethod_describe
   
      function FirstOrderMethod_RiemannFlux( self , edge ) result( Fstar )
         implicit none
         class(FirstOrderMethod_t)        :: self
         class(Edge_t), pointer           :: edge
         real(kind=RP), allocatable       :: Fstar(:,:)
         real(kind=RP), allocatable       :: QL(:) , QR(:) 
         real(kind=RP), pointer  :: T(:,:) , Tinv(:,:)
         integer                          :: iXi

         select type ( edge )

            type is (StraightBdryEdge_t)
               associate( N => edge % spA % N )
               allocate( Fstar ( 0 : N , NEC ) )

               if ( associated( edge % FB ) ) then
                  Fstar = edge % FB

               elseif ( associated ( edge % uB ) ) then
                  do iXi = 0 , N
                     allocate ( QL(NEC) , QR(NEC) )
   
                     QL = edge % Q(iXi , 1:NEC , 1)
                     QR = edge % uB(iXi, 1:NEC)
         
                     T  => edge % T(1:NEC , 1:NEC , iXi)
                     Tinv => edge % Tinv(1:NEC , 1:NEC , iXi)
                     Fstar(iXi , :) = self % RiemannSolver(QL , QR , T , Tinv)
         
                     deallocate( QL , QR )
   
                  end do
               end if

               end associate
               
            type is (CurvedBdryEdge_t)
               associate( N => edge % spA % N )

               allocate( Fstar ( 0 : N , NEC ) )

               if ( associated( edge % FB) ) then
                  Fstar = edge % FB

               elseif ( associated ( edge % uB) ) then
                  do iXi = 0 , N
                     allocate ( QL(NEC) , QR(NEC) )
   
                     QL = edge % Q(iXi , 1:NEC , 1)
                     QR = edge % uB(iXi, 1:NEC)
         
                     T  => edge % T(1:NEC , 1:NEC , iXi)
                     Tinv => edge % Tinv(1:NEC , 1:NEC , iXi)
      
                     Fstar(iXi , :) = self % RiemannSolver(QL , QR , T , Tinv)
         
                     deallocate( QL , QR )
   
                  end do
                  
               end if
      
               end associate
 
            type is (Edge_t)
      
               associate( N => edge % spA % N )

               allocate( Fstar ( 0 : N , NEC ) )

               do iXi = 0 , N

                  allocate( QL(NEC) , QR(NEC) )

                  QL    = edge % Q(iXi , 1:NEC , LEFT )
                  QR    = edge % Q(iXi , 1:NEC , RIGHT)
   
                  T     => edge % T(1:NEC , 1:NEC , iXi)
                  Tinv  => edge % Tinv(1:NEC , 1:NEC , iXi)

                  Fstar(iXi , : ) = self % RiemannSolver(QL , QR , T , Tinv)

                  deallocate( QL , QR )
  
               end do
               end associate
            class default
         end select
      end function FirstOrderMethod_RiemannFlux

end module DGFirstOrderMethods
