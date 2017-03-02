#ifdef NAVIER_STOKES

module DGViscousMethods
   use SMConstants
   use QuadMeshClass
   use QuadElementClass
   use Physics
   use Setup_class
   use QuadMeshDefinitions
   implicit none
!
!  *******************************************************************
   private
   public ViscousMethod_t , IPMethod_t , BR1Method_t , LDGMethod_t
   public ViscousMethod_Initialization
!  *******************************************************************
!
!                                *************************
   integer, parameter         :: STR_LEN_VISCOUS = 128
!                                *************************
!
!  *******************************************************************
   type ViscousMethod_t
      character(len=STR_LEN_VISCOUS)     :: method
      contains
         procedure ::   IntercellFlux      => BaseClass_IntercellFlux
         procedure ::   QDotFaceLoop       => BaseClass_QDotFaceLoop
         procedure ::   dQFaceLoop         => BaseClass_dQFaceLoop
         procedure ::   QDotVolumeLoop     => BaseClass_QDotVolumeLoop
         procedure ::   dQVolumeLoop       => BaseClass_dQVolumeLoop
         procedure ::   ComputeInnerFluxes => BaseClass_ComputeInnerFluxes
         procedure ::   ComputeFaceFluxes  => BaseClass_ComputeFaceFluxes
         procedure ::   Describe           => ViscousMethod_describe
   end type ViscousMethod_t
!  *******************************************************
!  ---------------- Interior penalty method --------------
!  *******************************************************
   type, extends(ViscousMethod_t) :: IPMethod_t
      character(len=STR_LEN_VISCOUS) :: subType
      real(kind=RP)                  :: sigma0
      real(kind=RP)                  :: sigma1
      real(kind=RP)                  :: epsilon
      contains
         procedure :: QDotFaceLoop   => IP_QDotFaceLoop
         procedure :: QDotVolumeLoop => IP_QDotVolumeLoop
         procedure :: dQFaceLoop     => IP_dQFaceLoop
         procedure :: dQVolumeLoop   => IP_dQVolumeLoop
   end type IPMethod_t
!  *******************************************************
!  ---------------- Bassy-Rebay 1 method -----------------
!  *******************************************************
   type, extends(ViscousMethod_t) ::  BR1Method_t
      contains
         procedure ::  QDotFaceLoop   => BR1_QDotFaceLoop
         procedure ::  QDotVolumeLoop => BR1_QDotVolumeLoop
         procedure ::  dQFaceLoop     => BR1_dQFaceLoop
         procedure ::  dQVolumeLoop   => BR1_dQVolumeLoop
         procedure ::  IntercellFlux  => BR1_IntercellFlux
   end type BR1Method_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(ViscousMethod_t) ::  LDGMethod_t

   end type LDGMethod_t
!  *******************************************************
!
!  ========
   contains
!  ========
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

         else

            write(STD_OUT , * ) "Method ",trim(Setup % viscous_discretization)," not implemented yet."
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

              else if ( trim ( Setup % IPMethod ) .eq. "NIPG") then
      
                  ViscousMethod % subType = "NIPG"
                  ViscousMethod % epsilon = 1.0_RP

              else if ( trim ( Setup % IPMethod ) .eq. "IIPG") then

                  ViscousMethod % subType = "IIPG"
                  ViscousMethod % epsilon = 0.0_RP

              else

                  write(STD_OUT , *) "Method ",trim(Setup % IPMethod)," not implemented yet."
                  stop "Stopped."

              end if
      
            type is (BR1Method_t)
!           ------------------------------------
!              No extra parameters are needed
!           ------------------------------------
            class default

                write(STD_OUT , *) "Second order method allocation went wrong."
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
      subroutine BaseClass_QDotFaceLoop( self , edge ) 
         implicit none
         class(ViscousMethod_t)          :: self
         class(Edge_t)                       :: edge
!
!        ------------------------------------
!           The base class does nothing.
!        ------------------------------------
!
      end subroutine BaseClass_QDotFaceLoop

      subroutine BaseClass_dQFaceLoop( self , edge )  
         implicit none
         class(ViscousMethod_t)          :: self
         class(Edge_t)                       :: edge
!
!        ------------------------------------
!           The base class does nothing.
!        ------------------------------------
!
      end subroutine BaseClass_dQFaceLoop
      
      subroutine BaseClass_QDotVolumeLoop( self , element ) 
         implicit none
         class(ViscousMethod_t)          :: self
         class(QuadElement_t)                       :: element
!
!        ------------------------------------
!           The base class does nothing.
!        ------------------------------------
!
      end subroutine BaseClass_QDotVolumeLoop

      subroutine BaseClass_dQVolumeLoop( self , element )  
         use MatrixOperations
         implicit none
         class(ViscousMethod_t)          :: self
         class(QuadElement_t)                       :: element
!
!        ------------------------------------
!           The base class does nothing.
!        ------------------------------------
!
      end subroutine BaseClass_dQVolumeLoop

      subroutine BaseClass_ComputeInnerFluxes( self , element ) 
!
!        **********************************************************************
!              This subroutine computes the contravariant fluxes of the element
!           The fluxes read:
!                 F <- F * Ja(1,1) + G * Ja(2,1)
!                 G <- F * Ja(1,2) + G * Ja(2,2)
!        **********************************************************************
!
         implicit none  
         class(ViscousMethod_t)    :: self
         class(QuadElement_t)       :: element
!        -------------------------------------------------------------
         real(kind=RP)              :: F(0:element % spA % N,0:element % spA % N,1:NCONS,1:NDIM)
         integer                    :: eq

         associate( N => element % spA % N )

         F = ViscousFlux( N , element % W , element % dQ )

         do eq = 1 , NCONS
!        
!           F flux (contravariant)
!           ----------------------
            element % F(0:N,0:N,eq,IX) = F(0:N,0:N,eq,IX) * element % Ja(0:N,0:N,1,1) + F(0:N,0:N,eq,IY) * element % Ja(0:N,0:N,2,1)
!        
!           G flux (contravariant)
!           ----------------------
            element % F(0:N,0:N,eq,IY) = F(0:N,0:N,eq,IX) * element % Ja(0:N,0:N,1,2) + F(0:N,0:N,eq,IY) * element % Ja(0:N,0:N,2,2)
         end do

         end associate
         
      end subroutine BaseClass_ComputeInnerFluxes

      subroutine BaseClass_ComputeFaceFluxes ( self , edge )
!
!        **********************************************************************
!              This subroutine computes the viscous normal fluxes at both
!           sides of the face
!        **********************************************************************
!
         implicit none
         class(ViscousMethod_t)        :: self
         class(Edge_t)          :: edge
!        ------------------------------------------------
         integer                       :: side

         associate ( N => edge % spA % N )

         select type ( edge )

            type is (Edge_t)

               do side = 1 , 2
                  edge % F ( 0:N , 1:NCONS , side ) = ViscousNormalFlux( N , edge % w(0:N , 1:NPRIM , side) , edge % dQ(0:N , 1:NDIM , 1:NGRAD , side) , edge % dS )

               end do

            type is (StraightBdryEdge_t)
               edge % F ( 0:N , 1:NCONS , 1 ) = ViscousNormalFlux( N , edge % w(0:N , 1:NPRIM , 1) , edge % dQ(0:N , 1:NDIM , 1:NGRAD , 1) , edge % dS )

            type is (CurvedBdryEdge_t)
               edge % F ( 0:N , 1:NCONS , 1 ) = ViscousNormalFlux( N , edge % w(0:N , 1:NPRIM , 1) , edge % dQ(0:N , 1:NDIM , 1:NGRAD , 1) , edge % dS)

         end select

         end associate


      end subroutine BaseClass_ComputeFaceFluxes

      function BaseClass_IntercellFlux ( self , edge ) result ( Fstar )
         implicit none
         class( ViscousMethod_t ) :: self
         class( Edge_t )          :: edge
         real(kind=RP)            :: Fstar(0:edge % spA % N,1:NCONS)


      end function BaseClass_IntercellFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           INTERIOR PENALTY PROCEDURES
!           ---------------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
     subroutine IP_dQFaceLoop( self , edge ) 
         use MatrixOperations
         implicit none
         class(IPMethod_t)             :: self
         class(Edge_t)                 :: edge

     end subroutine IP_dQFaceLoop

     subroutine IP_dQVolumeLoop( self , element )
         implicit none
         class(IPMethod_t)             :: self
         class(QuadElement_t)            :: element

     end subroutine IP_dQVolumeLoop

     subroutine IP_QDotFaceLoop( self , edge ) 
         use MatrixOperations
         implicit none
         class(IPMethod_t)             :: self
         class(Edge_t)                 :: edge

     end subroutine IP_QDotFaceLoop

     subroutine IP_QDotVolumeLoop( self , element ) 
         use MatrixOperations
         implicit none
         class(IPMethod_t)             :: self
         class(QuadElement_t)                 :: element

     end subroutine IP_QDotVolumeLoop
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           BASSY-REBAY 1 PROCEDURES
!           ------------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
     subroutine BR1_dQFaceLoop( self , edge ) 
         implicit none
         class(BR1Method_t)             :: self
         class(Edge_t)                  :: edge
         real(kind=RP)                  :: ustar(0:edge % spA % N,1:NGRAD)
         
         associate ( N => edge % spA % N ) 
        
         select type ( edge ) 

            type is ( Edge_t ) 

               uStar(0:N,IGU) = 0.5_RP * ( edge % W(0:N,IU,LEFT) + edge % W(0:N,IU,RIGHT) )
               uStar(0:N,IGV) = 0.5_RP * ( edge % W(0:N,IV,LEFT) + edge % W(0:N,IV,RIGHT) )
               uStar(0:N,IGT) = 0.5_RP * ( edge % W(0:N,IT,LEFT) + edge % W(0:N,IT,RIGHT) )
               
               associate ( dQ => edge % quads(LEFT) % e % dQ )
                  dQ = dQ + dQFaceContribution( edge , LEFT , uStar )
               end associate
               
               associate ( dQ => edge % quads(RIGHT) % e % dQ ) 
                  dQ = dQ - dQFaceContribution( edge , RIGHT , uStar )
               end associate

            type is ( StraightBdryEdge_t )

               uStar(0:N,IGU) = 0.5_RP * ( edge % W(0:N,IU,1) + edge % wB(0:N,IU) )
               uStar(0:N,IGV) = 0.5_RP * ( edge % W(0:N,IV,1) + edge % wB(0:N,IV) )    
               uStar(0:N,IGT) = 0.5_RP * ( edge % W(0:N,IT,1) + edge % wB(0:N,IT) )
         
               associate ( dQ => edge % quads(1) % e % dQ ) 

               if ( .not. edge % inverted ) then
!
!                 If the normal points towards the domain exterior
!                 ------------------------------------------------
                  dQ = dQ + dQFaceContribution( edge , 1 , uStar )

               else
!
!                 If the normal points towards the domain interior
!                 ------------------------------------------------
                  dQ = dQ - dQFaceContribution( edge , 1 , uStar )

               end if
               
               end associate
!
   
            type is ( CurvedBdryEdge_t )

               uStar(0:N,IGU) = 0.5_RP * ( edge % W(0:N,IU,1) + edge % wB(0:N,IU) )
               uStar(0:N,IGV) = 0.5_RP * ( edge % W(0:N,IV,1) + edge % wB(0:N,IV) )    
               uStar(0:N,IGT) = 0.5_RP * ( edge % W(0:N,IT,1) + edge % wB(0:N,IT) )

               associate ( dQ => edge % quads(1) % e % dQ )

               if ( .not. edge % inverted ) then
!
!                 If the normal points towards the domain exterior
!                 ------------------------------------------------
                  dQ = dQ + dQFaceContribution( edge , 1 , uStar )

               else
!
!                 If the normal points towards the domain interior
!                 ------------------------------------------------
                  dQ = dQ - dQFaceContribution( edge , 1 , uStar )

               end if
 
               end associate

            class default

         end select

         end associate
         
     end subroutine BR1_dQFaceLoop

     subroutine BR1_dQVolumeLoop( self , element ) 
         use MatrixOperations
         implicit none
         class(BR1Method_t)   :: self
         class(QuadElement_t) :: element
         integer              :: iDim
         integer              :: iVar
         integer              :: which(NDIM)
         integer, parameter   :: PrimVariable(3) = [IU,IV,IT]
         real(kind=RP)        :: JaTimesW(0:element % spA % N , 0:element % spA % N)

         associate( N => element % spA % N , W => element % W , dQ => element % dQ , MD => element % spA % MD , &
                    weights => element % spA % w , M => element % spA % M , gm1 => Thermodynamics % gm1)


         do iDim = 1 , NDIM
   
            do iVar = 1 , NGRAD
               JaTimesW = element % Ja(0:N,0:N,iDim,1) * W(0:N,0:N,PrimVariable(iVar))
               call Mat_x_Mat(A = -MD , &
                     B = MatrixByVectorInIndex_F( JaTimesW , weights , N+1 , N+1 , 2 ) , & 
                     C = dQ(0:N,0:N,iDim,iVar) , & 
                     trA = .true. , reset = .false. )

               JaTimesW = element % Ja(0:N,0:N,iDim,2) * W(0:N,0:N,PrimVariable(iVar))
               call Mat_x_Mat(A = MatrixByVectorInIndex_F( JaTimesW , weights , N+1 , N+1 , 1) , &
                     B = -MD , C = dQ(0:N,0:N,iDim,iVar) , &
                     reset = .false. )

            end do
         end do

         end associate


     end subroutine BR1_dQVolumeLoop

     subroutine BR1_QDotFaceLoop( self , edge )
!
!        **************************************************************
!              This routine computes the edge loop according to the
!           "Form I" formulation:
!                 QDot +-= \int_e F^* l_j l_i ds
!           in where F^* represents the F·n product, that is, the 
!           result is added to the LEFT element, and substracted to 
!           the RIGHT element. 
!        **************************************************************
!
         use MatrixOperations
         implicit none
         class(BR1Method_t)    :: self
         class(Edge_t)     :: edge
!        -------------------------------------------------------
         real(kind=RP)        :: Fstar(0:edge % spA % N,1:NCONS)

         associate ( N => edge % spA % N )
!
!        Compute face fluxes
!        -------------------
         call self % ComputeFaceFluxes ( edge )
!
!        Compute the edge Riemann Flux is proceeds, or uses the prescribed boundary flux
!        -------------------------------------------------------------------------------
         Fstar = self % IntercellFlux( edge )
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
!              The obtained term is added to the LEFT element
!              ----------------------------------------------
               associate ( QDot => edge % quads(LEFT) % e % QDot )
                  QDot = QDot + QDotFaceContribution( edge , LEFT , Fstar )
               end associate
!
!              The obtained term is substracted to the RIGHT element
!              -----------------------------------------------------
               associate ( QDot => edge % quads(RIGHT) % e % QDot ) 
                  QDot = QDot - QDotFaceContribution( edge , RIGHT , Fstar )
               end associate
!
!           --------------------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )

                  if ( .not. edge % inverted ) then
!
!                    If the normal points towards the domain exterior
!                    ------------------------------------------------
                     QDot = QDot + QDotFaceContribution( edge , 1 , Fstar )
                  else
!
!                    If the normal points towards the domain interior
!                    ------------------------------------------------
                     QDot = QDot - QDotFaceContribution( edge , 1 , Fstar )
                  end if
               
               end associate
!
!           --------------------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate ( QDot => edge % quads(1) % e % QDot )
!
!                 The normal for curved elements always points towards the domain exterior
!                 ------------------------------------------------------------------------
                  QDot = QDot + QDotFaceContribution( edge , 1 , Fstar ) 
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
 
     end subroutine BR1_QDotFaceLoop

     subroutine BR1_QDotVolumeLoop( self , element ) 
!
!        *******************************************************************************
!           This routine computes the standard DG volumetric terms according to:
!                 QDot -= tr(D) M F M + M G M D         
!           The details about this matricial form is shown in the doc HiODG2DTech.pdf
!        *******************************************************************************
!
         use MatrixOperations
         implicit none
         class(BR1Method_t) :: self
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
!           F Loop
!           ------
            call Mat_x_Mat(A = -MD ,B = MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IX) , w , N+1 , N+1 , 2 ) , C=QDot(0:N,0:N,eq) , &
                        trA = .true. , reset = .false. )

!
!           G Loop
!           ------
            call Mat_x_Mat(A = MatrixByVectorInIndex_F( element % F(0:N,0:N,eq,IY) , w , N+1 , N+1 , 1) , B = -MD , C=QDot(0:N,0:N,eq) , &
                         reset = .false. )

         end do

         end associate

     end subroutine BR1_QDotVolumeLoop

     function BR1_IntercellFlux ( self , edge ) result ( Fstar )
         implicit none
         class(BR1Method_t) :: self
         class(Edge_t)      :: edge
         real(kind=RP)      :: Fstar(0:edge % spA % N,1:NCONS)
!        -------------------------------------------------------
!
!        ********************
         select type ( edge )
!        ********************
!
!           -----------------------------------------------------------------
            type is (StraightBdryEdge_t)

               associate( N => edge % spA % N )

               Fstar ( 0:N , 1:NCONS ) = 0.5_RP * ( edge % F(0:N , 1:NCONS , 1 ) + viscousNormalFlux( N , edge % wB( 0:N , 1:NPRIM ) , edge % gB( 0:N , 1:NDIM , 1:NGRAD ) , edge % dS) ) 

               end associate
!               
!           -----------------------------------------------------------------
            type is (CurvedBdryEdge_t)

               associate( N => edge % spA % N )

               Fstar ( 0:N , 1:NCONS ) = 0.5_RP * ( edge % F(0:N , 1:NCONS , 1 ) + viscousNormalFlux( N , edge % wB( 0:N , 1:NPRIM ) , edge % gB( 0:N , 1:NDIM , 1:NGRAD ) , edge % dS ) ) 

               end associate
!
!           -----------------------------------------------------------------
            type is (Edge_t)
      
               associate( N => edge % spA % N )

               Fstar ( 0:N , 1:NCONS ) = 0.5_RP * ( edge % F(0:N , 1:NCONS , LEFT ) + edge % F( 0:N , 1:NCONS , RIGHT ) ) 

               end associate
!
!           ------------------------------------------------------------------
            class default
!
!        **********
         end select
!        **********
!
   end function BR1_IntercellFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           AUXILIAR SUBROUTINES
!           --------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      function dQFaceContribution( edge , loc , ustar ) result ( duF )
         use MatrixOperations
         implicit none
         class(Edge_t)              :: edge
         integer                    :: loc
         real(kind=RP)              :: ustar(0:,:)
         real(kind=RP)              :: duF(0:edge % spA % N,0:edge % spA % N,1:NDIM,1:NGRAD)
!        ---------------------------------------------------------------------
         real(kind=RP)          :: uTimesW(0:edge % spA % N,1:NGRAD)
         real(kind=RP)          :: auxdS(1:NDIM,0:edge % spA % N)
         real(kind=RP), pointer :: lj(:)
         integer                :: direction
         integer                :: pos
         integer                :: index
         integer                :: i , j , var , iDim
!
         associate( N=> edge % quads(loc) % e % spA % N, &
             e=> edge % quads(loc) % e)
 
         select case (edge % edgeLocation(loc))

            case (ERIGHT)
   
               direction = e % edgesDirection( edge % edgeLocation(loc) )
               pos = RIGHT                    ! Where to interpolate the test function
               index = iX                     ! The coordinate (xi/eta) in which the boundary is located
               
               if ( direction .eq. FORWARD ) then
                  uTimesW(0:N,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1 )
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1 )
                  auxdS(1:NDIM,N:0:-1)      = edge % dS(1:NDIM,0:N)

               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            case (ETOP)
      
               direction = - e % edgesDirection ( edge % edgeLocation(loc) ) 
               pos = RIGHT
               index = iY
   
               if ( direction .eq. FORWARD ) then
                  uTimesW(0:N,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1) 
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1 ) 
                  auxdS(1:NDIM,N:0:-1)      = edge % dS(1:NDIM,0:N)
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            case (ELEFT)
   
               direction = - e % edgesDirection ( edge % edgeLocation(loc) )
               pos = LEFT
               index = iX
   
               if ( direction .eq. FORWARD ) then
                  uTimesW(0:N,1:NGRAD) =  MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1 ) 
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1)  
                  auxdS(1:NDIM,N:0:-1)      = edge % dS(1:NDIM,0:N)
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            case (EBOTTOM)

               direction = e % edgesDirection ( edge % edgeLocation(loc) ) 
               pos = LEFT
               index = iY
   
               if ( direction .eq. FORWARD ) then
                  uTimesW(0:N,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1)  
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , N+1 , NGRAD , 1)  
                  auxdS(1:NDIM,N:0:-1)      = edge % dS(1:NDIM,0:N)
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            end select

             lj(0:N) => e % spA % lb(0:N,pos)

            select case (index)
            
               case (IY)
   
                  forall (var=1:NGRAD,iDim=1:NDIM,j=0:N,i=0:N)
                     duF(i,j,iDim,var) = uTimesW(i,var) * lj(j) * auxdS(iDim , i)
                  end forall

               case (IX)
   
                  forall (var=1:NGRAD,iDim=1:NDIM,j=0:N,i=0:N)
                     duF(i,j,iDim,var) = uTimesW(j,var) * lj(i) * auxdS(iDim , j)
                  end forall
            end select

            lj=>NULL()
   
         end associate

      end function dQFaceContribution

      function QDotFaceContribution( edge , loc , Fstar ) result ( dFJ )
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
         real(kind=RP)              :: dFJ(0:edge % spA % N,0:edge % spA % N,1:NCONS)
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
                  Fstar2D(1 , 0:N    , 1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS )

               elseif ( direction .eq. BACKWARD ) then
!
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(1 , N:0:-1 , 1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS ), rowC = N+1 , colC = NCONS )

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
                  Fstar2D(0:N,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(N:0:-1,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

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
                  Fstar2D(1,0:N,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(1,N:0:-1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

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
                  Fstar2D(0:N,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

               elseif ( direction .eq. BACKWARD ) then
!        
!                 Introduce the result in the opposite order
!                 ------------------------------------------
                  Fstar2D(N:0:-1,1,1:NCONS) = Mat_x_Mat_F( e % spA % M , Fstar(0:N,1:NCONS), rowC = N+1 , colC = NCONS ) 

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
         dFJ =  MatrixMultiplyInIndex_F( Fstar2D , lj2D , N+1 , N+1 , NCONS , index ) 
!
!        Free the variables
!        ------------------
         lj2D=>NULL()
         deallocate(Fstar2D)

         end associate

      end function QDotFaceContribution

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

#endif
