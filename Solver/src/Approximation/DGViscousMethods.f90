module DGViscousMethods
   use SMConstants
   use QuadMeshClass
   use QuadElementClass
   use Physics
   use Setup_class
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
         procedure ::   QDotFaceLoop => BaseClass_QDotFaceLoop
         procedure ::   dQFaceLoop => BaseClass_dQFaceLoop
         procedure ::   QDotVolumeLoop => BaseClass_QDotVolumeLoop
         procedure ::   dQVolumeLoop => BaseClass_dQVolumeLoop
         procedure ::   Describe => ViscousMethod_describe
   end type ViscousMethod_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(ViscousMethod_t) :: IPMethod_t
      character(len=STR_LEN_VISCOUS)        :: subType
      real(kind=RP)                             :: sigma0
      real(kind=RP)                             :: sigma1
      real(kind=RP)                             :: epsilon
      contains
         procedure :: QDotFaceLoop => IP_QDotFaceLoop
         procedure :: QDotVolumeLoop => IP_QDotVolumeLoop
         procedure :: dQVolumeLoop => IP_dQVolumeLoop
   end type IPMethod_t
!  *******************************************************
!  -------------------------------------------------------
!  *******************************************************
   type, extends(ViscousMethod_t) ::  BR1Method_t
      contains
         procedure ::  QDotFaceLoop => BR1_QDotFaceLoop
         procedure ::  QDotVolumeLoop => BR1_QDotVolumeLoop
         procedure ::  dQFaceLoop => BR1_dQFaceLoop
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

      subroutine BaseClass_QDotFaceLoop( self , face ) 
         implicit none
         class(ViscousMethod_t)          :: self
         class(Edge_t)                       :: face
!
!        ------------------------------------
!           The base class does nothing.
!        ------------------------------------
!
      end subroutine BaseClass_QDotFaceLoop

      subroutine BaseClass_dQFaceLoop( self , face )  
         implicit none
         class(ViscousMethod_t)          :: self
         class(Edge_t)                       :: face
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
!        ---------------------------------------
!           The base class consists in a 
!           standard DG volume loop for the
!           solution: -int_ l'j u dx
!              val = -tr(MD)*U^{el}
!        ---------------------------------------
!
!         associate( dQ => element % dQ , &
!                    MD => element % Interp % MD, &
!                    Q  => element % Q)
!
!         dQ = dQ - TransposeMat_x_NormalMat_F( MD , Q )
!         
!
!         end associate
      end subroutine BaseClass_dQVolumeLoop
 
      subroutine ViscousMethod_describe( self )
         implicit none
         class(ViscousMethod_t)          :: self
!
!         write(STD_OUT , *) "Second order method description: "
!         write(STD_OUT , '(20X,A,A)') "Method: ", trim(self % method)
!        
!         select type (self)
!            type is (IPMethod_t)
!               write(STD_OUT , '(20X,A,A)') "Sub method: ", trim(self % subType)
!               write(STD_OUT,'(20X,A,F10.4)') "Sigma0: " , self % sigma0
!               write(STD_OUT,'(20X,A,F10.4)') "Sigma1: " , self % sigma1
!               write(STD_OUT,'(20X,A,F10.4)') "Epsilon: " , self % epsilon
!            
!            type is (BR1Method_t)
!!           ---------------------
!!              Nothing to add 
!!           ---------------------
!            class default
!
!         end select
!
      end subroutine ViscousMethod_describe

      subroutine average_dQFaceLoop( self , face )
         use MatrixOperations
         implicit none
         class(ViscousMethod_t)             :: self
         class(Edge_t)                          :: face
         real(kind=RP), dimension(NEC)          :: ustar
!!
!!        Compute the averaged flux
!         ustar = average_uFlux( face ) 
!!
!!        Perform the loop in both elements
!         select type ( face )
!            type is (Edge_t)
!               associate(dQ => face % elements(LEFT) % e % dQ, &
!                   N=> face % elements(LEFT) % e % Interp % N, &
!                   e=> face % elements(LEFT) % e)
!
!                  dQ = dQ + vectorOuterProduct(e % Interp % lb(:,RIGHT) , ustar) * face % n
!
!               end associate
!
!               associate(dQ => face % elements(RIGHT) % e % dQ, &
!                   N=> face % elements(RIGHT) % e % Interp % N, &
!                   e=> face % elements(RIGHT) % e)
!
!                  dQ = dQ + vectorOuterProduct(e % Interp % lb(:,LEFT) , ustar) * (-1.0_RP * face % n)
!
!               end associate
!
!
!            type is (StraightBdryEdge_t)
!               associate(dQ => face % elements(1) % e % dQ , &
!                  N => face % elements(1) % e % Interp % N , &
!                  e => face % elements(1) % e)
!
!                  if ( face % BCLocation .eq. LEFT ) then
!                     dQ = dQ + vectorOuterProduct(e % Interp % lb(:,face % BCLocation) , ustar)*face % n
!                  elseif (face % BCLocation .eq. RIGHT) then
!                     dQ = dQ + vectorOuterProduct(e % Interp % lb(:,face % BCLocation) , ustar)* face % n
!                  end if
!               end associate
!            class default
!         end select
!
     end subroutine average_dQFaceLoop
!
!    ************************************************
!           INTERIOR PENALTY FLUXES 
!    ************************************************
!
     subroutine IP_dQVolumeLoop( self , element )
         implicit none
         class(IPMethod_t)             :: self
         class(QuadElement_t)            :: element

!         element % dQ = element % dQ + matmul( element % Interp % MD , element % Q )

     end subroutine IP_dQVolumeLoop

     subroutine IP_QDotFaceLoop( self , face ) 
         use MatrixOperations
         implicit none
         class(IPMethod_t)             :: self
         class(Edge_t)                 :: face
         real(kind=RP), allocatable    :: T1(:,:) , T2(:,:) , J0(:,:) , J1(:,:)
!!
!!        -------------------------------------------------------
!!           Computes the interior penalty face loop. Consists
!!         in four terms:
!!           T1= {du n}[[v]]
!!           T2= -e{dv n}[[u]]
!!           J0= [[u]][[v]]
!!           J1= [[du]][[dv]]
!!        -------------------------------------------------------
!         select type ( face )
!            type is (Edge_t)
!                 associate( uL => face % elements(LEFT)  % e % Qb(:,RIGHT) , &
!                            uR => face % elements(RIGHT) % e % Qb(:,LEFT) , &
!                           duL => face % elements(LEFT)  % e % dQb(:,RIGHT) , &
!                           duR => face % elements(RIGHT) % e % dQb(:,LEFT) )
!!
!!                   -------------------------------------------------
!!                    Terms added to the LEFT elements equation
!!                   -------------------------------------------------
!!
!                    associate( eL => face % elements(LEFT) % e , &
!                               eR => face % elements(RIGHT) % e )
!                       allocate( T1( 0:eL % Interp % N , NEC ) , T2( 0:eL % Interp % N , NEC ) , J0 ( 0:eL % Interp % N , NEC ) , J1( 0 : eL % Interp % N , NEC ) )
!                       
!                       T1 = 0.5_RP * (face % n) * vectorOuterProduct( eL % Interp % lb(:,RIGHT) , duL + duR )
!                       T2 = -0.5_RP * self % epsilon * (1.0_RP / eL % hdiv2) * (face % n) * vectorOuterProduct( matmul( eL % Interp % DT , eL % Interp % lb(:,RIGHT) ) , uL - uR)
!                       J0 = -0.5_RP * (self % sigma0)* (1.0_RP / el % hdiv2) * vectorOuterProduct( eL % Interp % lb(:,RIGHT) , uL - uR )
!                       J1 = -(self % sigma1)* vectorOuterProduct( matmul( eL % Interp % DT , eL % Interp % lb(:,RIGHT) ) , duL - duR )
!!                          TODO
!                       eL % QDot = eL % QDot! + viscousFlux(g = T1 + T2 + J0 + J1)
!
!                       deallocate( T1 , T2 , J0 , J1 )
!
!!
!!                      -------------------------------------------------
!!                         Terms added to the RIGHT element equation
!!                      -------------------------------------------------
!!
!                       allocate( T1( 0:eR % Interp % N , NEC ) , T2( 0:eR % Interp % N , NEC ) , J0 ( 0:eR % Interp % N , NEC ) , J1( 0 : eR % Interp % N , NEC ) )
!
!                       T1 = 0.5_RP * (-face % n) * vectorOuterProduct( eR % Interp % lb(:,LEFT) , duL + duR )
!                       T2 = -0.5_RP * self % epsilon * (1.0_RP / eR % hdiv2) * (face % n) * vectorOuterProduct( matmul( eR % Interp % DT , eR % Interp % lb(:,LEFT) ) , uL - uR)
!                       J0 = 0.5_RP * (self % sigma0)* (1.0_RP / el % hdiv2) * vectorOuterProduct( eR % Interp % lb(:,LEFT) , uL - uR )
!                       J1 = (self % sigma1)* vectorOuterProduct( matmul( eR % Interp % DT , eR % Interp % lb(:,LEFT) ) , duL - duR )
!!                                         TODO
!                       eR % QDot = eR % QDot !+ viscousFlux(g = T1 + T2 + J0 + J1)
!
!                       deallocate( T1 , T2 , J0 , J1 )
!
!                    end associate
!                    
!                 end associate
!
!            type is (StraightBdryEdge_t)
!               associate( uE => face % elements(1) % e % Qb(:,face % BCLocation) , &
!                          uB => face % uB , &
!                         duE => face % elements(1) % e % dQb(:,face % BCLocation) , &
!                         duB => face % gB , &
!                           e => face % elements(1) % e)
!
!                     allocate( T1( 0:e % Interp % N , NEC ) , T2( 0:e % Interp % N , NEC ) , J0 ( 0:e % Interp % N , NEC ) , J1( 0 : e % Interp % N , NEC ) )
!     
!                     T1 = 0.5_RP * (face % n) * vectorOuterProduct( e % Interp % lb(:,face % BCLocation ) , duE + duB )
!                     T2 = -0.5_RP * (self % epsilon) * (1.0_RP / e % hdiv2) * (face % n) * vectorOuterProduct( matmul( e % Interp % DT , e % Interp % lb(:,face % BCLocation) ) , uE - uB )
!                     J0 = -0.5_RP * (self % sigma0) * (1.0_RP / e % hdiv2) * vectorOuterProduct( e % Interp % lb(:,face % BCLocation) , uE - uB )
!                     J1 = -(self % sigma1) * vectorOuterProduct( matmul( e % Interp % DT , e % Interp % lb(:, face % BCLocation) ) , duE - duB )
!!                                      TODO 
!                     e % QDot = e % QDot !+ viscousFlux( T1 + T2 + J0 + J1 ) 
!
!                     deallocate( T1 , T2 , J0 , J1 )
!
!               end associate
!            class default
!         end select
     end subroutine IP_QDotFaceLoop

     subroutine IP_QDotVolumeLoop( self , element ) 
         use MatrixOperations
         implicit none
         class(IPMethod_t)             :: self
         class(QuadElement_t)                 :: element
!
!        -------------------------------------------
!           Computes the IP volume loop:
!              IPVol = - tr(D)*M*dQel 
!        -------------------------------------------
!                                      TODO
!         element % QDot = element % QDot! -  TransposeMat_x_NormalMat_F( element % Interp % MD , viscousFlux(g = element % dQ) )

     end subroutine IP_QDotVolumeLoop
!
!   ************************************************
!         BASSI-REBAY FLUXES
!   ************************************************
!
     subroutine BR1_dQFaceLoop( self , face ) 
         implicit none
         class(BR1Method_t)             :: self
         class(Edge_t)                 :: face

!         call average_dQFaceLoop( self , face )

     end subroutine BR1_dQFaceLoop

     subroutine BR1_QDotFaceLoop( self , face ) 
         use MatrixOperations
         implicit none
         class(BR1Method_t)                     :: self
         class(Edge_t)                          :: face
         real(kind=RP), dimension(NEC)          :: gstar
!!
!!        Compute the averaged flux
!         gstar = average_gFlux( face ) 
!!
!!        Perform the loop in both elements
!         select type ( face )
!            type is (Edge_t)
!               associate(QDot => face % elements(LEFT) % e % QDot, &
!                   N=> face % elements(LEFT) % e % Interp % N, &
!                   e=> face % elements(LEFT) % e)
!
!                  QDot = QDot + vectorOuterProduct(e % Interp % lb(:,RIGHT) , gstar) * face % n
!
!               end associate
!
!               associate(QDot => face % elements(RIGHT) % e % QDot, &
!                   N=> face % elements(RIGHT) % e % Interp % N, &
!                   e=> face % elements(RIGHT) % e)
!
!                  QDot = QDot + vectorOuterProduct(e % Interp % lb(:,LEFT) , gstar) * (-1.0_RP * face % n)
!
!               end associate
!
!
!            type is (StraightBdryEdge_t)
!               associate(QDot => face % elements(1) % e % QDot , &
!                  N => face % elements(1) % e % Interp % N , &
!                  e => face % elements(1) % e)
!
!                  QDot = QDot + vectorOuterProduct(e % Interp % lb(:,face % BCLocation) , gstar)*face % n
!               end associate
!            class default
!         end select
!

         
     end subroutine BR1_QDotFaceLoop

     subroutine BR1_QDotVolumeLoop( self , element ) 
         use MatrixOperations
         implicit none
         class(BR1Method_t)             :: self
         class(QuadElement_t)                 :: element

!!        Perform the matrix multiplication
!         associate( QDot => element % QDot , &
!                    MD => element % Interp % MD)
!!        TODO
!         QDot = QDot! - TransposeMat_x_NormalMat_F( MD , viscousFlux(g = element % dQ ))
!
!         end associate
!
!
!
     end subroutine BR1_QDotVolumeLoop


!
!    ****************************************************
!        Extra subroutines to compute fluxes
!    ****************************************************
!
     function average_uFlux( face ) result (val)
        implicit none
        class(Edge_t)                  :: face
        real(kind=RP), dimension(NEC)  :: val
        real(kind=RP), pointer         :: QL(:) , QR(:) , QBdry(:)
!!
!!       Compute the average of both states
!!
!        select type ( face )
!            type is (StraightBdryEdge_t)
!               
!                QBdry => face % elements(1) % e % Qb( : , face % BCLocation ) 
!               
!                val = 0.5_RP * (Qbdry + face % uB)
!
!                QBdry => NULL()
!            type is (Edge_t)
!                QL => face % elements(LEFT)  % e % Qb(: , RIGHT)
!                QR => face % elements(RIGHT) % e % Qb(: , LEFT )
!            
!                val = 0.5_RP * (QL + QR) 
!
!                QL => NULL()
!                QR => NULL()
!            class default
!        end select
!
     end function average_uFlux

     function average_gFlux( face ) result ( val )
         implicit none
         class(Edge_t)              :: face
         real(kind=RP), dimension(NEC)       :: val
         real(kind=RP), pointer        :: dQL(:) , dQR(:) , dQBdry(:)
!!
!!        Compute the average of both states
!!
!         select type ( face )
!            type is ( StraightBdryEdge_t )
!               
!               dQBdry => face % elements(1) % e % dQb( : , face % BCLocation )
!               ! TODO
!               val = 0.0_RP !0.5_RP * (viscousFlux(g=dQBdry) + viscousFlux(g=face % gB))
!
!               dQBdry => NULL()
!            
!            type is ( Edge_t )
!         
!               dQL => face % elements(LEFT) % e % dQb( : , RIGHT )
!               dQR => face % elements(RIGHT) % e %  dQb( : , LEFT ) 
!!              TODO
!               val = 0.0_RP !0.5_RP * (viscousFlux(g=dQL) + viscousFlux(g=dQR) )
!      
!               dQL => NULL()
!               dQR => NULL()
!
!            class default
!
!         end select
!
     end function average_gFlux
         




end module DGViscousMethods
