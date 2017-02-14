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
         procedure ::   QDotFaceLoop   => BaseClass_QDotFaceLoop
         procedure ::   dQFaceLoop     => BaseClass_dQFaceLoop
         procedure ::   QDotVolumeLoop => BaseClass_QDotVolumeLoop
         procedure ::   dQVolumeLoop   => BaseClass_dQVolumeLoop
         procedure ::   Describe       => ViscousMethod_describe
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
         real(kind=RP), allocatable     :: ustar(:,:)
         
         associate ( N => edge % spA % N ) 
        
         allocate( uStar ( 0 : N , NGRAD ) )

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

               uStar(0:N,IGU) = edge % W(0:N,IU,1)
               uStar(0:N,IGV) = edge % W(0:N,IV,1)
               uStar(0:N,IGT) = edge % W(0:N,IT,1)
         
               associate ( dQ => edge % quads(1) % e % dQ ) 
                  dQ = dQ + dQFaceContribution( edge , 1 , uStar )
               end associate
   
            type is ( CurvedBdryEdge_t )

               uStar(0:N,IGU) = edge % W(0:N,IU,1)
               uStar(0:N,IGV) = edge % W(0:N,IV,1)
               uStar(0:N,IGT) = edge % W(0:N,IT,1)

               associate ( dQ => edge % quads(1) % e % dQ )
                  dQ = dQ + dQFaceContribution( edge , 1 , uStar )
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

         associate( N => element % spA % N , W => element % W , dQ => element % dQ , MD => element % spA % MD , &
                    weights => element % spA % w , M => element % spA % M , gm1 => Thermodynamics % gm1)


         do iDim = 1 , NDIM
   
            do iVar = 1 , NGRAD

               which = [iDim , 1]
               call Mat_x_Mat(A = -MD , &
                     B = MatrixByVectorInIndex_F( element % Ja(which) * W(0:N,0:N,PrimVariable(iVar)) , weights , 2 ) , & 
                     C = dQ(0:N,0:N,iDim,iVar) , & 
                     trA = .true. , reset = .false. )

               which = [iDim , 2]
               call Mat_x_Mat(A = MatrixByVectorInIndex_F( element % Ja(which) * W(0:N,0:N,PrimVariable(iVar)) , weights , 1) , &
                     B = -MD , C = dQ(0:N,0:N,iDim,iVar) , &
                     reset = .false. )

            end do
         end do

         end associate


     end subroutine BR1_dQVolumeLoop

     subroutine BR1_QDotFaceLoop( self , edge ) 
         use MatrixOperations
         implicit none
         class(BR1Method_t)                     :: self
         class(Edge_t)                          :: edge
         
     end subroutine BR1_QDotFaceLoop

     subroutine BR1_QDotVolumeLoop( self , element ) 
         use MatrixOperations
         implicit none
         class(BR1Method_t)             :: self
         class(QuadElement_t)                 :: element

     end subroutine BR1_QDotVolumeLoop
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
         real(kind=RP), allocatable :: duF(:,:,:,:)
!        ---------------------------------------------------------------------
         real(kind=RP), allocatable         :: uTimesW(:,:)
         real(kind=RP), allocatable         :: auxdS(:,:)
         real(kind=RP), pointer             :: lj(:)
         integer                            :: direction
         integer                            :: pos
         integer                            :: index
         integer                            :: i , j , var , iDim
!

         associate( N=> edge % quads(loc) % e % spA % N, &
             e=> edge % quads(loc) % e)
 
         allocate(uTimesW(0:N,NGRAD))
         allocate(auxdS(NDIM,0:N))

         select case (edge % edgeLocation(loc))

            case (ERIGHT)
   
               direction = e % edgesDirection( edge % edgeLocation(loc) )
               pos = RIGHT                    ! Where to interpolate the test function
               index = iX                     ! The coordinate (xi/eta) in which the boundary is located

               
               if ( direction .eq. FORWARD ) then
                  uTimesW(0:N,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1 )
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1 )
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
                  uTimesW(0:N,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1) 
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1 ) 
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
                  uTimesW(0:N,1:NGRAD) =  MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1 ) 
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1)  
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
                  uTimesW(0:N,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1)  
                  auxdS(1:NDIM,0:N)      = edge % dS(1:NDIM,0:N)
               elseif ( direction .eq. BACKWARD ) then
                  uTimesW(N:0:-1,1:NGRAD) = MatrixByVectorInIndex_F ( ustar(0:N,1:NGRAD) , edge % spA % w , 1)  
                  auxdS(1:NDIM,N:0:-1)      = edge % dS(1:NDIM,0:N)
               else
                  print*, "Direction not forward nor backward"
                  stop "Stopped."
               end if

            end select

             lj(0:N) => e % spA % lb(0:N,pos)

             allocate ( duF(0:N,0:N,NDIM,NGRAD) )

            if ( index .eq. IY ) then
               do var = 1 , NGRAD
                  do iDim = 1 , NDIM
                     do j = 0 , N
                        do i = 0 , N
                           duF(i,j,iDim,var) = uTimesW(i,var) * lj(j) * auxdS(iDim , i)
                        end do
                     end do
                  end do
               end do

            elseif ( index .eq. IX ) then
               do var = 1 , NGRAD
                  do iDim = 1 , NDIM
                     do j = 0 , N
                        do i = 0 , N
                           duF(i,j,iDim,var) = uTimesW(j,var) * lj(i) * auxdS(iDim , j)
                        end do
                   end do
                  end do
               end do
            end if

            lj=>NULL()
   
            deallocate ( uTimesW )
         end associate

      end function dQFaceContribution

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
