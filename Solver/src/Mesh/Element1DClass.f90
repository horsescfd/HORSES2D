module Element1DClass
    use SMConstants
    use NodeClass
    use NodesAndWeights_class
    use Storage_module
    implicit none

    private
    public  Element1D_t , Element1D_p 

    type Element1D_t
        class(Node_p) , pointer           :: nodes(:)
        integer                           :: ID
        integer       , dimension(2)      :: facesID
        integer                           :: address
        !       Storage
        class(NodesAndWeights_t), pointer :: Interp
        class(NodesAndWeights_t), pointer :: spI
        real(kind=RP)                     :: hdiv2
        real(kind=RP), allocatable        :: x(:)
        real(kind=RP), pointer            :: Q(:,:) , QDot(:,:) , F(:,:) , dQ(:,:)
        real(kind=RP), pointer            :: Qb(:,:), dQb(:,:)
#ifdef ADVECTION ! -------------------------------------------------------------
        real(kind=RP), allocatable        :: A(:)
#endif ! -----------------------------------------------------------------------
        contains
            procedure   :: construct => constructElement
            procedure   :: SetStorage => Element1D_SetStorage
    end type Element1D_t

    type Element1D_p
        type(Element1D_t),  pointer     :: e
    end type Element1D_p

    contains
        subroutine ConstructElement(self , ID , leftNode , rightNode , faceLeftID , faceRightID , N , nodes , spA , address , storage , spI)
             use Setup_class
             use Physics
             class(Element1D_t)                :: self
             integer                           :: ID
             class(Node_t), pointer            :: leftNode
             class(Node_t), pointer            :: rightNode
             integer                           :: faceLeftID
             integer                           :: faceRightID
             integer                           :: nodes
             integer                           :: N
             class(NodalStorage)               :: spA
             integer                           :: address
             class(Storage_t)                  :: storage
             class(NodesAndWeights_t), pointer :: spI
!
             self % ID = ID
!
!            ************************
!            Point to neighbour nodes
!            ************************
!
             allocate( self % nodes(2) )
             self % nodes(LEFT) % n => leftNode
             self % nodes(RIGHT) % n => rightNode
!
!            **************
!            Get faces data
!            **************
!
             self % facesID(LEFT) = faceLeftID
             self % facesID(RIGHT) = faceRightID
!
!            *********************
!            Get NodalStorage data       
!            *********************
!
             call spA % add( N , nodes , self % Interp )
             self % spI => spI
!
!            ************
!            Linking data
!            ************
!            ---------------------------------------
!               These two arrays were before 
!               allocated inside each element.
!               Now they are just linked.
!                   allocate ( self % Q    ( 0:N , NEC )  ) 
!                   allocate ( self % QDot ( 0:N , NEC )  ) 
!                   allocate ( self % F    ( 0:N , NEC )  ) 
!                   allocate ( self % dQ   ( 0:N , NEC )  ) 
!            ---------------------------------------
!
             self % address = address
             call self % SetStorage( storage )
!
!            *************
!            Allocate data
!            *************
!
             allocate ( self % x    ( 0:N       )  ) 
             allocate ( self % Qb   ( NEC   , 2 )  ) 
             allocate ( self % dQb  ( NEC   , 2 )  ) 
!
!            **********
!            Set values
!            **********
!
             self % x = 0.5_RP * (self % nodes(LEFT) % n % x + self % nodes(RIGHT) % n % x + self % Interp % xi * (self % nodes(RIGHT) % n % x - self % nodes(LEFT) % n % x ) )
            ! TODO
             self % hdiv2 = 0.0_RP !0.5_RP*abs(self % nodes(LEFT) % n % x - self % nodes(RIGHT) % n % x)
             
#ifdef ADVECTION 
             allocate( self % A(0:N) )
#endif
        end subroutine ConstructElement

        subroutine Element1D_SetStorage( self , storage )
            use SMConstants
            use Storage_module
            use Setup_class
            use Physics
            implicit none
            class(Element1D_t)      :: self
            class(Storage_t)        :: storage
        

            associate ( N => self % Interp % N )
             self % Q    ( 0:N , 1:NEC )  => storage % Q    ( self % address: ) 
             self % QDot ( 0:N , 1:NEC )  => storage % QDot ( self % address: ) 

             if ( trim(Setup % inviscid_discretization) .eq. "Over-Integration" ) then
               self % F (0: self % spI % N , 1:NEC) => storage % F ( (self % ID -1)*(self % spI % N+1) + 1: self % ID * ( self % spI % N + 1) )

             else

               self % F    ( 0:N , 1:NEC )  => storage % F    ( self % address: ) 

             end if

             self % dQ   ( 0:N , 1:NEC )  => storage % dQ   ( self % address: ) 

            end associate

        end subroutine Element1D_SetStorage
            



end module Element1DClass
