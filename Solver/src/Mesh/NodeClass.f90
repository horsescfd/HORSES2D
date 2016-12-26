module nodeClass
     use SMConstants
     implicit none


     type node_t
          real(kind=RP)   :: x
          integer         :: ID
          contains
              procedure :: construct => constructNode
     end type node_t

     type node_p
          class(Node_t),  pointer       :: n
     end type Node_p

     contains
          subroutine constructNode(self , x , ID)
              implicit none
              class(Node_t)  :: self
              integer        :: ID
              real(kind=RP)  :: x
!
!             **************
!             Construct node 
!             **************
!
              self % ID = ID
              self % x  = x

          end subroutine constructNode



end module nodeClass
