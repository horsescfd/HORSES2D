module Storage_module
    use SMConstants
    implicit none

#include "Defines.h"

    private
    public  Storage_t , newStorage

    type Storage_t
        real(kind=RP), pointer      :: Q(:)
        real(kind=RP), pointer      :: QDot(:)
#ifdef NAVIER_STOKES
        real(kind=RP), pointer      :: dQ(:)
#endif
        contains
            procedure   :: AllocateMemory    => Storage_AllocateMemory
    end type Storage_t

    contains
        function newStorage() result(val)
            implicit none
            type(Storage_t)         :: val

            val % Q    => NULL()
            val % QDot => NULL()
#ifdef NAVIER_STOKES
            val % dQ   => NULL()
#endif

        end function newStorage

        subroutine Storage_AllocateMemory( self , totalPolynomialOrder )
            implicit none
            class(Storage_t)     :: self
            integer, intent(in)  :: totalPolynomialOrder
!
!           Allocate memory for Q , QDot , and dQ
!              The sizes are the following:
!                 Q    -> NCONS * (N+1) * (N+1) * no_of_elements
!                 QDot -> NCONS * (N+1) * (N+1) * no_of_elements
!                 dQ   -> NGRAD * NDIM * (N+1) * (N+1) * no_of_elements
!           -------------------------------------------------------------
            allocate ( self % Q    ( NCONS *         totalPolynomialOrder  )  ) 
            allocate ( self % QDot ( NCONS *         totalPolynomialOrder  )  ) 
#ifdef NAVIER_STOKES
            allocate ( self % dQ   ( NCONS  * NDIM * totalPolynomialOrder  )  ) 
#endif
        end subroutine Storage_AllocateMemory

end module Storage_module
