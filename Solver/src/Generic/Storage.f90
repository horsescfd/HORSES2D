module Storage_module
    use SMConstants
    implicit none

    private
    public  Storage_t , newStorage

    type Storage_t
        real(kind=RP), pointer, contiguous      :: Q(:)
        real(kind=RP), pointer, contiguous      :: QDot(:)
        real(kind=RP), pointer, contiguous      :: F(:)
        real(kind=RP), pointer, contiguous      :: dQ(:)
    end type Storage_t

    contains
        function newStorage() result(val)
            implicit none
            type(Storage_t)         :: val

            val % Q => NULL()
            val % QDot => NULL()
            val % F => NULL()
            val % dQ => NULL()

        end function newStorage

end module Storage_module
