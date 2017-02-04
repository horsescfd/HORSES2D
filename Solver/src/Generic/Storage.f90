module Storage_module
    use SMConstants
    implicit none

    private
    public  Storage_t , newStorage

    type Storage_t
        real(kind=RP), pointer      :: Q(:)
        real(kind=RP), pointer      :: QDot(:)
        real(kind=RP), pointer      :: F(:)
        real(kind=RP), pointer      :: dQ(:)
        real(kind=RP), pointer      :: W(:)
    end type Storage_t

    contains
        function newStorage() result(val)
            implicit none
            type(Storage_t)         :: val

            val % Q    => NULL()
            val % QDot => NULL()
            val % F    => NULL()
            val % dQ   => NULL()
            val % W    => NULL()

        end function newStorage

end module Storage_module
