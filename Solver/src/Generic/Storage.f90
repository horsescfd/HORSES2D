!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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
            self % Q = 0.0_RP 
            allocate ( self % QDot ( NCONS *         totalPolynomialOrder  )  ) 
            self % QDot = 0.0_RP
#ifdef NAVIER_STOKES
            allocate ( self % dQ   ( NCONS  * NDIM * totalPolynomialOrder  )  ) 
            self % dQ = 0.0_RP
#endif
        end subroutine Storage_AllocateMemory

end module Storage_module
