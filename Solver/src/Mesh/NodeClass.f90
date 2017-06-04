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
module nodeClass
     use SMConstants
     use Physics
     implicit none

#include "Defines.h"
     private
     public Node_t , Node_p

     type node_t
          real(kind=RP), dimension(NDIM)   :: x
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
              class(Node_t)                  :: self
              integer                        :: ID
              real(kind=RP), dimension(NDIM) :: x
!
!             **************
!             Construct node 
!             **************
!
              self % ID = ID
              self % x  = x

          end subroutine constructNode



end module nodeClass
