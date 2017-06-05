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
#include "Defines.h"

program main
    use SMConstants
    use DGSEM_Class
    use Setup_class
    use QuadMeshClass
    use QuadElementClass
    use Headers
    use ChecksModule
    implicit none
    type(DGSEM_t) :: sem
    integer       :: exit_code

#ifdef NAVIER_STOKES
    call Main_Header("2D Compressible Navier-Stokes equations")
#else
    call Main_Header("2D Compressible Euler equations")
#endif
!
!   Read the case file and load parameters on Setup class
!   -----------------------------------------------------
    call setup % Initialization()
!
!   Initialize the Physics module
!   -----------------------------
    call InitializePhysics
!
!   Initialize and build the DGSem structure
!   ----------------------------------------
    sem = DGSEM_Initialize()
    call sem % construct
!
!   Perform checks on the built framework
!   -------------------------------------
    call checks( sem )
!
!   Time integration
!   ----------------
    call sem % Integrate()
!
!   Check with an analytical condition
!   ----------------------------------
    exit_code = sem % Finalize()
!
!  Program finished
!  ----------------
   write(STD_OUT , '(/,/,30X,A)') "\x1B[1;32m ****************** \x1B[0m"
   write(STD_OUT , '(30X,A)' ) "\x1B[1;32m Program finished! \x1B[0m"
   write(STD_OUT , '(30X,A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

   call exit(exit_code)

end program main



