#include "Defines.h"

program main
    use SMConstants
    use DGSEM_Class
    use Setup_class
!    use QuadMeshClass
!    use QuadElementClass
    use Headers
    use ChecksModule
    implicit none
    type(DGSEM_t)          :: sem
    integer                :: eID , edID
   real(kind=RP)  :: tstart , tend

#ifdef NAVIER_STOKES
    call Main_Header("Compressible Navier-Stokes equations high-order discontinuous Galerkin CFD 2D Solver")
#else
    call Main_Header("Compressible Euler equations high-order discontinuous Galerkin CFD 2D Solver")
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

    call cpu_time(tstart)
    call sem % Integrate()
    call cpu_time(tend)
    print*, "Time: ",tend-tstart
!
!  Program finished
!  ----------------
   write(STD_OUT , '(/,/,30X,A)') "\x1B[1;32m ****************** \x1B[0m"
   write(STD_OUT , '(30X,A)' ) "\x1B[1;32m Program finished! \x1B[0m"
   write(STD_OUT , '(30X,A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

end program main



