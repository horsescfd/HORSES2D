program main
    use nodesAndWeights_class
    use SMConstants
    use Physics
    use MeshFileClass
    use DGSEM_Class
    use Setup_class
    use Headers
    implicit none
    type(MeshFile_t)            :: meshFile
    type(DGSEM_t)             :: sem
    interface
      subroutine checks( sem )
         use DGSEM_Class
         use SMConstants
         use Physics
         use NodesAndWeights_class
         use Mesh1DClass
         use MeshFileClass
         use DGSpatialDiscretizationMethods
         use Storage_module
         use DGBoundaryConditions  
         implicit none
         class(DGSEM_t) :: sem
      end subroutine checks
    end interface


!   =====================
!   The REAL main program
!   =====================
!
    call Main_Header("High-order discontinuous Galerkin CFD 2D Solver")

    call InitializePhysics

!
!   **********************************************************
!   Read the mesh
!   **********************************************************
!
!   
!   ****************************************
!   Initialize and build the DGSem structure
!   ****************************************
!
!    sem = DGSEM_Initialize()
!    call sem % construct(meshFile)
!   
!   ***********************************************
!   Set the initial condition to all flow variables
!   ***********************************************
!
!    call sem % SetInitialCondition()

    
!    call checks( sem ) 

!    call sem % Integrate()

    write(STD_OUT , '(/,/,A)') "\x1B[1;32m ****************** \x1B[0m"
    write(STD_OUT , '(A)' ) "\x1B[1;32m Program finished! \x1B[0m"
    write(STD_OUT , '(A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

end program main



