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
      integer           :: edge
    interface
      subroutine checks( sem )
         use DGSEM_Class
         use SMConstants
         use Physics
         use NodesAndWeights_class
         use QuadMeshClass
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
    call meshFile % read
!   
!   ****************************************
!   Initialize and build the DGSem structure
!   ****************************************
!
    sem = DGSEM_Initialize()
    call sem % construct(meshFile)
!   
!   ***********************************************
!   Set the initial condition to all flow variables
!   ***********************************************
!
    call sem % SetInitialCondition()

    do edge = 1 , sem % mesh % no_of_edges
      write(*,'(F24.16)') sem % mesh % edges(edge) % f % X(iY,0:5)
   end do
      
!    call checks( sem ) 

!    call sem % Integrate()

    write(STD_OUT , '(/,/,30X,A)') "\x1B[1;32m ****************** \x1B[0m"
    write(STD_OUT , '(30X,A)' ) "\x1B[1;32m Program finished! \x1B[0m"
    write(STD_OUT , '(30X,A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

end program main



