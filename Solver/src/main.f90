program main
    use nodesAndWeights_class
    use SMConstants
    use Physics
    use MeshFileClass
    use DGSEM_Class
    use Setup_class
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

    write(STD_OUT , * ) "Current solver: " , solver

!   =====================
!   The REAL main program
!   =====================
!

!
!   **********************************************************
!   Construct the mesh: This will be replaced by a mesh reader
!   **********************************************************
!
    call meshFile % construct( Setup % K , Setup % T )
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

    
    call checks( sem ) 

    call sem % Integrate()

    write(STD_OUT , '(/,/,A)') "\x1B[1;32m ****************** \x1B[0m"
    write(STD_OUT , '(A)' ) "\x1B[1;32m Program finished! \x1B[0m"
    write(STD_OUT , '(A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

end program main



