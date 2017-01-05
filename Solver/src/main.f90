program main
    use nodesAndWeights_class
    use SMConstants
    use Physics
    use MeshFileClass
    use DGSEM_Class
    use Setup_class
    use QuadMeshClass
    use QuadElementClass
    use Headers
    use ChecksModule
    use Tecplot
    implicit none
    type(MeshFile_t)       :: meshFile
    type(DGSEM_t)          :: sem
    real(kind=RP), pointer :: a(:) , b(:)
    integer                :: counter = 0
    integer                :: edge

   real(kind=RP), pointer                  :: rho(:)

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

    call ExportMeshToTecplot( sem % mesh , Setup % mesh_file ) 

    open (111 , FILE = "coords.dat" , status = "unknown" , action = "write" )
      do edge = 1 , sem % mesh % no_of_elements
         write(111, '(2F24.16)') sem % mesh % elements(edge) % X(iX:iY,:,:)

      end do
      close(111)

    call checks( sem ) 

    call ExportToTecplot( sem % mesh , './RESULTS/InitialCondition.plt')

    call sem % Integrate()

    write(STD_OUT , '(/,/,30X,A)') "\x1B[1;32m ****************** \x1B[0m"
    write(STD_OUT , '(30X,A)' ) "\x1B[1;32m Program finished! \x1B[0m"
    write(STD_OUT , '(30X,A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

end program main



