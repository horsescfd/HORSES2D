program main
    use SMConstants
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
    integer                :: eID , edID
   real(kind=RP)  :: tstart , tend

!   =====================
!   The REAL main program
!   =====================
!
    call Main_Header("High-order discontinuous Galerkin CFD 2D Solver")

    call setup % Initialization()

    call InitializePhysics

!
!   **********************************************************
!   Read the mesh TODO: Move into DGSEM
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
    call ExportMeshToTecplot( sem % mesh , Setup % mesh_file )       ! This one should be inside DGSEM or Mesh
!   
!   ***********************************************
!   Set the initial condition to all flow variables
!   ***********************************************
!
    call sem % SetInitialCondition()
    call ExportToTecplot( sem % mesh , './RESULTS/InitialCondition.plt')      ! This one should be inside DGSEM or Mesh

    call checks( sem )           ! This one should be inside DGSEM




call cpu_time(tstart)
    call sem % Integrate()
call cpu_time(tend)
   print*, "Time: ",tend-tstart

     
    write(STD_OUT , '(/,/,30X,A)') "\x1B[1;32m ****************** \x1B[0m"
    write(STD_OUT , '(30X,A)' ) "\x1B[1;32m Program finished! \x1B[0m"
    write(STD_OUT , '(30X,A,/,/)') "\x1B[1;32m ****************** \x1B[0m"

end program main



