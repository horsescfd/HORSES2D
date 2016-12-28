module MeshFileClass
    use SMConstants
    type MeshFile_t
       integer                        :: no_of_nodes
       integer                        :: no_of_elements
       integer                        :: no_of_bdryedges
       real(kind=RP), allocatable     :: nodes(:)
       integer, allocatable  :: elements(:,:)
       integer, allocatable  :: faceType(:) 
       integer, allocatable  :: polynomialOrder(:)
       integer, allocatable  :: cumulativePolynomialOrder(:)
       contains
         procedure      :: read => ReadMesh
    end type MeshFile_t

    private
    public  MeshFile_t

    contains

         subroutine ReadMesh( mesh )
            use Setup_class
            use NetCDFInterface
            implicit none
            class(MeshFile_t)          :: mesh

            mesh % no_of_nodes = NetCDF_getDimension( Setup % mesh_file , "no_of_nodes" )            

            print*, "no_of_nodes: " , mesh % no_of_nodes



         end subroutine ReadMesh

!         subroutine NewMesh(mesh,K,T)
!             use Setup_class
!             implicit none
!             class(MeshFile_t)        :: mesh
!             integer               :: K
!             real(kind=RP)               :: T
!             integer               :: p , el
!             real(kind=RP)         :: h
!
!             h = 2.0_RP * T / K
! 
!             mesh % Npoints = K+1
!             mesh % Nelements = K
!
!             allocate(mesh % nodes( mesh % Npoints ) )
!             mesh % nodes = reshape((/(-T + h*p,p=0,K)/),(/K+1/))
!
!             allocate(mesh % elements( mesh % Nelements , 2 ) )
!             mesh % elements = reshape((/((el + p,p=1,0,-1),el=1,K)/),(/K , 2/), ORDER=(/2,1/))
!
!             allocate(mesh % faceType( mesh % Npoints ) )
!             mesh % faceType = FACE_INTERIOR
!             mesh % faceType(1)              = 1
!             mesh % faceType(mesh % Npoints) = 2
!
!              
!
!             allocate(mesh % polynomialOrder ( mesh % Nelements ) )
!             mesh % polynomialOrder = setup % N
!
!!       
!!            ---------------------------------------------------------
!!                   The cumulativePolynomialOrder is an array that 
!!               goes from 0 to Nelements, and such that
!!                   cumul..(0) = 0
!!                   cumul..(i) = polynomialOrder(i) + cumul...(i-1)
!!           ---------------------------------------------------------
!             allocate(mesh % cumulativePolynomialOrder( 0 : mesh % Nelements ) )
!
!            mesh % cumulativePolynomialOrder(0) = 0
!
!            do el = 1 , mesh % Nelements
!                mesh % cumulativePolynomialOrder(el) = mesh % cumulativePolynomialOrder(el-1) + mesh % polynomialOrder(el)
!            end do
!
!         end subroutine NewMesh  
!
!

end module MeshFileClass
