module MeshFileClass
    use SMConstants


    integer, parameter           :: POINTS_PER_EDGE = 2
    integer, parameter           :: POINTS_PER_QUAD = 4

    type MeshFile_t
       integer                    :: no_of_nodes
       integer                    :: no_of_elements
       integer                    :: no_of_edges
       integer                    :: no_of_bdryedges
       integer, allocatable       :: no_of_curvedbdryedges
       integer, allocatable       :: curves_polynomialorder
       logical                    :: curvilinear = .false.
       real(kind=RP), allocatable :: points_coords(:,:)
       integer, allocatable       :: points_of_elements(:,:)
       integer, allocatable       :: points_of_edges(:,:)
       integer, allocatable       :: points_of_bdryedges(:,:)
       integer, allocatable       :: edges_of_elements(:,:)
       integer, allocatable       :: elements_of_edges(:,:)
       integer, allocatable       :: faceType(:)
       integer, allocatable       :: polynomialOrder(:)
       integer, allocatable       :: cumulativePolynomialOrder(:)
       real(Kind=RP), allocatable :: curvilinear_coords(:,:,:)
       contains
         procedure      :: Read => ReadMesh
         procedure      :: Compute => ComputeMesh
         procedure      :: Describe => DescribeMesh
    end type MeshFile_t

    private
    public  MeshFile_t

    contains

         subroutine ReadMesh( mesh )
            use Setup_class
            use Physics
            use NetCDFInterface
            implicit none
            class(MeshFile_t)          :: mesh
!           -----------------------------------------------------
            integer                    :: curved_bdryedges
            real(kind=RP), allocatable :: aux(:,:)
!
!           ----------------------------------------------------------------------------------------------------
!                 Read nodes, elements, and boundary edges
!           ----------------------------------------------------------------------------------------------------
!
!           Dimensions
            mesh % no_of_nodes     = NetCDF_getDimension( Setup % mesh_file , "no_of_nodes" )
            mesh % no_of_elements  = NetCDF_getDimension( Setup % mesh_file , "no_of_elements" )
            mesh % no_of_bdryedges = NetCDF_getDimension( Setup % mesh_file , "no_of_bdryedges" )
   
!           Allocate variables
            allocate( mesh % points_of_elements( POINTS_PER_QUAD , mesh % no_of_elements ) )
            allocate( mesh % points_coords( NDIM , mesh % no_of_nodes ) )
            allocate( mesh % points_of_bdryedges ( POINTS_PER_EDGE , mesh % no_of_bdryedges ) )

!           Gather variables
            call NetCDF_getVariable( Setup % mesh_file , "points_of_quads" , mesh % points_of_elements )
            call NetCDF_getVariable( Setup % mesh_file , "points" , mesh % points_coords)
            call NetCDF_getVariable( Setup % mesh_file , "points_of_bdryedges" , mesh % points_of_bdryedges )

!           Gather curved boundaries
            curved_bdryedges        = NetCDF_getDimension( Setup % mesh_file , "no_of_curvilinearedges" )
!
!           **********************************
!                 Curved boundaries
            if (curved_bdryedges .ne. -1) then
!           **********************************
!
!              ===========================
               mesh % curvilinear = .true.
!              ===========================
!
               allocate( mesh % no_of_curvedbdryedges )
               allocate( mesh % curves_polynomialorder )

               mesh % no_of_curvedbdryedges = curved_bdryedges
               mesh % curves_polynomialorder = NetCDF_getDimension( Setup % mesh_file , "Np1" ) - 1
!
!              -----------------------------------------------------------------------
!                 Get curved patches from file: An auxiliary variable "aux" is needed 
!              -----------------------------------------------------------------------
!
               allocate ( mesh % curvilinear_coords ( NDIM , mesh % curves_polynomialorder + 1 , mesh % no_of_curvedbdryedges )  ) 
               allocate ( aux                       (        mesh % curves_polynomialorder + 1 , mesh % no_of_curvedbdryedges )  ) 
!
!              Obtain x-coordinates
!              --------------------
               call NetCDF_getVariable( Setup % mesh_file , "x_curvilinear_edges" , aux )
               mesh % curvilinear_coords(1,:,:) = aux
!
!              Obtain y-coordinates
!              --------------------
               call NetCDF_getVariable( Setup % mesh_file , "y_curvilinear_edges" , aux )
               mesh % curvilinear_coords(2,:,:) = aux
!
!              ------------------------------------
!                 Free the auxiliary variable
!              ------------------------------------
!
               deallocate( aux )
            
            end if
!
!           ***************************
!              Compute the mesh       
!           ***************************
!
            call mesh % Compute
            call mesh % Describe

         end subroutine ReadMesh

         subroutine DescribeMesh( mesh )
            use Setup_class
            use Headers
            implicit none
            class(MeshFile_t)          :: mesh

            write(STD_OUT,'(/)')
            call Section_Header("Reading mesh")
            write(STD_OUT,'(/)')

            call SubSection_Header('Mesh file "' // trim(Setup % mesh_file) //'"')
            write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Number of nodes: ", mesh % no_of_nodes ,"."
            write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Number of elements: ", mesh % no_of_elements ,"."
            write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Number of edges: ", mesh % no_of_edges ,"."
            write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Number of boundary edges: ", mesh % no_of_bdryedges ,"."

            if (mesh % curvilinear) then

               write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Number of curved edges: ", mesh % no_of_curvedbdryedges ,"."
               write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Curved edges polynomial order: ", mesh % curves_polynomialorder ,"."
               
            end if

         end subroutine DescribeMesh

         subroutine ComputeMesh( self )
            implicit none
            class(MeshFile_t)          :: self
!           ----------------------------------------
            integer                    :: eID
            integer                    :: currentFace
            integer                    :: previousFaces

!
!           -------------------------------
!              Compute number of edges
!           -------------------------------
!
            if (mod(self % no_of_elements + self % no_of_bdryedges , 2) .eq. 0) then
               self % no_of_edges = (self % no_of_elements + self % no_of_bdryedges) / 2 
               allocate( self % points_of_edges ( POINTS_PER_EDGE , self % no_of_edges ) ) 
               allocate( self % edges_of_elements( POINTS_PER_QUAD , self % no_of_edges ) )
               allocate( self % elements_of_edges( POINTS_PER_EDGE , self % no_of_elements ) )
            else
               print*, "The mesh is not consistent."
               stop "Stopped."
            end if
!
!           -------------------------------
!              Obtain faces      
!           -------------------------------
!
            currentFace = 0
            do eID = 1 , self % no_of_elements
             
            end do

         end subroutine ComputeMesh
!
!        ****************************************************************************************
!           Auxiliar subroutines
!        ****************************************************************************************
!
         function facesEqual(face1 , face2) result(val)
            implicit none
            integer, dimension(2), intent(in)      :: face1
            integer, dimension(2), intent(in)      :: face2
            logical                                :: val

            if (((face1(1) .eq. face2(1)) .and. (face1(2) .eq. face2(2))) .or. (((face1(2) .eq. face2(1)) .and. (face1(1) .eq. face2(2))))) then
               val = .true.      ! They are equal
            else
               val = .false.     ! They are not equal
            end if

         end function facesEqual

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
