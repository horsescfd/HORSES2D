module MeshFileClass
    use SMConstants
    use ParamfileIO

#include "Defines.h"

    integer, parameter           :: STR_LEN_MESH = 128

    type MeshFile_t
       logical                                  :: curvilinear = .false.              ! Flag for curvilinear/not curvilinear meshes
       logical                                  :: has_pRefinement = .false.          ! Flag for p-Refinement
       integer                                  :: no_of_nodes                        ! Number of nodes in the mesh
       integer                                  :: no_of_elements                     ! Number of elements in the mesh
       integer                                  :: no_of_edges                        ! Number of edges in the mesh
       integer                                  :: no_of_bdryedges                    ! Number of edges which are boundaries
       integer                                  :: no_of_markers                      ! Number of markers
       integer                                  :: no_of_pRefinementZones             ! Number of zones in which pRefinement is performed
       integer, allocatable                     :: no_of_curvedbdryedges              ! Number of boundary edges which are curved
       integer, allocatable                     :: curves_polynomialorder             ! Curved edges polynomial order
       integer, allocatable                     :: points_of_elements(:,:)            ! Array with the points for each element ( # , element )
       integer, allocatable                     :: points_of_bdryedges(:,:)           ! Array with the points for each boundary edge ( # , edge ). Do not use it.
       integer, allocatable                     :: elements_of_edges(:,:)             ! Array with the elements which share an edge ( #el , edge ) 
       integer, allocatable                     :: polynomialOrder(:)               
       integer, allocatable                     :: cumulativePolynomialOrder(:) 
!                                                                                                           (-1 for the second point if boundary)
       integer, allocatable                     :: old_curved_bdryedges(:)                ! Which edges are curved in the mesh storage order
       integer, allocatable                     :: curved_bdryedges(:)                ! Which edges are curved
       integer, allocatable                     :: pRefinementZones(:)               ! The elements for the different pRefinement zones 
       integer, allocatable                     :: pRefinementOrder(:)               ! The elements for the different pRefinement zones 
       integer, allocatable                     :: edgeMarker(:)                      ! Array with the type of each edge ( interior, boundary, ...)
       real(kind=RP), allocatable               :: points_coords(:,:)                 ! Array with points_coordinates  (x/y , point)
       real(kind=RP), allocatable               :: curvilinear_coords(:,:,:)          ! Array with the coordinates of curvilinear edges (x/y , 0:N , edge)
       character(len=STR_LEN_MESH), allocatable :: bdryzones_names(:)
!
!      --------------------------------------------------------------------------------------------------------------------------------------
!              Intermediate arrays 
!      --------------------------------------------------------------------------------------------------------------------------------------
!
       integer, allocatable       :: points_of_edges(:,:)               ! Intermediate Array with the points for each edge ( # , edge ) 
       integer, allocatable       :: bdrymarker_of_edges(:)             ! Intermediate variable. Do not use it.
       contains
         procedure      :: Read => ReadMesh
         procedure      :: Compute => ComputeMesh
         procedure      :: Describe => DescribeMesh
         procedure      :: Destruct => MeshFile_Destruct
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
            integer                     :: curved_bdryedges
            integer                     :: marker
            integer                     :: el
            integer                     :: eID
            real(kind=RP), allocatable  :: aux(:,:)
            character(len=STR_LEN_MESH) :: name
            integer                     :: pRefZone
            integer, allocatable        :: zoneOrder
            character(len=STR_LEN_MESH) :: zoneID
!
!           ----------------------------------------------------------------------------------------------------
!                 Read nodes, elements, and boundary edges
!           ----------------------------------------------------------------------------------------------------
!
!           Dimensions
            mesh % no_of_nodes     = NetCDF_getDimension ( Setup % mesh_file , "no_of_nodes"     ) 
            mesh % no_of_elements  = NetCDF_getDimension ( Setup % mesh_file , "no_of_elements"  ) 
            mesh % no_of_bdryedges = NetCDF_getDimension ( Setup % mesh_file , "no_of_bdryedges" ) 
            mesh % no_of_markers   = NetCDF_getDimension ( Setup % mesh_file , "no_of_markers"   ) 

!           Allocate variables
            allocate ( mesh % points_of_elements  ( POINTS_PER_QUAD , mesh % no_of_elements  )  ) 
            allocate ( mesh % points_coords       ( NDIM , mesh % no_of_nodes                )  ) 
            allocate ( mesh % points_of_bdryedges ( POINTS_PER_EDGE , mesh % no_of_bdryedges )  ) 
            allocate ( mesh % bdrymarker_of_edges ( mesh % no_of_bdryedges                   )  ) 
            allocate ( mesh % bdryzones_names     ( mesh % no_of_markers                     )  )
            allocate ( mesh % polynomialOrder     ( mesh % no_of_elements                    )  )

!           Gather variables
            call NetCDF_getVariable ( Setup % mesh_file , "points_of_quads"     , mesh % points_of_elements  ) 
            call NetCDF_getVariable ( Setup % mesh_file , "points"              , mesh % points_coords       ) 
            call NetCDF_getVariable ( Setup % mesh_file , "points_of_bdryedges" , mesh % points_of_bdryedges ) 
            call NetCDF_getVariable ( Setup % mesh_file , "bdrymarker_of_edges" , mesh % bdrymarker_of_edges ) 
!
!           Assign the default polynomial order until replaced by the pRefinement polynomial order if proceeds
!           --------------------------------------------------------------------------------------------------
            mesh % polynomialOrder  = Setup % N

            mesh % points_coords = mesh % points_coords / RefValues % L       ! Scale with the reference length

            do marker = 1 , mesh % no_of_markers
               write(name , '(A,I0)') "marker" , marker
               call NetCDF_getVariable ( Setup % mesh_file , trim(name) , mesh % bdryzones_names(marker) )
            end do
!
!           Gather pRefinement zones
!           ------------------------
            mesh % no_of_pRefinementZones = NetCDF_getDimension( Setup % mesh_file , "no_of_pRefinementzones" )
!
!           *************************************************
!                 ----> The mesh has p-Refinement <----
            if ( mesh % no_of_pRefinementZones .ne. -1 ) then
!           *************************************************
!
!              ===============================
               mesh % has_pRefinement = .true.
!              ===============================
!
               allocate ( mesh % pRefinementZones ( mesh % no_of_elements ))
               allocate ( mesh % pRefinementOrder ( mesh % no_of_pRefinementZones ) )
               call NetCDF_getVariable ( Setup % mesh_file , "pRefinement_zones" , mesh % pRefinementZones )
!
!              Read from case file the polynomial order
!              ----------------------------------------
               do pRefZone = 1 , mesh % no_of_pRefinementZones 
                  write( zoneID , '(I0)') pRefZone
                  call ReadValueInRegion( trim(Setup % case_file) , trim(zoneID) , zoneOrder , "# define p-Refinement" , "# end" )
                  
                  if ( allocated ( zoneOrder ) ) then
                     mesh % pRefinementOrder(pRefZone) = zoneOrder
                     deallocate( zoneOrder )
                  end if
               end do
!
!              Assign the polynomial order to the elements
!              -------------------------------------------
               do eID = 1 , mesh % no_of_elements
                  mesh % polynomialOrder(eID) = mesh % pRefinementOrder ( mesh % pRefinementZones ( eID ) )  
               end do

            end if
            

!           Gather curved boundaries
            curved_bdryedges        = NetCDF_getDimension( Setup % mesh_file , "no_of_curvilinearedges" )
!
!           **********************************
!             ----> Curved boundaries <----
            if (curved_bdryedges .ne. -1) then
!           **********************************
!
!              ===========================
               mesh % curvilinear = .true.
!              ===========================
!
               allocate ( mesh % no_of_curvedbdryedges  ) 
               allocate ( mesh % curves_polynomialorder ) 

               mesh % no_of_curvedbdryedges  = curved_bdryedges
               mesh % curves_polynomialorder = NetCDF_getDimension( Setup % mesh_file , "Np1" ) - 1
!
!              --------------------------------------
!                 Get which edges are curved
!              --------------------------------------
!
               allocate ( mesh % old_curved_bdryedges ( mesh % no_of_curvedbdryedges ) ) 
               allocate ( mesh % curved_bdryedges ( mesh % no_of_curvedbdryedges ) ) 
               call NetCDF_getVariable ( Setup % mesh_file , "curvilinear_edges" , mesh % old_curved_bdryedges ) 
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
               mesh % curvilinear_coords(1,:,:) = aux / RefValues % L
!
!              Obtain y-coordinates
!              --------------------
               call NetCDF_getVariable( Setup % mesh_file , "y_curvilinear_edges" , aux )
               mesh % curvilinear_coords(2,:,:) = aux / RefValues % L
!
!              ------------------------------------
!                 Free the auxiliary variable
!              ------------------------------------
!
               deallocate( aux )
            
            end if
!       
!            ---------------------------------------------------------
!                   The cumulativePolynomialOrder is an array that 
!               goes from 0 to Nelements, and such that
!                   cumul..(0) = 0
!                   cumul..(i) = (polynomialOrder(i)+1)**2.0 + cumul...(i-1)
!           ---------------------------------------------------------
             allocate(mesh % cumulativePolynomialOrder( 0 : mesh % no_of_elements ) )

            mesh % cumulativePolynomialOrder(0) = 0
      
            do el = 1 , mesh % no_of_elements
               mesh % cumulativePolynomialOrder(el)   = mesh % cumulativePolynomialOrder(el-1) + (mesh % polynomialOrder(el) + 1) * (mesh % polynomialOrder(el) + 1)
            end do


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
            character(len=STR_LEN_MESH)   :: auxstr
            integer                    :: zone

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


            write(STD_OUT,'(/)')
            write(auxstr , '(I0,A)') mesh % no_of_markers , " boundary zones found"
            call SubSection_Header( trim(auxstr) )

            do zone = 1 , mesh % no_of_markers
               write(STD_OUT , '(30X,A,A,I0,A,A20)') "-> ", "Zone ",zone,": ", trim(mesh % bdryzones_names(zone))
            end do

         end subroutine DescribeMesh

         subroutine ComputeMesh( self )
            implicit none
            class(MeshFile_t)          :: self
!           -----------------------------------------
            integer                    :: eID
!
!           -------------------------------
!              Compute number of edges
!           -------------------------------
!
            if (mod(4*self % no_of_elements + self % no_of_bdryedges , 2) .eq. 0) then
               self % no_of_edges = (4*self % no_of_elements + self % no_of_bdryedges) / 2 
               allocate( self % points_of_edges ( POINTS_PER_EDGE , self % no_of_edges ) ) 
               allocate( self % elements_of_edges( POINTS_PER_EDGE , self % no_of_edges ) )
            else
               print*, "The mesh is not consistent."
               stop "Stopped."
            end if

            self % points_of_edges = -1
            self % elements_of_edges = -1

            call computeFaces( self )
            call computeFaceMarkers( self )

         end subroutine ComputeMesh

         subroutine computeFaces( mesh )
            use Headers
            implicit none
            class(MeshFile_t)          :: mesh
!           ----------------------------------------
            integer                    :: eID
            integer                    :: elFace
            integer                    :: currentFace
            integer                    :: previousFaces
            logical                    :: exists
            integer, dimension(2)      :: face
!
!           -------------------------------
!              Obtain faces      
!           -------------------------------
!
            currentFace = 1
            do eID = 1 , mesh % no_of_elements
               do elFace = 1 , POINTS_PER_QUAD
!
!                 ----------------------------------------
!                    Select face
!                 ----------------------------------------
!
                  if (elFace .ne. POINTS_PER_QUAD) then
                     face = mesh % points_of_elements([elFace , elFace + 1] , eID )
                  else
                     face = mesh % points_of_elements([elFace , 1] , eID ) 
                  end if
!
!                 ---------------------------------------
!                    Check if the face exists already
!                 ---------------------------------------
!
                  exists = .false.
                  do previousFaces = 1 , currentFace - 1             
                     if ( facesEqual( mesh % points_of_edges(:,previousFaces) , face )) then    ! There is an existing face
                        exists = .true.
                        mesh % elements_of_edges(2 , previousFaces) = eID
                        
                        exit
                     end if      ! Continue 
                  end do
!
!                 -------------------------------------------
!                    If does not exist, store it
!                 -------------------------------------------
!
                  if (.not. exists) then
                     mesh % points_of_edges( : , currentFace ) = face
                     mesh % elements_of_edges(1 , currentFace ) = eID

!                    Move to next face
!                    -----------------
                     currentFace = currentFace + 1
                  end if
               end do
            end do
         
         end subroutine computeFaces

         subroutine computeFaceMarkers( mesh )
            implicit none
            class(MeshFile_t)          :: mesh
            integer                    :: edge
            integer                    :: bdryface
!           ------------------------------------------------
            
!
!           --------------------------
!           Allocate face marker array
!           --------------------------
!
            allocate ( mesh % edgeMarker ( mesh % no_of_edges ) )

            do edge = 1 , mesh % no_of_edges

               if (mesh % elements_of_edges(2 , edge) .eq. -1) then   ! Is a boundary face
!
!                 ------------------------------------------------------------------------
!                    Search for a similar face in boundary faces structure
!                 ------------------------------------------------------------------------
!
                  do bdryface = 1 , mesh % no_of_bdryedges

                     if (facesEqual(mesh % points_of_bdryedges(:,bdryface) , mesh % points_of_edges(:,edge) ) ) then    ! Face found
!
!                       ----------------------------------------------
!
                        mesh % edgeMarker (edge) =  mesh % bdrymarker_of_edges(bdryface)        
!        
!                       --------------------------------------------------------------------
!                       Change the edge in curved_bdryedges to the new numeration
! 
                        if (mesh % curvilinear) then
                           call changeEntry( oldarray = mesh % old_curved_bdryedges , newarray = mesh % curved_bdryedges, old = bdryface , new = edge )
                        end if
!                       --------------------------------------------------------------------

                        exit        ! Exit from bdryface loop if the face is found

                     end if

                  end do

               else     ! Otherwise, an interior face is considered

                  mesh % edgeMarker (edge) = FACE_INTERIOR

               end if

            end do


         end subroutine computeFaceMarkers

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

         subroutine changeEntry( oldarray , newarray , old , new ) 
            implicit none
            integer,    intent(inout)        :: oldarray(:)
            integer,    intent(inout)        :: newarray(:)
            integer,    intent(in)           :: old
            integer,    intent(in)           :: new
!           -----------------------------------------------
            integer                          :: pos

            do pos = 1 , size(oldarray)

               if (oldarray(pos) .eq. old) then
                  newarray(pos) = new
                  return
               end if

            end do


         end subroutine

         subroutine MeshFile_Destruct ( self ) 
            implicit none
            class(MeshFile_t)          :: self
   
            self % no_of_nodes = 0
            self % no_of_elements = 0
            self % no_of_edges = 0
            self % no_of_bdryedges = 0
            self % no_of_markers = 0
            self % no_of_curvedBdryEdges = 0
            self % curves_polynomialorder = 0
            deallocate ( self % points_of_elements        ) 
            deallocate ( self % points_of_bdryedges       ) 
            deallocate ( self % elements_of_edges         ) 
            deallocate ( self % polynomialOrder           ) 
            deallocate ( self % cumulativePolynomialOrder ) 
            deallocate ( self % edgeMarker                ) 
            deallocate ( self % points_coords             ) 
            deallocate ( self % bdryzones_names           ) 
            deallocate ( self % points_of_edges           ) 
            deallocate ( self % bdrymarker_of_edges       ) 
            if ( self % curvilinear ) then
               deallocate ( self % old_curved_bdryedges ) 
               deallocate ( self % curved_bdryedges     ) 
               deallocate ( self % curvilinear_coords   ) 
            end if


         end subroutine MeshFile_Destruct

end module MeshFileClass
