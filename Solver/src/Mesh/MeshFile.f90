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
module MeshFileClass
    use SMConstants
    use ParamfileIO
    use IntegerArrayLinkedList

#include "Defines.h"

    integer, parameter           :: STR_LEN_MESH = 128
                                                                                      ! 
    type MeshFile_t                                                                   ! -------------------------------------------------------------------------
       logical                                  :: curvilinear = .false.              ! Flag for curvilinear/not curvilinear meshes
       logical                                  :: has_pRefinement = .false.          ! Flag for p-Refinement
       integer                                  :: no_of_nodes                        ! Number of nodes in the mesh
       integer                                  :: no_of_elements                     ! Number of elements in the mesh
       integer, pointer                         :: no_of_edges                        ! Number of edges in the mesh
       integer                                  :: no_of_bdryedges                    ! Number of edges which are boundaries
       integer                                  :: no_of_markers                      ! Number of markers
       integer                                  :: no_of_pRefinementZones             ! Number of zones in which pRefinement is performed
       integer, allocatable                     :: no_of_curvededges                  ! Number of edges which are curved
       integer, allocatable                     :: curves_polynomialorder             ! Curved edges polynomial order
       integer, allocatable                     :: points_of_elements(:,:)            ! Array with the points for each element ( # , element )
       type(IntegerArrayLinkedList_t)           :: points_of_edges                    ! Linked list with the points for each edge
       type(IntegerArrayLinkedList_t)           :: elements_of_edges                  ! Linked list with the elements which share an edge 
       integer, allocatable                     :: polynomialOrder(:)                 !
       integer, allocatable                     :: cumulativePolynomialOrder(:)       !
       integer, allocatable                     :: curved_edges_points(:,:)           ! Points of the edges which are curved
       integer, allocatable                     :: curved_edges(:)                    ! Which edges are curved
       integer, allocatable                     :: pRefinementZones(:)                ! The elements for the different pRefinement zones 
       integer, allocatable                     :: pRefinementOrder(:)                ! The elements for the different pRefinement zones 
       integer(kind=1), allocatable             :: edgeMarker(:)                      ! Array with the type of each edge ( interior, boundary, ...)
       real(kind=RP), allocatable               :: points_coords(:,:)                 ! Array with points_coordinates  (x/y , point)
       real(kind=RP), allocatable               :: curvilinear_coords(:,:,:)          ! Array with the coordinates of curvilinear edges (x/y , 0:N , edge)
       character(len=STR_LEN_MESH), allocatable :: bdryzones_names(:)                 ! 
!                                                                                     !
!      Intermediate arrays                                                            !
!      -------------------                                                            !
       integer, allocatable :: points_of_bdryedges(:,:)                               ! Array with the points for each boundary edge ( # , edge ). Do not use it.
       integer, allocatable :: bdrymarker_of_edges(:)                                 ! Intermediate variable. Do not use it.
       contains                                                                       !
         procedure      :: Read => ReadMesh                                           !
         procedure      :: Describe => DescribeMesh                                   !
         procedure      :: Destruct => MeshFile_Destruct                              !
    end type MeshFile_t                                                               ! --------------------------------------------------------------------------
                                                                                      !
    private
    public  MeshFile_t

    contains
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              READ FILE SUBROUTINES
!              ---------------------
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
         subroutine ReadMesh( mesh )
            use Setup_class
            use Physics
            use NetCDFInterface
            implicit none
            class(MeshFile_t)          :: mesh
!           -----------------------------------------------------
            integer                     :: no_of_curved_edges
            integer                     :: marker
            integer                     :: el
            integer                     :: eID
            real(kind=RP), allocatable  :: aux(:,:)
            character(len=STR_LEN_MESH) :: name
            integer                     :: pRefZone
            integer, allocatable        :: zoneOrder
            character(len=STR_LEN_MESH) :: zoneID
!$          integer                     :: omp_get_thread_num
!
!           ****************************************
!           Read nodes, elements, and boundary edges
!           ****************************************
!
!           Dimensions: Directly obtained from the NetCDF mesh file
!           -------------------------------------------------------
            mesh % no_of_nodes     = NetCDF_getDimension ( Setup % mesh_file , "no_of_nodes"     ) 
            mesh % no_of_elements  = NetCDF_getDimension ( Setup % mesh_file , "no_of_elements"  ) 
            mesh % no_of_bdryedges = NetCDF_getDimension ( Setup % mesh_file , "no_of_bdryedges" ) 
            mesh % no_of_markers   = NetCDF_getDimension ( Setup % mesh_file , "no_of_markers"   ) 
!
!           Allocate variables
!           ------------------
            allocate ( mesh % points_of_elements  ( POINTS_PER_QUAD , mesh % no_of_elements  )  ) 
            allocate ( mesh % points_coords       ( NDIM , mesh % no_of_nodes                )  ) 
            allocate ( mesh % points_of_bdryedges ( POINTS_PER_EDGE , mesh % no_of_bdryedges )  ) 
            allocate ( mesh % bdrymarker_of_edges ( mesh % no_of_bdryedges                   )  ) 
            allocate ( mesh % bdryzones_names     ( mesh % no_of_markers                     )  )
            allocate ( mesh % polynomialOrder     ( mesh % no_of_elements                    )  )
!
!           Gather variables from the NetCDF mesh file
!           ------------------------------------------
!$omp parallel num_threads(4)
!$ if ( omp_get_thread_num() .eq. 0 ) then
            call NetCDF_getVariable ( Setup % mesh_file , "points_of_quads"     , mesh % points_of_elements  ) 
!$ elseif ( omp_get_thread_num() .eq. 1 ) then
            call NetCDF_getVariable ( Setup % mesh_file , "points"              , mesh % points_coords       ) 
!$ elseif ( omp_get_thread_num() .eq. 2 ) then
            call NetCDF_getVariable ( Setup % mesh_file , "points_of_bdryedges" , mesh % points_of_bdryedges ) 
!$ elseif ( omp_get_thread_num() .eq. 3 ) then
            call NetCDF_getVariable ( Setup % mesh_file , "bdrymarker_of_edges" , mesh % bdrymarker_of_edges ) 
!$ end if 
!$omp barrier
!$omp end parallel
!
!           Assign the default polynomial order until replaced by the pRefinement polynomial order if proceeds
!           --------------------------------------------------------------------------------------------------
            mesh % polynomialOrder  = Setup % N
!
!           Fill each point coordinates (scaled with the reference length)
!           --------------------------------------------------------------
            mesh % points_coords = mesh % points_coords / RefValues % L

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
               allocate ( mesh % pRefinementOrder ( 0 : mesh % no_of_pRefinementZones ) )
               call NetCDF_getVariable ( Setup % mesh_file , "pRefinement_zones" , mesh % pRefinementZones )
!
!              Read from case file the polynomial order
!              ----------------------------------------
               mesh % pRefinementOrder(0) = Setup % N          ! Assign the default value to the "0" zones

               do pRefZone = 1 , mesh % no_of_pRefinementZones 
                  write( zoneID , '(I0)') pRefZone
                  call ReadValueInRegion( trim(Setup % case_file) , trim(zoneID) , zoneOrder , "# define p-Refinement" , "# end" )
                  
                  if ( allocated ( zoneOrder ) ) then
                     mesh % pRefinementOrder(pRefZone) = zoneOrder
                     deallocate( zoneOrder )
                  else
                     mesh % pRefinementOrder(pRefZone) = Setup % N
                  end if
               end do
!
!              Assign the polynomial order to the elements
!              -------------------------------------------
               do eID = 1 , mesh % no_of_elements
                  mesh % polynomialOrder(eID) = mesh % pRefinementOrder ( mesh % pRefinementZones ( eID ) )  
               end do

            end if
!
!           Gather curved boundaries
!           ------------------------
            no_of_curved_edges = NetCDF_getDimension( Setup % mesh_file , "no_of_curvilinearedges" )
!
!           **********************************
!             ----> Curved edges <----
            if (no_of_curved_edges .ne. -1) then
!           **********************************
!
!              ===========================
               mesh % curvilinear = .true.
!              ===========================
!
!              no_of_curvedEdges: How many edges are curved
!              --------------------------------------------
               allocate ( mesh % no_of_curvedEdges  ) 
               mesh % no_of_curvedEdges  = no_of_curved_edges
!
!              curves_polynomialorder: The polynomial order used for the edges
!              ---------------------------------------------------------------
               allocate ( mesh % curves_polynomialorder ) 
               mesh % curves_polynomialorder = NetCDF_getDimension( Setup % mesh_file , "Np1" ) - 1
!
!              curved_edges_points: The two nodes of each curved edge
!              ------------------------------------------------------
               allocate ( mesh % curved_edges_points ( NDIM , mesh % no_of_curvedEdges ) ) 
               call NetCDF_getVariable ( Setup % mesh_file , "curvilinear_edges" , mesh % curved_edges_points ) 

!
!              Get curved patches from file: An auxiliary variable "aux" is needed to store each coordinate
!              --------------------------------------------------------------------------------------------
               allocate ( mesh % curvilinear_coords ( NDIM , mesh % curves_polynomialorder + 1 , mesh % no_of_curvedEdges )  ) 
               allocate ( aux                       (        mesh % curves_polynomialorder + 1 , mesh % no_of_curvedEdges )  ) 
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
!              Free the auxiliary variable
!              ---------------------------
               deallocate( aux )
            
            end if
!       
!            -------------------------------------------------------------------------
!                   The cumulativePolynomialOrder is an array that 
!               goes from 0 to Nelements, and such that
!                   cumul..(0) = 0
!                   cumul..(i) = (polynomialOrder(i)+1)**2.0 + cumul...(i-1)
!           --------------------------------------------------------------------------
             allocate(mesh % cumulativePolynomialOrder( 0 : mesh % no_of_elements ) )

            mesh % cumulativePolynomialOrder(0) = 0
      
            do el = 1 , mesh % no_of_elements
               mesh % cumulativePolynomialOrder(el)   = mesh % cumulativePolynomialOrder(el-1) + (mesh % polynomialOrder(el) + 1) * (mesh % polynomialOrder(el) + 1)
            end do
!
!           ****************
!           Compute the mesh       
!           ****************
!
!           First, compute all the faces
!           ----------------------------
            call computeFaces( mesh )
!
!           Second, merge all divided faces
!           -------------------------------
            call mergeDividedFaces ( mesh ) 
            call mergeDividedFaces ( mesh ) 
!
!           Third, compute which elements belong to each edge
!           -------------------------------------------------
            call computeElementOfEdges ( mesh ) 
!
!           Fourth, assign each face a marker
!           ---------------------------------
            call computeFaceMarkers( mesh )
!
!           Fifth, if curvilinear, assign each edge a curve
!           -----------------------------------------------
            if ( mesh % curvilinear ) then
               call computeCurvedEdges( mesh ) 
            end if
!
!           Describe the mesh
!           -----------------
            call mesh % Describe

         end subroutine ReadMesh
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!                 COMPUTE SUBROUTINES
!                 -------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
         subroutine computeFaces( mesh )
            use Headers
            use Sorting
            implicit none
            class(MeshFile_t)          :: mesh
!
!           ---------------
!           Local variables
!           ---------------
!
            integer                             :: eID , edID
            integer                             :: elFace
            integer, dimension(2)               :: face
            integer                             :: pos ( mesh % no_of_bdryedges ) 
            integer                             :: counter
            class(IntegerArrayEntry_t), pointer :: current
            integer, parameter                  :: elFaceIDs(EDGES_PER_QUAD + 1) = [1,2,3,4,1]
!
!           Initialize the structure
!           ------------------------
            mesh % points_of_edges = ConstructIntegerArrayLinkedList()
!
!           no_of_edges points to the size of the linked list
!           -------------------------------------------------
            mesh % no_of_edges => mesh % points_of_edges % no_of_entries
!
!           Obtain the four faces of each element 
!           -------------------------------------
            do eID = 1 , mesh % no_of_elements
               do elFace = 1 , POINTS_PER_QUAD
!
!                 Select face
!                 -----------
                  face = mesh % points_of_elements( elFaceIDs( elFace : elFace + 1 ) , eID ) 
!
!                 Add it to points_of_edges
!                 -------------------------
                  call mesh % points_of_edges % Add ( POINTS_PER_EDGE , face )

               end do
            end do
!
!           Set which are boundaries in the Linked List attribute field (true/false)
!           ------------------------------------------------------------------------
            do edID = 1 , mesh % no_of_bdryedges
!
!              Look for the bdryedge points in the edges data
!              ----------------------------------------------
               pos(edID) = mesh % points_of_edges % Search( POINTS_PER_EDGE , mesh % points_of_bdryedges(:,edID) )

            end do
!
!           Set the value to the linked list attribute
!           ------------------------------------------
            current => mesh % points_of_edges % head 
            counter = 1 
!
!           Sort the positions to loop just once the list
!           ---------------------------------------------
            call Qsort(pos)         

            do edID = 1 , mesh % points_of_edges % no_of_entries

               if ( pos(counter) .eq. edID ) then
                  current % attribute = .true.
                  counter = counter + 1 

                  if ( counter .eq. mesh % no_of_bdryedges + 1) exit 

               end if

               current => current % next

            end do

         end subroutine computeFaces

         subroutine mergeDividedFaces( mesh ) 
!
!           ******************************************************************************************
!
!                 This subroutine looks for subdivided faces in the points_of_edges data
!              Basically, the process is:
!
!                 1/ Search for a point with just 3 faces (X in the sketch).
!
!                 2/ Check that this face is not a boundary face. 
!                                         (otherwise it is valid to have 3 edges)
!                 3/ Get the three points which are connected through its faces (O in the sketch).
!
!                 4/ Get the four points which each of the previous points are connected.
!
!                 5/ If one of those is of the three (O) points, they belong to a subdivided face.
!
!                          BEFORE                                     AFTER
!
!                 +-----------O----------+                   +----------2----------+
!                 |          /|          |                   |          |          |
!                 |         /=|          |           \       |          |          |
!                 |        /==|          |    ========\      |          |          |
!                 O-------X===|          |    =========\     +----------3          |
!                 |        \==|          |    =========/     |          |          |
!                 |         \=|          |    ========/      |          |          |
!                 |          \|          |           /       |          |          |
!                 +-----------O----------+                   +----------1----------+
!
!           ******************************************************************************************
!
            implicit none
            class(MeshFile_t)          :: mesh
!
!           ---------------
!           Local variables
!           ---------------
!
            integer                    :: edges(POINTS_PER_QUAD)
            integer                    :: neigh_edges(2*POINTS_PER_QUAD)
            logical                    :: atts(POINTS_PER_QUAD)
            integer                    :: points(POINTS_PER_SUBDIVIDED_EDGE)
            integer                    :: point
            integer                    :: newSubdividedEdge(POINTS_PER_SUBDIVIDED_EDGE)
            integer                    :: points_of_edge( POINTS_PER_EDGE )
            integer                    :: pos
            integer                    :: threeEdgesNode , j  , k
            integer                    :: oneDArray(1)
            integer                    :: edgeToRemove(2)
!
!           This subroutine removes the divided edges information
!           -----------------------------------------------------
!$omp parallel do default(private) shared(mesh)
mainloop:   do threeEdgesNode = 1 , mesh % no_of_nodes
               oneDArray = threeEdgesNode
               edges = mesh % points_of_edges % SearchIfPresent( ONE , oneDArray , POINTS_PER_QUAD , atts) 

               if ( edges     ( POINTS_PER_QUAD ) .ne. -1 )  cycle  ! It has four associated edges      -> Is not subdivided
               if ( edges     ( THREE           ) .eq. -1 )  cycle  ! It has already been merged
               if (       any ( atts            )         )  cycle  ! Not connected to boundary faces   -> Is not subdivided
!
!              Gather the three neighbouring points
!              ------------------------------------
               do j = 1 , POINTS_PER_SUBDIVIDED_EDGE
                  points_of_edge = mesh % points_of_edges % Get ( POINTS_PER_EDGE , edges(j) )
         
                  if ( points_of_edge(ONE) .ne. threeEdgesNode ) then
                     points(j) = points_of_edge(ONE)
               
                  elseif ( points_of_edge(TWO) .ne. threeEdgesNode ) then
                     points(j) = points_of_edge(TWO)
   
                  end if

               end do
!
!              Gather all the edges that emerge from each point
!              ------------------------------------------------
               do j = 1 , POINTS_PER_SUBDIVIDED_EDGE
                  oneDArray = points(j)
                  neigh_edges = mesh % points_of_edges % SearchIfPresent( ONE , oneDArray , 2 * POINTS_PER_QUAD ) 

                  do k = 1 , 2 * POINTS_PER_QUAD
                     if ( neigh_edges(k) .eq. -1 ) cycle 

                     points_of_edge = mesh % points_of_edges % Get ( POINTS_PER_EDGE , neigh_edges(k) )

                     if ( (points_of_edge(ONE) .ne. points(j)) .and. (points_of_edge(ONE) .ne. threeEdgesNode) ) then
                        point = points_of_edge(ONE)

                     elseif ( (points_of_edge(TWO) .ne. points(j)) .and. (points_of_edge(TWO) .ne. threeEdgesNode) ) then
                        point = points_of_edge(TWO)

                     else
                        point = -1

                     end if

                     if ( any( points .eq. point ) ) then       
!
!                       This is the found edge: It is arranged such that the middle node is the last (important!)
!                       ----------------------
                        newSubdividedEdge = [ points_of_edge , threeEdgesNode ] 
!
!                       Remove the three edges
!                       ----------------------
                        edgeToRemove = [newSubdividedEdge(ONE) , newSubdividedEdge(TWO)]
                        pos = mesh % points_of_edges % Search ( POINTS_PER_EDGE , edgeToRemove )
                        call  mesh % points_of_edges % Remove ( pos )

                        edgeToRemove = [newSubdividedEdge(ONE) , newSubdividedEdge(THREE)]
                        pos = mesh % points_of_edges % Search ( POINTS_PER_EDGE , edgeToRemove )
                        call  mesh % points_of_edges % Remove ( pos )

                        edgeToRemove = [newSubdividedEdge(TWO) , newSubdividedEdge(THREE)]
                        pos = mesh % points_of_edges % Search ( POINTS_PER_EDGE , edgeToRemove )
                        call  mesh % points_of_edges % Remove ( pos )
!
!                       Add the new subdivided edge
!                       ---------------------------
                        call mesh % points_of_edges % Add ( POINTS_PER_SUBDIVIDED_EDGE , newSubdividedEdge ) 
!
!                       Continue
!                       --------
                        cycle mainloop

                     end if
                  end do
               end do

            end do   mainloop
!$omp end parallel do

         end subroutine mergeDividedFaces

         subroutine computeElementOfEdges ( mesh ) 
            implicit none
            class(MeshFile_t)          :: mesh
!
!           ---------------
!           Local variables
!           ---------------
!
            integer                             :: edID
            class(IntegerArrayEntry_t), pointer :: current => NULL()
            integer                             :: simpleEdge            ( POINTS_PER_EDGE            ) 
            integer                             :: elementsOfSimpleEdge  ( QUADS_PER_EDGE             ) 
            integer                             :: dividedEdge           ( POINTS_PER_SUBDIVIDED_EDGE ) 
            integer                             :: elementsOfDividedEdge ( QUADS_PER_SUBDIVIDED_EDGE  ) 
!
!           Initialize the structure
!           ------------------------
            mesh % elements_of_edges = ConstructIntegerArrayLinkedList(type = WITH_REPEATING )
!
!           Get the list head
!           -----------------
            current => mesh % points_of_edges % head
!
!           Loop the list
!           -------------
            do edID = 1 , mesh % points_of_edges % no_of_entries

               if ( current % N .eq. POINTS_PER_EDGE ) then
!
!                 Get a simple edge
!                 -----------------
                  simpleEdge = current % val
!
!                 Compute the elements that share the edge
!                 ----------------------------------------
                  elementsOfSimpleEdge = searchSimpleEdgeInElements( mesh % no_of_elements , mesh % points_of_elements , simpleEdge ) 
!
!                 Append it to the elements_of_edges list
!                 ---------------------------------------
                  if ( elementsOfSimpleEdge(TWO) .eq. -1 ) then
!
!                    Boundary edge
!                    -------------
                     call mesh % elements_of_edges % Add( ONE , elementsOfSimpleEdge(ONE) )

                  else
!
!                    Interior edge
!                    -------------
                     call mesh % elements_of_edges % Add( QUADS_PER_EDGE , elementsOfSimpleEdge )

                  end if

               elseif ( current % N .eq. POINTS_PER_SUBDIVIDED_EDGE ) then
!
!                 Get a divided edge
!                 ------------------
                  dividedEdge = current % val 
!
!                 Compute the elements that share the edge
!                 ----------------------------------------
                  elementsOfDividedEdge = searchDividedEdgeInElements( mesh % no_of_elements , mesh % points_of_elements , dividedEdge ) 
!
!                 Append it to the elements_of_edges list
!                 ---------------------------------------
                  call mesh % elements_of_edges % Add( QUADS_PER_SUBDIVIDED_EDGE , elementsOfDividedEdge )

               end if
!
!              Move to the next edge
!              ---------------------
               current => current % next

            end do

         end subroutine computeElementOfEdges

         subroutine computeFaceMarkers( mesh )
!
!           ********************************************************************
!
!                   This subroutine compute the edgeMarker array, containing
!               the boundary zone to which each edge belongs to.
!
!           ********************************************************************
!
            implicit none
            class(MeshFile_t)          :: mesh
!
!           ---------------
!           Local variables
!           ---------------
!
            integer                    :: edge
            integer                    :: bdryface
            integer                    :: pos
!
!           Allocate face marker array
!           --------------------------
            allocate ( mesh % edgeMarker ( mesh % points_of_edges % no_of_entries ) )
!
!           Set all faces to interior by default
!           ------------------------------------
            mesh % edgeMarker = FACE_INTERIOR
!
!           Loop in boundary faces to assign markers
!           ----------------------------------------
            do bdryface = 1 , mesh % no_of_bdryedges
!
!               Find a boundary edge in the points_of_edges list
!               ------------------------------------------------            
                pos = mesh % points_of_edges % Search (POINTS_PER_EDGE , mesh % points_of_bdryedges(:,bdryface) )
!
!               Once found, assign it its marker
!               --------------------------------                
                mesh % edgeMarker(pos) = mesh % bdrymarker_of_edges(bdryface)

            end do
        
         end subroutine computeFaceMarkers

         subroutine computeCurvedEdges( mesh )
!
!           ********************************************************************
!
!                   This subroutine computes which edge is curved
!               If the edge is stored in the opposite way than the curve, it
!               is inverted.
!
!           ********************************************************************
!
            implicit none
            class(MeshFile_t)          :: mesh
!
!           ---------------
!           Local variables            
!           ---------------
!
            integer                    :: edID
            integer                    :: position , direction
            integer                    :: pos(ONE)
            integer                    :: edge(POINTS_PER_EDGE)
!
!           curved_edges: The ID of the edges which are curved
!           --------------------------------------------------
            allocate ( mesh % curved_edges ( mesh % no_of_curvedEdges ) ) 

            do edID = 1 , mesh % no_of_curvedEdges
!
!               Get the position in which the curved edge is stored
!               ---------------------------------------------------
                pos = mesh % points_of_edges % SearchIfContained ( POINTS_PER_EDGE , mesh % curved_edges_points(:,edID) , ONE  )
                mesh % curved_Edges(edID) = pos(ONE)
!
!               Get the curved edge linear mesh nodes
!               -------------------------------------            
                edge = mesh % points_of_edges % Get ( POINTS_PER_EDGE , pos(ONE) )
!
!               Check whether their orientation is consistent, and invert if necessary
!               ----------------------------------------------------------------------            
                if ( (edge(ONE) .eq. mesh % curved_edges_points(ONE,edID)) .and. (edge(TWO) .eq. mesh % curved_edges_points(TWO,edID) )) then
!
!                   Do nothing
!                   ----------
                elseif ( (edge(TWO) .eq. mesh % curved_edges_points(ONE,edID)) .and. (edge(ONE) .eq. mesh % curved_edges_points(TWO,edID) ) ) then
                     mesh % curvilinear_coords(:,:,edID) = mesh % curvilinear_coords(:,mesh % curves_polynomialorder + 1 : 1 : -1 , edID)

                else
                    errorMessage(STD_OUT)

                end if

            end do

         end subroutine computeCurvedEdges
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              AUXILIAR SUBROUTINES
!              --------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

         function searchFace( face , Nentries , faceList , position , direction ) result ( isPresent )
            implicit none
            integer,    intent(in)     :: face(2)
            integer,    intent(in)     :: Nentries
            integer,    intent(in)     :: faceList(2,Nentries)
            integer,    intent(out)    :: position
            integer,    intent(out)    :: direction
            logical                    :: isPresent
!
!           ---------------
!           Local variables
!           ---------------
!
            integer     :: edge

            do edge = 1 , Nentries
               
               if ( (face(1) .eq. faceList(1,edge)) .and. (face(2) .eq. faceList(2,edge)) ) then
                  position = edge
                  direction = FORWARD
                  isPresent = .true.
                  return

               elseif ( (face(2) .eq. faceList(1,edge)) .and. (face(1) .eq. faceList(2,edge))) then
                  position = edge
                  direction = BACKWARD
                  isPresent = .true.
                  return

               end if
               
            end do

            position = -1
            direction = 0
            isPresent = .false.
            
         end function searchFace

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

         end subroutine changeEntry

         function searchSimpleEdgeInElements( no_of_elements , points_of_elements , edge ) result ( elements ) 
            implicit none
            integer, intent(in)        :: no_of_elements
            integer, intent(in)        :: points_of_elements(POINTS_PER_QUAD , no_of_elements )
            integer, intent(in)        :: edge(POINTS_PER_EDGE)
            integer                    :: elements(QUADS_PER_EDGE)
!
!           ---------------
!           Local variables
!           ---------------
!
            integer     :: eID
            integer     :: counter

            elements = -1
            counter = 1 
            do eID = 1 , no_of_elements
               
               if ( any(points_of_elements(:,eID) .eq. edge(ONE)) .and. any(points_of_elements(:,eID) .eq. edge(TWO)) ) then
                  elements(counter) = eID
                  counter = counter + 1 

                  if ( counter .eq. QUADS_PER_EDGE + 1 ) return

               end if
            end do

            elements(counter+1:) = -1

         end function searchSimpleEdgeInElements

         function searchDividedEdgeInElements( no_of_elements , points_of_elements , dividedEdge ) result ( elements )
            implicit none
            integer, intent(in)        :: no_of_elements
            integer, intent(in)        :: points_of_elements(POINTS_PER_QUAD , no_of_elements)
            integer, intent(in)        :: dividedEdge(POINTS_PER_SUBDIVIDED_EDGE)
            integer                    :: elements(QUADS_PER_SUBDIVIDED_EDGE)
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: eID            
!
!           Assign a default value
!           ----------------------
            elements = -1 
!
!           Loop in all elements: If an element contains edge's nodes 1 and 2 -> First  position
!                                                                 ""  1 and 3 -> Second position                                                         
!                                                                 ""  2 and 3 -> Third  position                                                         
!           ------------------------------------------------------------------------------------
            do eID = 1 , no_of_elements

               if     ( (any(points_of_elements(:,eID) .eq. dividedEdge(ONE))) .and. (any(points_of_elements(:,eID) .eq. dividedEdge( TWO ))) ) then
                  elements( ONE )  = eID
                  if ( .not. any( elements .eq. -1 ) ) return

               elseif ( (any(points_of_elements(:,eID) .eq. dividedEdge(ONE))) .and. (any(points_of_elements(:,eID) .eq. dividedEdge(THREE))) ) then
                  elements( TWO ) = eID
                  if ( .not. any( elements .eq. -1 ) ) return

               elseif ( (any(points_of_elements(:,eID) .eq. dividedEdge(TWO))) .and. (any(points_of_elements(:,eID) .eq. dividedEdge(THREE))) ) then
                  elements(THREE) = eID 
                  if ( .not. any( elements .eq. -1 ) ) return

               end if

            end do
 
         end function searchDividedEdgeInElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              DESTRUCT SUBROUTINES
!              --------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
         subroutine MeshFile_Destruct ( self ) 
            implicit none
            class(MeshFile_t)          :: self
   
            self % no_of_nodes            = ZERO
            self % no_of_elements         = ZERO
            self % no_of_edges            = ZERO
            self % no_of_bdryedges        = ZERO
            self % no_of_markers          = ZERO
            self % no_of_curvedEdges      = ZERO
            self % curves_polynomialorder = ZERO

            deallocate ( self % points_of_elements        ) 
            deallocate ( self % points_of_bdryedges       ) 
            deallocate ( self % polynomialOrder           ) 
            deallocate ( self % cumulativePolynomialOrder ) 
            deallocate ( self % edgeMarker                ) 
            deallocate ( self % points_coords             ) 
            deallocate ( self % bdryzones_names           ) 
            deallocate ( self % bdrymarker_of_edges       ) 
            if ( self % curvilinear ) then
               deallocate ( self % curved_edges_points ) 
               deallocate ( self % curved_edges        ) 
               deallocate ( self % curvilinear_coords  ) 
            end if

            call self % points_of_edges   % Destruct
            call self % elements_of_edges % Destruct

         end subroutine MeshFile_Destruct
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              DESCRIBE SUBROUTINES
!              --------------------
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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

               write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Number of curved edges: ", mesh % no_of_curvedEdges ,"."
               write(STD_OUT,'(30X,A,A35,I10,A)') "-> ","Curved edges polynomial order: ", mesh % curves_polynomialorder ,"."
               
            end if


            write(STD_OUT,'(/)')
            write(auxstr , '(I0,A)') mesh % no_of_markers , " boundary zones found"
            call SubSection_Header( trim(auxstr) )

            do zone = 1 , mesh % no_of_markers
               write(STD_OUT , '(30X,A,A,I0,A,A20)') "-> ", "Zone ",zone,": ", trim(mesh % bdryzones_names(zone))
            end do

         end subroutine DescribeMesh

end module MeshFileClass
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
