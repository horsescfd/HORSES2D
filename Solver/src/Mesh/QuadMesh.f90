module QuadMeshClass
    use SMConstants
    use NodeClass
    use QuadElementClass
    use InitialConditions
    use Storage_module
    use QuadMeshDefinitions
    use DGBoundaryConditions

    private
    public Zone_t , QuadMesh_t , InitializeMesh

    integer, parameter        :: STR_LEN_MESH = 128

    type QuadMesh_t
         integer                               :: no_of_nodes
         integer                               :: no_of_edges
         integer                               :: no_of_elements
         class(Node_t),         pointer        :: nodes(:)
         class(Edge_p),         pointer        :: edges(:)           ! This is an array to pointers
         class(QuadElement_t) , pointer        :: elements(:)
         class(Zone_t)        , pointer        :: zones(:)
         procedure(ICFcn)   , pointer , NOPASS :: IC
         real(kind=RP)                         :: Volume
         contains
             procedure  :: ConstructFromFile
             procedure  :: ConstructZones            => Mesh_ConstructZones
             procedure  :: SetInitialCondition
             procedure  :: ApplyInitialCondition
             procedure  :: SetStorage                  => QuadMesh_SetStorage
             procedure  :: VolumeIntegral              => Compute_VolumeIntegral
             procedure  :: ScalarScalarSurfaceIntegral => Compute_ScalarScalarSurfaceIntegral
             procedure  :: ScalarVectorSurfaceIntegral => Compute_ScalarVectorSurfaceIntegral
             procedure  :: VectorVectorSurfaceIntegral => Compute_VectorVectorSurfaceIntegral
             procedure  :: TensorVectorSurfaceIntegral => Compute_TensorVectorSurfaceIntegral
             procedure  :: ComputeResiduals            => Mesh_ComputeResiduals
             procedure  :: ComputePrimitiveVariables   => Mesh_ComputePrimitiveVariables
             procedure  :: ComputeMaxJumps             => Mesh_ComputeMaxJumps
             procedure  :: FindElementWithCoords       => Mesh_FindElementWithCoords
    end type QuadMesh_t

    type Zone_t
       integer                             :: marker
       character(len=STR_LEN_MESH)         :: Name
       integer                             :: no_of_edges
       class(Edge_p), pointer              :: edges(:)
       class(BoundaryCondition_t), pointer :: BC
       contains
          procedure      :: Construct      => Zone_Construct
          procedure      :: UpdateSolution => Zone_UpdateSolution
#ifdef NAVIER_STOKES
          procedure      :: UpdateGradient => Zone_UpdateGradient
#endif
          procedure      :: Describe       => Zone_Describe
    end type Zone_t
 

    interface InitializeMesh
          module procedure newMesh
    end interface InitializeMesh

    contains

#include "ZoneProcedures.incf"
#include "QuadMeshIntegrals.incf"

         function newMesh()
             implicit none
             type(QuadMesh_t)    :: newMesh
!
!            **************************
!            Set to zero all dimensions            
!            **************************
!
             newMesh % no_of_nodes = 0
             newMesh % no_of_edges = 0
             newMesh % no_of_elements = 0
!
!            **************************
!            Point all contents to NULL
!            **************************
!
             newMesh % nodes    => NULL()
             newMesh % edges    => NULL()
             newMesh % elements => NULL() 

         end function newMesh
!
!        *********************************************************************
!           Subroutine to build the mesh structure from the MeshFile_t
!              Class already loaded.
!        *********************************************************************
!
         subroutine constructFromFile( self , meshFile , spA , Storage , spI)
             use MeshFileClass
             use Setup_class
             use Physics
             use NodesAndWeights_Class
             implicit none
             class(QuadMesh_t),                 intent (out)                 :: self
             class(MeshFile_t),                 intent (in )                 :: meshFile
             class(NodalStorage),               intent (in )                 :: spA
             class(Storage_t),                  intent (in )                 :: storage
             class(NodesAndWeights_t), pointer, intent (in )                 :: spI
!            ----------------------------------------------------------------------
             integer                                                         :: node

!            **************
!            Set dimensions
!            **************
!
             self % no_of_nodes    = meshFile % no_of_nodes
             self % no_of_edges    = meshFile % no_of_edges
             self % no_of_elements = meshFile % no_of_elements
!
!            *********************
!            Allocate the contents
!            *********************
!
             allocate ( self % nodes    ( self % no_of_nodes    )  ) 
             allocate ( self % edges    ( self % no_of_edges    )  ) 
             allocate ( self % elements ( self % no_of_elements )  ) 
!
!            ***********************
!            Construct the contents
!            ***********************
!
!
!            ================
!            Construct nodes
!            ================
!
             do node = 1 , self % no_of_nodes
                 call self % nodes(node) % construct( ID = node, x = meshFile % points_coords(1:NDIM,node))
             end do
!
!            ============================
!            Construct edges and elements
!            ============================
!
             call constructElementsAndEdges( self , meshFile , spA, Storage , spI )

             self % Volume = self % VolumeIntegral("One")            

         end subroutine constructFromFile

         subroutine constructElementsAndEdges( self  , meshFile , spA , Storage , spI)
             use MeshFileClass
             use Setup_class
             use Physics
             use NodesAndWeights_Class
             use MatrixOperations
             implicit none
             class(QuadMesh_t),                 intent (inout)               :: self
             class(MeshFile_t),                 intent (in )                 :: meshFile
             class(NodalStorage),               intent (in )                 :: spA
             class(Storage_t),                  intent (in )                 :: storage
             class(NodesAndWeights_t), pointer, intent (in )                 :: spI
!            ----------------------------------------------------------------------
             integer                                  :: address
             integer                                  :: node
             integer                                  :: edge
             integer                                  :: eID
             integer                                  :: el1 , el2 , elb
             integer                                  :: curve
             type(Node_p)                             :: nodes(POINTS_PER_QUAD)
             class(QuadElement_t), pointer            :: leftE , rightE , bdryE
             logical                                  :: curvilinear
!            ----------------------------------------------------------------------
!
!            ===================
!            Construct elements
!            ===================
!
             do eID = 1 , self % no_of_elements
                 
                 do node = 1 , POINTS_PER_QUAD
                    nodes(node) % n => self % nodes ( meshFile % points_of_elements(node , eID) ) 
                 end do


                 address = ( meshFile % cumulativePolynomialOrder(eID-1)  ) * NCONS + 1 
                 call self % elements(eID) % Construct( eID , nodes , meshFile % polynomialOrder(eID) , spA , address , storage , spI ) 

             end do

!
!            ================
!            Construct edges
!            ================
!
             do edge = 1 , self % no_of_edges

               do node = 1 , POINTS_PER_EDGE
                  nodes(node) % n => self % nodes ( meshFile % points_of_edges(node , edge) )
               end do

               if (meshFile % curvilinear) then
                  curvilinear = any(meshFile % curved_bdryedges == edge)
               else
                  curvilinear = .false.
               end if

               call self % edges(edge) % Construct( ID = edge , curvilinear = curvilinear , &
                                                nodes = nodes , edgeType = meshFile % edgeMarker(edge) , spA = spA , spI = spI )

               if (curvilinear) then
!
            
!              Add the curve to the edge
!              -------------------------
                  select type ( f => self % edges(edge) % f )

                     type is (Edge_t)
                        call self % edges(edge) % f % SetCurve()

                     type is (StraightBdryEdge_t) 
                        call self % edges(edge) % f % SetCurve()

                     type is (CurvedBdryEdge_t) 
                        curve = minloc(abs(meshFile % curved_bdryedges -  edge) , 1)
                        call self % edges(edge) % f % SetCurve( meshFile % curvilinear_coords(:,:,curve) , meshFile % curves_polynomialorder )  

                     class default

                  end select

               else
                  call self % edges(edge) % f % SetCurve()

               end if
             end do
!
!             ========================
!             Link elements with edges
!             ========================
!
              do edge = 1 , self % no_of_edges

                  if (self % edges(edge) % f % edgeType .eq. FACE_INTERIOR) then
                     el1 = meshFile % elements_of_edges( 1 , edge )
                     el2 = meshFile % elements_of_edges( 2 , edge )

                     call self % edges(edge) % linkWithElements( el1 = self % elements(el1) , el2 = self % elements(el2) )

                  else
                     
                     elb = meshFile % elements_of_edges( 1 , edge )
                     call self % edges(edge)  % linkWithElements( elb = self % elements(elb) )

                  end if
               end do

!
!              =================================
!              Compute the geometry of the quads
!              =================================
!
               do eID = 1 , self % no_of_elements
                  call self % elements(eID) % SetMappings
               end do
!
!              Compute areas and volumes
!              -------------------------
               do eID = 1 , self % no_of_elements
                  associate ( e => self % elements(eID) )
                  e % Volume = BilinearForm_F ( e % jac , e % spA % w , e % spA % w )
                  end associate
               end do

               do edge = 1 , self % no_of_edges
                  associate ( ed => self % edges(edge) % f )
                  ed % Area = dot_product( norm2(ed % dS , dim = 1) , ed % spA % w )
                  end associate
               end do

         end subroutine constructElementsAndEdges
           
         subroutine SetInitialCondition( self , which)
             use InitialConditions
             implicit none
             class(QuadMesh_t)            :: self
             character(len=*), optional   :: which
!
!            *******************************
!            Get Initial Condition procedure              
!            *******************************
!
             if (present(which)) then
               call InitialCondition( self % IC , which)

             else
               call InitialCondition( self % IC )

             end if
!
!            *******************************************
!            Apply the initial condition to the solution
!            *******************************************
!
             call self % ApplyInitialCondition()


          end subroutine setInitialCondition

          subroutine ApplyInitialCondition( self , argin)
             use Physics
             implicit none
             class(QuadMesh_t)        :: self
             real(kind=RP), optional  :: argin
             integer                  :: eID
             integer                  :: iXi
             integer                  :: iEta
             real(kind=RP)            :: X(NDIM)

             if ( present ( argin ) ) then
                do eID = 1 , self % no_of_elements
                  do iXi = 0 , self % elements(eID) % spA % N
                     do iEta = 0 , self % elements(eID) % spA % N
                        X = self % elements(eID) % X(iXi,iEta,IX:IY)
                        self % elements(eID) % Q(iXi,iEta,1:NCONS)  = self % IC( X , argin) 
                     end do
                  end do
                end do

            else
                do eID = 1 , self % no_of_elements
                  do iXi = 0 , self % elements(eID) % spA % N
                     do iEta = 0 , self % elements(eID) % spA % N
                        X = self % elements(eID) % X(iXi,iEta,IX:IY)
                        self % elements(eID) % Q(iXi,iEta,1:NCONS)  = self % IC( X ) 
                     end do
                  end do
                end do
   
            end if

          end subroutine ApplyInitialCondition

          subroutine QuadMesh_SetStorage( self , storage )
            use Storage_module
            implicit none
            class(QuadMesh_t)             :: self
            class(Storage_t)            :: storage
            integer                     :: eID

            do eID = 1 , self % no_of_elements
                call self % elements(eID) % SetStorage( storage )
            end do

         end subroutine QuadMesh_SetStorage

        subroutine Mesh_ConstructZones( self , meshFile  )
            use MeshFileClass
            use Headers
            implicit none
            class(QuadMesh_t)                        :: self
            class(MeshFile_t)                        :: meshFile
            character(len=STR_LEN_MESH), allocatable :: zoneNames(:)
            integer                                  :: zone

            write(STD_OUT,'(/)') 
            call Section_header("Boundary conditions overview")

            allocate( self % Zones( 0 : meshFile % no_of_markers ) )
            allocate( zoneNames( 0 :  meshFile % no_of_markers) )

            zoneNames(0) = "Interior"
            zoneNames(1 : meshFile % no_of_markers) = meshFile % bdryzones_names

            do zone = 0 , meshFile % no_of_markers
               call self % Zones(zone) % Construct( self , zone , zoneNames(zone) )
            end do
!
!           For periodic boundary conditions: It is neccesary to perform the linking
!           ------------------------------------------------------------------------
            do zone = 1 , meshFile % no_of_markers
               select type ( BC => self % Zones(zone) % BC ) 
                  type is (PeriodicBC_t)
                     if ( .not. BC % associated ) then
                        call Zone_LinkPeriodicZones( self % Zones(zone) , self % Zones( BC % connected_marker ) ) 
                     end if
                  class default
               end select
            end do

!
         end subroutine Mesh_ConstructZones

         subroutine Mesh_ComputePrimitiveVariables( self )
            implicit none
            class(QuadMesh_t)          :: self
            integer                    :: eID , edID

            do eID = 1 , self % no_of_elements

               call self % elements(eID) % ComputePrimitiveVariables

            end do

            do edID = 1 , self % no_of_edges

               call self % edges(edID) % f % ComputePrimitiveVariables

            end do

         end subroutine Mesh_ComputePrimitiveVariables
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              EXTRA SUBROUTINES
!              ----------------- 
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
         function Mesh_ComputeResiduals( self ) result ( residuals )
            use Physics
            implicit none
            class(QuadMesh_t )               :: self
            real(kind=RP)                    :: residuals(NCONS)
            integer                          :: eID

            residuals = 0.0_RP
         
            do eID = 1 , self % no_of_elements
               residuals(IRHO)  = max( residuals(IRHO) , maxval(abs(self % elements(eID) % QDot(:,:,IRHO))) )
               residuals(IRHOU) = max( residuals(IRHOU) , maxval(abs(self % elements(eID) % QDot(:,:,IRHOU))) )
               residuals(IRHOV) = max( residuals(IRHOV) , maxval(abs(self % elements(eID) % QDot(:,:,IRHOV))) )
               residuals(IRHOE) = max( residuals(IRHOE) , maxval(abs(self % elements(eID) % QDot(:,:,IRHOE))) )
            end do   

         end function Mesh_ComputeResiduals

         subroutine Mesh_FindElementWithCoords( self , x , elemID , xi , eta )
            use Physics
            implicit none
            class(QuadMesh_t) ,  intent (in)  :: self
            real(kind=RP)     ,  intent (in)  :: x(NDIM)
            integer           ,  intent (out) :: elemID
            real(kind=RP)     ,  intent (out)  :: xi
            real(kind=RP)     ,  intent (out)  :: eta
!           ------------------------------------------------------------------
            integer                           :: eID , edID
            logical                           :: isInside
            real(kind=RP)                     :: distance , minimumDistance

elloop:     do eID = 1 , self % no_of_elements

               isInside = self % elements(eID) % FindPointWithCoords( x , xi , eta )

               if ( isInside ) then
                  elemID = eID
                  exit elloop
      
               end if

            end do elloop

            if ( .not. isInside ) then
               print*, "Warning, the point probe was not found in the mesh."
               elemID = -1
               xi = huge(0.0_RP)
               eta = huge(0.0_RP)

            end if

         end subroutine Mesh_FindElementWithCoords

         function Mesh_ComputeMaxJumps( self ) result ( val ) 
            use Physics
            implicit none
            class(QuadMesh_t)          :: self
!           --------------------------------------------------------------
            integer                    :: edID
            real(kind=RP)              :: val
            real(kind=RP)              :: localJumps

            val = 0.0_RP 

            do edID = 1 , self % no_of_edges

               associate ( N => self % edges(edID) % f % spA % N )
               select type ( f => self % edges(edID) % f )


                  type is (Edge_t)
                     localJumps = maxval( abs ( f % Q(0:N,1:NCONS,LEFT) - f % Q(0:N,1:NCONS,RIGHT) ) )
   
                  type is (StraightBdryEdge_t)
                     localJumps = maxval( abs ( f % Q(0:N,1:NCONS,1) - f % uB(0:N,1:NCONS) ) ) 

                  type is (CurvedBdryEdge_t)
                     localJumps = maxval( abs ( f % Q(0:N,1:NCONS,1) - f % uB(0:N,1:NCONS) ) ) 
   
               end select
               end associate


               if ( localJumps .gt. val ) then
                  val = localJumps
      
               end if
            end do

         end function Mesh_ComputeMaxJumps

end module QuadMeshClass   
