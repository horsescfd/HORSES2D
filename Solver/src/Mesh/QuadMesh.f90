module QuadMeshClass
    use SMConstants
    use NodeClass
    use FaceClass
    use InitialConditions
    use Element1DClass
    use Storage_module
    use QuadMeshDefinitions

    private
    public QuadMesh_t , InitializeMesh

    type QuadMesh_t
         integer                               :: no_of_nodes
         integer                               :: no_of_edges
         integer                               :: no_of_elements
         class(Node_t),         pointer        :: nodes(:)
         class(Face_p),         pointer        :: edges(:)           ! This is an array to pointers
         class(QuadElement_t) , pointer        :: elements(:)
         procedure(ICFcn)   , pointer , NOPASS :: IC
         contains
             procedure  :: ConstructFromFile
             procedure  :: SetInitialCondition
             procedure  :: ApplyInitialCondition
             procedure  :: SetStorage => QuadMesh_SetStorage
    end type QuadMesh_t

    interface InitializeMesh
          module procedure newMesh
    end interface InitializeMesh

    contains
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

!            **************
!            Set dimensions
!            **************
!
             self % no_of_nodes    = meshFile % no_of_nodes
             self % no_of_edges    = meshFile % no_of_nodes
             self % no_of_elements = meshFile % no_of_elements
!
!            *********************
!            Allocate the contents
!            *********************
!
             allocate( self % nodes ( self % no_of_nodes ) )
             allocate( self % edges ( self % no_of_edges ) ) 
             allocate( self % elements ( self % no_of_elements ) )
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

         end subroutine constructFromFile

         subroutine constructElementsAndEdges( self  , meshFile , spA , Storage , spI)
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
             integer                                  :: address
             integer                                  :: node
             integer                                  :: edge
             integer                                  :: eID
             type(Node_p)                             :: nodes(POINTS_PER_QUAD)
             integer                                  :: nodesID(POINTS_PER_QUAD)
             class(QuadElement_t), pointer            :: leftE , rightE , bdryE
!
!
!            ===================
!            Construct elements
!            ===================
!
             do eID = 1 , self % no_of_elements
                 do node = 1 , POINTS_PER_QUAD
                    nodes(node) % n => self % nodes ( meshFile % points_of_elements(eID , node) ) 
                 end do

                 nodesID(:) = meshFile % points_of_elements(eID , :)

                 address = ( meshFile % cumulativePolynomialOrder(eID-1) + eID-1 ) * NEC + 1 
                 call self % elements(eID) % construct( eID , nodes , nodesID , &
                        meshFile % polynomialOrder(eID) , Setup % nodes , spA , address  , storage , spI)  
             end do

!
!            ================
!            Construct edges
!            ================
!
             do edge = 1 , self % no_of_edges


             end do
             do edge = 1 , self % no_of_edges
        
                    
                 if (edge .eq. 1) then
                    bdryE => self % elements(edge)
                    rightE => NULL()
                    call constructFace ( self=self % edges(edge) % f , ID=edge , faceType=meshFile % edgeMarker(edge) , bdryElement = bdryE)
                    self % edges(edge) % f % n = -1.0_RP       ! So that the normal points towards the outside of the domain
                 elseif (edge .eq. self % no_of_edges) then
                    bdryE => self % elements(edge-1)
                    rightE => NULL()
                    call constructFace ( self=self % edges(edge) % f , ID=edge , faceType=meshFile % edgeMarker(edge) , bdryElement = bdryE )
                 else
                    leftE => self % elements(edge-1)
                    rightE => self % elements(edge)
                    call constructFace ( self=self % edges(edge) % f , ID=edge , faceType=meshFile % edgeMarker(edge) , leftElement = leftE , rightElement = rightE)
                 end if
             end do



         end subroutine constructElementsAndEdges
           
         subroutine setInitialCondition( self )
             use InitialConditions
             implicit none
             class(QuadMesh_t)        :: self
!
!            *******************************
!            Get Initial Condition procedure              
!            *******************************
!
             call InitialCondition( self % IC )
!
!            *******************************************
!            Apply the initial condition to the solution
!            *******************************************
!
             call self % applyInitialCondition()


          end subroutine setInitialCondition

          subroutine applyInitialCondition( self )
             use Physics
             implicit none
             class(QuadMesh_t)        :: self
             integer                :: eID , j

             do eID = 1 , self % no_of_elements
               do j = 0 , self % elements(eID) % Interp % N
                  self % elements(eID) % Q(j,1:NEC)  = self % IC( self % elements(eID) % x(j) ) 
               end do
             end do
             
          end subroutine applyInitialCondition

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
            
end module QuadMeshClass
