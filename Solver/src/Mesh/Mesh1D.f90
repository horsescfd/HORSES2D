module mesh1DClass
    use SMConstants
    use NodeClass
    use FaceClass
    use InitialConditions
    use Element1DClass
    use Storage_module

    private
    public Mesh1D_t , InitializeMesh

    type Mesh1D_t
         integer                               :: no_of_nodes
         integer                               :: no_of_faces
         integer                               :: no_of_elements
         class(Node_t)      , pointer          :: nodes(:)
         class(Face_p)      , pointer          :: faces(:)           ! This is an array to pointers
         class(Element1D_t) , pointer          :: elements(:)
         procedure(ICFcn)   , pointer , NOPASS :: IC
         contains
             procedure  :: ConstructFromFile
             procedure  :: SetInitialCondition
             procedure  :: ApplyInitialCondition
             procedure  :: SetStorage => Mesh1D_SetStorage
    end type Mesh1D_t

    interface InitializeMesh
          module procedure newMesh
    end interface InitializeMesh

    contains
         function newMesh()
             implicit none
             type(Mesh1D_t)    :: newMesh
!
!            **************************
!            Set to zero all dimensions            
!            **************************
!
             newMesh % no_of_nodes = 0
             newMesh % no_of_faces = 0
             newMesh % no_of_elements = 0
!
!            **************************
!            Point all contents to NULL
!            **************************
!
             newMesh % nodes    => NULL()
             newMesh % faces    => NULL()
             newMesh % elements => NULL() 

         end function newMesh

         subroutine constructFromFile( self , meshFile , spA , Storage , spI)
             use MeshFileClass
             use Setup_class
             use Physics
             use NodesAndWeights_Class
             implicit none
             class(Mesh1D_t)                   :: self
             class(MeshFile_t)                 :: meshFile
             class(NodalStorage)               :: spA
             class(Storage_t)                  :: storage
             class(NodesAndWeights_t), pointer :: spI
!
             integer                :: address
             integer                :: node
             integer                :: face
             integer                :: eID
             class(Node_t), pointer :: leftNode , rightNode
             class(Element1D_t), pointer :: leftE , rightE , bdryE
!
!            **************
!            Set dimensions
!            **************
!
             self % no_of_nodes = meshFile % no_of_nodes
             self % no_of_faces = meshFile % no_of_nodes
             self % no_of_elements = meshFile % no_of_elements
!
!            *********************
!            Allocate the contents
!            *********************
!
             allocate( self % nodes ( self % no_of_nodes ) )
             allocate( self % faces ( self % no_of_faces ) ) 
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
                 call self % nodes(node) % construct( ID = node, x = meshFile % nodes(node))
             end do
!
!            ===================
!            Construct elements
!            ===================
!
             do eID = 1 , self % no_of_elements
                 leftNode => self % nodes( meshFile % elements(eID , LEFT) )
                 rightNode => self % nodes ( meshFile % elements( eID , RIGHT) )
                 address = ( meshFile % cumulativePolynomialOrder(eID-1) + eID-1 ) * NEC + 1 
                 call self % elements(eID) % construct( eID , leftNode , rightNode , leftNode % ID , rightNode % ID , &
                        meshFile % polynomialOrder(eID) , Setup % nodes , spA , address  , storage , spI)  
             end do

!
!            ================
!            Construct faces
!            ================
!
             do face = 1 , self % no_of_faces
        
                    
                 if (face .eq. 1) then
                    bdryE => self % elements(face)
                    rightE => NULL()
                    call constructFace ( self=self % faces(face) % f , ID=face , faceType=meshFile % faceType(face) , bdryElement = bdryE)
                    self % faces(face) % f % n = -1.0_RP       ! So that the normal points towards the outside of the domain
                 elseif (face .eq. self % no_of_faces) then
                    bdryE => self % elements(face-1)
                    rightE => NULL()
                    call constructFace ( self=self % faces(face) % f , ID=face , faceType=meshFile % faceType(face) , bdryElement = bdryE )
                 else
                    leftE => self % elements(face-1)
                    rightE => self % elements(face)
                    call constructFace ( self=self % faces(face) % f , ID=face , faceType=meshFile % faceType(face) , leftElement = leftE , rightElement = rightE)
                 end if
             end do

         end subroutine constructFromFile
           
         subroutine setInitialCondition( self )
             use InitialConditions
             implicit none
             class(Mesh1D_t)        :: self
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
             class(Mesh1D_t)        :: self
             integer                :: eID , j

             do eID = 1 , self % no_of_elements
               do j = 0 , self % elements(eID) % Interp % N
                  self % elements(eID) % Q(j,1:NEC)  = self % IC( self % elements(eID) % x(j) ) 
               end do
             end do
             
          end subroutine applyInitialCondition

          subroutine Mesh1D_SetStorage( self , storage )
            use Storage_module
            implicit none
            class(Mesh1D_t)             :: self
            class(Storage_t)            :: storage
            integer                     :: eID

            do eID = 1 , self % no_of_elements
                call self % elements(eID) % SetStorage( storage )
            end do

         end subroutine Mesh1D_SetStorage
            
end module mesh1DClass
