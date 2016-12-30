module EdgeClass
    use SMConstants
    use NodeClass
    use Element1DClass
    use QuadMeshDefinitions
    implicit none

    private
    public  Edge_t , StraightBdryEdge_t , CurvedBdryEdge_t , Edge_p , constructFace
!
!   *****************
!   Face_t definition
!   *****************
!
    type Edge_t
        integer                           :: ID
        class(QuadElement_p), pointer       :: elements(:)
        real(kind=RP)                     :: n = 1.0_RP       ! Normal direction, this is to prepare for 2/3D
        type(Node_t), pointer             :: node             !      In this case, it points from LEFT towards RIGHT in interior faces
        real(kind=RP)                     :: F                !      and towards the outside of the domain in boundary faces
        real(kind=RP)                     :: G
        integer                           :: FaceType
    end type Edge_t

    type, extends(Edge_t)  :: StraightBdryEdge_t
        real(kind=RP), pointer            :: uB(:)
        real(kind=RP), pointer            :: gB(:)
        integer                           :: BCLocation
    end type StraightBdryEdge_t 

    type, extends(Edge_t)  :: CurvedBdryEdge_t

    end type CurvedBdryEdge_t

    type Edge_p
        class(Edge_t),   pointer           :: f
        contains
            procedure         :: construct => constructFace
    end type Edge_p


!
!   ========
    contains
!   ========
!
        subroutine constructFace( self , ID , curvilinear , faceType , leftElement , rightElement , bdryElement)
            implicit none
            class(Edge_p)               :: self
            integer                     :: ID
            integer                     :: faceType
            logical                     :: curvilinear
            class(QuadElement_t), pointer, optional :: leftElement
            class(QuadElement_t), pointer, optional :: rightElement   
            class(QuadElement_t), pointer, optional :: bdryElement 
            
!
!           *************************************************
!              Allocate the edge depending on its type
!           *************************************************
!
            if (faceType .EQ. FACE_INTERIOR) then

                allocate(Face_t :: self % f)

!               Allocate its elements and point to the objects
                allocate( self % elements(QUADS_PER_EDGE) )
                if (present(leftElement) .and. present(rightElement) ) then
                  self % f % elements(LEFT) % e   => leftElement
                  self % f % elements(RIGHT) % e  => rightElement
                end if
                  self % faceType = FACE_INTERIOR

            elseif (faceType .NE. FACE_INTERIOR) then

               if ( .NOT. curvilinear ) then
                  allocate(StraightBdryEdge_t   :: self % f)

               else
                  allocate(CurvedBdryEdge_t     :: self % f)

               end if

               allocate( self % f % elements(1) )

                if (present(bdryElement)) then
                  self % elements(1) % e => bdryElement
                end if 
                self % faceType = faceType

                select type( self ) 
                  type is (BdryFace_t)
!
!                 ---------------------------------------------
!                  There is a need to indicate whether a 
!                    right or left face is for the element            
!                 ---------------------------------------------
!       
                     if ( self % elements(1) % e % facesID(LEFT) .eq. ID ) then
                        self % BCLocation = LEFT
                     elseif ( self % elements(1) % e % facesID(RIGHT) .eq. ID ) then
                        self % BCLocation = RIGHT
                     else
                        print*, "The face does not belong to the element."
                        stop "Stopped."
                     end if
                end select

            end if



        end subroutine constructFace            

        




end module EdgeClass
