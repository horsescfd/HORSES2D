module FaceClass
    use SMConstants
    use NodeClass
    use Element1DClass
    implicit none

    private
    public  Face_t , BdryFace_t , Face_p , constructFace
!
!   *****************
!   Face_t definition
!   *****************
!
    type Face_t
        integer                           :: ID
        class(QuadElement_p), pointer       :: elements(:)
        real(kind=RP)                     :: n = 1.0_RP       ! Normal direction, this is to prepare for 2/3D
        type(Node_t), pointer             :: node             !      In this case, it points from LEFT towards RIGHT in interior faces
        real(kind=RP)                     :: F                !      and towards the outside of the domain in boundary faces
        real(kind=RP)                     :: G
        integer                           :: FaceType
    end type Face_t

    type, extends(Face_t)  :: BdryFace_t
        real(kind=RP), pointer            :: uB(:)
        real(kind=RP), pointer            :: gB(:)
        integer                           :: BCLocation
    end type BdryFace_t 

    type Face_p
        class(Face_t),   pointer           :: f
    end type Face_p


!
!   ========
    contains
!   ========
!
        subroutine constructFace( self , ID , faceType , leftElement , rightElement , bdryElement)
            implicit none
            class(Face_t), pointer      :: self
            integer                     :: ID
            integer                     :: faceType
            class(QuadElement_t), pointer, optional :: leftElement
            class(QuadElement_t), pointer, optional :: rightElement   
            class(QuadElement_t), pointer, optional :: bdryElement 
            
!
!           It needs to be allocated
!
            if (faceType .EQ. FACE_INTERIOR) then

                allocate(Face_t :: self)

!               Allocate its elements and point to the objects
                allocate( self % elements(2) )
                if (present(leftElement) .and. present(rightElement) ) then
                  self % elements(LEFT) % e   => leftElement
                  self % elements(RIGHT) % e  => rightElement
                end if
                  self % faceType = FACE_INTERIOR
            elseif (faceType .NE. FACE_INTERIOR) then

                allocate(Bdryface_t :: self)

                allocate( self % elements(1) ) 
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

            self % ID = ID 


        end subroutine constructFace            

        




end module FaceClass
