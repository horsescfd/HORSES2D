module DGSpatialDiscretizationMethods
   use Mesh1DClass
   use DGSecondOrderMethods
   use DGFirstOrderMethods
   implicit none

   private
   public DGSpatial_Initialization , DGSpatial_computeTimeDerivative , DGSpatial_interpolateToBoundaries
   public DGSpatial_computeGradient

   class(SecondOrderMethod_t), pointer     :: SecondOrderMethod
   class(FirstOrderMethod_t), pointer      :: FirstOrderMethod
!
!  ========
   contains
!  ========
!
      subroutine DGSpatial_Initialization()
         use Setup_class
         implicit none

         FirstOrderMethod => FirstOrderMethod_Initialization()        
         SecondOrderMethod => SecondOrderMethod_Initialization()
  
      end subroutine DGSpatial_Initialization

      subroutine DGSpatial_computeTimeDerivative( mesh ) 
!        ------------------------------------------+
!           Subroutine that performs the spatial
!         discretisation and computes the time
!         derivative QDot
!        ------------------------------------------+
         implicit none
         class(Mesh1D_t)         :: mesh
!
!        -------------------------------------------
!           Prepare the mesh for a new iteration
!        -------------------------------------------
!
         call DGSpatial_newTimeStep( mesh )
!
!        -------------------------------------------
!           Compute QDot
!        -------------------------------------------
!
         call DGSpatial_computeQDot( mesh )

      end subroutine DGSpatial_computeTimeDerivative
      
      subroutine DGSpatial_newTimeStep( mesh )
!        --------------------------------------------------
!           This subroutine prepares the mesh struct
!          for a new time-step. 
!              1) Set QDot to zero
!              2) Compute the solution gradient
!              3) Interpolate the solution and gradient
!                    to boundaries.
!        --------------------------------------------------
         implicit none
         class(Mesh1D_t)         :: mesh
!
!        ----------
!        Reset QDot
!        ----------
!
         call DGSpatial_resetQDot( mesh )
!
!        ----------------------------------
!        Interpolate solution to boundaries
!        ----------------------------------
!
         call DGSpatial_interpolateToBoundaries( mesh , "Q" )

!
!        ----------------------------------
!        Compute the solution Q gradient dQ
!        ----------------------------------
!
         call DGSpatial_computeGradient( mesh )
!
!        ----------------------------------
!        Interpolate gradient to boundaries
!        ----------------------------------
!
         call DGSpatial_interpolateToBoundaries( mesh , "dQ" )

      end subroutine DGSpatial_newTimeStep
!
      subroutine DGSpatial_resetQDot( mesh )
         implicit none
         class(Mesh1D_t)         :: mesh
         integer                 :: eID

         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % QDot = 0.0_RP
         end do

      end subroutine DGSpatial_resetQDot

      subroutine DGSpatial_interpolateToBoundaries( mesh , var )
         use Physics
         use MatrixOperations
         implicit none
         class(Mesh1D_t)         :: mesh
         character(len=*)        :: var
         integer                 :: eID
         real(kind=RP), pointer  :: variable(:,:)     ! will point to both Q or dQ, (0:N , NEC)
         real(kind=RP), pointer  :: variable_b(:,:)     ! will point to both Qb or dQb, (2 , NEC)
!
         do eID = 1 , mesh % no_of_elements
            select case (trim(var))
               case ("Q")
                  variable => mesh % elements(eID) % Q
                  variable_b => mesh % elements(eID) % Qb
               case ("dQ")
                  variable => mesh % elements(eID) % dQ
                  variable_b => mesh % elements(eID) % dQb
            end select

            call TransposeMat_x_NormalMat( variable , mesh % elements(eID) % Interp % lb , variable_b )

         end do
            
      end subroutine DGSpatial_interpolateToBoundaries

      subroutine DGSpatial_computeGradient( mesh )
         use Element1DClass
         use FaceClass
         implicit none
         class(Mesh1D_t)         :: mesh
!        --------------------------------
         integer                 :: eID
         integer                 :: fID 
!        --------------------------------
!
!        ---------------------
!        Set gradients to zero
!        ---------------------
!
         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % dQ = 0.0_RP
         end do 
!
!        -------------------
!        Perform volume loop
!        -------------------
!
         do eID = 1 , mesh % no_of_elements
            call SecondOrderMethod % dQVolumeLoop(mesh % elements(eID))
         end do
!
!        -----------------
!        Perform face loop
!        -----------------
!
         do fID = 1 , mesh % no_of_faces
            call SecondOrderMethod % dQFaceLoop(mesh % faces(fID) % f)
         end do
            
!
!        ----------------------------------------
!        Perform the scaling with the mass matrix
!        ----------------------------------------
!
         do eID = 1 , mesh % no_of_elements
               mesh % elements(eID) % dQ = (1.0_RP / mesh % elements(eID) % hdiv2) * matmul(mesh % elements(eID) % Interp % Minv , mesh % elements(eID) % dQ)
         end do
      end subroutine DGSpatial_computeGradient

      subroutine DGSpatial_computeQDot( mesh )
         use Element1DClass
         use FaceClass
         implicit none
         class(Mesh1D_t)         :: mesh
!        -------------------------------
         integer                 :: eID
         integer                 :: fID
!        -------------------------------
!
!        Volume loops
!
         do eID = 1 , mesh % no_of_elements
            call FirstOrderMethod % QDotVolumeLoop( mesh % elements(eID) )
            call SecondOrderMethod % QDotVolumeLoop( mesh % elements(eID) )
         end do
!
!        Face loops
!
         do fID = 1 , mesh % no_of_faces
            call FirstOrderMethod % QDotFaceLoop( mesh % faces(fID) % f )
            call SecondOrderMethod % QDotFaceLoop( mesh % faces(fID) % f)
         end do
!
!        -------------------------------------------
!        Perform the scaling with the mass matrix
!        -------------------------------------------
!
         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % QDot = (1.0_RP / mesh % elements(eID) % hdiv2) * matmul(mesh % elements(eID) % Interp % Minv , mesh % elements(eID) % QDot)
         end do

      end subroutine DGSpatial_computeQDot

end module DGSpatialDiscretizationMethods

