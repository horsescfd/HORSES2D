module DGSpatialDiscretizationMethods
   use SMConstants
   use QuadMeshClass
   use QuadMeshDefinitions
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
         class(QuadMesh_t)         :: mesh
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
         class(QuadMesh_t)         :: mesh
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
!         call DGSpatial_computeGradient( mesh )
!
!        ----------------------------------
!        Interpolate gradient to boundaries
!        ----------------------------------
!
!         call DGSpatial_interpolateToBoundaries( mesh , "dQ" )

      end subroutine DGSpatial_newTimeStep
!
      subroutine DGSpatial_resetQDot( mesh )
         implicit none
         class(QuadMesh_t)         :: mesh
         integer                 :: eID

         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % QDot = 0.0_RP
         end do

      end subroutine DGSpatial_resetQDot

      subroutine DGSpatial_interpolateToBoundaries( mesh , var )
         use Physics
         use MatrixOperations
         implicit none
         class(QuadMesh_t)         :: mesh
         character(len=*)        :: var
         integer                 :: eID , edID , eq
         real(kind=RP), pointer  :: variable(:,:)     ! will point to both Q or dQ in the elements, (0:N,0:N)
         real(kind=RP), pointer  :: variable_b(:)     ! will point to both Q or dQ in the faces, (0:N)
         real(kind=RP)           :: direction

         do eID = 1 , mesh % no_of_elements
   
            associate ( N => mesh % elements(eID) % spA % N , e => mesh % elements(eID) )
            do edID = 1 , EDGES_PER_QUAD
               associate( ed => e % edges(edID) % f )

               allocate( variable_b ( 0 : N ) ) 

               if ( (edID .eq. EBOTTOM) .or. (edID .eq. ERIGHT) ) then     ! Same direction as stored
                  direction = e % edgesDirection(edID)
               else                                                        ! Change the direction
                  direction = - e % edgesDirection(edID)
               end if

               do eq = 1 , NEC
!
!                 Gather the variable
!                 -------------------
                  select case (trim(var))
                     case ("Q")
                        variable(0:, 0: )   => e % Q(0:,0:,eq)
               
                        if ( direction .eq. FORWARD ) then
                           variable_b = ed % Q(0:N , eq , e % quadPosition(edID) )
                        else
                           variable_b = ed % Q(N:0:-1 , eq , e % quadPosition(edID) )
                        end if

                     case ("dxiQ")
                        variable(0: , 0: )   => e % dQ(0:,0:,eq,iX)
               
                        if ( direction .eq. FORWARD ) then
                           variable_b = ed % dQ(0:N , eq , e % quadPosition(edID) , iX)
                        else
                           variable_b = ed % dQ(N:0:-1 , eq , e % quadPosition(edID) , iX)
                        end if

                     case ("detaQ")
                        variable(0: , 0: )   => e % dQ(0:,0:,eq,iY)
               
                        if ( direction .eq. FORWARD ) then
                           variable_b = ed % dQ(0:N , eq , e % quadPosition(edID) , iY)
                        else
                           variable_b = ed % dQ(N:0:-1 , eq , e % quadPosition(edID) , iY)
                        end if


                  end select
   
!                 Compute the interpolation
!                 -------------------------
                  if ( edID .eq. EBOTTOM ) then
                     variable_b = MatrixTimesVector_F( variable , e % spA % lb(:,LEFT) )
                  elseif ( edID .eq. ERIGHT ) then
                     variable_b = MatrixTimesVector_F( variable , e % spA % lb(:,RIGHT) , trA = .true. )
                  elseif ( edID .eq. ETOP ) then    
                     variable_b = MatrixTimesVector_F( variable , e % spA % lb(:,RIGHT) )
                  elseif ( edID .eq. ELEFT ) then 
                     variable_b = MatrixTimesVector_F( variable , e % spA % lb(:,LEFT) , trA = .true. )
                  end if

!                 Return its value
!                 ----------------
                   select case (trim(var))
                     case ("Q")
               
                        if ( direction .eq. FORWARD ) then
                           ed % Q(0:N , eq , e % quadPosition(edID)) = variable_b
                        else
                           ed % Q(0:N , eq , e % quadPosition(edID)) = variable_b(N:0:-1)
                        end if

                     case ("dxiQ")
               
                        if ( direction .eq. FORWARD ) then
                           ed % dQ(0:N , eq , e % quadPosition(edID),iX) = variable_b
                        else
                           ed % dQ(0:N , eq , e % quadPosition(edID),iX) = variable_b(N:0:-1)
                        end if

                     case ("detaQ")
               
                        if ( direction .eq. FORWARD ) then
                           ed % dQ(0:N , eq , e % quadPosition(edID),iY) = variable_b
                        else
                           ed % dQ(0:N , eq , e % quadPosition(edID),iY) = variable_b(N:0:-1)
                        end if

                  end select
               
               end do

               deallocate( variable_b )
   
               end associate
            end do
      
            end associate
         end do
            
      end subroutine DGSpatial_interpolateToBoundaries

      subroutine DGSpatial_computeGradient( mesh )
         use QuadElementClass
         implicit none
         class(QuadMesh_t)         :: mesh
!        --------------------------------
         integer                 :: eID
         integer                 :: fID 
!        --------------------------------
!
!        ---------------------
!        Set gradients to zero
!        ---------------------
!
!         do eID = 1 , mesh % no_of_elements
!            mesh % elements(eID) % dQ = 0.0_RP
!         end do 
!!
!!        -------------------
!!        Perform volume loop
!!        -------------------
!!
!         do eID = 1 , mesh % no_of_elements
!            call SecondOrderMethod % dQVolumeLoop(mesh % elements(eID))
!         end do
!!
!!        -----------------
!!        Perform face loop
!!        -----------------
!!
!         do fID = 1 , mesh % no_of_edges
!            call SecondOrderMethod % dQFaceLoop(mesh % edges(fID) % f)
!         end do
!            
!!
!!        ----------------------------------------
!!        Perform the scaling with the mass matrix
!!        ----------------------------------------
!!
!         do eID = 1 , mesh % no_of_elements
!               mesh % elements(eID) % dQ = (1.0_RP / mesh % elements(eID) % hdiv2) * matmul(mesh % elements(eID) % Interp % Minv , mesh % elements(eID) % dQ)
!         end do
      end subroutine DGSpatial_computeGradient

      subroutine DGSpatial_computeQDot( mesh )
         use QuadElementClass
         implicit none
         class(QuadMesh_t)         :: mesh
!        -------------------------------
         integer                 :: eID
         integer                 :: fID
!        -------------------------------
!
!        Volume loops
!
         do eID = 1 , mesh % no_of_elements
            call FirstOrderMethod % QDotVolumeLoop( mesh % elements(eID) )
!            call SecondOrderMethod % QDotVolumeLoop( mesh % elements(eID) )
         end do
!
!        Face loops
!
         do fID = 1 , mesh % no_of_edges
            call FirstOrderMethod % QDotFaceLoop( mesh % edges(fID) % f )
!            call SecondOrderMethod % QDotFaceLoop( mesh % edges(fID) % f)
         end do
!
!        -------------------------------------------
!        Perform the scaling with the mass matrix
!        -------------------------------------------
!
!         do eID = 1 , mesh % no_of_elements
!            mesh % elements(eID) % QDot = (1.0_RP / mesh % elements(eID) % hdiv2) * matmul(mesh % elements(eID) % Interp % Minv , mesh % elements(eID) % QDot)
!         end do
!
      end subroutine DGSpatial_computeQDot

end module DGSpatialDiscretizationMethods

