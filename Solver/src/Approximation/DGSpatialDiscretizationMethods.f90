module DGSpatialDiscretizationMethods
   use SMConstants
   use Physics
   use QuadMeshClass
   use QuadMeshDefinitions
   use DGViscousMethods
   use DGInviscidMethods
   implicit none

   private
   public DGSpatial_Initialization , DGSpatial_computeTimeDerivative , DGSpatial_interpolateToBoundaries
   public DGSpatial_computeGradient , DGSpatial_newTimeStep

   class(ViscousMethod_t), pointer     :: ViscousMethod
   class(InviscidMethod_t), pointer      :: InviscidMethod


   integer, parameter         :: IQ = 1
   integer, parameter         :: IDXIQ = 2
   integer, parameter         :: IDETAQ = 3
   integer, parameter         :: IFLUXES = 4
!
!  ========
   contains
!  ========
!
      subroutine DGSpatial_Initialization()
         use Setup_class
         use Headers
         implicit none

         write(STD_OUT , '(/)')
         call Section_header("Spatial discretization overview")
         InviscidMethod => InviscidMethod_Initialization()        
         ViscousMethod => ViscousMethod_Initialization()
  
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
!              3) Compute the primitive variables
!              4) Interpolate the solution and gradient
!                    to boundaries.
!        --------------------------------------------------
         implicit none
         class(QuadMesh_t)         :: mesh
         real(kind=RP)              :: tstart , tend
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
!        ---------------------------
!        Compute primitive variables
!        ---------------------------
!
         call mesh % computePrimitiveVariables
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

      recursive subroutine DGSpatial_interpolateToBoundaries( mesh , var )
         use Physics
         use MatrixOperations
         implicit none
         class(QuadMesh_t)         :: mesh
         character(len=*)        :: var
         integer                 :: eID , edID , eq
         real(kind=RP), pointer  :: variable(:,:,:)     ! will point to both Q or dQ in the elements, (0:N,0:N)
         real(kind=RP), pointer  :: variable_b(:,:)     ! will point to both Q or dQ in the faces, (0:N)
         integer                 :: varID
         real(kind=RP)           :: edgeSign
         
         select case ( trim(var) )
            case ("Q")
               varID = IQ
            case ("dxiQ")
               varID = IDXIQ
            case ("detaQ")
               varID = IDETAQ
            case ("Fluxes")
               varID = IFLUXES
            case default
               varID = -1
         end select

         do eID = 1 , mesh % no_of_elements
   
            associate ( N => mesh % elements(eID) % spA % N , e => mesh % elements(eID) )
            do edID = 1 , EDGES_PER_QUAD
               associate( ed => e % edges(edID) % f )

               allocate( variable_b ( 0 : N , NCONS ) ) 

               if ( ( edID .eq. ETOP ) .or. (edID .eq. EBOTTOM) ) then
                  edgeSign = -1.0_RP          ! Outside-pointing edges
               else
                  edgeSign = 1.0_RP         ! Inside-pointing edges
               end if
!
!                 Gather the variable
!                 -------------------
                  select case ( varID )
                     case (IQ)
                        variable(0:, 0:,1: )   => e % Q(0:,0:,1:)
                     case (IDXIQ)
                        variable(0: , 0:,1: )   => e % dQ(0:,0:,1:,iX)
                     case (IDETAQ)
                        variable(0: , 0:,1: )   => e % dQ(0:,0:,1:,iY)
                     case (IFLUXES)
                        if ( (edID .eq. EBOTTOM) .or. (edID .eq. ETOP) ) then
                           variable(0: , 0:,1: )   => e % F(0:,0:,1:,iY)
                        elseif ( (edID .eq. ELEFT) .or. (edID .eq. ERIGHT) ) then
                           variable(0: , 0:,1: )   => e % F(0:,0:,1:,iX)
                        end if

                  end select

               do eq = 1 , NCONS

                  if ( e % spA % nodes .eq. LG ) then   
!
!                    Compute the interpolation
!                    -------------------------
                     if ( edID .eq. EBOTTOM ) then
                        variable_b(:,eq) = MatrixTimesVector_F( variable(:,:,eq) , e % spA % lb(:,LEFT) )
                     elseif ( edID .eq. ERIGHT ) then
                        variable_b(:,eq) = MatrixTimesVector_F( variable(:,:,eq) , e % spA % lb(:,RIGHT) , trA = .true. )
                     elseif ( edID .eq. ETOP ) then    
                        variable_b(:,eq) = MatrixTimesVector_F( variable(:,:,eq) , e % spA % lb(:,RIGHT) )
                     elseif ( edID .eq. ELEFT ) then 
                        variable_b(:,eq) = MatrixTimesVector_F( variable(:,:,eq) , e % spA % lb(:,LEFT) , trA = .true. )
                     end if

                  end if
               end do
               if ( e % spA % nodes .eq. LGL ) then
!
!                    Just associate with its value
!                    -----------------------------
                     if ( edID .eq. EBOTTOM ) then
                        variable_b = variable(0:N,0,1:NCONS)
                     elseif ( edID .eq. ERIGHT ) then
                        variable_b = variable(N , 0:N,1:NCONS) 
                     elseif ( edID .eq. ETOP ) then
                        variable_b = variable(0:N,N,1:NCONS)
                     elseif ( edID .eq. ELEFT ) then
                        variable_b = variable(0,0:N,1:NCONS)
                     end if

               end if
              
!                 Return its value
!                 ----------------
                select case (varID)
                  case (IQ)
            
                     if ( e % edgesAssemblyDir(edID) .eq. FORWARD ) then
                        ed % Q(0:N , 1:NCONS , e % quadPosition(edID)) = variable_b
                     else
                        ed % Q(0:N , 1:NCONS , e % quadPosition(edID)) = variable_b(N:0:-1,1:NCONS)
                     end if

                  case (IDXIQ)
            
                     if ( e % edgesAssemblyDir(edID) .eq. FORWARD ) then
                        ed % dQ(0:N , 1:NCONS , e % quadPosition(edID),iX) = variable_b
                     else
                        ed % dQ(0:N , 1:NCONS , e % quadPosition(edID),iX) = variable_b(N:0:-1,1:NCONS)
                     end if

                  case (IDETAQ)
            
                     if ( e % edgesAssemblyDir(edID) .eq. FORWARD ) then
                        ed % dQ(0:N , 1:NCONS , e % quadPosition(edID),iY) = variable_b
                     else
                        ed % dQ(0:N , 1:NCONS , e % quadPosition(edID),iY) = variable_b(N:0:-1,1:NCONS)
                     end if

                  case (IFLUXES)
   
                     if ( e % edgesAssemblyDir(edID) .eq. FORWARD ) then
                        ed % F (0:N , 1:NCONS , e % quadPosition(edID)) = edgeSign * variable_b
                     elseif ( e % edgesAssemblyDir(edID) .eq. BACKWARD ) then
                        ed % F (0:N , 1:NCONS , e % quadPosition(edID)) = -edgeSign * variable_b(N:0:-1,1:NCONS)     ! To ensure that is consistent with the edge normal
                     end if

               end select
 


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
         integer                 :: iDim , eq
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
            call ViscousMethod % dQVolumeLoop(mesh % elements(eID))
         end do
!
!        -----------------
!        Perform face loop
!        -----------------
!
         do fID = 1 , mesh % no_of_edges
            call ViscousMethod % dQFaceLoop(mesh % edges(fID) % f)
         end do
!
!        ----------------------------------------
!        Perform the scaling with the mass matrix
!        ----------------------------------------
         do eID = 1 , mesh % no_of_elements
            associate( N => mesh % elements(eID) % spA % N )
            do eq = 1 , NGRAD
               do iDim = 1 , NDIM
                  mesh % elements(eID) % dQ(0:N,0:N,iDim,eq) = mesh % elements(eID) % dQ(0:N,0:N,iDim,eq) * mesh % elements(eID) %  invM2Djac
               end do
            end do
            end associate
        end do


      end subroutine DGSpatial_computeGradient

      subroutine DGSpatial_computeQDot( mesh )
         use QuadElementClass
         use Setup_class
         implicit none
         class(QuadMesh_t)         :: mesh
!        -------------------------------
         integer                 :: eID
         integer                 :: fID
         integer                 :: zoneID
         integer                 :: eq
!        -------------------------------
!
!        Volume loops
!
         do eID = 1 , mesh % no_of_elements
            call InviscidMethod % QDotVolumeLoop( mesh % elements(eID) )
!            call ViscousMethod % QDotVolumeLoop( mesh % elements(eID) )
         end do

!
!        **********
!        Face loops
!        **********
!
!
!        Update the contents
!        -------------------
         do zoneID = 1 , size(mesh % zones) - 1
            call mesh % zones(zoneID) % Update
         end do 
!
!        Interpolate the fluxes and perform the face loops
!        -------------------------------------------------
         if ( Setup % inviscid_formulation .eq. FORMI ) then
            do fID = 1 , mesh % no_of_edges
               call InviscidMethod % QDotFaceLoopFormI( mesh % edges(fID) % f )
!               call ViscousMethod % QDotFaceLoop( mesh % edges(fID) % f)
            end do
        
         elseif ( Setup % inviscid_formulation .eq. FORMII ) then

            call DGSpatial_InterpolateToBoundaries( mesh , "Fluxes" )

            do fID = 1 , mesh % no_of_edges
               call InviscidMethod % QDotFaceLoopFormII( mesh % edges(fID) % f )
!               call ViscousMethod % QDotFaceLoop( mesh % edges(fID) % f)
            end do

         end if
!
!        -------------------------------------------
!        Perform the scaling with the mass matrix
!        -------------------------------------------
!
         do eID = 1 , mesh % no_of_elements
            associate( N => mesh % elements(eID) % spA % N )
            do eq = 1 , NCONS
               mesh % elements(eID) % QDot(0:N,0:N,eq) = mesh % elements(eID) % QDot(0:N,0:N,eq) * mesh % elements(eID) %  invM2Djac
            end do
            end associate
         end do

      end subroutine DGSpatial_computeQDot

end module DGSpatialDiscretizationMethods

