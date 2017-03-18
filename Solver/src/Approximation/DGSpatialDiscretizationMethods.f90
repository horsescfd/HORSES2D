module DGSpatialDiscretizationMethods
   use SMConstants
   use Physics
   use QuadMeshClass
   use QuadMeshDefinitions
   use DGInviscidMethods
#ifdef NAVIER_STOKES
   use DGViscousMethods
#endif
   implicit none
!
   private
   public DGSpatial_Initialization  , DGSpatial_computeTimeDerivative , DGSpatial_interpolateSolutionToBoundaries
   public DGSpatial_newTimeStep
#ifdef NAVIER_STOKES
   public DGSpatial_computeGradient
#endif
!
!  ************************************
!  Inviscid and Viscous methods objects
!  ************************************
!
   class(InviscidMethod_t), pointer  :: InviscidMethod
#ifdef NAVIER_STOKES
   class(ViscousMethod_t),  pointer  :: ViscousMethod
#endif
!
!                                ***********                             
   integer, parameter         :: IQ      = 1
   integer, parameter         :: IDXIQ   = 2
   integer, parameter         :: IDETAQ  = 3
   integer, parameter         :: IFLUXES = 4
!                                ***********                             
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_Initialization()
         use Setup_class
         use Headers
         implicit none

         write(STD_OUT , '(/)')
         call Section_header("Spatial discretization overview")
!
!        Initialize Inviscid method
!        --------------------------
         InviscidMethod => InviscidMethod_Initialization()        
!
!        Initialize Viscous method
!        -------------------------
#ifdef NAVIER_STOKES
         ViscousMethod  => ViscousMethod_Initialization()
#endif
  
      end subroutine DGSpatial_Initialization

      subroutine DGSpatial_computeTimeDerivative( mesh ) 
!
!        ***************************************************
!           Subroutine that performs the spatial
!         discretization and computes the time
!         derivative QDot.
!        ***************************************************
!
         implicit none
         class(QuadMesh_t)         :: mesh
!
!        Prepare the mesh for a new iteration
!        ------------------------------------
         call DGSpatial_newTimeStep( mesh )
!
!        Compute QDot
!        ------------
         call DGSpatial_computeQDot( mesh )

      end subroutine DGSpatial_computeTimeDerivative
      
      subroutine DGSpatial_newTimeStep( mesh )
!
!        *************************************************************
!           This subroutine prepares the mesh struct
!          for a new time-step. 
!              1) Set QDot to zero
!              2) Interpolate the solution to boundaries
!              3) Compute the primitive variables
!              4) Compute the solution gradient
!              5) Interpolate the gradient to boundaries
!        *************************************************************
!
         implicit none
         class(QuadMesh_t)         :: mesh
         integer                   :: zoneID
!
!        Reset QDot
!        ----------
         call DGSpatial_resetQDot( mesh )
!
!        Interpolate solution to boundaries
!        ----------------------------------
         call DGSpatial_interpolateSolutionToBoundaries( mesh )
!
!        Update the zones solution
!        -------------------------
         do zoneID = 1 , size(mesh % zones) - 1
            call mesh % zones(zoneID) % UpdateSolution
         end do 
#ifdef NAVIER_STOKES
!
!        Compute the solution Q gradient dQ
!        ----------------------------------
         call DGSpatial_computeGradient( mesh )
!
!        Interpolate gradient to boundaries
!        ----------------------------------
         call DGSpatial_interpolateGradientsToBoundaries( mesh )
#endif

      end subroutine DGSpatial_newTimeStep
!
      subroutine DGSpatial_resetQDot( mesh )
         implicit none
         class(QuadMesh_t)         :: mesh
!        --------------------------------------
         integer                   :: eID

         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % QDot = 0.0_RP
         end do

      end subroutine DGSpatial_resetQDot

      subroutine DGSpatial_interpolateSolutionToBoundaries( mesh )
         use Physics
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(QuadMesh_t) :: mesh
!        --------------------------------------------------------------------
         integer                       :: eID , edID , eq
         class(QuadElement_t), pointer :: e
         class(Edge_t), pointer        :: ed
         integer, pointer              :: N

         do eID = 1 , mesh % no_of_elements

            e => mesh % elements(eID) 
            N => e % spA % N
!
!           Prolong the BOTTOM edge
!           -----------------------   
            ed => e % edges(EBOTTOM) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(EBOTTOM)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,LEFT) , N+1 )
            end do
!
!           Prolong the RIGHT edge. TODOO: implement MatrixTimesVectorInIndex to avoid transposes.
!           -----------------------   
            ed => e % edges(ERIGHT) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(ERIGHT)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,RIGHT) , N+1 , trA = .true.)
            end do
!
!           Prolong the TOP edge
!           -----------------------   
            ed => e % edges(ETOP) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(ETOP)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,RIGHT) , N+1 )
            end do
!
!           Prolong the LEFT edge. TODOO: implement MatrixTimesVectorInIndex to avoid transposes.
!           -----------------------   
            ed => e % edges(ELEFT) % f
            do eq = 1 , NCONS
               ed % storage(e % quadPosition(ELEFT)) % Q(0:N,eq) = MatrixTimesVector_F( e % Q(0:N,0:N,eq) , e % spA % lb(0:N,LEFT) , N+1 , trA = .true.)
            end do
        
         end do
            
      end subroutine DGSpatial_interpolateSolutionToBoundaries

#ifdef NAVIER_STOKES
      subroutine DGSpatial_interpolateGradientsToBoundaries( mesh )
!
!        ***************************************************************************
!              This routine interpolates to the edges the gradients.

!        ***************************************************************************
!
         use Physics
         use MatrixOperations
         implicit none
         class(QuadMesh_t) :: mesh
!        --------------------------------------------------------------------
         integer                 :: eID , edID , eq , iDim
         real(kind=RP), pointer  :: variable(:,:,:,:)     ! will point to dQ in the elements, (0:N,0:N,NDIM,NGRAD)
         real(kind=RP), pointer  :: variable_b(:,:,:)     ! will point to dQ in the faces, (0:N,NDIM,NGRAD,SIDE)

         do eID = 1 , mesh % no_of_elements
   
            associate ( N => mesh % elements(eID) % spA % N , e => mesh % elements(eID) )

            do edID = 1 , EDGES_PER_QUAD

               associate( ed => e % edges(edID) % f )

               allocate( variable_b ( 0 : N , 1 : NDIM , 1 : NGRAD ) ) 
!
!              Gather the variable
!              -------------------
               variable(0:,0:,1:,1:) => e % dQ(0:,0:,1:,1:)

               if ( e % spA % nodes .eq. LG ) then   
                  do eq = 1 , NGRAD
                     do iDim = 1 , NDIM
!
!                       Compute the interpolation
!                       -------------------------
                        if ( edID .eq. EBOTTOM ) then
                           variable_b(:,iDim,eq) = MatrixTimesVector_F( variable(:,:,iDim,eq) , e % spA % lb(:,LEFT) , N+1)
                        elseif ( edID .eq. ERIGHT ) then
                           variable_b(:,iDim,eq) = MatrixTimesVector_F( variable(:,:,iDim,eq) , e % spA % lb(:,RIGHT) , N+1 , trA = .true. )
                        elseif ( edID .eq. ETOP ) then    
                           variable_b(:,iDim,eq) = MatrixTimesVector_F( variable(:,:,iDim,eq) , e % spA % lb(:,RIGHT) , N+1)
                        elseif ( edID .eq. ELEFT ) then 
                           variable_b(:,iDim,eq) = MatrixTimesVector_F( variable(:,:,iDim,eq) , e % spA % lb(:,LEFT) , N+1 , trA = .true. )
                        end if
                     end do
                  end do

               elseif ( e % spA % nodes .eq. LGL ) then
!
!                    Just associate with its value
!                    -----------------------------
                     if ( edID .eq. EBOTTOM ) then
                        variable_b = variable(0:N,0,1:NDIM,1:NGRAD)
                     elseif ( edID .eq. ERIGHT ) then
                        variable_b = variable(N , 0:N,1:NDIM,1:NGRAD) 
                     elseif ( edID .eq. ETOP ) then
                        variable_b = variable(0:N,N,1:NDIM,1:NGRAD)
                     elseif ( edID .eq. ELEFT ) then
                        variable_b = variable(0,0:N,1:NDIM,1:NGRAD)
                     end if

               end if
              
!              Return its value
!              ----------------
             !  if ( e % edgesAssemblyDir(edID) .eq. FORWARD ) then
             !     ed % storage(e % quadPosition(edID)) % dQ ( 0:N , 1:NDIM , 1:NGRAD ) = variable_b
!
 !              else
  !                ed % storage(e % quadPosition(edID)) % dQ ( 0:N , 1:NDIM , 1:NGRAD ) = variable_b( N:0:-1 , 1:NDIM , 1:NGRAD )
!
 !              end if

               deallocate( variable_b )
   
               end associate
            end do
      
            end associate
         end do
            
      end subroutine DGSpatial_interpolateGradientsToBoundaries

      subroutine DGSpatial_computeGradient( mesh )
!
!        ***************************************************
!              This routine computes the solution gradient.
!           This is performed in four stages:
!              1) Initialization
!              2) Volume loops
!              3) Face loops
!              4) Scaling
!        ***************************************************
!
         use QuadElementClass
         implicit none
         class(QuadMesh_t)         :: mesh
!        --------------------------------
         integer                 :: eID
         integer                 :: fID 
         integer                 :: iDim , eq
         integer                 :: zoneID
!        --------------------------------
!
!        Set gradients to zero
!        ---------------------
         do eID = 1 , mesh % no_of_elements
            mesh % elements(eID) % dQ = 0.0_RP
         end do 
!
!        Perform volume loop
!        -------------------
         do eID = 1 , mesh % no_of_elements
            call ViscousMethod % dQVolumeLoop(mesh % elements(eID) )
         end do
!
!        Perform face loop
!        -----------------
         do fID = 1 , mesh % no_of_edges
            call ViscousMethod % dQFaceLoop(mesh % edges(fID) % f)
         end do
!
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
#endif
      subroutine DGSpatial_computeQDot( mesh )
!
!        *************************************************************
!              Once the mesh has been prepared for the Time Derivative
!           calculations, it is performed in this routine.
!           This is performed in three stages:
!              1) Volume loop ( Inviscid & Viscous )
!              2) Face loop ( Inviscid & Viscous )
!              3) Scaling
!        *************************************************************
!  
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
#ifdef NAVIER_STOKES
!
!        Update the zones gradients
!        --------------------------
         do zoneID = 1 , size(mesh % zones) - 1
            call mesh % zones(zoneID) % UpdateGradient
         end do 
#endif

!
!        **************
!        Inviscid terms
!        **************
!
!
!        Volume loops
!        ------------
         do eID = 1 , mesh % no_of_elements
            call InviscidMethod % QDotVolumeLoop( mesh % elements(eID) )
         end do
!
!        Face loops ( Perform the fluxes interpolation if needed )
!        ---------------------------------------------------------
         if ( Setup % inviscid_formulation .eq. FORMI ) then
            do fID = 1 , mesh % no_of_edges
               call InviscidMethod % QDotFaceLoopFormI( mesh % edges(fID) % f )

            end do
        
         elseif ( Setup % inviscid_formulation .eq. FORMII ) then
!            call DGSpatial_InterpolateToBoundaries( mesh , "Fluxes" )

            do fID = 1 , mesh % no_of_edges
               call InviscidMethod % QDotFaceLoopFormII( mesh % edges(fID) % f )

            end do

         end if
!
!        *************
!        Viscous terms
!        *************
!
#ifdef NAVIER_STOKES
!
!        Volume loops
!        ------------
         do eID = 1 , mesh % no_of_elements
            call ViscousMethod % QDotVolumeLoop( mesh % elements(eID) )
         end do
!
!        Face loops
!        ----------
!
!        If the method is IP, the gradients must be replaced by the (\nabla u) version
         select type (ViscousMethod)

            type is (IPMethod_t)

               do eID = 1 , mesh % no_of_elements
                  call mesh % elements(eID) % ComputeInteriorGradient
               end do
               call DGSpatial_interpolateGradientsToBoundaries( mesh )

         end select 

         do fID = 1 , mesh % no_of_edges
            call ViscousMethod % QDotFaceLoop( mesh % edges(fID) % f ) 

         end do
#endif
!
!        Perform the scaling with the mass matrix
!        ----------------------------------------
         do eID = 1 , mesh % no_of_elements
            associate( N => mesh % elements(eID) % spA % N )

            do eq = 1 , NCONS
               mesh % elements(eID) % QDot(0:N,0:N,eq) = mesh % elements(eID) % QDot(0:N,0:N,eq) * mesh % elements(eID) %  invM2Djac

            end do

            end associate

         end do

      end subroutine DGSpatial_computeQDot

end module DGSpatialDiscretizationMethods

