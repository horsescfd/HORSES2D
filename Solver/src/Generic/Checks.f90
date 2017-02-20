module ChecksModule
   use SMConstants
   use Physics
   implicit none
!
!  ========   
   contains
!  ========   
!
!
      subroutine checks( sem )
          use DGSEM_Class
          use Headers
          use NodesAndWeights_class
          use QuadMeshClass
          use MeshFileClass
          use Setup_class
          use DGSpatialDiscretizationMethods
          use Storage_module
          use DGBoundaryConditions  
          implicit none
          class(DGSEM_t)                :: sem
          integer, parameter            :: STR_LEN_CHECKS = 128
          integer                       :: fID
          integer                       :: eID
          integer                       :: command_argument_count
          integer                       :: nArgs
          character(len=STR_LEN_CHECKS) :: arg
          integer                       :: iArg
          logical                       :: perform_tests = .false.
         
          nArgs = command_argument_count()
          do iArg = 1 , nArgs
            call get_command_argument(iArg , arg)
            if (trim(arg) .eq. "-check") then
               perform_tests = .true.
               exit
            end if
          end do
      
          if (perform_tests) then
      
            write(STD_OUT , '(/)')
            call Section_Header("Performing tests on the built framework") 
            write(STD_OUT , '(/)')
   
            call CheckMappings( sem % mesh )

            if ( trim ( Setup % IC ) .ne. "Restart" ) then
               call CheckInterpolationToBoundaries( sem % mesh ) 
            end if

            call CheckMetricIdentities( sem % mesh )

            call checkQDot( sem )

            call checkGradients ( sem ) 

         end if

        end subroutine checks
      
        subroutine CheckMappings( mesh )
            use QuadElementClass
            use Headers
            use QuadMeshClass
            use QuadMeshDefinitions
            use MatrixOperations
            implicit none
            type(QuadMesh_t)           :: mesh
!           ---------------------------------------
            integer                    :: eID , edID
            real(kind=RP), allocatable :: dxiX(:,:,:)
            real(kind=RP), allocatable :: detaX(:,:,:)
            real(kind=RP)              :: error
            integer                    :: current , location
            integer                    :: zone
            real(kind=RP), allocatable :: dSx(:) , dSy(:)
            real(kind=RP), allocatable :: dSe(:,:)
            integer                    :: which(NDIM)

            call SubSection_Header("Testing the mappings")
            
!           This is to test the elements mappings derivatives formula
!           ---------------------------------------------------------
            do eID = 1 , mesh % no_of_elements
               associate (e => mesh % elements(eID) )

               if ( allocated (dxiX ) ) deallocate ( dxiX  ) 
               if ( allocated (detaX) ) deallocate ( detaX ) 
          
               allocate(dxiX  ( 0 : e % spA % N , 0 : e % spA % N , NDIM) ) 
               allocate(detaX ( 0 : e % spA % N , 0 : e % spA % N , NDIM) ) 

               dxiX         = MatrixMultiplyInIndex_F( e % X , e % spA % DT  , 1)
               detaX        = MatrixMultiplyInIndex_F( e % X , e % spA % DT , 2)

               if (eID .eq. 1) then
                  current = eID
                  error = maxval(abs(dxiX - e % dX(:,:,:,IX))) 
                   
               else

                  if ( maxval(abs(dxiX - e % dX(:,:,:,iX) ) ) .gt. error ) then
                     error = maxval(abs(dxiX - e % dX(:,:,:,IX) ) )
                     current = eID
                  end if
               end if
                     
               if ( maxval(abs(detaX - e % dX(:,:,:,iY)) ) .gt. error) then
                     error = maxval(abs(detaX - e % dX(:,:,:,iY) ) ) 
                     current = eID
               end if
               end associate

            end do

            write(STD_OUT , '(30X,A,A,ES10.3,A,I0,A)') "-> ", "Maximum error found in elements mapping: ",error,"  (Cell ",current,")."

!
!           This is to test the extrapolation to boundaries of the normal vectors from elements
!           -----------------------------------------------------------------------------------
            error = 0.0_RP
            current = 0

            do eID = 1 , mesh % no_of_elements
               if ( allocated(dSx) ) deallocate(dSx) 
               if ( allocated(dSy) ) deallocate(dSy) 
               if ( allocated(dSe) ) deallocate(dSe) 
               associate( e => mesh % elements(eID) )
               
               allocate( dSx ( 0 : e % spA % N ) )
               allocate( dSy ( 0 : e % spA % N ) )
               allocate( dSe ( NDIM , 0 : e % spA % N ) )

!              BOTTOM Edge
!              -----------
               which = [1,2]
               dSx = -MatrixTimesVector_F( e % Ja(which) , e % spA % lj(0.0_RP)  )
               which = [2,2]
               dSy = -MatrixTimesVector_F( e % Ja(which) , e % spA % lj(0.0_RP)  )

               if ( e % edgesDirection(EBOTTOM) .eq. FORWARD ) then
                  dSe = e % edges(EBOTTOM) % f % dS 
               else
                  dSe = - e % edges(EBOTTOM) % f % dS(iX:iY , e % spA % N : 0 : -1 )
               end if
               
               if ( maxval(abs(dSx - dSe(iX,:) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,:) ) )
                  current = eID
                  location = EBOTTOM
               end if
               if (  maxval(abs(dSy - dSe(iY,:) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,:) ) )
                  current = eID 
                  location = EBOTTOM
               end if

!              RIGHT Edge
!              -----------
               which = [1,1]
               dSx = MatrixTimesVector_F( e % Ja(which) , e % spA % lj(1.0_RP) , trA = .true.  )
               which = [2,1]
               dSy = MatrixTimesVector_F( e % Ja(which) , e % spA % lj(1.0_RP) , trA = .true.  )
               
               if ( e % edgesDirection(ERIGHT) .eq. FORWARD ) then
                  dSe = e % edges(ERIGHT) % f % dS 
               else
                  dSe = - e % edges(ERIGHT) % f % dS(iX:iY , e % spA % N : 0 : -1 )
               end if
               
               if ( maxval(abs(dSx - dSe(iX,:) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,:) ) )
                  current = eID
                  location = ERIGHT
               end if
               if (  maxval(abs(dSy - dSe(iY,:) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,:) ) )
                  current = eID 
                  location = ERIGHT
               end if              

!              TOP Edge
!              -----------
               which = [1,2]
               dSx = MatrixTimesVector_F( e % Ja(which) , e % spA % lj(1.0_RP)  )
               which = [2,2]
               dSy = MatrixTimesVector_F( e % Ja(which) , e % spA % lj(1.0_RP)  )

               if ( e % edgesDirection(ETOP) .eq. FORWARD ) then
                  dSe = e % edges(ETOP) % f % dS 
               else
                  dSe = - e % edges(ETOP) % f % dS(iX:iY , e % spA % N : 0 : -1 )
               end if
               
               if ( maxval(abs(dSx - dSe(iX,e % spA % N : 0 : -1) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,e % spA % N : 0 : -1) ) )
                  current = eID
                  location = ETOP
               end if
               if (  maxval(abs(dSy - dSe(iY,e % spA % N:0:-1) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,e % spA % N : 0 : -1) ) )
                  current = eID 
                  location = ETOP
               end if              

!              LEFT Edge
!              -----------
               which = [1,1]
               dSx = -MatrixTimesVector_F( e % Ja(which) , e % spA % lj(0.0_RP)  , trA = .true.)
               which = [2,1]
               dSy = -MatrixTimesVector_F( e % Ja(which) , e % spA % lj(0.0_RP)  , trA = .true.)

               if ( e % edgesDirection(ELEFT) .eq. FORWARD ) then
                  dSe = e % edges(ELEFT) % f % dS 
               else
                  dSe = - e % edges(ELEFT) % f % dS(iX:iY , e % spA % N : 0 : -1 )
               end if
               
               if ( maxval(abs(dSx - dSe(iX,e % spA % N : 0 : -1) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,e % spA % N : 0 : -1) )) 
                  current = eID
                  location = ELEFT
               end if
               if (  maxval(abs(dSy - dSe(iY,e % spA % N : 0 : -1) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,e % spA % N : 0 : -1) ) )
                  current = eID 
                  location = ELEFT
               end if              

              
               end associate
            end do

            write(STD_OUT , '(30X,A,A,ES10.3,A,I0,A,I0,A)') "-> ", "Maximum error found in edges mapping: ",error,"  (Cell ",current,", in " , location , ")."

!           Compute the volume of the domain
            write(STD_OUT , '(30X,A,A35,F16.10,A)') "-> ", "Computed domain volume: " , mesh % VolumeIntegral("One") * RefValues % L**2.0_RP,"."

!           Compute faces surface            
            do zone = 1 , size(mesh % Zones) - 1
               write(STD_OUT,'(30X,A,A35,F16.10,A)') "-> ", "Computed surface in zone " // trim(mesh % Zones(zone) % Name) // ": ",mesh % SurfaceIntegral("One",zone) * RefValues % L ,"." 
            end do

        end subroutine CheckMappings

        subroutine CheckInterpolationToBoundaries( mesh ) 
          use HEaders
          use DGSpatialDiscretizationMethods
          use QuadMeshClass
          use QuadElementClass
          implicit none
          class(QuadMesh_t)            :: mesh
          integer                      :: eID , edID , quad
          real(kind=RP)                :: error = 0.0_RP
          real(kind=RP)                :: currentError = 0.0_RP
          integer                      :: iXi , iEta 
 
          write(STD_OUT,'(/)')
          call SubSection_Header("Checking the interpolation to boundaries")

          call DGSpatial_interpolateToBoundaries( mesh ,"Q")

          do eID = 1 , mesh % no_of_elements
            do iXi = 0 , mesh % elements(eID) % spA % N
               do iEta = 0 , mesh % elements(eiD) % spA % N
               
                  currentError = norm2( mesh % elements(eID) % Q(iXi,iEta,:) - mesh % IC(mesh % elements(eID) % x(:,iXi,iEta) ) ) 
                  if (currentError .gt. error) then
                     error = currentError
                  end if

               end do
            end do
            
          end do

          write(STD_OUT , '(30X,A,A50,ES16.10,A)') "-> ", "Initial condition interpolation error in quads: " , error,"."

          error = 0.0_RP
      
          do edID = 1 , mesh % no_of_edges
            do quad = 1 , size(mesh % edges(edID) % f % quads)
               do iXi = 0 , mesh % edges(edID) % f % spA % N
                  currentError = norm2( mesh % edges(edID) % f % Q(iXi,:,quad) - mesh % IC(mesh % edges(edID) % f % x(:,iXi) ) )

                  if (currentError .gt. error) then
                     error = currentError
                  end if
               end do
            end do
          end do
      
          write(STD_OUT , '(30X,A,A50,ES16.10,A)') "-> ", "Initial condition interpolation error in edges: " , error,"."

        end subroutine CheckInterpolationToBoundaries

        subroutine CheckMetricIdentities( mesh ) 
         use SMConstants
         use Headers
         use Physics
         use QuadMeshClass
         use QuadElementClass
         use Setup_class
         use MatrixOperations
         implicit none
         class(QuadMesh_t)          :: mesh
!        --------------------------------------------
         integer                    :: eID
         integer                    :: coord
         integer                    :: which(NDIM)
         real(kind=RP)              :: error = 0.0_RP , currenterror = 0.0_RP
         real(kind=RP), allocatable :: Ja1(:,:) , Ja2(:,:)
         real(kind=RP), allocatable :: metricID(:,:)

     
         write(STD_OUT,'(/)')
         call Subsection_Header("Checking discrete metric identities")

         do eID = 1 , mesh % no_of_elements

            associate( e => mesh % elements(eID) )

            allocate( Ja1(0:e % spA % N , 0 : e % spA % N ) )
            allocate( Ja2(0:e % spA % N , 0 : e % spA % N ) )
            allocate( metricID(0:e % spA % N , 0 : e % spA % N ) )
            do coord = 1 , NDIM
               
               which = [coord,1]
               Ja1 = e % Ja(which) 
               which = [coord,2]
               Ja2 = e % Ja(which)

               metricID = Mat_x_Mat_F( A = e % spA % D , B = Ja1 ) + Mat_x_Mat_F( A = Ja2 , B = e % spA % DT )
               
               currenterror = abs(BilinearForm_F( A = metricID , X = e % spA % w , Y = e % spA % w ))
               
               if ( currenterror .gt. error ) then
                  error = currenterror
               end if

            end do

            deallocate(Ja1 , Ja2 , metricID)

            end associate
         end do

         write(STD_OUT , '(30X,A,A50,F16.10,A)') "-> ", "Maximum discrete metric identities residual: " , error,"."

        end subroutine CheckMetricIdentities
   
        subroutine Integration_checks( sem ) 
          use DGSEM_Class
          use SMConstants
          use Physics
          use NodesAndWeights_class
          use QuadMeshClass
          use MeshFileClass
          use Setup_class
          use DGSpatialDiscretizationMethods
          use Storage_module
          use DGBoundaryConditions  
          implicit none
          class(DGSEM_t) :: sem
          real(kind=RP), allocatable      :: Manalytical(:,:) , Mquadrature(:,:)
          real(kind=RP) , allocatable     :: a(:)
          real(kind=RP) , allocatable     :: b(:)
          real(kind=RP) , allocatable     :: aM(:)
          real(kind=RP) , allocatable     :: bM(:)
          integer                         :: i , j
          write(STD_OUT ,'(/)')
          write(STD_OUT , * ) "This subroutine performs checks on the numerical quadratures enforced. 3N polynomials are computed and stored into a matrix which entries are:"
         write(STD_OUT ,'(40X , A)') "M_{ij} = x^i (x^j)^2"
      
         allocate( Manalytical( 0 : sem % spA % head % N , 0 : sem % sPA % head % N ) )
         allocate( Mquadrature( 0 : sem % spA % head % N , 0 : sem % sPA % head % N ) )
         allocate( a( 0 : sem % spA % head % N ) )
         allocate( b( 0 : sem % spA % head % N ) )
      
         do i = 0 , sem % spA % head % N
            do j = 0 , sem % spA % head % N
               Manalytical(i,j) = ((1.0_RP+0.1_RP)**(2.0_RP*j+i+1.0_RP)- (-1.0_RP+0.1_RP)**(2.0_RP*j+i+1.0_RP))/(2.0_RP*j + i + 1.0_RP )
            end do
         end do 
      
         if ( Setup % inviscid_discretization .eq. "Over-Integration" ) then
                 allocate( aM( 0 : sem % spI % N ) )
                 allocate( bM( 0 : sem % spI % N ) )
         end if
      
         do i = 0 , sem % spA % head % N
            do j = 0 , sem % spA % head % N
               a = (sem % spA % head % xi + 0.1_RP)**(1.0_RP * i)
               b = (sem % spA % head % xi + 0.1_RP)**(1.0_RP * j)
            
               if ( Setup % inviscid_discretization .eq. "Standard") then
                Mquadrature(i,j) = sum( sem % spA % head % w * a * b**2.0_RP)
               elseif ( Setup % inviscid_discretization .eq. "Over-Integration") then
                  
                  aM = matmul( sem % spA % head % T , a)
                  bM = matmul( sem % spA % head % T , b)
      
                  Mquadrature(i,j) = sum( sem % spI % w * aM * bM**2.0_RP)
      
               end if
            end do
         end do
      
         write(STD_OUT , * ) "Quadratures analysis:"
         write(STD_OUT , '(20X , A , I0)' ) "Polynomial order: " , sem % spA % head % N
         if ( Setup % nodes .eq. LG) then
         write(STD_OUT , '(20X , A , I0)' ) "Expected number of quadrature points: ", ceiling((3.0_RP * sem % spA % head % N - 1.0_RP )/2.0_RP)
         elseif ( Setup % nodes .eq. LGL) then
         write(STD_OUT , '(20X , A , I0)' ) "Expected number of quadrature points: ", ceiling((3.0_RP * sem % spA % head % N + 1.0_RP )/2.0_RP)
         end if 
      
         if ( associated (sem % spI ) ) then
         write(STD_OUT , '(20X , A , I0)' ) "Quadrature points: " , sem % spI % N
         else
         write(STD_OUT , '(20X , A , I0)' ) "Quadrature points: " , sem % spA % head % N
         endif 
         print*, "Maximum error found in quadratures: " , maxval( abs( (Manalytical - Mquadrature) / (Manalytical + Mquadrature )) ) 

        end subroutine Integration_checks

        subroutine CheckQDot( sem )
         use DGSEM_Class
         use DGSpatialDiscretizationMethods
         use Headers
         use QuadMeshDefinitions
         implicit none
         class(DGSem_t)          :: sem
         integer                 :: iXi , iEta
         integer                 :: eID , elem = -1 , zoneID
         real(kind=RP)           :: error = 0.0_RP , localerror
         real(kind=RP)           :: x(NDIM)
         real(kind=RP)           :: L 
         real(kind=RP)           :: QDot(NCONS)
         logical                 :: elementIsInterior

         write(STD_OUT,'(/)')
         call SubSection_Header("Testing QDot calculation")

!
!        Apply the "ChecksPolynomic" initial condition
!        ------------------------------------
         L = sqrt( sem % mesh % VolumeIntegral("One") ) / 4.0_RP
         call sem % mesh % SetInitialCondition("ChecksPolynomic")
         call sem % mesh % ApplyInitialCondition( L )    

         call DGSpatial_ComputeTimeDerivative( sem % mesh )

         error = 0.0_RP
         do eID = 1 , sem % mesh % no_of_elements
            associate( e => sem % mesh % elements(eID) ,  N => sem % mesh % elements(eID) % spA % N )

            do iXi = 0 , N
               do iEta = 0 , N

                  x = e % X(iXi , iEta , IX:IY)

                  QDot = QDotPolynomicFCN( x , L )

                  localerror = maxval(abs([e % QDot(iXi,iEta,1:4) - QDot(1:4)]))

                  elementIsInterior = .true.

                  if ( e % edges(EBOTTOM) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.
                  if ( e % edges(ETOP) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.
                  if ( e % edges(ELEFT) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.
                  if ( e % edges(ERIGHT) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.

                  if ( (localerror .gt. error) .and. elementIsInterior ) then
                     error = localerror
                     elem = eID
                  end if

               end do
            end do

            end associate
         end do
         write(STD_OUT , '(30X,A,A60,ES16.10,A,I0,A)') "-> ", "Polynomic initial condition time derivative error: " , error," (cell  ",elem,")."
!
!        Apply the "ChecksTrigonometric" initial condition
!        ------------------------------------
         L = sqrt( sem % mesh % VolumeIntegral("One") ) / 4.0_RP
         call sem % mesh % SetInitialCondition("ChecksTrigonometric")
         call sem % mesh % ApplyInitialCondition( L )


         call DGSpatial_ComputeTimeDerivative( sem % mesh )

         error = 0.0_RP
         do eID = 1 , sem % mesh % no_of_elements
            associate( e => sem % mesh % elements(eID) ,  N => sem % mesh % elements(eID) % spA % N )

            do iXi = 0 , N
               do iEta = 0 , N

                  x = e % X(iXi , iEta , IX:IY)
                  QDot = QDotTrigonometricFCN( x , L )

                  localerror = maxval(abs([e % QDot(iXi,iEta,1:4) - QDot(1:4)]))
                  
                  elementIsInterior = .true.
                  if ( e % edges(EBOTTOM) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.
                  if ( e % edges(ETOP) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.
                  if ( e % edges(ELEFT) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.
                  if ( e % edges(ERIGHT) % f % edgeType .ne. FACE_INTERIOR ) elementIsInterior = .false.

                  if ( (localerror .gt. error) .and. elementIsInterior ) then
                     error = localerror
                     elem = eID
                  end if

               end do
            end do

            end associate
         end do
         write(STD_OUT , '(30X,A,A60,ES16.10,A,I0,A)') "-> ", "Trigonometric initial condition time derivative error: " , error," (cell  ",elem,")."
!
!        Return to the problem initial condition and RiemannSolvers
!        ----------------------------------------------------------
         call sem % SetInitialCondition ( verbose = .false. )

        end subroutine CheckQDot
   
        subroutine CheckGradients( sem ) 
          use DGSEM_Class
          use SMConstants
          use Physics
          use NodesAndWeights_class
          use QuadMeshClass
          use MeshFileClass
          use Setup_class
          use DGSpatialDiscretizationMethods
          use Storage_module
          use DGBoundaryConditions  
          use Headers
          implicit none
          class(DGSem_t)         :: sem
          integer                :: eID
          real(kind=RP)          :: error = 0.0_RP , localerror = 0.0_RP
          integer                :: elem = -1
          integer                :: iXi, iEta
          real(kind=RP)          :: x(NDIM)
          real(kind=RP)          :: dQ(NDIM , NGRAD)
          real(kind=RP)          :: L 

         write(STD_OUT,'(/)')
         call SubSection_Header("Testing Gradients")
!
!        Set the polynomic initial condition
!        -----------------------------------
         call sem % mesh % SetInitialCondition("ChecksPolynomic")
         call sem % mesh % ApplyInitialCondition

          call DGSpatial_newTimeStep( sem % mesh )


          do eID = 1 , sem % mesh % no_of_elements
            do iXi = 0 , sem % mesh % elements(eID) % spA % N
               do iEta = 0 , sem % mesh % elements(eID) % spA % N

                  x = sem % mesh % elements(eID) % x(iXi,iEta,IX:IY) 
                  dQ = dQPolynomicFcn(x)
                  localerror = maxval(abs(sem % mesh % elements(eID) % dQ(iXi,iEta,1:NDIM,1:NGRAD) - dQ ) )

                  if ( localerror .gt. error ) then
                     error = localerror
                     elem = eID
                  end if

               end do
            end do
         end do

         write(STD_OUT , '(30X,A,A50,ES16.10,A,I0,A)') "-> ", "Error in gradients for polynomic state: " , error," (cell  ",elem,")."
!
!        Set the trigonometric initial condition
!        -----------------------------------
         L = sqrt( sem % mesh % VolumeIntegral("One") ) / 4.0_RP
         call sem % mesh % SetInitialCondition("ChecksTrigonometric")
         call sem % mesh % ApplyInitialCondition( L )

         call DGSpatial_newTimeStep( sem % mesh )

         error = 0.0_RP
         elem = -1

          do eID = 1 , sem % mesh % no_of_elements
            do iXi = 0 , sem % mesh % elements(eID) % spA % N
               do iEta = 0 , sem % mesh % elements(eID) % spA % N

                  x = sem % mesh % elements(eID) % x(iXi,iEta,IX:IY) 
                  dQ = dQTrigonometricFcn(x , L )
                  localerror = maxval(abs(sem % mesh % elements(eID) % dQ(iXi,iEta,1:NDIM,1:NGRAD) - dQ ) ) * L

                  if ( localerror .gt. error ) then
                     error = localerror
                     elem = eID
                  end if

               end do
            end do
         end do

         write(STD_OUT , '(30X,A,A50,ES16.10,A,I0,A)') "-> ", "Error in gradients for trigonometric state: " , error," (cell  ",elem,")."

!
!        Return to the problem initial condition and RiemannSolvers
!        ----------------------------------------------------------
         call sem % SetInitialCondition( verbose = .false. )

        end subroutine CheckGradients
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           CHECK AUXILIAR ROUTINES
!           -----------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        function dQPolynomicFcn ( x ) result ( val )
         use Physics
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: val(NDIM,NGRAD)
         
         associate ( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )

         val(IX:IY , IGU) = [sqrt(gamma)*Mach , 0.0_RP]
         val(IX:IY , IGV) = [0.0_RP , sqrt(gamma)*Mach]
         val(IX:IY , IGT) = 2.0_RP * gamma * Mach * Mach * [ x(IX) , x(IY) ]

         end associate

        end function dQPolynomicFcn

        function dQTrigonometricFcn ( x , L ) result ( val )
         use Physics
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)        :: L 
         real(kind=RP)        :: val(NDIM,NGRAD)
         
         associate ( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )

         val(IX:IY , IGU) = sqrt(gamma)*Mach*PI * [ cos(PI*x(IX)/L)*cos(PI*x(IY)/L) , -sin(PI*x(IX)/L)*sin(PI*x(IY)/L)] / L
         val(IX:IY , IGV) = sqrt(gamma)*Mach*PI * [ sin(PI*x(IX)/L)*sin(PI*x(IY)/L) , -cos(PI*x(IX)/L)*cos(PI*x(IY)/L)] / L
         val(IX:IY , IGT) = -0.25_RP * gamma * Mach * Mach * PI * [ sin(2.0_RP * PI * x(IX)/L) , sin(2.0_RP * PI * x(IY)/L) ] / L 

         end associate

        end function dQTrigonometricFcn

        function QDotTrigonometricFCN( x , L ) result( val )
         use Physics
         implicit none
         real(kind=RP)           :: x(NDIM)
         real(kind=RP)           :: L 
         real(kind=RP)           :: val(NCONS)
         real(kind=RP)           :: u , v , p
         real(kind=RP)           :: ux , vy , H , uy , vx , px , py , Hx , Hy


         associate( gamma => Thermodynamics % gamma , Mach => dimensionless % Mach , cp => Dimensionless % cp)

         u = sqrt(gamma) * Mach * sin(PI * x(IX) / L ) * cos(PI * x(IY) / L)
         v = -sqrt(gamma)* Mach * cos(PI * x(IX) / L ) * sin(PI * x(IY) / L)
         p = 1.0_RP + (gamma * Mach**2.0_RP / 8.0_RP) * ( cos(2.0_RP * PI * x(IX) / L ) + cos(2.0_RP * PI * x(IY) / L) )
         H = cp * p + 0.5_RP * u**2.0_RP + 0.5_RP * v**2.0_RP


         ux = sqrt(gamma) * Mach * PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L
         uy = -sqrt(gamma)* Mach * PI * sin(PI * x(IX) / L) * sin(PI * x(IY) / L ) / L

         vx = sqrt(gamma) * Mach * PI * sin(PI * x(IX) / L) * sin(PI * x(IY) / L ) / L 
         vy = -sqrt(gamma) * Mach * PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L 

         px = -0.25_RP * gamma * Mach**2.0_RP * PI * sin(2.0_RP * PI * x(IX) / L ) / L 
         py = -0.25_RP * gamma * Mach**2.0_RP * PI * sin(2.0_RP * PI * x(IY) / L ) / L 

         Hx = cp * px + u * ux + v * vx
         Hy = cp * py + u * uy + v * vy

         
   
         val(IRHO)      = 0.0_RP
         val(IRHOU)     = 2.0_RP * u * ux + u*vy + uy * v + px
         val(IRHOV)     = 2.0_RP * v * vy + py + u*vx + ux*v
         val(IRHOE) = H * (ux + vy) + Hx * u + Hy * v
!
         val = -val / (sqrt(gamma) * Mach)
         
         end associate

        end function QDotTrigonometricFCN

        function QDotPolynomicFCN( x , L ) result( val )
         use Physics
         implicit none
         real(kind=RP)           :: x(NDIM) 
         real(kind=RP)           :: L
         real(kind=RP)           :: val(NCONS)
         real(kind=RP)           :: u , v , p
         real(kind=RP)           :: ux , vy , H , uy , vx , px , py , Hx , Hy

         associate( gamma => Thermodynamics % gamma , Mach => dimensionless % Mach , cp => Dimensionless % cp)

         u = sqrt(gamma) * Mach * x(IX) / L
         v = sqrt(gamma) * Mach * x(IY) / L
         p = 1.0_RP + gamma * Mach * Mach * ( x(IX) * x(IX) + x(IY) * x(IY)) / (L **2.0_RP)
         H = cp * p + 0.5_RP * u**2.0_RP + 0.5_RP * v**2.0_RP

         ux = sqrt(gamma) * Mach / L
         uy = 0.0_RP

         vx = 0.0_RP
         vy = sqrt(gamma) * Mach / L 

         px = gamma * Mach * Mach * 2 * x(IX) / L**2.0_RP
         py = gamma * Mach * Mach * 2 * x(IY) / L**2.0_RP

         Hx = cp * px + u * ux + v * vx
         Hy = cp * py + u * uy + v * vy

         val(IRHO)      = ux + vy
         val(IRHOU)     = 2.0_RP * u * ux + u*vy + uy * v + px
         val(IRHOV)     = 2.0_RP * v * vy + py + u*vx + ux*v
         val(IRHOE)     = H * (ux + vy) + Hx * u + Hy * v
!
         val = -val / (sqrt(gamma) * Mach)

         end associate

         associate ( gamma => thermodynamics % gamma , Mach => dimensionless % Mach , mu => dimensionless % mu , kappa => dimensionless % kappa )

         val(IRHOE) = gamma * Mach * Mach / (L*L) * ( 4.0_RP/3.0_RP * mu + 4.0_RP * kappa ) 


         end associate

        end function QDotPolynomicFCN        
      
end module
