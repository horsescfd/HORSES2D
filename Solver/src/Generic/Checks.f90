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

            call CheckMesh( sem % mesh )
   
            call CheckMappings( sem % mesh )

            if ( trim ( Setup % IC ) .ne. "Restart" ) then
               call CheckInterpolationToBoundaries( sem % mesh ) 
            end if

            call CheckMetricIdentities( sem % mesh )

            call checkQDot( sem )
#ifdef NAVIER_STOKES
            call checkGradients ( sem ) 
#endif

         end if

        end subroutine checks

        subroutine CheckMesh( mesh )
            use QuadElementClass
            use Headers
            use QuadMeshClass
            use QuadMeshDefinitions
            use MatrixOperations
            implicit none
            class(QuadMesh_t)          :: mesh
!           ---------------------------------------------
            integer                :: eID
            integer                :: edID
            integer                :: nodesID(EDGES_PER_QUAD,QUADS_PER_EDGE)
            real(kind=RP)          :: xL(NDIM) , xR(NDIM)
            real(kind=RP)          :: dx(NDIM)
            class(Edge_t), pointer :: edge
            class(QuadElement_t), pointer       :: eL , eR

            nodesID ( EBOTTOM , : )  = [1 , 2]
            nodesID ( ERIGHT  , : )  = [2 , 3]
            nodesID ( ETOP    , : )  = [4 , 3]
            nodesID ( ELEFT   , : )  = [1 , 4]

            call SubSection_Header("Checking the Mesh")

!
!           Checking that in each mesh the element and edges connectivities are correct
!           ---------------------------------------------------------------------------
            do eID = 1 , mesh % no_of_elements

               do edID = 1 , EDGES_PER_QUAD
                  if ( mesh % elements(eID) % edgesDirection(edID) .eq. FORWARD ) then
                     if ( mesh % elements(eID) % nodes(nodesID(edID,1)) % n % ID .ne. mesh % elements(eID) % edges(edID) % f % nodes(1) % n % ID ) then
                        print*, "Connectivities are not correct."
                        stop "Stopped"
                     elseif ( mesh % elements(eID) % nodes(nodesID(edID,2)) % n % ID .ne. mesh % elements(eID) % edges(edID) % f % nodes(2) % n % ID ) then
                        print*, "Connectivities are not correct."
                        stop "Stopped"

                     end if

                  elseif ( mesh % elements(eID) % edgesDirection(edID) .eq. BACKWARD ) then
                     if ( mesh % elements(eID) % nodes(nodesID(edID,1)) % n % ID .ne. mesh % elements(eID) % edges(edID) % f % nodes(2) % n % ID ) then
                        print*, "Connectivities are not correct."
                        stop "Stopped"
                     elseif ( mesh % elements(eID) % nodes(nodesID(edID,2)) % n % ID .ne. mesh % elements(eID) % edges(edID) % f % nodes(1) % n % ID ) then
                        print*, "Connectivities are not correct."
                        stop "Stopped"

                     end if
               
                  end if

                  if ( minval( mesh % elements(eID) % jac ) .le. 0.0_RP ) then
                     print*, "Negative volume in element " , eID
                     stop "Stopped"
                  end if
               end do
            end do
!
!           Checking that the normal direction points towards the LEFT element, or towards the outside for boundaries
!           ---------------------------------------------------------------------------------------------------------
            do edID = 1 , mesh % no_of_edges

               edge => mesh % edges(edID) % f

               select type ( edge )
   
                  type is (Edge_t)

                     eL => mesh % edges(edID) % f % quads(LEFT) % e
                     eR => mesh % edges(edID) % f % quads(RIGHT) % e

                     xL = sum( sum( eL % X , dim = 1) , dim = 1 ) / ( eL % spA % N + 1.0_RP) ** 2.0_RP
                     xR = sum( sum( eR % X , dim = 1) , dim = 1 ) / ( eR % spA % N + 1.0_RP) ** 2.0_RP

                     if ( dot_product( xR-xL , edge % n(IX:IY,0) ) .lt. 0.0_RP  ) then
                        print*, "Edges normal orientation is not correct. "
                        stop "Stopped"

                     elseif ( edge % dS(0) .lt. 0.0_RP ) then
                        print*, "Negative surface in edge " , edID

                     end if
                  
                  type is (StraightBdryEdge_t)
      
                     eL => mesh % edges(edID) % f % quads(1) % e
                     eR => NULL()

                     xL = sum( sum( eL % X , dim = 1) , dim = 1 ) / ( eL % spA % N + 1.0_RP) ** 2.0_RP
                     xR = sum( edge % X , dim = 2 ) / (edge % spA % N + 1.0_RP)

                     if ( dot_product( xR - xL , edge % n(IX:IY,0) ) .lt. 0.0_RP ) then
                        print*, "Edges normal orientation is not correct. "
                        stop "Stopped"
   
                     elseif ( edge % dS(0) .lt. 0.0_RP ) then
                        print*, "Negative surface in edge " , edID

                     end if
               
                  type is (CurvedBdryEdge_t)
      
                     eL => mesh % edges(edID) % f % quads(1) % e
                     eR => NULL()

                     xL = sum( sum( eL % X , dim = 1) , dim = 1 ) / ( eL % spA % N + 1.0_RP) ** 2.0_RP
                     xR = sum( edge % X , dim = 2 ) / (edge % spA % N + 1.0_RP)

                     if ( dot_product( xR - xL , sum( edge % n , dim = 2) / (edge % spA % N + 1.0_RP)  ) .lt. 0.0_RP ) then
                        print*, "Edges normal orientation is not correct. "
                        stop "Stopped"
   
                     elseif ( minval(edge % dS) .lt. 0.0_RP ) then
                        print*, "Negative surface in edge " , edID

                     end if
               end select
            end do
!
!           Checking that the tangent vector follows the direction of the element
!           ---------------------------------------------------------------------
            do eID = 1 , mesh % no_of_elements
   
               eL => mesh % elements(eID)
               do edID = 1 , EDGES_PER_QUAD

                  dx = eL % nodes(nodesID(edID,2)) % n % x - eL % nodes(nodesID(edID,1)) % n % x 

                  if (( eL % edgesDirection(edID) .eq. FORWARD ) .and. ( dot_product( dx , eL % edges(edID) % f % dX(IX:IY,0) ) .le. 0.0_RP )) then
                     print*, "Edges tangent vector orientation is not correct."
                     stop "Stopped"

                  elseif (( eL % edgesDirection(edID) .eq. BACKWARD) .and. ( dot_product( dx , eL % edges(edID) % f % dX(IX:IY,0) ) .ge. 0.0_RP) ) then
                     print*, "Edges tangent vector orientation is not correct."
                     stop "Stopped"

                  end if
            
               end do
            end do
            write(STD_OUT , '(30X,A,A)') "-> ", "All tests succeeded."

        end subroutine CheckMesh
      
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
            real(kind=RP)              :: error , localerror
            integer                    :: current , location
            integer                    :: zone
            real(kind=RP), allocatable :: dSx(:) , dSy(:)
            real(kind=RP), allocatable :: dSe(:,:)
            real(kind=RP)              :: dX(NDIM,NDIM)
            integer                    :: i , j 

            call SubSection_Header("Testing the mappings")
            
!           This is to test the elements mappings derivatives formula
!           ---------------------------------------------------------
            error = 0.0_RP
   
            do eID = 1 , mesh % no_of_elements
               associate (e => mesh % elements(eID) )

               if ( allocated (dxiX ) ) deallocate ( dxiX  ) 
               if ( allocated (detaX) ) deallocate ( detaX ) 
          
               allocate(dxiX  ( 0 : e % spA % N , 0 : e % spA % N , NDIM) ) 
               allocate(detaX ( 0 : e % spA % N , 0 : e % spA % N , NDIM) ) 

               associate ( N => e % spA % N )
               dxiX         = MatrixMultiplyInIndex_F ( e % X , e % spA % DT , N+1 , N+1 , NDIM , IX ) 
               detaX        = MatrixMultiplyInIndex_F ( e % X , e % spA % DT , N+1 , N+1 , NDIM , IY ) 

               do i = 0 , N
                  do j = 0 , N
!
!                    Recover the mapping derivatives from the metric matrix
!                    ------------------------------------------------------   
                     dX(IX,IX) = e % Ja(i,j,IY,IY)
                     dX(IX,IY) = -e % Ja(i,j,IY,IX)
                     dX(IY,IX) = -e % Ja(i,j,IX,IY)
                     dX(IY,IY) = e % Ja(i,j,IX,IX)

                     localerror = maxval(abs(dxiX(i,j,1:NDIM) - dX(1:NDIM,IX)))
      
                     if ( localerror .gt. error ) then
                        error = localerror
                        current = eID
                     end if

                     localerror = maxval(abs(detaX(i,j,1:NDIM) - dX(1:NDIM,IY)))
                           
                     if ( localerror .gt. error ) then
                        error = localerror
                        current = eID
                     end if

                  end do
               end do 
               end associate
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
               dSx = -MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,1,2) , e % spA % lj(0.0_RP) , e % spA % N + 1 )
               dSy = -MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,2,2) , e % spA % lj(0.0_RP) , e % spA % N + 1 )

               if (( e % quadPosition(EBOTTOM) .eq. LEFT )) then
                  select type (ed => e % edges(EBOTTOM) % f)
                     type is (Edge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               elseif (( e % quadPosition(EBOTTOM) .eq. RIGHT) .and. (e % edgesDirection(EBOTTOM) .eq. FORWARD) ) then
                  select type (ed => e % edges(EBOTTOM) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               else
                  select type (ed => e % edges(EBOTTOM) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,e % spA % N : 0 : -1)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,e % spA % N : 0 : -1)
                  end select

               end if
               
               if ( maxval(abs(dSx - dSe(iX,0:e % spA % N) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,0:e % spA % N) ) )
                  current = eID
                  location = EBOTTOM
               end if
               if (  maxval(abs(dSy - dSe(iY,0:e % spA % N) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,0:e % spA % N) ) )
                  current = eID 
                  location = EBOTTOM
               end if

!              RIGHT Edge
!              -----------
               dSx = MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,1,1) , e % spA % lj(1.0_RP) , e % spA % N + 1 , trA = .true.)
               dSy = MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,2,1) , e % spA % lj(1.0_RP) , e % spA % N + 1 , trA = .true.)
               
               if (( e % quadPosition(ERIGHT) .eq. LEFT )) then
                  select type (ed => e % edges(ERIGHT) % f)
                     type is (Edge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               elseif (( e % quadPosition(ERIGHT) .eq. RIGHT) .and. (e % edgesDirection(ERIGHT) .eq. FORWARD) ) then
                  select type (ed => e % edges(ERIGHT) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               else
                  select type (ed => e % edges(ERIGHT) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,e % spA % N : 0 : -1)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,e % spA % N : 0 : -1)
                  end select

               end if
 
               if ( maxval(abs(dSx - dSe(iX,0:e % spA % N) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,0:e % spA % N) ) )
                  current = eID
                  location = ERIGHT
               end if
               if (  maxval(abs(dSy - dSe(iY,0:e % spA % N) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,:) ) )
                  current = eID 
                  location = ERIGHT
               end if              

!              TOP Edge
!              -----------
               dSx = MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,1,2) , e % spA % lj(1.0_RP) , e % spA % N + 1 )
               dSy = MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,2,2) , e % spA % lj(1.0_RP) , e % spA % N + 1 )

               if (( e % quadPosition(ETOP) .eq. LEFT )) then
                  select type (ed => e % edges(ETOP) % f)
                     type is (Edge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               elseif (( e % quadPosition(ETOP) .eq. RIGHT) .and. (e % edgesDirection(ETOP) .eq. FORWARD) ) then
                  select type (ed => e % edges(ETOP) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               else
                  select type (ed => e % edges(ETOP) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,e % spA % N : 0 : -1)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,e % spA % N : 0 : -1)
                  end select

               end if
 
               if ( maxval(abs(dSx - dSe(iX,:) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(iX,:) ) )
                  current = eID
                  location = ETOP
               end if
               if (  maxval(abs(dSy - dSe(iY,:) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,:) ) )
                  current = eID 
                  location = ETOP
               end if              

!              LEFT Edge
!              -----------
               dSx = -MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,1,1) , e % spA % lj(0.0_RP) , e % spA % N + 1 , trA = .true.)
               dSy = -MatrixTimesVector_F( e % Ja(0:e % spA % N,0:e % spA % N,2,1) , e % spA % lj(0.0_RP) , e % spA % N + 1 , trA = .true.)

               if (( e % quadPosition(ELEFT) .eq. LEFT )) then
                  select type (ed => e % edges(ELEFT) % f)
                     type is (Edge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               elseif (( e % quadPosition(ELEFT) .eq. RIGHT) .and. (e % edgesDirection(ELEFT) .eq. FORWARD) ) then
                  select type (ed => e % edges(ELEFT) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,0 : e % spA % N)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,0 : e % spA % N)
                  end select

               else
                  select type (ed => e % edges(ELEFT) % f)
                     type is (Edge_t)
                        dSe = -spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (StraightBdryEdge_t)
                        dSe = spread( ed % dS(0) * ed % n(IX:IY,0) , ncopies = e % spA % N + 1 , dim = 2 ) 
                     type is (CurvedBdryEdge_t)
                        dSe(IX,0 : e % spA % N) = ed % dS * ed % n(IX,e % spA % N : 0 : -1)
                        dSe(IY,0 : e % spA % N) = ed % dS * ed % n(IY,e % spA % N : 0 : -1)
                  end select

               end if
               
               if ( maxval(abs(dSx - dSe(IX,:) )) .gt. error ) then
                  error =  maxval(abs(dSx - dSe(IX,:) )) 
                  current = eID
                  location = ELEFT
               end if
               if (  maxval(abs(dSy - dSe(iY,:) ))  .gt. error ) then
                  error =  maxval(abs(dSy - dSe(iY,:) ) )
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
               write(STD_OUT,'(30X,A,A35,F16.10,A)') "-> ", "Computed surface in zone " // trim(mesh % Zones(zone) % Name) // ": ",mesh % ScalarScalarSurfaceIntegral("Surface",zone) * RefValues % L ,"." 
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
          integer                      :: direction
 
          write(STD_OUT,'(/)')
          call SubSection_Header("Checking the interpolation to boundaries")

          call DGSpatial_interpolateSolutionToBoundaries( mesh )

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

                  direction = mesh % edges(edID) % f % quads(quad) % e % edgesDirection( mesh % edges(edID) % f % edgeLocation(quad) )

                  if ( mesh % edges(edID) % f % transform(quad) ) then 
                     if ( direction .eq. FORWARD ) then
                        currentError = norm2( matmul( mesh % edges(edID) % f % T_forward(iXi,0:mesh % edges(edID) % f % NLow) , &
                                                      mesh % edges(edID) % f % storage(quad) % Q(0:,1:NCONS)) - mesh % IC(mesh % edges(edID) % f % x(:,iXi) ) )
                     else
                        currentError = norm2( matmul( mesh % edges(edID) % f % T_forward(iXi,0:mesh % edges(edID) % f % NLow) , mesh % edges(edID) % f % storage(quad) % Q(0:,1:NCONS)) - mesh % IC(mesh % edges(edID) % f % x(:,mesh % edges(edID) % f % spA % N - iXi) ) )
                     end if
                  else
                     if ( direction .eq. FORWARD ) then
                        currentError = norm2( mesh % edges(edID) % f % storage(quad) % Q(iXi,1:NCONS) - mesh % IC(mesh % edges(edID) % f % x(:,iXi) ) )
                     else
                        currentError = norm2( mesh % edges(edID) % f % storage(quad) % Q(iXi,1:NCONS) - mesh % IC(mesh % edges(edID) % f % x(:,mesh % edges(edID) % f % spA % N - iXi) ) )
                     end if

                  end if

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
         real(kind=RP)              :: error = 0.0_RP , currenterror = 0.0_RP
         real(kind=RP), allocatable :: Ja1(:,:) , Ja2(:,:)
         real(kind=RP), allocatable :: metricID(:,:)
         integer                    :: elem = 0

     
         write(STD_OUT,'(/)')
         call Subsection_Header("Checking discrete metric identities")

         do eID = 1 , mesh % no_of_elements

            associate( e => mesh % elements(eID) )

            allocate( Ja1(0:e % spA % N , 0 : e % spA % N ) )
            allocate( Ja2(0:e % spA % N , 0 : e % spA % N ) )
            allocate( metricID(0:e % spA % N , 0 : e % spA % N ) )
            do coord = 1 , NDIM
               
               Ja1 = e % Ja(0:e % spA % N,0:e % spA % N,coord,1) 
               Ja2 = e % Ja(0:e % spA % N,0:e % spA % N,coord,2) 

               associate ( N => e % spA % N )
               metricID = Mat_x_Mat_F( A = e % spA % D , B = Ja1 , rowC = N+1, colC = N+1 ) + Mat_x_Mat_F( A = Ja2 , B = e % spA % DT , rowC = N+1 , colC = N+1)
               end associate
               
               currenterror = abs(BilinearForm_F( A = metricID , X = e % spA % w , Y = e % spA % w ))
               
               if ( currenterror .gt. error ) then
                  error = currenterror
                  elem = e % ID
               end if

            end do

            deallocate(Ja1 , Ja2 , metricID)

            end associate
         end do

         write(STD_OUT , '(30X,A,A50,F16.10,A,I0,A)') "-> ", "Maximum discrete metric identities residual: " , error," (element " , elem, "). "

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
#ifdef NAVIER_STOKES   
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
#endif
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           CHECK AUXILIAR ROUTINES
!           -----------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#ifdef NAVIER_STOKES
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
#endif
        function QDotTrigonometricFCN( x , L ) result( val )
         use Physics
         implicit none
         real(kind=RP)           :: x(NDIM)
         real(kind=RP)           :: L 
         real(kind=RP)           :: val(NCONS)
         real(kind=RP)           :: u , v , p
         real(kind=RP)           :: ux , vy , H , uy , vx , px , py , Hx , Hy
         real(kind=RP)           :: tauxx , tauxy , tauyy , tauxx_x , tauyy_y
         real(kind=RP)           :: T_xx  , T_yy


         associate( gamma => Thermodynamics % gamma , Mach => dimensionless % Mach , cp => Dimensionless % cp)

#ifdef _DIMENSIONLESS_TAU
         u = sqrt(gamma) * Mach * sin(PI * x(IX) / L ) * cos(PI * x(IY) / L)
         v = -sqrt(gamma) * Mach * cos(PI * x(IX) / L ) * sin(PI * x(IY) / L)
         p = 1.0_RP + gamma * Mach * Mach * (1.0_RP / 8.0_RP) * ( cos(2.0_RP * PI * x(IX) / L ) + cos(2.0_RP * PI * x(IY) / L) )
         H = cp * p + 0.5_RP * u**2.0_RP + 0.5_RP * v**2.0_RP


         ux = sqrt(gamma) * Mach * PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L
         uy = -sqrt(gamma) * Mach * PI * sin(PI * x(IX) / L) * sin(PI * x(IY) / L ) / L

         vx = sqrt(gamma) * Mach * PI * sin(PI * x(IX) / L) * sin(PI * x(IY) / L ) / L 
         vy = -sqrt(gamma) * Mach * PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L 

         px = -0.25_RP * gamma * Mach * Mach * PI * sin(2.0_RP * PI * x(IX) / L ) / L 
         py = -0.25_RP * gamma * Mach * Mach * PI * sin(2.0_RP * PI * x(IY) / L ) / L 
#else
         u = sin(PI * x(IX) / L ) * cos(PI * x(IY) / L)
         v = -cos(PI * x(IX) / L ) * sin(PI * x(IY) / L)
         p = 1.0_RP + (1.0_RP / 8.0_RP) * ( cos(2.0_RP * PI * x(IX) / L ) + cos(2.0_RP * PI * x(IY) / L) )
         H = cp * p + 0.5_RP * u**2.0_RP + 0.5_RP * v**2.0_RP


         ux = PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L
         uy = -PI * sin(PI * x(IX) / L) * sin(PI * x(IY) / L ) / L

         vx = PI * sin(PI * x(IX) / L) * sin(PI * x(IY) / L ) / L 
         vy = -PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L 

         px = -0.25_RP * PI * sin(2.0_RP * PI * x(IX) / L ) / L 
         py = -0.25_RP * PI * sin(2.0_RP * PI * x(IY) / L ) / L 
#endif

         Hx = cp * px + u * ux + v * vx
         Hy = cp * py + u * uy + v * vy

         
   
         val(IRHO)      = 0.0_RP
         val(IRHOU)     = 2.0_RP * u * ux + u*vy + uy * v + px
         val(IRHOV)     = 2.0_RP * v * vy + py + u*vx + ux*v
         val(IRHOE) = H * (ux + vy) + Hx * u + Hy * v

#ifdef _DIMENSIONLESS_TAU
         val = -val * dimensionless % invSqrtGammaMach
#else
         val = -val
#endif
         
#ifdef NAVIER_STOKES
         associate ( mu => dimensionless % mu , kappa => dimensionless % kappa )

         tauxx = 2.0_RP * mu *  sqrt(gamma) * Mach * PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L
         tauyy = -2.0_RP * mu *  sqrt(gamma) * Mach * PI * cos(PI * x(IX) / L ) * cos(PI * x(IY) / L) / L
         tauxy = 0.0_RP

         tauxx_x = -2.0_RP * mu * sqrt(gamma) * Mach * PI * PI * sin(PI * x(IX) / L) * cos(PI * x(IY) / L ) / ( L * L )
         tauyy_y = 2.0_RP * mu * sqrt(gamma) * Mach * PI * PI * cos(PI * x(IX) / L) * sin(PI * x(IY) / L ) / ( L * L )

         T_xx = -0.5_RP * gamma * Mach * Mach * PI * PI * cos(2.0_RP * PI * x(IX) / L ) / (L * L)
         T_yy = -0.5_RP * gamma * Mach * Mach * PI * PI * cos(2.0_RP * PI * x(IY) / L ) / (L * L)

         val(IRHOU) = val(IRHOU) + tauxx_x
         val(IRHOV) = val(IRHOV) + tauyy_y
         val(IRHOE) = val(IRHOE) + ux * tauxx + u * tauxx_x + vy * tauyy + v*tauyy_y + kappa * T_xx + kappa * T_yy


         end associate
#endif
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

#ifdef _DIMENSIONLESS_TAU
         u = sqrt(gamma) * Mach * x(IX) / L
         v = sqrt(gamma) * Mach * x(IY) / L
         p = 1.0_RP + gamma * Mach * Mach * ( x(IX) * x(IX) + x(IY) * x(IY)) / (L **2.0_RP)
         H = cp * p + 0.5_RP * u**2.0_RP + 0.5_RP * v**2.0_RP

         ux = sqrt(gamma) * Mach * 1.0_RP/ L
         uy = 0.0_RP

         vx = 0.0_RP
         vy = sqrt(gamma) * Mach * 1.0_RP / L 

         px = 2.0_RP * gamma * Mach * Mach * x(IX) / L**2.0_RP
         py = 2.0_RP * gamma * Mach * Mach * x(IY) / L**2.0_RP
#else
         u = x(IX) / L
         v = x(IY) / L
         p = 1.0_RP + ( x(IX) * x(IX) + x(IY) * x(IY)) / (L **2.0_RP)
         H = cp * p + 0.5_RP * u**2.0_RP + 0.5_RP * v**2.0_RP

         ux = 1.0_RP/ L
         uy = 0.0_RP

         vx = 0.0_RP
         vy = 1.0_RP / L 

         px = 2.0_RP * x(IX) / L**2.0_RP
         py = 2.0_RP * x(IY) / L**2.0_RP
#endif

         Hx = cp * px + u * ux + v * vx
         Hy = cp * py + u * uy + v * vy

         val(IRHO)      = ux + vy
         val(IRHOU)     = 2.0_RP * u * ux + u*vy + uy * v + px
         val(IRHOV)     = 2.0_RP * v * vy + py + u*vx + ux*v
         val(IRHOE)     = H * (ux + vy) + Hx * u + Hy * v

         end associate
#ifdef NAVIER_STOKES
         associate ( gamma => thermodynamics % gamma , Mach => dimensionless % Mach , mu => dimensionless % mu , kappa => dimensionless % kappa )

         val(IRHOE) = val(IRHOE) + gamma * Mach * Mach / (L*L) * ( 4.0_RP/3.0_RP * mu + 4.0_RP * kappa ) 


         end associate
#endif

#ifdef _DIMENSIONLESS_TAU
         val = -val * dimensionless % invSqrtGammaMach
#else
         val = -val
#endif

        end function QDotPolynomicFCN        
      
end module
