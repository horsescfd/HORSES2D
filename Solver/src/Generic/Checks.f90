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
          !use Physics
          !use NodesAndWeights_class
          !use QuadMeshClass
          !use MeshFileClass
          !use DGSpatialDiscretizationMethods
          !use Storage_module
          !use DGBoundaryConditions  
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
      !    do eID = 1 , sem % mesh % no_of_elements
      !      write(STD_OUT , '(6F24.16)') sem % mesh % elements(eID) % x
      !    end do
       !           write(STD_OUT , *)  "IC"
       !            do eID = 1 , sem % mesh  % no_of_elements
                   !write(STD_OUT , '(6F24.16)') sem % mesh % elements(eID) % Q(:,1)
       !            end do
      
      
      !    call DGSpatial_interpolateToBoundaries( sem % mesh ,"Q")
      
      !    call DGSpatial_computeGradient( sem % mesh )
      
      !    call DGSpatial_computeTimeDerivative ( sem % mesh ) 
      
      !    print*, "array"
      !     write(STD_OUT , '(6F24.16)') sem % Storage % Q
      
      
      !    print*, "Lets check the interpolation to boundaries"
      !    do eID = 1 , sem % mesh % no_of_elements
      !       write(STD_OUT , '(2F24.16)') sem % mesh % elements(eID) % Qb
      !    enddo
      
      !            print*, "Checking boundary conditions"
      
      
       !           print*, "normal faces"
       !           do fID = 1 , sem % mesh % no_of_faces
       !              write(STD_OUT , '(F10.3)') sem % mesh % faces(fID) % f % n
       !           end do
               
       !           print*, "markers"
       !           do fID = 1 , sem % mesh % no_of_faces
       !              write(STD_OUT , '(I10)') sem % mesh % faces(fID) % f % faceType
       !           end do
         
       !           print*, "BCLocations"
       !           do fID = 1 , sem % mesh % no_of_faces
       !              select type (f=>sem % mesh % faces(fID) % f)
       !                 type is (BdryFace_t)
       !                    print*, "Boundary face no ", fID ,"."
       !                    print*, "      BCLocation: " , f % BCLocation
       !                    print*, "Boundary value: " , f % uB
       !              end select
       !           end do
                        
       !           print*, "Elements ID , LEFT Face , RIGHT Face"
       !           do fID = 1 , sem % mesh % no_of_elements
       !              print*, fID, sem % mesh % elements(fID) % facesID(LEFT) , sem % mesh % elements(fID) % facesID(RIGHT)
       !           end do
       !           print*, "Faces ID"
       !           do fID = 1 , sem % mesh % no_of_faces
       !              print*, sem % mesh % faces(fID) % f % ID
       !           end do
      
      
       !    print*, "Checking gradients"
       !    do eID = 1 , sem % mesh % no_of_elements
       !       write(STD_OUT , '(6F24.16)') sem % mesh % elements(eID) % dQ
       !    end do
      
      !      print*, "Computing QDot......."
      
      !     do eID = 1 , sem % mesh % no_of_elements
      !        write(STD_OUT , '(6F24.16)') sem % mesh % elements(eID) % QDot  
      !     end do
      
      !     call Integration_checks(sem)
      
      
            
         end if
        end subroutine checks
      
        subroutine CheckMappings( mesh )
            use QuadElementClass
            use Headers
            use QuadMeshClass
            use MatrixOperations
            implicit none
            type(QuadMesh_t)           :: mesh
!           ---------------------------------------
            integer                    :: eID
            real(kind=RP), allocatable :: dxiX(:,:,:)
            real(kind=RP), allocatable :: detaX(:,:,:)
            real(kind=RP)              :: error
            integer                    :: current

            call SubSection_Header("Testing the mappings")
            
            do eID = 1 , mesh % no_of_elements
               associate (e => mesh % elements(eID) )

               if ( allocated (dxiX ) ) deallocate ( dxiX  ) 
               if ( allocated (detaX) ) deallocate ( detaX ) 
          
               allocate(dxiX  ( NDIM , 0 : e % spA % N , 0 : e % spA % N ) ) 
               allocate(detaX ( NDIM , 0 : e % spA % N , 0 : e % spA % N ) ) 

               dxiX         = MatrixMultiplyInIndex_F( e % X , transpose( e % spA % D ) , 2)
               detaX        = MatrixMultiplyInIndex_F( e % X , transpose( e % spA % D ) , 3)

               if (eID .eq. 1) then
                  current = eID
                  error = maxval(abs(dxiX - e % dX(:,:,:,iX)))
                   
               else

                  if ( maxval(abs(dxiX - e % dX(:,:,:,iX) ) ) .gt. error ) then
                     error = maxval(abs(dxiX - e % dX(:,:,:,iX) ) )
                     current = eID
                  end if
               end if
                     
               if ( maxval(abs(detaX - e % dX(:,:,:,iY)) ) .gt. error) then
                     error = maxval(abs(detaX - e % dX(:,:,:,iY) ) )
                     current = eID
               end if
               end associate

            end do

            write(STD_OUT , '(30X,A,A,E10.3,A,I0,A)') "-> ", "Maximum error found: ",error,"  (Cell ",current,")."

!           Compute the volume of the domain
            write(STD_OUT , '(30X,A,A,F16.10,A)') "-> ", "Domain volume: " , mesh % VolumeIntegral("One"),"."
        end subroutine CheckMappings
      
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
end module
