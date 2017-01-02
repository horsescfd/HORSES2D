module MatrixOperations
   use SMConstants   

   interface InnerProduct
      module procedure InnerProduct2D
   end interface

   contains
      function vectorOuterProduct(a,b) result(val)
         implicit none  
         real(kind=RP)              :: a(:)
         real(kind=RP)              :: b(:)
         real(kind=RP), allocatable :: val(:,:)
         integer                    :: N , M

         N = size(a,1)
         M = size(b,1)

         allocate(val(1:N , 1:M) )

         val = spread(a(1:N),dim=2,ncopies=M)*spread(b(1:M),dim=1,ncopies=N)

      end function vectorOuterProduct

!      function MatrixTimesVector_F( A , V ) result( C )
!!     -----------------------------
!!        Computes the product
!!           C = A * B
!!     -----------------------------
!         implicit none
!         real(kind=RP), intent(in)        :: A(:,:)
!         real(kind=RP), intent(in)        :: V(:,:)
!         real(kind=RP), allocatable       :: C(:,:)
!!
!!        ----------------------------------------
!!           Variables for lapack dgemm
!!        ----------------------------------------
!!
!#ifdef _USE_LAPACK  ! ---------------------------------------
!         integer                 :: N , M , K 
!         integer                 :: LDA , LDB , LDC
!#endif ! ----------------------------------------------------
!!
!
!         allocate(C(size(A,1) , size(B,2)) )
!
!#ifdef _USE_LAPACK  ! ------------------------------------------------------------------------------------------
!
!!           Set dimensions
!            M = size(A,1)
!            K = size(A,2)
!            N = size(B,2)
!            LDA = M
!            LDB = K
!            LDC = M
!
!            call dgemm( "N" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
!
!#else ! ------------------------------------------------------------------------------------------------------------
!            C = MATMUL(  A  , B ) 
!#endif ! -------------------------------------------------------------------------------------------------------------
!
!      end function MatrixTimesVector_F
!
      subroutine Mat_x_Mat( trA , trB , A , B , C )
!     -----------------------------
!        Computes the product
!           C = tr(A) * B
!     -----------------------------
         implicit none
         logical                          :: trA
         logical                          :: trB
         real(kind=RP), intent(in)        :: A(:,:)
         real(kind=RP), intent(in)        :: B(:,:)
         real(kind=RP), intent(out)       :: C(:,:)
!
!        ----------------------------------------
!           Variables for lapack dgemm
!        ----------------------------------------
!
#ifdef _USE_LAPACK  ! ---------------------------------------
         integer                 :: N , M , K 
         integer                 :: LDA , LDB , LDC
#endif ! ----------------------------------------------------
!

#ifdef _USE_LAPACK  ! ------------------------------------------------------------------------------------------

!           Set dimensions
            if (trA) then
               M = size(A,2)
            else 
               M = size(A,1)
            end if

            if (trB) then
               N = size(B,1)
            else
               N = size(B,2)
            end if

            if (trA) then
               K = size(A,1)
            else
               K = size(A,2)
            end if

            LDA = K
            LDB = K
            LDC = M

            if ( (.not. trA) .and. (.not. trB) ) then
               call dgemm( "N" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            elseif ( (.not. trA) .and. ( trB ) ) then
               call dgemm( "N" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            elseif ( ( trA ) .and. (.not. trB) ) then
               call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            elseif ( ( trA ) .and. ( trB ) ) then
               call dgemm( "T" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            end if
#else ! ------------------------------------------------------------------------------------------------------------
            C = MATMUL( TRANSPOSE( A ) , B ) 
            if ( (.not. trA) .and. (.not. trB) ) then
               C = matmul(A,B)
            elseif ( (.not. trA) .and. ( trB ) ) then
               C = matmul(A,transpose(B))
            elseif ( ( trA ) .and. (.not. trB) ) then
               call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
               C = matmul(transpose(A),B)
            elseif ( ( trA ) .and. ( trB ) ) then
               C = matmul(transpose(A),transpose(B))
            end if
#endif ! -------------------------------------------------------------------------------------------------------------

      end subroutine Mat_x_Mat

      function Mat_x_Mat_F( trA , trB , A , B ) result ( C )
!     -----------------------------
!        Computes the product
!           C = tr(A) * B
!     -----------------------------
         implicit none
         logical                          :: trA
         logical                          :: trB
         real(kind=RP), intent(in)        :: A(:,:)
         real(kind=RP), intent(in)        :: B(:,:)
         real(kind=RP), allocatable       :: C(:,:)
!
!        ----------------------------------------
!           Variables for lapack dgemm
!        ----------------------------------------
!
         integer                 :: N , M , K 
#ifdef _USE_LAPACK  ! ---------------------------------------
         integer                 :: LDA , LDB , LDC
#endif ! ----------------------------------------------------
!


!           Set dimensions
            if (trA) then
               M = size(A,2)
            else 
               M = size(A,1)
            end if

            if (trB) then
               N = size(B,1)
            else
               N = size(B,2)
            end if

            if (trA) then
               K = size(A,1)
            else
               K = size(A,2)
            end if


            allocate( C(M,N) )

#ifdef _USE_LAPACK  ! ------------------------------------------------------------------------------------------
            LDA = K
            LDB = K
            LDC = M

            if ( (.not. trA) .and. (.not. trB) ) then
               call dgemm( "N" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            elseif ( (.not. trA) .and. ( trB ) ) then
               call dgemm( "N" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            elseif ( ( trA ) .and. (.not. trB) ) then
               call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            elseif ( ( trA ) .and. ( trB ) ) then
               call dgemm( "T" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
            end if
#else ! ------------------------------------------------------------------------------------------------------------
            C = MATMUL( TRANSPOSE( A ) , B ) 
            if ( (.not. trA) .and. (.not. trB) ) then
               C = matmul(A,B)
            elseif ( (.not. trA) .and. ( trB ) ) then
               C = matmul(A,transpose(B))
            elseif ( ( trA ) .and. (.not. trB) ) then
               call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )
               C = matmul(transpose(A),B)
            elseif ( ( trA ) .and. ( trB ) ) then
               C = matmul(transpose(A),transpose(B))
            end if
#endif ! -------------------------------------------------------------------------------------------------------------

      end function Mat_x_Mat_F

      function MatrixMultiplyInIndex_F( A , B , index) result( C )
!     -----------------------------
!        Computes the product
!           C = A * B
!     -----------------------------
         implicit none
         real(kind=RP), intent(in)        :: A(:,:,:)
         real(kind=RP), intent(in)        :: B(:,:)
         integer                          :: index
         real(kind=RP), allocatable       :: C(:,:,:)
         integer                          :: I1 , I2 , I3
         integer                          :: i , j
         
         I1 = size(A,1)
         I2 = size(A,2)
         I3 = size(A,3)
         
         if (index .eq. 1) then
            I1 = size(B,2)
         elseif ( index .eq. 2) then
            I2 = size(B,2)
         elseif ( index .eq. 3) then
            I3 = size(B,2)
         end if

         
         allocate(C(I1,I2,I3))

         if (index .eq. 1) then
            do i = 1 , I2
               do j = 1 , I3
                  C(:,i,j) = matmul(A(:,i,j) , B) 
               end do
            end do
         
         elseif (index .eq. 2) then
            do i = 1 , I1
               do j = 1 , I3
                  C(i,:,j) = matmul(A(i,:,j) , B) 
               end do
            end do
         elseif (index .eq. 3) then
            do i = 1 , I1
               do j = 1 , I2
                  C(i,j,:) = matmul(A(i,j,:) , B) 
               end do
            end do
         end if

      end function MatrixMultiplyInIndex_F

      subroutine innerProduct2D( A , M , val )
!     
!     ******************************************************
!        Computes the product
!           val = tr(A) M A
!     ******************************************************
!
         implicit none
         real(kind=RP), intent(in)     :: A(:,:)
         real(kind=RP), intent(in)     :: M(:,:)
         real(kind=RP), intent(out)     :: val(:,:)

         val = matmul(matmul(transpose(A) , M ) , A)

      end subroutine innerProduct2D

      function inv(A) result(Ainv)
!
!--------------------------------------------------------------
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
!     Author: Gonzalo Rubio (g.rubio@upm.es)
!--------------------------------------------------------------

        real(KIND = RP), dimension(:,:), intent(in) :: A
        real(KIND = RP), dimension(size(A,1),size(A,2)) :: Ainv
      
        real(KIND = RP), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
      
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
      
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
      
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        !PRINT*, "Call DGETRF"
        call DGETRF(n, n, Ainv, n, ipiv, info)
        !PRINT*, "info", info
        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
      
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
      
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
end function inv
!
!////////////////////////////////////////////////////////////////////////
!


end module MatrixOperations
