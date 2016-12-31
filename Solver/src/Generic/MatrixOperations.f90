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

      subroutine TransposeMat_x_NormalMat( A , B , C )
!     -----------------------------
!        Computes the product
!           C = tr(A) * B
!     -----------------------------
         implicit none
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
            M = size(A,2)
            K = size(A,1)
            N = size(B,2)
            LDA = K
            LDB = K
            LDC = M

            call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )

#else ! ------------------------------------------------------------------------------------------------------------
            C = MATMUL( TRANSPOSE( A ) , B ) 
#endif ! -------------------------------------------------------------------------------------------------------------

      end subroutine TransposeMat_x_NormalMat

      function TransposeMat_x_NormalMat_F( A , B ) result( C )
!     -----------------------------
!        Computes the product
!           C = tr(A) * B
!     -----------------------------
         implicit none
         real(kind=RP), intent(in)        :: A(:,:)
         real(kind=RP), intent(in)        :: B(:,:)
         real(kind=RP), allocatable       :: C(:,:)
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

         allocate(C(size(A,2) , size(B,2)) )

#ifdef _USE_LAPACK  ! ------------------------------------------------------------------------------------------

!           Set dimensions
            M = size(A,2)
            K = size(A,1)
            N = size(B,2)
            LDA = K
            LDB = K
            LDC = M

            call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )

#else ! ------------------------------------------------------------------------------------------------------------
            C = MATMUL( TRANSPOSE( A ) , B ) 
#endif ! -------------------------------------------------------------------------------------------------------------

      end function TransposeMat_x_NormalMat_F

      function NormalMat_x_TransposeMat_F( A , B ) result( C )
!     -----------------------------
!        Computes the product
!           C = A * tr(B)
!     -----------------------------
         implicit none
         real(kind=RP), intent(in)        :: A(:,:)
         real(kind=RP), intent(in)        :: B(:,:)
         real(kind=RP), allocatable       :: C(:,:)
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

         allocate(C(size(A,1) , size(B,1)) )

#ifdef _USE_LAPACK  ! ------------------------------------------------------------------------------------------

!           Set dimensions
            M = size(A,1)
            K = size(A,2)
            N = size(B,1)
            LDA = M
            LDB = N
            LDC = M

            call dgemm( "N" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , 0.0_RP , C , LDC )

#else ! ------------------------------------------------------------------------------------------------------------
            C = MATMUL(  A  , transpose(B) ) 
#endif ! -------------------------------------------------------------------------------------------------------------

      end function NormalMat_x_TransposeMat_F

      subroutine TripleMatrixProduct( A , B , C , val )
!
!        ***********************************
!           Computes the product 
!              val = A B C
!        ***********************************
!
         implicit none
         real(kind=RP), intent(in)        :: A(:,:)
         real(kind=RP), intent(in)        :: B(:,:)
         real(kind=RP), intent(in)        :: C(:,:)
         real(kind=RP), intent(out)       :: val(:,:)

         val = matmul(matmul(A,B),C)

      end subroutine TripleMatrixProduct

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
