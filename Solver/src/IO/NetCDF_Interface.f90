module NetCDFInterface
   use SMConstants
   use NetCDF
   
   private
   public NetCDF_getDimension , NetCDF_getVariable


   integer, parameter            :: STR_LEN_NETCDF    = 128

   interface NetCDF_getDimension
      module procedure getDimension
   end interface NetCDF_getDimension
   
   interface NetCDF_getVariable
      module procedure getDouble1DVariable , getInteger1DVariable , getDouble2DVariable , getInteger2DVariable
   end interface NetCDF_getVariable

   contains
      function getDimension(file , name) result(value)
         implicit none
         character(len=*), intent(in)                 :: file
         character(len=*), intent(in)                 :: name
         integer                                      :: value
!        ------------------------------------------------------------------
         integer                                      :: ncid
         integer                                      :: ndim , nvar
         integer                                      :: dim
         integer                                      :: currentVal
         character(len=STR_LEN_NETCDF)                :: currentName
         
!        Default value
         value = -1   

!        Open file
         call check ( NF90_OPEN( trim(file) , NF90_NOWRITE , ncid ) )

!        Inquire number of variables and dimensions
         call check ( NF90_INQUIRE( ncid , ndim , nvar ) )

!        Search for the desired dimension
         do dim = 1 , ndim
            call check ( NF90_INQUIRE_DIMENSION( ncid , dim , currentName , currentVal ) )
      
            if (trim(currentName) .eq. trim(name)) then
               value = currentVal
               return
            end if

         end do

!        Close file
         call check ( NF90_CLOSE( ncid ) )

      end function getDimension
      subroutine getDouble1DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                  :: file
         character(len=*), intent(in)                  :: name
         real(kind=RP),    allocatable,    intent(out) :: var(:)

      end subroutine getDouble1DVariable
         
      subroutine getDouble2DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                  :: file
         character(len=*), intent(in)                  :: name
         real(kind=RP),    allocatable,    intent(out) :: var(:,:)

      end subroutine getDouble2DVariable

      subroutine getInteger1DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                  :: file
         character(len=*), intent(in)                  :: name
         integer,          allocatable,    intent(out) :: var(:)

      end subroutine getInteger1DVariable

      subroutine getInteger2DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                  :: file
         character(len=*), intent(in)                  :: name
         integer      ,    allocatable,    intent(out) :: var(:,:)

      end subroutine getInteger2DVariable
!
!     *************************************
!        Check subroutine
!     *************************************
!
      subroutine check( status )
         integer, intent(in)        :: status
         
         if (status .ne. NF90_NOERR) then
            print*, trim(NF90_STRERROR( status ) ) 
            stop "Stopped."
         end if
      end subroutine check

end module NetCDFInterface
