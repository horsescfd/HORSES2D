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
         real(kind=RP),    allocatable,    intent(inout) :: var(:)
!        -------------------------------------------------------------------------
         integer                                       :: ncid
         integer                                       :: nDim , nVar
         integer                                       :: varID
         integer, dimension(2)                         :: dimids
         integer                                       :: varType , varNDim
         integer, dimension(1)                         :: dims
         character(len=STR_LEN_NETCDF)                 :: varName , dimName

!        Open file
         call check ( NF90_OPEN( trim(file) , NF90_NOWRITE , ncid ) )

!        Inquire number of variables and dimensions
         call check ( NF90_INQUIRE( ncid , ndim , nvar ) )

!
!        **************************************
!           Gather the variable dimensions
!        **************************************
!
         do varID = 1 , nVar
            call check ( NF90_INQUIRE_VARIABLE ( ncid , varID , varName , varType , varNDim , dimids ) ) 

            if (trim(varName) .eq. trim(name) ) then        ! Is the desired variable

!              Obtain dimensions
               call check ( NF90_INQUIRE_DIMENSION( ncid , dimids(1) , dimName ,  dims(1) ) )
              
!              If allocated, check if dimensions are consistent
               if (allocated(var)) then

                  if (dims(1) .eq. size(var,1) ) then
!                    Read the variable
                     call check ( NF90_GET_VAR( ncid , varID , var ) )
                  else
                     print*, "Preallocated matrix does not match the specified size"
                     stop "Stopped."
                  end if

               else
!                 Allocate and read
                  allocate( var ( dims(1) ) )
                  call check ( NF90_GET_VAR( ncid , varID , var ) )

               end if

               return

             end if

         end do

      end subroutine getDouble1DVariable
         
      subroutine getDouble2DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                  :: file
         character(len=*), intent(in)                  :: name
         real(kind=RP),    allocatable,    intent(inout) :: var(:,:)
!        -------------------------------------------------------------------------
         integer                                       :: ncid
         integer                                       :: nDim , nVar
         integer                                       :: varID
         integer, dimension(2)                         :: dimids
         integer                                       :: varType , varNDim
         integer, dimension(2)                         :: dims
         character(len=STR_LEN_NETCDF)                 :: varName , dimName
         real(kind=RP), allocatable                          :: aux(:,:)

!        Open file
         call check ( NF90_OPEN( trim(file) , NF90_NOWRITE , ncid ) )

!        Inquire number of variables and dimensions
         call check ( NF90_INQUIRE( ncid , ndim , nvar ) )

!
!        **************************************
!           Gather the variable dimensions
!        **************************************
!
         do varID = 1 , nVar
            call check ( NF90_INQUIRE_VARIABLE ( ncid , varID , varName , varType , varNDim , dimids ) ) 

            if (trim(varName) .eq. trim(name) ) then        ! Is the desired variable

!              Obtain dimensions
               call check ( NF90_INQUIRE_DIMENSION( ncid , dimids(1) , dimName ,  dims(1) ) )
               call check ( NF90_INQUIRE_DIMENSION( ncid , dimids(2) , dimName ,  dims(2) ) )
              
!              If allocated, check if dimensions are consistent
               if (allocated(var)) then

                  if ((dims(1) .eq. size(var,1)).and.(dims(2) .eq. size(var,2))) then
!                    Read the variable
                     call check ( NF90_GET_VAR( ncid , varID , var ) )
                  elseif ((dims(1) .eq. size(var,2)) .and. (dims(2) .eq. size(var,1)))  then
!                    Read the variable on an auxiliary variable, and transpose
                     print*, "Warning: dimensions are incorrect"
                     allocate( aux( dims(1) , dims(2) ) ) 
                     call check ( NF90_GET_VAR( ncid , varID , aux ) )
                     var = transpose(aux)
                     deallocate(aux)
                  else
                     print*, "Preallocated matrix does not match the specified size"
                     stop "Stopped."
                  end if
      
               else
!                 Allocate and read
                  allocate( var ( dims(1) , dims(2) ) )
                  call check ( NF90_GET_VAR( ncid , varID , var ) )

               end if

               return

             end if

         end do

      end subroutine getDouble2DVariable

      subroutine getInteger1DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                    :: file
         character(len=*), intent(in)                    :: name
         integer,          allocatable,    intent(inout) :: var(:)
!        -------------------------------------------------------------------------
         integer                                       :: ncid
         integer                                       :: nDim , nVar
         integer                                       :: varID
         integer, dimension(2)                         :: dimids
         integer                                       :: varType , varNDim
         integer, dimension(1)                         :: dims
         character(len=STR_LEN_NETCDF)                 :: varName , dimName
!        Open file
         call check ( NF90_OPEN( trim(file) , NF90_NOWRITE , ncid ) )

!        Inquire number of variables and dimensions
         call check ( NF90_INQUIRE( ncid , ndim , nvar ) )

!
!        **************************************
!           Gather the variable dimensions
!        **************************************
!
         do varID = 1 , nVar
            call check ( NF90_INQUIRE_VARIABLE ( ncid , varID , varName , varType , varNDim , dimids ) ) 

            if (trim(varName) .eq. trim(name) ) then        ! Is the desired variable

!              Obtain dimensions
               call check ( NF90_INQUIRE_DIMENSION( ncid , dimids(1) , dimName ,  dims(1) ) )
              
!              If allocated, check if dimensions are consistent
               if (allocated(var)) then

                  if (dims(1) .eq. size(var,1) ) then
!                    Read the variable
                     call check ( NF90_GET_VAR( ncid , varID , var ) )
                  else
                     print*, "Preallocated matrix does not match the specified size"
                     stop "Stopped."
                  end if

               else
!                 Allocate and read
                  allocate( var ( dims(1) ) )
                  call check ( NF90_GET_VAR( ncid , varID , var ) )

               end if

               return

             end if

         end do

      end subroutine getInteger1DVariable

      subroutine getInteger2DVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                  :: file
         character(len=*), intent(in)                  :: name
         integer      ,    allocatable,    intent(inout) :: var(:,:)
!        -------------------------------------------------------------------------
         integer                                       :: ncid
         integer                                       :: nDim , nVar
         integer                                       :: varID
         integer, dimension(2)                         :: dimids
         integer                                       :: varType , varNDim
         integer, dimension(2)                         :: dims
         character(len=STR_LEN_NETCDF)                 :: varName , dimName
         integer, allocatable                          :: aux(:,:)

!        Open file
         call check ( NF90_OPEN( trim(file) , NF90_NOWRITE , ncid ) )

!        Inquire number of variables and dimensions
         call check ( NF90_INQUIRE( ncid , ndim , nvar ) )

!
!        **************************************
!           Gather the variable dimensions
!        **************************************
!
         do varID = 1 , nVar
            call check ( NF90_INQUIRE_VARIABLE ( ncid , varID , varName , varType , varNDim , dimids ) ) 

            if (trim(varName) .eq. trim(name) ) then        ! Is the desired variable

!              Obtain dimensions
               call check ( NF90_INQUIRE_DIMENSION( ncid , dimids(1) , dimName ,  dims(1) ) )
               call check ( NF90_INQUIRE_DIMENSION( ncid , dimids(2) , dimName ,  dims(2) ) )
              
!              If allocated, check if dimensions are consistent
               if (allocated(var)) then

                  if ((dims(1) .eq. size(var,1)).and.(dims(2) .eq. size(var,2))) then
!                    Read the variable
                     call check ( NF90_GET_VAR( ncid , varID , var ) )
                  elseif ((dims(1) .eq. size(var,2)) .and. (dims(2) .eq. size(var,1)))  then
!                    Read the variable on an auxiliary variable, and transpose
                     allocate( aux( dims(1) , dims(2) ) ) 
                     call check ( NF90_GET_VAR( ncid , varID , aux ) )
                     var = transpose(aux)
                     deallocate(aux)
                  else
                     print*, "Preallocated matrix does not match the specified size"
                     stop "Stopped."
                  end if
      
               else
!                 Allocate and read
                  allocate( var ( dims(1) , dims(2) ) )
                  call check ( NF90_GET_VAR( ncid , varID , var ) )

               end if

               return

             end if

         end do

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
