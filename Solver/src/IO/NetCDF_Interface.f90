!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////
!
module NetCDFInterface
   use SMConstants
   use NetCDF
   
   private
   public NetCDF_CreateFile 
   public NetCDF_putDimension , NetCDF_getDimension 
   public NetCDF_putVariable  , NetCDF_getVariable


   integer, parameter            :: STR_LEN_NETCDF    = 128

   interface NetCDF_CreateFile
      module procedure CreateFile
   end interface NetCDF_CreateFile

   interface NetCDF_putDimension
      module procedure putDimension
   end interface NetCDF_putDimension

   interface NetCDF_getDimension
      module procedure getDimension
   end interface NetCDF_getDimension
   
   interface NetCDF_putVariable
      module procedure putDouble1DVariable , putInteger1DVariable !, putDouble2DVariable , putInteger2DVariable
!      module procedure putCharacterVariable
   end interface NetCDF_putVariable

   interface NetCDF_getVariable
      module procedure getDouble1DVariable , getInteger1DVariable , getDouble2DVariable , getInteger2DVariable
      module procedure getCharacterVariable
   end interface NetCDF_getVariable


   contains

      subroutine createFile( fileName ) 
         implicit none
         character(len=*), intent(in)                 :: fileName
!        -------------------------------------------------------
         integer                                      :: ncid

         call check ( NF90_CREATE ( trim(fileName) , NF90_CLOBBER , ncid ) ) 

         call check ( NF90_ENDDEF ( ncid ) )

         call check ( NF90_CLOSE ( ncid ) )

      end subroutine createFile
   
      subroutine putDimension( fileName , Name , Value )
         implicit none
         character(len=*), intent(in)                 :: fileName
         character(len=*), intent(in)                 :: Name
         integer, intent(in)                          :: Value
!        ---------------------------------------------------------------------
         integer                                      :: ncid
         integer                                      :: dimid

!
!        Open file
!        ---------
         call check ( NF90_OPEN ( trim(fileName) , NF90_WRITE , ncid ) ) 
!
!        Put it into define mode
!        -----------------------   
         call check ( NF90_REDEF ( ncid ) ) 
!
!        Insert the dimension
!        --------------------
         call check ( NF90_DEF_DIM ( ncid , trim(Name) , Value , dimid ) )
!
!        Put it back into variable mode
!        ------------------------------
         call check ( NF90_ENDDEF ( ncid ) )
!
!        Close file
!        ----------
         call check ( NF90_CLOSE( ncid ) )
          
      end subroutine putDimension

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

      subroutine putDouble1DVariable( fileName , varName , dimNames , var )
         implicit none
         character(len=*), intent(in)                 :: fileName
         character(len=*), intent(in)                 :: varName
         character(len=*), intent(in)                 :: dimNames(:)
         real(kind=RP), intent(in)                    :: var(:)
!        ----------------------------------------------------------------
         integer                                      :: ncid
         integer, allocatable                         :: dimids(:)
         integer                                      :: nDim , nVar
         integer                                      :: iDim
         integer                                      :: currentDim
         character(len=STR_LEN_NETCDF)                :: currentName
         integer                                      :: currentVal
         integer                                      :: varID
!
!        Open file
!        ---------
         call check ( NF90_OPEN ( trim(fileName) , NF90_WRITE , ncid ) )
!
!        Get number of variables and dimensions
!        --------------------------------------
         call check ( NF90_INQUIRE( ncid , nDim , nVar) ) 
!
!        Look for the variable dimensions
!        --------------------------------
         allocate( dimids( size(dimNames) ) )

         currentDim = 1
         do iDim = 1 , nDim
            call check ( NF90_INQUIRE_DIMENSION ( ncid , iDim , currentName , currentVal ) )
            
            if ( trim(currentName) .eq. trim(dimNames(currentDim)) )  then
!              
!              The dimension is found
!              ----------------------
               dimids(currentDim)   = iDim

               if ( currentDim .eq. size(dimNames) ) then
                  exit
               end if

               currentDim = currentDim + 1
        
            end if
         end do
!
!        Once the dimensions are found, put the variable
!        -----------------------------------------------
         call check ( NF90_REDEF ( ncid ) ) 

         call check ( NF90_DEF_VAR ( ncid , trim(varName) , NF90_DOUBLE , dimids , varID ) ) 

         call check ( NF90_ENDDEF ( ncid ) )

         call check ( NF90_PUT_VAR ( ncid , varID , var ) ) 

         call check ( NF90_CLOSE ( ncid ) )

      end subroutine

      subroutine putInteger1DVariable( fileName , varName , dimNames , var )
         implicit none
         character(len=*), intent(in)                 :: fileName
         character(len=*), intent(in)                 :: varName
         character(len=*), intent(in)                 :: dimNames(:)
         integer, intent(in)                          :: var(:)
!        ----------------------------------------------------------------
         integer                                      :: ncid
         integer, allocatable                         :: dimids(:)
         integer                                      :: nDim , nVar
         integer                                      :: iDim
         integer                                      :: currentDim
         character(len=STR_LEN_NETCDF)                :: currentName
         integer                                      :: currentVal
         integer                                      :: varID
!
!        Open file
!        ---------
         call check ( NF90_OPEN ( trim(fileName) , NF90_WRITE , ncid ) )
!
!        Get number of variables and dimensions
!        --------------------------------------
         call check ( NF90_INQUIRE( ncid , nDim , nVar) ) 
!
!        Look for the variable dimensions
!        --------------------------------
         allocate( dimids( size(dimNames) ) )

         currentDim = 1
         do iDim = 1 , nDim
            call check ( NF90_INQUIRE_DIMENSION ( ncid , iDim , currentName , currentVal ) )
            
            if ( trim(currentName) .eq. trim(dimNames(currentDim)) )  then
!              
!              The dimension is found
!              ----------------------
               dimids(currentDim)   = iDim

               if ( currentDim .eq. size(dimNames) ) then
                  exit
               end if

               currentDim = currentDim + 1
        
            end if
         end do
!
!        Once the dimensions are found, put the variable
!        -----------------------------------------------
         call check ( NF90_REDEF ( ncid ) ) 

         call check ( NF90_DEF_VAR ( ncid , trim(varName) , NF90_INT , dimids , varID ) ) 

         call check ( NF90_ENDDEF ( ncid ) )

         call check ( NF90_PUT_VAR ( ncid , varID , var ) ) 

         call check ( NF90_CLOSE ( ncid ) )

      end subroutine putInteger1DVariable

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

               call check ( NF90_CLOSE ( ncid ) )
               return

             end if

         end do
!
!        Close file
         call check ( NF90_CLOSE ( ncid ) )
         

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

               call check ( NF90_CLOSE ( ncid ) )
               return

             end if

         end do
!
!        Close file
         call check ( NF90_CLOSE ( ncid ) )

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

               call check ( NF90_CLOSE ( ncid ) )
               return

             end if

         end do
!
!        Close file
         call check ( NF90_CLOSE ( ncid ) )

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

               call check ( NF90_CLOSE ( ncid ) )
               return

             end if

         end do
!
!        Close file
         call check ( NF90_CLOSE ( ncid ) )


      end subroutine getInteger2DVariable

      subroutine getCharacterVariable( file , name , var )
         implicit none
         character(len=*), intent(in)                    :: file
         character(len=*), intent(in)                    :: name
         character(len=*), intent(inout)                 :: var
!        -------------------------------------------------------------------------
         integer                                       :: ncid
         integer                                       :: nDim , nVar
         integer                                       :: varID
         integer, dimension(2)                         :: dimids
         integer                                       :: varType , varNDim
         integer, dimension(1)                         :: dims
         character(len=STR_LEN_NETCDF)                 :: varName , dimName
!
!        Initialize variable
         var = ""

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
              
               call check ( NF90_GET_VAR( ncid , varID , var(1:dims(1)) ) )

               call check ( NF90_CLOSE ( ncid ) )
               return

             end if

         end do

         call check ( NF90_CLOSE ( ncid ) )

      end subroutine getCharacterVariable

 
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
