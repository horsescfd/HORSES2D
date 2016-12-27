module FileWriting
   use SMConstants
   use NetCDF
   use Physics
   use Mesh1DClass
   use Storage_module
   use Setup_class

   integer, parameter         :: STR_FILENAME_LENGTH = 128

   private
   public   FileWriting_SaveSolution

   type NetCDF_File
      integer                            :: fID
      character(len=STR_FILENAME_LENGTH) :: fileName
      class(Dimension_t), pointer        :: dimension_head => NULL()
      integer                            :: ndim
      contains
         procedure  :: create          => FileWriting_createFile
         procedure  :: writeDimensions => FileWriting_writeDimensions
   end type NetCDF_File


   type Dimension_t
      integer                            :: value
      integer                            :: dimID
      character(len=STR_FILENAME_LENGTH) :: name
      class(Dimension_t), pointer        :: next => NULL()
   end type Dimension_t

   contains

      subroutine FileWriting_SaveSolution( mesh , fileName , t , Storage)
         use NetCDF
         implicit none
         class(Mesh1D_t)  :: mesh
         character(len=*) :: fileName
         real(kind=RP)    :: t
         class(Storage_t) :: Storage
         type(NetCDF_File)    :: file
      
         file % fileName = fileName
         

         call file % create() 
         call file % writeDimensions(mesh)
      
      end subroutine FileWriting_SaveSolution

      subroutine FileWriting_createFile( self ) 
         use NetCDF
         implicit none
         class(NetCDF_File)         :: self

         call check ( NF90_OPEN ( trim( self % fileName) , NF90_CLOBBER , self % fID ) )

      end subroutine FileWriting_createFile

      subroutine FileWriting_writeDimensions( self , mesh )
         use NetCDF
         use Physics
         implicit none
         class(NetCDF_File)                     :: self 
         class(Mesh1D_t)                        :: mesh
         class(Dimension_t), pointer            :: current
!
!        **********************************************
!           Allocate all dimensions
!        **********************************************
!
         current => self % dimension_head

         allocate( current ) 
!
!        Number of equations
!  
         current % value = NEC
         current % Name = "NEC"
         call check ( NF90_DEF_DIM ( self % fID , trim(current % Name) , current % value , current % dimID ) )

         allocate( current % next )
         current => current % next
!
!        Number of degrees of freedom
!
         current % value = 1  ! TODO
         current % Name = "nDOF"
         call check ( NF90_DEF_DIM ( self % fID , trim(current % Name) , current % Value , current % dimID ) )

         allocate( current % next ) 
         current => current % next
!
!        One
!  
         current % Value = 1
         current % Name  = "one"
         call check ( NF90_DEF_DIM ( self % fID , trim(current % Name) , current % Value , current % dimID ) )         
         

      end subroutine FileWriting_writeDimensions
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

end module FileWriting
