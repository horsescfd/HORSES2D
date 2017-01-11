!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
module ParamfileIO
   use SMConstants
   implicit none
!
!  ******************
!  File configuration
!  ******************
!
   private
   public   readValue
   character, parameter       :: comment = '!'
   character, parameter       :: equal = '='
   integer,   parameter       :: STR_LEN_PARAM = 512

   interface readValue
      module procedure readCharacterValue , readIntegerValue , readRealValue , readLogicalValue
   end interface readValue
!
!   interface readValueInRegion
!      module procedure readCharacterValueInRegion , readIntegerValueInRegion , readLogicalValueInRegion
!   end interface readValueInRegion
!
   
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
      subroutine readCharacterValue(fileName , label , var)
         implicit none
         character(len=*), intent(in)        :: fileName
         character(len=*), intent(in)        :: label
         character(len=*)       :: var
!        -------------------------------------------------------------
         integer               :: fID
         character(len=STR_LEN_PARAM)            :: auxstr
         integer               :: io
         integer               :: position
         
!
!        Open file
!        ---------
         open ( newunit = fID , file = trim(fileName) , status = 'old' , action = 'read' ) 

         do 
            read( fID , '(A)' , iostat = io ) auxstr
      
            if (io .lt. 0) then
               print*, 'The value "',trim(label),'" is not present in the file ',trim(fileName),'.'
               exit

            elseif (io .gt. 0) then
               print*, "IOSTAT returned a positive number"
               stop "Stopped."
         
            else
               
!              Removed commented part of the string
!              ------------------------------------
               position = index(trim(auxstr) , comment)
               if ( position .gt. 0 ) then
                  auxstr = auxstr(1:position)
               end if
!
!              Look for the label
!              ------------------
               position = index(trim(auxstr) , trim(label) )
               if ( position .eq. 0) then
                  cycle
               end if

!
!              The label is present. The value will be everything from the equal
!              -----------------------------------------------------------------
               position = index(trim(auxstr) , equal)
               if ( position .eq. 0) then
                  cycle
               else
                  auxstr = auxstr(position+1:)
               end if
!
!              Get the value
!              -------------
               var = trim(adjustl(auxstr))
               exit

            end if
         end do

!
!        Close file
!        ----------
         close ( fID ) 
      
      end subroutine readCharacterValue

      subroutine readIntegerValue(fileName , label , var)
         implicit none
         character(len=*), intent(in)        :: fileName
         character(len=*), intent(in) :: label
         integer         :: var
!        -------------------------------------------------
         character(len=STR_LEN_PARAM)        :: auxstr

         call readCharacterValue(fileName , label , auxstr)

         auxstr = trim(adjustl(auxstr))
         read(auxstr,*) var
      
      end subroutine readIntegerValue

      subroutine readRealValue(fileName , label , var)
         implicit none
         character(len=*), intent(in)        :: fileName
         character(len=*), intent(in) :: label
         real(kind=RP)   :: var
!        -------------------------------------------------
         character(len=STR_LEN_PARAM)        :: auxstr

         call readCharacterValue(fileName , label , auxstr)

         auxstr = trim(adjustl(auxstr))
         read(auxstr,*) var
      
      end subroutine readRealValue

      subroutine readLogicalValue(fileName , label , var)
         implicit none
         character(len=*), intent(in)        :: fileName
         character(len=*), intent(in) :: label
         logical         :: var
!        -------------------------------------------------
         character(len=STR_LEN_PARAM)        :: auxstr

         call readCharacterValue(fileName , label , auxstr)

         auxstr = trim(adjustl(auxstr))
         read(auxstr,*) var
      
      end subroutine readLogicalValue

      subroutine readCharacterValueInRegion(file , label , var)
         implicit none
         character(len=*), intent(in)        :: file
         character(len=*), intent(in)        :: label
         character(len=STR_LEN_PARAM)       :: var
      
         var = ''

      end subroutine readCharacterValueInRegion

      subroutine readIntegerValueInRegion(file , label , var)
         implicit none
         character(len=*), intent(in)        :: file
         character(len=*), intent(in) :: label
         integer         :: var
      
         var = 0
      end subroutine readIntegerValueInRegion

      subroutine readRealValueInRegion(file , label , var)
         implicit none
         character(len=*), intent(in)        :: file
         character(len=*), intent(in) :: label
         real(kind=RP)   :: var
      
         var = 0.0_RP
      
      end subroutine readRealValueInRegion

      subroutine readLogicalValueInRegion(file , label , var)
         implicit none
         character(len=*), intent(in)        :: file
         character(len=*), intent(in) :: label
         logical         :: var
      
         var = .false.
      end subroutine readLogicalValueInRegion


end module ParamfileIO
