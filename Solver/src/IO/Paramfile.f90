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
   public   readValue , readValueInRegion , getSquashedLine
   character, parameter       :: comment = '!'
   character, parameter       :: equal(2) = ['=',':'] 
   integer,   parameter       :: STR_LEN_PARAM = 512

   interface readValue
      module procedure readCharacterValue , readIntegerValue , readRealValue , readLogicalValue
   end interface readValue

   interface readValueInRegion
      module procedure readCharacterValueInRegion , readIntegerValueInRegion , readRealValueInRegion , readLogicalValueInRegion
   end interface readValueInRegion
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
               var = ""
               exit

            elseif (io .gt. 0) then
               print*, "IOSTAT returned a positive number"
               stop "Stopped."
         
            else
               
!              Removed commented part of the string
!              ------------------------------------
               position = index(trim(auxstr) , comment)
               if ( position .gt. 0 ) then
                  auxstr = auxstr(1:position-1)
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
               position = max(index(trim(auxstr) , equal(1)) , index(trim(auxstr) , equal(2) ) )
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
         integer, allocatable         :: var
!        -------------------------------------------------
         character(len=STR_LEN_PARAM)        :: auxstr
         integer                             :: io

         call readCharacterValue(fileName , label , auxstr)

         auxstr = trim(adjustl(auxstr))

         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if
      
      end subroutine readIntegerValue

      subroutine readRealValue(fileName , label , var)
         implicit none
         character(len=*), intent(in)        :: fileName
         character(len=*), intent(in) :: label
         real(kind=RP), allocatable   :: var
!        -------------------------------------------------
         character(len=STR_LEN_PARAM)        :: auxstr
         integer                             :: io

         call readCharacterValue(fileName , label , auxstr)

         auxstr = trim(adjustl(auxstr))
      
         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if

      end subroutine readRealValue

      subroutine readLogicalValue(fileName , label , var)
         implicit none
         character(len=*), intent(in)        :: fileName
         character(len=*), intent(in) :: label
         logical, allocatable         :: var
!        -------------------------------------------------
         character(len=STR_LEN_PARAM)        :: auxstr
         integer                             :: io

         call readCharacterValue(fileName , label , auxstr)

         auxstr = trim(adjustl(auxstr))

         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if
      
      end subroutine readLogicalValue

      subroutine readCharacterValueInRegion(fileName , label , var , in_label , out_label)
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         character(len=*)             :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
!        -------------------------------------------------------------
         integer               :: fID
         character(len=STR_LEN_PARAM)            :: auxstr
         integer               :: io
         integer               :: position
         logical               :: inside
         
         inside = .false.
!
!        Open file
!        ---------
         open ( newunit = fID , file = trim(fileName) , status = 'old' , action = 'read' ) 

         do 
            read( fID , '(A)' , iostat = io ) auxstr
      
            if (io .lt. 0) then
               var = ""
               exit

            elseif (io .gt. 0) then
               print*, "IOSTAT returned a positive number"
               stop "Stopped."
         
            else
               
!
!              Check if inside a zone
!              ----------------------
               if ( getSquashedLine(auxstr) .eq. getSquashedLine(in_label) ) then
                  inside = .true.
                  cycle
               elseif( getSquashedLine(auxstr) .eq. getSquashedLine(out_label) ) then
                  inside = .false.
                  cycle
               end if

!              Removed commented part of the string
!              ------------------------------------
               if (inside) then
                  position = index(trim(auxstr) , comment)
                  if ( position .gt. 0 ) then
                     auxstr = auxstr(1:position-1)
                  end if
!   
!                 Look for the label
!                 ------------------
                  position = index(trim(auxstr) , trim(label) )
                  if ( position .eq. 0) then
                     cycle
                  end if

   
!   
!                 The label is present. The value will be everything from the equal
!                 -----------------------------------------------------------------
                  auxstr = trim(auxstr(position + len_trim(label) :))
                  position = max(index(trim(auxstr) , equal(1)) , index(trim(auxstr) , equal(2) ) )
                  if ( position .eq. 0) then
                     cycle
                  else
                     auxstr = auxstr(position+1:)
                  end if
!   
!                 Get the value
!                 -------------
                  var = trim(adjustl(auxstr))
                  exit
               end if

            end if
         end do

!
!        Close file
!        ----------
         close ( fID ) 
      

      end subroutine readCharacterValueInRegion

      subroutine readIntegerValueInRegion(fileName , label , var , in_label , out_label )
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         integer, allocatable         :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
         character(len=STR_LEN_PARAM) :: auxstr
         integer                      :: io

         call readCharacterValueInRegion(fileName , label, auxstr , in_label , out_label)

         auxstr = trim(adjustl(auxstr))
         
         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if
 

      end subroutine readIntegerValueInRegion

      subroutine readRealValueInRegion(fileName , label , var , in_label , out_label)
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         real(kind=RP), allocatable   :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
         character(len=STR_LEN_PARAM) :: auxstr
         integer                      :: io
      
         call readCharacterValueInRegion(fileName , label, auxstr , in_label , out_label)

         auxstr = trim(adjustl(auxstr))
         
         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if
      
      end subroutine readRealValueInRegion

      subroutine readLogicalValueInRegion(fileName , label , var , in_label , out_label)
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=*), intent(in) :: label
         logical, allocatable         :: var
         character(len=*), intent(in) :: in_label
         character(len=*), intent(in) :: out_label
         character(len=STR_LEN_PARAM) :: auxstr
         integer                      :: io
      
         call readCharacterValueInRegion(fileName , label, auxstr , in_label , out_label)

         auxstr = trim(adjustl(auxstr))

         if( .not. allocated(var)) allocate(var)
         read(auxstr,*,iostat=io) var
         if (io .lt. 0) then
            deallocate(var)
         end if

      end subroutine readLogicalValueInRegion

      function getSquashedLine(line) result (squashed)
         implicit none
         character(len=*), intent(in)     :: line
         character(len=STR_LEN_PARAM)     :: squashed
         integer                          :: pos

         pos = index(trim(adjustl(line)),' ')

         if (pos .eq. 0) then
            squashed = trim(adjustl(line))
         else
            squashed = line(1:pos-1) // line(pos+1:)
         end if

         do
            pos = index(trim(adjustl(squashed)),' ')

            if (pos .eq. 0) then
               squashed = trim(adjustl(squashed))
               return
            else
               squashed = squashed(1:pos-1) // squashed(pos+1:)
            end if
         end do
            

      end function getSquashedLine


end module ParamfileIO
