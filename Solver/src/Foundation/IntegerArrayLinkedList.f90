!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     File: IntegerArrayLinkedList.f90
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
module IntegerArrayLinkedList
   use Sorting
!
   private 
   public ConstructIntegerArrayLinkedList , IntegerArrayLinkedList_t , IntegerArrayEntry_t
   public WITHOUT_REPEATING , WITH_REPEATING

   integer, parameter      :: WITHOUT_REPEATING = 1
   integer, parameter      :: WITH_REPEATING    = 2
!
!//////////////////////////////////////////////////////////////////////////////////
!
!  Linked list class
!  -----------------
   type IntegerArrayLinkedList_t
      class(IntegerArrayEntry_t), pointer    :: head => NULL()
      class(IntegerArrayEntry_t), pointer    :: tail => NULL() 
      integer, pointer                       :: no_of_entries
      integer, private                       :: type = WITHOUT_REPEATING
      contains
         procedure   :: Get               => IntegerArrayLinkedList_Get
         procedure   :: Add               => IntegerArrayLinkedList_Add
         procedure   :: Search            => IntegerArrayLinkedList_Search
         procedure   :: SearchIfContained => IntegerArrayLinkedList_SearchIfContained
         procedure   :: SearchIfPresent   => IntegerArrayLinkedList_SearchIfPresent
         procedure   :: List              => IntegerArrayLinkedList_List
         procedure   :: Remove            => IntegerArrayLinkedList_Remove
         procedure   :: CheckConsistency  => IntegerArrayLinkedList_CheckConsistency
         procedure   :: Destruct          => IntegerArrayLinkedList_Destruct
   end type IntegerArrayLinkedList_t
!
!  Linked list entry class
!  -----------------------
   type IntegerArrayEntry_t
      integer                             :: N
      logical                             :: attribute = .false.
      integer, allocatable                :: val(:)
      class(IntegerArrayEntry_t), pointer :: next => NULL()
      class(IntegerArrayEntry_t), pointer :: prev => NULL()
   end type IntegerArrayEntry_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           BUILDING SUBROUTINES
!           --------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function ConstructIntegerArrayLinkedList(type) result (LL)
         implicit none
         integer, intent(in), optional       :: type
         type(IntegerArrayLinkedList_t)      :: LL

         LL % tail => NULL()
         LL % head => NULL()

         allocate ( LL % no_of_entries ) 

         if ( present (type) ) then
            LL % type = type
         else
            LL % type = WITHOUT_REPEATING
         end if
         
         LL % no_of_entries = 0

      end function ConstructIntegerArrayLinkedList

      subroutine IntegerArrayLinkedList_Add( self , N , array ) 
         implicit none
         class(IntegerArrayLinkedList_t)        :: self
         integer, intent(in)                    :: N
         integer, intent(in)                    :: array(N)
         class(IntegerArrayEntry_t), pointer    :: current

         if ( (self % type .eq. WITHOUT_REPEATING) .and. (self % Search(N,array) .ne. -1) ) then
!
!           It is already stored
!           --------------------
            return

         end if

         if ( self % no_of_entries .eq. 0 ) then
            if ( associated(self % head) ) stop "Fatal error"
            allocate ( self % head )
            
            current => self % head
            self % tail    => self % head

         else
!
!           Point to the end of the list
!           ----------------------------
            current => self % tail
!
!           Create a new slave and link it in the list
!           ------------------------------------------
            allocate ( self % tail % next ) 
            self % tail        => self % tail % next
            self % tail % prev => current
            current % next     => self % tail
            current            => self % tail
            

         end if
!
!        Save the data
!        -------------
         allocate ( current % val ( N ) ) 
         current % N             = N
         current % val           = array
         current % attribute     = .false.
         self    % no_of_entries = self % no_of_entries + 1

      end subroutine IntegerArrayLinkedList_Add

      subroutine IntegerArrayLinkedList_Remove( self , which ) 
         implicit none
         class(IntegerArrayLinkedList_t)     :: self
         integer, intent(in)                 :: which
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                             :: i
         class(IntegerArrayEntry_t), pointer :: current

         if ( which .gt. self % no_of_entries ) then
            return

         elseif ( which .le. 0 ) then
            return

         elseif ( self % no_of_entries .eq. 1 ) then
!
!           Reset the linked list
!           ---------------------
            self % no_of_entries = 0
            deallocate( self % head )
            self % head => NULL()
            self % tail => NULL()

         elseif ( which .eq. 1 ) then
!
!           Remove the head
!           ---------------
            current => self % head

            self % head        => current % next
            self % head % prev => NULL()
            deallocate(current)
            self % no_of_entries = self % no_of_entries - 1 
            return

         elseif ( which .eq. self % no_of_entries ) then
!
!           Remove the tail
!           ---------------
            current => self % tail

            self % tail        => current % prev
            self % tail % next => NULL()
            deallocate(current)
            self % no_of_entries = self % no_of_entries - 1 
            return

         else
         
            current => self % head    

            do i = 1 , which - 1
               current => current % next
            end do
!
!           Current item is the one to remove
!           ---------------------------------
            current % next % prev => current % prev 
            current % prev % next => current % next
            deallocate(current)
            self % no_of_entries = self % no_of_entries - 1 
            return

         end if
            
   
      end subroutine IntegerArrayLinkedList_Remove

      subroutine IntegerArrayLinkedList_Destruct( self ) 
         implicit none
         class(IntegerArrayLinkedList_t)     :: self
!
!        ---------------
!        Local variables         
!        ---------------
!
         integer     :: i

         do i = 1 , self % no_of_entries
            call self % Remove(1)
         
         end do

         self % head => NULL()
         self % tail => NULL()
         deallocate ( self % no_of_entries ) 

      end subroutine IntegerArrayLinkedList_Destruct


      function IntegerArrayLinkedList_CheckConsistency(self) result ( code )
         implicit none
         class(IntegerArrayLinkedList_t)        :: self
         integer                                :: code
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i
         class(IntegerArrayEntry_t), pointer    :: current , previous

         if ( associated ( self % head % prev ) ) then
            code = -1
            return

         elseif ( associated ( self % tail % next ) ) then
            code = -1   
            return

         elseif ( self % no_of_entries .eq. 1 ) then
            if ( .not. associated ( self % head , self % tail ) ) then
               code = -1 
               return
            elseif ( associated( self % head % next ) ) then
               code = -1
               return
            elseif ( associated( self % tail % prev ) ) then
               code = -1 
               return
            end if

         else
            current  => self % head

            do i = 1 , self % no_of_entries - 1
               previous => current
               current  => current % next

               if ( .not. associated( current % prev , previous ) ) then
                  code = -1 
                  return
               elseif ( .not. associated ( previous % next , current ) ) then
                  code = -1 
                  return

               end if

            end do

            if ( .not. associated ( current , self % tail ) ) then
               code = -1 
               return

            elseif ( .not. associated ( previous , self % tail % prev ) ) then
               code = -1 
               return

            end if
            
         end if


         code = 1 

      end function IntegerArrayLinkedList_CheckConsistency
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           SEARCHING SUBROUTINES
!           ---------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function IntegerArrayLinkedList_Get( self , N , which ) result ( val ) 
         implicit none
         class(IntegerArrayLinkedList_t) ,  intent (in) :: self
         integer                         ,  intent (in) :: N
         integer                         ,  intent (in) :: which
         integer                                        :: val (N)
!
!        ---------------
!        Local variables
!        ---------------
!
         class(IntegerArrayEntry_t), pointer       :: current
         integer                                   :: i 
!
!        Get the LinkedList head
!        -----------------------
         current => self % head
!
!        Loop the LinkedList until the desired position is found
!        -------------------------------------------------------
         do i = 1 , which - 1 
            current => current % next
         end do
!
!        Return its value
!        ----------------
         val = current % val(1:N)

      end function IntegerArrayLinkedList_Get

      function IntegerArrayLinkedList_Search( self , N , array ) result ( pos )
!
!        **********************************************************************
!
!              This subroutine checks if an entry has the value in "array".
!           There is just one possible solution, since the LinkedList does
!           not include repeated solution.
!           Recall that the LinkedList does not account for the order of the
!           stored values in the array. E.g.:
!
!             >   array = [1,2] can be found in entry % val = [2,1]
!
!        **********************************************************************
!
         implicit none
         class(IntegerArrayLinkedList_t)     :: self
         integer, intent(in)                 :: N
         integer, intent(in)                 :: array(N)
         integer                             :: pos
!
!        ---------------
!        Local variables
!        ---------------
!
         class(IntegerArrayEntry_t), pointer :: current
         integer                             :: array_in(N)
         integer                             :: array_stored(N)
         integer                             :: i , j
!
!        Sort the array
!        --------------
         array_in = array
         call Qsort(array_in)
!
!        Loop the list to find the item
!        ------------------------------
         current => self % head

list:    do i = 1 , self % no_of_entries 
            if ( current % N .eq. N ) then
!
!              The item has the same length than the requested
!              -----------------------------------------------
               array_stored = current % val
               call Qsort(array_stored)

               do j = 1 , N
                  if ( array_in(j) .ne. array_stored(j) ) then
                     current => current % next
                     cycle list
                  end if
               end do
!
!              If arrives here, it has been found
!              ----------------------------------
               pos = i
               return

            end if

         end do   list
!
!        Default value: -1
!        -------------
         pos = -1

      end function IntegerArrayLinkedList_Search

      function IntegerArrayLinkedList_SearchIfContained( self , N , array , howMany ) result ( pos )
!
!        *******************************************************************************************
!
!              This subroutine Searches if a subarray is entirely contained in val. It looks for 
!           "howMay" solutions (if no_of_solutions < howMany, they are filled with -1). 
!           It does not account for the order. E.g.:
!
!                 > [1,2] is found in [4,2,5,1]    (since both 1 and 2 are contained)
!                 > [1,2] is not found in [1,3]    (since 2 is missing)
!
!           It returns the positions in the LinkedList of these entries.
!
!        *******************************************************************************************
!
         implicit none
         class(IntegerArrayLinkedList_t)     :: self
         integer, intent(in)                 :: N
         integer, intent(in)                 :: array(N)
         integer, intent(in)                 :: howMany
         integer                             :: pos(howMany)
!
!        ---------------
!        Local variables
!        ---------------
!
         class(IntegerArrayEntry_t), pointer :: current
         integer                             :: array_in(N)
         integer                             :: i , j , counter
!
!        Sort the array
!        --------------
         array_in = array
         call Qsort(array_in)
!
!        Loop the list to find the item
!        ------------------------------
         current => self % head
   
         counter = 0
list:    do i = 1 , self % no_of_entries 
            if ( current % N .ge. N ) then
               do j = 1 , N
                  if ( .not. any( current % val .eq. array_in(j) ) ) then
                     current => current % next
                     cycle list
                  end if
               end do
!
!              If arrives here, it has been found
!              ----------------------------------
               counter = counter + 1 
               pos(counter) = i
               if ( counter .eq. howMany ) then
                 return
               end if
!
!              Advance in the list
!              -------------------
               current => current % next
               cycle list

            else
               current => current % next
               cycle list
            end if

         end do   list
!
!        Default value: -1
!        -------------
         pos(counter+1:) = -1

      end function IntegerArrayLinkedList_SearchIfContained

      function IntegerArrayLinkedList_SearchIfPresent( self , N , array , howMany , attributes ) result ( pos )
!
!        *******************************************************************************************
!
!              This subroutine Searches whether an entry contains a certain entry in its value array 
!           "howMay" solutions are seeked (if no_of_solutions < howMany, they are filled with -1). 
!           E.g.:
!                 > [1,2] is found in [4,5,1]      (since 1 is present)
!                 > [1,2] is found in [2,4,3]      (since 2 is present)
!                 > [1,2] is found in [1,2,3]      (since both 1 and 2 are present)
!                 > [1,2] is not found in [3,4]    (since none of 1 and 2 are present)
!
!           It returns the positions in the LinkedList of these entries.
!
!        *******************************************************************************************
!

         implicit none
         class(IntegerArrayLinkedList_t)     :: self
         integer, intent(in)                 :: N
         integer, intent(in)                 :: array(N)
         integer, intent(in)                 :: howMany
         logical, intent(out), optional      :: attributes(howMany)
         integer                             :: pos(howMany)
!
!        ---------------
!        Local variables
!        --------------
!
         class(IntegerArrayEntry_t), pointer :: current
         integer                             :: i , j , counter
!
!        Loop the list to find the item
!        ------------------------------
         current => self % head
   
         counter = 0
list:    do i = 1 , self % no_of_entries 

            do j = 1 , N
               if (  any( current % val .eq. array(j) ) ) then
!
!                 Found
!                 -----
                  counter = counter + 1 
                  pos(counter) = i

                  if ( present ( attributes ) ) then
                     attributes(counter) = current % attribute
                  end if

                  if ( counter .eq. howMany ) then
                     return
                  end if
!
!                 Advance in the list
!                 -------------------
                  current => current % next
                  cycle list
               end if
            end do

!
!           Advance in the list
!           -------------------
            current => current % next

         end do   list
!
!        Default value: -1
!        -------------
         pos(counter+1:) = -1

         if ( present ( attributes ) ) then
            attributes(counter+1:) = .false.

         end if

      end function IntegerArrayLinkedList_SearchIfPresent
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!              DESCRIBE SUBROUTINES
!              --------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IntegerArrayLinkedList_List( self , UNIT ) 
         implicit none
         class(IntegerArrayLinkedList_t)     :: self
         integer, intent(in)                 :: UNIT
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i , j
         class(IntegerArrayEntry_t), pointer    :: current

         current => self % head

         do i = 1 , self % no_of_entries

            write(UNIT,'(A,I0,A)', advance = "no" ) "Position " , i , ": ["

            do j = 1 , current % N - 1 
               write(UNIT,'(I0,A)', advance = "no" ) current % val(j) ," , "
            end do

            write(UNIT,'(I0,A)'  ) current % val(current % N) , "] "

            current => current % next

         end do
         

      end subroutine IntegerArrayLinkedList_List

!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
end module IntegerArrayLinkedList
