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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     File: Stopwatch.f90
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
module StopwatchClass
   use SMConstants
   implicit none

#define IYEAR 1
#define IMONTH 2
#define IDAY 3
#define IHOUR 5
#define IMINUTES 6
#define ISECONDS 7
#define IMILLISECONDS 8

   private
   public   Stopwatch

   integer, parameter      :: STR_LEN_STOPWATCH = 128
   integer, parameter      :: reference_date(8) = [2017,6,1,0,0,0,0,0]

   type Event_t
      character(len=STR_LEN_STOPWATCH)    :: name
      logical                             :: running
      real(kind=RP)                       :: tstart
      real(kind=RP)                       :: elapsedTime
      class(Event_t), pointer, private    :: next  => NULL()
      contains
         procedure      :: Start          => Event_Start
         procedure      :: Pause          => Event_Pause
         procedure      :: Reset          => Event_Reset
         procedure      :: GetElapsedTime => Event_GetElapsedTime
         procedure      :: Destruct       => Event_Destruct
   end type Event_t

   type Stopwatch_t
      class(Event_t),   pointer, private  :: head => NULL()
      integer, private                    :: no_of_events = 0
      contains
         procedure   :: CreateNewEvent => Stopwatch_CreateNewEvent
         procedure   :: Start          => Stopwatch_Start
         procedure   :: Pause          => Stopwatch_Pause
         procedure   :: Reset          => Stopwatch_Reset
         procedure   :: ElapsedTime    => Stopwatch_ElapsedTime
         procedure   :: Destruct       => Stopwatch_Destruct
   end type Stopwatch_t

   type(Stopwatch_t)    :: Stopwatch 

   contains
!
!/////////////////////////////////////////////////////////////////////
!
!           STOPWATCH PROCEDURES
!           --------------------
!/////////////////////////////////////////////////////////////////////
!
   function NewStopwatch()
      implicit none
      type(Stopwatch_t)       :: NewStopwatch

      NewStopwatch % head => NULL()
      NewStopwatch % no_of_events = 0

   end function NewStopwatch

   subroutine Stopwatch_CreateNewEvent( self , Name )
      implicit none
      class(Stopwatch_t)      :: self
      character(len=*)    :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                 :: i
      class(Event_t), pointer :: event

      if ( self % no_of_events .eq. 0 ) then
         self % no_of_events = 1 
         allocate(self % head)

         select type (ev => self % head)
            type is (Event_t)
               ev = NewEvent(trim(Name))
         end select

      else

         event => self % head
         do i = 1 , self % no_of_events-1
            if ( trim(event % name) .eq. trim(name) ) then
               print*, "Warning: The event already exists."
               return
            end if
    
            event => event % next

         end do

         if ( trim(event % name) .eq. trim(name) ) then
            print*, "Warning: The event already exists."
         end if
        
         self % no_of_events = self % no_of_events + 1 
         allocate(event % next)
         event => event % next
   
         select type (ev => event)
            type is (Event_t)
               ev = NewEvent(trim(Name))
         end select

      end if
      
   end subroutine Stopwatch_CreateNewEvent

   subroutine Stopwatch_Start(self , Name)
      implicit none
      class(Stopwatch_t)               :: self
      character(len=*) :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 

      event => self % head
      do i = 1 , self % no_of_events

         if ( trim(event % Name) .eq. trim(Name) ) then
            call event % Start()
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

   end subroutine Stopwatch_Start

   subroutine Stopwatch_Pause(self , Name)
      implicit none
      class(Stopwatch_t)               :: self
      character(len=*) :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd

      event => self % head
      do i = 1 , self % no_of_events

         if ( trim(event % Name) .eq. trim(Name) ) then
            call event % Pause
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

   end subroutine Stopwatch_Pause

   function Stopwatch_ElapsedTime(self,Name,units) result(elapsedTime)
      implicit none
      class(Stopwatch_t)                     :: self
      character(len=*)                       :: Name
      character(len=*), intent(in), optional :: units
      real(kind=RP)                          :: elapsedTime
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd
      character(len=STR_LEN_STOPWATCH)       :: chosenUnits

      if ( present(units) ) then
         chosenUnits = units

      else
         chosenUnits = "sec"

      end if

      event => self % head
      do i = 1 , self % no_of_events 

         if ( trim(event % Name) .eq. trim(Name) ) then
            elapsedTime = event % GetElapsedTime(chosenUnits)
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

      elapsedTime = -1.0_RP

   end function Stopwatch_ElapsedTime

   subroutine Stopwatch_Reset(self,Name) 
      implicit none
      class(Stopwatch_t)               :: self
      character(len=*) :: Name
!
!     ---------------
!     Local variables
!     ---------------
!
      class(Event_t), pointer    :: event => NULL()
      integer                    :: i 
      real(kind=RP)              :: tEnd

      event => self % head
      do i = 1 , self % no_of_events

         if ( trim(event % Name) .eq. trim(Name) ) then
            call event % Reset
            return
         else
            event => event % next
         end if
      end do

      print*, "Warning: Stopwatch event ",trim(Name)," was not found."

   end subroutine Stopwatch_Reset

   subroutine Stopwatch_Destruct(self)
      implicit none
      class(Stopwatch_t)      :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i
      class(Event_t), pointer    :: event , nextEvent

      event => self % head
      do i = 1 , self % no_of_events 
         nextEvent => event % next
         call event % Destruct  
         deallocate(event)
         event => nextEvent

      end do

      self % no_of_events = 0
      self % head => NULL()

   end subroutine Stopwatch_Destruct
!
!/////////////////////////////////////////////////////////////////////
!
!           EVENT PROCEDURES
!           ----------------
!/////////////////////////////////////////////////////////////////////
!
   function NewEvent(Name)
      implicit none
      character(len=*), intent(in) :: Name
      type(Event_t)                :: NewEvent
      
      NewEvent % Name        = trim(Name)
      NewEvent % running     = .false.
      NewEvent % tstart      = 0.0_RP
      NewEvent % elapsedTime = 0.0_RP
      NewEvent % next        => NULL()

   end function NewEvent

   subroutine Event_Start(self)
      implicit none
      class(Event_t)    :: self

      self % running = .true.
      self % tstart = GetCurrentSeconds()

   end subroutine Event_Start

   subroutine Event_Pause(self)
      implicit none
      class(Event_t)    :: self
      real(kind=RP)     :: tEnd 

      tEnd = GetCurrentSeconds()
      self % elapsedTime = tEnd - self % tStart + self % elapsedTime
      self % tStart = 0.0_RP
      self % running = .false.

   end subroutine Event_Pause

   function Event_GetElapsedTime(self,units) result (time)
      implicit none
      class(Event_t)               :: self
      character(len=*), intent(in) :: units
      real(kind=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP)     :: tActual

      if ( self % running ) then
         tActual = GetCurrentSeconds()
         time = tActual - self % tStart + self % elapsedTime

      else
         time = self % elapsedTime

      end if

      if ( (trim(units) .eq. "minutes") .or. (trim(units) .eq. "Minutes") .or. (trim(units) .eq. "min") ) then
         time = time / 60.0_RP
   
      else if ( (trim(units) .eq. "hours") .or. (trim(units) .eq. "Hours") .or. (trim(units) .eq. "h") ) then
         time = time / 3600.0_RP

      else if ( (trim(units) .eq. "days") .or. (trim(units) .eq. "Days") .or. (trim(units) .eq. "d") ) then
         time = time / (24.0_RP * 3600.0_RP)
         
      end if

   end function Event_GetElapsedTime

   subroutine Event_Reset(self)
      implicit none
      class(Event_t)          :: self
      
      self % tStart = 0.0_RP
      self % elapsedTime = 0.0_RP
      self % running = .false.

   end subroutine Event_Reset

   subroutine Event_Destruct(self)
      implicit none
      class(Event_t)       :: self

      self % tStart = 0.0_RP
      self % elapsedTime = 0.0_RP
      self % running = .false.
      self % next => NULL()

   end subroutine Event_Destruct
!
!/////////////////////////////////////////////////////////////////////
!
!           AUXILIAR PROCEDURES
!           -------------------
!////////////////////////////////////////////////////////////////////
!
   function GetCurrentSeconds() result(sec)
      implicit none
      real(kind=RP)         :: sec
!
!     ---------------
!     Local variables
!     ---------------
!
      integer, dimension(8) :: currentDate
      integer               :: no_of_days

      call date_and_time (values = currentDate)

      
      no_of_days = DaysBetweenTwoDates( reference_date(IDAY) , reference_date(IMONTH) , reference_date(IYEAR) , &
                                        currentDate(IDAY)    , currentDate(IMONTH)    , currentDate(IYEAR) )

      sec = no_of_days * 24.0_RP * 3600.0_RP + SecondsToMidnight(currentDate(IHOUR),&
                                                                 currentDate(IMINUTES),&
                                                                 currentDate(ISECONDS),&
                                                                 currentDate(IMILLISECONDS) )
   end function GetCurrentSeconds
   
   function SecondsToMidnight(hour,minute,seconds,millisecs) result (secs)
      implicit none
      integer, intent(in)     :: hour
      integer, intent(in)     :: minute
      integer, intent(in)     :: seconds
      integer, intent(in)     :: millisecs
      real(kind=RP)           :: secs

      secs = seconds + millisecs/1000.0_RP + minute * 60.0_RP + hour * 3600.0_RP

   end function SecondsToMidnight

   function DaysBetweenTwoDates(d1,m1,y1,d2,m2,y2) result (days)
      implicit none
      integer, intent(in)        :: d1
      integer, intent(in)        :: m1
      integer, intent(in)        :: y1
      integer, intent(in)        :: d2
      integer, intent(in)        :: m2
      integer, intent(in)        :: y2
      integer                    :: days

      days = DaysBetweenTwoYears(y1,y2) + DaysToNewYearsDay(d2,m2,y2) - DaysToNewYearsDay(d1,m1,y1)

   end function DaysBetweenTwoDates

   function DaysBetweenTwoYears(year1,year2) result(days)
      implicit none
      integer, intent(in)        :: year1 
      integer, intent(in)        :: year2
      integer                    :: days
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: no_of_leapYears
      integer        :: i 

      no_of_leapYears = 0

      do i = year1 , year2 - 1
         if ( checkIfLeapYear(i) ) then
            no_of_leapYears = no_of_leapYears + 1 
         end if
      end do

      days = 365 * (year2-year1) + no_of_leapYears

   end function DaysBetweenTwoYears

   function DaysToNewYearsDay(day,month,year) result(days)
      implicit none
      integer, intent(in)        :: day
      integer, intent(in)        :: month
      integer, intent(in)        :: year
      integer                    :: days
      integer, parameter         :: days_per_month(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: i

      days = day

      do i = 1 , month-1
         days = days + days_per_month(i)
      end do

      if ( checkIfLeapYear(year) .and. (month .gt. 2)) then
         days = days + 1 
      end if 

   end function DaysToNewYearsDay

   function checkIfLeapYear(year) result (leap)
      implicit none
      integer, intent(in)     :: year
      logical                 :: leap

      if ( mod(year,4) .ne. 0 ) then
         leap = .false.
         return
      end if

      if ( (mod(year,100) .eq. 0) .and. (mod(year,400) .ne. 0 )) then
         leap = .false.
         return

      else
         leap = .true.
         return
      end if

   end function checkIfLeapYear
      

end module StopwatchClass
