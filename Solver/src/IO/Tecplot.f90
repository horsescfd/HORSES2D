module Tecplot
   use SMConstants

   private
   public         :: ExportToTecplot

   integer, parameter         :: STR_LEN_TECPLOT = 128

   type Tecplot_t
      integer        :: no_of_variables = 0
      character(len=STR_LEN_TECPLOT), allocatable   :: variables(:)
      contains
         procedure      :: gatherVariables      => Tecplot_GatherVariables
   end type Tecplot_t

   type LinkedList_t
      integer        :: no_of_entries = 0
      class(Charlist), pointer    :: HEAD => NULL()
   end type LinkedList_t

   type Charlist
      character(len=STR_LEN_TECPLOT)      :: str
      class(Charlist), pointer            :: next => NULL()
   end type Charlist

   interface ExportToTecplot
      module procedure Tecplot_Save
   end interface ExportToTecplot

!
!  ========
   contains
!  ========
!
      subroutine TecPlot_Save( mesh ) 
         use QuadMeshClass
         implicit none
         class(QuadMesh_t)       :: mesh
         type(Tecplot_t)         :: tec 
         integer                 :: var


         call tec % GatherVariables

         do var = 1 , tec % no_of_variables
            print*, trim( tec % variables(var) )
         end do

      end subroutine

      subroutine TecPlot_GatherVariables( self ) 
         use Setup_Class
         implicit none
         class(Tecplot_t)               :: self
         logical                        :: flag = .true.
         integer                        :: pos
         character(len=STR_LEN_TECPLOT) :: auxstr
         type(LinkedList_t)             :: entries
         class(CharList), pointer       :: current
         integer                        :: i

         auxstr = Setup % saveVariables

!        Prepare the linked list
!        -----------------------

         do 

            pos = index(auxstr , "_")

            if ( pos .eq. 0 ) then        ! Is not present: All the string is a variable
!
!              Prepare a new entry in the list
!              -------------------------------
               if ( entries % no_of_entries .eq. 0) then
                  allocate( entries % HEAD ) 
                  current => entries % HEAD
               else
                  allocate(current % next)
                  current => current % next
               end if 

               entries % no_of_entries = entries % no_of_entries + 1 
            
               current % str = auxstr
               auxstr        = auxstr
               
               exit

            else
!
!              Prepare a new entry in the list
!              -------------------------------
               if ( entries % no_of_entries .eq. 0) then
                  allocate( entries % HEAD ) 
                  current => entries % HEAD
               else
                  allocate(current % next)
                  current => current % next
               end if 

               entries % no_of_entries = entries % no_of_entries + 1 
            
               current % str = auxstr(1:pos-1)
               auxstr        = auxstr(pos+1:)
            end if
               
         end do

!
!        Store the results in the tecplot typedef
!        ----------------------------------------
         allocate( self % variables ( entries % no_of_entries ) )
         current => entries % HEAD

         self % no_of_variables = entries % no_of_entries
         do i = 1 , entries % no_of_entries
            self % variables(i)  = current % str 
            current => current % next
         end do

      end subroutine TecPlot_gatherVariables

end module Tecplot
