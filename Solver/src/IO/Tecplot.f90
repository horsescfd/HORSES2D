module Tecplot
   use SMConstants

   private
   public         :: ExportToTecplot , ExportMeshToTecplot

   integer, parameter         :: STR_LEN_TECPLOT = 128

   type Tecplot_t
      integer        :: no_of_variables = 0
      integer        :: fID
      character(len=STR_LEN_TECPLOT), allocatable   :: variables(:)
      character(len=STR_LEN_TECPLOT)                :: Name
      
      contains
         procedure      :: gatherVariables      => Tecplot_GatherVariables
         procedure      :: Open                 => Tecplot_OpenFile
         procedure      :: NewZone              => Tecplot_NewZone
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

   interface ExportMeshToTecplot
      module procedure Tecplot_SaveMesh
   end interface ExportMeshToTecplot

!
!  ========
   contains
!  ========
!
      subroutine Tecplot_SaveMesh( mesh , Name ) 
         use QuadMeshClass
         implicit none
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
         type(Tecplot_t)         :: tec
         integer                 :: eID
         character(len=STR_LEN_TECPLOT)      :: auxname

         auxname = Name(1: len_trim(Name) - len(".HiOMesh")) // ".plt"

         tec % Name = trim(auxname)


         call tec % Open
      
         do eID = 1 , mesh % no_of_elements
            call tec % NewZone( mesh , eID ) 
         end do
   
         close ( tec % fID )

      end subroutine Tecplot_SaveMesh

      subroutine TecPlot_Save( mesh , Name) 
         use QuadMeshClass
         implicit none
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
         type(Tecplot_t)         :: tec 
         integer                 :: var
         integer                 :: eID

         tec % Name = trim(Name)

         call tec % GatherVariables

         call tec % Open

         do eID = 1 , mesh % no_of_elements
            call tec % NewZone( mesh , eID )
         end do

         close( tec % fID )

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

      subroutine Tecplot_OpenFile( self ) 
         implicit none
         class(Tecplot_t)        :: self
         integer                 :: var

         open( newunit = self % fID , file = trim(self % Name) , status = "unknown" , action = "write" ) 

!
!        Print the header into the file
!        ------------------------------
         write( self % fID , '(A,A,A)') 'TITLE = "',trim(self % Name),'"'
         write( self % fID , '(A)' , advance="no") 'VARIABLES = "X" "Y" "Z" '

         do var = 1 , self % no_of_variables
            write( self % fID , '(A,A,A)' , advance="no" ) '"',trim(self % variables(var)),'" '
         end do
         write( self % fID , * )

      end subroutine Tecplot_OpenFile

      subroutine Tecplot_NewZone( self , mesh , eID) 
         use QuadMeshClass
         use Physics
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: eID 
         real(kind=RP), pointer  :: rho (:,:) , rhou (:,:) , rhov (:,:) , rhoe (:,:)
         real(kind=RP), pointer  :: rhot(:,:) , rhout(:,:) , rhovt(:,:) , rhoet(:,:)
         integer                 :: iXi , iEta
         integer                 :: var


         associate ( N => mesh % elements(eID) % spA % N )
         
!        New header
!        ----------
         write( self % fID , '(A,I0,A)' , advance="no" ) "ZONE N=",(N+1)*(N+1),", "
         write( self % fID , '(A,I0,A)' , advance="no" ) "E=",(N)*(N),", "
         write( self % fID , '(A)'                     ) "DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"

!
!        Point to the quantities
!        -----------------------
         rho(0:,0:)  => mesh % elements(eID) % Q(0:,0:,IRHO) 
         rhou(0:,0:) => mesh % elements(eID) % Q(0:,0:,IRHOU)
         rhov(0:,0:) => mesh % elements(eID) % Q(0:,0:,IRHOV)
         rhoe(0:,0:) => mesh % elements(eID) % Q(0:,0:,IRHOE)
         rhot(0:,0:)  => mesh % elements(eID) % QDot(0:,0:,IRHO) 
         rhout(0:,0:) => mesh % elements(eID) % QDot(0:,0:,IRHOU)
         rhovt(0:,0:) => mesh % elements(eID) % QDot(0:,0:,IRHOV)
         rhoet(0:,0:) => mesh % elements(eID) % QDot(0:,0:,IRHOE)

         do iEta = 0 , N
            do iXi = 0 , N
               write( self % fID , '(E16.10,1X,E16.10,1X,E16.10)',advance="no") mesh % elements(eID) % x(iX,iXi,iEta) , mesh % elements(eID) % x(iY,iXi,iEta) , 0.0_RP  
!
!              Save quantities
!              ---------------
               do var = 1 , self % no_of_variables

                  select case ( trim( self % variables(var) ) )
                     case ("rho")
                        write(self % fID,'(1X,E16.10)',advance="no") rho(iXi,iEta) * refValues % rho

                     case ("rhou")
                        write(self % fID,'(1X,E16.10)',advance="no") rhou(iXi,iEta) * refValues % rho * refValues % V

                     case ("rhov")
                        write(self % fID,'(1X,E16.10)',advance="no") rhov(iXi,iEta) * refValues % rho * refValues % V

                     case ("rhoe")
                        write(self % fID,'(1X,E16.10)',advance="no") rhoe(iXi,iEta) * refValues % rho * refValues % p

                     case ("rhot")
                        write(self % fID,'(1X,E16.10)',advance="no") rhot(iXi,iEta) * refValues % rho / refValues % tc

                     case ("rhout")
                        write(self % fID,'(1X,E16.10)',advance="no") rhout(iXi,iEta) * refValues % rho * refValues % V / refValues % tc

                     case ("rhovt")
                        write(self % fID,'(1X,E16.10)',advance="no") rhovt(iXi,iEta) * refValues % rho * refValues % V / refValues % tc

                     case ("rhoet")
                        write(self % fID,'(1X,E16.10)',advance="no") rhoet(iXi,iEta) * refValues % rho * refValues % p / refValues % tc

                     case ("u")
                        write(self % fID,'(1X,E16.10)',advance="no") rhou(iXi,iEta)/rho(iXi,iEta) * refValues % V

                     case ("v")
                        write(self % fID,'(1X,E16.10)',advance="no") rhov(iXi,iEta)/rho(iXi,iEta) * refValues % V
   
                     case ("p")
                        write(self % fID,'(1X,E16.10)',advance="no") Thermodynamics % gm1 * ( rhoe(iXi,iEta) - 0.5*rhou(iXi,iEta)*rhou(iXi,iEta)/rho(iXi,iEta) - 0.5*rhov(iXi,iEta)*rhov(iXi,iEta)/rho(iXi,iEta) ) * refValues % p
      
                     case ("Mach")
                        write(self % fID,'(1X,E16.10)',advance="no") sqrt(rhou(iXi,iEta)*rhou(iXi,iEta)+rhov(iXi,iEta)*rhov(iXi,iEta))/rho(iXi,iEta)

                  end select                        

               end do

!              Jump to next line
!              -----------------
               write( self % fID , *)

            end do
         end do

         write( self % fID , * )    ! One blank line

         do iEta = 1 , N
            do iXi = 1 , N
               write(self % fID , '(I0,1X,I0,1X,I0,1X,I0)')  pointPosition(iXi,iEta,N)
            end do
         end do
         end associate


      end subroutine  Tecplot_NewZone

      function pointPosition(iXi , iEta , N) result( val )
         use QuadMeshDefinitions
         implicit none
         integer        :: iXi
         integer        :: iEta
         integer        :: N
         integer        :: val(POINTS_PER_QUAD)

         val(1) = (N+1)*(iEta-1) + iXi + 1
         val(2) = (N+1)*(iEta-1) + iXi 
         val(3) = (N+1)*iEta + iXi
         val(4) = (N+1)*iEta + iXi + 1
      end function pointPosition

end module Tecplot  
