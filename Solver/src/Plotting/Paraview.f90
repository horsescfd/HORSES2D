submodule (Plotter) Paraview
   use SMConstants

#include "Defines.h"

   type LinkedList_t
      integer        :: no_of_entries = 0
      class(Charlist), pointer    :: HEAD => NULL()
      contains
         procedure   :: Destruct => LinkedList_Destruct
   end type LinkedList_t

   type Charlist
      character(len=STR_LEN_PLOTTER)      :: str
      class(Charlist), pointer            :: next => NULL()
   end type Charlist

   integer     :: point_position = 0 
!
!
!  ========
   contains
!  ========
!
      module subroutine Paraview_Initialization ( self )
         use QuadMeshClass
         class(Paraview_t)     :: self

         call Paraview_GatherVariables( self )

      end subroutine Paraview_Initialization
         
      module subroutine Paraview_ExportMesh( self , mesh , Name ) 
         use QuadMeshClass
         implicit none
         class(Paraview_t)        :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
         integer                 :: eID
         character(len=STR_LEN_PLOTTER)      :: auxname

         auxname = Name(1: len_trim(Name) - len(".HiOMesh")) // ".vtk"

         self % Name = trim(auxname)

         call Paraview_OpenFile ( self , isMesh = .true. ) 
!
!        Compute number of cells and points
!        ----------------------------------
         self % npoints = 0
         self % ncells  = 0
         do eID = 1 , mesh % no_of_elements
            self % npoints = self % npoints + (mesh % elements(eID) % spA % N + 1)**2
            self % ncells = self % ncells + (mesh % elements(eID) % spA % N)**2
         end do
!
         
         write ( self % fID , '(A)') "DATASET UNSTRUCTURED_GRID"

         call Paraview_WriteMesh( self , mesh )
  
         close ( self % fID )

      end subroutine Paraview_ExportMesh

      module subroutine Paraview_Export( self , mesh , Name) 
         use QuadMeshClass
         implicit none
         class(Paraview_t)         :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
         integer                 :: eID

         self % Name = trim(Name) // ".vtk"

         call Paraview_OpenFile ( self , isMesh = .false. ) 
!
!        Compute number of cells and points
!        ----------------------------------
         self % npoints = 0
         self % ncells  = 0
         do eID = 1 , mesh % no_of_elements
            self % npoints = self % npoints + (mesh % elements(eID) % spA % N + 1)**2
            self % ncells = self % ncells + (mesh % elements(eID) % spA % N)**2
         end do
!
!        Create the zone         
!        ---------------
         write ( self % fID , '(A)') "DATASET UNSTRUCTURED_GRID"
!
!        Write the mesh
!        --------------
         call Paraview_WriteMesh( self , mesh )
!
!        Write the variables
!        -------------------
         call Paraview_WriteVariables(self , mesh)
!
!        Close file
!        ----------
         close ( self % fID )
!
      end subroutine Paraview_Export

      subroutine Paraview_GatherVariables( self ) 
         use Setup_Class
         implicit none
         class(Paraview_t)               :: self
         integer                        :: pos
         character(len=STR_LEN_PLOTTER) :: auxstr
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

      end subroutine Paraview_gatherVariables

      subroutine Paraview_OpenFile( self , IsMesh ) 
         implicit none
         class(Paraview_t)        :: self
         logical                 :: IsMesh
         integer                 :: var

         open( newunit = self % fID , file = trim(self % Name) , status = "unknown" , action = "write" ) 

!
!        Print the header into the file
!        ------------------------------
         write( self % fID , '(A)') "# vtk DataFile Version 2.0"
         write( self % fID , '(A)') trim(self % Name)
         write( self % fID , '(A)') "ASCII"

      end subroutine Paraview_OpenFile

      subroutine Paraview_WriteMesh(self , mesh ) 
         use QuadMeshClass
         use Physics
         implicit none
         class(Paraview_t)        :: self
         class(QuadMesh_t)       :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: eID
         integer       :: iXi , iEta
         integer       :: var
         real(kind=RP) :: Q(1:NCONS)

         write ( self % fID , *)
         write ( self % fID , '(A,I0,A)') "POINTS " , self % npoints , " float"
   
                  
         do eID = 1 , mesh % no_of_elements
            associate ( N => mesh % elements(eID) % spA % N )
            do iEta = 0 , N
               do iXi = 0 , N
                  write( self % fID , '(ES17.10,1X,ES17.10,1X,ES17.10)') mesh % elements(eID) % x(iXi,iEta,IX) * RefValues % L &
                                                                                 , mesh % elements(eID) % x(iXi,iEta,IY) * RefValues % L &
                                                                                 , 0.0_RP  
               end do
            end do
            end associate
         end do

         write( self % fID , * )    ! One blank line
   
         write( self % fID , '(A,I0,1X,I0)' ) "CELLS ", self % ncells,5*self % ncells

         point_position = -1
         do eID = 1 , mesh % no_of_elements
            associate ( N => mesh % elements(eID) % spA % N )
            do iEta = 1 , N
               do iXi = 1 , N
                  write(self % fID , '(I0,1X,I0,1X,I0,1X,I0,1X,I0)')  4,pointPosition(iXi,iEta,N) + point_position
               end do
            end do
            point_position = point_position + (N+1)*(N+1)
            end associate
         end do

         write( self % fID , * )    ! One blank line
         write( self % fID , '(A,I0)' ) "CELL_TYPES ", self % ncells
         do eID = 1 , mesh % no_of_elements
            associate ( N => mesh % elements(eID) % spA % N )
            do iEta = 1 , N
               do iXi = 1 , N
                  write(self % fID , '(I0)')  9
               end do
            end do
            end associate
         end do

      end subroutine Paraview_WriteMesh

      subroutine Paraview_WriteVariables( self , mesh )
         use QuadMeshClass
         use Physics
         implicit none
         class(Paraview_t),   intent(in)     :: self
         class(QuadMesh_t),   intent(in)     :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iVar , eID , iXi , iEta
         real(kind=RP)  :: Q(NCONS)
         real(kind=RP)  :: dQ(NDIM,NCONS)

         write(self % fID,*)
         write(self % fID , '(A,I0)') "POINT_DATA " , self % npoints
         do iVar = 1 , self % no_of_variables

            select case ( trim( self % variables(iVar) ) )

               case ("rho")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rho float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % Q(iXi,iEta,IRHO) * refValues % rho
                     end do            ; end do
                     end associate
                  end do

               case ("rhou")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhou float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % Q(iXi,iEta,IRHOU) * refValues % rho * refValues % a
                     end do            ; end do
                     end associate
                  end do

               case ("rhov")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhov float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % Q(iXi,iEta,IRHOV) * refValues % rho * refValues % a
                     end do            ; end do
                     end associate
                  end do
 
               case ("rhoe")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhoe float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % Q(iXi,iEta,IRHOE) * refValues % p
                     end do            ; end do
                     end associate
                  end do

               case ("rhot")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhot float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % QDot(iXi,iEta,IRHO) * refValues % rho / refValues % tc
                     end do            ; end do
                     end associate
                  end do

               case ("rhout")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhout float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % QDot(iXi,iEta,IRHOU) * refValues % rho * refValues % a / refValues % tc
                     end do            ; end do
                     end associate
                  end do

               case ("rhovt")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhovt float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % QDot(iXi,iEta,IRHOV) * refValues % rho * refValues % a / refValues % tc
                     end do            ; end do
                     end associate
                  end do
 
               case ("rhoet")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS rhoet float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % QDot(iXi,iEta,IRHOE) * refValues % p / refValues % tc
                     end do            ; end do
                     end associate
                  end do

                case ("u")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS u float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % Q(iXi,iEta,IRHOU) / mesh % elements(eID) % Q(iXi,iEta,IRHO) * refValues % a
                     end do            ; end do
                     end associate
                  end do

               case ("v")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS v float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        write(self % fID , '(ES17.10)') mesh % elements(eID) % Q(iXi,iEta,IRHOV) / mesh % elements(eID) % Q(iXi,iEta,IRHO) * refValues % a
                     end do            ; end do
                     end associate
                  end do
 
               case ("p")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS p float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        write(self % fID , '(ES17.10)') getPressure(Q) * refValues % p
                     end do            ; end do
                     end associate
                  end do

               case ("Mach")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS Mach float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        write(self % fID , '(ES17.10)') ( sqrt(Q(IRHOU)*Q(IRHOU)+Q(IRHOV)*Q(IRHOV))/ Q(IRHO)) / getSoundSpeed(Q)
                     end do            ; end do
                     end associate
                  end do
 
               case ("s")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS s float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        write(self % fID , '(ES17.10)') getPressure(Q) * refValues % p / (Q(IRHO) * refValues % rho) ** ( Thermodynamics % gamma) 
                     end do            ; end do
                     end associate
                  end do
 
#ifdef NAVIER_STOKES

               case ("ux")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS ux float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        dQ = mesh % elements(eID) % dQ(iXi,iEta,:,:)
                        write(self % fID , '(ES17.10)') (- Q(IRHOU) / Q(IRHO) * dQ(IX,IRHO) + dQ(IX,IRHOU))/Q(IRHO)
                     end do            ; end do
                     end associate
                  end do

               case ("uy")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS uy float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        dQ = mesh % elements(eID) % dQ(iXi,iEta,:,:)
                        write(self % fID , '(ES17.10)') (- Q(IRHOU) / Q(IRHO) * dQ(IY,IRHO) + dQ(IY,IRHOU))/Q(IRHO)
                     end do            ; end do
                     end associate
                  end do
 
               case ("vx")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS vx float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        dQ = mesh % elements(eID) % dQ(iXi,iEta,:,:)
                        write(self % fID , '(ES17.10)') (- Q(IRHOV) / Q(IRHO) * dQ(IX,IRHO) + dQ(IX,IRHOV))/Q(IRHO)
                     end do            ; end do
                     end associate
                  end do

               case ("vy")
                  write(self % fID,*)
                  write(self % fID , '(A)') "SCALARS vy float"
                  write(self % fID , '(A)') "LOOKUP_TABLE default"
                  do eID = 1 , mesh % no_of_elements
                     associate ( N => mesh % elements(eID) % spA % N ) 
                     do iEta = 0 , N   ; do iXi = 0 , N
                        Q = mesh % elements(eID) % Q(iXi,iEta,:)
                        dQ = mesh % elements(eID) % dQ(iXi,iEta,:,:)
                        write(self % fID , '(ES17.10)') (- Q(IRHOV) / Q(IRHO) * dQ(IY,IRHO) + dQ(IY,IRHOV))/Q(IRHO)
                     end do            ; end do
                     end associate
                  end do

#endif

                  
               end select
            end do

      end subroutine Paraview_WriteVariables
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!                 AUXILIAR SUBROUTINES
!                 --------------------
!///////////////////////////////////////////////////////////////////////////////////////////
!
      function pointPosition(iXi , iEta , N) result( val )
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

      subroutine LinkedList_Destruct( self ) 
         implicit none
         class(LinkedList_t)      :: self
         class(Charlist), pointer :: current
         class(Charlist), pointer :: next
         integer                  :: i

         current => self % head

         do i = 1 , self % no_of_entries
            next => current % next
            deallocate( current )
            current => next

         end do


      end subroutine LinkedList_Destruct

end submodule Paraview  
