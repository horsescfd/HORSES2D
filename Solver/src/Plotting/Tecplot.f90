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
submodule (Plotter) Tecplot
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
!
!
!  ========
   contains
!  ========
!
      module subroutine Tecplot_Initialization ( self )
         use QuadMeshClass
         class(Tecplot_t)     :: self

         call Tecplot_GatherVariables( self )

      end subroutine Tecplot_Initialization
         
      module subroutine Tecplot_ExportMesh( self , mesh , Name ) 
         use QuadMeshClass
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
         integer                 :: eID
         character(len=STR_LEN_PLOTTER)      :: auxname

         auxname = Name(1: len_trim(Name) - len(".HiOMesh")) // ".plt"

         self % Name = trim(auxname)

         call Tecplot_OpenFile ( self , isMesh = .true. ) 
      
         do eID = 1 , mesh % no_of_elements
            call Tecplot_NewMeshZone( self , mesh , eID ) 
         end do
   
         close ( self % fID )

      end subroutine Tecplot_ExportMesh

      module subroutine TecPlot_Export( self , mesh , Name) 
         use QuadMeshClass
         implicit none
         class(Tecplot_t)         :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
         integer                 :: eID

         self % Name = trim(Name) // ".plt"

         call TecPlot_OpenFile ( self , isMesh = .false. )

         do eID = 1 , mesh % no_of_elements
            call Tecplot_NewZone( self , mesh , eID )
         end do

         close( self % fID )

      end subroutine TecPlot_Export

      subroutine TecPlot_GatherVariables( self ) 
         use Setup_Class
         implicit none
         class(Tecplot_t)               :: self
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

      end subroutine TecPlot_gatherVariables

      subroutine Tecplot_OpenFile( self , IsMesh ) 
         implicit none
         class(Tecplot_t)        :: self
         logical                 :: IsMesh
         integer                 :: var

         open( newunit = self % fID , file = trim(self % Name) , status = "unknown" , action = "write" ) 

!
!        Print the header into the file
!        ------------------------------
         write( self % fID , '(A,A,A)') 'TITLE = "',trim(self % Name),'"'
         write( self % fID , '(A)' , advance="no") 'VARIABLES = "X" "Y" "Z" '

         if ( .not. IsMesh ) then
            do var = 1 , self % no_of_variables
               write( self % fID , '(A,A,A)' , advance="no" ) '"',trim(self % variables(var)),'" '
            end do
         end if
         write( self % fID , * )

      end subroutine Tecplot_OpenFile

      subroutine Tecplot_NewMeshZone(self , mesh , eID ) 
         use QuadMeshClass
         use Physics
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: eID 
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                 :: iXi , iEta
         integer                 :: var
         real(kind=RP)              :: Q(1:NCONS)

         associate ( N => mesh % elements(eID) % spA % N )
!         
!        New header
!        ----------
         write( self % fID , '(A,I0,A)' , advance="no" ) "ZONE N=",(N+1)*(N+1),", "
         write( self % fID , '(A,I0,A)' , advance="no" ) "E=",(N)*(N),", "
         write( self % fID , '(A)'                     ) "DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"

         do iEta = 0 , N
            do iXi = 0 , N
               write( self % fID , '(ES17.10,1X,ES17.10,1X,ES17.10)',advance="no") mesh % elements(eID) % x(iXi,iEta,IX) * RefValues % L &
                                                                              , mesh % elements(eID) % x(iXi,iEta,IY) * RefValues % L &
                                                                              , 0.0_RP  
!
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

      end subroutine Tecplot_NewMeshZone

      subroutine Tecplot_NewZone( self , mesh , eID , zoneType) 
         use QuadMeshClass
         use Physics
         use Setup_Class
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: eID 
         character(len=*), optional :: zoneType
         character(len=STR_LEN_PLOTTER)           :: zType

         associate ( N => mesh % elements(eID) % spA % N )

         if ( present ( zoneType ) ) then
            zType = zoneType

         else
            zType = trim( Setup % outputType ) 
   
         end if

         if ( trim(zType) .eq. "DGSEM" ) then
            call Tecplot_StandardZone( self , mesh , eID )
         elseif ( trim(zType) .eq. "Interpolated") then
            call Tecplot_InterpolatedZone( self , mesh , eID ) 
         else
            print*, "Unknown output type " , trim(Setup % outputType)
            print*, "Options available are: "
            print*, "   * DGSEM"
            print*, "   * Interpolated"
            stop "Stopped."
         end if

         end associate

      end subroutine  Tecplot_NewZone

      subroutine Tecplot_StandardZone(self , mesh , eID ) 
         use QuadMeshClass
         use Physics
         use DGSpatialDiscretizationMethods
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: eID 
!        --------------------------------------------------------------------------
         real(kind=RP), pointer  :: rho (:,:) , rhou (:,:) , rhov (:,:) , rhoe (:,:)
         real(kind=RP), pointer  :: rhot(:,:) , rhout(:,:) , rhovt(:,:) , rhoet(:,:)
#ifdef NAVIER_STOKES
         real(kind=RP), target   :: du(0:mesh % elements(eID) % spA % N , 0:mesh % elements(eID) % spA % N , 1:NDIM , 1:NDIM)
         real(kind=RP), pointer  :: ux(:,:) , uy(:,:) , vx(:,:) , vy(:,:)
#endif
         integer                 :: iXi , iEta
         integer                 :: var
         real(kind=RP)              :: Q(1:NCONS)

         associate ( N => mesh % elements(eID) % spA % N )
!         
!        New header
!        ----------
         write( self % fID , '(A,I0,A)' , advance="no" ) "ZONE N=",(N+1)*(N+1),", "
         write( self % fID , '(A,I0,A)' , advance="no" ) "E=",(N)*(N),", "
         write( self % fID , '(A)'                     ) "DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"
!
!        Point to the quantities
!        -----------------------
         rho(0:,0:)  => mesh % elements(eID) % Q(IRHO,0:,0:) 
         rhou(0:,0:) => mesh % elements(eID) % Q(IRHOU,0:,0:)
         rhov(0:,0:) => mesh % elements(eID) % Q(IRHOV,0:,0:)
         rhoe(0:,0:) => mesh % elements(eID) % Q(IRHOE,0:,0:)
         rhot(0:,0:)  => mesh % elements(eID) % QDot(IRHO,0:,0:) 
         rhout(0:,0:) => mesh % elements(eID) % QDot(IRHOU,0:,0:)
         rhovt(0:,0:) => mesh % elements(eID) % QDot(IRHOV,0:,0:)
         rhoet(0:,0:) => mesh % elements(eID) % QDot(IRHOE,0:,0:)

#ifdef NAVIER_STOKES
   if ( self % no_of_variables .ne. 0 ) then
         du = getStrainTensor( N , mesh % elements(eID) % Q , mesh % elements(eID) % dQ )
         ux(0:,0:) => du(0:,0:,IX,IX)
         uy(0:,0:) => du(0:,0:,IY,IX)
         vx(0:,0:) => du(0:,0:,IX,IY)
         vy(0:,0:) => du(0:,0:,IY,IY)
   end if
#endif


         do iEta = 0 , N
            do iXi = 0 , N
               write( self % fID , '(ES17.10,1X,ES17.10,1X,ES17.10)',advance="no") mesh % elements(eID) % x(iXi,iEta,IX) * RefValues % L &
                                                                              , mesh % elements(eID) % x(iXi,iEta,IY) * RefValues % L &
                                                                              , 0.0_RP  
!
!              Save quantities
!              ---------------
               do var = 1 , self % no_of_variables

                  select case ( trim( self % variables(var) ) )
                     case ("rho")
                        write(self % fID,'(1X,ES17.10)',advance="no") rho(iXi,iEta) * refValues % rho

                     case ("rhou")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhou(iXi,iEta) * refValues % rho * refValues % a

                     case ("rhov")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhov(iXi,iEta) * refValues % rho * refValues % a

                     case ("rhoe")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhoe(iXi,iEta) * refValues % p

                     case ("rhot")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhot(iXi,iEta) * refValues % rho / refValues % tc

                     case ("rhout")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhout(iXi,iEta) * refValues % rho * refValues % a / refValues % tc

                     case ("rhovt")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhovt(iXi,iEta) * refValues % rho * refValues % a / refValues % tc

                     case ("rhoet")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhoet(iXi,iEta) * refValues % p / refValues % tc

                     case ("u")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhou(iXi,iEta)/rho(iXi,iEta) * refValues % a

                     case ("v")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhov(iXi,iEta)/rho(iXi,iEta) * refValues % a
   
                     case ("p")
                        Q = [rho(iXi,iEta) , rhou(iXi,iEta) , rhov(iXi,iEta) , rhoe(iXi,iEta) ]
                        write(self % fID,'(1X,ES17.10)',advance="no") getPressure( Q ) * refValues % p
      
                     case ("Mach")
                        Q = [rho(iXi,iEta) , rhou(iXi,iEta) , rhov(iXi,iEta) , rhoe(iXi,iEta) ]
                        write(self % fID,'(1X,ES17.10)',advance="no") sqrt(rhou(iXi,iEta)*rhou(iXi,iEta)+rhov(iXi,iEta)*rhov(iXi,iEta))/rho(iXi,iEta)/ getSoundSpeed( Q )

                     case ("s")
                        write(self % fID,'(1X,ES17.10)',advance="no") Thermodynamics % gm1 * ( rhoe(iXi,iEta) - &
                                                  0.5*rhou(iXi,iEta)*rhou(iXi,iEta)/rho(iXi,iEta) - 0.5*rhov(iXi,iEta)*rhov(iXi,iEta)/rho(iXi,iEta) ) &
                                                    / (rho(iXi,iEta) * refValues % rho)**(Thermodynamics % gamma) * refValues % p
#ifdef NAVIER_STOKES
                      case ( "ux" ) 
                        write(self % fID,'(1X,ES17.10)',advance="no") ux(iXi,iEta) * refValues % a

                      case ( "uy" ) 
                        write(self % fID,'(1X,ES17.10)',advance="no") uy(iXi,iEta) * refValues % a
   
                      case ( "vx" ) 
                        write(self % fID,'(1X,ES17.10)',advance="no") ux(iXi,iEta) * refValues % a

                      case ( "vy" )
                        write(self % fID,'(1X,ES17.10)',advance="no") vy(iXi,iEta) * refValues % a

                     case ("vort")
                        write(self % fID,'(1X,ES17.10)',advance="no") ( vx(iXi,iEta) - uy(iXi,iEta) ) * refValues % a / refValues % L
            
                     case ("muart")
                        write(self % fID,'(1X,ES17.10)',advance="no") ArtificialDissipation % ElementViscosity( ArtificialDissipation , mesh % elements(eID) )
#endif

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

      end subroutine Tecplot_StandardZone

      subroutine Tecplot_InterpolatedZone(self , mesh , eID ) 
         use QuadMeshClass
         use Physics
         use Setup_Class
         use MatrixOperations
         use InterpolationAndDerivatives
         use DGSpatialDiscretizationMethods
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         integer                 :: eID 
!        --------------------------------------------------------------------------
         real(kind=RP), pointer     :: rhoDG  (:,:) ,  rhouDG  (:,:) ,  rhovDG  (:,:) ,  rhoeDG  (:,:) 
         real(kind=RP), pointer     :: rhotDG (:,:) ,  rhoutDG (:,:) ,  rhovtDG (:,:) ,  rhoetDG (:,:) 
         real(kind=RP), pointer     :: rho    (:,:) ,  rhou    (:,:) ,  rhov    (:,:) ,  rhoe    (:,:) 
         real(kind=RP), pointer     :: rhot   (:,:) ,  rhout   (:,:) ,  rhovt   (:,:) ,  rhoet   (:,:) 
#ifdef NAVIER_STOKES
         real(kind=RP), target      :: du(0:mesh % elements(eID) % spA % N , 0:mesh % elements(eID) % spA % N , 1:NDIM , 1:NDIM)
         real(kind=RP), pointer     :: uxDG   (:,:) ,  uyDG    (:,:) ,  vxDG    (:,:) ,  vyDG    (:,:)
         real(kind=RP), pointer     :: ux     (:,:) ,  uy      (:,:) ,  vx      (:,:) ,  vy      (:,:)
#endif
         real(kind=RP), allocatable :: xi(:) , T(:,:) , x(:)
         integer                    :: iXi , iEta
         integer                    :: Nout
         integer                    :: var
         real(kind=RP)              :: Q(1:NCONS)

         associate ( N => mesh % elements(eID) % spA % N , spA => mesh % elements(eID) % spA)
!
!        Construct the interpolation framework
!        -------------------------------------
         Nout = Setup % no_of_plotPoints - 1

!         
!        New header
!        ----------
         write( self % fID , '(A,I0,A)' , advance="no" ) "ZONE N=",(Nout+1)*(Nout+1),", "
         write( self % fID , '(A,I0,A)' , advance="no" ) "E=",(Nout)*(Nout),", "
         write( self % fID , '(A)'                     ) "DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"

         allocate( xi ( 0 : Nout ) , T (0 : Nout , 0 : N ) , x(NDIM)) 

         xi = reshape( (/( (1.0_RP * iXi) / Nout , iXi = 0 , Nout )/) , (/Nout + 1/) )
         call PolynomialInterpolationMatrix( N , Nout , spA % xi , spA % wb , xi , T )
!
!        Point to the quantities
!        -----------------------
         rhoDG(0:,0:)  => mesh % elements(eID) % Q(IRHO,0:,0:) 
         rhouDG(0:,0:) => mesh % elements(eID) % Q(IRHOU,0:,0:)
         rhovDG(0:,0:) => mesh % elements(eID) % Q(IRHOV,0:,0:)
         rhoeDG(0:,0:) => mesh % elements(eID) % Q(IRHOE,0:,0:)
         rhotDG(0:,0:)  => mesh % elements(eID) % QDot(IRHO,0:,0:) 
         rhoutDG(0:,0:) => mesh % elements(eID) % QDot(IRHOU,0:,0:)
         rhovtDG(0:,0:) => mesh % elements(eID) % QDot(IRHOV,0:,0:)
         rhoetDG(0:,0:) => mesh % elements(eID) % QDot(IRHOE,0:,0:)

#ifdef NAVIER_STOKES
      if ( self % no_of_variables .ne. 0 ) then
         du = getStrainTensor( N , mesh % elements(eID) % Q , mesh % elements(eID) % dQ )
         uxDG(0:,0:) => du(0:,0:,IX,IX)
         uyDG(0:,0:) => du(0:,0:,IY,IX)
         vxDG(0:,0:) => du(0:,0:,IX,IY)
         vyDG(0:,0:) => du(0:,0:,IY,IY)
      end if
#endif
!
!        Obtain the interpolation to a new set of equispaced points
!        ----------------------------------------------------------
         allocate  (  rho(0:Nout  , 0:Nout) , rhou(0:Nout  , 0:Nout) , rhov(0:Nout  , 0:Nout) , rhoe(0:Nout  , 0:Nout) )
         allocate  (  rhot(0:Nout , 0:Nout) , rhout(0:Nout , 0:Nout) , rhovt(0:Nout , 0:Nout) , rhoet(0:Nout , 0:Nout) )
#ifdef NAVIER_STOKES
         allocate  (  ux(0:Nout , 0:Nout) , uy(0:Nout , 0:Nout) , vx(0:Nout , 0:Nout) , vy(0:Nout,0:Nout) )
#endif
         call TripleMatrixProduct ( T , rhoDG   , T , rho   , trC = .true. ) 
         call TripleMatrixProduct ( T , rhouDG  , T , rhou  , trC = .true. ) 
         call TripleMatrixProduct ( T , rhovDG  , T , rhov  , trC = .true. ) 
         call TripleMatrixProduct ( T , rhoeDG  , T , rhoe  , trC = .true. ) 
         call TripleMatrixProduct ( T , rhotDG  , T , rhot  , trC = .true. ) 
         call TripleMatrixProduct ( T , rhoutDG , T , rhout , trC = .true. ) 
         call TripleMatrixProduct ( T , rhovtDG , T , rhovt , trC = .true. ) 
         call TripleMatrixProduct ( T , rhoetDG , T , rhoet , trC = .true. ) 
#ifdef NAVIER_STOKES
         call TripleMatrixProduct ( T , uxDG    , T , ux    , trC = .true. ) 
         call TripleMatrixProduct ( T , uyDG    , T , uy    , trC = .true. ) 
         call TripleMatrixProduct ( T , vxDG    , T , vx    , trC = .true. ) 
         call TripleMatrixProduct ( T , vyDG    , T , vy    , trC = .true. ) 
#endif

         do iEta = 0 , Nout
            do iXi = 0 , Nout

               x = mesh % elements(eID) % compute_X ( xi(iXi) , xi(iEta) )
               write( self % fID , '(ES17.10,1X,ES17.10,1X,ES17.10)',advance="no") x(IX)  * RefValues % L &
                                                                              , x(IY) * RefValues % L &
                                                                              , 0.0_RP  
!
!              Save quantities
!              ---------------
               do var = 1 , self % no_of_variables

                  select case ( trim( self % variables(var) ) )
                     case ("rho")
                        write(self % fID,'(1X,ES17.10)',advance="no") rho(iXi,iEta) * refValues % rho

                     case ("rhou")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhou(iXi,iEta) * refValues % rho * refValues % a

                     case ("rhov")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhov(iXi,iEta) * refValues % rho * refValues % a

                     case ("rhoe")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhoe(iXi,iEta) * refValues % p

                     case ("rhot")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhot(iXi,iEta) * refValues % rho / refValues % tc

                     case ("rhout")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhout(iXi,iEta) * refValues % rho * refValues % a / refValues % tc

                     case ("rhovt")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhovt(iXi,iEta) * refValues % rho * refValues % a / refValues % tc

                     case ("rhoet")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhoet(iXi,iEta) * refValues % p / refValues % tc

                     case ("u")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhou(iXi,iEta)/rho(iXi,iEta) * refValues % a

                     case ("v")
                        write(self % fID,'(1X,ES17.10)',advance="no") rhov(iXi,iEta)/rho(iXi,iEta) * refValues % a
   
                     case ("p")
                        Q = [rho(iXi,iEta) , rhou(iXi,iEta) , rhov(iXi,iEta) , rhoe(iXi,iEta) ]
                        write(self % fID,'(1X,ES17.10)',advance="no") getPressure( Q ) * refValues % p
      
                     case ("Mach")
                        Q = [rho(iXi,iEta) , rhou(iXi,iEta) , rhov(iXi,iEta) , rhoe(iXi,iEta) ]
                        write(self % fID,'(1X,ES17.10)',advance="no") sqrt(rhou(iXi,iEta)*rhou(iXi,iEta)+rhov(iXi,iEta)*rhov(iXi,iEta))/rho(iXi,iEta)/ getSoundSpeed( Q )

                     case ("s")
                        write(self % fID,'(1X,ES17.10)',advance="no") Thermodynamics % gm1 * ( rhoe(iXi,iEta) - &
                                                  0.5*rhou(iXi,iEta)*rhou(iXi,iEta)/rho(iXi,iEta) - 0.5*rhov(iXi,iEta)*rhov(iXi,iEta)/rho(iXi,iEta) ) &
                                                    / (rho(iXi,iEta) * refValues % rho)**(Thermodynamics % gamma) * refValues % p
#ifdef NAVIER_STOKES
                      case ( "ux" ) 
                        write(self % fID,'(1X,ES17.10)',advance="no") ux(iXi,iEta) * refValues % a

                      case ( "uy" ) 
                        write(self % fID,'(1X,ES17.10)',advance="no") uy(iXi,iEta) * refValues % a
   
                      case ( "vx" ) 
                        write(self % fID,'(1X,ES17.10)',advance="no") ux(iXi,iEta) * refValues % a

                      case ( "vy" )
                        write(self % fID,'(1X,ES17.10)',advance="no") vy(iXi,iEta) * refValues % a

                     case ("vort")
                        write(self % fID,'(1X,ES17.10)',advance="no") ( vx(iXi,iEta) - uy(iXi,iEta) ) * refValues % a / refValues % L

                     case ("muart")
                        write(self % fID,'(1X,ES17.10)',advance="no") ArtificialDissipation % ElementViscosity( ArtificialDissipation , mesh % elements(eID) )
#endif

                  end select                        

               end do

!              Jump to next line
!              -----------------
               write( self % fID , *)

            end do
         end do

         write( self % fID , * )    ! One blank line

         do iEta = 1 , Nout
            do iXi = 1 , Nout
               write(self % fID , '(I0,1X,I0,1X,I0,1X,I0)')  pointPosition(iXi,iEta,Nout)
            end do
         end do

         end associate

      end subroutine Tecplot_InterpolatedZone

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

end submodule Tecplot  
