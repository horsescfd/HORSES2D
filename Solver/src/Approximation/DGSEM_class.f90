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
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      DGSEM_class.f90
!      Created: 2015-09-26 12:05:17
!      REV:     2016-02-24 21:16:23
!      By: Juan MANZANERO
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!   *****************
    module DGSEM_class
!   *****************

!   ------------
!   Modules used
!   ------------
    use SMConstants
    use Physics
    use NodesAndWeights_class
    use QuadMeshClass
    use MeshFileClass
    use DGSpatialDiscretizationMethods
    use DGTimeIntegrator
    use Storage_module
    use DGBoundaryConditions
    use Plotter
    implicit none
!
#include "Defines.h"
!
    private
    public DGSEM_t , DGSEM_Initialize
!
!                                   *******************
    integer, parameter           :: STR_LEN_DGSEM = 128
!                                   *******************
!
!   ---------------------
!   DGSEM type DEFINITION
!   ---------------------
!
    type DGSEM_t
        type(QuadMesh_t)                    :: mesh
        type(NodalStorage)                  :: spA             ! Interpolation nodes and weights structure
        class(NodesAndWeights_t),   pointer :: spI => NULL()   ! Integration nodes and weights structure
        type(Storage_t)                     :: Storage
        class(BoundaryCondition_t), pointer :: BoundaryConditions(:)
        type(TimeIntegrator_t)              :: Integrator
        class(Plotter_t), allocatable       :: Plotter
        contains
            procedure :: Construct           => DGSEM_construct
            procedure :: SetInitialCondition => DGSEM_SetInitialCondition
            procedure :: Integrate           => DGSEM_Integrate
            procedure :: LoadRestartFile     => DGSEM_LoadRestartFile
            procedure :: Finalize            => DGSEM_Finalize
    end type DGSEM_t

!
!   ========
    contains
!   ========
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        function DGSEM_Initialize() result(sem)
            type(DGSEM_t)         :: sem
!
!           Initialize mesh
!           ---------------
            sem % mesh = InitializeMesh()
!
!           Initialize Spectral Approximation object
!           ----------------------------------------
            sem % spA  = newNodalStorage()
!
!           Initialize Storage
!           ------------------ 
            sem % Storage = newStorage()
             
        end function DGSEM_Initialize
         
        subroutine DGSEM_construct( self )
            use Setup_class
            use QuadElementClass
            use StopwatchClass
            implicit none
            class(DGSEM_t)   :: self
!
!           ---------------
!           Local variables
!           ---------------
!
            type(MeshFile_t) :: meshFile
            integer          :: totalPolynomialOrder
!
!           Create an event to measure the preprocessing time
!           -------------------------------------------------
            call Stopwatch % CreateNewEvent("Preprocessing")
            call Stopwatch % Start("Preprocessing")
!
!           Read the mesh file
!           ------------------
            call meshFile % Read
!
!           Allocate memory for the solution, time derivative, and gradients
!           ----------------------------------------------------------------
            totalPolynomialOrder = meshFile % cumulativePolynomialOrder( meshFile % no_of_elements )
            call self % Storage % AllocateMemory( totalPolynomialOrder )
!
!           Construct the spectral Integration class if Over-Integration is selected
!           ------------------------------------------------------------------------
            if (Setup % inviscid_discretization .eq. "Over-Integration") then
               allocate( self % spI )
               call self % spI % init( Setup % integration_points , Setup % nodes )
            else
               self % spI => NULL()
            end if
!
!           Construct the spectral element mesh object
!           ------------------------------------------
            call self % mesh % constructFromFile(meshFile , self % spA , self % Storage , self % spI)
!
!           Prepare the spectral Approximation structures generated for Over-Integration
!           ----------------------------------------------------------------------------
            if (Setup % inviscid_discretization .eq. "Over-Integration") then
!          
!              Set the interpolation matrices
!              ------------------------------
               call self % spA % computeInterpolationMatrices( self % spI )

            end if
!
!           Construct the domain zones
!           --------------------------              
            call self % mesh % ConstructZones( meshFile )
!
!           Initialize Inviscid and Viscous discretization methods
!           ------------------------------------------------------
            call DGSpatial_Initialization() 
!
!           Construct plotter and Export the mesh            
!           -------------------------------------
            call ConstructPlotter( self % Plotter )
            call self % Plotter % ExportMesh( self % mesh , Setup % mesh_file )   
!
!           Destruct the mesh file object
!           -----------------------------
            call meshFile % Destruct
!
!           Describe the problem file
!           -------------------------
            call DescribeProblemFile()
!
!           Set the initial condition to all flow variables
!           -----------------------------------------------
            call self % SetInitialCondition()
            call self % Plotter % Export( self % mesh , './RESULTS/InitialCondition')     
!
!           Stop the preprocessing event time
!           ---------------------------------
            call Stopwatch % Pause("Preprocessing")
            
        end subroutine DGSEM_construct
            
        subroutine DGSEM_SetInitialCondition( self , verbose )
!
!           *********************************************************
!                 This routine sets the initial condition from one 
!              previously selected from the library, or from the
!              Restart file *.HiORst
!           *********************************************************
!              
            use InitialConditions
            use Setup_class
            implicit none
            class(DGSEM_t)                   :: self
            logical, optional                :: verbose

            if ( Setup % IC .eq. "Restart" ) then
               call self % loadRestartFile( trim ( Setup % Restart_file ) ) 
            else
               call self % mesh % SetInitialCondition ()
            end if
!
!           Describe the initial condition
!           ------------------------------
            if ( present ( verbose ) ) then
               if ( verbose ) then
                  call InitialCondition_Describe
   
               end if
            else
               call InitialCondition_Describe

            end if

            call DGSpatial_ComputeTimeDerivative( self % mesh , Setup % initialTime )

        end subroutine DGSEM_SetInitialCondition
      
        subroutine DGSEM_Integrate( self )
!
!           *********************************************************
!                 This routine constructs the Time Integrator,
!              and performs the time integration of the problem.
!           *********************************************************
!
            use Setup_Class
            implicit none
            class(DGSEM_t)                   :: self
            character(len=STR_LEN_DGSEM)     :: solutionpltName
            integer                          :: pos
!
!           Construct Time Integrator
!           -------------------------
            self % Integrator = NewTimeIntegrator( self % mesh )
            call self % Integrator % Describe()
!
!           Integrate
!           ---------
            call self % Integrator % Integrate( self % mesh , self % Storage)
!
!           Save the solution file
!           ----------------------
            pos = index( Setup % solution_file , '.HiORst')

            if ( pos .gt. 0 ) then
               solutionpltName = Setup % solution_file(1:pos-1) 
            end if
   
            call self % Plotter % Export ( self % mesh , trim(solutionpltname))  

        end subroutine DGSEM_Integrate
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!            
!              AUXILIAR FUNCTIONS
!              ------------------ 
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
        subroutine DGSEM_loadRestartFile( self , fileName ) 
!     
!           ************************************************************************
!                 Loads a previous state from a Restart *.HiORst file
!           ************************************************************************
!
            use NetCDFInterface
            use Setup_class
            implicit none
            class(DGSEM_t)                :: self
            character(len=*)              :: fileName
            real(kind=RP), allocatable    :: t(:)
            integer, allocatable          :: iter(:)
            real(kind=RP), allocatable    :: Q(:)
!
!           Get the current time and iteration
!           ----------------------------------
            call NetCDF_getVariable ( trim ( fileName )  , "t"    , t    ) 
            call NetCDF_getVariable ( trim ( fileName )  , "iter" , iter ) 
            call Setup % SetInitialTime ( t(1) , iter(1)) 
!
!           Load the state vector
!           ---------------------
            call NetCDF_getVariable ( trim ( fileName )  , "Q"    , Q    ) 
            self % Storage % Q = Q
      
            deallocate( t , Q ) 
   
        end subroutine DGSEM_loadRestartFile

        function DGSEM_Finalize( self ) result ( exit_code )
            use SMConstants
            use Physics
            use Setup_class
            use DGTimeIntegrator
            use MonitorsClass
            use Headers
            use StopwatchClass
            implicit none
            class(DGSEM_t)       :: self
            integer              :: exit_code
!
!           ----------
!           Interfaces            
!           ----------
!
            interface
               function Finalize( sem_ , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ , Monitors_ ) result(exit_code)
                  use SMConstants
                  use Setup_class
                  use QuadMeshClass
                  use QuadElementClass
                  use Physics
                  use MonitorsClass
                  import DGSEM_t
                  implicit none
                  class(DGSEM_t)                      :: sem_
                  class(Thermodynamics_t), intent(in) :: thermodynamics_
                  class(Setup_t),          intent(in) :: Setup_
                  class(RefValues_t),      intent(in) :: refValues_
                  class(Dimensionless_t),  intent(in) :: dimensionless_
                  class(Monitor_t),        intent(in) :: Monitors_
                  integer                             :: exit_code
               end function Finalize
            end interface

            write(STD_OUT,*)
            call Section_header("Solution analysis")

            exit_code = Finalize ( self , Thermodynamics , Setup , refValues , dimensionless , Monitors) 

            write(STD_OUT,*)
            call Section_header("Simulation statistics")
            write(STD_OUT,*)
            write(STD_OUT,'(30X,A,A20,ES10.3,A)') "-> ","Preprocessing time: ", Stopwatch % ElapsedTime("Preprocessing") , " (s)."
            write(STD_OUT,'(30X,A,A20,ES10.3,A)') "-> ","Simulation time: ", Stopwatch % ElapsedTime("Simulation") , " (s)."
            write(STD_OUT,'(30X,A,A20,ES10.3,A)') "-> ","Solver perfomance: ", &
                  1.0e06_RP * Stopwatch % ElapsedTime("Simulation") / self % Integrator % no_of_iterations / size(self % Storage % Q), " (s/iter 1MDOF)."

        end function DGSEM_Finalize

        subroutine DescribeProblemFile()
            use Headers
            implicit none
            interface
               function getProblemFileName() result (name)
                  implicit none
                  character(len = LINE_LENGTH)  :: name
               end function getProblemFileName
            end interface

            write(STD_OUT,'(/)') 
            call Section_header("Problem file")
            write(STD_OUT,*)
            write(STD_OUT , '(30X,A,A,A)') "-> ", "Linked problem file: " , trim(getProblemFileName())

         end subroutine DescribeProblemFile

   end module DGSEM_class
