module Plotter

   use SMConstants

#include "Defines.h"

   private
   public   Plotter_t , ConstructPlotter

   integer, parameter         :: STR_LEN_PLOTTER = 128

   type Plotter_t
      contains
         procedure   :: Construct   => Plotter_Initialization
         procedure   :: Export      => Plotter_Export
         procedure   :: ExportMesh  => Plotter_ExportMesh
   end type Plotter_t

   type, extends(Plotter_t)   ::  Tecplot_t
      integer        :: no_of_variables = 0
      integer        :: fID
      character(len=STR_LEN_PLOTTER), allocatable   :: variables(:)
      character(len=STR_LEN_PLOTTER)                :: Name
      contains
         procedure      :: Construct            => Tecplot_Initialization
         procedure      :: Export               => Tecplot_Export
         procedure      :: ExportMesh           => Tecplot_ExportMesh
   end type Tecplot_t

   type, extends(Plotter_t)   ::  Paraview_t
      integer        :: no_of_variables = 0
      integer        :: npoints
      integer        :: ncells
      integer        :: fID
      character(len=STR_LEN_PLOTTER), allocatable   :: variables(:)
      character(len=STR_LEN_PLOTTER)                :: Name
      contains
         procedure      :: Construct            => Paraview_Initialization
         procedure      :: Export               => Paraview_Export
         procedure      :: ExportMesh           => Paraview_ExportMesh
   end type Paraview_t

!
!  ---------------
!  Interface block 
!  ---------------
!
   interface
      module subroutine Tecplot_Initialization ( self )
         use QuadMeshClass
         class(Tecplot_t)     :: self
      end subroutine Tecplot_Initialization
         
      module subroutine Tecplot_ExportMesh( self , mesh , Name ) 
         use QuadMeshClass
         implicit none
         class(Tecplot_t)        :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
      end subroutine Tecplot_ExportMesh

      module subroutine TecPlot_Export( self , mesh , Name) 
         use QuadMeshClass
         implicit none
         class(Tecplot_t)         :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
      end subroutine TecPlot_Export
   end interface

   interface
      module subroutine Paraview_Initialization ( self )
         use QuadMeshClass
         class(Paraview_t)     :: self
      end subroutine Paraview_Initialization
         
      module subroutine Paraview_ExportMesh( self , mesh , Name ) 
         use QuadMeshClass
         implicit none
         class(Paraview_t)        :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
      end subroutine Paraview_ExportMesh

      module subroutine Paraview_Export( self , mesh , Name) 
         use QuadMeshClass
         implicit none
         class(Paraview_t)         :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
      end subroutine Paraview_Export
   end interface
!
!  ========  
   contains
!  ========  
!
      subroutine ConstructPlotter( self )
         use Setup_class
         implicit none
         class(Plotter_t), allocatable       :: self

!
!        TODO: Select format in case file
!        --------------------------------
         if ( (trim(Setup % exportFormat) .eq. "Tecplot") .or. (trim(Setup % exportFormat) .eq. "plt") ) then
            allocate ( Tecplot_t    :: self )

         elseif (( trim(Setup % exportFormat) .eq. "Paraview") .or. (trim(Setup % exportFormat) .eq. "vtk" ) ) then
            allocate ( Paraview_t    :: self )

         else
            print*, "Plotter was not loaded"
            errorMessage(STD_OUT)
            stop "Stopped"

         end if

         call self % Construct

      end subroutine ConstructPlotter

      subroutine Plotter_Initialization ( self )
         use QuadMeshClass
         class(Plotter_t)     :: self
      end subroutine Plotter_Initialization
         
      subroutine Plotter_ExportMesh( self , mesh , Name ) 
         use QuadMeshClass
         implicit none
         class(Plotter_t)        :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
      end subroutine Plotter_ExportMesh

      subroutine Plotter_Export( self , mesh , Name) 
         use QuadMeshClass
         implicit none
         class(Plotter_t)         :: self
         class(QuadMesh_t)       :: mesh
         character(len=*)        :: Name
      end subroutine Plotter_Export

end module Plotter
