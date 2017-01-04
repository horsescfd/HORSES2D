module Tecplot
   use SMConstants

   private
   public 

   integer, parameter         :: STR_LEN_TECPLOT = 128

   type Tecplot_t
      integer        :: no_of_variables
      character(len=STR_LEN_TECPLOT)   :: variables(:)
      contains
         procedure      :: gatherVariables      => Tecplot_GatherVariables
   end type Tecplot_t

!
!  ========
   contains
!  ========
!
      subroutine TecPlot_Save( mesh ) 
         implicit none
         type(Tecplot_t)         :: tec 


      end subroutine


      subroutine TecPlot_GatherVariables( self ) 
         implicit none
         class(Tecplot_t)     :: self

      

      end subroutine TecPlot_gatherVariables

end module Tecplot
