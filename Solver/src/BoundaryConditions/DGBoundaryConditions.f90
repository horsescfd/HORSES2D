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
!///////////////////////////////////////////////////////////////////////////////////////////
!
!  Boundary conditions file
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
module DGBoundaryConditions
   use SMConstants
   use Physics
   use QuadElementClass
   use ParamfileIO
   use InitialConditions
   implicit none

#include "Defines.h"

   private
   public BoundaryCondition_t , Construct
   public PeriodicBC_t , DirichletBC_t , FarfieldBC_t , EulerWall_t , PressureOutletBC_t , PressureInletBC_t
   public RiemannBC_t , newDirichletBC_t
   public BoundaryConditions_SetInitialCondition
!
!  ********************************************************************
   integer, parameter         :: STR_LEN_BC                       = 128
   integer, parameter         :: SPECIFY_SPEED                    = 1
   integer, parameter         :: SPECIFY_TOTAL_PRESSURE           = 2
   integer, parameter         :: ISOTHERMAL_WALL                  = 1
   integer, parameter         :: ADIABATIC_WALL                   = 2
   integer, parameter         :: REFLECTIVE_OUTFLOW               = 1
   integer, parameter         :: PARTIALLY_NON_REFLECTIVE_OUTFLOW = 2
   integer, parameter         :: STATIC_PRESSURE                  = 1
   integer, parameter         :: TOTAL_PRESSURE                   = 2
!  ********************************************************************
!
!
!  **********************************************
!  Boundary condition function abstract procedure
!  **********************************************
!
   abstract interface 
      function BoundaryStateFunction(x , time ) result (state)
         use SMConstants
         implicit none
         real(kind=RP), intent(in)           :: x(NDIM)
         real(kind=RP), intent(in)           :: time
         real(kind=RP)                       :: state(NCONS)
      end function BoundaryStateFunction
   end interface
!
!  **********************************
!  Base class for boundary conditions
!  **********************************
!
   type BoundaryCondition_t
      character(len=STR_LEN_BC)                         :: Name
      integer                                           :: marker
      integer                                           :: BCType
      integer(kind=1)                                   :: WeakType
      character(len=STR_LEN_BC)                         :: RiemannSolverName
      character(len=STR_LEN_BC)                         :: BoundaryStateName
      procedure(RiemannSolverFunction), pointer, nopass :: RiemannSolver => NULL()
      procedure(BoundaryStateFunction), pointer, nopass :: BoundaryState => NULL()
      contains
         procedure :: Construct        => BaseClass_Construct
         procedure :: SetRiemannSolver => BoundaryConditions_SetRiemannSolver
         procedure :: Associate        => BaseClass_Associate
         procedure :: UpdateSolution   => BaseClass_UpdateSolution
         procedure :: Describe         => BaseClass_Describe
   end type BoundaryCondition_t
!
!  *********************************
!  Periodic boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: PeriodicBC_t
      logical                          :: associated
      logical                          :: invert
      integer                          :: direction
      integer                          :: connected_marker
      contains
         procedure ::     Construct => PeriodicBC_Construct
         procedure ::     Describe  => PeriodicBC_Describe
   end type PeriodicBC_t
!
!  **********************************
!  Dirichlet boundary condition class
!  **********************************
!
   type, extends(BoundaryCondition_t)           :: DirichletBC_t
      real(kind=RP), dimension(NCONS)       :: q
      contains
         procedure ::      Construct => DirichletBC_Construct
         procedure ::      Associate => DirichletBC_Associate
         procedure ::      Describe  => DirichletBC_Describe
         procedure :: UpdateSolution => DirichletBC_UpdateSolution
   end type DirichletBC_t
!
!  *********************************
!  Farfield boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: FarfieldBC_t
      real(kind=RP), dimension(NCONS)       :: q
      contains
         procedure ::      Construct => FarfieldBC_Construct
         procedure ::      Associate => FarfieldBC_Associate
         procedure ::      UpdateSolution    => FarfieldBC_UpdateSolution
         procedure ::      Describe  => FarfieldBC_Describe
   end type FarFieldBC_t
!
!  ********************************
!  Outflow boundary condition class
!  ********************************
!
   type, extends(BoundaryCondition_t)           :: PressureOutletBC_t
      integer                            :: outflowType
      real(kind=RP), dimension(NCONS)    :: q
      contains
         procedure ::      Construct      => PressureOutletBC_Construct
         procedure ::      Associate      => PressureOutletBC_Associate
         procedure ::      UpdateSolution => PressureOutletBC_UpdateSolution
         procedure ::      Describe       => PressureOutletBC_Describe
   end type PressureOutletBC_t
!
!  ***************************************
!  Inflow/Outflow boundary condition class
!  ***************************************
!
   type, extends(BoundaryCondition_t)           :: PressureInletBC_t
      real(kind=RP), dimension(NCONS) :: q
      contains
         procedure ::      Construct      => PressureInletBC_Construct
         procedure ::      Associate      => PressureInletBC_Associate
         procedure ::      UpdateSolution => PressureInletBC_UpdateSolution
         procedure ::      Describe       => PressureInletBC_Describe
   end type PressureInletBC_t
!
!  ********************************
!  Riemann boundary condition class
!  ********************************
!
   type, extends(BoundaryCondition_t)           :: RiemannBC_t
      real(kind=RP), dimension(NCONS) :: q
      integer                         :: mode
      contains
         procedure ::      Construct      => RiemannBC_Construct
         procedure ::      Associate      => RiemannBC_Associate
         procedure ::      UpdateSolution => RiemannBC_UpdateSolution
         procedure ::      Describe       => RiemannBC_Describe
   end type RiemannBC_t
!
!  **************************************
!  New Dirichlet boundary condition class
!  **************************************
!
   type, extends(BoundaryCondition_t)           :: newDirichletBC_t
      integer                         :: outflowPressure
      integer                         :: mode
      real(kind=RP), dimension(NCONS) :: q
      contains
         procedure ::      Construct      => newDirichletBC_Construct
         procedure ::      Associate      => newDirichletBC_Associate
         procedure ::      UpdateSolution => newDirichletBC_UpdateSolution
         procedure ::      Describe       => newDirichletBC_Describe
   end type newDirichletBC_t
!
!  ***********************************
!  Euler wall boundary condition class
!  ***********************************
!
   type, extends(BoundaryCondition_t)           :: EulerWall_t
      contains
         procedure   ::    Associate      => EulerWall_Associate
         procedure   ::    UpdateSolution => EulerWall_UpdateSolution
         procedure   ::    Describe       => EulerWall_Describe
   end type EulerWall_t
!
!  *************************************
!  Viscous wall boundary condition class
!  *************************************
!
   type, extends(BoundaryCondition_t)           :: ViscousWall_t
      real(kind=RP)                 :: Tw
      integer                       :: wall_type
      real(kind=RP)                 :: v(NDIM)
      contains
         procedure   ::    Construct      => ViscousWall_Construct
         procedure   ::    Associate      => ViscousWall_Associate
         procedure   ::    UpdateSolution => ViscousWall_UpdateSolution
         procedure   ::    Describe       => ViscousWall_Describe
   end type ViscousWall_t
!
!  *******************
!  Construct procedure
!  *******************
!
   interface Construct
      module procedure BoundaryConditions_construct
   end interface Construct

   procedure(ICFcn), pointer     :: IC => NULL()
!
!///////////////////////////////////////////////////////////////////////////////////
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////
!
#include "./Dirichlet.incf"
#include "./EulerWall.incf"
#include "./Farfield.incf"
#include "./Periodic.incf"
#include "./PressureInlet.incf"
#include "./PressureOutlet.incf"
#include "./newDirichlet.incf"
#include "./Riemann.incf"
#include "./ViscousWall.incf"
#include "./BoundaryConditionFunctionDrivers.incf"

      subroutine BoundaryConditions_construct( self , marker)
         use Setup_class
         implicit none
         class(BoundaryCondition_t), pointer :: self
         integer                             :: marker
         character(len=STR_LEN_BC)           :: BCType
         character(len=STR_LEN_BC)           :: in_label
!
!        Get boundary type from case file
!        --------------------------------
         write(in_label,'(A,I0)') "# define zone ",marker
         call ReadValueInRegion( trim(Setup % bdry_file) , "Type" , BCType , in_label , "# end" )
!
!        ============================
         select case ( trim(BCType) ) 
!        ============================
!
!           ----------------------------------------
            case ( "Dirichlet" )
               allocate( DirichletBC_t    :: self )
               self % BCType = DIRICHLET_BC
               self % WeakType = WEAK_RIEMANN
!
!           ----------------------------------------   
            case ( "Farfield" )
               allocate( FarfieldBC_t     :: self )
               self % BCType = FARFIELD_BC
               self % WeakType = WEAK_RIEMANN
!
!           ----------------------------------------   
            case ( "Pressure outlet" )
               allocate( PressureOutletBC_t     :: self )
               self % BCType = PRESSUREOUTLET_BC
               self % WeakType = WEAK_RIEMANN
!
!           ----------------------------------------   
            case ( "Pressure inlet" )
               allocate( PressureInletBC_t     :: self )
               self % BCType = PRESSUREINLET_BC
               self % WeakType = WEAK_PRESCRIBED
!
!           ----------------------------------------   
            case ( "Periodic" )
               allocate( PeriodicBC_t     :: self )
               self % BCType = PERIODIC_BC
               self % WeakType = WEAK_RIEMANN
!
!           ----------------------------------------   
            case ( "Riemann" )
               allocate( RiemannBC_t     :: self )
               self % BCType = RIEMANN_BC
               self % WeakType = WEAK_PRESCRIBED

!           ----------------------------------------   
            case ("newDirichlet" )
               allocate( newDirichletBC_t    :: self )
               self % BCType = NEWDIRICHLET_BC
               self % WeakType = WEAK_RIEMANN

!           ----------------------------------------   
            case ( "Euler wall" )
               allocate( EulerWall_t      :: self )
               self % BCType = EULERWALL_BC
               self % WeakType = WEAK_PRESCRIBED
!
!           ----------------------------------------   
            case ( "Viscous wall" )
               allocate( ViscousWall_t    :: self ) 
               self % BCType = VISCOUSWALL_BC
               self % WeakType = WEAK_RIEMANN
!
!           ----------------------------------------   
            case default

               if ( len_trim (BCType) .eq. 0 ) then
                  print*, "Boundary condition type was not specified for marker " , marker , "."
               else
                  print*, 'Boundary condition "',trim(BCType),'" in marker ', marker, ' not implemented yet.'
                  print*, "Options available are:"
                  print*, "   * Dirichlet"
                  print*, "   * Farfield"
                  print*, "   * Pressure outlet"
                  print*, "   * Pressure inlet"
                  print*, "   * Riemann"
                  print*, "   * Periodic"
                  print*, "   * newDirichlet"
                  print*, "   * Euler wall"
                  print*, "   * Viscous wall"
                  stop "Stopped."

               end if
!         
!        **********
         end select
!        **********
!
!
!        Set marker
!        ----------   
         self % marker = marker
!
!        Obtain zone tag
!        ---------------
         call readValueInRegion( trim(Setup % bdry_file) , "Name" , self % Name , in_label , "# end" )
!
!        Set the Riemann solver
!        ----------------------
         call self % SetRiemannSolver()
!
!        Construct the zone
!        ------------------
         call self % Construct( marker , in_label)

      end subroutine BoundaryConditions_construct

      subroutine BoundaryConditions_SetRiemannSolver( self  )
!
!        ********************************************************************
!              This function sets the zone Riemann solver from the 
!           case file. 
!
!              #define zone X
!                 Riemann solver: Exact/Roe/HLL/Interior
!              #end
!        ********************************************************************
!
         use Setup_class
         implicit none
         class(BoundaryCondition_t)                :: self
         character(len=STR_LEN_BC)                 :: in_label
!
!        Get Riemann solver from case file
!        ---------------------------------         
         write(in_label,'(A,I0)') "# define zone ", self % marker
         call readValueInRegion( trim(Setup % bdry_file) , "Riemann solver" , self % RiemannSolverName ,  in_label , "# end" )
!
!        ************************************************
         select case ( trim( self % RiemannSolverName ) ) 
!        ************************************************
!         
!           ---------------------------------------------
            case ( "Exact" )
               self % RiemannSolver => ExactRiemannSolver
!         
!           ---------------------------------------------
            case ( "Roe" ) 
               self % RiemannSolver => RoeFlux
!         
!           ---------------------------------------------
            case ( "HLL" )
               self % RiemannSolver => HLLFlux
!         
!           ---------------------------------------------
            case ( "Interior" )
               self % RiemannSolver => NULL()
!         
!           ---------------------------------------------
            case default
               self % RiemannSolver => NULL()
               self % RiemannSolverName = "Interior"
!
!        **********
         end select
!        **********
!
      end subroutine BoundaryConditions_SetRiemannSolver

      subroutine BoundaryConditions_SetInitialCondition( IC_in )
         implicit none
         procedure(ICFcn), pointer      :: IC_in

         IC => IC_in

      end subroutine BoundaryConditions_SetInitialCondition
   
      subroutine BaseClass_Construct( self , marker , in_label)
         implicit none
         class(BoundaryCondition_t)                :: self
         integer                                   :: marker
         character(len=*)                          :: in_label
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Construct

      subroutine BaseClass_Associate( self , edge )
         implicit none
         class(BoundaryCondition_t)          :: self
         class(Edge_t)                       :: edge
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Associate

      subroutine BaseClass_UpdateSolution( self , edge , time)
         implicit none
         class(BoundaryCondition_t)          :: self
         class(Edge_t)                       :: edge
         real(kind=RP)                       :: time
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_UpdateSolution
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              Describe subroutines
!              --------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_Describe( self ) 
         implicit none
         class(BoundaryCondition_t)          :: self
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Describe

      subroutine DirichletBC_Describe( self )
         use Physics
         implicit none
         class(DirichletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )

         write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary condition type: " , "Dirichlet."
         write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         if ( .not. associated(self % BoundaryState) ) then
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , getPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , getTotalPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Temperature: "             , getTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total Temperature: "       , getTotalTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , sqrt( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP ) / ( getSoundSpeed(self % q) * self % q(IRHO) )
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Angle of attack: "         , atan2(self % q(IRHOV),self % q(IRHOU)) * 180.0_RP / PI
         else
            write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary state loaded: " , trim(self % BoundaryStateName)
         end if
         end associate

      end subroutine DirichletBC_Describe

      subroutine FarfieldBC_Describe( self )
         use Physics
         implicit none
         class(FarfieldBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )

         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Boundary condition type: " , "Farfield."
         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         if ( .not. associated(self % BoundaryState) ) then
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , getPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , getTotalPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Temperature: "             , getTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total Temperature: "       , getTotalTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , sqrt( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP ) / ( getSoundSpeed(self % q) * self % q(IRHO) )
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Angle of attack: "         , atan2(self % q(IRHOV),self % q(IRHOU)) * 180.0_RP / PI
         else
            write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary state loaded: " , trim(self % BoundaryStateName)
         end if

         end associate

      end subroutine FarfieldBC_Describe

      subroutine PressureOutletBC_Describe( self )
         use Physics
         implicit none
         class(PressureOutletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )

         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Boundary condition type: " , "Pressure outlet."
         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         if ( self % outflowType .eq. REFLECTIVE_OUTFLOW ) then
            write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Outflow type: "          , "Reflective"
         elseif ( self % outflowType .eq. PARTIALLY_NON_REFLECTIVE_OUTFLOW ) then
            write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Outflow type: "          , "Partially non-reflective"
         end if
 
         if ( .not. associated(self % BoundaryState) ) then
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , getPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , getTotalPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Temperature: "             , getTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total Temperature: "       , getTotalTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , sqrt( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP ) / ( getSoundSpeed(self % q) * self % q(IRHO) )
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Angle of attack: "         , atan2(self % q(IRHOV),self % q(IRHOU)) * 180.0_RP / PI
         else
            write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary state loaded: " , trim(self % BoundaryStateName)
         end if

         end associate

      end subroutine PressureOutletBC_Describe

      subroutine PressureInletBC_Describe( self )
         use Physics
         implicit none
         class(PressureInletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )

         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Boundary condition type: " , "Pressure inlet."
         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         if ( .not. associated(self % BoundaryState) ) then
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , getPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , getTotalPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Temperature: "             , getTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total Temperature: "       , getTotalTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , sqrt( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP ) / ( getSoundSpeed(self % q) * self % q(IRHO) )
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Angle of attack: "         , atan2(self % q(IRHOV),self % q(IRHOU)) * 180.0_RP / PI
         else
            write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary state loaded: " , trim(self % BoundaryStateName)
         end if

         end associate

      end subroutine PressureInletBC_Describe

      subroutine newDirichletBC_Describe( self )
         use Physics
         implicit none
         class(newDirichletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )

         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Boundary condition type: " , "new Dirichlet."
         if ( self % outflowPressure .eq. STATIC_PRESSURE ) then
            write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Outflow pressure: " , "Static."

         elseif ( self % outflowPressure .eq. TOTAL_PRESSURE ) then
            write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Outflow pressure: " , "Total."

         end if
   
         if ( self % mode .eq. SPECIFY_SPEED ) then
            write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Speed: " , "Specified."

         elseif ( self % mode .eq. SPECIFY_TOTAL_PRESSURE ) then
            write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Speed: " , "Computed from total pressure."

         end if


         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         if ( .not. associated(self % BoundaryState) ) then
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , getPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , getTotalPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Temperature: "             , getTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total Temperature: "       , getTotalTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , sqrt( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP ) / ( getSoundSpeed(self % q) * self % q(IRHO) )
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Angle of attack: "         , atan2(self % q(IRHOV),self % q(IRHOU)) * 180.0_RP / PI
         else
            write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary state loaded: " , trim(self % BoundaryStateName)
         end if


!         write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
!         write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , self % pt * refValues % p
!         write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
!         write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
!         write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
!         write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , self % Mach
!         write(STD_OUT , '(30X , A , A25 , I10   )') "-> " , "Angle of attack: "         , nint(self % AngleOfAttack * 180.0_RP / PI)
!
         end associate

      end subroutine newDirichletBC_Describe

      subroutine RiemannBC_Describe( self )
         use Physics
         implicit none
         class(RiemannBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )

         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Boundary condition type: " , "Pressure inlet."
         write(STD_OUT , '(30X , A , A25 , A     )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)

         if ( .not. associated(self % BoundaryState) ) then
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Pressure: "                , getPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total pressure: "          , getTotalPressure(self % q) * refValues % p
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Temperature: "             , getTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.2 )') "-> " , "Total Temperature: "       , getTotalTemperature(self % q) * refValues % T
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Density: "                 , self % q(IRHO) * refValues % rho
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "              , self % q(IRHOU) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "              , self % q(IRHOV) / self % q(IRHO) * refValues % a
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Mach number: "             , sqrt( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP ) / ( getSoundSpeed(self % q) * self % q(IRHO) )
            write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Angle of attack: "         , atan2(self % q(IRHOV),self % q(IRHOU)) * 180.0_RP / PI
         else
            write(STD_OUT , '(30X , A , A25 , A)') "-> " , "Boundary state loaded: " , trim(self % BoundaryStateName)
         end if

         end associate

      end subroutine RiemannBC_Describe

      subroutine PeriodicBC_Describe( self ) 
         implicit none
         class(PeriodicBC_t)          :: self

         write(STD_OUT , '(30X , A , A25 , A      )') "-> " , "Boundary condition type: " , "Periodic."
         write(STD_OUT , '(30X , A , A25 , A      )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X , A , A25 , I0 , A )') "-> " , "Connecting zone: "         , self % connected_marker
!
!        ********************************
         select case ( self % direction )
!        ********************************
!
!           --------------------------------------------------------------
            case (IX)
               write(STD_OUT , '(30X,A,A25,A)') "-> " , "Direction: " ,"x"
!
!           --------------------------------------------------------------
            case (IY)
               write(STD_OUT , '(30X,A,A25,A)') "-> " , "Direction: " ,"y"
!
!        **********
         end select
!        **********
!         
      end subroutine PeriodicBC_Describe

      subroutine EulerWall_Describe( self )
         implicit none
         class(EulerWall_t)                :: self

         write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Boundary condition type: " , "Euler wall."
!         if ( self % wall_type .eq. ADIABATIC_WALL ) then
!            write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Wall type: "          , "Adiabatic"
!         elseif ( self % wall_type .eq. ISOTHERMAL_WALL ) then
!            write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Wall type: "          , "Adiabatic"
!            write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Wall temperature: "   , self % Tw * refValues % T 
!         end if


         write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)

      end subroutine EulerWall_Describe

      subroutine ViscousWall_Describe( self )
         implicit none
         class(ViscousWall_t)                :: self

         write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Boundary condition type: " , "Viscous wall."
         if ( self % wall_type .eq. ADIABATIC_WALL ) then
            write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Wall type: "          , "Adiabatic"
         elseif ( self % wall_type .eq. ISOTHERMAL_WALL ) then
            write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Wall type: "          , "Adiabatic"
            write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Wall temperature: "   , self % Tw * refValues % T 
         end if
         write(STD_OUT , '( 30X , A , A25 , A )') "-> " , "Riemann solver: "          , trim(self % RiemannSolverName)
         write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "X-Velocity: "           , self % v(IX) * refValues % a
         write(STD_OUT , '(30X , A , A25 , F10.4 )') "-> " , "Y-Velocity: "           , self % v(IY) * refValues % a

      end subroutine ViscousWall_Describe

      function getTotalPressure(q) result (pt)
         implicit none
         real(kind=RP)     :: q(NCONS)
         real(kind=RP)     :: pt
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: p , Msq

         p = getPressure(q)
         Msq = ( q(IRHOU) * q(IRHOU) + q(IRHOV) * q(IRHOV) ) / ( thermodynamics % gamma * p * q(IRHO) )
         pt = p * (1.0_RP + 0.5_RP * thermodynamics % gm1 * Msq) ** (dimensionless % cp)
         
      end function getTotalPressure

      function getTotalTemperature(q) result (Tt)
         implicit none
         real(kind=RP)     :: q(NCONS)
         real(kind=RP)     :: Tt
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: T , Msq

         T = getTemperature(q)
         Msq = ( q(IRHOU) * q(IRHOU) + q(IRHOV) * q(IRHOV) ) / ( thermodynamics % gamma * getPressure(q) * q(IRHO) )
         Tt = T * (1.0_RP + 0.5_RP * thermodynamics % gm1 * Msq) 
         
      end function getTotalTemperature

end module DGBoundaryConditions
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
