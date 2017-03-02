!
! ///////////////////////////////////////////////////////////////////////////////////////
!
!     SMConstants.F
!
!!
!!     Modification History:
!!       version 0.0 August 10, 2005 David A. Kopriva
!
!     MODULE SMConstants
!
!!        Defines constants for use by the spectral demonstaration
!!        routines, including precision definitions. 
!
!!    @author David A. Kopriva
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      module SMConstants
!
!     *************************************************************************
!           Floating point parameters 
!     *************************************************************************
!
           integer       , parameter, private :: DOUBLE_DIGITS   =  15                                  ! # of desired digits
           integer       , parameter, private :: SINGLE_DIGITS   =  6                                   ! # of desired digits
           integer       , parameter          :: RP              =  SELECTED_REAL_KIND( DOUBLE_DIGITS ) ! Real Kind
           integer       , parameter          :: SP              =  SELECTED_REAL_KIND( SINGLE_DIGITS ) ! Single Real Kind
           integer       , parameter          :: CP              =  SELECTED_REAL_KIND( DOUBLE_DIGITS ) ! Complex Kind
!
!     *************************************************************************
!           Constants
!     *************************************************************************
!
           integer       , parameter          :: FORWARD         = 1
           integer       , parameter          :: BACKWARD       = -1
           integer       , parameter          :: LEFT            = 2
           integer       , parameter          :: RIGHT           = 1
           real(kind=RP) , parameter          :: PI              =  3.141592653589793238462643_RP
           complex(kind=CP)                   :: ImgI            =  ( 0.0_RP, 1.0_RP)                   !                         =  SQRT(-1.0_RP)
!
!     *************************************************************************
!           Interpolation node type aliases              
!     *************************************************************************
!
           integer, parameter                 :: LG              =  1               ! Parameter for Legendre-Gauss nodes
           integer, parameter                 :: LGL             =  2               ! Parameter for Legendre-Gauss-Lobatto nodes
!
!     *************************************************************************
!           Equation type aliases 
!     *************************************************************************
!
           integer, parameter                 :: FORMI           =  1               ! Green form
           integer, parameter                 :: FORMII          =  2               ! Divergence form
!
!
!     *************************************************************************
!           Parameters for I/O
!     *************************************************************************
!
           integer, parameter                 :: STD_OUT         =  6
           integer, parameter                 :: STD_IN          =  5
           integer, parameter                 :: LINE_LENGTH     =  132
!
!     *************************************************************************
!           Boundary conditions and faces classification
!     *************************************************************************
!
           integer, parameter                 :: FACE_INTERIOR     = 0
           integer, parameter                 :: PERIODIC_BC       = 1
           integer, parameter                 :: DIRICHLET_BC      = 2
           integer, parameter                 :: EULERWALL_BC      = 3
           integer, parameter                 :: VISCOUSWALL_BC    = 4
           integer, parameter                 :: FARFIELD_BC       = 5
           integer, parameter                 :: PRESSUREOUTLET_BC = 6
           integer, parameter                 :: PRESSUREINLET_BC  = 7
           integer, parameter                 :: RIEMANN_BC        = 8
!
!     *************************************************************************
!           Time integration mode
!     *************************************************************************
!             
           integer, parameter                 :: STEADY          =  0
           integer, parameter                 :: TRANSIENT       =  1
     
      end module SMConstants
