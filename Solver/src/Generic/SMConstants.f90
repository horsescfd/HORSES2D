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
           integer, parameter                 :: FACE_INTERIOR   =  0
           integer, parameter                 :: FACE_BOUNDARY   =  1
     
           integer, parameter                 :: PERIODIC_BC     =  0
           integer, parameter                 :: DIRICHLET_BC    =  1
           integer, parameter                 :: EULERWALL_BC    =  2
           integer, parameter                 :: VISCOUSWALL_BC  =  3
           integer, parameter                 :: INFLOW_BC       =  4
           integer, parameter                 :: OUTFLOW_BC      =  5
!
!     *************************************************************************
!           Time integration mode
!     *************************************************************************
!             
           integer, parameter                 :: STEADY          =  0
           integer, parameter                 :: TRANSIENT       =  1
     
      end module SMConstants
