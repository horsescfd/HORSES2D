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
    
      end module SMConstants
