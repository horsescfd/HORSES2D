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
! ///////////////////////////////////////////////////////////////////////////////////////
!
!     SMConstants.F
!
!
!     Modification History:
!       version 0.0 August 10, 2005 David A. Kopriva
!
!     MODULE SMConstants
!
!        Defines constants for use by the spectral demonstaration
!        routines, including precision definitions. 
!
!    @author David A. Kopriva
!
! ////////////////////////////////////////////////////////////////////////////////////////
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
