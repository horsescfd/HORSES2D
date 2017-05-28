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
module InitialConditions
   use SMConstants
   use Physics
   implicit none

#include "Defines.h"
!
!  *******
   private
   public ICFcn , getInitialCondition , InitialCondition_Describe , SetLCheck
!  *******
!

   integer,            parameter :: STR_LEN_IC     = 128
   character(len = *), parameter :: ConstantIC     = "Uniform"
   character(len = *), parameter :: SteadyIC       = "Steady"
   character(len = *), parameter :: TurbulenceIC   = "Turbulence2D"
   character(len = *), parameter :: UserDefinedIC  = "User-defined"
   character(len = *), parameter :: RestartIC      = "Restart"
   character(len = *), parameter :: ChecksPolyIC   = "ChecksPolynomic"
   character(len = *), parameter :: ChecksTrigIC   = "ChecksTrigonometric"
   real(kind=RP)                 :: Lcheck

! 
!  ******************
   abstract interface
!  ******************
!
      function ICFcn (x) result (val)
         use SMConstants
         use Physics
         implicit none
         real(kind=RP)           :: x(NDIM)
         real(kind=RP)           :: val(NCONS)
      end function ICFcn
!
!  *************
   end interface
!  *************
!
!
!  ========   
   contains
!  ========   
!

      subroutine getInitialCondition( fcn , which )
         use SMConstants
         use Setup_class
         implicit none
         procedure(ICFcn), pointer     :: fcn
         character(len=*), optional    :: which
         character(len=STR_LEN_IC)     :: ICName
         
         if (present(which)) then
            ICName = which
         else
            ICName = Setup % IC
         end if

!
!        ********************************
         select case ( trim(ICName) )
!        ********************************
!
!           =========================
            case ( trim(ConstantIC) ) 
!           =========================
!
               fcn => UniformInitialCondition
!
!           ============================
            case ( trim(SteadyIC) )
!           ============================
!
               fcn => SteadyInitialCondition
!
!           ============================
            case ( trim(UserDefinedIC) )
!           ============================
!
               fcn => UserDefinedInitialConditionDriver
!
!           ============================
            case ( trim(ChecksPolyIC) )
!           ============================
!
               fcn => ChecksPolynomicInitialCondition
!
!           ============================
            case ( trim(ChecksTrigIC) )
!           ============================
!
               fcn => ChecksTrigonometricInitialCondition

!
!           ============================
            case ( trim(TurbulenceIC) )
!           ============================
!
               fcn => Turbulence2DInitialCondition
!
!           ============================
            case ( trim(RestartIC) )
!           ============================
!
               fcn => NULL()
!
!
!
!           ============               
            case default
!           ============               
!
               print*, 'Error: Unknown initial condition "',trim(ICName),'".'
               print*, "   Options available:"
               print*, "      * ",trim(ConstantIC)
               print*, "      * ",trim(SteadyIC)
               print*, "      * ",trim(ChecksPolyIC)
               print*, "      * ",trim(ChecksTrigIC)
               print*, "      * ",trim(RestartIC) 
               print*, "      * ",trim(UserDefinedIC) 
               STOP "Stopped." 
!
!        **********
         end select
!        **********
!


      end subroutine getInitialCondition
            
      function UniformInitialCondition(x) result(val)
!        ***************************************************************
!           Loads an uniform initial condition from reference values
!        ***************************************************************
         use SMConstants
         use Setup_Class
         implicit none
         real(kind=RP)            :: x(NDIM)
         real(kind=RP)            :: val(NCONS)
         real(kind=RP), parameter :: AngleOfAttack = 0.0_RP

         val(IRHO)  = refValues % rho
         val(IRHOU) = refValues % rho * refValues % V * cos(AngleOfAttack)
         val(IRHOV) = refValues % rho * refValues % V * sin(AngleOfAttack)
         val(IRHOE) = thermodynamics % cv * refValues % rho * refValues % T + 0.5_RP * refValues % rho * refValues % V * refValues % V

      end function UniformInitialCondition

      function UserDefinedInitialConditionDriver( x ) result ( val )
         use SMConstants
         use Setup_Class
         use Physics
         implicit none
         real(kind=RP)        :: x(NDIM)
         real(kind=RP)           :: val(NCONS)
         interface
            function UserDefinedInitialCondition(x , Thermodynamics_ , Setup_ , refValues_ , dimensionless_ ) result (val)
               use SMConstants
               use Setup_class
               use Physics
               implicit none
               real(kind=RP)                        :: x(NDIM)
               class(Thermodynamics_t), intent(in)  :: thermodynamics_
               class(Setup_t),          intent(in)  :: Setup_
               class(RefValues_t),      intent(in)  :: refValues_
               class(Dimensionless_t),  intent(in)  :: dimensionless_
               real(kind=RP)                        :: val(NCONS)
            end function UserDefinedInitialCOndition
         end interface 

         val = UserDefinedInitialCondition( x , Thermodynamics , Setup , refValues , dimensionless )

      end function UserDefinedInitialConditionDriver

      function SteadyInitialCondition(x) result(val)
!        ***************************************************************
!           Loads a steady flow from the reference values
!        ***************************************************************
         use SMConstants
         use Setup_Class
         implicit none
         real(kind=RP)           :: x(NDIM)
         real(kind=RP)           :: val(NCONS)
         real(kind=RP)           :: pressure

         val(IRHO)  = refValues % rho
         val(IRHOU) = 0.0_RP
         val(IRHOV) = 0.0_RP
         val(IRHOE) = thermodynamics % cv * refValues % rho * refValues % T

      end function SteadyInitialCondition
         
      function TaylorVortexInitialCondition(x) result(val)
!        ****************************************************************
!           Loads an initial condition consistent with the Taylor Green 
!          vortex problem
!        ****************************************************************
         use SMConstants
         implicit none
         real(kind=RP)           :: x(NDIM)
         real(kind=RP)           :: val(NCONS)
         real(kind=RP)           :: rho , u , v , p

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         rho = 1.0_RP
         u = sqrt(gamma) * Mach * sin(PI*x(1)) * cos(PI*x(2))
         v = -sqrt(gamma) * Mach *  cos(PI*x(1)) * sin(PI*x(2))
         p = 1.0_RP + 0.125_RP * gamma * Mach * Mach * (cos(2.0_RP * PI * x(1)) + cos(2.0_RP * PI * x(2)) )

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
            
      end function TaylorVortexInitialCondition

      function Turbulence2DInitialCondition(x) result(val)
!        ****************************************************************
!           Loads an initial condition to generate 2D Turbulence
!        ****************************************************************
         use SMConstants
         implicit none
         real(kind=RP)            :: x(NDIM)
         real(kind=RP)            :: val(NCONS)
         real(kind=RP)            :: rho , u , v , p
         real(kind=RP), parameter :: k = 10.0_RP

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         rho = 1.0_RP
         u = -sqrt(gamma) * Mach * cos(2.0_RP * k * PI*x(1)) * sin(2.0_RP * k * PI*x(2))
         v = sqrt(gamma) * Mach *  sin(2.0_RP * k * PI*x(1)) * cos(2.0_RP * k * PI*x(2))
         p = 1.0_RP 

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
            
      end function Turbulence2DInitialCondition

      function ChecksPolynomicInitialCondition (x) result (val)
         use SMConstants
         use Physics
         implicit none
         real(kind=RP)              :: x(NDIM)
         real(kind=RP)              :: val(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)              :: rho , u , v , p

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => dimensionless % cv )

         rho = refValues % rho
         u = refValues % V * x(IX) / LCheck
         v = refValues % V * x(IY) / LCheck
         p = ( 1.0_RP + gamma * Mach * Mach * ( x(IX) * x(IX)  + x(IY) * x(IY) ) / (LCheck **2.0_RP) ) * refValues % p

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
         
      end function ChecksPolynomicInitialCondition

      function ChecksTrigonometricInitialCondition (x) result (val)
         use SMConstants
         use Physics
         implicit none
         real(kind=RP)              :: x(NDIM)
         real(kind=RP)              :: val(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: rho , u , v , p

         associate ( gamma => Thermodynamics % Gamma , Mach => Dimensionless % Mach , cv => Dimensionless % cv )

         rho = refValues % rho
         u = refValues % V *  sin(PI * x(IX) / LCheck) * cos(PI*x(IY) / LCheck)
         v = -refValues % V * cos(PI*x(IX)/LCheck) * sin(PI*x(IY) / LCheck)
         p = refValues % p * ( 1.0_RP + 0.125_RP * gamma * Mach * Mach * ( cos(2.0_RP * PI * x(IX) / LCheck) + cos(2.0_RP * PI * x(IY) / LCheck ) ) )

         val(IRHO)  = rho
         val(IRHOU) = rho * u
         val(IRHOV) = rho * v
         val(IRHOE) = cv * p + 0.5_RP * rho * ( u*u + v*v )

         end associate
         
      end function ChecksTrigonometricInitialCondition

      subroutine SetLCheck( val ) 
         implicit none
         real(kind=RP),    intent(in)     :: val
         
         LCheck = val

      end subroutine SetLCheck
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        DESCRIBE
!        --------
!/////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine InitialCondition_Describe()
         use Headers
         use Setup_class
         implicit none
   
         write(STD_OUT ,'(/)')
         call Section_header("Initial condition description")
         write(STD_OUT ,*)
         
         write(STD_OUT , '(30X,A,A25,A,A)') "-> ", "Initial condition type: " , trim(Setup % IC) , "."
   
         if ( trim(Setup % IC) .eq. "Restart" ) then
            write(STD_OUT , '(30X,A,A25,A,A)') "-> ", "Data restart from file: " , trim(Setup % restart_file) , "."
         end if

         write(STD_OUT , '(30X,A,A25,ES10.4,A)' ) "-> " , "Initial time: " , Setup % initialTime

      end subroutine InitialCondition_Describe

end module InitialConditions
