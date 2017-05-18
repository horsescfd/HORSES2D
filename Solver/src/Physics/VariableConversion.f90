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
submodule(PhysicsNS)   VariableConversion
   implicit none

#include "Defines.h"

   contains
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!           Compute pressure
!
!
      module pure function getPressure0D(Q) result (p)
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)                       :: p

         p = thermodynamics % gm1 * (Q(IRHOE) - 0.5_RP * ( Q(IRHOU)*Q(IRHOU) + Q(IRHOV)*Q(IRHOV) ) / Q(IRHO) )

      end function getPressure0D

      module pure function getPressure1D(N,Q) result (p)
         implicit none
         integer,       intent (in)  :: N
         real(kind=RP), intent (in)  :: Q(0:N,1:NCONS)
         real(kind=RP)              :: p(0:N)

         p = thermodynamics % gm1 * (Q(:,IRHOE) - 0.5_RP * ( Q(:,IRHOU)*Q(:,IRHOU) + Q(:,IRHOV)*Q(:,IRHOV) ) / Q(:,IRHO) )

      end function getPressure1D

      module pure function getPressure2D(N,Q) result (p)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)              :: p(0:N,0:N)

         p = thermodynamics % gm1 * (Q(:,:,IRHOE) - 0.5_RP * ( Q(:,:,IRHOU)*Q(:,:,IRHOU) + Q(:,:,IRHOV)*Q(:,:,IRHOV) ) / Q(:,:,IRHO) )

      end function getPressure2D
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!           
!        Compute Temperature
!
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      module pure function getTemperature0D(Q) result (T)
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)                      :: T

#ifdef _DIMENSIONLESS_TAU
         T = getPressure0D(Q) / Q(IRHO)
#else
         T = dimensionless % gammaMach2 * getPressure0D(Q) / Q(IRHO)
#endif

      end function getTemperature0D

      module pure function getTemperature1D(N,Q) result (T)
         implicit none
         integer,       intent (in)  :: N
         real(kind=RP), intent (in)  :: Q(0:N,1:NCONS)
         real(kind=RP)              :: T(0:N)

#ifdef _DIMENSIONLESS_TAU
         T = getPressure1D(N,Q) / Q(:,IRHO)
#else
         T = dimensionless % gammaMach2 * getPressure1D(N,Q) / Q(:,IRHO)
#endif

      end function getTemperature1D

      module pure function getTemperature2D(N,Q) result (T)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)              :: T(0:N,0:N)

#ifdef _DIMENSIONLESS_TAU
         T = getPressure2D(N,Q) / Q(:,:,IRHO) 
#else
         T = dimensionless % gammaMach2 * getPressure2D(N,Q) / Q(:,:,IRHO)
#endif

      end function getTemperature2D
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!           
!        Compute Temperature
!
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      module pure function getSoundSpeed0D(Q) result (a)
         implicit none
         real(kind=RP), intent(in)        :: Q(1:NCONS)
         real(kind=RP)                    :: a

         a = sqrt( thermodynamics % gamma * getPressure0D(Q) / Q(IRHO) )

      end function getSoundSpeed0D

      module pure function getSoundSpeed1D(N,Q) result (a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(0:N,1:NCONS)
         real(kind=RP)                    :: a(0:N)

         a = sqrt( thermodynamics % gamma * getPressure1D(N,Q) / Q(:,IRHO) )

      end function getSoundSpeed1D

      module pure function getSoundSpeed2D(N,Q) result(a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)                    :: a(0:N,0:N)

         a = sqrt( thermodynamics % gamma * getPressure2D(N,Q) / Q(:,:,IRHO) )

      end function getSoundSpeed2D
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        Compute velocity gradients
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#ifdef NAVIER_STOKES
      module pure function getStrainTensor0D ( Q , dQ ) result ( du )
!
!        ***************************************************
!                 Computes the velocities gradients, which
!           are stored in the form:
!                 u(IX:IY , U:V)
!        ***************************************************
!
         implicit none
         real(kind=RP), intent(in)     :: Q(1:NCONS)
         real(kind=RP), intent(in)     :: dQ(1:NDIM,1:NCONS)
         real(kind=RP)                 :: du(1:NDIM , 1:NDIM)     
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v 

         invRho = 1.0_RP / Q(IRHO)
      
         u = Q(IRHOU) * invRho
         v = Q(IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         du(IX,IX) = -uDivRho * dQ(IX,IRHO) + invRho * dQ(IX,IRHOU) 
         du(IY,IX) = -uDivRho * dQ(IY,IRHO) + invRho * dQ(IY,IRHOU) 
         du(IX,IY) = -vDivRho * dQ(IX,IRHO) + invRho * dQ(IX,IRHOV)
         du(IY,IY) = -vDivRho * dQ(IY,IRHO) + invRho * dQ(IY,IRHOV)

      end function getStrainTensor0D
 
      module pure function getStrainTensor1D ( N ,  Q , dQ ) result ( du )
!
!        ***************************************************
!                 Computes the velocities gradients, which
!           are stored in the form:
!                 u(IX:IY , U:V)
!        ***************************************************
!
         implicit none
         integer,       intent(in)     :: N
         real(kind=RP), intent(in)     :: Q(0:N,1:NCONS)
         real(kind=RP), intent(in)     :: dQ(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                 :: du(0:N,1:NDIM , 1:NDIM)     
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N) , uDivRho(0:N) , vDivRho(0:N) , u(0:N) , v(0:N) 


         invRho = 1.0_RP / Q(:,IRHO)
      
         u = Q(:,IRHOU) * invRho
         v = Q(:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         du(:,IX,IX) = -uDivRho * dQ(:,IX,IRHO) + invRho * dQ(:,IX,IRHOU) 
         du(:,IY,IX) = -uDivRho * dQ(:,IY,IRHO) + invRho * dQ(:,IY,IRHOU) 
         du(:,IX,IY) = -vDivRho * dQ(:,IX,IRHO) + invRho * dQ(:,IX,IRHOV)
         du(:,IY,IY) = -vDivRho * dQ(:,IY,IRHO) + invRho * dQ(:,IY,IRHOV)

      end function getStrainTensor1D
 
      module pure function getStrainTensor2D ( N ,  Q , dQ ) result ( du )
!
!        ***************************************************
!                 Computes the velocities gradients, which
!           are stored in the form:
!                 u(IX:IY , U:V)
!        ***************************************************
!
         implicit none
         integer,       intent(in)     :: N
         real(kind=RP), intent(in)     :: Q(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)     :: dQ(0:N,0:N,1:NDIM,1:NCONS)
         real(kind=RP)                 :: du(0:N,0:N,1:NDIM , 1:NDIM)     
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N,0:N) , uDivRho(0:N,0:N) , vDivRho(0:N,0:N) , u(0:N,0:N) , v(0:N,0:N) 


         invRho = 1.0_RP / Q(:,:,IRHO)
      
         u = Q(:,:,IRHOU) * invRho
         v = Q(:,:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         du(:,:,IX,IX) = -uDivRho * dQ(:,:,IX,IRHO) + invRho * dQ(:,:,IX,IRHOU) 
         du(:,:,IY,IX) = -uDivRho * dQ(:,:,IY,IRHO) + invRho * dQ(:,:,IY,IRHOU) 
         du(:,:,IX,IY) = -vDivRho * dQ(:,:,IX,IRHO) + invRho * dQ(:,:,IX,IRHOV)
         du(:,:,IY,IY) = -vDivRho * dQ(:,:,IY,IRHO) + invRho * dQ(:,:,IY,IRHOV)

      end function getStrainTensor2D

      module pure function getTemperatureGradient0D( Q , dQ ) result ( gradT )
         implicit none
         real(kind=RP),    intent(in)        ::  Q(1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(1:NDIM , 1:NCONS)
         real(kind=RP)                       :: gradT(1:NDIM) 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: u , v , invRho 

         invRho = 1.0_RP / Q(IRHO)
         u  = q(IRHOU) * invRho
         v  = q(IRHOV) * invRho

         associate ( gm1 => thermodynamics % gm1 )

         gradT(IX) = gm1 * invRho * ( -Q(IRHOE) * invRho * dQ(IX,IRHO) + dQ(IX,IRHOE) + u * ( u * dQ(IX,IRHO) - dQ(IX,IRHOU)) + v * (v * dQ(IX,IRHO) - dQ(IX,IRHOV)) ) 
         gradT(IY) = gm1 * invRho * ( -Q(IRHOE) * invRho * dQ(IY,IRHO) + dQ(IY,IRHOE) + u * ( u * dQ(IY,IRHO) - dQ(IY,IRHOU)) + v * (v * dQ(IY,IRHO) - dQ(IY,IRHOV)) ) 
         
         end associate

      end function getTemperatureGradient0D

      module pure function getTemperatureGradient1D( N , Q , dQ ) result ( gradT )
         implicit none
         integer      ,    intent(in)        :: N 
         real(kind=RP),    intent(in)        ::  Q(0:N,1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(0:N,1:NDIM , 1:NCONS)
         real(kind=RP)                       :: gradT(0:N,1:NDIM) 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: u(0:N) , v(0:N) , invRho(0:N) 

         invRho = 1.0_RP / Q(:,IRHO)
         u  = q(:,IRHOU) * invRho
         v  = q(:,IRHOV) * invRho

         associate ( gm1 => thermodynamics % gm1 )

         gradT(:,IX) = gm1 * invRho * ( -Q(:,IRHOE) * invRho * dQ(:,IX,IRHO) + dQ(:,IX,IRHOE) &
                              + u * ( u * dQ(:,IX,IRHO) - dQ(:,IX,IRHOU)) + v * (v * dQ(:,IX,IRHO) - dQ(:,IX,IRHOV)) ) 
         gradT(:,IY) = gm1 * invRho * ( -Q(:,IRHOE) * invRho * dQ(:,IY,IRHO) + dQ(:,IY,IRHOE) &
                              + u * ( u * dQ(:,IY,IRHO) - dQ(:,IY,IRHOU)) + v * (v * dQ(:,IY,IRHO) - dQ(:,IY,IRHOV)) ) 
         
         end associate

      end function getTemperatureGradient1D

      module pure function getTemperatureGradient2D( N , Q , dQ ) result ( gradT )
         implicit none
         integer      ,    intent(in)        :: N 
         real(kind=RP),    intent(in)        ::  Q(0:N,0:N,1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(0:N,0:N,1:NDIM , 1:NCONS)
         real(kind=RP)                       :: gradT(0:N,0:N,1:NDIM) 
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: u(0:N,0:N) , v(0:N,0:N) , invRho(0:N,0:N) 

         invRho = 1.0_RP / Q(:,:,IRHO)
         u  = q(:,:,IRHOU) * invRho
         v  = q(:,:,IRHOV) * invRho

         associate ( gm1 => thermodynamics % gm1 )

         gradT(:,:,IX) = gm1 * invRho * ( -Q(:,:,IRHOE) * invRho * dQ(:,:,IX,IRHO) + dQ(:,:,IX,IRHOE) &
                              + u * ( u * dQ(:,:,IX,IRHO) - dQ(:,:,IX,IRHOU)) + v * (v * dQ(:,:,IX,IRHO) - dQ(:,:,IX,IRHOV)) ) 
         gradT(:,:,IY) = gm1 * invRho * ( -Q(:,:,IRHOE) * invRho * dQ(:,:,IY,IRHO) + dQ(:,:,IY,IRHOE) &
                              + u * ( u * dQ(:,:,IY,IRHO) - dQ(:,:,IY,IRHOU)) + v * (v * dQ(:,:,IY,IRHO) - dQ(:,:,IY,IRHOV)) ) 
         
         end associate

      end function getTemperatureGradient2D

#endif 

end submodule VariableConversion
