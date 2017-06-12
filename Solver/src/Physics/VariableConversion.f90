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
!           Compute dimensionless variables
!
      module pure function getDimensionlessVariables(QWithDim) result ( Q )
         implicit none
         real(kind=RP), intent(in)           :: QWithDim(1:NCONS)
         real(kind=RP)                       :: Q(1:NCONS)

         Q(IRHO)  = QWithDim(IRHO)  / refValues % rho
         Q(IRHOU) = QWithDim(IRHOU) / ( refValues % rho * refValues % a )
         Q(IRHOV) = QWithDim(IRHOV) / ( refValues % rho * refValues % a )
         Q(IRHOE) = QWithDim(IRHOE) / refValues % p

      end function getDimensionlessVariables
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
         real(kind=RP), intent (in)  :: Q(1:NCONS,0:N)
         real(kind=RP)              :: p(0:N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i 

         do i = 0 , N 
            p(i) = thermodynamics % gm1 * (Q(IRHOE,i) - 0.5_RP * ( Q(IRHOU,i)*Q(IRHOU,i) + Q(IRHOV,i)*Q(IRHOV,i) ) / Q(IRHO,i) )
         end do

      end function getPressure1D

      module pure function getPressure2D(N,Q) result (p)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(1:NCONS,0:N,0:N)
         real(kind=RP)              :: p(0:N,0:N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i , j

         do j = 0 , N   ; do i = 0 , N 
            p(i,j) = thermodynamics % gm1 * (Q(IRHOE,i,j) - 0.5_RP * ( Q(IRHOU,i,j)*Q(IRHOU,i,j) + Q(IRHOV,i,j)*Q(IRHOV,i,j) ) / Q(IRHO,i,j) )
         end do         ; end do

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
         real(kind=RP), intent (in)  :: Q(1:NCONS,0:N)
         real(kind=RP)              :: T(0:N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i 

         T = getPressure1D(N,Q) 

         do i = 0 , N
#ifdef _DIMENSIONLESS_TAU
            T(i) = T(i) / Q(IRHO,i)
#else
            T(i) = dimensionless % gammaMach2 * T(i) / Q(IRHO,i)
#endif
         end do


      end function getTemperature1D

      module pure function getTemperature2D(N,Q) result (T)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(1:NCONS,0:N,0:N)
         real(kind=RP)              :: T(0:N,0:N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i , j

         T = getPressure2D(N,Q)

         do j = 0 , N   ; do i = 0 , N 
#ifdef _DIMENSIONLESS_TAU
            T(i,j) = T(i,j) / Q(IRHO,i,j) 
#else
            T(i,j) = dimensionless % gammaMach2 * T(i,j) / Q(IRHO,i,j)
#endif
         end do         ; end do

      end function getTemperature2D
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!           
!        Compute sound speed
!        -------------------
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
         real(kind=RP), intent(in)        :: Q(1:NCONS,0:N)
         real(kind=RP)                    :: a(0:N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i 

         a = getPressure1D(N,Q)

         do i = 0 , N 
            a(i) = sqrt(thermodynamics % gamma * a(i) / Q(IRHO,i) )
         end do

      end function getSoundSpeed1D

      module pure function getSoundSpeed2D(N,Q) result(a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(1:NCONS,0:N,0:N)
         real(kind=RP)                    :: a(0:N,0:N)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i , j 

         a = getPressure2D(N,Q) 

         do j = 0 , N   ; do i = 0 , N
            a(i,j) = sqrt( thermodynamics % gamma * a(i,j) / Q(IRHO,i,j) ) 
         end do         ; end do

      end function getSoundSpeed2D
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        Compute velocity gradients
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#ifdef NAVIER_STOKES
      module pure function getStrainTensor0D ( Q , dQ ) result ( S )
!
!        ***************************************************
!                 Computes the velocities gradients, which
!           are stored in the form:
!                 u(IX:IY , U:V)
!        ***************************************************
!
         implicit none
         real(kind=RP), intent(in)     :: Q(1:NCONS)
         real(kind=RP), intent(in)     :: dQ(1:NCONS,1:NDIM)
         real(kind=RP)                 :: S(1:NDIM , 1:NDIM)     
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

         S(IX,IX) = -uDivRho * dQ(IRHO,IX) + invRho * dQ(IRHOU,IX) 
         S(IY,IX) = -uDivRho * dQ(IRHO,IY) + invRho * dQ(IRHOU,IY) 
         S(IX,IY) = -vDivRho * dQ(IRHO,IX) + invRho * dQ(IRHOV,IX)
         S(IY,IY) = -vDivRho * dQ(IRHO,IY) + invRho * dQ(IRHOV,IY)

      end function getStrainTensor0D
 
      module pure function getStrainTensor1D ( N ,  Q , dQ ) result ( S )
!
!        ***************************************************
!                 Computes the velocities gradients, which
!           are stored in the form:
!                 u(IX:IY , U:V)
!        ***************************************************
!
         implicit none
         integer,       intent(in)         :: N
         real(kind=RP), intent(in)         :: Q(1:NCONS,0:N)
         real(kind=RP), intent(in)         :: dQ(1:NCONS,0:N,1:NDIM)
         real(kind=RP), target             :: S(0:N,1:NDIM , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i 
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v
         real(kind=RP), pointer :: ux(:) , uy(:) , vx(:) , vy(:)

         ux(0:) => S(0:,IX,IX)
         uy(0:) => S(0:,IY,IX)
         vx(0:) => S(0:,IX,IY)
         vy(0:) => S(0:,IY,IY)

         do i = 0 , N
            invRho = 1.0_RP / Q(IRHO,i)
      
            u = Q(IRHOU,i) * invRho
            v = Q(IRHOV,i) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho

            ux(i) = -uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX) 
            uy(i) = -uDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOU,i,IY) 
            vx(i) = -vDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOV,i,IX)
            vy(i) = -vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)
         end do

      end function getStrainTensor1D
 
      module pure function getStrainTensor2D ( N ,  Q , dQ ) result ( S )
!
!        ***************************************************
!                 Computes the velocities gradients, which
!           are stored in the form:
!                 u(IX:IY , U:V)
!        ***************************************************
!
         implicit none
         integer,       intent(in)         :: N
         real(kind=RP), intent(in)         :: Q(1:NCONS,0:N,0:N)
         real(kind=RP), intent(in)         :: dQ(1:NCONS,0:N,0:N,1:NDIM)
         real(kind=RP), target             :: S(0:N,0:N,1:NDIM , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i , j
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v
         real(kind=RP), pointer :: ux(:,:) , uy(:,:) , vx(:,:) , vy(:,:) 

         ux(0:,0:) => S(0:,0:,IX,IX)
         uy(0:,0:) => S(0:,0:,IY,IX)
         vx(0:,0:) => S(0:,0:,IX,IY)
         vy(0:,0:) => S(0:,0:,IY,IY)

         do j = 0 , N   ; do i = 0 , N
            invRho = 1.0_RP / Q(IRHO,i,j)
      
            u = Q(IRHOU,i,j) * invRho
            v = Q(IRHOV,i,j) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho

            ux(i,j) = -uDivRho * dq(IRHO,i,j,IX) + invRho * dq(IRHOU,i,j,IX) 
            uy(i,j) = -uDivRho * dq(IRHO,i,j,IY) + invRho * dq(IRHOU,i,j,IY) 
            vx(i,j) = -vDivRho * dq(IRHO,i,j,IX) + invRho * dq(IRHOV,i,j,IX)
            vy(i,j) = -vDivRho * dq(IRHO,i,j,IY) + invRho * dq(IRHOV,i,j,IY)
         end do      ; end do

      end function getStrainTensor2D

      module pure function getTemperatureGradient0D( Q , dQ ) result ( gradT )
         implicit none
         real(kind=RP),    intent(in)        ::  Q(1:NCONS)
         real(kind=RP),    intent(in)        :: dQ(1:NCONS,1:NDIM)
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

         gradT(IX) = gm1 * invRho * ( -Q(IRHOE) * invRho * dQ(IRHO,IX) + dQ(IRHOE,IX) + u * ( u * dQ(IRHO,IX) - dQ(IRHOU,IX)) + v * (v * dQ(IRHO,IX) - dQ(IRHOV,IX)) ) 
         gradT(IY) = gm1 * invRho * ( -Q(IRHOE) * invRho * dQ(IRHO,IY) + dQ(IRHOE,IY) + u * ( u * dQ(IRHO,IY) - dQ(IRHOU,IY)) + v * (v * dQ(IRHO,IY) - dQ(IRHOV,IY)) ) 
         
         end associate

      end function getTemperatureGradient0D

      module pure function getTemperatureGradient1D( N , Q , dQ ) result ( gradT )
         implicit none
         integer      ,    intent(in)         :: N
         real(kind=RP),    intent(in)         :: Q(1:NCONS,0:N)
         real(kind=RP),    intent(in)         :: dQ(1:NCONS,0:N,1:NDIM)
         real(kind=RP),    target             :: gradT(0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i
         real(kind=RP)          :: u , v , invRho
         real(kind=RP), pointer :: Tx(:) , Ty(:) 

         Tx(0:)   => gradT(0:,IX)
         Ty(0:)   => gradT(0:,IY)

         do i = 0 , N
            invRho = 1.0_RP / Q(IRHO,i)
            u  = q(IRHOU,i) * invRho
            v  = q(IRHOV,i) * invRho

            associate ( gm1 => thermodynamics % gm1 )

            Tx(i) = gm1 * invRho * ( -Q(IRHOE,i) * invRho * dq(IRHO,i,IX) + dq(IRHOE,i,IX) &
                                 + u * ( u * dq(IRHO,i,IX) - dq(IRHOU,i,IX)) + v * (v * dq(IRHO,i,IX) - dq(IRHOV,i,IX)) ) 
            Ty(i) = gm1 * invRho * ( -Q(IRHOE,i) * invRho * dq(IRHO,i,IY)+ dq(IRHOE,i,IY) &
                                 + u * ( u * dq(IRHO,i,IY) - dq(IRHOU,i,IY)) + v * (v * dq(IRHO,i,IY) - dq(IRHOV,i,IY)) ) 
            
            end associate
         end do

      end function getTemperatureGradient1D

      module pure function getTemperatureGradient2D( N , Q , dQ ) result ( gradT )
         implicit none
         integer      ,    intent(in)         :: N
         real(kind=RP),    intent(in)         :: Q(1:NCONS,0:N,0:N)
         real(kind=RP),    intent(in)         :: dQ(1:NCONS,0:N,0:N,1:NDIM)
         real(kind=RP),    target             :: gradT(0:N,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i , j 
         real(kind=RP)          :: u , v , invRho
         real(kind=RP), pointer :: Tx(:,:) , Ty(:,:) 

         Tx(0:,0:)   => gradT(0:,0:,IX)
         Ty(0:,0:)   => gradT(0:,0:,IY)

         do j = 0 , N   ; do i = 0 , N 
            invRho = 1.0_RP / Q(IRHO,i,j)
            u  = q(IRHOU,i,j) * invRho
            v  = q(IRHOV,i,j) * invRho

            associate ( gm1 => thermodynamics % gm1 )

            Tx(i,j) = gm1 * invRho * ( -Q(IRHOE,i,j) * invRho * dq(IRHO,i,j,IX) + dq(IRHOE,i,j,IX) &
                                 + u * ( u * dq(IRHO,i,j,IX) - dq(IRHOU,i,j,IX)) + v * (v * dq(IRHO,i,j,IX) - dq(IRHOV,i,j,IX)) ) 
            Ty(i,j) = gm1 * invRho * ( -Q(IRHOE,i,j) * invRho * dq(IRHO,i,j,IY)+ dq(IRHOE,i,j,IY) &
                                 + u * ( u * dq(IRHO,i,j,IY) - dq(IRHOU,i,j,IY)) + v * (v * dq(IRHO,i,j,IY) - dq(IRHOV,i,j,IY)) ) 
            
            end associate
         end do       ; end do

      end function getTemperatureGradient2D

#endif 

end submodule VariableConversion
