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
submodule (PhysicsNS)  InviscidFluxes
   use SMConstants

#include "Defines.h"

   contains

      module pure function inviscidFlux0D(u) result(val)
         implicit none
         real(kind=RP), intent(in) :: u(NCONS)
         real(kind=RP)             :: val(NCONS,NDIM)
         real(kind=RP)             :: vx , vy  , p

         associate ( Gamma => Thermodynamics % Gamma , gm1 => Thermodynamics % gm1 , Mach => Dimensionless % Mach ) 

         vx = u(IRHOU) / u(IRHO)
         vy = u(IRHOV) / u(IRHO)
         p  = gm1 * ( u(IRHOE) - 0.5_RP * u(IRHOU) * vx - 0.5_RP * u(IRHOV) * vy )
         
         val(IRHO,IX)  = u(IRHOU)
         val(IRHOU,IX) = u(IRHOU) * vx + p 
         val(IRHOV,IX) = u(IRHOU) * vy
         val(IRHOE,IX) = (u(IRHOE) + p) * vx

         val(IRHO,IY)  = u(IRHOV)
         val(IRHOU,IY) = val(IRHOV,IX)
         val(IRHOV,IY) = u(IRHOV) * vy + p
         val(IRHOE,IY) = (u(IRHOE) + p) * vy

#ifdef _DIMENSIONLESS_TAU
         val = val * dimensionless % invSqrtGammaMach
#endif

         end associate

      end function inviscidFlux0D

      module pure function inviscidFlux1D(N,u) result(val)
         implicit none
         integer, intent(in)       :: N
         real(kind=RP), intent(in) :: u(1:NCONS,0:N)
         real(kind=RP), target     :: val(1:NCONS,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                   :: i 
         real(kind=RP)             :: vx , vy , p , invRho
         real(kind=RP), pointer    :: f(:,:) , g(:,:)

         associate ( Gamma => Thermodynamics % Gamma , gm1 => Thermodynamics % gm1 , Mach => Dimensionless % Mach ) 
    
         f(1:,0:) => val(1:,0:,IX)
         g(1:,0:) => val(1:,0:,IY)

         do i = 0 , N
            invRho = 1.0_RP / u(IRHO,i)
            vx = u(IRHOU,i) * invRho
            vy = u(IRHOV,i) * invRho
            p  = gm1 * ( u(IRHOE,i) - 0.5_RP * u(IRHOU,i) * vx - 0.5_RP * u(IRHOV,i) * vy )
         
            f(IRHO ,i) = u(IRHOU,i)
            f(IRHOU,i) = u(IRHOU,i) * vx + p
            f(IRHOV,i) = u(IRHOU,i) * vy
            f(IRHOE,i) = ( u(IRHOE,i) + p ) * vx

            g(IRHO ,i) = u(IRHOV,i)  
            g(IRHOU,i) = f(IRHOV,i)
            g(IRHOV,i) = u(IRHOV,i) * vy + p
            g(IRHOE,i) = (u(IRHOE,i) + p) * vy 

         end do

#ifdef _DIMENSIONLESS_TAU
         val = val * dimensionless % invSqrtGammaMach
#endif

         end associate

      end function inviscidFlux1D

      module pure function InviscidFlux2D(N,u) result(val)
         implicit none
         integer, intent(in)       :: N
         real(kind=RP), intent(in) :: u(1:NCONS,0:N,0:N)
         real(kind=RP), target     :: val(1:NCONS,0:N,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                   :: i , j
         real(kind=RP)             :: vx , vy , p , invRho
         real(kind=RP), pointer    :: f(:,:,:) , g(:,:,:)

         associate ( Gamma => Thermodynamics % Gamma , gm1 => Thermodynamics % gm1 , Mach => Dimensionless % Mach ) 
    
         f(1:,0:,0:) => val(1:,0:,0:,IX)
         g(1:,0:,0:) => val(1:,0:,0:,IY)

         do j = 0 , N ;    do i = 0 , N
            invRho = 1.0_RP / u(IRHO,i,j)
            vx = u(IRHOU,i,j) * invRho
            vy = u(IRHOV,i,j) * invRho
            p  = gm1 * ( u(IRHOE,i,j) - 0.5_RP * u(IRHOU,i,j) * vx - 0.5_RP * u(IRHOV,i,j) * vy )
         
            f(IRHO ,i,j) = u(IRHOU,i,j)
            f(IRHOU,i,j) = u(IRHOU,i,j) * vx + p
            f(IRHOV,i,j) = u(IRHOU,i,j) * vy
            f(IRHOE,i,j) = ( u(IRHOE,i,j) + p ) * vx

            g(IRHO ,i,j) = u(IRHOV,i,j)  
            g(IRHOU,i,j) = f(IRHOV,i,j)
            g(IRHOV,i,j) = u(IRHOV,i,j) * vy + p
            g(IRHOE,i,j) = (u(IRHOE,i,j) + p) * vy 

         end do ;          end do

#ifdef _DIMENSIONLESS_TAU
         val = val * dimensionless % invSqrtGammaMach
#endif

         end associate

      end function inviscidFlux2D

      module pure function F_inviscidFlux(rho,u,v,p,H) result(F)
         implicit none
         real(kind=RP), intent(in) :: rho
         real(kind=RP), intent(in) :: u
         real(kind=RP), intent(in) :: v
         real(kind=RP), intent(in) :: p
         real(kind=RP), intent(in) :: H
         real(kind=RP)             :: F(NCONS)
   
         F(IRHO)  = rho * u
         F(IRHOU) = F(IRHO) * u + p
         F(IRHOV) = F(IRHO) * v
         F(IRHOE) = H * F(IRHO)

      end function F_inviscidFlux

end submodule InviscidFluxes
