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
!///////////////////////////////////////////////////////////////////////////////////////////////
!
!     File: ViscousFluxes.f90
!     Description:   This is a submodule of PhysicsNS.f90 containing the viscous fluxes 
!                  procedures.
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
submodule (PhysicsNS)   ViscousFluxes
   use SMConstants
   implicit none

#include "Defines.h"

   contains

      module pure function viscousFlux0D( q , dq) result(F)
         implicit none
         real(kind=RP), intent(in)  :: q(NCONS)
         real(kind=RP), intent(in)  :: dq(NCONS,NDIM)
         real(kind=RP)              :: F(1:NCONS,NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV , dxT , dyT

         invRho  = 1.0_RP / q(IRHO)
         u       = q(IRHOU) * invRho
         v       = q(IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(IRHO,IX) + invRho * dq(IRHOU,IX) - vDivRho * dq(IRHO,IY) + invRho * dq(IRHOV,IY)

         dxT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IRHO,IX) + dq(IRHOE,IX) + u * ( u * dq(IRHO,IX) - dq(IRHOU,IX)) + v * (v * dq(IRHO,IX) - dq(IRHOV,IX)) ) 
         dyT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IRHO,IY) + dq(IRHOE,IY) + u * ( u * dq(IRHO,IY) - dq(IRHOU,IY)) + v * (v * dq(IRHO,IY) - dq(IRHOV,IY)) ) 
   
         F(IRHO,IX) = 0.0_RP
         F(IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(IRHO,IX) + invRho * dq(IRHOU,IX)) + lambda * divV
         F(IRHOV,IX) = - vDivRho * dq(IRHO,IX)  &
                       - uDivRho * dq(IRHO,IY)  &
                       + invRho  * dq(IRHOU,IY) &
                       + invRho  * dq(IRHOV,IX) 
         F(IRHOE,IX) = u * F(IRHOU,IX) + v * F(IRHOV,IX) + ( dimensionless % cp / dimensionless % Pr ) * dxT

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX)
         F(IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(IRHO,IY) + invRho * dq(IRHOV,IY)) + lambda * divV
         F(IRHOE,IY) = u * F(IRHOU,IY) + v * F(IRHOV,IY) + ( dimensionless % cp / dimensionless % Pr ) * dyT

         end associate

      end function viscousFlux0D

      module pure function viscousFlux1D( N , q , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(1:NCONS,0:N)
         real(kind=RP), intent(in)          :: dq(1:NCONS,0:N,1:NDIM)
         real(kind=RP), target              :: F(1:NCONS,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i 
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV , dxT , dyT
         real(kind=RP), pointer  :: fv(:,:) , gv(:,:)

         fv(1:,0:)    => f(1:,0:,IX)
         gv(1:,0:)    => f(1:,0:,IY)

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 

         do i = 0 , N
            invRho  = 1.0_RP / q(IRHO,i)
            u       = q(IRHOU,i) * invRho
            v       = q(IRHOV,i) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho
   
            divV = - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX) - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)
   
            dxT = gm1 * invRho * ( -q(IRHOE,i) * invRho * dq(IRHO,i,IX) + dq(IRHOE,i,IX) + u * ( u * dq(IRHO,i,IX) - dq(IRHOU,i,IX)) + v * (v * dq(IRHO,i,IX) - dq(IRHOV,i,IX)) ) 
            dyT = gm1 * invRho * ( -q(IRHOE,i) * invRho * dq(IRHO,i,IY) + dq(IRHOE,i,IY) + u * ( u * dq(IRHO,i,IY) - dq(IRHOU,i,IY)) + v * (v * dq(IRHO,i,IY) - dq(IRHOV,i,IY)) ) 
      
            fv(IRHO,i)  = 0.0_RP
            fv(IRHOU,i) = 2.0_RP * ( - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX)) + lambda * divV
            fv(IRHOV,i) =  - vDivRho * dq(IRHO,i,IX)  &
                            - uDivRho * dq(IRHO,i,IY)  &
                            + invRho  * dq(IRHOU,i,IY) &
                            + invRho  * dq(IRHOV,i,IX) 
            fv(IRHOE,i) =  u * fv(IRHOU,i) + v * fv(IRHOV,i) + ( dimensionless % cp / dimensionless % Pr ) * dxT

            gv(IRHO,i)  = 0.0_RP
            gv(IRHOU,i) = fv(IRHOV,i)
            gv(IRHOV,i) =  2.0_RP * ( - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)) + lambda * divV
            gv(IRHOE,i) = u * gv(IRHOU,i) + v * gv(IRHOV,i) + ( dimensionless % cp / dimensionless % Pr ) * dyT
 
         end do

         end associate
      end function viscousFlux1D

      module pure function viscousFlux2D( N , q , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(1:NCONS,0:N,0:N)
         real(kind=RP), intent(in)          :: dq(1:NCONS,0:N,0:N,1:NDIM)
         real(kind=RP), target              :: F(1:NCONS,0:N,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i , j
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV , dxT , dyT
         real(kind=RP), pointer  :: fv(:,:,:) , gv(:,:,:)

         fv(1:,0:,0:)    => f(1:,0:,0:,IX)
         gv(1:,0:,0:)    => f(1:,0:,0:,IY)

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 

         do j = 0 , N   ; do i = 0 , N 
            invRho  = 1.0_RP / q(IRHO,i,j)
            u       = q(IRHOU,i,j) * invRho
            v       = q(IRHOV,i,j) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho
   
            divV = - uDivRho * dq(IRHO,i,j,IX) + invRho * dq(IRHOU,i,j,IX) - vDivRho * dq(IRHO,i,j,IY) + invRho * dq(IRHOV,i,j,IY)
   
            dxT = gm1 * invRho * ( -q(IRHOE,i,j) * invRho * dq(IRHO,i,j,IX) + dq(IRHOE,i,j,IX) + u * ( u * dq(IRHO,i,j,IX) - dq(IRHOU,i,j,IX)) + v * (v * dq(IRHO,i,j,IX) - dq(IRHOV,i,j,IX)) ) 
            dyT = gm1 * invRho * ( -q(IRHOE,i,j) * invRho * dq(IRHO,i,j,IY) + dq(IRHOE,i,j,IY) + u * ( u * dq(IRHO,i,j,IY) - dq(IRHOU,i,j,IY)) + v * (v * dq(IRHO,i,j,IY) - dq(IRHOV,i,j,IY)) ) 
      
            fv(IRHO,i,j)  = 0.0_RP
            fv(IRHOU,i,j) = 2.0_RP * ( - uDivRho * dq(IRHO,i,j,IX) + invRho * dq(IRHOU,i,j,IX)) + lambda * divV
            fv(IRHOV,i,j) =  - vDivRho * dq(IRHO,i,j,IX)  &
                            - uDivRho * dq(IRHO,i,j,IY)  &
                            + invRho  * dq(IRHOU,i,j,IY) &
                            + invRho  * dq(IRHOV,i,j,IX) 
            fv(IRHOE,i,j) =  u * fv(IRHOU,i,j) + v * fv(IRHOV,i,j) + ( dimensionless % cp / dimensionless % Pr ) * dxT

            gv(IRHO,i,j)  = 0.0_RP
            gv(IRHOU,i,j) = fv(IRHOV,i,j)
            gv(IRHOV,i,j) =  2.0_RP * ( - vDivRho * dq(IRHO,i,j,IY) + invRho * dq(IRHOV,i,j,IY)) + lambda * divV
            gv(IRHOE,i,j) = u * gv(IRHOU,i,j) + v * gv(IRHOV,i,j) + ( dimensionless % cp / dimensionless % Pr ) * dyT
 
         end do      ; end do

         end associate

      end function viscousFlux2D

      module pure function viscousFluxBC0D( q , qB , dq) result(F)
         implicit none
         real(kind=RP), intent(in)  :: q(NCONS)
         real(kind=RP), intent(in)  :: qB(NCONS)
         real(kind=RP), intent(in)  :: dq(NCONS,NDIM)
         real(kind=RP)              :: F(1:NCONS,NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV , dxT , dyT

         invRho  = 1.0_RP / q(IRHO)
         u       = q(IRHOU) * invRho
         v       = q(IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(IRHO,IX) + invRho * dq(IRHOU,IX) - vDivRho * dq(IRHO,IY) + invRho * dq(IRHOV,IY)

         dxT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IRHO,IX) + dq(IRHOE,IX) + u * ( u * dq(IRHO,IX) - dq(IRHOU,IX)) + v * (v * dq(IRHO,IX) - dq(IRHOV,IX)) ) 
         dyT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IRHO,IY) + dq(IRHOE,IY) + u * ( u * dq(IRHO,IY) - dq(IRHOU,IY)) + v * (v * dq(IRHO,IY) - dq(IRHOV,IY)) ) 
   
         F(IRHO,IX) = 0.0_RP
         F(IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(IRHO,IX) + invRho * dq(IRHOU,IX)) + lambda * divV
         F(IRHOV,IX) = - vDivRho * dq(IRHO,IX)  &
                       - uDivRho * dq(IRHO,IY)  &
                       + invRho  * dq(IRHOU,IY) &
                       + invRho  * dq(IRHOV,IX) 
         F(IRHOE,IX) = (qB(IRHOU) * F(IRHOU,IX) + qB(IRHOV) * F(IRHOU,IY))/qB(IRHO) + ( dimensionless % cp / dimensionless % Pr ) * dxT

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX)
         F(IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(IRHO,IY) + invRho * dq(IRHOV,IY)) + lambda * divV
         F(IRHOE,IY) = (qB(IRHOU) * F(IRHOV,IX) + qB(IRHOV) * F(IRHOV,IY))/qB(IRHO) + ( dimensionless % cp / dimensionless % Pr ) * dyT

         end associate

      end function viscousFluxBC0D

      module pure function viscousFluxBC1D( N , q , qb , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(1:NCONS,0:N)
         real(kind=RP), intent(in)          :: qB(1:NCONS,0:N)
         real(kind=RP), intent(in)          :: dq(1:NCONS,0:N,1:NDIM)
         real(kind=RP), target              :: F(1:NCONS,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i 
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV , dxT , dyT
         real(kind=RP), pointer  :: fv(:,:) , gv(:,:)

         fv(1:,0:)    => f(1:,0:,IX)
         gv(1:,0:)    => f(1:,0:,IY)

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 

         do i = 0 , N
            invRho  = 1.0_RP / q(IRHO,i)
            u       = q(IRHOU,i) * invRho
            v       = q(IRHOV,i) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho
   
            divV = - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX) - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)
   
            dxT = gm1 * invRho * ( -q(IRHOE,i) * invRho * dq(IRHO,i,IX) + dq(IRHOE,i,IX) + u * ( u * dq(IRHO,i,IX) - dq(IRHOU,i,IX)) + v * (v * dq(IRHO,i,IX) - dq(IRHOV,i,IX)) ) 
            dyT = gm1 * invRho * ( -q(IRHOE,i) * invRho * dq(IRHO,i,IY) + dq(IRHOE,i,IY) + u * ( u * dq(IRHO,i,IY) - dq(IRHOU,i,IY)) + v * (v * dq(IRHO,i,IY) - dq(IRHOV,i,IY)) ) 
      
            fv(IRHO,i)  = 0.0_RP
            fv(IRHOU,i) = 2.0_RP * ( - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX)) + lambda * divV
            fv(IRHOV,i) =  - vDivRho * dq(IRHO,i,IX)  &
                            - uDivRho * dq(IRHO,i,IY)  &
                            + invRho  * dq(IRHOU,i,IY) &
                            + invRho  * dq(IRHOV,i,IX) 
            fv(IRHOE,i) =  (qB(IRHOU,i) * fv(IRHOU,i) + qB(IRHOV,i) * fv(IRHOV,i))/qB(IRHO,i) + ( dimensionless % cp / dimensionless % Pr ) * dxT

            gv(IRHO,i)  = 0.0_RP
            gv(IRHOU,i) = fv(IRHOV,i)
            gv(IRHOV,i) =  2.0_RP * ( - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)) + lambda * divV
            gv(IRHOE,i) = (qB(IRHOU,i) * gv(IRHOU,i) + qB(IRHOV,i) * gv(IRHOV,i))/qB(IRHO,i) + ( dimensionless % cp / dimensionless % Pr ) * dyT
 
         end do

         end associate

      end function viscousFluxBC1D

      module pure function AdiabaticViscousFlux0D( q , qB , dq) result(F)
         implicit none
         real(kind=RP), intent(in)  :: q(NCONS)
         real(kind=RP), intent(in)  :: qB(NCONS)
         real(kind=RP), intent(in)  :: dq(NCONS,NDIM)
         real(kind=RP)              :: F(1:NCONS,NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV 

         invRho  = 1.0_RP / q(IRHO)
         u       = q(IRHOU) * invRho
         v       = q(IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(IRHO,IX) + invRho * dq(IRHOU,IX) - vDivRho * dq(IRHO,IY) + invRho * dq(IRHOV,IY)

         F(IRHO,IX) = 0.0_RP
         F(IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(IRHO,IX) + invRho * dq(IRHOU,IX)) + lambda * divV
         F(IRHOV,IX) = - vDivRho * dq(IRHO,IX)  &
                       - uDivRho * dq(IRHO,IY)  &
                       + invRho  * dq(IRHOU,IY) &
                       + invRho  * dq(IRHOV,IX) 
         F(IRHOE,IX) = (qB(IRHOU) * F(IRHOU,IX) + qB(IRHOV) * F(IRHOU,IY))/qB(IRHO)

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX)
         F(IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(IRHO,IY) + invRho * dq(IRHOV,IY)) + lambda * divV
         F(IRHOE,IY) = (qB(IRHOU) * F(IRHOV,IX) + qB(IRHOV) * F(IRHOV,IY))/qB(IRHO)

         end associate

      end function AdiabaticViscousFlux0D

      module pure function AdiabaticViscousFlux1D( N , q , qB , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(1:NCONS,0:N)
         real(kind=RP), intent(in)          :: qB(1:NCONS,0:N)
         real(kind=RP), intent(in)          :: dq(1:NCONS,0:N,1:NDIM)
         real(kind=RP), target              :: F(1:NCONS,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i 
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV 
         real(kind=RP), pointer  :: fv(:,:) , gv(:,:)

         fv(1:,0:)    => f(1:,0:,IX)
         gv(1:,0:)    => f(1:,0:,IY)

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 

         do i = 0 , N
            invRho  = 1.0_RP / q(IRHO,i)
            u       = q(IRHOU,i) * invRho
            v       = q(IRHOV,i) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho
   
            divV = - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX) - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)
   
            fv(IRHO,i)  = 0.0_RP
            fv(IRHOU,i) = 2.0_RP * ( - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX)) + lambda * divV
            fv(IRHOV,i) =  - vDivRho * dq(IRHO,i,IX)  &
                            - uDivRho * dq(IRHO,i,IY)  &
                            + invRho  * dq(IRHOU,i,IY) &
                            + invRho  * dq(IRHOV,i,IX) 
            fv(IRHOE,i) =  (qB(IRHOU,i) * fv(IRHOU,i) + qB(IRHOV,i) * fv(IRHOV,i))/qB(IRHO,i) 

            gv(IRHO,i)  = 0.0_RP
            gv(IRHOU,i) = fv(IRHOV,i)
            gv(IRHOV,i) =  2.0_RP * ( - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)) + lambda * divV
            gv(IRHOE,i) = (qB(IRHOU,i) * gv(IRHOU,i) + qB(IRHOV,i) * gv(IRHOV,i))/qB(IRHO,i)  
         end do

         end associate

      end function AdiabaticViscousFlux1D

      module pure function AdiabaticViscousFlux2D( N , q , qB , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(1:NCONS,0:N,0:N)
         real(kind=RP), intent(in)          :: qB(1:NCONS,0:N,0:N)
         real(kind=RP), intent(in)          :: dq(1:NCONS,0:N,0:N,1:NDIM)
         real(kind=RP), target              :: F(1:NCONS,0:N,0:N,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i , j
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV 
         real(kind=RP), pointer  :: fv(:,:,:) , gv(:,:,:)

         fv(1:,0:,0:)    => f(1:,0:,0:,IX)
         gv(1:,0:,0:)    => f(1:,0:,0:,IY)

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 

         do j = 0 , N   ; do i = 0 , N
            invRho  = 1.0_RP / q(IRHO,i,j)
            u       = q(IRHOU,i,j) * invRho
            v       = q(IRHOV,i,j) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho
   
            divV = - uDivRho * dq(IRHO,i,j,IX) + invRho * dq(IRHOU,i,j,IX) - vDivRho * dq(IRHO,i,j,IY) + invRho * dq(IRHOV,i,j,IY)
   
            fv(IRHO,i,j)  = 0.0_RP
            fv(IRHOU,i,j) = 2.0_RP * ( - uDivRho * dq(IRHO,i,j,IX) + invRho * dq(IRHOU,i,j,IX)) + lambda * divV
            fv(IRHOV,i,j) =  - vDivRho * dq(IRHO,i,j,IX)  &
                            - uDivRho * dq(IRHO,i,j,IY)  &
                            + invRho  * dq(IRHOU,i,j,IY) &
                            + invRho  * dq(IRHOV,i,j,IX) 
            fv(IRHOE,i,j) =  (qB(IRHOU,i,j) * fv(IRHOU,i,j) + qB(IRHOV,i,j) * fv(IRHOV,i,j))/qB(IRHO,i,j) 

            gv(IRHO,i,j)  = 0.0_RP
            gv(IRHOU,i,j) = fv(IRHOV,i,j)
            gv(IRHOV,i,j) =  2.0_RP * ( - vDivRho * dq(IRHO,i,j,IY) + invRho * dq(IRHOV,i,j,IY)) + lambda * divV
            gv(IRHOE,i,j) = (qB(IRHOU,i,j) * gv(IRHOU,i,j) + qB(IRHOV,i,j) * gv(IRHOV,i,j))/qB(IRHO,i,j)  
         end do      ; end do

         end associate

      end function AdiabaticViscousFlux2D

      module pure function ComputeViscousTensor ( N , Q , dQ ) result ( tau )
         implicit none
         integer,          intent(in)     :: N
         real(kind=RP),    intent(in)     ::  Q(1:NCONS,0:N)
         real(kind=RP),    intent(in)     :: dQ(1:NCONS,0:N,1:NDIM)
         real(kind=RP),    target         :: tau(0:N,1:NDIM,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i 
         real(kind=RP)          :: invRho , uDivRho , vDivRho , u , v , divV 
         real(kind=RP), pointer :: tauxx(:) , tauyy(:) , tauxy(:) , tauyx(:)

         tauxx(0:)   => tau(0:,IX,IX)
         tauxy(0:)   => tau(0:,IX,IY)
         tauyx(0:)   => tau(0:,IY,IX)
         tauyy(0:)   => tau(0:,IY,IY)

         associate ( mu => dimensionless % mu , lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 

         do i = 0 , N

            invRho  = 1.0_RP / q(IRHO,i)
            u       = q(IRHOU,i) * invRho
            v       = q(IRHOV,i) * invRho
            uDivRho = u * invRho
            vDivRho = v * invRho
   
            divV = - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX) - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)
      
            tauxx(i) = 2.0_RP * mu * ( - uDivRho * dq(IRHO,i,IX) + invRho * dq(IRHOU,i,IX)) + lambda * mu * divV
            tauyy(i) = 2.0_RP * mu * ( - vDivRho * dq(IRHO,i,IY) + invRho * dq(IRHOV,i,IY)) + lambda * mu * divV
   
            tauxy(i) = - mu * vDivRho * dq(IRHO,i,IX)  &
                       - mu * uDivRho * dq(IRHO,i,IY)  &
                       + mu * invRho  * dq(IRHOU,i,IY) &
                       + mu * invRho  * dq(IRHOV,i,IX) 
   
            tauyx(i) = tauxy(i)
   
         end do
         end associate

      end function ComputeViscousTensor

end submodule ViscousFluxes
