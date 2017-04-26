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
         real(kind=RP), intent(in)  :: dq(NDIM , NCONS)
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


         divV = - uDivRho * dq(IX,IRHO) + invRho * dq(IX,IRHOU) - vDivRho * dq(IY,IRHO) + invRho * dq(IY,IRHOV)

         dxT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IX,IRHO) + dq(IX,IRHOE) + u * ( u * dq(IX,IRHO) - dq(IX,IRHOU)) + v * (v * dq(IX,IRHO) - dq(IX,IRHOV)) ) 
         dyT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IY,IRHO) + dq(IY,IRHOE) + u * ( u * dq(IY,IRHO) - dq(IY,IRHOU)) + v * (v * dq(IY,IRHO) - dq(IY,IRHOV)) ) 
   
         F(IRHO,:) = 0.0_RP
         F(IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(IX,IRHO) + invRho * dq(IX,IRHOU)) + lambda * divV
         F(IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(IY,IRHO) + invRho * dq(IY,IRHOV)) + lambda * divV

         F(IRHOU,IY) = - vDivRho * dq(IX,IRHO)  &
                       - uDivRho * dq(IY,IRHO)  &
                       + invRho  * dq(IY,IRHOU) &
                       + invRho  * dq(IX,IRHOV) 


         F(IRHOV,IX) = F(IRHOU,IY)

         F(IRHOE,IX) = u * F(IRHOU,IX) + v * F(IRHOU,IY) + ( dimensionless % cp / dimensionless % Pr ) * dxT
         F(IRHOE,IY) = u * F(IRHOV,IX) + v * F(IRHOV,IY) + ( dimensionless % cp / dimensionless % Pr ) * dyT

         end associate

      end function viscousFlux0D

      module pure function viscousFlux1D( N , q , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,1:NCONS,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N) , uDivRho(0:N) , vDivRho(0:N) , u(0:N) , v(0:N) , divV(0:N) , dxT(0:N) , dyT(0:N)

         invRho  = 1.0_RP / q(:,IRHO)
         u       = q(:,IRHOU) * invRho
         v       = q(:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU) - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)


         dxT = gm1 * invRho * ( -q(:,IRHOE) * invRho * dq(:,IX,IRHO) + dq(:,IX,IRHOE) + u * ( u * dq(:,IX,IRHO) - dq(:,IX,IRHOU)) + v * (v * dq(:,IX,IRHO) - dq(:,IX,IRHOV)) ) 
         dyT = gm1 * invRho * ( -q(:,IRHOE) * invRho * dq(:,IY,IRHO) + dq(:,IY,IRHOE) + u * ( u * dq(:,IY,IRHO) - dq(:,IY,IRHOU)) + v * (v * dq(:,IY,IRHO) - dq(:,IY,IRHOV)) ) 
   
         F(:,IRHO,:)   = 0.0_RP
         F(:,IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU)) + lambda * divV
         F(:,IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)) + lambda * divV

         F(:,IRHOU,IY) = - vDivRho * dq(:,IX,IRHO)  &
                         - uDivRho * dq(:,IY,IRHO)  &
                         + invRho  * dq(:,IY,IRHOU) &
                         + invRho  * dq(:,IX,IRHOV) 


         F(:,IRHOV,IX) = F(:,IRHOU,IY)

         F(:,IRHOE,IX) = u * F(:,IRHOU,IX) + v * F(:,IRHOU,IY) + ( dimensionless % cp / dimensionless % Pr ) * dxT
         F(:,IRHOE,IY) = u * F(:,IRHOV,IX) + v * F(:,IRHOV,IY) + ( dimensionless % cp / dimensionless % Pr ) * dyT

         end associate

      end function viscousFlux1D

      module pure function viscousFlux2D( N , q , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,0:N,1:NCONS,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N,0:N) , uDivRho(0:N,0:N) , vDivRho(0:N,0:N) , u(0:N,0:N) , v(0:N,0:N) , divV(0:N,0:N) , dxT(0:N,0:N) , dyT(0:N,0:N)

         invRho  = 1.0_RP / q(:,:,IRHO)
         u       = q(:,:,IRHOU) * invRho
         v       = q(:,:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(:,:,IX,IRHO) + invRho * dq(:,:,IX,IRHOU) - vDivRho * dq(:,:,IY,IRHO) + invRho * dq(:,:,IY,IRHOV)

         dxT = gm1 * invRho * ( -q(:,:,IRHOE) * invRho * dq(:,:,IX,IRHO) + dq(:,:,IX,IRHOE) + u * ( u * dq(:,:,IX,IRHO) - dq(:,:,IX,IRHOU)) + v * (v * dq(:,:,IX,IRHO) - dq(:,:,IX,IRHOV)) ) 
         dyT = gm1 * invRho * ( -q(:,:,IRHOE) * invRho * dq(:,:,IY,IRHO) + dq(:,:,IY,IRHOE) + u * ( u * dq(:,:,IY,IRHO) - dq(:,:,IY,IRHOU)) + v * (v * dq(:,:,IY,IRHO) - dq(:,:,IY,IRHOV)) ) 

         F(:,:,IRHO,:)   = 0.0_RP
         F(:,:,IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(:,:,IX,IRHO) + invRho * dq(:,:,IX,IRHOU)) + lambda * divV
         F(:,:,IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(:,:,IY,IRHO) + invRho * dq(:,:,IY,IRHOV)) + lambda * divV

         F(:,:,IRHOU,IY) = - vDivRho * dq(:,:,IX,IRHO)  &
                           - uDivRho * dq(:,:,IY,IRHO)  &
                           + invRho  * dq(:,:,IY,IRHOU) &
                           + invRho  * dq(:,:,IX,IRHOV) 
 

         F(:,:,IRHOV,IX) = F(:,:,IRHOU,IY)

         F(:,:,IRHOE,IX) = u * F(:,:,IRHOU,IX) + v * F(:,:,IRHOU,IY) + ( dimensionless % cp / dimensionless % Pr ) * dxT
         F(:,:,IRHOE,IY) = u * F(:,:,IRHOV,IX) + v * F(:,:,IRHOV,IY) + ( dimensionless % cp / dimensionless % Pr ) * dyT

         end associate

      end function viscousFlux2D

      module pure function viscousFluxBC0D( q , qB , dq) result(F)
         implicit none
         real(kind=RP), intent(in)  :: q(NCONS)
         real(kind=RP), intent(in)  :: qB(NCONS)
         real(kind=RP), intent(in)  :: dq(NDIM , NCONS)
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


         divV = - uDivRho * dq(IX,IRHO) + invRho * dq(IX,IRHOU) - vDivRho * dq(IY,IRHO) + invRho * dq(IY,IRHOV)

         dxT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IX,IRHO) + dq(IX,IRHOE) + u * ( u * dq(IX,IRHO) - dq(IX,IRHOU)) + v * (v * dq(IX,IRHO) - dq(IX,IRHOV)) ) 
         dyT = gm1 * invRho * ( -q(IRHOE) * invRho * dq(IY,IRHO) + dq(IY,IRHOE) + u * ( u * dq(IY,IRHO) - dq(IY,IRHOU)) + v * (v * dq(IY,IRHO) - dq(IY,IRHOV)) ) 
   

         F(IRHO,:)   = 0.0_RP
         F(IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(IX,IRHO) + invRho * dq(IX,IRHOU)) + lambda * divV
         F(IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(IY,IRHO) + invRho * dq(IY,IRHOV)) + lambda * divV

         F(IRHOU,IY) = - vDivRho * dq(IX,IRHO)  &
                       - uDivRho * dq(IY,IRHO)  &
                       + invRho  * dq(IY,IRHOU) &
                       + invRho  * dq(IX,IRHOV) 


         F(IRHOV,IX) = F(IRHOU,IY)

         F(IRHOE,IX) = ( qB(IRHOU) * F(IRHOU,IX) + qB(IRHOV) * F(IRHOU,IY) ) / qB(IRHO) + ( dimensionless % cp / dimensionless % Pr ) * dxT
         F(IRHOE,IY) = ( qB(IRHOU) * F(IRHOV,IX) + qB(IRHOV) * F(IRHOV,IY) ) / qB(IRHO) + ( dimensionless % cp / dimensionless % Pr )  * dyT

         end associate

      end function viscousFluxBC0D

      module pure function viscousFluxBC1D( N , q , qb , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: qB(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,1:NCONS,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N) , uDivRho(0:N) , vDivRho(0:N) , u(0:N) , v(0:N) , divV(0:N) , dxT(0:N) , dyT(0:N)

         invRho  = 1.0_RP / q(:,IRHO)
         u       = q(:,IRHOU) * invRho
         v       = q(:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU) - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)


         dxT = gm1 * invRho * ( -q(:,IRHOE) * invRho * dq(:,IX,IRHO) + dq(:,IX,IRHOE) + u * ( u * dq(:,IX,IRHO) - dq(:,IX,IRHOU)) + v * (v * dq(:,IX,IRHO) - dq(:,IX,IRHOV)) ) 
         dyT = gm1 * invRho * ( -q(:,IRHOE) * invRho * dq(:,IY,IRHO) + dq(:,IY,IRHOE) + u * ( u * dq(:,IY,IRHO) - dq(:,IY,IRHOU)) + v * (v * dq(:,IY,IRHO) - dq(:,IY,IRHOV)) ) 
   
         F(:,IRHO,:)   = 0.0_RP
         F(:,IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU)) + lambda * divV
         F(:,IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)) + lambda * divV

         F(:,IRHOU,IY) = - vDivRho * dq(:,IX,IRHO)  &
                         - uDivRho * dq(:,IY,IRHO)  &
                         + invRho  * dq(:,IY,IRHOU) &
                         + invRho  * dq(:,IX,IRHOV) 


         F(:,IRHOV,IX) = F(:,IRHOU,IY)

         F(:,IRHOE,IX) = (qB(:,IRHOU) * F(:,IRHOU,IX) + qB(:,IRHOV) * F(:,IRHOU,IY)) / qB(:,IRHO) + ( dimensionless % cp / dimensionless % Pr ) * dxT
         F(:,IRHOE,IY) = (qB(:,IRHOU) * F(:,IRHOV,IX) + qB(:,IRHOV) * F(:,IRHOV,IY)) / qB(:,IRHO) + ( dimensionless % cp / dimensionless % Pr ) * dyT

         end associate

      end function viscousFluxBC1D

      module pure function AdiabaticViscousFlux0D( q , qB , dq) result(F)
         implicit none
         real(kind=RP), intent(in)  :: q(NCONS)
         real(kind=RP), intent(in)  :: qB(NCONS)
         real(kind=RP), intent(in)  :: dq(NDIM , NCONS)
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


         divV = - uDivRho * dq(IX,IRHO) + invRho * dq(IX,IRHOU) - vDivRho * dq(IY,IRHO) + invRho * dq(IY,IRHOV)

         F(IRHO,:)   = 0.0_RP
         F(IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(IX,IRHO) + invRho * dq(IX,IRHOU)) + lambda * divV
         F(IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(IY,IRHO) + invRho * dq(IY,IRHOV)) + lambda * divV

                                                               !                      -                                             -
         F(IRHOU,IY) = - vDivRho * dq(IX,IRHO)  &         !        du   dv      |    u  drho    1  drhou    v  drho    1  drhov |
                       - uDivRho * dq(IY,IRHO)  &         !     mu -- + -- = mu | - --- ---- + --- ----- - --- ---- + --- ----- |
                       + invRho  * dq(IY,IRHOU) &         !        dy   dx      |   rho  dy    rho  dy     rho  dx    rho  dx   |
                       + invRho  * dq(IX,IRHOV)           !                      -                                             -


         F(IRHOV,IX) = F(IRHOU,IY)

         F(IRHOE,IX) = (qB(IRHOU) * F(IRHOU,IX) + qB(IRHOV) * F(IRHOU,IY)) / qB(IRHO) 
         F(IRHOE,IY) = (qB(IRHOU) * F(IRHOV,IX) + qB(IRHOV) * F(IRHOV,IY)) / qB(IRHO) 

         end associate

      end function AdiabaticViscousFlux0D

      module pure function AdiabaticViscousFlux1D( N , q , qB , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: qB(0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,1:NCONS,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N) , uDivRho(0:N) , vDivRho(0:N) , u(0:N) , v(0:N) , divV(0:N)

         invRho  = 1.0_RP / q(:,IRHO)
         u       = q(:,IRHOU) * invRho
         v       = q(:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU) - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)

         F(:,IRHO ,: ) = 0.0_RP
         F(:,IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU)) + lambda * divV
         F(:,IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)) + lambda * divV

         F(:,IRHOU,IY) = - vDivRho * dq(:,IX,IRHO)  &
                         - uDivRho * dq(:,IY,IRHO)  &
                         + invRho  * dq(:,IY,IRHOU) &
                         + invRho  * dq(:,IX,IRHOV) 


         F(:,IRHOV,IX) = F(:,IRHOU,IY)

         F(:,IRHOE,IX) = (qB(:,IRHOU) * F(:,IRHOU,IX) + qB(:,IRHOV) * F(:,IRHOU,IY) ) / qB(:,IRHO)
         F(:,IRHOE,IY) = (qB(:,IRHOU) * F(:,IRHOV,IX) + qB(:,IRHOV) * F(:,IRHOV,IY) ) / qB(:,IRHO)

         end associate

      end function AdiabaticViscousFlux1D

      module pure function AdiabaticViscousFlux2D( N , q , qB , dq ) result ( F )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP), intent(in)          :: q(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)          :: qB(0:N,0:N,1:NCONS)
         real(kind=RP), intent(in)          :: dq(0:N,0:N,1:NDIM,1:NCONS)
         real(kind=RP)                      :: F(0:N,0:N,1:NCONS,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N,0:N) , uDivRho(0:N,0:N) , vDivRho(0:N,0:N) , u(0:N,0:N) , v(0:N,0:N) , divV(0:N,0:N)

         invRho  = 1.0_RP / q(:,:,IRHO)
         u       = q(:,:,IRHOU) * invRho
         v       = q(:,:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(:,:,IX,IRHO) + invRho * dq(:,:,IX,IRHOU) - vDivRho * dq(:,:,IY,IRHO) + invRho * dq(:,:,IY,IRHOV)

         F(:,:,IRHO , :) = 0.0_RP
         F(:,:,IRHOU,IX) = 2.0_RP * ( - uDivRho * dq(:,:,IX,IRHO) + invRho * dq(:,:,IX,IRHOU)) + lambda * divV
         F(:,:,IRHOV,IY) = 2.0_RP * ( - vDivRho * dq(:,:,IY,IRHO) + invRho * dq(:,:,IY,IRHOV)) + lambda * divV

         F(:,:,IRHOU,IY) = - vDivRho * dq(:,:,IX,IRHO)  &
                           - uDivRho * dq(:,:,IY,IRHO)  &
                           + invRho  * dq(:,:,IY,IRHOU) &
                           + invRho  * dq(:,:,IX,IRHOV) 


         F(:,:,IRHOV,IX) = F(:,:,IRHOU,IY)

         F(:,:,IRHOE,IX) = ( qB(:,:,IRHOU) * F(:,:,IRHOU,IX) + qB(:,:,IRHOV) * F(:,:,IRHOU,IY) ) / qB(:,:,IRHO)
         F(:,:,IRHOE,IY) = ( qB(:,:,IRHOU) * F(:,:,IRHOV,IX) + qB(:,:,IRHOV) * F(:,:,IRHOV,IY) ) / qB(:,:,IRHO) 

         end associate

      end function AdiabaticViscousFlux2D

      module pure function ComputeViscousTensor ( N , Q , dQ ) result ( tau )
         implicit none
         integer,          intent(in)     :: N
         real(kind=RP),    intent(in)     ::  Q(0:N,1:NCONS)
         real(kind=RP),    intent(in)     :: dQ(0:N,1:NDIM,1:NCONS)
         real(kind=RP)                    :: tau(0:N,1:NDIM,1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: invRho(0:N) , uDivRho(0:N) , vDivRho(0:N) , u(0:N) , v(0:N) , divV(0:N) 

         invRho  = 1.0_RP / q(:,IRHO)
         u       = q(:,IRHOU) * invRho
         v       = q(:,IRHOV) * invRho
         uDivRho = u * invRho
         vDivRho = v * invRho

         associate ( mu => dimensionless % mu , lambda => thermodynamics % lambda , gm1 => thermodynamics % gm1 ) 


         divV = - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU) - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)
   
         tau(:,IX,IX) = 2.0_RP * mu * ( - uDivRho * dq(:,IX,IRHO) + invRho * dq(:,IX,IRHOU)) + lambda * mu * divV
         tau(:,IY,IY) = 2.0_RP * mu * ( - vDivRho * dq(:,IY,IRHO) + invRho * dq(:,IY,IRHOV)) + lambda * mu * divV

         tau(:,IX,IY) = - mu * vDivRho * dq(:,IX,IRHO)  &
                         - mu * uDivRho * dq(:,IY,IRHO)  &
                         + mu * invRho  * dq(:,IY,IRHOU) &
                         + mu * invRho  * dq(:,IX,IRHOV) 

         tau(:,IY,IX) = tau(:,IX,IY)

         end associate

      end function ComputeViscousTensor

end submodule ViscousFluxes
