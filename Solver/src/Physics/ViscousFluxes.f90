submodule (PhysicsNS)   ViscousFluxes
   use SMConstants
   implicit none

   contains


      module function viscousFlux0D( w , dq) result(val)
         implicit none
         real(kind=RP)          :: w(NPRIM)
         real(kind=RP)          :: dq(NDIM , NGRAD)
         real(kind=RP), target  :: val(NCONS,NDIM)
         real(kind=RP), pointer :: F(:) , G(:)

         F(1:NCONS)    => val(1:NCONS,IX)
         G(1:NCONS)    => val(1:NCONS,IY)

         associate ( mu => dimensionless % mu , lambda => thermodynamics % lambda , kappa => dimensionless % kappa ) 
           
         F(IRHO)  = 0.0_RP
         F(IRHOU) = mu * ( 2.0_RP * dq(IX , IGU) - lambda * ( dq(IX,IGU) + dq(IY,IGV) ) )
         F(IRHOV) = mu * ( dq(IY,IGU) + dq(IX,IGV) )
         F(IRHOE) = F(IRHOU) * w(IU) + F(IRHOV) * w(IV) + kappa * dq(IX,IGT) 

         G(IRHO)  = 0.0_RP
         G(IRHOU) = F(IRHOV)
         G(IRHOV) = mu * ( 2.0_RP * dq(IY , IGV) - lambda * ( dq(IX,IGU) + dq(IY,IGV) ) )
         G(IRHOE) = G(IRHOU) * w(IU) + G(IRHOV) * w(IV) + kappa * dq(IY,IGT) 

         end associate

      end function viscousFlux0D

      module function viscousFlux1D( N , w , dq ) result ( val )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP)                      :: w(0:N,1:NPRIM)
         real(kind=RP)                      :: dq(0:N,1:NDIM,1:NGRAD)
         real(kind=RP), target              :: val(0:N,1:NCONS,1:NDIM)
         real(kind=RP), pointer             :: F(:,:) , G(:,:)

         F(0:,1:)    => val(0:,1:,IX)
         G(0:,1:)    => val(0:,1:,IY)

         associate ( mu => dimensionless % mu , lambda => thermodynamics % lambda , kappa => dimensionless % kappa ) 

         F(:,IRHO)  = 0.0_RP
         F(:,IRHOU) = mu * ( 2.0_RP * dq(:,IX , IGU) - lambda * ( dq(:,IX,IGU) + dq(:,IY,IGV) ) )
         F(:,IRHOV) = mu * ( dq(:,IY,IGU) + dq(:,IX,IGV) )
         F(:,IRHOE) = F(:,IRHOU) * w(:,IU) + F(:,IRHOV) * w(:,IV) + kappa * dq(:,IX,IGT) 

         G(:,IRHO)  = 0.0_RP
         G(:,IRHOU) = F(:,IRHOV)
         G(:,IRHOV) = mu * ( 2.0_RP * dq(:,IY , IGV) - lambda * ( dq(:,IX,IGU) + dq(:,IY,IGV) ) )
         G(:,IRHOE) = G(:,IRHOU) * w(:,IU) + G(:,IRHOV) * w(:,IV) + kappa * dq(:,IY,IGT) 

         end associate 

      end function viscousFlux1D

      module function viscousFlux2D( N , w , dq ) result ( val )
         implicit none
         integer, intent(in)                :: N 
         real(kind=RP)                      :: w(0:N,0:N,1:NPRIM)
         real(kind=RP)                      :: dq(0:N,0:N,1:NDIM,1:NGRAD)
         real(kind=RP), target              :: val(0:N,0:N,1:NCONS,1:NDIM)
         real(kind=RP), pointer             :: F(:,:,:) , G(:,:,:)

         F(0:,0:,1:)    => val(0:,0:,1:,iX)
         G(0:,0:,1:)    => val(0:,0:,1:,iY)

         associate ( mu => dimensionless % mu , lambda => thermodynamics % lambda , kappa => dimensionless % kappa ) 

         F(:,:,IRHO)  = 0.0_RP
         F(:,:,IRHOU) = mu * ( 2.0_RP * dq(:,:,IX , IGU) - lambda * ( dq(:,:,IX,IGU) + dq(:,:,IY,IGV) ) )
         F(:,:,IRHOV) = mu * ( dq(:,:,IY,IGU) + dq(:,:,IX,IGV) )
         F(:,:,IRHOE) = F(:,:,IRHOU) * w(:,:,IU) + F(:,:,IRHOV) * w(:,:,IV) + kappa * dq(:,:,IX,IGT) 

         G(:,:,IRHO)  = 0.0_RP
         G(:,:,IRHOU) = F(:,:,IRHOV)
         G(:,:,IRHOV) = mu * ( 2.0_RP * dq(:,:,IY , IGV) - lambda * ( dq(:,:,IX,IGU) + dq(:,:,IY,IGV) ) )
         G(:,:,IRHOE) = G(:,:,IRHOU) * w(:,:,IU) + G(:,:,IRHOV) * w(:,:,IV) + kappa * dq(:,:,IY,IGT) 

         end associate 

      end function viscousFlux2D

      module function viscousNormalFlux0D( w , dq , dS ) result ( val )
         implicit none
         real(kind=RP), intent(in)  :: w(NPRIM)
         real(kind=RP), intent(in)  :: dq(NDIM,NGRAD)
         real(kind=RP), intent(in)  :: dS(NDIM) 
         real(kind=RP)              :: val(NCONS)
!        ----------------------------------------------------------
         real(kind=RP)              :: Fv(NCONS , NDIM)

         Fv = viscousFlux0D ( w , dq )

         val = Fv(1:NCONS,IX) * dS(IX) + Fv(1:NCONS,IY) * dS(IY)

      end function viscousNormalFlux0D

      module function viscousNormalFlux1D ( N , w , dq , dS ) result ( val )
         implicit none
         integer, intent(in)        :: N
         real(kind=RP), intent(in)  :: w(0:N ,1:NPRIM)
         real(kind=RP), intent(in)  :: dq(0:N , 1:NDIM , 1:NGRAD )
         real(kind=RP), intent(in)  :: dS(1:NDIM , 0:N )
         real(kind=RP)              :: val(0:N , 1:NCONS)
!        ---------------------------------------------------------
         real(kind=RP)              :: Fv(0:N,1:NCONS,1:NDIM)
         integer                    :: eq

         Fv = viscousFlux1D ( N , w , dq )
   
         do eq = 1 , NCONS
            val(:,eq) = Fv(:,eq,IX) * dS(IX,:) + Fv(:,eq,IY) * dS(IY,:)
         end do

      end function viscousNormalFlux1D

      module function viscousNormalFlux2D ( N , w , dq , dS ) result ( val )
         implicit none
         integer, intent(in)        :: N 
         real(kind=RP), intent(in)  :: w(0:N , 0:N , 1:NPRIM)
         real(kind=RP), intent(in)  :: dq(0:N , 0:N , 1:NDIM , 1:NGRAD )
         real(kind=RP), intent(in)  :: dS(1:NDIM , 0:N , 0:N )
         real(kind=RP)              :: val(0:N , 0:N , 1:NCONS)
!        ---------------------------------------------------------
         real(kind=RP)              :: Fv(0:N,0:N,1:NCONS,1:NDIM)
         integer                    :: eq

         Fv = viscousFlux2D ( N , w , dq )
   
         do eq = 1 , NCONS
            val(:,:,eq) = Fv(:,:,eq,IX) * dS(IX,:,:) + Fv(:,:,eq,IY) * dS(IY,:,:)
         end do

      end function viscousNormalFlux2D

      module function ComputeViscousTensor ( N , dQ ) result ( tau )
         implicit none
         integer,          intent(in)     :: N
         real(kind=RP),    intent(in)     :: dQ(0:N,1:NDIM,1:NGRAD)
         real(kind=RP)                    :: tau(0:N,1:NDIM,1:NDIM)

         associate ( mu => dimensionless % mu , lambda => thermodynamics % lambda )

         tau(0:N , IX , IX ) = mu * ( 2.0_RP * dQ(0:N , IX , IGU ) - lambda * ( dQ(0:N,IX,IGU) + dQ(0:N,IY,IGV) ) )
         tau(0:N , IY , IY ) = mu * ( 2.0_RP * dQ(0:N , IY , IGV ) - lambda * ( dQ(0:N,IX,IGU) + dQ(0:N,IY,IGV) ) )
         tau(0:N , IX , IY ) = mu * ( dQ(0:N,IY,IGU) + dQ(0:N,IX,IGV) )
         tau(0:N , IY , IX ) = tau(0:N , IX , IY )

         end associate

      end function ComputeViscousTensor














end submodule ViscousFluxes
