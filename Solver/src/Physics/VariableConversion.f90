submodule(PhysicsNS)   VariableConversion


   contains
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!           Compute pressure
!
!
      module function getPressure0D(Q) result (p)
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)                       :: p

         p = thermodynamics % gm1 * (Q(IRHOE) - 0.5_RP * ( Q(IRHOU)*Q(IRHOU) + Q(IRHOV)*Q(IRHOV) ) / Q(IRHO) )

      end function getPressure0D

      module function getPressure1D(N,Q) result (p)
         implicit none
         integer,       intent (in)  :: N
         real(kind=RP), intent (in)  :: Q(0:N,1:NCONS)
         real(kind=RP)              :: p(0:N)

         p = thermodynamics % gm1 * (Q(:,IRHOE) - 0.5_RP * ( Q(:,IRHOU)*Q(:,IRHOU) + Q(:,IRHOV)*Q(:,IRHOV) ) / Q(:,IRHO) )

      end function getPressure1D

      module function getPressure2D(N,Q) result (p)
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
      module function getTemperature0D(Q) result (T)
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)                      :: T

         T = dimensionless % gammaMach2 * getPressure0D(Q) / Q(IRHO)

      end function getTemperature0D

      module function getTemperature1D(N,Q) result (T)
         implicit none
         integer,       intent (in)  :: N
         real(kind=RP), intent (in)  :: Q(0:N,1:NCONS)
         real(kind=RP)              :: T(0:N)

         T = dimensionless % gammaMach2 * getPressure0D(Q) / Q(:,IRHO)

      end function getTemperature1D

      module function getTemperature2D(N,Q) result (T)
         implicit none
         integer,       intent(in)     :: N 
         real(kind=RP), intent (in)  :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)              :: T(0:N,0:N)

         T = dimensionless % gammaMach2 * getPressure0D(Q) / Q(:,:,IRHO)

      end function getTemperature2D

end submodule VariableConversion
