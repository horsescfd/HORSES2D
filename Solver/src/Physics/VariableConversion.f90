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

#ifdef _DIMENSIONLESS_TAU
         T = getPressure0D(Q) / Q(IRHO)
#else
         T = dimensionless % gammaMach2 * getPressure0D(Q) / Q(IRHO)
#endif

      end function getTemperature0D

      module function getTemperature1D(N,Q) result (T)
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

      module function getTemperature2D(N,Q) result (T)
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
      module function getSoundSpeed0D(Q) result (a)
         implicit none
         real(kind=RP), intent(in)        :: Q(1:NCONS)
         real(kind=RP)                    :: a

         a = sqrt( thermodynamics % gamma * getPressure0D(Q) / Q(IRHO) )

      end function getSoundSpeed0D

      module function getSoundSpeed1D(N,Q) result (a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(0:N,1:NCONS)
         real(kind=RP)                    :: a(0:N)

         a = sqrt( thermodynamics % gamma * getPressure1D(N,Q) / Q(:,IRHO) )

      end function getSoundSpeed1D

      module function getSoundSpeed2D(N,Q) result(a)
         implicit none
         integer,       intent(in)        :: N
         real(kind=RP), intent(in)        :: Q(0:N,0:N,1:NCONS)
         real(kind=RP)                    :: a(0:N,0:N)

         a = sqrt( thermodynamics % gamma * getPressure2D(N,Q) / Q(:,:,IRHO) )

      end function getSoundSpeed2D
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
 

end submodule VariableConversion
