submodule (PhysicsNS)  RiemannSolvers
   use SMConstants
   contains

#include "Defines.h"

!
!     ****************************************************
!        Riemann solvers
!     ****************************************************
!
      module pure function ExactRiemannSolver(qL , qR , n) result (Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS), intent(in) :: qL
         real(kind=RP), dimension(NCONS), intent(in) :: qR
         real(kind=RP), dimension(NDIM) , intent(in) :: n
         real(kind=RP), dimension(NCONS) :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP), dimension(NCONS) :: qnL , qnR
         real(kind=RP), dimension(NPRIM) :: wL , wR
         real(kind=RP)                   :: pstar , ustar
         real(kind=RP)                   :: rhostar , uFan , pFan
         real(kind=RP)                   :: Ft

!        0/ Gather variables
!           ----------------
            qnL(IRHO)  = qL(IRHO)
            qnL(IRHOU) = qL(IRHOU) * n(IX) + qL(IRHOV) * n(IY)
            qnL(IRHOV) = -qL(IRHOU) * n(IY) + qL(IRHOV) * n(IX)
            qnL(IRHOE) = qL(IRHOE) 

            qnR(IRHO)  = qR(IRHO)
            qnR(IRHOU) = qR(IRHOU) * n(IX) + qR(IRHOV) * n(IY)
            qnR(IRHOV) = -qR(IRHOU) * n(IY) + qR(IRHOV) * n(IX)
            qnR(IRHOE) = qR(IRHOE) 


            associate( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )
            wL(IRHO) = qnL(IRHO)
            wL(IU)   = qnL(IRHOU) / qnL(IRHO)
            wL(IV)   = qnL(IRHOV) / qnL(IRHO)
            wL(IP)   = gm1 * (qnL(IRHOE) - 0.5_RP * (qnL(IRHOU) * wL(IU) + qnL(IRHOV) * wL(IV) ) )
            wL(IA)   = sqrt(gamma * wL(IP) / wL(IRHO))

            wR(IRHO) = qnR(IRHO)
            wR(IU)   = qnR(IRHOU) / qnR(IRHO)
            wR(IV)   = qnR(IRHOV) / qnR(IRHO)
            wR(IP)   = gm1 * (qnR(IRHOE) - 0.5_RP * (qnR(IRHOU) * wR(IU) + qnR(IRHOV) * wR(IV) ) )
            wR(IA)   = sqrt(gamma * wR(IP) / wR(IRHO))
            end associate

!        1/ Compute the star region
!           -----------------------
            call ExactRiemann_ComputePStar (wL , wR , pstar , ustar )

!        2/ Check to which region belongs the solution
!           ------------------------------------------
            associate ( cv => Dimensionless % cv , cp => Dimensionless % cp , gamma => Thermodynamics % gamma ) 
            if ( ustar .ge. 0.0_RP ) then

               if ( ( pstar .le. wL(IP) ) .and. ( wL(IU) .ge. wL(IA) ) ) then
                  Fstar(IRHO)  = qL(IRHOU)
                  Fstar(IRHOU) = qL(IRHOU) * wL(IU) + wL(IP)
                  Fstar(IRHOV) = qL(IRHOV) * wL(IU)
                  Fstar(IRHOE) = wL(IU)*( cp * wL(IP) + 0.5_RP * (qL(IRHOU)*wL(IU) + qL(IRHOV)*wL(IV)) )

               elseif ( ( pstar .le. wL(IP) ) .and. ( wL(IU) .lt. wL(IA)) .and. (ustar .lt. wL(IA) * (pstar/wL(IP))**(0.5_RP / cp ) ) ) then
                  rhostar = wL(IRHO) * ( pstar / wL(IP) ) ** ( 1.0_RP / gamma )

                  Fstar(IRHO) = rhostar * ustar 
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = Fstar(IRHO) * wL(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + wL(IV) * wL(IV)) ) 

               elseif ( ( pstar .le. wL(IP) ) .and. ( wL(IU) .lt. wL(IA) ) .and. ( ustar .ge. wL(IA) * (pstar / wL(IP))**(0.5_RP / cp))) then
                  
                  uFan    = 2.0_RP * wL(IA) / ( gamma + 1.0_RP)  + (gamma-1.0_RP)/(gamma+1.0_RP) * wL(IU)
                  rhostar = wL(IRHO) * (uFan / wL(IA)) ** (2.0_RP * cv)
                  pFan    = wL(IP) * (rhostar / wL(IRHO)) ** (gamma)

                  Fstar(IRHO) = rhostar * uFan
                  Fstar(IRHOU) = rhostar * uFan * uFan + pFan
                  Fstar(IRHOV) = rhostar * uFan * wL(IV)
                  Fstar(IRHOE) = uFan * ( cp * pFan + 0.5_RP * rhostar * ( uFan*uFan + wL(IV)*wL(IV) ) )
 
               elseif ( ( pstar .gt. wL(IP) ) .and. ( wL(IU) .ge. wL(IA) * sqrt( (gamma+1.0_RP)/(2.0_RP * gamma)*pstar/wL(IP) + 0.5_RP/cp) ) ) then
                  Fstar(IRHO)  = qL(IRHOU)
                  Fstar(IRHOU) = qL(IRHOU) * wL(IU) + wL(IP)
                  Fstar(IRHOV) = qL(IRHOV) * wL(IU)
                  Fstar(IRHOE) = wL(IU)*( cp * wL(IP) + 0.5_RP * (qL(IRHOU)*wL(IU) + qL(IRHOV)*wL(IV)) )

               elseif ( ( pstar .gt. wL(IP) ) .and. ( wL(IU) .lt.  wL(IA) * sqrt( (gamma+1.0_RP)/(2.0_RP * gamma)*pstar/wL(IP) + 0.5_RP/cp) ) ) then
                  rhostar = wL(IRHO) * ( (pstar/wL(IP) + (gamma-1.0_RP)/(gamma+1.0_RP))/(pstar/wL(IP)*(gamma-1.0_RP)/(gamma+1.0_RP) + 1.0_RP))
                  Fstar(IRHO) = rhostar * ustar
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = rhostar * ustar * wL(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + wL(IV)*wL(IV) ) )
               
               end if

            else     ! ( ustar .lt. 0.0_RP )
               if ( ( pstar .le. wR(IP) ) .and. ( wR(IU) + wR(IA) .le. 0.0_RP ) ) then
                  Fstar(IRHO)  = qR(IRHOU)
                  Fstar(IRHOU) = qR(IRHOU) * wR(IU) + wR(IP)
                  Fstar(IRHOV) = qR(IRHOV) * wR(IU)
                  Fstar(IRHOE) = wR(IU)*( cp * wR(IP) + 0.5_RP * (qR(IRHOU)*wR(IU) + qR(IRHOV)*wR(IV)) )

               elseif ( ( pstar .le. wR(IP) ) .and. (wR(IU) + wR(IA) .ge. 0.0_RP) .and. (ustar + wR(IA)*(pstar/wR(IP))**(0.5_RP / cp) .ge. 0.0_RP) ) then
                  rhostar = wR(IRHO) * ( pstar / wR(IP) ) ** ( 1.0_RP / gamma )

                  Fstar(IRHO) = rhostar * ustar 
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = Fstar(IRHO) * wR(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + wR(IV) * wR(IV)) ) 

               elseif ( ( pstar .le. wR(IP) ) .and. (wR(IU) + wR(IA) .gt. 0.0_RP) .and. (ustar + wR(IA) * (pstar/wR(IP))**(0.5_RP / cp) .lt. 0.0_RP ) ) then
                  uFan    = -2.0_RP * wR(IA) / ( gamma + 1.0_RP)  + (gamma-1.0_RP)/(gamma+1.0_RP) * wR(IU)
                  rhostar = wR(IRHO) * (-uFan / wR(IA)) ** (2.0_RP * cv)
                  pFan    = wR(IP) * (rhostar / wR(IRHO)) ** (gamma)

                  Fstar(IRHO) = rhostar * uFan
                  Fstar(IRHOU) = rhostar * uFan * uFan + pFan
                  Fstar(IRHOV) = rhostar * uFan * wR(IV)
                  Fstar(IRHOE) = uFan * ( cp * pFan + 0.5_RP * rhostar * ( uFan*uFan + wR(IV)*wR(IV) ) )

               elseif ( (pstar .gt. wR(IP)) .and. (wR(IU) + wR(IA)*( 0.5_RP*(gamma+1.0_RP)/gamma*pstar/wR(IP) + 0.5_RP/cp ) .le. 0.0_RP) ) then
                  Fstar(IRHO)  = qR(IRHOU)
                  Fstar(IRHOU) = qR(IRHOU) * wR(IU) + wR(IP)
                  Fstar(IRHOV) = qR(IRHOV) * wR(IU)
                  Fstar(IRHOE) = wR(IU)*( cp * wR(IP) + 0.5_RP * (qR(IRHOU)*wR(IU) + qR(IRHOV)*wR(IV)) )

               elseif ( (pstar .gt. wR(IP)) .and. (wR(IU) + wR(IA)*( (gamma+1.0_RP)/(2.0_RP*gamma)*pstar/wR(IP) + 0.5_RP/cp) .gt. 0.0_RP) ) then

                  rhostar = wR(IRHO) * ( (pstar/wR(IP) + (gamma-1.0_RP)/(gamma+1.0_RP))/(pstar/wR(IP)*(gamma-1.0_RP)/(gamma+1.0_RP) + 1.0_RP))
                  Fstar(IRHO) = rhostar * ustar
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = rhostar * ustar * wR(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + wR(IV)*wR(IV) ) )
               end if

            end if
            end associate

!        3/ Return to the 3D Space
!           ----------------------
            Ft = Fstar(IRHOU) * n(IX) - Fstar(IRHOV) * n(IY)
            Fstar(IRHOV) = Fstar(IRHOU) * n(IY) + Fstar(IRHOV) * n(IX)
            Fstar(IRHOU) = Ft

#ifdef _DIMENSIONLESS_TAU
            Fstar = Fstar * dimensionless % invSqrtGammaMach
#endif

      end function ExactRiemannSolver
         
      module pure function RoeFlux(qL, qR , n) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS), intent(in)     :: qL
         real(kind=RP), dimension(NCONS), intent(in)     :: qR
         real(kind=RP), dimension(NDIM) , intent(in)     :: n
         real(kind=RP), dimension(NCONS)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP)                 :: sqrtRhoL , invRhoL , rhoL , uL , vL , HL , TL , pL
         real(kind=RP)                 :: sqrtRhoR , invRhoR , rhoR , uR , vR , HR , TR , pR
         real(kind=RP)                 :: invrho , u , v , H , a
         real(kind=RP)                 :: dq(NCONS)
         real(kind=RP)                 :: lambda(NCONS)
         real(kind=RP)                 :: K(NCONS,NCONS)
         real(kind=RP)                 :: alpha(NCONS)
         real(kind=RP)                 :: Ft
         integer                       :: eq
         integer                       :: negativeWaves
         integer                       :: wave
        

!        0/ Gather variables
!           ----------------
            rhoL = qL(IRHO)
            invRhoL = 1.0_RP / rhoL
            sqrtRhoL = sqrt(rhoL)
            uL   = (  qL(IRHOU) * n(IX) + qL(IRHOV) * n(IY) ) * invRhoL
            vL   = ( -qL(IRHOU) * n(IY) + qL(IRHOV) * n(IX) ) * invRhoL
            pL   = thermodynamics % gm1 * (qL(IRHOE) - 0.5_RP * rhoL * (uL * uL + vL * vL) )
            TL   = dimensionless % gammaMach2 * pL * invRhoL
            HL   = dimensionless % cp * pL * invRhoL + 0.5_RP * ( uL*uL + vL*vL )

            rhoR = qR(IRHO)
            invRhoR = 1.0_RP / rhoR
            sqrtRhoR = sqrt(rhoR)
            uR   = (  qR(IRHOU) * n(IX) + qR(IRHOV) * n(IY) ) * invRhoR
            vR   = ( -qR(IRHOU) * n(IY) + qR(IRHOV) * n(IX) ) * invRhoR
            pR   = thermodynamics % gm1 * (qR(IRHOE) - 0.5_RP * rhoR * (uR * uR + vR * vR) )
            TR   = dimensionless % gammaMach2 * pR * invRhoR
            HR   = dimensionless % cp * pR * invRhoR + 0.5_RP * ( uR*uR + vR*vR )
!
!        1/ Compute Roe averages
!           --------------------
            invrho = 1.0_RP / (sqrtRhoL + sqrtRhoR)
            u      = (sqrtRhoL*uL + sqrtRhoR*uR) * invrho
            v      = (sqrtRhoL*vL + sqrtRhoR*vR) * invrho
            H      = (sqrtRhoL*HL + sqrtRhoR*HR) * invrho
            associate( gm1 => Thermodynamics % gm1 ) 
            a   = sqrt(gm1*(H - 0.5_RP*(u*u + v*v) ) )
            end associate

!
!        2/ Compute Roe matrix eigenvalues
!           ------------------------------
            lambda(1) = u - a
            lambda(2) = u
            lambda(3) = u
            lambda(4) = u + a 

            if ( lambda(1) .gt. 0.0_RP ) then
               negativeWaves = 0
            elseif ( lambda(2) .gt. 0.0_RP ) then
               negativeWaves = 1
            elseif ( lambda(4) .gt. 0.0_RP ) then
               negativeWaves = 3
            else
               negativeWaves = 4
            end if

!
!        3/ Compute the averaged right eigenvectors
!           ---------------------------------------
            K(1:NCONS , 1)  = reshape( (/ 1.0_RP , u-a    , v      , H-u*a                /)  ,  (/ NCONS /) )
            K(1:NCONS , 2)  = reshape( (/ 1.0_RP , u      , v      , 0.5_RP * (u*u + v*v) /)  ,  (/ NCONS /) )
            K(1:NCONS , 3)  = reshape( (/ 0.0_RP , 0.0_RP , 1.0_RP , v                    /)  ,  (/ NCONS /) )
            K(1:NCONS , 4)  = reshape( (/ 1.0_RP , u + a  , v      , H + u*a              /)  ,  (/ NCONS /) )
!
!        4/ Compute the wave strengths
!           --------------------------
            associate( gm1 => Thermodynamics % gm1 ) 
            dq(IRHO) = rhoR - rhoL
            dq(IRHOU) = rhoR*uR - rhoL*uL
            dq(IRHOV) = rhoR*vR - rhoL*vL
            dq(IRHOE) = qR(IRHOE) - qL(IRHOE)

            alpha(3) = dq(IRHOV) - v*dq(IRHO)
            alpha(2) = gm1 * (dq(IRHO) * (H - u*u) + u*dq(IRHOU) - dq(IRHOE) + (dq(IRHOV) - v*dq(IRHO))*v) / ( a*a )
            alpha(1) = (dq(IRHO) * (u + a) - dq(IRHOU) - a * alpha(2)) / (2.0_RP*a)
            alpha(4) = dq(IRHO) - (alpha(1) + alpha(2)) 
            end associate
!
!        5/ Compute the flux
!           ----------------
            Fstar = F_inviscidFlux(rhoL,uL,vL,pL,HL)
               
            do wave = 1 , negativeWaves
               Fstar = Fstar + alpha(wave) * lambda(wave) * K(1:NCONS , wave)
            end do
!
!        6/ Return to the 3D space
!           ----------------------
            Ft = Fstar(IRHOU) * n(IX) - Fstar(IRHOV) * n(IY)
            Fstar(IRHOV) = Fstar(IRHOU) * n(IY) + Fstar(IRHOV) * n(IX)
            Fstar(IRHOU) = Ft

#ifdef _DIMENSIONLESS_TAU
            Fstar = Fstar * dimensionless % invSqrtGammaMach
#endif
  
      end function RoeFlux

      module pure function HLLFlux(qL , qR , n) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS), intent(in)     :: qL
         real(kind=RP), dimension(NCONS), intent(in)     :: qR
         real(kind=RP), dimension(NDIM) , intent(in)     :: n
         real(kind=RP), dimension(NCONS)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP)                 :: rhoL , invRhoL , uL , vL , HL , aL , pL , rhoeL , sqrtRhoL
         real(kind=RP)                 :: rhoR , invRhoR , uR , vR , HR , aR , pR , rhoeR , sqrtRhoR
         real(kind=RP)                 :: invrho , u , v , H , a
         real(kind=RP)                 :: SL , SR , invdS
         real(kind=RP)                 :: sqrtS
         real(kind=RP)                 :: Ft
        

!        0/ Gather variables
!           ----------------
            rhoL     = qL(IRHO)
            invRhoL  = 1.0_RP / rhoL
            sqrtRhoL = sqrt(rhoL)
            uL       = (  qL(IRHOU) * n(IX) + qL(IRHOV) * n(IY) ) * invRhoL
            vL       = ( -qL(IRHOU) * n(IY) + qL(IRHOV) * n(IX) ) * invRhoL
            pL       = thermodynamics % gm1 * (qL(IRHOE) - 0.5_RP * rhoL * (uL * uL + vL * vL) )
            HL       = dimensionless % cp * pL * invRhoL + 0.5_RP * ( uL*uL + vL*vL )
            aL       = sqrt( thermodynamics % gamma * pL * invRhoL )
            rhoeL    = qL(IRHOE)

            rhoR     = qR(IRHO)
            invRhoR  = 1.0_RP / rhoR
            sqrtRhoR = sqrt(rhoR)
            uR       = (  qR(IRHOU) * n(IX) + qR(IRHOV) * n(IY) ) * invRhoR
            vR       = ( -qR(IRHOU) * n(IY) + qR(IRHOV) * n(IX) ) * invRhoR
            pR       = thermodynamics % gm1 * (qR(IRHOE) - 0.5_RP * rhoR * (uR * uR + vR * vR) )
            HR       = dimensionless % cp * pR * invRhoR + 0.5_RP * ( uR*uR + vR*vR )
            aR       = sqrt( thermodynamics % gamma * pR * invRhoR )
            rhoeR    = qR(IRHOE)

!        1/ Compute Roe averages
!           --------------------
            sqrtRhoL = sqrt(rhoL)
            sqrtRhoR = sqrt(rhoR)
            invrho = 1.0_RP / (sqrtRhoL + sqrtRhoR)
            u      = (sqrtRhoL*uL + sqrtRhoR*uR) * invrho
            v      = (sqrtRhoL*vL + sqrtRhoR*vR) * invrho
            H      = (sqrtRhoL*HL + sqrtRhoR*HR) * invrho
            associate( gm1 => Thermodynamics % gm1 ) 
            a   = sqrt(gm1*(H - 0.5_RP*(u*u + v*v) ) )
            end associate

!
!        2/ Compute wave speeds
!           -------------------
            SL = min(uL - aL , u - a)
            SR = max(uR + aR , u + a)

!
!        3/ Compute the fluxes depending on the speeds 
!           ------------------------------------------
            if ( SL .ge. 0.0_RP ) then
               Fstar = F_inviscidFlux(rhoL,uL,vL,pL,HL)

            elseif ( SR .le. 0.0_RP ) then
               Fstar = F_inviscidFlux(rhoR,uR,vR,pR,HR)

            elseif ( (SL .lt. 0.0_RP) .or. (SR .gt. 0.0_RP) ) then
               invdS = 1.0_RP / (SR - SL)
               Fstar(IRHO) = (SR*rhoL*uL - SL*rhoR*uR + SL*SR*(rhoR-rhoL)) * invdS
               Fstar(IRHOU) = (SR*(rhoL*uL*uL+pL) - SL*(rhoR*uR*uR+pR) + SL*SR*(rhoR*uR-rhoL*uL)) * invdS
               Fstar(IRHOV) = (SR*rhoL*uL*vL - SL*rhoR*uR*vR + SR*SL*(rhoR*vR-rhoL*vL)) * invdS
               Fstar(IRHOE) = (SR*uL*rhoL*HL - SL*uR*rhoR*HR + SR*SL*(rhoeR-rhoeL)) * invdS

            end if
!         
!        6/ Return to the 3D space
!           ----------------------
            Ft = Fstar(IRHOU) * n(IX) - Fstar(IRHOV) * n(IY)
            Fstar(IRHOV) = Fstar(IRHOU) * n(IY) + Fstar(IRHOV) * n(IX)
            Fstar(IRHOU) = Ft

#ifdef _DIMENSIONLESS_TAU
            Fstar = Fstar * dimensionless % invSqrtGammaMach
#endif

      end function HLLFlux
!
!     ****************************************************************************
!           Auxiliar module functions
!     ****************************************************************************
!
      pure subroutine ExactRiemann_ComputePStar(WL , WR , pstar , ustar )
         implicit none
         real(kind=RP), intent(in)           :: WL(NPRIM)         !  Left and right primitive variable
         real(kind=RP), intent(in)           :: WR(NPRIM)         !  sets.
         real(kind=RP), intent(out)           :: pstar
         real(kind=RP), intent(out)           :: ustar
!        -------------------------------------------------------------
         real(kind=RP), parameter :: TOL = 1.0e-6_RP
         integer, parameter       :: max_no_of_iterations = 50
         integer                  :: iter
         real(kind=RP)            :: FL , dFL , FR , dFR
         real(kind=RP)            :: F , dF
         real(kind=RP)            :: pold
         real(kind=RP)            :: cha
         
         

!        1/ Initial value for pstar: "Two-rarefaction approximation"
!           --------------------------------------------------------
            associate ( gm1 => Thermodynamics % gm1 , cp => Dimensionless % cp )
            pstar = ( (WL(IA) + WR(IA) - 0.5_RP * gm1 * (WR(IU) - WL(IU))) / ( WL(IA)/WL(IP)**(0.5_RP / cp) + WR(IA)/WR(IP)**(0.5_RP / cp) ) ) **(2.0_RP * cp)
            end associate

!
!        2/ Perform Newton iterations until the tolerance is reached
!           --------------------------------------------------------
            do iter = 1 , max_no_of_iterations

               call ExactRiemann_F(pstar , WL , FL , dFL)
               call ExactRiemann_F(pstar , WR , FR , dFR)
               F = FL + FR + WR(IU) - WL(IU)
               dF = dFL + dFR

               pold  = pstar
               pstar = pold - F / dF

!              Check for convergence
!              ---------------------
               cha = 0.5_RP * abs(pstar - pold) / (pstar + pold)

               if ( abs(cha) .lt. TOL ) then
!
!                 Compute ustar with pstar
!                 ------------------------
                  call ExactRiemann_F(pstar , WL , FL , dFL)
                  call ExactRiemann_F(pstar , WR , FR , dFR)
            
                  ustar = 0.5_RP * (WL(IU) + WR(IU) + FR - FL)
                  
                  return
               end if

            end do

      end subroutine ExactRiemann_ComputePStar

      pure subroutine ExactRiemann_F(p , W , F , dF )
         implicit none
         real(kind=RP), intent(IN)           :: p
         real(kind=RP), intent(IN)           :: W(NPRIM)
         real(kind=RP), intent(OUT)          :: F
         real(kind=RP), intent(OUT)          :: dF
!        ------------------------------------------
         real(kind=RP)           :: A , B
         real(kind=RP)           :: sqrtADivpPlusB

         associate ( gamma => Thermodynamics % gamma )
         if ( p .gt. W(IP) ) then
            A = 2.0_RP / (W(IRHO) * ( gamma + 1.0_RP ) )
            B = Thermodynamics % gm1 * W(IP) / (gamma + 1.0_RP)
            sqrtADivpPlusB = sqrt(A / (p+B) )
            F = (p - W(IP)) * sqrtADivpPlusB
            dF = sqrtADivpPlusB * (1.0_RP - 0.5_RP * (p-W(IP)) / (B + p) )
         else
            F = 2.0_RP * Dimensionless % cv * W(IA) * ( ( p/W(IP) ) ** (0.5_RP /  Dimensionless % cp)  - 1  )
            dF = 1.0_RP / (W(IRHO) * W(IA)) * (p / W(IP)) ** ( -0.5_RP * (gamma + 1.0_RP)/gamma) 
         end if
         end associate

      end subroutine ExactRiemann_F

end submodule RiemannSolvers
