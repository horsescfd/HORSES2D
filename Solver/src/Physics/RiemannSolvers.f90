submodule (PhysicsNS)  RiemannSolvers
   use SMConstants
   contains


!
!     ****************************************************
!        Riemann solvers
!     ****************************************************
!
      module function ExactRiemannSolver(qL3D , qR3D , wL3D , wR3D , T , Tinv ) result (Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS)     :: qL3D
         real(kind=RP), dimension(NCONS)     :: qR3D
         real(kind=RP), dimension(NPRIM)     :: wL3D
         real(kind=RP), dimension(NPRIM)     :: wR3D
         real(kind=RP), dimension(NCONS,NCONS) :: T
         real(kind=RP), dimension(NCONS,NCONS) :: Tinv
         real(kind=RP), dimension(NCONS)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP), dimension(NCONS)   :: qL , qR
         real(kind=RP), dimension(NPRIM) :: WL , WR
         real(kind=RP)                   :: pstar , ustar
         real(kind=RP)                   :: rhostar , uFan , pFan

!        0/ Gather variables
!           ----------------
            qL(IRHO) = qL3D(IRHO)
            qL(IRHOU) = qL3D(IRHOU) * T(IRHOU,IRHOU) + qL3D(IRHOV) * T(IRHOU,IRHOV)
            qL(IRHOV) = qL3D(IRHOU) * T(IRHOV,IRHOU) + qL3D(IRHOV) * T(IRHOV,IRHOV)
            qL(IRHOE) = qL3D(IRHOE) * T(IRHOE,IRHOE)

            qR(IRHO) = qR3D(IRHO)
            qR(IRHOU) = qR3D(IRHOU) * T(IRHOU,IRHOU) + qR3D(IRHOV) * T(IRHOU,IRHOV)
            qR(IRHOV) = qR3D(IRHOU) * T(IRHOV,IRHOU) + qR3D(IRHOV) * T(IRHOV,IRHOV)
            qR(IRHOE) = qR3D(IRHOE) * T(IRHOE,IRHOE)


            associate( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )
            WL(IRHO) = qL(IRHO)
            WL(IU)   = qL(IRHOU) / qL(IRHO)
            WL(IV)   = qL(IRHOV) / qL(IRHO)
            WL(IP)   = gm1 * (qL(IRHOE) - 0.5_RP * (qL(IRHOU) * WL(IU) + qL(IRHOV) * WL(IV) ) )
            WL(IA)   = sqrt(gamma * WL(IP) / WL(IRHO))

            WR(IRHO) = qR(IRHO)
            WR(IU)   = qR(IRHOU) / qR(IRHO)
            WR(IV)   = qR(IRHOV) / qR(IRHO)
            WR(IP)   = gm1 * (qR(IRHOE) - 0.5_RP * (qR(IRHOU) * WR(IU) + qR(IRHOV) * WR(IV) ) )
            WR(IA)   = sqrt(gamma * WR(IP) / WR(IRHO))
            end associate

!        1/ Compute the star region
!           -----------------------
            call ExactRiemann_ComputePStar (WL , WR , pstar , ustar )

!        2/ Check to which region belongs the solution
!           ------------------------------------------
            associate ( cv => Dimensionless % cv , cp => Dimensionless % cp , gamma => Thermodynamics % gamma ) 
            if ( ustar .ge. 0.0_RP ) then

               if ( ( pstar .le. WL(IP) ) .and. ( WL(IU) .ge. WL(IA) ) ) then
                  Fstar(IRHO)  = qL(IRHOU)
                  Fstar(IRHOU) = qL(IRHOU) * WL(IU) + WL(IP)
                  Fstar(IRHOV) = qL(IRHOV) * WL(IU)
                  Fstar(IRHOE) = WL(IU)*( cp * WL(IP) + 0.5_RP * (qL(IRHOU)*WL(IU) + qL(IRHOV)*WL(IV)) )

               elseif ( ( pstar .le. WL(IP) ) .and. ( WL(IU) .lt. WL(IA)) .and. (ustar .lt. WL(IA) * (pstar/WL(IP))**(0.5_RP / cp ) ) ) then
                  rhostar = WL(IRHO) * ( pstar / WL(IP) ) ** ( 1.0_RP / gamma )

                  Fstar(IRHO) = rhostar * ustar 
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = Fstar(IRHO) * WL(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + WL(IV) * WL(IV)) ) 

               elseif ( ( pstar .le. WL(IP) ) .and. ( WL(IU) .lt. WL(IA) ) .and. ( ustar .ge. WL(IA) * (pstar / WL(IP))**(0.5_RP / cp))) then
                  
                  uFan    = 2.0_RP * WL(IA) / ( gamma + 1.0_RP)  + (gamma-1.0_RP)/(gamma+1.0_RP) * WL(IU)
                  rhostar = WL(IRHO) * (uFan / WL(IA)) ** (2.0_RP * cv)
                  pFan    = WL(IP) * (rhostar / WL(IRHO)) ** (gamma)

                  Fstar(IRHO) = rhostar * uFan
                  Fstar(IRHOU) = rhostar * uFan * uFan + pFan
                  Fstar(IRHOV) = rhostar * uFan * WL(IV)
                  Fstar(IRHOE) = uFan * ( cp * pFan + 0.5_RP * rhostar * ( uFan*uFan + WL(IV)*WL(IV) ) )
 
               elseif ( ( pstar .gt. WL(IP) ) .and. ( WL(IU) .ge. WL(IA) * sqrt( (gamma+1.0_RP)/(2.0_RP * gamma)*pstar/WL(IP) + 0.5_RP/cp) ) ) then
                  Fstar(IRHO)  = qL(IRHOU)
                  Fstar(IRHOU) = qL(IRHOU) * WL(IU) + WL(IP)
                  Fstar(IRHOV) = qL(IRHOV) * WL(IU)
                  Fstar(IRHOE) = WL(IU)*( cp * WL(IP) + 0.5_RP * (qL(IRHOU)*WL(IU) + qL(IRHOV)*WL(IV)) )

               elseif ( ( pstar .gt. WL(IP) ) .and. ( WL(IU) .lt.  WL(IA) * sqrt( (gamma+1.0_RP)/(2.0_RP * gamma)*pstar/WL(IP) + 0.5_RP/cp) ) ) then
                  rhostar = WL(IRHO) * ( (pstar/WL(IP) + (gamma-1.0_RP)/(gamma+1.0_RP))/(pstar/WL(IP)*(gamma-1.0_RP)/(gamma+1.0_RP) + 1.0_RP))
                  Fstar(IRHO) = rhostar * ustar
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = rhostar * ustar * WL(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + WL(IV)*WL(IV) ) )
               
               end if

            else     ! ( ustar .lt. 0.0_RP )
               if ( ( pstar .le. WR(IP) ) .and. ( WR(IU) + WR(IA) .le. 0.0_RP ) ) then
                  Fstar(IRHO)  = qR(IRHOU)
                  Fstar(IRHOU) = qR(IRHOU) * WR(IU) + WR(IP)
                  Fstar(IRHOV) = qR(IRHOV) * WR(IU)
                  Fstar(IRHOE) = WR(IU)*( cp * WR(IP) + 0.5_RP * (qR(IRHOU)*WR(IU) + qR(IRHOV)*WR(IV)) )

               elseif ( ( pstar .le. WR(IP) ) .and. (WR(IU) + WR(IA) .ge. 0.0_RP) .and. (ustar + WR(IA)*(pstar/WR(IP))**(0.5_RP / cp) .ge. 0.0_RP) ) then
                  rhostar = WR(IRHO) * ( pstar / WR(IP) ) ** ( 1.0_RP / gamma )

                  Fstar(IRHO) = rhostar * ustar 
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = Fstar(IRHO) * WR(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + WR(IV) * WR(IV)) ) 

               elseif ( ( pstar .le. WR(IP) ) .and. (WR(IU) + WR(IA) .gt. 0.0_RP) .and. (ustar + WR(IA) * (pstar/WR(IP))**(0.5_RP / cp) .lt. 0.0_RP ) ) then
                  uFan    = -2.0_RP * WR(IA) / ( gamma + 1.0_RP)  + (gamma-1.0_RP)/(gamma+1.0_RP) * WR(IU)
                  rhostar = WR(IRHO) * (-uFan / WR(IA)) ** (2.0_RP * cv)
                  pFan    = WR(IP) * (rhostar / WR(IRHO)) ** (gamma)

                  Fstar(IRHO) = rhostar * uFan
                  Fstar(IRHOU) = rhostar * uFan * uFan + pFan
                  Fstar(IRHOV) = rhostar * uFan * WR(IV)
                  Fstar(IRHOE) = uFan * ( cp * pFan + 0.5_RP * rhostar * ( uFan*uFan + WR(IV)*WR(IV) ) )

               elseif ( (pstar .gt. WR(IP)) .and. (WR(IU) + WR(IA)*( 0.5_RP*(gamma+1.0_RP)/gamma*pstar/WR(IP) + 0.5_RP/cp ) .le. 0.0_RP) ) then
                  Fstar(IRHO)  = qR(IRHOU)
                  Fstar(IRHOU) = qR(IRHOU) * WR(IU) + WR(IP)
                  Fstar(IRHOV) = qR(IRHOV) * WR(IU)
                  Fstar(IRHOE) = WR(IU)*( cp * WR(IP) + 0.5_RP * (qR(IRHOU)*WR(IU) + qR(IRHOV)*WR(IV)) )

               elseif ( (pstar .gt. WR(IP)) .and. (WR(IU) + WR(IA)*( (gamma+1.0_RP)/(2.0_RP*gamma)*pstar/WR(IP) + 0.5_RP/cp) .gt. 0.0_RP) ) then

                  rhostar = WR(IRHO) * ( (pstar/WR(IP) + (gamma-1.0_RP)/(gamma+1.0_RP))/(pstar/WR(IP)*(gamma-1.0_RP)/(gamma+1.0_RP) + 1.0_RP))
                  Fstar(IRHO) = rhostar * ustar
                  Fstar(IRHOU) = rhostar * ustar * ustar + pstar
                  Fstar(IRHOV) = rhostar * ustar * WR(IV)
                  Fstar(IRHOE) = ustar * ( cp * pstar + 0.5_RP * rhostar * (ustar * ustar + WR(IV)*WR(IV) ) )
               end if

            end if
            end associate

!        3/ Return to the 3D Space
!           ----------------------
            Fstar = MatrixTimesVector_F( A = Tinv , X = Fstar , Nout = NCONS)

!        4/ Scale it with the Mach number
!           -----------------------------
            associate( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
            Fstar = Fstar * dimensionless % invSqrtGammaMach
            end associate


      end function ExactRiemannSolver
         
      module function RoeFlux(qL3D , qR3D , wL3D , wR3D , T , Tinv) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS)     :: qL3D
         real(kind=RP), dimension(NCONS)     :: qR3D
         real(kind=RP), dimension(NPRIM)     :: wL3D
         real(kind=RP), dimension(NPRIM)     :: wR3D
         real(kind=RP), dimension(NCONS,NCONS) :: T
         real(kind=RP), dimension(NCONS,NCONS) :: Tinv
         real(kind=RP), dimension(NCONS)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP), dimension(NCONS) :: qL , qR
         real(kind=RP)                 :: sqrtRhoL , rhoL , uL , vL , HL , TL , pL
         real(kind=RP)                 :: sqrtRhoR , rhoR , uR , vR , HR , TR , pR
         real(kind=RP)                 :: invrho , u , v , H , a
         real(kind=RP)                 :: dq(NCONS)
         real(kind=RP)                 :: lambda(NCONS)
         real(kind=RP)                 :: K(NCONS,NCONS)
         real(kind=RP)                 :: alpha(NCONS)
         integer                       :: eq
         integer                       :: negativeWaves
         integer                       :: wave
        

!        0/ Gather variables
!           ----------------
            rhoL = wL3D(IRHO)
            sqrtRhoL = sqrt(rhoL)
            uL   = wL3D(IU) * T(IRHOU,IRHOU) + wL3D(IV) * T(IRHOU,IRHOV)
            vL   = wL3D(IU) * T(IRHOV,IRHOU) + wL3D(IV) * T(IRHOV,IRHOV)
            pL   = wL3D(IP) * T(IRHOE,IRHOE)
            TL   = wL3D(IT) * T(IRHOE,IRHOE)
            HL   = dimensionless % cp * TL + 0.5_RP * ( uL*uL + vL*vL )

            rhoR = wR3D(IRHO)
            sqrtRhoR = sqrt(rhoR)
            uR   = wR3D(IU) * T(IRHOU,IRHOU) + wR3D(IV) * T(IRHOU,IRHOV)
            vR   = wR3D(IU) * T(IRHOV,IRHOU) + wR3D(IV) * T(IRHOV,IRHOV)
            pR   = wR3D(IP) * T(IRHOE,IRHOE)
            TR   = wR3D(IT) * T(IRHOE,IRHOE)
            HR   = dimensionless % cp * TR + 0.5_RP * ( uR*uR + vR*vR )

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
            dq(IRHOE) = (qR3D(IRHOE) - qL3D(IRHOE)) * T(IRHOE,IRHOE)

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
  
!        6/ Return to the 3D Space
!           ----------------------
            Fstar = MatrixTimesVector_F( A = Tinv , X = Fstar , Nout = NCONS)

!        7/ Scale it with the Mach number
!           -----------------------------
            associate( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
            Fstar = Fstar * dimensionless % invSqrtGammaMach
            end associate

!
      end function RoeFlux

      module function HLLFlux(qL3D , qR3D , wL3D , wR3D , T , Tinv) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS)     :: qL3D
         real(kind=RP), dimension(NCONS)     :: qR3D
         real(kind=RP), dimension(NPRIM)     :: wL3D
         real(kind=RP), dimension(NPRIM)     :: wR3D
         real(kind=RP), dimension(NCONS,NCONS) :: T
         real(kind=RP), dimension(NCONS,NCONS) :: Tinv
         real(kind=RP), dimension(NCONS)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP)                 :: rhoL , uL , vL , HL , aL , pL , rhoeL , TL , sqrtRhoL
         real(kind=RP)                 :: rhoR , uR , vR , HR , aR , pR , rhoeR , TR , sqrtRhoR
         real(kind=RP)                 :: invrho , u , v , H , a
         real(kind=RP)                 :: SL , SR , invdS
         real(kind=RP)                 :: sqrtS
        

!        0/ Gather variables
!           ----------------
            associate( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )
            sqrtS = sqrt(T(IRHOE,IRHOE))

            rhoL = wL3D(IRHO)
            uL   = wL3D(IU) * T(IRHOU,IRHOU) + wL3D(IV) * T(IRHOU,IRHOV)
            vL   = wL3D(IU) * T(IRHOV,IRHOU) + wL3D(IV) * T(IRHOV,IRHOV)
            pL   = wL3D(IP) * T(IRHOE,IRHOE)
            TL   = wL3D(IT) * T(IRHOE,IRHOE)
            rhoeL = dimensionless % cv * pL + 0.5_RP * rhoL * (uL * uL + vL * vL)
            HL   = dimensionless % cp * TL + 0.5_RP * ( uL * uL + vL * vL )
            aL   = wL3D(IA) * sqrtS

            rhoR = wR3D(IRHO)
            uR   = wR3D(IU) * T(IRHOU,IRHOU) + wR3D(IV) * T(IRHOU,IRHOV)
            vR   = wR3D(IU) * T(IRHOV,IRHOU) + wR3D(IV) * T(IRHOV,IRHOV)
            pR   = wR3D(IP) * T(IRHOE,IRHOE)
            TR   = wR3D(IT) * T(IRHOE,IRHOE)
            rhoeR = dimensionless % cv * pR + 0.5_RP * rhoR * (uR * uR + vR * vR)
            HR   = dimensionless % cp * TR + 0.5_RP * ( uR * uR + vR * vR )
            aR   = wR3D(IA) * sqrtS

            end associate

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
         
!        4/ Return to the 3D Space
!           ----------------------
            Fstar = MatrixTimesVector_F( A = Tinv , X = Fstar , Nout = NCONS)

!        5/ Scale it with the Mach number
!           -----------------------------
            associate( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
            Fstar = Fstar * dimensionless % invSqrtGammaMach
            end associate
!        
      end function HLLFlux

      module function AUSMFlux(qL3D , qR3D , wL3D , wR3D , T , Tinv) result(Fstar)
         use MatrixOperations
         implicit none
         real(kind=RP), dimension(NCONS)     :: qL3D
         real(kind=RP), dimension(NCONS)     :: qR3D
         real(kind=RP), dimension(NPRIM)     :: wL3D
         real(kind=RP), dimension(NPRIM)     :: wR3D
         real(kind=RP), dimension(NCONS,NCONS) :: T
         real(kind=RP), dimension(NCONS,NCONS) :: Tinv
         real(kind=RP), dimension(NCONS)     :: Fstar
!        ---------------------------------------------------------------
         real(kind=RP), dimension(NCONS) :: qL , qR
         real(kind=RP)                 :: rhoL , uL , vL , pL , aL , ML
         real(kind=RP)                 :: rhoR , uR , vR , pR , aR , MR
         real(kind=RP)                 :: MplusL , MminusR , pplusL , pminusR 
         real(kind=RP)                 :: M , p 
        

!        0/ Gather variables
!           ----------------
            qL(IRHO) = qL3D(IRHO)
            qL(IRHOU) = qL3D(IRHOU) * T(IRHOU,IRHOU) + qL3D(IRHOV) * T(IRHOU,IRHOV)
            qL(IRHOV) = qL3D(IRHOU) * T(IRHOV,IRHOU) + qL3D(IRHOV) * T(IRHOV,IRHOV)
            qL(IRHOE) = qL3D(IRHOE) * T(IRHOE,IRHOE)

            qR(IRHO) = qR3D(IRHO)
            qR(IRHOU) = qR3D(IRHOU) * T(IRHOU,IRHOU) + qR3D(IRHOV) * T(IRHOU,IRHOV)
            qR(IRHOV) = qR3D(IRHOU) * T(IRHOV,IRHOU) + qR3D(IRHOV) * T(IRHOV,IRHOV)
            qR(IRHOE) = qR3D(IRHOE) * T(IRHOE,IRHOE)


            associate( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 )
            rhoL = sqrt(qL(IRHO))
            uL   = qL(IRHOU) / qL(IRHO)
            vL   = qL(IRHOV) / qL(IRHO)
            pL   = gm1 * (qL(IRHOE) - 0.5_RP  * ( qL(IRHOU)*uL + qL(IRHOV)*vL ) )  
            aL   = sqrt(gamma * pL / rhoL )
            ML   = uL / aL

            rhoR = sqrt(qR(IRHO))
            uR   = qR(IRHOU) / qR(IRHO)
            vR   = qR(IRHOV) / qR(IRHO)
            pR   = gm1 * (qR(IRHOE) - 0.5_RP  * ( qR(IRHOU)*uR + qR(IRHOV)*vR ) ) 
            aR   = sqrt(gamma * pR / rhoR )
            MR   = uR / aR
            end associate

!        1/ Compute splitted Mach numbers and pressures
!           -------------------------------------------
            if ( abs(ML) .le. 1.0_RP ) then
               MplusL = 0.25_RP * (ML + 1.0_RP)*(ML + 1.0_RP)
               pplusL = 0.5_RP * pL * (1.0_RP + ML)
            else
               MplusL = 0.5_RP * (ML + abs(ML) )
               pplusL = 0.5_RP * pL * ( sign(1.0_RP,ML) + 1.0_RP )
            end if

            if ( abs(MR) .le. 1.0_RP) then
               MminusR = -0.25_RP * ( MR - 1.0_RP) * ( MR - 1.0_RP )
               pminusR = 0.5_RP * p * (1.0_RP - MR)
            else
               MminusR = 0.5_RP * (MR - abs(MR))
               pminusR = 0.5_RP * p * (1.0_RP - sign(1.0_RP , MR))
            end if

!        2/ Compute intercell pressure and Mach number
!           ------------------------------------------
            p = pplusL + pminusR
            M = MplusL + MminusR

!        3/ Compute intercell flux
!           ---------------------------------
            if ( M .ge. 0.0_RP ) then
               Fstar(IRHO) = M * rhoL * aL
               Fstar(IRHOU) = M * rhoL * aL * uL + p
               Fstar(IRHOV) = M * rhoL * aL * vL
               Fstar(IRHOE) = M * aL * ( Dimensionless % cp * pL + 0.5_RP * rhoL * (uL * uL + vL * vL))
            else
               Fstar(IRHO) = M * rhoR * aR
               Fstar(IRHOU) = M * rhoR * aR * uR + p
               Fstar(IRHOV) = M * rhoR * aR * vR
               Fstar(IRHOE) = M * aR * ( Dimensionless % cp * pR + 0.5_RP * rhoR * (uR * uR + vR * vR))
            end if
               
!        4/ Return to the 3D Space
!           ----------------------
            Fstar = MatrixTimesVector_F( A = Tinv , X = Fstar  , Nout = NCONS )

!        5/ Scale it with the Mach number
!           -----------------------------
            associate( gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
            Fstar = Fstar * dimensionless % invSqrtGammaMach
            end associate
!
      end function AUSMFlux
!
!     ****************************************************************************
!           Auxiliar module functions
!     ****************************************************************************
!
      subroutine ExactRiemann_ComputePStar(WL , WR , pstar , ustar )
         implicit none
         real(kind=RP)           :: WL(NPRIM)         !  Left and right primitive variable
         real(kind=RP)           :: WR(NPRIM)         !  sets.
         real(kind=RP)           :: pstar
         real(kind=RP)           :: ustar
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

            print*, "Warning: Exact Riemann solver exceeded the maximum number of iterations."

      end subroutine ExactRiemann_ComputePStar

      subroutine ExactRiemann_F(p , W , F , dF )
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
