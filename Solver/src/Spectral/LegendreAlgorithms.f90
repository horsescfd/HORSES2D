!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
!    Copyright (C) 2009  David A. Kopriva 
!    Copyright (C) 2011  Gonzalo Rubio Calzado (g.rubio@upm.es)
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
!////////////////////////////////////////////////////////////////////////
!
!      LegendreAlgorithms.f90
!      Created: 2009-12-08 13:47:42 -0500 
!      By: David Kopriva
!
!
!      Modifications: 25-08-2011
!      By: Gonzalo Rubio      
!
!      Contains:
!            ALGORITHM 22: SUBROUTINE LegendrePolyAndDerivative( N, x, L_N, LPrime_N )
!            ALGORITHM 23: SUBROUTINE GaussLegendreNodesAndWeights( N, x, w )
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE LegendreAlgorithms
      USE SMConstants
      IMPLICIT NONE

#include "Defines.h"
      
      PUBLIC  :: GaussLegendreNodesAndWeights
      PRIVATE :: LegendrePolyAndDerivative
!
!     ========
      CONTAINS
!     ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LegendrePolyAndDerivative( N, x, L_N, LPrime_N )
!
!     Compute the Legendre Polynomial of degree k and its derivative

      IMPLICIT NONE 
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER      , INTENT(IN) :: N
      REAL(KIND=RP), INTENT(IN) :: x
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP), INTENT(OUT) :: L_N, LPrime_N
!
!     ----------------
!     Local Variables:
!     ----------------
!
      INTEGER       :: k
      REAL(KIND=RP) :: L_NM1, L_NM2, LPrime_NM2, LPrime_NM1
      
      IF( N == 0 )     THEN
         L_N      = 1.0_rp
         LPrime_N = 0.0_rp
      ELSE IF ( N == 1 )     THEN
         L_N      = x
         LPrime_N = 1.0_rp
      ELSE
         L_NM2 = 1.0_rp
         LPrime_NM2 = 0.0_rp
         L_NM1 = x
         LPrime_NM1 = 1.0_rp
         DO k = 2, N
            L_N        = ((2.0_RP*k-1.0_RP)*x*L_NM1 - (k-1.0_RP)*L_NM2)/k
            LPrime_N   = LPrime_NM2 + (2.0_RP*k-1.0_RP)*L_NM1
            L_NM2      = L_NM1
            L_NM1      = L_N
            LPrime_NM2 = LPrime_NM1
            LPrime_NM1 = LPrime_N
         END DO
      END IF
      
      END SUBROUTINE LegendrePolyAndDerivative
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE GaussLegendreNodesAndWeights( N, x, w )
      
      IMPLICIT NONE 
!
!     Compute the Gauss-legendre quadrature nodes and
!     weights
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER, INTENT(IN) :: N
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: x, w
!
!     ----------------
!     Local Variables:
!     ----------------
!
      REAL(KIND=RP) :: xj, L_NP1, LPrime_NP1, delta, tolerance
      INTEGER       :: j, k, NDiv2
!
!     ----------
!     Constants:
!     ----------
!
      INTEGER, PARAMETER       :: noNewtonIterations = 100
      REAL(KIND=RP), PARAMETER :: toleranceFactor    = 4.0_RP
      
      tolerance = toleranceFactor*EPSILON(L_NP1)
      IF( N == 0 )     THEN
         x(0) = 0.0_RP
         w(0) = 2.0_RP
         RETURN
      ELSE IF( N == 1 )     THEN
         x(0) = -SQRT(1.0_RP/3.0_RP)
         w(0) =  1.0_RP
         x(1) = -x(0)
         w(1) =  w(0)
      ELSE
!
!        ----------------------------------
!        Iterate on half the interior nodes
!        ----------------------------------
!
         NDiv2 = (N+1)/2
         DO j = 0, NDiv2-1
            xj = -COS( (2*j+1)*PI/(2*N+2) )
            DO k = 0, noNewtonIterations
               CALL LegendrePolyAndDerivative( N+1, xj, L_NP1, LPrime_NP1 )
               delta = -L_NP1/LPrime_NP1
               xj = xj + delta
               IF( ABS(delta) <=  tolerance*ABS(xj) )     EXIT
            END DO
            CALL LegendrePolyAndDerivative( N+1, xj, L_NP1, LPrime_NP1 )
            x(j)   = xj
            w(j)   = 2.0_RP/( (1.0_RP - xj*xj)*LPrime_NP1*LPrime_NP1 )
            x(N-j) = -xj
            w(N-j) = w(j)
         END DO
      END IF
!
!     ---------------------------
!     Fill in middle if necessary
!     ---------------------------
!
      IF( MOD(N,2) == 0 )     THEN
         CALL LegendrePolyAndDerivative( N+1, 0.0_RP, L_NP1, LPrime_NP1 )
         x(N/2) = 0.0_RP
         w(N/2) = 2.0_RP/(LPrime_NP1*LPrime_NP1)
      END IF
      
      END SUBROUTINE GaussLegendreNodesAndWeights


    SUBROUTINE LegendreGaussLobattoNodesAndWeights(N,x,w)
        ! Algorithm to compute N-th order LegendreGaussLobatto quadrature nodes and weights
        implicit none
        integer,intent(IN)::N
        double precision,intent(OUT)::x(0:),w(0:)
        integer:: j=0,ITERS=10,k=0
        double precision:: pi=acos(-1d0),delta=0d0, TOL=1d-10,q=0d0,dq=0d0,LN=0d0
        if (N.EQ.1) then
            x(0)=-1d0
            x(1)=1d0
            w(0)=1d0
            w(1)=1d0
        else
            x(0)=-1d0
            w(0)=2d0/(1d0*N*(1d0*N+1d0))
            x(N)=1d0
            w(N)=w(0)
            do j=1,(((N+1)/2)-1)
                x(j)=-cos((j+0.25d0)*pi/N-3d0/(8d0*N*pi*(j+0.25d0)))

                do k=0,ITERS
                    call qAndLEvaluation(N,x(j),q,dq,LN)
                    delta=-q/dq
                    x(j)=x(j)+delta
                    if (abs(delta) .LT. TOL*abs(x(j))) then
                        exit
                    endif

                end do
                call qAndLEvaluation(N,x(j),q,dq,LN)
                x(N-j)=-x(j)
                w(j)=2d0/(1d0*N*(N+1d0)*LN**2d0)
                w(N-j)=w(j)
            end do
        endif

        if (mod(N,2) .EQ. 0) then
            call qAndLEvaluation(N,0d0,q,dq,LN)
            x(N/2)=0d0
            w(N/2)=2d0/(1d0*N*(N+1d0)*LN**2d0)
        endif
    END SUBROUTINE LegendreGaussLobattoNodesAndWeights

    subroutine qAndLEvaluation(N,x,q,dq,LN)
        implicit none
        ! Combined algorithm to compute q=L_{N+1}-L_{N-1}, q', and L_N
        integer,intent(IN)::N
        real(kind=RP),intent(IN)::x
        real(kind=RP),intent(OUT)::q,dq,LN
        integer::k=2
        real(kind=RP)         ::dLN=0.0_RP ,LNm2=0.0_RP ,LNm1=0.0_RP ,dLNm2=0.0_RP ,dLNm1=0.0_RP, LNp1=0.0_RP ,dLNp1=0.0_RP
        LNm2=1.0_RP
        LNm1=x
        dLNm2=0.0_RP
        dLNm1=1.0_RP

        do k=2,N
            LN=(2.0_RP *k-1.0_RP)*x*LNm1/k-(1.0_RP*k-1.0_RP)*LNm2/k
            dLN=dLNm2+(2.0_RP*k-1.0_RP)*LNm1
            LNm2=LNm1
            LNm1=LN
            dLNm2=dLNm1
            dLNm1=dLN
        end do

        k=N+1
        LNp1=(2.0_RP*k-1.0_RP)*x*LN/k-(1.0_RP*k-1.0_RP)*LNm2/k
        dLNp1=dLNm2+(2.0_RP*k-1.0_RP)*LN
        q=LNp1-LNm2
        dq=dLNp1-dLNm2

    end subroutine qAndLEvaluation
 
!     ==========      
      END MODULE LegendreAlgorithms
!     ==========      
