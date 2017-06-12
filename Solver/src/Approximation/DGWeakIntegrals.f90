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
!        File:    DGWeakIntegrals.f90
!
!        This file contains the Weak Integrals involved in the DG Formulations. There are two
!     main classes:
!
!           >> Scalar weak integrals: Those whose test function is scalar
!
!                 StdVolumeGreen  :   \int \vec{F} \cdot \nabla v dx
!                 StdFace         :   \int \vec{F} \cdot \vec{n} v ds
!                 StdGradientFace :   \int \vec{F} \cdot \nabla v ds
!
!              These weak integrals are computed for 0:N - 0:N - 1:NCONS test functions
!                    (i,j,eq) < l_i^eq l_j^eq 
!
!                  ---------------
!              And THEY ARE SCALED by the Mass Matrix: 1/(w_i w_j)
!                  ---------------
!
!           >> Vector weak integrals: Those whose test function is vectorial
!
!                 StdVolumeGreen:   \int u \nabla \cdot \vec{tau} dx
!                 StdFace       :   \int u \vec{tau} \cdot \vec{n} ds
!
!              These weak integrals are computed for 0:N - 0:N - 1:NDIM - 1:NCONS test functions
!                    (i,j,dim,eq) < l_i^{dim,eq} l_j^{dim,eq} 
!
!                  ---------------
!              And THEY ARE SCALED by the Mass Matrix: 1/(w_i w_j)
!                  ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
module DGWeakIntegrals
   use SMConstants
   use Physics
   use QuadElementClass
!
#include "Defines.h"
!
!  *******
   private
   public   ScalarWeakIntegrals_t , VectorWeakIntegrals_t
   public   ScalarWeakIntegrals   , VectorWeakIntegrals
!  *******
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!                 --------------------
!                 | TYPE DEFINITIONS |
!                 --------------------
!
!  ----------------------------------------------------
!        SCALAR WEAK INTEGRALS TYPE
!  ----------------------------------------------------
!
   type     ScalarWeakIntegrals_t
      contains
         procedure, nopass      :: StdVolumeGreen  => Scalar_StdVolumeGreen
         procedure, nopass      :: StdFace         => Scalar_StdFace
         procedure, nopass      :: StdGradientFace => Scalar_StdGradientFace
   end type ScalarWeakIntegrals_t
!
!  ----------------------------------------------------
!        VECTORIAL WEAK INTEGRALS TYPE
!  ----------------------------------------------------
!
   type     VectorWeakIntegrals_t
      contains
         procedure, nopass      :: StdVolumeGreen  => Vector_StdVolumeGreen
         procedure, nopass      :: StdFace         => Vector_StdFace
   end type VectorWeakIntegrals_t
!
!
   type(ScalarWeakIntegrals_t)      :: ScalarWeakIntegrals
   type(VectorWeakIntegrals_t)      :: VectorWeakIntegrals
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with scalar test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure function Scalar_StdVolumeGreen( e , F) result ( volInt )
!
!     ***************************************************************************
!           Computes a weak volume integrals with:
!              -> Scalar test function
!              -> Vector fluxes
!              -> Green form
!
!                 volInt = \int_{\Omega} \vec{F} \cdot \nabla v dx
!
!           Its discrete version:
!        
!            volInt(i,j,eq) = \frac{1}{w_i w_j} \int_{element} \vec{F} \cdot \nabla \phi_{ij}^{eq} dx   
!
!        This is performed by means of the volume derivative matrix:
!
!           volInt = MatMultiplyInIndex(F,hatD,IX) + MatMultiplyInIndex(G,hatD,IY)
!
!        Since F and G stand for the CONTRAVARIANT fluxes.
!                                    -------------
!     ***************************************************************************
!
         use MatrixOperations
         implicit none
         class(QuadElement_t), intent(in) :: e
         real(kind=RP), intent(in)        :: F(1:NCONS , 0:e % spA % N , 0:e % spA % N , 1:NDIM)
         real(kind=RP)                    :: volInt(1:NCONS , 0:e % spA % N, 0:e % spA % N )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: contravariant_F( 1:NCONS , 0 : e % spA % N , 0 : e % spA % N )
         real(kind=RP)                    :: contravariant_G( 1:NCONS , 0 : e % spA % N , 0 : e % spA % N )
         integer                          :: eq , i , j , l , N , dimID

         N = e % spA % N 
!
!        Initialization
!        -------------- 
         contravariant_F = 0.0_RP
         contravariant_G = 0.0_RP
         volInt = 0.0_RP
!
!        Contravariant fluxes
!        --------------------
         do dimID = 1 , NDIM  ; do j = 0 , e % spA % N ; do i = 0 , e % spA % N 
!           
!           F flux (contravariant)
!           ----------------------
            contravariant_F(:,i,j) = contravariant_F(:,i,j) + F(:,i,j,dimID) * e % Ja(i,j,dimID,IX) 

         end do               ; end do                 ; end do

         do dimID = 1 , NDIM  ; do j = 0 , e % spA % N ; do i = 0 , e % spA % N 
!           
!           G flux (contravariant)
!           ----------------------
            contravariant_G(:,i,j) = contravariant_G(:,i,j) + F(:,i,j,dimID) * e % Ja(i,j,dimID,IY) 

         end do               ; end do                 ; end do
!
!        Compute the integral
!        --------------------   
         do l = 0 , e % spA % N
            do j = 0 , e % spA % N ; do i = 0 , e % spA % N  
               volInt(:,i,j) = volInt(:,i,j) + e % spA % hatD(i,l) * contravariant_F(:,l,j) + e % spA % hatD(j,l) * contravariant_G(:,i,l)
            end do                 ; end do
         end do
     
      end function Scalar_StdVolumeGreen

      pure function Scalar_StdFace( ed , loc , F ) result ( faceInt )
!
!     *****************************************************************************
!           This computes the following weak integral with scalar test function:
!           
!              faceInt = \frac{1}{w_i w_j} \int_0^1 F \phi_{ij} ds
!
!           The integral is computed at the face, which is located at the
!           "loc" (e.g. LEFT/RIGHT) element.
!
!           First, the position referred to the element is computed, to apply
!           the appropriate test function \phi. Then, it is integrated
!           along the face.
!
!           Important: To this subroutine, F = \vec{F} \cdot \vec{n} dS, enters such 
!           that a scalar integral is computed for F in the computational domain!
!        
!     *****************************************************************************
!
         implicit none
         class(Edge_t), intent(in) :: ed
         integer,       intent(in) :: loc
         real(kind=RP), intent(in) :: F  (1:NCONS , 0 : ed % storage(loc) % spA % N )
         real(kind=RP)             :: faceInt(1:NCONS , 0 : ed % storage(loc) % spA % N , 0 : ed % storage(loc) % spA % N )
!        ------------------------------------------------------------------------------------------------------------------------
         integer     :: iXi , iEta , eq

         associate ( e => ed % quads(loc) % e )
!
!        ------------------------------------------------------------
!>       The location of the edge referred to the element is gathered
!        ------------------------------------------------------------
!
         select case ( ed % edgeLocation(loc) ) 

            case (EBOTTOM)
!
!              -----------------------------------------------------------            
!>             Bottom edge: Variables are defined for Xi-coordinate, 
!                       while the test function is written in Eta=-1
!              -----------------------------------------------------------            
!
               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = F(:,iXi) * e % spA % lbw(iEta,LEFT)
               end do ; end do

            case (ERIGHT)
!
!              -----------------------------------------------------------            
!>             Right edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=+1
!              -----------------------------------------------------------            
!
               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = F(:,iEta) * e % spA % lbw(iXi,RIGHT)
               end do ; end do
                  
            case (ETOP)
!
!              -----------------------------------------------------------            
!>             Top edge: Variables are defined for Xi-coordinate, 
!                       while the test function is written in Eta=+1
!              -----------------------------------------------------------            
!
               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = F(:,iXi) * e % spA % lbw(iEta,RIGHT)
               end do ; end do
   
            case (ELEFT)
!
!              -----------------------------------------------------------            
!>             Left edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=-1
!              -----------------------------------------------------------            
!
               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = F(:,iEta) * e % spA % lbw(iXi,LEFT)
               end do ; end do

         end select
      
         end associate

      end function Scalar_StdFace

      pure function Scalar_StdGradientFace( ed , loc , F ) result ( faceInt )
!
!     *****************************************************************************
!           This computes the following weak integral with scalar test function:
!           
!              faceInt = \frac{1}{w_i w_j} \int_0^1 \vec{F}\cdot\nabla  \phi_{ij} dxi
!
!           The integral is computed at the face, which is located at the
!           "loc" (e.g. LEFT/RIGHT) element.
!
!           First, the position referred to the element is computed, to apply
!           the appropriate test function \phi. Then, it is integrated
!           along the face.
!
!           Important: To this subroutine, \vec{F} <= \vec{F} dS, enters such 
!           that a scalar integral is computed for \vec{F} in the computational domain!
!        
!     *****************************************************************************
!
         use MatrixOperations
         implicit none
         class(Edge_t), intent(in)     :: ed
         integer,       intent(in)     :: loc
         real(kind=RP), intent(in)     :: F(1:NCONS,0:ed % storage(loc) % spA % N,1:NDIM)
         real(kind=RP)                 :: faceInt(1:NCONS,0:ed % storage(loc) % spA % N , 0:ed % storage(loc) % spA % N)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Ja_xi(0:ed % storage(loc) % spA % N , 1:NDIM)
         real(kind=RP)     :: Ja_eta(0:ed % storage(loc) % spA % N , 1:NDIM)
         real(kind=RP)     :: jac(0:ed % storage(loc) % spA % N)
         real(kind=RP)     :: F_contravariant( 1:NCONS , 0 : ed % storage(loc) % spA % N )
         real(kind=RP)     :: G_contravariant( 1:NCONS , 0 : ed % storage(loc) % spA % N )
         real(kind=RP)     :: dF             ( 1:NCONS , 0 : ed % storage(loc) % spA % N )
         real(kind=RP)     :: dG             ( 1:NCONS , 0 : ed % storage(loc) % spA % N )
         integer     :: iXi , iEta , eq , dimID , i , l

         associate ( e => ed % quads(loc) % e )
!
!        ------------------------------------------------------------
!>       The location of the edge referred to the element is gathered
!        ------------------------------------------------------------
!
         select case ( ed % edgeLocation(loc) ) 

            case (EBOTTOM)
!
!              -----------------------------------------------------------            
!>             Bottom edge: Variables are defined for the Xi-coordinate, 
!                       while the test function is written in Eta=-1
!              -----------------------------------------------------------            
!
               do dimID = 1 , NDIM         ! TODO   : improve this
                  Ja_xi (:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IX) , e % spA % lb(:,LEFT) , e % spA % N + 1 )
                  Ja_eta(:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IY) , e % spA % lb(:,LEFT) , e % spA % N + 1 )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,LEFT) , e % spA % N + 1 )

               F_contravariant = 0.0_RP
               G_contravariant = 0.0_RP
               
               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  F_contravariant(:,i) = F_contravariant(:,i) + F(:,i,dimID) * Ja_xi (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  F_contravariant(:,i) = F_contravariant(:,i) / jac(i)
               end do             

               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  G_contravariant(:,i) = G_contravariant(:,i) + F(:,i,dimID) * Ja_eta (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  G_contravariant(:,i) = G_contravariant(:,i) / jac(i)
               end do             
!
!              dF(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dF = 0.0_RP
               do l = 0 , ed % storage(loc) % spA % N    ; do i = 0 , ed % storage(loc) % spA % N 
                  dF(:,i)  = dF(:,i) + e % spA % hatD(i,l) * F_contravariant(:,l)
               end do                                    ; end do
               
               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = dF(:,iXi) * e % spA % lbw(iEta,LEFT) + G_contravariant(:,iXi) * e % spA % dlbw(iEta,LEFT)
               end do ;                    end do

            case (ERIGHT)
!
!              -----------------------------------------------------------            
!>             Right edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=+1
!              -----------------------------------------------------------            
!
               do dimID = 1 , NDIM
                  Ja_xi (:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IX) , e % spA % lb(:,RIGHT) , e % spA % N + 1 , trA = .true. )
                  Ja_eta(:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IY) , e % spA % lb(:,RIGHT) , e % spA % N + 1 , trA = .true. )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,RIGHT) , e % spA % N + 1 , trA = .true.)

               F_contravariant = 0.0_RP
               G_contravariant = 0.0_RP
               
               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  F_contravariant(:,i) = F_contravariant(:,i) + F(:,i,dimID) * Ja_xi (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  F_contravariant(:,i) = F_contravariant(:,i) / jac(i)
               end do

               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  G_contravariant(:,i) = G_contravariant(:,i) + F(:,i,dimID) * Ja_eta (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  G_contravariant(:,i) = G_contravariant(:,i) / jac(i)
               end do
!
!              dG(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dG = 0.0_RP
               do l = 0 , ed % storage(loc) % spA % N    ; do i = 0 , ed % storage(loc) % spA % N 
                  dG(:,i)  = dG(:,i) + e % spA % hatD(i,l) * G_contravariant(:,l)
               end do                                    ; end do

               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = F_contravariant(:,iEta) * e % spA % dlbw(iXi,RIGHT) + dG(:,iEta) * e % spA % lbw(iXi,RIGHT)
               end do ;                    end do

            case (ETOP)
!
!              -----------------------------------------------------------            
!>             Top edge: Variables are defined for Xi-coordinate, 
!                       while the test function is written in Eta=1
!              -----------------------------------------------------------            
!
               do dimID = 1 , NDIM
                  Ja_xi (:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IX) , e % spA % lb(:,RIGHT) , e % spA % N + 1 )
                  Ja_eta(:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IY) , e % spA % lb(:,RIGHT) , e % spA % N + 1 )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,RIGHT) , e % spA % N + 1 )

               F_contravariant = 0.0_RP
               G_contravariant = 0.0_RP
               
               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  F_contravariant(:,i) = F_contravariant(:,i) + F(:,i,dimID) * Ja_xi (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  F_contravariant(:,i) = F_contravariant(:,i) / jac(i)
               end do

               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  G_contravariant(:,i) = G_contravariant(:,i) + F(:,i,dimID) * Ja_eta (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  G_contravariant(:,i) = G_contravariant(:,i) / jac(i)
               end do
!
!              dF(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dF = 0.0_RP
               do l = 0 , ed % storage(loc) % spA % N    ; do i = 0 , ed % storage(loc) % spA % N 
                  dF(:,i)  = dF(:,i) + e % spA % hatD(i,l) * F_contravariant(:,l)
               end do                                    ; end do
               
               do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = dF(:,iXi) * e % spA % lbw(iEta,RIGHT) +  G_contravariant(:,iXi) * e % spA % dlbw(iEta,RIGHT)
               end do ;                    end do  

            case (ELEFT)
!
!              -----------------------------------------------------------            
!>             Left edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=+0
!              -----------------------------------------------------------            
!
               do dimID = 1 , NDIM
                  Ja_xi (:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IX) , e % spA % lb(:,LEFT) , e % spA % N + 1 , trA = .true. )
                  Ja_eta(:,dimID) = MatrixTimesVector_F( e % Ja(:,:,dimID,IY) , e % spA % lb(:,LEFT) , e % spA % N + 1 , trA = .true. )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,LEFT) , e % spA % N + 1 , trA = .true.)

               F_contravariant = 0.0_RP
               G_contravariant = 0.0_RP
               
               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  F_contravariant(:,i) = F_contravariant(:,i) + F(:,i,dimID) * Ja_xi (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  F_contravariant(:,i) = F_contravariant(:,i) / jac(i)
               end do

               do dimID = 1 , NDIM  ; do i = 0 , ed % storage(loc) % spA % N 
                  G_contravariant(:,i) = G_contravariant(:,i) + F(:,i,dimID) * Ja_eta (i,dimID) 
               end do               ; end do

               do i = 0 , ed % storage(loc) % spA % N
                  G_contravariant(:,i) = G_contravariant(:,i) / jac(i)
               end do
!
!              dG(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dG = 0.0_RP
               do l = 0 , ed % storage(loc) % spA % N    ; do i = 0 , ed % storage(loc) % spA % N 
                  dG(:,i)  = dG(:,i) + e % spA % hatD(i,l) * G_contravariant(:,l)
               end do               ;                      end do

               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta) = F_contravariant(:,iEta) * e % spA % dlbw(iXi,LEFT) + dG(:,iEta) * e % spA % lbw(iXi,LEFT)
               end do ;            end do ;                    end do

         end select
      
         end associate

      end function Scalar_StdGradientFace
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with vector test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure function Vector_StdVolumeGreen( e , u ) result (volInt)
!
!     ***************************************************************************
!           Computes a weak volume integrals with:
!              -> Vector test function (NDIM)
!              -> Scalar argument (u)
!              -> Green form
!
!                 volInt = \int_{\Omega} u (\nabla \cdot \vec{v}) dx
!
!           Its discrete version:
!        
!            volInt(i,j,dim,eq) = \frac{1}{w_i w_j} \int_{element} u  \partial_{dim} \phi_{ij}^{eq} dx   
!
!        This is performed by means of the volume derivative matrix:
!
!           volInt(i,j,dim,eq) =    MatMultiplyInIndex(u^{eq}*Ja^{xi }_dim,hatD,IX) 
!                                 + MatMultiplyInIndex(u^{eq}*Ja^{eta}_dim,hatD,IY)
!     ******************************************************************************************************************
!
         use MatrixOperations
         implicit none
         class(QuadElement_t), intent (in) :: e
         real(kind=RP),        intent (in) :: u(1:NCONS , 0:e % spA % N , 0:e % spA % N)
         real(kind=RP)                     :: volInt(1:NCONS , 0:e % spA % N, 0:e % spA % N , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: eq , i , j , l , dimID

         volInt = 0.0_RP
!
!        Derivative in the local Xi coordinate
!        -------------------------------------
         do dimID = 1 , NDIM ; 
            do j = 0 , e % spA % N ;  do l = 0 , e % spA % N ; do i = 0 , e % spA % N 
               volInt(:,i,j,dimID) = volInt(:,i,j,dimID) + u(:,l,j) * e % Ja(l,j,dimID,IX) * e % spA % hatD(i,l)
            end do                 ;  end do                 ; end do
         end do
!
!        Derivative in the local Eta coordinate
!        --------------------------------------
         do dimID = 1 , NDIM ; 
            do l = 0 , e % spA % N ;  do j = 0 , e % spA % N ; do i = 0 , e % spA % N 
               volInt(:,i,j,dimID) = volInt(:,i,j,dimID) + u(:,i,l) * e % Ja(i,l,dimID,IY) * e % spA % hatD(j,l)
            end do                 ;  end do                 ; end do
         end do

      end function Vector_StdVolumeGreen

      pure function Vector_StdFace( ed , loc , Fu ) result ( faceInt )
!
!     ***************************************************************************************
!           This computes the following weak integral with vectorial test function:
!           
!           faceInt = \frac{1}{w_i w_j} \int_0^1 u^{eq} \vec{\tau}_{ij} \cdot \hat{n} ds
!
!           The integral is computed at the face, which is located at the
!           "loc" (e.g. LEFT/RIGHT) element.
!
!           First, the position referred to the element is computed, to apply
!           the appropriate test function \tau. Then, it is integrated
!           along the face.
!
!           The obtained result is:
!
!              faceInt(i,j,dim,eq) = \frac{1}{w_i w_j} \int_0^1 u^{eq} \vec{\tau}_{ij} n_{dim} ds

!        
!     ***************************************************************************************
!
         use NodesAndWeights_Class
         implicit none
         class(Edge_t), intent(in) :: ed
         integer,       intent(in) :: loc
         real(kind=RP), intent(in) :: Fu (1:NCONS , 0 : ed % storage(loc) % spA % N , 1:NDIM)
         real(kind=RP)             :: faceInt(1:NCONS , 0 : ed % storage(loc) % spA % N , 0 : ed % storage(loc) % spA % N , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: iXi , iEta , eq , i , dimID

         associate ( e => ed % quads(loc) % e ) 

         select case ( ed % edgeLocation(loc) ) 

            case (EBOTTOM)
            
               do dimID = 1 , NDIM ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta,dimID) = Fu(:,iXi,dimID) * e % spA % lbw(iEta,LEFT)
               end do ;              end do ;                   end do
      
            case (ERIGHT)

               do dimID = 1 , NDIM ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta,dimID) = Fu(:,iEta,dimID) * e % spA % lbw(iXi,RIGHT) 
               end do ; end do ; end do
                  
            case (ETOP)

               do dimID = 1 , NDIM ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta,dimID) = Fu(:,iXi,dimID) * e % spA % lbw(iEta,RIGHT)
               end do ; end do ; end do
   
            case (ELEFT)

               do dimID = 1 , NDIM ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(:,iXi,iEta,dimID) = Fu(:,iEta,dimID) * e % spA % lbw(iXi,LEFT)
               end do ; end do ; end do

            case default 
         end select

         end associate

      end function Vector_StdFace

end module DGWeakIntegrals
