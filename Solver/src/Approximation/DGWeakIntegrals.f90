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
      function Scalar_StdVolumeGreen( e , F) result ( volInt )
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
         class(QuadElement_t)             :: e
         real(kind=RP), intent(in)        :: F(0:e % spA % N , 0:e % spA % N, 1:NCONS , 1:NDIM)
         real(kind=RP)                    :: volInt(0:e % spA % N, 0:e % spA % N , 1:NCONS)
!        ----------------------------------------------------------------------------------------------
         real(kind=RP)                    :: contravariant_F( 0 : e % spA % N , 0 : e % spA % N , 1:NCONS )
         real(kind=RP)                    :: contravariant_G( 0 : e % spA % N , 0 : e % spA % N , 1:NCONS )
         integer, pointer                 :: N
         integer                          :: eq
   
         N => e % spA % N
 
         do eq = 1 , NCONS
!           
!           F flux (contravariant)
!           ----------------------
            contravariant_F(0:N,0:N,eq) = F(0:N,0:N,eq,IX) * e % Ja(0:N,0:N,1,1) + F(0:N,0:N,eq,IY) * e % Ja(0:N,0:N,2,1)

         end do
   
         do eq = 1 , NCONS
!           
!           G flux (contravariant)
!           ----------------------
            contravariant_G(0:N,0:N,eq) = F(0:N,0:N,eq,IX) * e % Ja(0:N,0:N,1,2) + F(0:N,0:N,eq,IY) * e % Ja(0:N,0:N,2,2)

         end do

  

         volInt =    MatrixMultiplyInIndex_F(contravariant_F , e % spA % hatD , N+1 , N+1 , NCONS , IX)          &
                   + MatrixMultiplyInIndex_F(contravariant_G , e % spA % hatD , N+1 , N+1 , NCONS , IY)
   
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
         real(kind=RP), intent(in) :: F  (0 : ed % storage(loc) % spA % N , 1:NCONS)
         real(kind=RP)             :: faceInt(0 : ed % storage(loc) % spA % N , 0 : ed % storage(loc) % spA % N , 1:NCONS)
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
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = F(iXi,eq) * e % spA % lbw(iEta,LEFT)
               end do ; end do ; end do

            case (ERIGHT)
!
!              -----------------------------------------------------------            
!>             Right edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=+1
!              -----------------------------------------------------------            
!
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = F(iEta,eq) * e % spA % lbw(iXi,RIGHT)
               end do ; end do ; end do
                  
            case (ETOP)
!
!              -----------------------------------------------------------            
!>             Top edge: Variables are defined for Xi-coordinate, 
!                       while the test function is written in Eta=+1
!              -----------------------------------------------------------            
!
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = F(iXi,eq) * e % spA % lbw(iEta,RIGHT)
               end do ; end do ; end do
   
            case (ELEFT)
!
!              -----------------------------------------------------------            
!>             Left edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=-1
!              -----------------------------------------------------------            
!
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = F(iEta,eq) * e % spA % lbw(iXi,LEFT)
               end do ; end do ; end do

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
         real(kind=RP), intent(in)     :: F(0:ed % storage(loc) % spA % N,1:NCONS,1:NDIM)
         real(kind=RP)                 :: faceInt(0:ed % storage(loc) % spA % N , 0:ed % storage(loc) % spA % N , 1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Ja_xi(0:ed % storage(loc) % spA % N , 1:NDIM)
         real(kind=RP)     :: Ja_eta(0:ed % storage(loc) % spA % N , 1:NDIM)
         real(kind=RP)     :: jac(0:ed % storage(loc) % spA % N)
         real(kind=RP)     :: F_contravariant( 0 : ed % storage(loc) % spA % N , 1:NCONS )
         real(kind=RP)     :: G_contravariant( 0 : ed % storage(loc) % spA % N , 1:NCONS )
         real(kind=RP)     :: dF             ( 0 : ed % storage(loc) % spA % N , 1:NCONS )
         real(kind=RP)     :: dG             ( 0 : ed % storage(loc) % spA % N , 1:NCONS )
         integer     :: iXi , iEta , eq , iDim

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
               do iDim = 1 , NDIM
                  Ja_xi (:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IX) , e % spA % lb(:,LEFT) , e % spA % N + 1 )
                  Ja_eta(:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IY) , e % spA % lb(:,LEFT) , e % spA % N + 1 )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,LEFT) , e % spA % N + 1 )

               do eq = 1 , NCONS
                  F_contravariant(:,eq) = ( F(:,eq,IX) * Ja_xi (:,IX) + F(:,eq,IY) * Ja_xi (:,IY) ) / jac
                  G_contravariant(:,eq) = ( F(:,eq,IX) * Ja_eta(:,IX) + F(:,eq,IY) * Ja_eta(:,IY) ) / jac
               end do
!
!              dF(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dF = Mat_x_Mat_F( e % spA % hatD , F_contravariant , e % spA % N + 1 , NCONS ,  trA = .true. ) 
               
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = dF(iXi,eq) * e % spA % lbw(iEta,LEFT) +  G_contravariant(iXi,eq) * e % spA % dlbw(iEta,LEFT)
               end do ;            end do ;                    end do

            case (ERIGHT)
!
!              -----------------------------------------------------------            
!>             Right edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=+1
!              -----------------------------------------------------------            
!
               do iDim = 1 , NDIM
                  Ja_xi (:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IX) , e % spA % lb(:,RIGHT) , e % spA % N + 1 , trA = .true. )
                  Ja_eta(:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IY) , e % spA % lb(:,RIGHT) , e % spA % N + 1 , trA = .true. )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,RIGHT) , e % spA % N + 1 , trA = .true.)

               do eq = 1 , NCONS
                  F_contravariant(:,eq) = ( F(:,eq,IX) * Ja_xi (:,IX) + F(:,eq,IY) * Ja_xi (:,IY) ) / jac
                  G_contravariant(:,eq) = ( F(:,eq,IX) * Ja_eta(:,IX) + F(:,eq,IY) * Ja_eta(:,IY) ) / jac
               end do
!
!              dG(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dG = Mat_x_Mat_F( e % spA % hatD , G_contravariant , e % spA % N + 1 , NCONS ,  trA = .true. ) 
               
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = F_contravariant(iEta,eq) * e % spA % dlbw(iXi,RIGHT) + dG(iEta,eq) * e % spA % lbw(iXi,RIGHT)
               end do ;            end do ;                    end do

            case (ETOP)
!
!              -----------------------------------------------------------            
!>             Top edge: Variables are defined for Xi-coordinate, 
!                       while the test function is written in Eta=1
!              -----------------------------------------------------------            
!
               do iDim = 1 , NDIM
                  Ja_xi (:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IX) , e % spA % lb(:,RIGHT) , e % spA % N + 1 )
                  Ja_eta(:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IY) , e % spA % lb(:,RIGHT) , e % spA % N + 1 )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,RIGHT) , e % spA % N + 1 )

               do eq = 1 , NCONS
                  F_contravariant(:,eq) = ( F(:,eq,IX) * Ja_xi (:,IX) + F(:,eq,IY) * Ja_xi (:,IY) ) / jac
                  G_contravariant(:,eq) = ( F(:,eq,IX) * Ja_eta(:,IX) + F(:,eq,IY) * Ja_eta(:,IY) ) / jac
               end do
!
!              dF(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dF = Mat_x_Mat_F( e % spA % hatD , F_contravariant , e % spA % N + 1 , NCONS ,  trA = .true. ) 
               
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = dF(iXi,eq) * e % spA % lbw(iEta,RIGHT) +  G_contravariant(iXi,eq) * e % spA % dlbw(iEta,RIGHT)
               end do ;            end do ;                    end do

            case (ELEFT)
!
!              -----------------------------------------------------------            
!>             Left edge: Variables are defined for Eta-coordinate, 
!                       while the test function is written in Xi=+0
!              -----------------------------------------------------------            
!
               do iDim = 1 , NDIM
                  Ja_xi (:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IX) , e % spA % lb(:,LEFT) , e % spA % N + 1 , trA = .true. )
                  Ja_eta(:,iDim) = MatrixTimesVector_F( e % Ja(:,:,iDim,IY) , e % spA % lb(:,LEFT) , e % spA % N + 1 , trA = .true. )
               end do
               
               jac = MatrixTimesVector_F( e % jac , e % spA % lb(:,LEFT) , e % spA % N + 1 , trA = .true.)

               do eq = 1 , NCONS
                  F_contravariant(:,eq) = ( F(:,eq,IX) * Ja_xi (:,IX) + F(:,eq,IY) * Ja_xi (:,IY) ) / jac
                  G_contravariant(:,eq) = ( F(:,eq,IX) * Ja_eta(:,IX) + F(:,eq,IY) * Ja_eta(:,IY) ) / jac
               end do
!
!              dG(i,eq) = hatD(m,i) * F(m,eq)
!              ------------------------------
               dG = Mat_x_Mat_F( e % spA % hatD , G_contravariant , e % spA % N + 1 , NCONS ,  trA = .true. ) 
               
               do eq = 1 , NCONS ; do iEta = 0 , e % spA % N ; do iXi = 0 , e % spA % N
                  faceInt(iXi,iEta,eq) = F_contravariant(iEta,eq) * e % spA % dlbw(iXi,LEFT) + dG(iEta,eq) * e % spA % lbw(iXi,LEFT)
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
      function Vector_StdVolumeGreen( e , u ) result (volInt)
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
         class(QuadElement_t), target     :: e
         real(kind=RP), intent(in)        :: u(0:e % spA % N , 0:e % spA % N, 1:NCONS)
         real(kind=RP)                    :: volInt(0:e % spA % N, 0:e % spA % N , 1:NDIM , 1:NCONS)
!        ----------------------------------------------------------------------------------------------
         integer, pointer                 :: N
         real(kind=RP), pointer           :: JaXi(:,:,:) , JaEta(:,:,:)
         real(kind=RP)                    :: contravariantU(0:e % spA % N, 0:e % spA % N, 1:NDIM)
         integer                          :: eq
   
         N => e % spA % N
         JaXi (0:,0:,1:)    => e % Ja(0:,0:,1:,IX)
         JaEta(0:,0:,1:)    => e % Ja(0:,0:,1:,IY)
!
!>       Perform the xi-component of the test function divergence
!        --------------------------------------------------------
         do eq = 1 , NCONS
!
!>          contravariantU^xi = U * Ja^xi
!           -----------------------------
            contravariantU(:,:,IX) = u(:,:,eq) * JaXi(:,:,IX) 
            contravariantU(:,:,IY) = u(:,:,eq) * JaXi(:,:,IY) 
!
!>          Compute the volume integral
!           ---------------------------           
            volInt(:,:,:,eq) = MatrixMultiplyInIndex_F( contravariantU , e % spA % hatD , N+1 , N+1 , NDIM , IX ) 
         end do
!
!>       Perform the eta-component of the test function divergence
!        ---------------------------------------------------------
         do eq = 1 , NCONS
!
!>          contravariantU^eta = U * Ja^eta
!           -------------------------------
            contravariantU(:,:,IX) = u(:,:,eq) * JaEta(:,:,IX) 
            contravariantU(:,:,IY) = u(:,:,eq) * JaEta(:,:,IY) 
!
!>          Compute the volume integral
!           ---------------------------           
            volInt(:,:,:,eq) =   volInt(:,:,:,eq)  &
                               + MatrixMultiplyInIndex_F( contravariantU , e % spA % hatD , N+1 , N+1 , NDIM , IY ) 
         end do
!
      end function Vector_StdVolumeGreen

      pure function Vector_StdFace( ed , loc , u ) result ( faceInt )
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
         implicit none
         class(Edge_t), intent(in) :: ed
         integer,       intent(in) :: loc
         real(kind=RP), intent(in) :: u (0 : ed % storage(loc) % spA % N , 1:NCONS)
         real(kind=RP)             :: faceInt(0 : ed % storage(loc) % spA % N , 0 : ed % storage(loc) % spA % N , 1:NDIM , 1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: iXi , iEta , eq
         real(kind=RP)                 :: loc_sign
!
!        -------------------------------------------------------------
!>       If the edge is interior, it is precise to change the normal
!>       sign so that it points towards the element outside
!        -------------------------------------------------------------
!
         if ( loc .eq. LEFT ) then
!
!           Already points towards the outside
!           ----------------------------------
            loc_sign = 1.0_RP
      
         elseif ( loc .eq. RIGHT ) then
!
!           Change its sign
!           ---------------
            loc_sign = -1.0_RP

         end if

         associate ( e => ed % quads(loc) % e ) 

         select case ( ed % edgeLocation(loc) ) 

            case (EBOTTOM)
            
               select type ( ed ) 
                  type is (Edge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX:IY,eq) = loc_sign * u(iXi,eq) * e % spA % lbw(iEta,LEFT) * ed % n(IX:IY,0) * ed % dS(0)
                     end do ; end do ; end do

                  type is (CurvedEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX:IY,eq) = loc_sign * u(iXi,eq) * e % spA % lbw(iEta,LEFT) * ed % n(IX:IY,iXi) * ed % dS(iXi)
                     end do ; end do ; end do
                  
                  type is (StraightBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iXi,eq) * e % spA % lbw(iEta,LEFT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = u(iXi,eq) * e % spA % lbw(iEta,LEFT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do
   
                  type is (CurvedBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iXi,eq) * e % spA % lbw(iEta,LEFT) * ed % n(IX,iXi) * ed % dS(iXi)
                        faceInt(iXi,iEta,IY,eq) = u(iXi,eq) * e % spA % lbw(iEta,LEFT) * ed % n(IY,iXi) * ed % dS(iXi)
                     end do ; end do ; end do

               end select
      
            case (ERIGHT)

               select type ( ed ) 
                  type is (Edge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do
                  
                  type is (CurvedEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IX,iEta) * ed % dS(iEta)
                        faceInt(iXi,iEta,IY,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IY,iEta) * ed % dS(iEta)
                     end do ; end do ; end do

                  type is (StraightBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do
   
                  type is (CurvedBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IX,iEta) * ed % dS(iEta)
                        faceInt(iXi,iEta,IY,eq) = u(iEta,eq) * e % spA % lbw(iXi,RIGHT) * ed % n(IY,iEta) * ed % dS(iEta)
                     end do ; end do ; end do

               end select
                  
            case (ETOP)

               select type ( ed ) 
                  type is (Edge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = loc_sign * u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = loc_sign * u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do

                  type is (CurvedEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = loc_sign * u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IX,iXi) * ed % dS(iXi)
                        faceInt(iXi,iEta,IY,eq) = loc_sign * u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IY,iXi) * ed % dS(iXi)
                     end do ; end do ; end do
                  
                  type is (StraightBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do
   
                  type is (CurvedBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IX,iXi) * ed % dS(iXi)
                        faceInt(iXi,iEta,IY,eq) = u(iXi,eq) * e % spA % lbw(iEta,RIGHT) * ed % n(IY,iXi) * ed % dS(iXi)
                     end do ; end do ; end do

               end select
   
            case (ELEFT)

               select type ( ed ) 

                  type is (Edge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do

                  type is (CurvedEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IX,iEta) * ed % dS(iEta)
                        faceInt(iXi,iEta,IY,eq) = loc_sign * u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IY,iEta) * ed % dS(iEta)
                     end do ; end do ; end do
                  
                  type is (StraightBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IX,0) * ed % dS(0)
                        faceInt(iXi,iEta,IY,eq) = u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IY,0) * ed % dS(0)
                     end do ; end do ; end do
   
                  type is (CurvedBdryEdge_t)
                     do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                        faceInt(iXi,iEta,IX,eq) = u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IX,iEta) * ed % dS(iEta)
                        faceInt(iXi,iEta,IY,eq) = u(iEta,eq) * e % spA % lbw(iXi,LEFT) * ed % n(IY,iEta) * ed % dS(iEta)
                     end do ; end do ; end do

                  

               end select
            case default 
         end select

         end associate

      end function Vector_StdFace

end module DGWeakIntegrals
