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
!                 StdVolumeGreen:   \int \vec{F} \cdot \nabla v dx
!                 StdFace       :   \int \vec{F} \cdot \vec{n} v ds
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
         integer, pointer                 :: N
   
         N => e % spA % N
   
         volInt =    MatrixMultiplyInIndex_F(F(:,:,:,IX) , e % spA % hatD , N+1 , N+1 , NCONS , IX)          &
                   + MatrixMultiplyInIndex_F(F(:,:,:,IY) , e % spA % hatD , N+1 , N+1 , NCONS , IY)
   
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
