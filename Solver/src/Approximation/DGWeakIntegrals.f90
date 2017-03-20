module DGWeakIntegrals
   use SMConstants
   use Physics
   use QuadElementClass
   use QuadMeshDefinitions


   private
   public   WeakIntegrals_t

   type WeakIntegrals_t
      contains
         procedure, nopass      :: StdVolumeFormI
         procedure, nopass      :: StdFace
   end type WeakIntegrals_t


   contains

      function StdVolumeFormI( e , F) result ( val )
!
!     ***************************************************************************
!           Computes a weak volume integral, which consists in:
!        
!            Val = \frac{1}{w_i w_j} \int_{element} F \nabla \phi_{ij} dx   
!
!        This is performed by means of the volume derivative matrix:
!
!           Val = MatMultiplyInIndex(F,hatD,IX) + MatMultiplyInIndex(G,hatD,IY)
!     ***************************************************************************
!
         use MatrixOperations
         implicit none
         class(QuadElement_t)             :: e
         real(kind=RP), intent(in)        :: F(0:e % spA % N , 0:e % spA % N, 1:NCONS , 1:NDIM)
         real(kind=RP)                    :: val(0:e % spA % N, 0:e % spA % N , 1:NCONS)
!        ----------------------------------------------------------------------------------------------
         integer, pointer                 :: N
   
         N => e % spA % N
   
         val =    MatrixMultiplyInIndex_F(F(:,:,:,IX) , e % spA % hatD , N+1 , N+1 , NCONS , IX)          &
                + MatrixMultiplyInIndex_F(F(:,:,:,IY) , e % spA % hatD , N+1 , N+1 , NCONS , IY)
   
      end function StdVolumeFormI

      function StdFace( ed , loc , F ) result ( val )
!
!     *******************************************************************
!           This computes the following weak integral:
!           
!              Val = \frac{1}{w_i w_j} \int_0^1 F \phi_{ij} ds
!        
!     *******************************************************************
!
         implicit none
         class(Edge_t)                            :: ed
         integer,       intent(in)                :: loc
         real(kind=RP), intent(in)                 :: F  (0 : ed % storage(loc) % spA % N , 1:NCONS)
         real(kind=RP)                 :: val(0 : ed % storage(loc) % spA % N , 0 : ed % storage(loc) % spA % N , 1:NCONS)
!        ------------------------------------------------------------------------------------------------------------------------
         class(QuadElement_t), pointer :: e
         integer                       :: iXi , iEta , eq

         e => ed % quads(loc) % e

         select case ( ed % edgeLocation(loc) ) 

            case (EBOTTOM)
            
               do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                  val(iXi,iEta,eq) = F(iXi,eq) * e % spA % lbw(iEta,LEFT)
               end do ; end do ; end do

            case (ERIGHT)

               do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                  val(iXi,iEta,eq) = F(iEta,eq) * e % spA % lbw(iXi,RIGHT)
               end do ; end do ; end do
                  
            case (ETOP)

               do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                  val(iXi,iEta,eq) = F(iXi,eq) * e % spA % lbw(iEta,RIGHT)
               end do ; end do ; end do
   
            case (ELEFT)

               do eq = 1 , NCONS ; do iXi = 0 , e % spA % N ; do iEta = 0 , e % spA % N
                  val(iXi,iEta,eq) = F(iEta,eq) * e % spA % lbw(iXi,LEFT)
               end do ; end do ; end do

         end select

      end function StdFace

end module DGWeakIntegrals
