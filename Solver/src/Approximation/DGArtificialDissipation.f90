#ifdef NAVIER_STOKES

module DGArtificialDissipation
   use SMConstants
   use QuadMeshClass
   use QuadElementClass
   implicit none

#include "Defines.h"


   private
   public   ArtificialDissipation_t    , ArtificialDissipation_Initialization
!
!////////////////////////////////////////////////////////////////////////////
!
!           DISSIPATION MATRIX CLASSES
!           --------------------------
!////////////////////////////////////////////////////////////////////////////
!
   type DissipationMatrix_t
      contains
         procedure      :: Compute     => DissipationMatrix_Compute
   end type DissipationMatrix_t

   type, extends(DissipationMatrix_t)  :: ResidualsBasedDM_t
      real(kind=RP)        :: maxViscosity
      contains
         procedure      :: Compute     => ResidualsBasedDM_Compute
   end type ResidualsBasedDM_t

   type, extends(DissipationMatrix_t)  :: JumpsBasedDM_t
      real(kind=RP)        :: maxViscosity
      contains
         procedure      :: Compute  => JumpsBasedDM_Compute
   end type JumpsBasedDM_t
!
!////////////////////////////////////////////////////////////////////////////
!
!           ARTIFICIAL DISSIPATION CLASSES
!           ------------------------------
!////////////////////////////////////////////////////////////////////////////
!
   type ArtificialDissipation_t
      class(DissipationMatrix_t),   allocatable    :: DissipationMatrix
      contains
         procedure   :: ComputeVolumeFluxes  => ArtificialDissipation_ComputeVolumeFluxes
   end type ArtificialDissipation_t

   type, extends(ArtificialDissipation_t)    :: GreenVolumeAD_t
      contains
         procedure   :: ComputeVolumeFluxes  => GreenVolumeAD_ComputeVolumeFluxes
   end type GreenVolumeAD_t
!
!  ========
   contains
!  ========
!
!
!//////////////////////////////////////////////////////////////////////////////////////
!
!           INITIALIZATION
!           --------------
!//////////////////////////////////////////////////////////////////////////////////////
!
      function ArtificialDissipation_Initialization() result( ArtificialDissipation )
         use Setup_class
         implicit none
         class(ArtificialDissipation_t), pointer       :: ArtificialDissipation


         if ( Setup % artificialDissipation ) then
            allocate    ( GreenVolumeAD_t    :: ArtificialDissipation )
   
            if ( trim(Setup % artificialDissipationIndicator) .eq. "Residuals-based" ) then
               allocate    ( ResidualsBasedDM_t    :: ArtificialDissipation % DissipationMatrix )
   
            elseif ( trim(Setup % artificialDissipationIndicator) .eq. "Jumps-based" ) then
               allocate    ( JumpsBasedDM_t    :: ArtificialDissipation % DissipationMatrix )

            end if

            select type ( DM => ArtificialDissipation % DissipationMatrix )
               type is (ResidualsBasedDM_t)
                  DM % maxViscosity = Setup % artificialDissipation_maxViscosity
               
               type is (JumpsBasedDM_t)
                  DM % maxViscosity = Setup % artificialDissipation_maxViscosity

            end select
   
         else
            allocate ( ArtificialDissipation_t  :: ArtificialDissipation )

         end if

      end function ArtificialDissipation_Initialization
!
!//////////////////////////////////////////////////////////////////////////////////////
!
!           ARTIFICIAL DISSIPATION PROCEDURES
!           ---------------------------------
!//////////////////////////////////////////////////////////////////////////////////////
!
!TODO: pure
      function ArtificialDissipation_ComputeVolumeFluxes( self , e ) result ( F )
         implicit none
         class(ArtificialDissipation_t), intent(in)       :: self
         class(QuadElement_t)      , intent(in)       :: e
         real(kind=RP)                                :: F( 0 : e % spA % N , 0 : e % spA % N , 1 : NCONS , 1:NDIM)
!
!        ---------------------------
!        The base class does nothing.
!        ---------------------------
!
         F = 0.0_RP
      end function ArtificialDissipation_ComputeVolumeFluxes
! TODO: pure
      function GreenVolumeAD_ComputeVolumeFluxes ( self , e ) result ( F )
         implicit none
         class(GreenVolumeAD_t), intent(in)       :: self
         class(QuadElement_t)      , intent(in)       :: e
         real(kind=RP)                                :: F( 0 : e % spA % N , 0 : e % spA % N , 1 : NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!  
         real(kind=RP)                                :: eps( 0 : e % spA % N , 0 : e % spA % N , 1 : NCONS , 1:NDIM)

         eps = self % DissipationMatrix % Compute( e ) 

         F(:,:,:,IX) = e % dQ(:,:,IX,:) * eps(:,:,:,IX)
         F(:,:,:,IY) = e % dQ(:,:,IY,:) * eps(:,:,:,IY)
         

      end function GreenVolumeAD_ComputeVolumeFluxes

!
!//////////////////////////////////////////////////////////////////////////////////////
!
!           DISSIPATION MATRIX PROCEDURES
!           -----------------------------
!//////////////////////////////////////////////////////////////////////////////////////
!
!TODO: pure
      function DissipationMatrix_Compute( self , e ) result ( eps )
         implicit none
         class(DissipationMatrix_t), intent(in)       :: self
         class(QuadElement_t)      , intent(in)       :: e
         real(kind=RP)                                :: eps( 0 : e % spA % N , 0 : e % spA % N , 1 : NCONS , 1:NDIM)
!
!        ---------------------------
!        The base class does nothing.
!        ---------------------------
!
      end function DissipationMatrix_Compute
! TODO:pure
      function ResidualsBasedDM_Compute( self , e ) result ( eps )
         implicit none
         class(ResidualsBasedDM_t), intent(in)       :: self
         class(QuadElement_t)      , intent(in)       :: e
         real(kind=RP)                                :: eps( 0 : e % spA % N , 0 : e % spA % N , 1 : NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: hx , hy
         integer           :: eq

!         hx = self % C_eps * ( e % Volume / max( e % edges(ELEFT) % f % Area , e % edges(ERIGHT) % f % Area ) ) ** ( 2.0_RP - self % beta )
!         hy = self % C_eps * ( e % Volume / max( e % edges(ETOP) % f % Area , e % edges(EBOTTOM) % f % Area ) ) ** ( 2.0_RP - self % beta )

         do eq = 1 , NCONS
            eps(:,:,eq,IX) =  abs(e % QDot(:,:,IRHO)) * self % maxViscosity 
            eps(:,:,eq,IY) =  abs(e % QDot(:,:,IRHO)) * self % maxViscosity
         end do
      
      end function ResidualsBasedDM_Compute
!TODO: pure
      function JumpsBasedDM_Compute( self , e ) result ( eps )
         implicit none
         class(JumpsBasedDM_t), intent(in)      :: self
         class(QuadElement_t)      , intent(in) :: e
         real(kind=RP)                          :: eps( 0 : e % spA % N , 0 : e % spA % N , 1 : NCONS , 1:NDIM)

         real(kind=RP)     :: gk
         real(kind=RP)     :: capitalGk
         real(kind=RP)     :: jumpsIntegral
         integer           :: ed
         real(kind=RP)     :: xiL = 1.0e-5_RP
         real(kind=RP)     :: xiR = 1.0e-4_RP
!
!        The jumps indicator is computed squaring all interface jumps
!        ------------------------------------------------------------
         gk =   e % edges(EBOTTOM) % f % computeJumps(IRHO) * e % edges(EBOTTOM) % f % invh   &
              + e % edges(ERIGHT ) % f % computeJumps(IRHO) * e % edges(ERIGHT ) % f % invh   & 
              + e % edges(ETOP   ) % f % computeJumps(IRHO) * e % edges(ETOP   ) % f % invh   & 
              + e % edges(ELEFT  ) % f % computeJumps(IRHO) * e % edges(ELEFT  ) % f % invh  

         gk = gk / e % Volume
!
!        Compute the smooth discrete jumps indicator
!        -------------------------------------------
         if ( gk .lt. xiL ) then
            capitalGk = 0.0_RP

         elseif ( gk .lt. xiR ) then
            capitalGk = 0.5_RP * sin( 0.5_RP * PI * ( gk - (xiR-xiL) )/(xiR-xiL) ) + 0.5_RP

         else
            capitalGk = 1.0_RP

         end if

         eps = self % maxViscosity * capitalGk

   if ( e % ID .eq. 426 ) then
      print*, gk , capitalGk
end if

      end function JumpsBasedDM_Compute

!
!//////////////////////////////////////////////////////////////////////////////////////
!
end module DGArtificialDissipation

#endif
