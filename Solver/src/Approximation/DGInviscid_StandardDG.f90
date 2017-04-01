submodule (DGInviscidMethods)       DGInviscid_StandardDG
   use Physics

#include "Defines.h"

   contains
      module subroutine StdDG_ComputeRiemannFluxes_Interior( self , ed , FL , FR )
         use MatrixOperations
         use QuadElementClass
         implicit none
         class(InviscidMethod_t)    :: self
         type(Edge_t)               :: ed
         real ( kind=RP )           :: FL ( 0 : ed % storage ( LEFT  ) % spA % N , 1 : NCONS )
         real ( kind=RP )           :: FR ( 0 : ed % storage ( RIGHT ) % spA % N , 1 : NCONS )
!        ---------------------------------------------------------------------------------------------
         real ( kind=RP ), target   :: QL ( 0 : ed                     % spA % N , 1 : NCONS )
         real ( kind=RP ), target   :: QR ( 0 : ed                     % spA % N , 1 : NCONS )
         real ( kind=RP )           :: QL1D(1:NCONS) , QR1D(1:NCONS)
         integer                    :: iXi

         if ( ed % transform(LEFT) .and. ed % inverse ) then
!
!           Transform the LEFT edge
!           -----------------------            
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
!
!           Invert the RIGHT edge 
!           ---------------------
            QR = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS )
!
!           Compute the Riemann solver
!           --------------------------
            do iXi = 0 , ed % spA % N
               QL1D = QL(iXi,:)
               QR1D = QR(iXi,:)
               FR(iXi,1:NCONS) = self % RiemannSolver(QL1D , QR1D , ed % n(IX:IY,0) ) * ed % dS(0)
            end do
!
!           Transform the LEFT 
!           -----------------------
            call Mat_x_Mat( ed % T_backward , FR , FL )
!
!           Invert the RIGHT 
!           ---------------------
            FR(0:ed % spA % N , 1:NCONS) = FR(ed % spA % N : 0 : -1 , 1:NCONS)

         elseif ( ed % transform(LEFT) ) then
!
!           Transform the LEFT 
!           -----------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(LEFT) % Q , QL)
!
!           Get the RIGHT  state
!           ------------------------
            QR = ed % storage(RIGHT) % Q
!
!           Compute the Riemann solver
!           --------------------------
            do iXi = 0 , ed % spA % N
               QL1D = QL(iXi,:)
               QR1D = QR(iXi,:)
               FR(iXi,1:NCONS) = self % RiemannSolver(QL1D , QR1D , ed % n(IX:IY,0) ) * ed % dS(0)
            end do
!
!           Transform the LEFT 
!           -----------------------
            call Mat_x_Mat( ed % T_backward , FR , FL )


         elseif ( ed % transform(RIGHT) .and.  ed % inverse ) then
!
!           Get the LEFT 
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Transform the RIGHT 
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
!
!           Invert the RIGHT 
!           ---------------------
            QR(0:ed % spA % N,1:NCONS) = QR(ed % spA % N : 0 : -1 , 1:NCONS)
!
!           Compute the Riemann solver
!           --------------------------
            do iXi = 0 , ed % spA % N
               QL1D = QL(iXi,:)
               QR1D = QR(iXi,:)
               FL(iXi,1:NCONS) = self % RiemannSolver(QL1D , QR1D , ed % n(IX:IY,0) ) * ed % dS(0)
            end do
!
!           Undo the transformation for the RIGHT 
!           ------------------------------------------ 
            call Mat_x_Mat( ed % T_backward , FL , FR )
!
!           Invert the 
!           ---------------
            FR(0:ed % NLow,1:NCONS) = FR(ed % NLow:0:-1 , 1:NCONS)

         elseif ( ed % transform(RIGHT) ) then
!
!           Get the LEFT 
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Transform the RIGHT 
!           ------------------------
            call Mat_x_Mat( ed % T_forward , ed % storage(RIGHT) % Q , QR)
!
!           Compute the Riemann solver
!           --------------------------
            do iXi = 0 , ed % spA % N
               QL1D = QL(iXi,:)
               QR1D = QR(iXi,:)
               FL(iXi,1:NCONS) = self % RiemannSolver(QL1D , QR1D , ed % n(IX:IY,0) ) * ed % dS(0)
            end do
!
!
!           Undo the transformation for the RIGHT 
!           ------------------------------------------
            call Mat_x_Mat( ed % T_backward , FL , FR )
         
         elseif ( ed % inverse ) then
!
!           Get the LEFT 
!           -----------------
            QL = ed % storage(LEFT) % Q
!
!           Invert the RIGHT 
!           ---------------------
            QR(0:ed % spA % N , 1:NCONS ) = ed % storage(RIGHT) % Q(ed % spA % N : 0 : -1 , 1:NCONS)
!
!           Compute the Riemann solver
!           --------------------------         
            do iXi = 0 , ed % spA % N
               QL1D = QL(iXi,:)
               QR1D = QR(iXi,:)
               FL(iXi,1:NCONS) = self % RiemannSolver(QL1D , QR1D , ed % n(IX:IY,0) ) * ed % dS(0)
            end do
!
!           Invert the RIGHT 
!           ---------------------
            FR( 0 : ed % spA % N , 1:NCONS ) = FL ( ed % spA % N : 0 : -1 , 1:NCONS )

         else
!
!           Get the LEFT 
!           -----------------   
            QL = ed % storage(LEFT) % Q
!
!           Get the RIGHT 
!           ------------------
            QR = ed % storage(RIGHT) % Q
!
!           Compute the Riemann solver
!           --------------------------
            do iXi = 0 , ed % spA % N
               QL1D = QL(iXi,:)
               QR1D = QR(iXi,:)
               FL(iXi,1:NCONS) = self % RiemannSolver(QL1D , QR1D , ed % n(IX:IY,0) ) * ed % dS(0)
               FR(iXi,1:NCONS) = FL(iXi,1:NCONS)
            end do
         
         end if
!
!        Change the sign of FR so that it points towards the element outside
!        -------------------------------------------------------------------
         FR = -1.0_RP * FR

      end subroutine StdDG_ComputeRiemannFluxes_Interior

      module subroutine StdDG_ComputeRiemannFluxes_StraightBdry( self , ed , F )
         use QuadElementClass
         implicit none
         class(InviscidMethod_t)       :: self
         type(StraightBdryEdge_t)      :: ed
         real(kind=RP)                 :: F( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: QL(1:NCONS) , QR(1:NCONS)
         integer, pointer              :: N
         integer                       :: iXi
         procedure(RiemannSolverFunction), pointer    :: RiemannSolver

         N => ed % spA % N
!
!        Select boundary condition type
!        ------------------------------
         select case ( ed % BCWeakType ) 

            case ( WEAK_PRESCRIBED )
!
!              Prescribed boundary conditions: just grab the value
!              ---------------------------------------------------
               F = ed % FB

            case ( WEAK_RIEMANN )
!
!              Weak boundary conditions
!              ------------------------
               if ( associated ( ed % RiemannSolver ) ) then
                  RiemannSolver => ed % RiemannSolver

               else
                  RiemannSolver => self % RiemannSolver

               end if
!
!              Select LEFT and RIGHT states
!              ----------------------------
               if ( .not. ed % inverse ) then       ! This just occurs in periodic BCs
                   do iXi = 0 , N
                     QL = ed % storage(1) % Q(iXi , 1:NCONS)
                     QR = ed % uB(iXi , 1:NCONS)

                     F(iXi , 1:NCONS) = RiemannSolver( QL , QR , ed % n(IX:IY,0) ) * ed % dS(0)
                  end do

               else
                  do iXi = 0 , N
                     QL = ed % storage(1) % Q(iXi , 1:NCONS)
                     QR = ed % uB(N - iXi , 1:NCONS)

                     F(iXi , 1:NCONS) = RiemannSolver( QL , QR , ed % n(IX:IY,0) ) * ed % dS(0)
                  end do

               end if
         end select

      end subroutine StdDG_ComputeRiemannFluxes_StraightBdry

      module subroutine StdDG_ComputeRiemannFluxes_CurvedBdry( self , ed , F)
         use QuadElementClass
         implicit none
         class(InviscidMethod_t)       :: self
         type(CurvedBdryEdge_t)        :: ed
         real(kind=RP)                 :: F( 0 : ed % spA % N , 1 : NCONS )
         real(kind=RP)                 :: QL(1:NCONS) , QR(1:NCONS)
         integer, pointer              :: N 
         integer                       :: iXi
         procedure(RiemannSolverFunction), pointer    :: RiemannSolver

         N => ed % spA % N
!
!        Select boundary condition type
!        ------------------------------
         select case ( ed % BCWeakType ) 

            case ( WEAK_PRESCRIBED )
!
!              Prescribed boundary conditions: just grab the value
!              ---------------------------------------------------
               F = ed % FB

            case ( WEAK_RIEMANN )
!
!              Weak boundary conditions
!              ------------------------
               if ( associated ( ed % RiemannSolver ) ) then
                  RiemannSolver => ed % RiemannSolver

               else
                  RiemannSolver => self % RiemannSolver

               end if
!
!              Select LEFT and RIGHT states
!              ----------------------------
               do iXi = 0 , N
                  QL = ed % storage(1) % Q(iXi , 1:NCONS)
                  QR = ed % uB(iXi , 1:NCONS)

                  F(iXi , 1:NCONS) = RiemannSolver( QL , QR , ed % n(IX:IY,iXi) ) * ed % dS(iXi)
               end do

         end select

      end subroutine StdDG_ComputeRiemannFluxes_CurvedBdry

end submodule DGInviscid_StandardDG
