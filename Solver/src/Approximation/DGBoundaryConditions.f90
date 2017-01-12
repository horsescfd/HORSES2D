!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
module DGBoundaryConditions
   use SMConstants
   use Physics
   use QuadElementClass
   use ParamfileIO
   implicit none

   private
   public BoundaryCondition_t , Construct
   public PeriodicBC_t , DirichletBC_t , Farfield_t , EulerWall_t

   integer, parameter         :: STR_LEN_BC = 128
!
!  **********************************
!  Base class for boundary conditions
!  **********************************
!
   type BoundaryCondition_t
      character(len=STR_LEN_BC)        :: Name
      integer                          :: BCType
      contains
         procedure ::     Construct => BaseClass_Construct
         procedure ::     Associate => BaseClass_Associate
         procedure ::     Update    => BaseClass_Update
   end type BoundaryCondition_t
!
!  *********************************
!  Periodic boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: PeriodicBC_t
   end type PeriodicBC_t
!
!  **********************************
!  Dirichlet boundary condition class
!  **********************************
!
   type, extends(BoundaryCondition_t)           :: DirichletBC_t
      real(kind=RP), dimension(NEC)       :: q
      contains
         procedure ::      Construct => DirichletBC_Construct
         procedure ::      Associate => DirichletBC_Associate
   end type DirichletBC_t
!
!  *********************************
!  Farfield boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: Farfield_t
      real(kind=RP)           :: u
      real(kind=RP)           :: v
      real(kind=RP)           :: p
      real(kind=RP)           :: T
   end type FarField_t
!
!  ***********************************
!  Euler wall boundary condition class
!  ***********************************
!
   type, extends(BoundaryCondition_t)           :: EulerWall_t
      contains
         procedure   ::    Associate => EulerWall_Associate
         procedure   ::    Update    => EulerWall_Update
   end type EulerWall_t

!
!  *******************
!  Construct procedure
!  *******************
!
   interface Construct
      module procedure BoundaryConditions_construct
   end interface Construct
!
!///////////////////////////////////////////////////////////////////////////////////
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////
!

      subroutine BoundaryConditions_construct( self , marker)
         use Setup_class
         implicit none
         class(BoundaryCondition_t), pointer :: self
         integer                             :: marker
         character(len=STR_LEN_BC)           :: BCType
         character(len=STR_LEN_BC)           :: in_label

         write(in_label,'(A,I0)') "# define zone ",marker

         call ReadValueInRegion( trim(Setup % bdry_file) , "Type" , BCType , in_label , "# end" )

         if (BCType .eq. "Dirichlet") then
            allocate( DirichletBC_t    :: self )
            self % BCType = DIRICHLET_BC
      
         elseif ( BCType .eq. "EulerWall") then
            allocate( EulerWall_t      :: self )
            self % BCType = EULERWALL_BC

         end if

         call self % Construct( marker , in_label)

      end subroutine BoundaryConditions_construct
   
      subroutine BaseClass_Construct( self , marker , in_label)
         implicit none
         class(BoundaryCondition_t)                :: self
         integer                                   :: marker
         character(len=*)                          :: in_label
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Construct

      subroutine BaseClass_Associate( self , edge )
         implicit none
         class(BoundaryCondition_t)          :: self
         class(Edge_t)                       :: edge
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Associate

      subroutine BaseClass_Update( self , edge )
         implicit none
         class(BoundaryCondition_t)          :: self
         class(Edge_t)                       :: edge
!
!        *****************************************
!           The base class does nothing
!        *****************************************
!
      end subroutine BaseClass_Update
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           DIRICHLET BC
!           ------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine DirichletBC_Construct( self , marker , in_label)
         use Setup_class
         implicit none
         class(DirichletBC_t)      :: self
         integer                   :: marker
         character(len=*)          :: in_label
         real(kind=RP), allocatable             :: pressure
         real(kind=RP), allocatable         :: Temperature
         real(kind=RP), allocatable         :: Mach
         real(kind=RP), allocatable         :: AngleOfAttack
         real(kind=RP)                          :: rho

         call readValueInRegion( trim(Setup % bdry_file) , "Name" , self % Name , in_label , "# end" )
         call readValueInRegion( trim(Setup % bdry_file) , "pressure" , pressure , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Temperature", Temperature , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Mach" , Mach , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Angle of attack" , AngleOfAttack , in_label , "# end")
         
         if ( allocated(pressure) ) then
            pressure = pressure / refValues % p
         else
            allocate(pressure)
            pressure = 1.0_RP
         end if

         if ( allocated(Temperature) ) then
            Temperature = Temperature / refValues % T
         else
            allocate(Temperature)
            Temperature = 1.0_RP
         end if

         if ( .not. allocated(Mach) ) then
            allocate(Mach)
            Mach = Dimensionless % Mach
         end if

         if ( allocated(AngleOfAttack) ) then
            AngleOfAttack = AngleOfAttack * PI / 180.0_RP
         else
            allocate(AngleOfAttack)
            AngleOfAttack = 0.0_RP
         end if
!
!        Construct the state vector
!        --------------------------
         associate ( gamma => Thermodynamics % Gamma , cv => Dimensionless % cv)
         rho = pressure / Temperature
         self % q(IRHO) = rho
         self % q(IRHOU) = rho * sqrt(gamma) * Mach * cos(AngleOfAttack)
         self % q(IRHOV) = rho * sqrt(gamma) * Mach * sin(AngleOfAttack)
         self % q(IRHOE) = cv * pressure + 0.5_RP * rho * gamma * Mach * Mach
         end associate

      end subroutine DirichletBC_Construct

      subroutine DirichletBC_Associate(self , edge)
         implicit none
         class(DirichletBC_t)          :: self
         class(Edge_t)                 :: edge
         integer                       :: i

         associate ( N => edge % spA % N )

         select type ( edge )
         
            type is (Edge_t)
               print*, "Only boundary edges are expected."
               stop "Stopped"
      
            type is (StraightBdryEdge_t)
               allocate( edge % uB(0:N,NEC) )
               allocate( edge % gB(0:N,NEC,NDIM) )   ! Normal gradients
               
               do i = 0 , N
                  edge % uB(i , 1:NEC) = self % q
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP
               
   
            type is (CurvedBdryEdge_t)
               allocate( edge % uB(0:N,NEC) )
               allocate( edge % gB(0:N,NEC,NDIM) )

               do i = 0 , N
                  edge % uB(i , 1:NEC) = self % q
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP

         end select
         end associate
      end subroutine DirichletBC_Associate
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           EULER WALL
!           ----------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine EulerWall_Associate( self , edge ) 
         implicit none  
         class(EulerWall_t)                  :: self
         class(Edge_t)                       :: edge

         associate( N => edge % spA % N )
         select type ( edge )
         
            type is (Edge_t)
               print*, "Only boundary edges are expected."
               stop "Stopped"

            type is (StraightBdryEdge_t)
               allocate( edge % FB(0:N,NEC) )
         
            type is (CurvedBdryEdge_t)
               allocate( edge % FB(0:N,NEC) )

            class default
         end select
         end associate
                
      end subroutine EulerWall_Associate


      subroutine EulerWall_Update( self , edge )
         implicit none  
         class(EulerWall_t)                  :: self
         class(Edge_t)                       :: edge
         real(kind=RP), allocatable          :: p(:)
!
!        ***************************************************
!           For each edge the flux is computed as:
!              F(IRHO)  = 0
!              F(IRHOU) = p dSx
!              F(IRHOV) = p dSy
!              F(IRHOE) = 0
!        ***************************************************
!
         associate( N => edge % spA % N , gm1 => Thermodynamics % gm1 , gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
   
         allocate( p(0:N) ) 
         
         p(0:N) = gm1 * ( edge % Q(0:N,IRHOE,1) - 0.5_RP * (edge % Q(0:N,IRHOU,1) * edge % Q(0:N,IRHOU,1) +  &
                                                            edge % Q(0:N,IRHOV,1) * edge % Q(0:N,IRHOV,1) / edge % Q(0:N,IRHO,1)) ) 

         select type ( edge ) 
            type is (StraightBdryEdge_t)

               edge % FB(0:N,IRHO)  = 0.0_RP
               edge % FB(0:N,IRHOU) = p(0:N) * edge % dS(iX,0:N) / (sqrt(gamma) * Mach)
               edge % FB(0:N,IRHOV) = p(0:N) * edge % dS(iY,0:N) / (sqrt(gamma) * Mach)
               edge % FB(0:N,IRHOE) = 0.0_RP
         
            type is (CurvedBdryEdge_t)

               edge % FB(0:N,IRHO)  = 0.0_RP
               edge % FB(0:N,IRHOU) = p(0:N) * edge % dS(iX,0:N) / (sqrt(gamma) * Mach)
               edge % FB(0:N,IRHOV) = p(0:N) * edge % dS(iY,0:N) / (sqrt(gamma) * Mach)
               edge % FB(0:N,IRHOE) = 0.0_RP
      
            class default
         end select

         deallocate( p ) 

         end associate

      end subroutine EulerWall_Update


!
!      subroutine DGBoundaryConditions_setFace( self  )
!         use Physics
!         use Setup_class
!         implicit none
!         class(BoundaryCondition_t)       :: self
!
!         select type (f1=>self % edge)
!            type is (StraightBdryEdge_t)
!               if (self % type  .eq. PERIODIC_BC) then
!!                 -------------------------------
!!                    Look for its pairing
!!                 -------------------------------
!                  select type (f2=>self % periodicPair % edge)
!                     type is (StraightBdryEdge_t)
!                        f1 % uB => NULL()    ! TODO f2 % quads(1) % e % Qb( : , f2 % BCLocation )  
!                        f1 % gB => NULL()    ! TODO f2 % quads(1) % e % dQb( : , f2 % BCLocation )  
!                  end select
!               elseif ( self % type .eq. DIRICHLET_BC) then
!!                 -------------------------------
!!                    Allocate data and set values
!!                 -------------------------------
!                        allocate( f1 % uB (NEC , 0 : f1 % spA % N) )
!                        f1 % uB = Setup % dirichletBC( f1 % edgeType )
!                        f1 % gB => NULL()    ! TODO f1 % quads(1) % e % dQb( : , f1 % BCLocation )
!
!               end if
!         end select
!
!         
!
!
!      end subroutine DGBoundaryConditions_setFace   
!
!
end module DGBoundaryConditions
