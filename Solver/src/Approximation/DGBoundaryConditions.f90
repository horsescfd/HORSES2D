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
         procedure ::     Describe  => BaseClass_Describe
   end type BoundaryCondition_t
!
!  *********************************
!  Periodic boundary condition class
!  *********************************
!
   type, extends(BoundaryCondition_t)           :: PeriodicBC_t
      logical                          :: associated
      integer                          :: direction
      integer                          :: connected_marker
      contains
         procedure ::     Construct => PeriodicBC_Construct
         procedure ::     Describe  => PeriodicBC_Describe
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
         procedure ::      Describe  => DirichletBC_Describe
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
         procedure   ::    Describe  => EulerWall_Describe
   end type EulerWall_t
!
!  *************************************
!  Viscous wall boundary condition class
!  *************************************
!
   type, extends(BoundaryCondition_t)           :: ViscousWall_t
      contains
         procedure   ::    Associate => ViscousWall_Associate
         procedure   ::    Update    => ViscousWall_Update
         procedure   ::    Describe  => ViscousWall_Describe
   end type ViscousWall_t


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
      
         elseif ( BCType .eq. "Periodic") then
            allocate( PeriodicBC_t     :: self )
            self % BCType = PERIODIC_BC

         elseif ( BCType .eq. "EulerWall") then
            allocate( EulerWall_t      :: self )
            self % BCType = EULERWALL_BC

         elseif ( BCType .eq. "ViscousWall") then
            allocate( ViscousWall_t    :: self ) 
            self % BCType = VISCOUSWALL_BC

         end if
   
         call readValueInRegion( trim(Setup % bdry_file) , "Name" , self % Name , in_label , "# end" )

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
!           PERIODIC BC
!           -----------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine PeriodicBC_Construct ( self , marker , in_label) 
         use Setup_class
         implicit none
         class(PeriodicBC_t)        :: self
         integer                   :: marker
         character(len=*)          :: in_label
         character(len=STR_LEN_BC)  :: direction
         integer, allocatable          :: connected_marker

         call readValueInRegion( trim(Setup % bdry_file) , "Direction" , direction , in_label , "# end")
         call readValueInRegion( trim(Setup % bdry_file) , "Marker" , connected_marker , in_label , "# end")

         if ( allocated( connected_marker ) ) then
            self % connected_marker = connected_marker
         else
            print*, "You need to specify a marker to connect zone " , marker , " with."
            stop "Stopped."
         end if         

         if ( trim(direction) .eq. "x" ) then
            self % direction = IX
         elseif ( trim(direction) .eq. "y") then
            self % direction = IY
         else
            print*, "Direction in zone ", marker, " must be 'x' or 'y'."
            stop "Stopped."
         end if
         
         self % associated = .false.

      end subroutine PeriodicBC_Construct

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
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              VISCOUS WALL BC
!              ---------------
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ViscousWall_Associate(self , edge)
         implicit none
         class(ViscousWall_t)          :: self
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
   
            type is (CurvedBdryEdge_t)
               allocate( edge % uB(0:N,NEC) )
               allocate( edge % gB(0:N,NEC,NDIM) )

         end select
         end associate
      end subroutine ViscousWall_Associate

      subroutine ViscousWall_Update( self , edge )
         implicit none  
         class(ViscousWall_t)                  :: self
         class(Edge_t)                       :: edge
!
!        ***************************************************
!           For each edge the ghost cell state is computed as
!              UB(IRHO) = Q(IRHO)
!              UB(IRHOU) = -Q(IRHOU)
!              UB(IRHOV) = -Q(IRHOV)
!              UB(IRHOE) = Q(IRHOE)
!        ***************************************************
!
         associate( N => edge % spA % N , gm1 => Thermodynamics % gm1 , gamma => Thermodynamics % gamma , Mach => Dimensionless % Mach )
   
         select type ( edge ) 
            type is (StraightBdryEdge_t)
   
               edge % uB(0 : N , IRHO)  =  edge % Q(0 : N , IRHO  , 1)
               edge % uB(0 : N , IRHOU) = -edge % Q(0 : N , IRHOU , 1)
               edge % uB(0 : N , IRHOV) = -edge % Q(0 : N , IRHOV , 1)
               edge % uB(0 : N , IRHOE) =  edge % Q(0 : N , IRHOE , 1)

            type is (CurvedBdryEdge_t)

               edge % uB(0 : N , IRHO)  =  edge % Q(0 : N , IRHO  , 1)
               edge % uB(0 : N , IRHOU) = -edge % Q(0 : N , IRHOU , 1)
               edge % uB(0 : N , IRHOV) = -edge % Q(0 : N , IRHOV , 1)
               edge % uB(0 : N , IRHOE) =  edge % Q(0 : N , IRHOE , 1)
      
            class default
         end select

         end associate

      end subroutine ViscousWall_Update
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              Describe subroutines
!              --------------------
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_Describe( self ) 
         implicit none
         class(BoundaryCondition_t)          :: self

      end subroutine BaseClass_Describe

      subroutine DirichletBC_Describe( self )
         use Physics
         implicit none
         class(DirichletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Dirichlet."
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         end associate

      end subroutine DirichletBC_Describe

      subroutine PeriodicBC_Describe( self ) 
         implicit none
         class(PeriodicBC_t)          :: self

         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Periodic."
         write(STD_OUT , '(30X,A,A25,I0,A)') "-> ", "Connecting zone: ",self % connected_marker
         if ( self % direction .eq. IX) then
            write(STD_OUT , '(30X,A,A25,A)') "-> " , "Direction: " ,"x"
         elseif ( self % direction .eq. IY) then
            write(STD_OUT , '(30X,A,A25,A)') "-> " , "Direction: " ,"y"
         end if
         
      end subroutine PeriodicBC_Describe

      subroutine EulerWall_Describe( self )
         implicit none
         class(EulerWall_t)                :: self

         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Euler wall."

      end subroutine EulerWall_Describe

      subroutine ViscousWall_Describe( self )
         implicit none
         class(ViscousWall_t)                :: self

         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Viscous wall."

      end subroutine ViscousWall_Describe

end module DGBoundaryConditions
