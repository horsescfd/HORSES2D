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
   public PeriodicBC_t , DirichletBC_t , FarfieldBC_t , EulerWall_t , PressureOutletBC_t , PressureInletBC_t

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
   type, extends(BoundaryCondition_t)           :: FarfieldBC_t
      real(kind=RP), dimension(NEC)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      contains
         procedure ::      Construct => FarfieldBC_Construct
         procedure ::      Associate => FarfieldBC_Associate
         procedure ::      Update    => FarfieldBC_Update
         procedure ::      Describe  => FarfieldBC_Describe
   end type FarFieldBC_t
!
!  ********************************
!  Outflow boundary condition class
!  ********************************
!
   type, extends(BoundaryCondition_t)           :: PressureOutletBC_t
      real(kind=RP), dimension(NEC)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      contains
         procedure ::      Construct => PressureOutletBC_Construct
         procedure ::      Associate => PressureOutletBC_Associate
         procedure ::      Update    => PressureOutletBC_Update
         procedure ::      Describe  => PressureOutletBC_Describe
   end type PressureOutletBC_t
!
!  ***************************************
!  Inflow/Outflow boundary condition class
!  ***************************************
!
   type, extends(BoundaryCondition_t)           :: PressureInletBC_t
      real(kind=RP), dimension(NEC)       :: q
      real(kind=RP)                       :: AngleOfAttack
      real(kind=RP)                       :: Tt
      real(kind=RP)                       :: pt
      real(kind=RP)                       :: Ht
      contains
         procedure ::      Construct => PressureInletBC_Construct
         procedure ::      Associate => PressureInletBC_Associate
         procedure ::      Update    => PressureInletBC_Update
         procedure ::      Describe  => PressureInletBC_Describe
   end type PressureInletBC_t

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
   
         elseif (BCType .eq. "Farfield") then
            allocate( FarfieldBC_t     :: self )
            self % BCType = FARFIELD_BC

         elseif (BCType .eq. "Pressure outlet") then
            allocate( PressureOutletBC_t     :: self )
            self % BCType = OUTFLOW_BC
      
         elseif (BCType .eq. "Pressure inlet") then
            allocate( PressureInletBC_t     :: self )
            self % BCType = INFLOWOUTFLOW_BC

         elseif ( BCType .eq. "Periodic") then
            allocate( PeriodicBC_t     :: self )
            self % BCType = PERIODIC_BC

         elseif ( BCType .eq. "Euler wall") then
            allocate( EulerWall_t      :: self )
            self % BCType = EULERWALL_BC

         elseif ( BCType .eq. "Viscous wall") then
            allocate( ViscousWall_t    :: self ) 
            self % BCType = VISCOUSWALL_BC

         else
            print*, 'Boundary condition "',trim(BCType),'" not implemented yet.'
            print*, "Options available are:"
            print*, "   * Dirichlet"
            print*, "   * Farfield"
            print*, "   * Pressure outlet"
            print*, "   * Pressure inlet"
            print*, "   * Periodic"
            print*, "   * Euler wall"
            print*, "   * Viscous wall"
            stop "Stopped."

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
!           FARFIELD BC
!           -----------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine FarfieldBC_Construct( self , marker , in_label)
         use Setup_class
         implicit none
         class(FarfieldBC_t)      :: self
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
         associate ( gamma => Thermodynamics % Gamma , cv => Dimensionless % cv , gm1 => Thermodynamics % gm1 , gogm1 => Thermodynamics % gogm1 )
         rho = pressure / Temperature
         self % q(IRHO) = rho
         self % q(IRHOU) = rho * sqrt(gamma) * Mach * cos(AngleOfAttack)
         self % q(IRHOV) = rho * sqrt(gamma) * Mach * sin(AngleOfAttack)
         self % q(IRHOE) = cv * pressure + 0.5_RP * rho * gamma * Mach * Mach
         self % AngleOfAttack = AngleOfAttack
         self % Tt            = Temperature * ( 1.0_RP + 0.5_RP * gm1 * Mach * Mach)
         self % pt            = pressure * ( self % Tt / Temperature ) ** gogm1 
         end associate

      end subroutine FarfieldBC_Construct

      subroutine FarfieldBC_Associate(self , edge)
         implicit none
         class(FarfieldBC_t)          :: self
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
                  edge % uB(i , 1:NEC) = 0.0_RP      ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP
               
   
            type is (CurvedBdryEdge_t)
               allocate( edge % uB(0:N,NEC) )
               allocate( edge % gB(0:N,NEC,NDIM) )

               do i = 0 , N
                  edge % uB(i , 1:NEC) = 0.0_RP    ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP

         end select
         end associate
      end subroutine FarfieldBC_Associate

      subroutine FarfieldBC_Update( self , edge )
         implicit none
         class(FarfieldBC_t)          :: self
         class(Edge_t)                       :: edge
         integer                       :: iXi
         integer                    :: N 
         real(kind=RP)              :: rhoL , vnL , uL , vL , pL , ML , aL 
         real(kind=RP)              :: rhoR , uR , vR , pR , VtR , TR , aR , MR
         real(kind=RP)              :: rhoInfty, uInfty, vInfty , pInfty , ptInfty , TtInfty
         real(kind=RP)              :: acoef , bcoef , ccoef , Ht
         real(kind=RP)              :: nInfty(NDIM)
         real(kind=RP)              :: Rminus
!
!        *********************************************************************
!           This routine computes the "Right" state of a Farfield boundary
!          condition. Once this is done, the Riemann flux is computed
!          from the computed state, and the real boundary state. This
!          ficticial state is computed by means of the characteristics method
!          which yields in four cases:
!              * Supersonic inflow (Four entering characteristics)
!              * Subsonic inflow (Three entering characteristics)
!              * Subsonic outflow (One entering characteristics)
!              * Supersonic outflow (All characteristics leave the domain)
!        *********************************************************************
!
         associate ( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 , cp => Dimensionless % cp , cv => Dimensionless % cv)

         rhoInfty = self % q(IRHO)
         uInfty   = self % q(IRHOU) / rhoInfty
         vInfty   = self % q(IRHOV) / rhoInfty
         pInfty   = gm1 * ( self % q(IRHOE) - 0.5_RP * self % q(IRHOU) * uInfty - 0.5_RP * self % q(IRHOV) * vInfty )
         nInfty   = [uInfty,vInfty] / norm2([uInfty,vInfty])

         N = edge % spA % N

         select type ( edge )
            type is (StraightBdryEdge_t) 
               do iXi = 0 , N
!
!                 First stage: Determine the boundary flow character
!                 --------------------------------------------------
                  rhoL = edge % Q(iXi , IRHO , 1)
                  uL  = edge % Q(iXi , IRHOU , 1) / rhoL
                  vL  = edge % Q(iXi , IRHOV , 1) / rhoL
                  vnL = uL * edge % n (IX , iXi) + vL * edge % n(IY, iXi)
                  pL  = gm1 * ( edge % Q(iXi , IRHOE , 1) - 0.5_RP * edge % Q(iXi,IRHOU , 1) * uL - 0.5_RP * edge % Q(iXi,IRHOV , 1) * vL )
                  aL  = sqrt( gamma * pL / rhoL ) 
                  ML  = vnL / aL
!
!                 Second stage: Compute the "Right" state depending on the result
!                 ---------------------------------------------------------------
                  if ( ML .lt. -1.0_RP ) then      ! Supersonic inflow
                     edge % uB(iXi , 1:NEC)  = self % q
               
                  elseif ( ML .lt. 0.0_RP ) then   ! Subsonic inflow
                     pR = 0.5_RP * ( pL + pInfty - rhoL * aL * ( edge % n(IX,iXi) * (uL-uInfty) + edge % n(IY,iXi) * (vL - vInfty) ) )
                     rhoR = rhoInfty + (pR - pInfty) / (aL * aL)
                     uR = uInfty + edge % n(IX,iXi) * (pR-pInfty) / (aL * rhoL)
                     vR = vInfty + edge % n(IY,iXi) * (pR-pInfty) / (aL * rhoL) 

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )
                     
                  elseif ( ML .lt. 1.0_RP ) then   ! Subsonic outflow
!
!                 ****************************************************************************************************
!                    According to TAU Technical Doc.
!                 ****************************************************************************************************
!
                     pR = pInfty
                     rhoR = rhoL + (pInfty - pL) / (aL * aL)
                     uR = uL - edge % n(IX , iXi) * (pInfty - pL) / (rhoL * aL)
                     vR = vL - edge % n(IY , iXi) * (pInfty - pL) / (rhoL * aL)

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )


                  elseif ( ML .ge. 1.0_RP ) then   ! Supersonic outflow
                     edge % uB(iXi , 1:NEC) = edge % Q(iXi, 1:NEC , 1)
      
                  end if
               
      
               end do 

            type is (CurvedBdryEdge_t) 
               do iXi = 0 , N
!
!                 First stage: Determine the boundary flow character
!                 --------------------------------------------------
                  rhoL = edge % Q(iXi , IRHO , 1)
                  uL  = edge % Q(iXi , IRHOU , 1) / rhoL
                  vL  = edge % Q(iXi , IRHOV , 1) / rhoL
                  vnL = uL * edge % n (IX , iXi) + vL * edge % n(IY, iXi)
                  pL  = gm1 * ( edge % Q(iXi , IRHOE , 1) - 0.5_RP * edge % Q(iXi,IRHOU , 1) * uL - 0.5_RP * edge % Q(iXi,IRHOV , 1) * vL )
                  aL  = sqrt( gamma * pL / rhoL ) 
                  ML  = vnL / aL
!
!                 Second stage: Compute the "Right" state depending on the result
!                 ---------------------------------------------------------------
                  if ( ML .lt. -1.0_RP ) then      ! Supersonic inflow
                     edge % uB(iXi , 1:NEC)  = self % q
               
                  elseif ( ML .lt. 0.0_RP ) then   ! Subsonic inflow
                     pR = 0.5_RP * ( pL + pInfty - rhoL * aL * ( edge % n(IX,iXi) * (uL-uInfty) + edge % n(IY,iXi) * (vL - vInfty) ) )
                     rhoR = rhoInfty + (pR - pInfty) / (aL * aL)
                     uR = uInfty + edge % n(IX,iXi) * (pR-pInfty) / (aL * rhoL)
                     vR = vInfty + edge % n(IY,iXi) * (pR-pInfty) / (aL * rhoL) 

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )
                     
                  elseif ( ML .lt. 1.0_RP ) then   ! Subsonic outflow
!
!                 ****************************************************************************************************
!                    According to TAU Technical Doc.
!                 ****************************************************************************************************
!
                     pR = pInfty
                     rhoR = rhoL + (pInfty - pL) / (aL * aL)
                     uR = uL - edge % n(IX , iXi) * (pInfty - pL) / (rhoL * aL)
                     vR = vL - edge % n(IY , iXi) * (pInfty - pL) / (rhoL * aL)

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )


                  elseif ( ML .ge. 1.0_RP ) then   ! Supersonic outflow
                     edge % uB(iXi , 1:NEC) = edge % Q(iXi, 1:NEC , 1)
      
                  end if
               
      
               end do 

            class default
         end select



         end associate
         
      end subroutine FarfieldBC_Update
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           PRESSURE OUTLET OUTFLOW BC
!           ----------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine PressureOutletBC_Construct( self , marker , in_label)
         use Setup_class
         implicit none
         class(PressureOutletBC_t)      :: self
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
         associate ( gamma => Thermodynamics % Gamma , cv => Dimensionless % cv , gm1 => Thermodynamics % gm1 , gogm1 => Thermodynamics % gogm1 )
         rho = pressure / Temperature
         self % q(IRHO) = rho
         self % q(IRHOU) = rho * sqrt(gamma) * Mach * cos(AngleOfAttack)
         self % q(IRHOV) = rho * sqrt(gamma) * Mach * sin(AngleOfAttack)
         self % q(IRHOE) = cv * pressure + 0.5_RP * rho * gamma * Mach * Mach
         self % AngleOfAttack = AngleOfAttack
         self % Tt            = Temperature * ( 1.0_RP + 0.5_RP * gm1 * Mach * Mach)
         self % pt            = pressure * ( self % Tt / Temperature ) ** gogm1 
         end associate

      end subroutine PressureOutletBC_Construct

      subroutine PressureOutletBC_Associate(self , edge)
         implicit none
         class(PressureOutletBC_t)          :: self
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
                  edge % uB(i , 1:NEC) = 0.0_RP      ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP
               
   
            type is (CurvedBdryEdge_t)
               allocate( edge % uB(0:N,NEC) )
               allocate( edge % gB(0:N,NEC,NDIM) )

               do i = 0 , N
                  edge % uB(i , 1:NEC) = 0.0_RP    ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP

         end select
         end associate
      end subroutine PressureOutletBC_Associate

      subroutine PressureOutletBC_Update( self , edge )
         implicit none
         class(PressureOutletBC_t)          :: self
         class(Edge_t)                       :: edge
         integer                       :: iXi
         integer                    :: N 
         real(kind=RP)              :: rhoL , vnL , uL , vL , pL , ML , aL  , TL
         real(kind=RP)              :: rhoR , uR , vR , pR , VtR , TR , aR , MR
         real(kind=RP)              :: rhoInfty, uInfty, vInfty , pInfty , ptInfty , TtInfty
         real(kind=RP)              :: acoef , bcoef , ccoef , Ht
         real(kind=RP)              :: nInfty(NDIM)
         real(kind=RP)              :: Rminus
!
!        *********************************************************************
!           This routine computes the "Right" state of a Outflow boundary
!          condition. Once this is done, the Riemann flux is computed
!          from the computed state, and the real boundary state. This
!          ficticial state is computed by means of the characteristics method
!          which yields in four cases:
!              * Supersonic inflow (Four entering characteristics)
!              * Subsonic inflow (Three entering characteristics)
!              * Subsonic outflow (One entering characteristics)
!              * Supersonic outflow (All characteristics leave the domain)
!        *********************************************************************
!
         associate ( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 , cp => Dimensionless % cp , cv => Dimensionless % cv)

         rhoInfty = self % q(IRHO)
         uInfty   = self % q(IRHOU) / rhoInfty
         vInfty   = self % q(IRHOV) / rhoInfty
         pInfty   = gm1 * ( self % q(IRHOE) - 0.5_RP * self % q(IRHOU) * uInfty - 0.5_RP * self % q(IRHOV) * vInfty )
         nInfty   = [cos(self % AngleOfAttack) , sin(self % AngleOfAttack)]
         

         N = edge % spA % N

         select type ( edge )
            type is (StraightBdryEdge_t) 
               do iXi = 0 , N
!
!                 First stage: Determine the boundary flow character
!                 --------------------------------------------------
                  rhoL = edge % Q(iXi , IRHO , 1)
                  uL  = edge % Q(iXi , IRHOU , 1) / rhoL
                  vL  = edge % Q(iXi , IRHOV , 1) / rhoL
                  vnL = uL * edge % n (IX , iXi) + vL * edge % n(IY, iXi)
                  pL  = gm1 * ( edge % Q(iXi , IRHOE , 1) - 0.5_RP * edge % Q(iXi,IRHOU , 1) * uL - 0.5_RP * edge % Q(iXi,IRHOV , 1) * vL )
                  aL  = sqrt( gamma * pL / rhoL ) 
                  TL  = pL / rhoL
                  ML  = vnL / aL
!
!                 Second stage: Compute the "Right" state depending on the result
!                 ---------------------------------------------------------------
                  if ( ML .lt. -1.0_RP ) then      ! Supersonic inflow
                     edge % uB(iXi , 1:NEC)  = self % q
               
                  elseif ( ML .lt. 0.0_RP ) then   ! Subsonic inflow
                     pR = pInfty
                     uR = abs(vnL) * nInfty(IX)
                     vR = abs(vnL) * nInfty(IY)
                     rhoR = pR / TL 

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )
                     
                  elseif ( ML .lt. 1.0_RP ) then   ! Subsonic outflow
                     pR = pInfty
                     uR = uL + edge % n(IX,iXi) * (pL - pR) / (rhoL * aL)
                     vR = vL + edge % n(IY,iXi) * (pL - pR) / (rhoL * aL)
                     rhoR = rhoL + (pR - pL)/(aL*aL)

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )

                  elseif ( ML .ge. 1.0_RP ) then   ! Supersonic outflow
                     edge % uB(iXi , 1:NEC) = edge % Q(iXi, 1:NEC , 1)
      
                  end if
               
      
               end do 

            type is (CurvedBdryEdge_t) 
               do iXi = 0 , N
!
!                 First stage: Determine the boundary flow character
!                 --------------------------------------------------
                  rhoL = edge % Q(iXi , IRHO , 1)
                  uL  = edge % Q(iXi , IRHOU , 1) / rhoL
                  vL  = edge % Q(iXi , IRHOV , 1) / rhoL
                  vnL = uL * edge % n (IX , iXi) + vL * edge % n(IY, iXi)
                  pL  = gm1 * ( edge % Q(iXi , IRHOE , 1) - 0.5_RP * edge % Q(iXi,IRHOU , 1) * uL - 0.5_RP * edge % Q(iXi,IRHOV , 1) * vL )
                  aL  = sqrt( gamma * pL / rhoL ) 
                  TL  = pL / rhoL
                  ML  = vnL / aL
!
!                 Second stage: Compute the "Right" state depending on the result
!                 ---------------------------------------------------------------
                  if ( ML .lt. -1.0_RP ) then      ! Supersonic inflow
                     edge % uB(iXi , 1:NEC)  = self % q
               
                  elseif ( ML .lt. 0.0_RP ) then   ! Subsonic inflow
                     pR = pInfty
                     uR = abs(vnL) * nInfty(IX)
                     vR = abs(vnL) * nInfty(IY)
                     rhoR = pR / TL 

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )
                     
                  elseif ( ML .lt. 1.0_RP ) then   ! Subsonic outflow
                     pR = pInfty
                     uR = uL + edge % n(IX,iXi) * (pL - pR) / (rhoL * aL)
                     vR = vL + edge % n(IY,iXi) * (pL - pR) / (rhoL * aL)
                     rhoR = rhoL + (pR - pL)/(aL*aL)

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )

                  elseif ( ML .ge. 1.0_RP ) then   ! Supersonic outflow
                     edge % uB(iXi , 1:NEC) = edge % Q(iXi, 1:NEC , 1)
      
                  end if
               
      
               end do 

            class default
         end select



         end associate

      end subroutine PressureOutletBC_Update
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           INFLOW OUTFLOW BC
!           -----------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine PressureInletBC_Construct( self , marker , in_label)
         use Setup_class
         implicit none
         class(PressureInletBC_t)      :: self
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
         associate ( gamma => Thermodynamics % Gamma , cv => Dimensionless % cv , gm1 => Thermodynamics % gm1 , gogm1 => Thermodynamics % gogm1 )
         rho = pressure / Temperature
         self % q(IRHO) = rho
         self % q(IRHOU) = rho * sqrt(gamma) * Mach * cos(AngleOfAttack)
         self % q(IRHOV) = rho * sqrt(gamma) * Mach * sin(AngleOfAttack)
         self % q(IRHOE) = cv * pressure + 0.5_RP * rho * gamma * Mach * Mach
         self % AngleOfAttack = AngleOfAttack
         self % Tt            = Temperature * ( 1.0_RP + 0.5_RP * gm1 * Mach * Mach)
         self % pt            = pressure * ( self % Tt / Temperature ) ** gogm1 
         self % Ht            = Temperature * Dimensionless % cp + 0.5_RP * gamma * Mach * Mach
         end associate

      end subroutine PressureInletBC_Construct

      subroutine PressureInletBC_Associate(self , edge)
         implicit none
         class(PressureInletBC_t)          :: self
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
                  edge % uB(i , 1:NEC) = 0.0_RP      ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP
               
   
            type is (CurvedBdryEdge_t)
               allocate( edge % uB(0:N,NEC) )
               allocate( edge % gB(0:N,NEC,NDIM) )

               do i = 0 , N
                  edge % uB(i , 1:NEC) = 0.0_RP    ! Its value is not given until the update routine is invoked
               end do

               edge % gB(0:N , 1:NEC , 1:NDIM ) = 0.0_RP

         end select
         end associate
      end subroutine PressureInletBC_Associate

      subroutine PressureInletBC_Update( self , edge )
         implicit none
         class(PressureInletBC_t)          :: self
         class(Edge_t)                       :: edge
         integer                       :: iXi
         integer                    :: N 
         real(kind=RP)              :: rhoL , vnL , uL , vL , pL , ML , aL  , TL
         real(kind=RP)              :: rhoR , vnR , uR , vR , pR , VtR , TR , aR , MR
         real(kind=RP)              :: rhoInfty, uInfty, vInfty , pInfty 
         real(kind=RP)              :: a1,b1,c1,d1,a2,b2,c2 
         real(kind=RP)              :: alpha , beta
         real(kind=RP)              :: nInfty(NDIM)
!
!        *********************************************************************
!           This routine computes the "Right" state of a InflowOutflow boundary
!          condition. Once this is done, the Riemann flux is computed
!          from the computed state, and the real boundary state. This
!          ficticial state is computed by means of the characteristics method
!          which yields in four cases:
!              * Supersonic inflow (Four entering characteristics)
!              * Subsonic inflow (Three entering characteristics)
!              * Subsonic outflow (One entering characteristics)
!              * Supersonic outflow (All characteristics leave the domain)
!        *********************************************************************
!
         associate ( gamma => Thermodynamics % gamma , gm1 => Thermodynamics % gm1 , cp => Dimensionless % cp , cv => Dimensionless % cv)

         rhoInfty = self % q(IRHO)
         uInfty   = self % q(IRHOU) / rhoInfty
         vInfty   = self % q(IRHOV) / rhoInfty
         pInfty   = gm1 * ( self % q(IRHOE) - 0.5_RP * self % q(IRHOU) * uInfty - 0.5_RP * self % q(IRHOV) * vInfty )
         nInfty   = [cos(self % AngleOfAttack) , sin(self % AngleOfAttack)]
         

         N = edge % spA % N

         select type ( edge )
            type is (StraightBdryEdge_t) 
               do iXi = 0 , N
!
!                 First stage: Determine the boundary flow character
!                 --------------------------------------------------
                  rhoL = edge % Q(iXi , IRHO , 1)
                  uL  = edge % Q(iXi , IRHOU , 1) / rhoL
                  vL  = edge % Q(iXi , IRHOV , 1) / rhoL
                  vnL = uL * edge % n (IX , iXi) + vL * edge % n(IY, iXi)
                  pL  = gm1 * ( edge % Q(iXi , IRHOE , 1) - 0.5_RP * edge % Q(iXi,IRHOU , 1) * uL - 0.5_RP * edge % Q(iXi,IRHOV , 1) * vL )
                  aL  = sqrt( gamma * pL / rhoL ) 
                  TL  = pL / rhoL
                  ML  = vnL / aL
!
!                 Second stage: Compute the "Right" state depending on the result
!                 ---------------------------------------------------------------
                  if ( ML .lt. -1.0_RP ) then      ! Supersonic inflow
                     edge % uB(iXi , 1:NEC)  = self % q
               
                  elseif ( ML .le. 0.0_RP ) then   ! Subsonic inflow
                     alpha = rhoInfty - pInfty / (aL*aL)
                     beta  = vnL - pL / (rhoL * aL) 

                     a1 = 0.5_RP/(rhoL*rhoL * aL*aL*aL*aL);
                     b1 = beta / (rhoL * aL*aL*aL) + alpha / (2*rhoL*rhoL*aL*aL);
                     c1 = beta * beta / ( 2 * aL * aL) + alpha * beta / ( rhoL * aL ) + Dimensionless % cp - self % Ht / ( aL * aL);
                     d1= alpha * beta * beta / 2 + self % Ht * pInfty / ( aL * aL ) - rhoInfty * self % Ht;

                     a2 = b1 / a1
                     b2 = c1 / a1
                     c2 = d1 / a1
                     
                     pR = ThirdDegreeRoots(a2,b2,c2)

                     rhoR = rhoInfty + (pR - pInfty) / (aL*aL)
                     vnR = vnL + (pR - pInfty) / (rhoL*aL)
                     uR = abs(vnR) * nInfty(IX)
                     vR = abs(vnR) * nInfty(IY)

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = edge % uB(iXi , IRHO) * uR
                     edge % uB(iXi , IRHOV) = edge % uB(iXi , IRHO) * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )
                     
                  elseif ( ML .lt. 1.0_RP ) then   ! Subsonic outflow
                     pR = pInfty
                     uR = uL
                     vR = vL
                     rhoR = pR / TL

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )

                  elseif ( ML .ge. 1.0_RP ) then   ! Supersonic outflow
                     edge % uB(iXi , 1:NEC) = edge % Q(iXi, 1:NEC , 1)
      
                  end if
               
      
               end do 

            type is (CurvedBdryEdge_t) 
               do iXi = 0 , N
!
!                 First stage: Determine the boundary flow character
!                 --------------------------------------------------
                  rhoL = edge % Q(iXi , IRHO , 1)
                  uL  = edge % Q(iXi , IRHOU , 1) / rhoL
                  vL  = edge % Q(iXi , IRHOV , 1) / rhoL
                  vnL = uL * edge % n (IX , iXi) + vL * edge % n(IY, iXi)
                  pL  = gm1 * ( edge % Q(iXi , IRHOE , 1) - 0.5_RP * edge % Q(iXi,IRHOU , 1) * uL - 0.5_RP * edge % Q(iXi,IRHOV , 1) * vL )
                  aL  = sqrt( gamma * pL / rhoL ) 
                  TL  = pL / rhoL
                  ML  = vnL / aL
!
!                 Second stage: Compute the "Right" state depending on the result
!                 ---------------------------------------------------------------
                  if ( ML .lt. -1.0_RP ) then      ! Supersonic inflow
                     edge % uB(iXi , 1:NEC)  = self % q
               
                  elseif ( ML .le. 0.0_RP ) then   ! Subsonic inflow
                     alpha = rhoInfty - pInfty / (aL*aL)
                     beta  = vnL - pL / (rhoL * aL) 

                     a1 = 0.5_RP/(rhoL*rhoL * aL*aL*aL*aL);
                     b1 = beta / (rhoL * aL*aL*aL) + alpha / (2*rhoL*rhoL*aL*aL);
                     c1 = beta * beta / ( 2 * aL * aL) + alpha * beta / ( rhoL * aL ) + Dimensionless % cp - self % Ht / ( aL * aL);
                     d1= alpha * beta * beta / 2 + self % Ht * pInfty / ( aL * aL ) - rhoInfty * self % Ht;

                     a2 = b1 / a1
                     b2 = c1 / a1
                     c2 = d1 / a1
                     
                     pR = ThirdDegreeRoots(a2,b2,c2)

                     rhoR = rhoInfty + (pR - pInfty) / (aL*aL)
                     vnR = vnL + (pR - pInfty) / (rhoL*aL)
                     uR = abs(vnR) * nInfty(IX)
                     vR = abs(vnR) * nInfty(IY)

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = edge % uB(iXi , IRHO) * uR
                     edge % uB(iXi , IRHOV) = edge % uB(iXi , IRHO) * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )
                     
                  elseif ( ML .lt. 1.0_RP ) then   ! Subsonic outflow
                     pR = pInfty
                     uR = uL
                     vR = vL
                     rhoR = pR / TL

                     edge % uB(iXi , IRHO) = rhoR
                     edge % uB(iXi , IRHOU) = rhoR * uR
                     edge % uB(iXi , IRHOV) = rhoR * vR
                     edge % uB(iXi , IRHOE) = cv * pR + 0.5_RP * (edge % uB(iXi,IRHOU) * uR + edge % uB(iXi,IRHOV) * vR )

                  elseif ( ML .ge. 1.0_RP ) then   ! Supersonic outflow
                     edge % uB(iXi , 1:NEC) = edge % Q(iXi, 1:NEC , 1)
      
                  end if
               
      
               end do 

            class default
         end select



         end associate

      end subroutine PressureInletBC_Update

!
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

      subroutine FarfieldBC_Describe( self )
         use Physics
         implicit none
         class(FarfieldBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Farfield."
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine FarfieldBC_Describe

      subroutine PressureOutletBC_Describe( self )
         use Physics
         implicit none
         class(PressureOutletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Pressure outlet."
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine PressureOutletBC_Describe

      subroutine PressureInletBC_Describe( self )
         use Physics
         implicit none
         class(PressureInletBC_t)                :: self

         associate ( gm1 => Thermodynamics % gm1 )
         write(STD_OUT , '(30X,A,A25,A)') "-> " , "Boundary condition type: " , "Pressure inlet."
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Pressure: " , gm1*(self % q(IRHOE) - 0.5_RP * ( self % q(IRHOU)**2.0_RP + self % q(IRHOV)**2.0_RP) / self % q(IRHO) ) * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total pressure: " , self % pt * refValues % p 
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total temperature: " , self % Tt * refValues % T
         write(STD_OUT , '(30X,A,A25,F10.2)') "-> " , "Total enthalpy: " , self % Ht * refValues % p
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Density: " , self % q(IRHO) * refValues % rho 
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "X-Velocity: " , self % q(IRHOU) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,F10.4)') "-> " , "Y-Velocity: " , self % q(IRHOV) / self % q(IRHO) * refValues % a
         write(STD_OUT , '(30X,A,A25,I10)') "-> " , "Angle of attack: " , nint(self % AngleOfAttack * 180.0_RP / PI)
         end associate

      end subroutine PressureInletBC_Describe

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
