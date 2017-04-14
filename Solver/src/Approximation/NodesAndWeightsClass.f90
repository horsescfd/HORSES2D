module nodesAndWeights_class
    use Physics
    use SMConstants
    use InterpolationAndDerivatives
    use LegendreAlgorithms
    implicit none

#include "Defines.h" 

    private
    public NodalStorage , newNodalStorage
    public NodesAndWeights_t
!  
!   **************************************************************
!   NodesAndWeights entry type:
!       Contains all the data involving the spatial discretisation
!   **************************************************************
!
    type  NodesAndWeights_t
        integer                   :: N
        integer                   :: nodes
        real(kind=RP), pointer    :: xi(:)
        real(kind=RP), pointer    :: w(:)
        real(kind=RP), pointer    :: wb(:)
        real(kind=RP), pointer    :: M(:,:)
        real(kind=RP), pointer    :: Minv(:,:)
        real(kind=RP), pointer    :: invM2D(:,:)
        real(kind=RP), pointer    :: D(:,:)
        real(kind=RP), pointer    :: DT(:,:)
        real(kind=RP), pointer    :: MD(:,:)
        real(kind=RP), pointer    :: hatD(:,:)
        real(kind=RP), pointer    :: trMD(:,:)
        real(kind=RP), pointer    :: tildeMTD(:,:)
        real(kind=RP), pointer    :: lb(:,:)
        real(kind=RP), pointer    :: lbw(:,:)
        real(kind=RP), pointer    :: T(:,:)     ! Interpolation matrix to the TAIL of the linked list 
!                                                  (integration points are placed there)
        class(NodesAndWeights_t), pointer :: next => NULL()
        contains
            procedure            :: lj   => polyevaluation
            procedure            :: dlj  => polyDerivativeEvaluation
            procedure            :: init => initWithDegree
    end type NodesAndWeights_t
!  
!   **************************************************************
!   NodalStorage  type:
!       LinkedList that contains all the stored entries
!   **************************************************************
!
    type NodalStorage
        class(NodesAndWeights_t), pointer   :: head => NULL()
        class(NodesAndWeights_t), pointer   :: tail => NULL()
        integer                     :: nrecords = 0
        contains
            procedure :: add => addNodesAndWeights
            procedure :: find => findEntryWithDegree
            procedure :: list => listNodesAndWeights
            procedure :: computeInterpolationMatrices
    end type NodalStorage

!
!   ========
    contains
!   ========
!
!
!       **************************
!       Initialization subroutines
!       **************************
!
        function newNodalStorage()
            implicit none
            type(NodalStorage)     :: newNodalStorage

            newNodalStorage % head => NULL()
            newNodalStorage % tail => NULL()
            newNodalStorage % nrecords = 0
        
        end function newNodalStorage

!
!       *****************************
!       Memory management subroutines
!       *****************************
!
        subroutine addNodesAndWeights(self,N,nodes,NAW)
            implicit none
            class(NodalStorage)        :: self
            integer                            :: N
            integer                            :: nodes
            class(NodesAndWeights_t), pointer, optional :: NAW
            integer                            :: i
            class(NodesAndWeights_t), pointer  :: current


            if (self % nrecords .GT. 0) then

                current => self % head

                do i = 1 , self % nrecords

                    if ( (N .EQ. current % N).AND.(nodes .EQ. current % nodes) ) then

                        if (present(NAW)) NAW => current
                        return
                    end if
            
                    current => current % next
         
                end do
    
                if (associated(current)) then
                    write(STD_OUT , * )     "Error in the linked list"
                    STOP "Stopped."
                end if

                current => NULL()

            end if 

!
!           Code arrives here just if the requested degree/nodes is not found
!
            if ( .NOT. associated(self % head) ) then
!
!               Initiate the list
!   
                allocate(self % head)

                self % tail => self % head

                self % nrecords = 1

            else
!           
!               Append a new entry to the list
!   
                allocate( self % tail % next )
                self % tail => self % tail % next

                self % nrecords = self % nrecords + 1
                

            end if

            call self % tail % init(N,nodes)
          
            if (present(NAW)) NAW => self % tail


        end subroutine addNodesAndWeights

        subroutine initWithDegree(self , N , nodes)
            implicit none
            class(NodesAndWeights_t)        :: self
            integer                         :: N
            integer                         :: nodes

            self % N = N
            self % nodes = nodes

            allocate ( self % xi     ( 0:N     )  ) 
            allocate ( self % w      ( 0:N     )  ) 
            allocate ( self % wb     ( 0:N     )  ) 
            allocate ( self % M      ( 0:N,0:N )  ) 
            allocate ( self % Minv   ( 0:N,0:N )  ) 
            allocate ( self % invM2D ( 0:N,0:N )  ) 
            allocate ( self % D      ( 0:N,0:N )  ) 
            allocate ( self % DT     ( 0:N,0:N )  ) 
            allocate ( self % hatD   ( 0:N,0:N )  ) 
            allocate ( self % MD     ( 0:N,0:N )  ) 
            allocate ( self % trMD   ( 0:N,0:N )  ) 
            allocate ( self % lb     ( 0:N,2   )  ) 
            allocate ( self % lbw    ( 0:N,2   )  ) 

            self % T => NULL()
            self % tildeMTD => NULL()


            call computeNodesAndWeights(self)
                

        end subroutine initWithDegree

        function findEntryWithDegree(self,N,nodes)
            implicit none
            class(NodalStorage)        :: self
            integer                         :: N
            integer                         :: nodes
            class(NodesAndWeights_t), pointer       :: findEntryWithDegree
            integer                         :: i

            if (self % nrecords .EQ. 0) THEN
!
!               The list contains no items
!   
                findEntryWithDegree => NULL()

            else
                
                findEntryWithDegree => self % head
                do i = 1 , self % nrecords
    
                    if((findEntryWithDegree % N .EQ. N).AND.(findEntryWithDegree % nodes .EQ. nodes)) then
                        return
                    end if

                    findEntryWithDegree => findEntryWithDegree % next

                end do
            
            end if

!          
!           If the code arrives here, the desired N is not found
!
            findEntryWithDegree => NULL()

        end function findEntryWithDegree

!
!       *********************
!       Computing subroutines
!       *********************
!
        subroutine computeNodesAndWeights(self)
            implicit none
            class(NodesAndWeights_t)        :: self
            integer                         :: N
            integer                         :: i
        
            N = self % N

!
!           -------------------------
!           Compute nodes and weights
!           -------------------------
!
            if (self % nodes .EQ. LG) then
!               Legendre-Gauss nodes

                call GaussLegendreNodesAndWeights( N, self %  xi, self % w )
            
            elseif (self % nodes .EQ. LGL) then
!               Legendre-Gauss-Lobatto nodes

                call LegendreGaussLobattoNodesAndWeights( N, self % xi, self % w )

            end if
!                
!           -------------------------------
!           Adation to (0,1) domain
!           -------------------------------
!
            self % xi = 0.5_RP + 0.5_RP * self % xi
            self % w  = 0.5_RP * self % w
            
!
!           ---------------------------
!           Compute barycentric weights
!           ---------------------------
!
            call BarycentricWeights( N, self % xi, self % wb )
!
!           -----------------------------------
!           Compute mass matrix and its inverse
!           -----------------------------------
!
            self % M = 0.0_RP
            self % Minv = 0.0_RP

            do i = 0 , self % N
                 self % M(i,i)    =          self % w(i)
                 self % Minv(i,i) = 1.0_RP / self % w(i)
                 self % invM2D(i,:)  = 1.0_RP / (self % w(i) * self % w)
            end do 
!
!           -----------------------------------------
!           Compute the first order derivative matrix
!           -----------------------------------------
!
            call PolynomialDerivativeMatrix( N , self % xi , self % D )
            self % DT   = transpose( self % D )
            self % MD   = matmul( self % M , self % D )
            self % trMD = transpose( self % MD )
            self % hatD = matmul( self % MD , self % Minv )

!       
!           --------------------------------------------
!           Compute the boundary lagrange interpolations
!           --------------------------------------------
!
            call LagrangeInterpolatingPolynomialBarycentric(  0.0_RP, N, self % xi, self % wb, self % lb(: ,LEFT) )
            call LagrangeInterpolatingPolynomialBarycentric(  1.0_RP, N, self % xi, self % wb, self % lb(:,RIGHT) )
   
            self % lbw ( :,LEFT  ) = self % lb ( :,LEFT  ) / self % w
            self % lbw ( :,RIGHT ) = self % lb ( :,RIGHT ) / self % w
            

        end subroutine computeNodesAndWeights
            
        pure function polyevaluation(self,x)   result(val)
            implicit none
            class(NodesAndWeights_t), intent(in)   :: self
            real(kind=RP)           , intent(in)   :: x
            real(kind=RP)                          :: val( 0 : self % N )
    
            call LagrangeInterpolatingPolynomialBarycentric(x , self % N , self % xi , self % wb , val)

        end function polyevaluation

        pure function polyDerivativeEvaluation( self , x ) result(val)
            implicit none
            class(NodesAndWeights_t), intent(in) :: self
            real(kind=RP), intent(in)            :: x
            real(kind=RP)                        :: val( 0 : self % N )
!           -------------------------------------------------------------
            integer                          :: i

            do i = 0 , self % N
               val(i) = EvaluateLagrangePolyDerivative( i, x, self % N , self % xi)
            end do

        end function polyDerivativeEvaluation

        subroutine computeInterpolationMatrices( self , spI )
!           
!           ********************************************
!              If the code enters here, is because
!              an over-integration mode is 
!              considered. There are several things 
!              performed here:
!                 -> The interpolation matrix Tnj=lj(xin)
!                 -> The derivative matrix tilde(M)TD
!                 -> The mass matrix tr(T)tilde(M)T
!                 -> The inverse of the mass matrix
!           ********************************************
!  
            use InterpolationAndDerivatives
            use MatrixOperations
            implicit none
            class(NodalStorage)         :: self
            class(NodesAndWeights_t)      :: spI
            class(NodesAndWeights_t), pointer   :: current

            current => self % head
            do 
!              Compute interpolation matrix
               allocate( current % T(0: spI % N ,0: current % N ) )
               call PolynomialInterpolationMatrix( current % N , spI % N , current % xi , current % wb , spI % xi , current % T)

!              Compute the derivative matrix tilde(M) T D
               allocate( current % tildeMTD( 0 : spI % N , 0 : current % N ) )
               call TripleMatrixProduct(spI % M , current % T , current % D , current % tildeMTD )

!              Compute the mass matrix M
               call InnerProduct( current % T , spI % M , size(current % T,1) , size(current % T,2) , current % M )

!              Compute the inverse of the mass matrix
               current % Minv = inv( current % M )
               
!              
               current => current % next
               if (.not. associated(current)) exit

            end do 
         end subroutine computeInterpolationMatrices


!
!       ***************
!       I/O subroutines
!       ***************
!
        subroutine listNodesAndWeights(self,unit)
            implicit none
            class(NodalStorage)         :: self
            integer                             :: unit
            integer                             :: i
            class(NodesAndWeights_t), pointer   :: current
            
            write(unit , '(/)' ) 
            write(unit , '(10X , A)') "* Listing all entries added to the Nodes And Weights list:"
            write(unit , '(/)' ) 

            if (self % nrecords .GT. 0) then
                current => self % head
                do i = 1 , self % nrecords
                    write(unit , '(20X,A,I0,A)') "Position ",i,": " 
                    write(unit , '(40X,A,A10,I0,A)') "->" ,"N: ", current % N,"."
                    if (current % nodes .EQ. LG) then
                        write(unit , '(40X,A,A10,A)') "->","Nodes: ", "Legendre-Gauss."
                    else if (current % nodes .EQ. LGL) then
                        write(unit , '(40X,A,A10,A)') "->", "Nodes: ", "Legendre-Gauss-Lobatto."
                    end if
                    write(unit , '(18X,A)') "_________________________________________________________________"
                    write(unit , * ) 
                    current => current % next
                end do

                current => NULL()
            end if

        end subroutine listNodesAndWeights





end module nodesAndWeights_class
