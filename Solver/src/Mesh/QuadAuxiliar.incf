
!
!        **********************************************************************************
!                 Auxiliar subroutines
!        **********************************************************************************
!
         subroutine searchEdge( nodesEl , nodesEdge , edgePosition , quadPosition , edgeDirection )
            implicit none
            integer, intent(in)        :: nodesEl(:)
            integer, intent(in)        :: nodesEdge(:)
            integer, intent(out)       :: edgePosition
            integer, intent(out)       :: quadPosition
            integer, intent(out)       :: edgeDirection
!           -----------------------------------------------------------
            integer                    :: currentEdge(POINTS_PER_EDGE)
            integer                    :: edge
            logical                    :: edges_are_equal

            do edge = 1 , EDGES_PER_QUAD
!
!              Obtain current edge of the element
!              ----------------------------------
               currentEdge(1) = nodesEl(edge)
               if (edge .eq. EDGES_PER_QUAD) then
                  currentEdge(2) = nodesEl(1)
               else
                  currentEdge(2) = nodesEl(edge+1)
               end if
!
!              Compare the element edge with the original edge
!              -----------------------------------------------
               call compareEdges( edge1 = currentEdge , edge2 = nodesEdge , edges_are_equal = edges_are_equal , edgeDirection = edgeDirection )

               if (edges_are_equal) then
                  edgePosition = edge
                  if ( edgeDirection .eq. FORWARD ) then
                     quadPosition = LEFT
                  else
                     quadPosition = RIGHT
                  end if

                  return

               end if

            end do   

         end subroutine searchEdge

         subroutine compareEdges ( edge1 , edge2 , edges_are_equal , edgeDirection )
            implicit none
            integer, intent(in)        :: edge1(:)
            integer, intent(in)        :: edge2(:)
            logical, intent(out)       :: edges_are_equal
            integer, intent(out)       :: edgeDirection

            if ((edge1(1) .eq. edge2(1)) .and. (edge1(2) .eq. edge2(2))) then
               edges_are_equal = .true.
               edgeDirection = FORWARD
            
            elseif ((edge1(1) .eq. edge2(2)) .and. (edge1(2) .eq. edge2(1))) then
               edges_are_equal = .true.
               edgeDirection = BACKWARD

            else
               edges_are_equal = .false.
               edgeDirection = 0
            end if

         end subroutine compareEdges


