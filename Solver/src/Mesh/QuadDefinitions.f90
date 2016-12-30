module QuadMeshDefinitions
   use SMConstants
   implicit none  


   integer, parameter         :: POINTS_PER_QUAD = 4
   integer, parameter         :: POINTS_PER_EDGE = 2
   integer, parameter         :: EDGES_PER_QUAD  = 4
   integer, parameter         :: QUADS_PER_EDGE  = 2

   integer, parameter         :: EBOTTOM = 1
   integer, parameter         :: ERIGHT = 2
   integer, parameter         :: ETOP = 3
   integer, parameter         :: ELEFT = 4


end module QuadMeshDefinitions
