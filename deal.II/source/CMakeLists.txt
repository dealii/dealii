ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(distributed)
ADD_SUBDIRECTORY(dofs)
ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(grid)
ADD_SUBDIRECTORY(hp)
ADD_SUBDIRECTORY(lac)
ADD_SUBDIRECTORY(multigrid)
ADD_SUBDIRECTORY(numerics)

ADD_LIBRARY(deal_II SHARED
 $<TARGET_OBJECTS:obj_base>
 $<TARGET_OBJECTS:obj_distributed>
 $<TARGET_OBJECTS:obj_dofs>
 $<TARGET_OBJECTS:obj_fe>
 $<TARGET_OBJECTS:obj_grid>
 $<TARGET_OBJECTS:obj_hp>
 $<TARGET_OBJECTS:obj_lac>
 $<TARGET_OBJECTS:obj_multigrid>
 $<TARGET_OBJECTS:obj_numerics>
 ${deal_ii_additional_object_files}
 )

 message( "${deal_ii_additional_object_files}" )
