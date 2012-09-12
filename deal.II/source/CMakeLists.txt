#
# Compile the deal.II library
#

#
# Each subdirectory specifies an OBJECT target:
#
ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(distributed)
ADD_SUBDIRECTORY(dofs)
ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(grid)
ADD_SUBDIRECTORY(hp)
ADD_SUBDIRECTORY(lac)
ADD_SUBDIRECTORY(multigrid)
ADD_SUBDIRECTORY(numerics)

#
# Combine all OBJECT targets to a shared library:
#

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

 # for contrib stuff:
 ${deal_ii_additional_object_files}
 )

#
# Link all external libraries into the target:
# TODO: toggle for as-needed
#

IF (CMAKE_BUILD_TYPE MATCHES "debug")
  TARGET_LINK_LIBRARIES(deal_II ${deal_ii_required_linker_flags} ${deal_ii_external_debug_libraries})
ELSE()
  TARGET_LINK_LIBRARIES(deal_II ${deal_ii_required_linker_flags} ${deal_ii_external_libraries})
ENDIF()
