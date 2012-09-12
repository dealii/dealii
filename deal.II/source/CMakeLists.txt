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
# Sonaming: Well... we just use the version number. No point to wrack
# one's brain over the question whether a C++ library is still ABI
# backwards compatible.
#
SET_TARGET_PROPERTIES(deal_II PROPERTIES
  VERSION ${VERSION}
  SOVERSION ${VERSION}
  )

TARGET_LINK_LIBRARIES(deal_II ${deal_ii_required_linker_flags} ${deal_ii_external_libraries}) # TODO

