#
# Compile the deal.II library
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

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
# Combine all OBJECT targets to a library:
#
ADD_LIBRARY(deal_II
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
# Sonaming: Well... we just use the version number.
# No point to wrack one's brain over the question whether a new version of
# a C++ library is still ABI backwards compatible :-]
#
SET_TARGET_PROPERTIES(deal_II PROPERTIES
  VERSION ${VERSION}
  SOVERSION ${VERSION}
  DEBUG_POSTFIX ".g"
  )

SET_TARGET_PROPERTIES(deal_II PROPERTIES
  LINK_FLAGS "${deal_ii_required_linker_flags}"
  )

TARGET_LINK_LIBRARIES(deal_II ${deal_ii_external_libraries})


#
# Finally, add a rule how to install the library:
#
INSTALL(TARGETS deal_II
  LIBRARY
  DESTINATION ${DEAL_II_LIBRARY_RELDIR}
  COMPONENT library
  )


IF(DEAL_II_INSTALL_COMPAT_FILES)
  #
  # If ${LIB_SUFFIX} is not empty, create a symlink from lib to
  # lib${LIB_SUFFIX} to ensure compatibility with the old build system:
  #
  IF(NOT "${LIB_SUFFIX}" STREQUAL "")
    INSTALL(
      CODE
        "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
          ${CMAKE_INSTALL_PREFIX}/lib
          )"
      COMPONENT compat_files
      )
  ENDIF()
ENDIF()

