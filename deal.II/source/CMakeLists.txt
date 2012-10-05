#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Compile the deal.II library
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(distributed)
ADD_SUBDIRECTORY(dofs)
ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(grid)
ADD_SUBDIRECTORY(hp)
ADD_SUBDIRECTORY(lac)
ADD_SUBDIRECTORY(multigrid)
ADD_SUBDIRECTORY(numerics)


FOREACH(build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${build} build_lowercase)

  #
  # Combine all ${build} OBJECT targets to a ${build} library:
  #
  FILE(STRINGS
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${build_lowercase}
    deal_ii_objects_${build_lowercase}
    )
  ADD_LIBRARY(${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    ${deal_ii_objects_${build_lowercase}}
    )

  SET_TARGET_PROPERTIES(${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    PROPERTIES
    VERSION ${VERSION}
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION ${VERSION}
    LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_${build}}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${build}}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${build}}"
    )

  TARGET_LINK_LIBRARIES(${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
    )

  INSTALL(TARGETS ${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    COMPONENT library
    )

ENDFOREACH()

INSTALL(EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

