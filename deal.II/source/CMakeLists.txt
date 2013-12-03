## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

MESSAGE(STATUS "Setup library")

#
# Compile the deal.II library
#

INCLUDE_DIRECTORIES(
  BEFORE # Ensure deal.II include directories come first
  #
  # Reverse order due to BEFORE:
  #
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

#
# List the directories where we have source files. the ones with the longest
# compile jobs come first so that 'make -j N' saturates many processors also
# towards the end of compiling rather than having to wait for one long
# compilation that, because it has been listed last, is started towards the
# end of everything (e.g. numerics/vectors.cc takes several minutes to
# compile...)
#
ADD_SUBDIRECTORY(numerics)
ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(dofs)
ADD_SUBDIRECTORY(lac)
ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(grid)
ADD_SUBDIRECTORY(hp)
ADD_SUBDIRECTORY(multigrid)
ADD_SUBDIRECTORY(distributed)
ADD_SUBDIRECTORY(algorithms)
ADD_SUBDIRECTORY(integrators)
ADD_SUBDIRECTORY(matrix_free)
ADD_SUBDIRECTORY(meshworker)

# A custom target to build the library:
ADD_CUSTOM_TARGET(build_library)

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
    dummy.cc # Workaround for a bug in the Xcode generator
    ${deal_ii_objects_${build_lowercase}}
    )

  ADD_DEPENDENCIES(build_library ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})
  ADD_DEPENDENCIES(library ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})

  SET_TARGET_PROPERTIES(${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    PROPERTIES
    VERSION "${DEAL_II_PACKAGE_VERSION}"
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION "${DEAL_II_PACKAGE_VERSION}"
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${build}}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${build}}"
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    )

  TARGET_LINK_LIBRARIES(${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
    ${DEAL_II_EXTERNAL_LIBRARIES}
    )

  FILE(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}
    )
  EXPORT(TARGETS ${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )

  INSTALL(TARGETS ${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    COMPONENT library
    )
ENDFOREACH()

INSTALL(EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

MESSAGE(STATUS "Setup library - Done")
