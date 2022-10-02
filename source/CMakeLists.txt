## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

message(STATUS "Setting up library")

#
# This happens surprisingly often:
#
if("${DEAL_II_PACKAGE_NAME}" STREQUAL  "")
  message(FATAL_ERROR
    "\nWrong source directory specified.\n"
    "You pointed cmake to the subdirectory \"source\" instead of the base "
    "directory (the one containing the file \"README.md\"). Make sure that "
    "the path in the cmake invocation does not end in \".../source\". If you "
    "have created a \"build\" subdirectory, the cmake invocation should simply "
    "be of the form:\n"
    "  $ cmake  [...options...]  ../\n"
    "(path not ending in \".../source\")"
    )
endif()

#
# Compile the deal.II library
#

include_directories(
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  SYSTEM
  ${DEAL_II_INCLUDE_DIRS}
  )

#
# List the directories where we have source files. the ones with the longest
# compile jobs come first so that 'make -j N' saturates many processors also
# towards the end of compiling rather than having to wait for one long
# compilation that, because it has been listed last, is started towards the
# end of everything (e.g. numerics/vectors.cc takes several minutes to
# compile...)
#
add_subdirectory(numerics)
add_subdirectory(fe)
add_subdirectory(dofs)
add_subdirectory(lac)
add_subdirectory(base)
add_subdirectory(cgal)
add_subdirectory(gmsh)
add_subdirectory(grid)
add_subdirectory(hp)
add_subdirectory(multigrid)
add_subdirectory(distributed)
add_subdirectory(algorithms)
add_subdirectory(integrators)
add_subdirectory(matrix_free)
add_subdirectory(meshworker)
add_subdirectory(opencascade)
add_subdirectory(particles)
add_subdirectory(differentiation)
add_subdirectory(physics)
add_subdirectory(optimization/rol)
add_subdirectory(non_matching)
add_subdirectory(sundials)
add_subdirectory(trilinos)
add_subdirectory(arborx)

foreach(build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${build} build_lowercase)

  #
  # Combine all ${build} OBJECT targets to a ${build} library:
  #

  get_property(_objects GLOBAL PROPERTY DEAL_II_OBJECTS_${build})

  add_library(${DEAL_II_NAMESPACE}_${build_lowercase}
    dummy.cc # Workaround for a bug in the Xcode generator
    ${_objects}
    )
  add_dependencies(library ${DEAL_II_NAMESPACE}_${build_lowercase})

  set_target_properties(${DEAL_II_NAMESPACE}_${build_lowercase}
    PROPERTIES
    VERSION "${DEAL_II_PACKAGE_VERSION}"
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION "${DEAL_II_PACKAGE_VERSION}"
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}"
    LINKER_LANGUAGE "CXX"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${build}}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}}"
    ARCHIVE_OUTPUT_NAME "${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}"
    LIBRARY_OUTPUT_NAME "${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}"
    RUNTIME_OUTPUT_NAME "${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    )

  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set_target_properties(${DEAL_II_NAMESPACE}_${build_lowercase}
      PROPERTIES
      MACOSX_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}"
      )
  endif()

  # Under Windows (MSVC) cmake will always generate multi-configuration
  # projects. When building on the command line with 'cmake --build .',
  # release and debug builds of the library are done with the default 'Debug'
  # configuration. This causes the debug and release .lib to be built inside
  # ./lib/Debug/. This is not very pretty and confuses example/test projects,
  # so we just hard-wire the location here. We only really need to set static
  # lib locations for _DEBUG (no support for dynamic linking, _RELEASE will be
  # ignored), but we do it anyhow.
  if (DEAL_II_MSVC)
    set_target_properties(${DEAL_II_NAMESPACE}_${build_lowercase}
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
      ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
      )
  endif()


  target_link_libraries(${DEAL_II_NAMESPACE}_${build_lowercase}
    ${DEAL_II_LIBRARIES_${build}}
    ${DEAL_II_LIBRARIES}
    )

  file(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}
    )
  export(TARGETS ${DEAL_II_NAMESPACE}_${build_lowercase}
    NAMESPACE "${DEAL_II_NAMESPACE}::"
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )

  install(TARGETS ${DEAL_II_NAMESPACE}_${build_lowercase}
    COMPONENT library
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    )
endforeach()

install(EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
  NAMESPACE "${DEAL_II_NAMESPACE}::"
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

message(STATUS "Setting up library - Done")
