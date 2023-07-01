## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2023 by the deal.II authors
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

#
# Automatically add the current source and binary directory as include
# directory. The binary directory is needed in order to find the expanded
# .inst include files.
#
set(CMAKE_INCLUDE_CURRENT_DIR true)

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
add_subdirectory(matrix_free)
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

#
# Define a "dealii", i.e., "${DEAL_II_TARGET_NAME}", target that contains
# DEAL_II_CXX_FLAGS[_DEBUG|_RELEASE] and otherwise aliases dealii_debug or
# dealii_release depending on the downstream CMAKE_BUILD_TYPE.
#

add_library(${DEAL_II_TARGET_NAME} INTERFACE)
target_compile_flags(${DEAL_II_TARGET_NAME} INTERFACE
  "$<COMPILE_LANGUAGE:CXX>" "${DEAL_II_CXX_FLAGS}"
  )
target_link_flags(${DEAL_II_TARGET_NAME} INTERFACE "${DEAL_II_LINKER_FLAGS}")

foreach(build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${build} build_lowercase)
  if("${build}" MATCHES "DEBUG")
    set(build_camelcase "Debug")
  elseif("${build}" MATCHES "RELEASE")
    set(build_camelcase "Release")
  endif()

  #
  # Combine all ${build} OBJECT targets to a ${build} library:
  #

  set(_object_files "")
  if(NOT BUILD_SHARED_LIBS)
    get_property(_object_targets GLOBAL PROPERTY DEAL_II_OBJECT_TARGETS_${build})
    foreach(_target ${_object_targets})
      list(APPEND _object_files "$<TARGET_OBJECTS:${_target}>")
    endforeach()
  endif()

  add_library(${DEAL_II_TARGET_NAME}_${build_lowercase} dummy.cc ${_object_files})
  add_dependencies(library ${DEAL_II_TARGET_NAME}_${build_lowercase})

  #
  # Add compile flags, definitions and (public) feature (recorded in
  # DEAL_II_TARGETS(|_DEBUG|_RELEASE)).
  #
  populate_target_properties(${DEAL_II_TARGET_NAME}_${build_lowercase} ${build})

  #
  # Record the expected C++ standard as a compile feature. This target
  # property ensures that support for our expected C++ standard is always
  # enabled in client user code irrespective of what compile flags/options
  # they have set.
  #
  target_compile_features(${DEAL_II_TARGET_NAME}_${build_lowercase}
    INTERFACE cxx_std_${CMAKE_CXX_STANDARD}
    )

  if(BUILD_SHARED_LIBS)
    #
    # Add all object targets as private link targets
    #
    get_property(_object_targets GLOBAL PROPERTY DEAL_II_OBJECT_TARGETS_${build})
    target_link_libraries(${DEAL_II_TARGET_NAME}_${build_lowercase}
      PRIVATE ${_object_targets}
      )
  endif()

  set_target_properties(${DEAL_II_TARGET_NAME}_${build_lowercase}
    PROPERTIES
    LINKER_LANGUAGE "CXX"
    VERSION "${DEAL_II_PACKAGE_VERSION}"
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION "${DEAL_II_PACKAGE_VERSION}"
    ARCHIVE_OUTPUT_NAME "${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}"
    LIBRARY_OUTPUT_NAME "${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}"
    RUNTIME_OUTPUT_NAME "${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    )

  target_compile_flags(${DEAL_II_TARGET_NAME} INTERFACE
    "$<AND:$<CONFIG:${build_camelcase}>,$<COMPILE_LANGUAGE:CXX>>"
    "${DEAL_II_CXX_FLAGS_${build}}"
    )
  target_link_flags(${DEAL_II_TARGET_NAME} INTERFACE
    "$<CONFIG:${build_camelcase}>" "${DEAL_II_LINKER_FLAGS_${build}}"
    )

  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set_target_properties(${DEAL_II_TARGET_NAME}_${build_lowercase}
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
    set_target_properties(${DEAL_II_TARGET_NAME}_${build_lowercase}
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
      ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
      )
  endif()

  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR})
  export(TARGETS ${DEAL_II_TARGET_NAME}_${build_lowercase}
    NAMESPACE "${DEAL_II_TARGET_NAME}::"
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )

  install(TARGETS ${DEAL_II_TARGET_NAME}_${build_lowercase}
    COMPONENT library
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    )

  target_link_libraries(${DEAL_II_TARGET_NAME} INTERFACE
    "$<$<CONFIG:${build_camelcase}>:${DEAL_II_TARGET_NAME}::${DEAL_II_TARGET_NAME}_${build_lowercase}>"
    )
endforeach()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR})
export(TARGETS ${DEAL_II_TARGET_NAME}
  NAMESPACE "${DEAL_II_TARGET_NAME}::"
  FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
  APPEND
)

install(TARGETS ${DEAL_II_TARGET_NAME}
  COMPONENT library
  EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
  RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
  ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
)

install(EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
  NAMESPACE "${DEAL_II_TARGET_NAME}::"
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

message(STATUS "Setting up library - Done")
