## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
# List the directories where we have source files. the ones with the longest
# compile jobs come first so that 'make -j N' saturates many processors also
# towards the end of compiling rather than having to wait for one long
# compilation that, because it has been listed last, is started towards the
# end of everything (e.g. some of the numerics/vector_tools_*.cc files can
# take quite a long time to compile).
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
# Set up and configure the two library targets dealii_debug and
# dealii_release, i.e., "${DEAL_II_TARGET_NAME}_${build}."
#
# Both library targets are composed of the object and bundled object
# targets that have been set up in the various subdirectories under source/
# and bundled/ and that we recorded in the DEAL_II_OBJECT_TARGETS_${build}
# global property.
#

foreach(build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${build} build_lowercase)

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

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.30")
    #
    # Create alias targets "dealii::dealii_release" and "dealii::dealii_debug"
    # to have the exported target names available when populating the link
    # interface for our "dealii::dealii" convenience target:
    #
    add_library(${DEAL_II_TARGET_NAME}::${DEAL_II_TARGET_NAME}_${build_lowercase}
      ALIAS ${DEAL_II_TARGET_NAME}_${build_lowercase}
      )
  endif()

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.28")
    get_property(_header_files GLOBAL PROPERTY DEAL_II_HEADER_FILES)
    target_sources(${DEAL_II_TARGET_NAME}_${build_lowercase}
      PUBLIC
      FILE_SET headers
      TYPE HEADERS
      FILES
        ${_header_files}
        ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
        ${CMAKE_BINARY_DIR}/include/deal.II/base/revision.h
      BASE_DIRS
        ${CMAKE_SOURCE_DIR}/include/
        ${CMAKE_BINARY_DIR}/include/
      )
  endif()
endforeach()

#
# Now, define a "dealii", i.e., "${DEAL_II_TARGET_NAME}", target that can
# be used by user projects to automatically switch between the release and
# debug versions of the library depending on the build type.
#
# For convenience, we also add our compiler and linker flags to the target.
# This is useful for user projects that simply want to use our compiler and
# linker flags and not set up their own.
#
# If a downstream project does not want us to override their compiler and
# linker flags, then they will need to use the "dealii_debug" and
# "dealii_release" targets directly.
#

add_library(${DEAL_II_TARGET_NAME} INTERFACE)
add_dependencies(library ${DEAL_II_TARGET_NAME})

#
# Record the libraries, as well as our compiler and linker flags as a
# *public* interface property.
#
# This way all targets linking against dealii::dealii will be compiled and
# linked with our compiler and linker flags.
#
# We start with the configuration independent bits...
#

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
  # ... and then add debug and release configuration dependent flags
  # guarded with a "generator expression" [1]. Note that we further guard
  # the compiler flags with $<COMPILE_LANGUAGE:CXX> so that they are only
  # used for compiling C++ code.
  #

  target_compile_flags(${DEAL_II_TARGET_NAME} INTERFACE
    "$<AND:$<CONFIG:${build_camelcase}>,$<COMPILE_LANGUAGE:CXX>>"
    "${DEAL_II_CXX_FLAGS_${build}}"
    )
  target_link_flags(${DEAL_II_TARGET_NAME} INTERFACE
    "$<CONFIG:${build_camelcase}>" "${DEAL_II_LINKER_FLAGS_${build}}"
    )

  #
  # Attach the debug and release libraries:
  #

  target_link_libraries(${DEAL_II_TARGET_NAME} INTERFACE
    "$<$<CONFIG:${build_camelcase}>:${DEAL_II_TARGET_NAME}::${DEAL_II_TARGET_NAME}_${build_lowercase}>"
    )

endforeach()

message(STATUS "Setting up library - Done")
