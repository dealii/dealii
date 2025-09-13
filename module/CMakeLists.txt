## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


if(DEAL_II_WITH_CXX20_MODULE)
  message(STATUS "Setting up the C++20 module")

  # First version with proper module support
  cmake_minimum_required(VERSION 3.28)

  #
  # The scripts used below are all written in Python. Ensure that we have
  # it available.
  #

  find_package (Python3 COMPONENTS Interpreter)
  if (NOT Python3_FOUND)
    message(FATAL_ERROR
      "To build a C++20 module, you need to have the Python3 interpreter available"
      )
  endif()

  #
  # Create deal.II/macros.h:
  #
  # Modules cannot export macros. So we need to still provide a way to
  # export the macros that are part of deal.II, which includes
  # specifically the configuration stuff in config.h and the exception
  # stuff in exception_macros.h, along with macros we have cloned from
  # external libraries.
  #
  # To make things easier for users, we combine all of these into a
  # file deal.II/macros.h that we can rely on in all of the module
  # files. Define how to build this file first, so that we can let
  # other things depend on it later on.
  #

  set(_dealii_macros_header "${CMAKE_CURRENT_BINARY_DIR}/../include/deal.II/macros.h")
  add_custom_command(
    OUTPUT ${_dealii_macros_header}
    DEPENDS
      "${CMAKE_CURRENT_BINARY_DIR}/../include/deal.II/base/config.h"
      "${CMAKE_CURRENT_BINARY_DIR}/../include/deal.II/base/revision.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/../include/deal.II/base/exception_macros.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/deal.II/boost_macros.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/deal.II/petsc_macros.h"
    COMMAND
      ${CMAKE_COMMAND}
      ARGS -E cat
      "${CMAKE_CURRENT_BINARY_DIR}/../include/deal.II/base/config.h"
      "${CMAKE_CURRENT_BINARY_DIR}/../include/deal.II/base/revision.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/../include/deal.II/base/exception_macros.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/deal.II/boost_macros.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/deal.II/petsc_macros.h"
      ">"
      "${_dealii_macros_header}"
    )
  add_custom_target(build_macros_header DEPENDS ${_dealii_macros_header})


  ########################################################################
  #                                                                      #
  #                        interface module units:                       #
  #                                                                      #
  ########################################################################

  #
  # First, convert header files into *.ccm module partitions:
  #
  # Next, we need to set up interface module partitions for each of
  # deal.II's header files. We query the list of all header files from the
  # DEAL_II_HEADER_FILES property (which includes everything we picked up
  # from the include/ directory, but not config.h which resided in the
  # build directory). From the name of the header we construct the name of
  # an interface module partition file. These files, by convention, have
  # suffix .ccm (or .cppm).
  #

  set(_interface_module_partition_units)
  get_property(_header_files GLOBAL PROPERTY DEAL_II_HEADER_FILES)
  foreach (_header_file ${_header_files})

    # Exclude preprocessor #defines only headers already inlucded in macros.h:
    if (${_header_file} STREQUAL "${CMAKE_SOURCE_DIR}/include/base/exception_macros.h")
      continue()
    endif()

    # Check whether the file actually exports anything into namespace dealii.
    # Only those files should be part of the module, but we have some header
    # files that are only there for backward compatibility, or that provide
    # wrappers for external projects, and we need to exlude those.
    file(READ "${_header_file}" _header_file_as_string)
    string(FIND "${_header_file_as_string}" "DEAL_II_NAMESPACE_OPEN" _match)
    if(${_match} EQUAL -1)
      continue()
    endif()

    string(REGEX REPLACE
      ".*/deal.II/(.*)\.h$" "${CMAKE_CURRENT_BINARY_DIR}/interface_module_partitions/\\1.ccm"
      _module_file ${_header_file}
      )

    # Create an interface module partition file from the header file:
    add_custom_command(
      OUTPUT ${_module_file}
      DEPENDS
        "${_header_file}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/convert_header_file_to_interface_module_unit.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/convert_to_module_units_common.py"

      COMMAND ${Python3_EXECUTABLE}
      ARGS "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/convert_header_file_to_interface_module_unit.py"
           "${_header_file}"
           "${_module_file}"
      )

    # Compiling this file will require access to the deal.II/macros.h header:
    set_property(SOURCE _module_file APPEND PROPERTY
      OBJECT_DEPENDS ${_dealii_macros_header}
      )

    list(APPEND _interface_module_partition_units "${_module_file}")
  endforeach()

  #
  # Append interface units for external libraries:
  #
  # We have some other interface units that we need to add to the list,
  # namely the ones that were written by hand instead of created via
  # script above. We need to add these to the list of interface units
  # as well.
  #
  # Specifically: It is very inefficient in the module system to have
  # repeated #includes in many module partition files because when you
  # 'import' those partitions, you also have to load everything they
  # #included. In other words, you get the same content *many times*,
  # once from each imported partition, rather than only once via the
  # old-style #include system. We deal with this by wrapping all of
  # our external dependencies into partitions that we can 'import'
  # wherever we need.
  #
  set(_external_interface_module_partition_units
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/adolc.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/boost.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/cgal.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/hdf5.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/kokkos.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/metis.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/mumps.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/mpi.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/muparser.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/opencascade.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/p4est.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/petsc.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/slepc.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/std.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/sundials.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/taskflow.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/tbb.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/trilinos.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/umfpack.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/vtk.ccm"
    "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions/zlib.ccm"
    )

  #
  # And finally, generate the primary module interface unit.
  #
  # This unit (the "primary module interface unit") simply imports all of
  # the partitions and then re-exports them.
  #
  # Build this file via a script from the list of all of the other module
  # input files:
  #
  set(_primary_module_interface_unit
    "${CMAKE_CURRENT_BINARY_DIR}/interface_module_partitions/dealii.ccm"
    )
  add_custom_command(
    OUTPUT ${_primary_module_interface_unit}
    DEPENDS
      "${_interface_module_partition_units}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/build_primary_interface_unit.py"

    COMMAND ${Python3_EXECUTABLE}
    ARGS "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/build_primary_interface_unit.py"
         ${_interface_module_partition_units}
         ">" "${_primary_module_interface_unit}"
    )


  ########################################################################
  #                                                                      #
  #                   Implementation module partitions:                  #
  #                                                                      #
  ########################################################################

  #
  # Next, we also need to set up implementation module partitions for
  # each of deal.II's source files.
  #

  set(_implementation_module_partition_units)
  get_property(_source_files GLOBAL PROPERTY DEAL_II_SOURCE_FILES)
  foreach (_source_file ${_source_files})
    # Check whether the file actually exports anything into namespace
    # dealii. Unlike for header files, this realy shouldn't be the case, so
    # check here and error out if we find such a file -- that's because if
    # that happened, the script would produce an implementation unit that
    # does not actually have a partition and the compiler will error
    # anyway.
    file(READ "${_source_file}" _source_file_as_string)
    string(FIND "${_source_file_as_string}" "DEAL_II_NAMESPACE_OPEN" _match)
    if(${_match} EQUAL -1)
      message(FATAL_ERROR "${_source_file} does not appear to implement anything in namespace dealii.")
      continue()
    endif()

    # From the name of the source, construct the name of an implementation module
    # partition file:
    string(REGEX REPLACE
      ".*/source/(.*)\.cc$" "${CMAKE_CURRENT_BINARY_DIR}/implementation_module_partitions/\\1.cc"
      _module_file ${_source_file}
      )

    # Create an interface module partition file from the source file:
    add_custom_command(
      OUTPUT ${_module_file}
      DEPENDS
        "${_source_file}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/convert_source_file_to_implementation_module_unit.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/convert_to_module_units_common.py"

      COMMAND ${Python3_EXECUTABLE}
      ARGS "${CMAKE_CURRENT_SOURCE_DIR}/../contrib/utilities/convert_source_file_to_implementation_module_unit.py"
           "${_source_file}"
           "${_module_file}"
      )

    # Compiling this file will require access to the deal.II/macros.h header:
    set_property(SOURCE _module_file APPEND PROPERTY
      OBJECT_DEPENDS ${_dealii_macros_header}
      )

    list(APPEND _implementation_module_partition_units "${_module_file}")
  endforeach()

  #
  # Generate two top-level targets for compiling the implementation module
  # partitions into a shared library:
  #

  foreach(build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${build} build_lowercase)
    if("${build}" MATCHES "DEBUG")
      set(build_camelcase "Debug")
    elseif("${build}" MATCHES "RELEASE")
      set(build_camelcase "Release")
    endif()

    #
    # Define the library. Compile it with the usual flags we also use for
    # all other targets, but do set the DEAL_II_BUILDING_CXX20_MODULE
    # preprocessor variable for a small number of cases where we need to
    # work around things that work differently between module and
    # non-module builds.
    #

    add_library(${DEAL_II_TARGET_NAME}_module_${build_lowercase})

    populate_target_properties(${DEAL_II_TARGET_NAME}_module_${build_lowercase} ${build})
    target_compile_definitions(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PRIVATE "DEAL_II_BUILDING_CXX20_MODULE"
      )

    #
    # Say what the source files for the library are. These are all of
    # the interface and implementation partition units we have created
    # above (including the ones that wrap external libraries), plus
    # the primary interface unit.
    #

    target_sources(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PUBLIC
      FILE_SET module_interface_units
      TYPE CXX_MODULES
      FILES
        ${_interface_module_partition_units}
        ${_primary_module_interface_unit}
      BASE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/interface_module_partitions"
      )
    # Work around an issue with CMake 3.28.3 ... 4.0.3 (and possible later) where the
    # export() logic loses one of the BASE_DIRS. Thus, record the following
    # files separately
    target_sources(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PUBLIC
      FILE_SET external_module_interface_units
      TYPE CXX_MODULES
      FILES ${_external_interface_module_partition_units}
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/interface_module_partitions"
      )

    target_sources(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PRIVATE
      FILE_SET module_implementation_units
      TYPE CXX_MODULES
      FILES
        ${_implementation_module_partition_units}
      BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
      )

    target_sources(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PUBLIC
      FILE_SET headers
      TYPE HEADERS
      FILES ${_dealii_macros_header}
      BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/../include/
      )

    #
    # In addition, link in all objects from the bundled object targets
    # (listed in the DEAL_II_BUNDLED_TARGETS_* global property).
    #
    get_property(_bundled_object_targets GLOBAL PROPERTY DEAL_II_BUNDLED_TARGETS_${build})
    set(_object_files "")
    foreach(_target ${_bundled_object_targets})
      list(APPEND _object_files "$<TARGET_OBJECTS:${_target}>")
    endforeach()
    target_sources(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PRIVATE ${_object_files}
      )

    # Ensure that all .inst files and the macro headers are in place
    add_dependencies(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      expand_all_instantiations
      build_macros_header
      )

    #
    # Record the expected C++ standard as a compile feature. This target
    # property ensures that support for our expected C++ standard is always
    # enabled in client user code irrespective of what compile flags/options
    # they have set.
    #
    target_compile_features(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      INTERFACE cxx_std_${CMAKE_CXX_STANDARD}
      )

    set_target_properties(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
      PROPERTIES
      LINKER_LANGUAGE "CXX"
      VERSION "${DEAL_II_PACKAGE_VERSION}"
      #
      # Sonaming: Well... we just use the version number.
      # No point to wrack one's brain over the question whether a new version of
      # a C++ library is still ABI backwards compatible :-]
      #
      SOVERSION "${DEAL_II_PACKAGE_VERSION}"
      ARCHIVE_OUTPUT_NAME "${DEAL_II_BASE_NAME}_module${DEAL_II_${build}_SUFFIX}"
      LIBRARY_OUTPUT_NAME "${DEAL_II_BASE_NAME}_module${DEAL_II_${build}_SUFFIX}"
      RUNTIME_OUTPUT_NAME "${DEAL_II_BASE_NAME}_module${DEAL_II_${build}_SUFFIX}"
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
      )

    # Use the same logic to determine which compile and link flags to use
    # as for the non-module libraries:
    target_compile_flags(${DEAL_II_TARGET_NAME}_module_${build_lowercase} INTERFACE
      "$<AND:$<CONFIG:${build_camelcase}>,$<COMPILE_LANGUAGE:CXX>>"
      "${DEAL_II_CXX_FLAGS_${build}}"
      )
    target_link_flags(${DEAL_II_TARGET_NAME}_module_${build_lowercase} INTERFACE
      "$<CONFIG:${build_camelcase}>" "${DEAL_II_LINKER_FLAGS_${build}}"
      )

    if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
      set_target_properties(${DEAL_II_TARGET_NAME}_module_${build_lowercase}
        PROPERTIES
        MACOSX_RPATH OFF
        BUILD_WITH_INSTALL_RPATH OFF
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}"
        )
    endif()

    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.30")
      #
      # Create alias targets "dealii::dealii_module_release" and
      # "dealii::dealii_module_debug" to have the exported target
      # names available when populating the link interface of
      # downstream targets:
      #
      add_library(${DEAL_II_TARGET_NAME}::${DEAL_II_TARGET_NAME}_module_${build_lowercase}
        ALIAS ${DEAL_II_TARGET_NAME}_module_${build_lowercase}
        )
    endif()

    # Finally, let the top-level 'library' target depend on this module:
    add_dependencies(library ${DEAL_II_TARGET_NAME}_module_${build_lowercase})
  endforeach()

  message(STATUS "Setting up the C++20 module - Done")
endif()
