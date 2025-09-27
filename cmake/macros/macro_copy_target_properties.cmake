## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2023 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# copy_target_properties(<destination target> [<source targets>])
#
# Copies relevant (imported) target properties from the source targets and
# appends the properties to the destination target.
#

function(copy_target_properties _destination_target)

  #
  # Make sure that we expand all targets recursively for the sake of a
  # cleaner output we check whether we have already included a given target
  #

  set(_source_targets ${ARGN})
  set(_processed_targets)

  #
  # Collect values over all targets found in ${_source_targets}:
  #

  set(_libraries)
  set(_include_directories)
  set(_compile_definitions)
  set(_compile_options)
  set(_link_options)

  while(NOT "${_source_targets}" STREQUAL "")
    # Future FIXME: replace by POP_FRONT (CMake 3.15)
    list(GET _source_targets 0 _entry)
    list(REMOVE_AT _source_targets 0)

    if(${_entry} IN_LIST _processed_targets)
      continue()
    endif()

    list(APPEND _processed_targets ${_entry})
    message(STATUS "    copying ${_entry} into ${_destination_target} ...")

    #
    # Only query the LOCATION from non-interface libraries. An interface
    # library is a CMake target that only consists of "INTERFACE" target
    # properties and does not by itself refer to an executable or shared
    # object/static archive. CMake prior to 3.19 will throw an error if we
    # query for the LOCATION. Newer CMake versions simply return
    # "-NOTFOUND".
    #
    set(_location)
    get_target_property(_test ${_entry} TYPE)
    if(NOT "${_test}" STREQUAL "INTERFACE_LIBRARY")
      get_target_property(_location ${_entry} LOCATION)
      if("${_location}" MATCHES "-NOTFOUND")
        set(_location)
      endif()
    endif()

    get_target_property(_values ${_entry} INTERFACE_LINK_LIBRARIES)
    if("${_values}" MATCHES "-NOTFOUND")
      set(_values)
    endif()

    foreach(_lib ${_location} ${_values})
      strip_known_generator_expressions(_lib)
      if(TARGET ${_lib})
        list(APPEND _source_targets ${_lib})
      else()
        #
        # Complain loudly if we encounter an undefined target:
        #
        if("${_lib}" MATCHES "::")
          message(FATAL_ERROR
            "Undefined imported target name \"${_lib}\" present in interface "
            "of target \"${_entry}\"."
            )
        endif()
        list(APPEND _libraries ${_lib})
      endif()
    endforeach()

    get_target_property(_values ${_entry} INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT "${_values}" MATCHES "-NOTFOUND")
      list(APPEND _include_directories ${_values})
    endif()

    get_target_property(_values ${_entry} INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    if(NOT "${_values}" MATCHES "-NOTFOUND")
      list(APPEND _include_directories ${_values})
    endif()

    get_target_property(_values ${_entry} INTERFACE_COMPILE_DEFINITIONS)
    if(NOT "${_values}" MATCHES "-NOTFOUND")
      strip_known_generator_expressions(_values)
      list(APPEND _compile_definitions ${_values})
    endif()

    get_target_property(_values ${_entry} INTERFACE_COMPILE_OPTIONS)
    if(NOT "${_values}" MATCHES "-NOTFOUND")
      strip_known_generator_expressions(_values)
      list(APPEND _compile_options ${_values})
    endif()

    get_target_property(_values ${_entry} INTERFACE_LINK_OPTIONS)
    if(NOT "${_values}" MATCHES "-NOTFOUND")
      strip_known_generator_expressions(_values)
      list(APPEND _link_options ${_values})
    endif()
  endwhile()

  #
  # Populate destination target:
  #

  if(NOT "${_libraries}" STREQUAL "")
    remove_duplicates(_libraries REVERSE)
    message(STATUS "    LINK_LIBRARIES:      ${_libraries}")
    target_link_libraries(${_destination_target} INTERFACE ${_libraries})
  endif()

  if(NOT "${_include_directories}" STREQUAL "")
    remove_duplicates(_include_directories)
    message(STATUS "    INCLUDE_DIRECTORIES: ${_include_directories}")
    target_include_directories(${_destination_target} SYSTEM INTERFACE ${_include_directories})
  endif()

  if(NOT "${_compile_definitions}" STREQUAL "")
    remove_duplicates(_compile_definitions REVERSE)
    message(STATUS "    COMPILE_DEFINITIONS: ${_compile_definitions}")
    target_compile_definitions(${_destination_target} INTERFACE ${_compile_definitions})
  endif()

  if(NOT "${_compile_options}" STREQUAL "")
    remove_duplicates(_compile_options REVERSE)
    message(STATUS "    COMPILE_OPTIONS: ${_compile_options}")
    target_compile_options(${_destination_target} INTERFACE ${_compile_options})
  endif()

  if(NOT "${_link_options}" STREQUAL "")
    remove_duplicates(_link_options REVERSE)
    message(STATUS "    LINK_OPTIONS: ${_link_options}")
    target_link_options(${_destination_target} INTERFACE ${_link_options})
  endif()
endfunction()
