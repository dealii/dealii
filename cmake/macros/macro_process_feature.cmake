## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2025 by the deal.II authors
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
# process_feature(<feature>
#  {<conf. variable> {(REQUIRED|OPTIONAL) <variables>}}
#  [CLEAR <variables>]
#  )
#
# This macro processes a configuration statement similar to
# FIND_PACKAGE_HANDLE_STANDARD_ARGS() macro shipped with CMake - but in
# contrast to the CMake macro it populates corresponding FEATURE_*
# variables.
#
# Its usage is best explained with an example:
#
#   process_feature(PETSC
#     LIBRARIES
#       REQUIRED PETSC_LIBRARY
#       OPTIONAL _petsc_libraries
#     INCLUDE_DIRS
#       REQUIRED PETSC_INCLUDE_DIR_COMMON PETSC_INCLUDE_DIR_ARCH
#       OPTIONAL _petsc_includes
#     CLEAR PETSC_LIBRARY PETSC_INCLUDE_DIR_COMMON PETSC_INCLUDE_DIR_ARCH
#     )
#
# This will check whether all REQUIRED variables are non-empty and
# different from "-NOTFOUND". If so, the PETSC_LIBRARIES and
# PETSC_INCLUDE_DIRS will be defined and populated with the contents of all
# specified variables. Optional variables with no content or whose content
# is "-NOTFOUND" are filtered out.
#
# After the 'CLEAR' statement all internally cached variables should be
# listed - this is used to provide a possibility to undo a feature
# search.
#
# Valid suffixes are
#   TARGETS      TARGETS_RELEASE      TARGETS_DEBUG
#   LIBRARIES    LIBRARIES_RELEASE    LIBRARIES_DEBUG
#   INCLUDE_DIRS
#   DEFINITIONS  DEFINITIONS_RELEASE  DEFINITIONS_DEBUG
#   CXX_FLAGS    CXX_FLAGS_RELEASE    CXX_FLAGS_DEBUG
#   LINKER_FLAGS LINKER_FLAGS_RELEASE LINKER_FLAGS_DEBUG
#   EXECUTABLE
#
# Note: The contents of all <feature>_<suffix> variables will be cleared
#
# In addition the macro defines the booleans
#
#   <feature>_FOUND                -  true if processing was successful
#   <feature>_SPLIT_CONFIGURATION  -  true if separate debug and release
#                                     values have been defined.
#

macro(process_feature _feature)

  message(STATUS "Processing ${_feature} variables and targets")

  #
  # Respect a possible ${_feature}_FOUND variable that is set to a truth
  # value. We need this for modern™ MPI detection where CMake's
  # FindMPI.cmake might only set MPI_FOUND to true and nothing else.
  #
  if(NOT DEFINED ${_feature}_FOUND)
    set(${_feature}_FOUND TRUE)
  endif()

  #
  # Clear temporary variables
  #
  foreach(_suffix ${DEAL_II_LIST_SUFFIXES} ${DEAL_II_STRING_SUFFIXES})
    set(_temp_${_suffix} "")
  endforeach()

  #
  # State variables for parsing keywords and arguments. We store the
  # currently encountered keyword in ${_current_suffix} and store the fact
  # whether we encountered an "OPTIONAL" or "REQUIRED" keyword in
  # ${_required}
  #
  set(_current_suffix "")
  set(_required TRUE)

  #
  # Record whether we have encountered split *_DEBUG/*_RELEASE variables
  #
  set(_split_configuration FALSE)

  #
  # A temporary list accumulating all variables that should be "cleared"
  # when the feature gets disabled.
  #
  set(_clear_variables_list "")

  foreach(_arg ${ARGN})
    if(("${_arg}" IN_LIST DEAL_II_LIST_SUFFIXES) OR
       ("${_arg}" IN_LIST DEAL_II_STRING_SUFFIXES) OR
       ("${_arg}" STREQUAL "CLEAR"))
      #
      # We encountered a new keyword.
      #
      set(_current_suffix "${_arg}")

      if(_current_suffix MATCHES "_DEBUG" OR _current_suffix MATCHES "_RELEASE")
        set(_split_configuration TRUE)
      endif()

    elseif("${_arg}" STREQUAL "REQUIRED")
      set(_required TRUE)

    elseif("${_arg}" STREQUAL "OPTIONAL")
      set(_required FALSE)

    else()
      if ("${_current_suffix}" STREQUAL "")
        message(FATAL_ERROR
          "Internal configuration error: the second "
          "argument to process_feature must be a keyword"
          )
      elseif(_arg MATCHES "^(optimized|debug|general)$")
        message(FATAL_ERROR
          "Internal configuration error: process_feature() does not support "
          "\"debug«, »optimized«, or »general\" library identifiers, use the "
          "appropriate keyword instead."
          )
      endif()

      mark_as_advanced(${_arg})

      if("${_current_suffix}" STREQUAL "CLEAR")
        list(APPEND _clear_variables_list ${_arg})

      else()

        if("${${_arg}}" MATCHES "^\\s*$" OR "${${_arg}}" MATCHES "-NOTFOUND")
          if(_required)
            if("${${_arg}}" MATCHES "^\\s*$")
              message(STATUS
                "  ${_feature}_${_current_suffix}: *** Required variable \"${_arg}\" empty ***"
                )
            else()
              message(STATUS
                "  ${_feature}_${_current_suffix}: *** Required variable \"${_arg}\" set to NOTFOUND ***"
                )
            endif()
            set(${_feature}_FOUND FALSE)
          endif()
        else()
          list(APPEND _temp_${_current_suffix} ${${_arg}})
        endif()
      endif()
    endif()
  endforeach()

  set(${_feature}_CLEAR_VARIABLES ${_clear} CACHE INTERNAL "")

  if(${_feature}_FOUND)
    #
    # Take care of "optimized", "debug", and "general" keywords in
    # LIBRARIES and TARGETS variables by distributing affected entries into
    # the respective LIBRARIES_(DEBUG|RELEASE) and TARGETS_(DEBUG|RELEASE)
    # variables.
    #
    foreach(_suffix LIBRARIES TARGETS)
      set(_temp ${_temp_${_suffix}})
      set(_temp_${_suffix} "")
      set(_switch "")
      foreach(_entry ${_temp})
        if("${_entry}" STREQUAL "optimized")
          set(_split_configuration TRUE)
          set(_switch "_RELEASE")
        elseif("${_entry}" STREQUAL "debug")
          set(_split_configuration TRUE)
          set(_switch "_DEBUG")
        elseif("${_entry}" STREQUAL "general")
          set(_switch "")
        else()
          list(APPEND _temp_${_suffix}${_switch} ${_entry})
        endif()
      endforeach()
    endforeach()

    #
    # Deduplicate and stringify entries:
    #
    foreach(_suffix ${DEAL_II_LIST_SUFFIXES})
      if(_suffix MATCHES "INCLUDE_DIRS$")
        remove_duplicates(_temp_${_suffix})
      else()
        remove_duplicates(_temp_${_suffix} REVERSE)
      endif()
    endforeach()
    foreach(_suffix ${_DEAL_II_STRING_SUFFIXES})
      to_string(_temp_${_suffix} ${_temp_${_suffix}})
    endforeach()
    set(${_feature}_SPLIT_CONFIGURATION ${_split_configuration})

    #
    # Remove certain system libraries from the link interface. This is
    # purely cosmetic (we always implicitly link against the C library, and
    # we always set up threading by linking against libpthread.so if
    # necessary).
    #
    foreach(_suffix LIBRARIES LIBRARIES_DEBUG LIBRARIES_RELEASE)
      if(NOT "${_temp_${_suffix}}" STREQUAL "")
        list(REMOVE_ITEM _temp_${_suffix}
          "pthread" "-pthread" "-lpthread" "c" "-lc"
          )
      endif()
    endforeach()

    #
    # Write back into global variables:
    #
    clear_feature(${_feature})
    foreach(_suffix ${DEAL_II_LIST_SUFFIXES} ${DEAL_II_STRING_SUFFIXES})
      if(NOT "${_temp_${_suffix}}" STREQUAL "")
        set(${_feature}_${_suffix} "${_temp_${_suffix}}")
      endif()
    endforeach()

    mark_as_advanced(${_feature}_DIR ${_feature}_ARCH)

    message(STATUS "Processing ${_feature} variables and targets - Done")

  else()

    clear_feature(${_feature})
    message(STATUS "Unable to process ${_feature}")
  endif()
endmacro()
