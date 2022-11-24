## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2022 by the deal.II authors
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

#
# deal_ii_package_handle(<feature>
#  {<conf. variable> {(REQUIRED|OPTIONAL) <variables>}}
#  [CLEAR <variables>]
#  )
#
# This macro is an alternative implementation of the
# FIND_PACKAGE_HANDLE_STANDARD_ARGS macro shipped with CMake - aka do
# everything that was expected from CMake in the first place *sigh*
#
# Its usage is best explained with an example:
#
#   deal_ii_package_handle(PETSC
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
# different from "-NOTFOUND". If so, PETSC_LIBRARIES and PETSC_INCLUDE_DIRS
# is defined and populated with the contents of all specified variables.
# Optional variables with no content or whose content is "-NOTFOUND" are
# filtered out.
# After the 'CLEAR' statement all internally cached variables should be
# listed - this is used to provide a possibility to undo a feature
# search.
#

macro(DEAL_II_PACKAGE_HANDLE _feature)

  if(DEFINED ${_feature}_VERSION)
    message(STATUS "  ${_feature}_VERSION: ${${_feature}_VERSION}")
  endif()

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

    elseif("${_arg}" STREQUAL "REQUIRED")
      set(_required TRUE)

    elseif("${_arg}" STREQUAL "OPTIONAL")
      set(_required FALSE)

    elseif(_arg MATCHES "^(optimized|debug|general)$"
            AND "${_current_suffix}" STREQUAL "LIBRARIES")
      list(APPEND _temp_${_current_suffix} ${_arg})

    else()
      if ("${_current_suffix}" STREQUAL "")
        message(FATAL_ERROR
          "Internal configuration error: the second "
          "argument to DEAL_II_PACKAGE_HANDLE must be a keyword"
          )
      endif()

      mark_as_advanced(${_arg})

      if("${_current_suffix}" STREQUAL "CLEAR")
        if(NOT _arg MATCHES "^(optimized|debug|general)$")
          list(APPEND _clear_variables_list ${_arg})
        endif()

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

    #
    # Write back into global variables:
    #
    clear_feature(${_feature})
    foreach(_suffix ${DEAL_II_LIST_SUFFIXES} ${DEAL_II_STRING_SUFFIXES})
      if(NOT "${_temp_${_suffix}}" STREQUAL "")
        set(${_feature}_${_suffix} "${_temp_${_suffix}}")
        message(STATUS "  ${_feature}_${_suffix}: ${${_feature}_${_suffix}}")
      endif()
    endforeach()

    #
    # Remove certain system libraries from the link interface. This is
    # purely cosmetic (we always implicitly link against the C library, and
    # we always set up threading by linking against libpthread.so if
    # necessary).
    #
    foreach(_suffix LIBRARIES LIBRARIES_DEBUG LIBRARIES_RELEASE)
      if(NOT "${${_feature}_${_suffix}}" STREQUAL "")
        list(REMOVE_ITEM ${_feature}_${_suffix}
          "pthread" "-pthread" "-lpthread" "c" "-lc"
          )
      endif()
    endforeach()

    message(STATUS "Found ${_feature}")

    mark_as_advanced(${_feature}_DIR ${_feature}_ARCH)

  else()

    message(STATUS "Could NOT find ${_feature}")
  endif()
endmacro()
