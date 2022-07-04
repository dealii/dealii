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
# DEAL_II_PACKAGE_HANDLE(<feature>
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
#   DEAL_II_PACKAGE_HANDLE(PETSC
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

MACRO(DEAL_II_PACKAGE_HANDLE _feature)

  IF(DEFINED ${_feature}_VERSION)
    MESSAGE(STATUS "  ${_feature}_VERSION: ${${_feature}_VERSION}")
  ENDIF()

  #
  # Respect a possible ${_feature}_FOUND variable that is set to a truth
  # value. We need this for modern™ MPI detection where CMake's
  # FindMPI.cmake might only set MPI_FOUND to true and nothing else.
  #
  IF(NOT DEFINED ${_feature}_FOUND)
    SET(${_feature}_FOUND TRUE)
  ENDIF()

  #
  # Clear temporary variables
  #
  FOREACH(_suffix ${DEAL_II_LIST_SUFFIXES} ${DEAL_II_STRING_SUFFIXES})
    set(_temp_${_suffix} "")
  ENDFOREACH()

  #
  # State variables for parsing keywords and arguments. We store the
  # currently encountered keyword in ${_current_suffix} and store the fact
  # whether we encountered an "OPTIONAL" or "REQUIRED" keyword in
  # ${_required}
  #
  SET(_current_suffix "")
  SET(_required TRUE)

  #
  # A temporary list accumulating all variables that should be "cleared"
  # when the feature gets disabled.
  #
  SET(_clear_variables_list "")

  FOREACH(_arg ${ARGN})
    IF(("${_arg}" IN_LIST DEAL_II_LIST_SUFFIXES) OR
       ("${_arg}" IN_LIST DEAL_II_STRING_SUFFIXES) OR
       ("${_arg}" STREQUAL "CLEAR"))
      #
      # We encountered a new keyword.
      #
      SET(_current_suffix "${_arg}")

    ELSEIF("${_arg}" STREQUAL "REQUIRED")
      SET(_required TRUE)

    ELSEIF("${_arg}" STREQUAL "OPTIONAL")
      SET(_required FALSE)

    ELSEIF(_arg MATCHES "^(optimized|debug|general)$"
            AND "${_current_suffix}" STREQUAL "LIBRARIES")
      LIST(APPEND _temp_${_current_suffix} ${_arg})

    ELSE()
      IF ("${_current_suffix}" STREQUAL "")
        MESSAGE(FATAL_ERROR
          "Internal configuration error: the second "
          "argument to DEAL_II_PACKAGE_HANDLE must be a keyword"
          )
      ENDIF()

      MARK_AS_ADVANCED(${_arg})

      IF("${_current_suffix}" STREQUAL "CLEAR")
        IF(NOT _arg MATCHES "^(optimized|debug|general)$")
          LIST(APPEND _clear_variables_list ${_arg})
        ENDIF()

      ELSE()

        IF("${${_arg}}" MATCHES "^\\s*$" OR "${${_arg}}" MATCHES "-NOTFOUND")
          IF(_required)
            IF("${${_arg}}" MATCHES "^\\s*$")
              MESSAGE(STATUS
                "  ${_feature}_${_current_suffix}: *** Required variable \"${_arg}\" empty ***"
                )
            ELSE()
              MESSAGE(STATUS
                "  ${_feature}_${_current_suffix}: *** Required variable \"${_arg}\" set to NOTFOUND ***"
                )
            ENDIF()
            SET(${_feature}_FOUND FALSE)
          ENDIF()
        ELSE()
          LIST(APPEND _temp_${_current_suffix} ${${_arg}})
        ENDIF()
      ENDIF()
    ENDIF()
  ENDFOREACH()

  SET(${_feature}_CLEAR_VARIABLES ${_clear} CACHE INTERNAL "")

  IF(${_feature}_FOUND)
    #
    # Deduplicate and stringify entries:
    #
    FOREACH(_suffix ${DEAL_II_LIST_SUFFIXES})
      IF(_suffix MATCHES "INCLUDE_DIRS$")
        REMOVE_DUPLICATES(_temp_${_suffix})
      ELSE()
        REMOVE_DUPLICATES(_temp_${_suffix} REVERSE)
      ENDIF()
    ENDFOREACH()
    FOREACH(_suffix ${_DEAL_II_STRING_SUFFIXES})
      TO_STRING(_temp_${_suffix} ${_temp_${_suffix}})
    ENDFOREACH()

    #
    # Write back into global variables:
    #
    CLEAR_FEATURE(${_feature})
    FOREACH(_suffix ${DEAL_II_LIST_SUFFIXES} ${DEAL_II_STRING_SUFFIXES})
      IF(NOT "${_temp_${_suffix}}" STREQUAL "")
        set(${_feature}_${_suffix} "${_temp_${_suffix}}")
        MESSAGE(STATUS "  ${_feature}_${_suffix}: ${${_feature}_${_suffix}}")
      ENDIF()
    ENDFOREACH()

    #
    # Remove certain system libraries from the link interface. This is
    # purely cosmetic (we always implicitly link against the C library, and
    # we always set up threading by linking against libpthread.so if
    # necessary).
    #
    FOREACH(_suffix LIBRARIES LIBRARIES_DEBUG LIBRARIES_RELEASE)
      IF(NOT "${${_feature}_${_suffix}}" STREQUAL "")
        LIST(REMOVE_ITEM ${_feature}_${_suffix}
          "pthread" "-pthread" "-lpthread" "c" "-lc"
          )
      ENDIF()
    ENDFOREACH()

    MESSAGE(STATUS "Found ${_feature}")

    MARK_AS_ADVANCED(${_feature}_DIR ${_feature}_ARCH)

  ELSE()

    MESSAGE(STATUS "Could NOT find ${_feature}")
  ENDIF()
ENDMACRO()
