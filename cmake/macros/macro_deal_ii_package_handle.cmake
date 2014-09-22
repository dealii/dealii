## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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

MACRO(DEAL_II_PACKAGE_HANDLE _feature _var)

  IF(DEFINED ${_feature}_VERSION)
    MESSAGE(STATUS "  ${_feature}_VERSION: ${${_feature}_VERSION}")
  ENDIF()

  SET(${_feature}_FOUND TRUE)

  SET(_variable ${_var})
  SET(${_feature}_${_variable} "")
  SET(_required TRUE)
  SET(_fine TRUE)
  SET(_fill_clear FALSE)
  SET(_clear "")

  FOREACH(_arg ${ARGN})
    IF(_arg MATCHES "^LIBRARIES(|_DEBUG|_RELEASE)$"
       OR _arg MATCHES "^(|BUNDLED_|USER_)INCLUDE_DIRS$"
       OR _arg MATCHES "^(|USER_)DEFINITIONS(|_DEBUG|_RELEASE)$"
       OR _arg MATCHES "^CXX_FLAGS(|_DEBUG|_RELEASE)"
       OR _arg MATCHES "^LINKER_FLAGS(|_DEBUG|_RELEASE)")

      IF(_fine)
        IF(_variable MATCHES "^CXX_FLAGS(|_DEBUG|_RELEASE)"
           OR _variable MATCHES "^LINKER_FLAGS(|_DEBUG|_RELEASE)")
          TO_STRING(${_feature}_${_variable} ${${_feature}_${_variable}})
        ENDIF()
        MESSAGE(STATUS "  ${_feature}_${_variable}: ${${_feature}_${_variable}}")
      ENDIF()

      #
      # *Yay* a new keyword.
      #
      SET(_variable ${_arg})
      SET(${_feature}_${_variable} "")
      SET(_required TRUE)
      SET(_fine TRUE)

    ELSEIF("${_arg}" STREQUAL "REQUIRED")
      SET(_required TRUE)
    ELSEIF("${_arg}" STREQUAL "OPTIONAL")
      SET(_required FALSE)
    ELSEIF(_arg MATCHES "^(optimized|debug|general)$"
            AND "${_variable}" STREQUAL "LIBRARIES")
      #
      # Keywords are special...
      #
      LIST(APPEND ${_feature}_${_variable} ${_arg})
    ELSEIF("${_arg}" STREQUAL "CLEAR")
      SET(_fill_clear TRUE)
    ELSE()
      MARK_AS_ADVANCED(${_arg})
      IF(_fill_clear)
        IF(NOT _arg MATCHES "^(optimized|debug|general)$")
          LIST(APPEND _clear ${_arg})
        ENDIF()
      ELSE()
        IF("${${_arg}}" MATCHES "^\\s*$" OR "${${_arg}}" MATCHES "-NOTFOUND")
          IF(_required AND _fine)
            IF("${${_arg}}" MATCHES "^\\s*$")
              MESSAGE(STATUS
                "  ${_feature}_${_variable}: *** Required variable \"${_arg}\" empty ***"
                )
            ELSE()
              MESSAGE(STATUS
                "  ${_feature}_${_variable}: *** Required variable \"${_arg}\" set to NOTFOUND ***"
                )
            ENDIF()
            SET(${_feature}_FOUND FALSE)
            SET(_fine FALSE)
          ENDIF()
        ELSE()
          LIST(APPEND ${_feature}_${_variable} ${${_arg}})
        ENDIF()
      ENDIF()
    ENDIF()
  ENDFOREACH()

  SET(${_feature}_CLEAR_VARIABLES ${_clear} CACHE INTERNAL "")

  IF(_fine)
    IF(_variable MATCHES "^CXX_FLAGS(|_DEBUG|_RELEASE)"
       OR _variable MATCHES "^LINKER_FLAGS(|_DEBUG|_RELEASE)")
      TO_STRING(${_feature}_${_variable} ${${_feature}_${_variable}})
    ENDIF()
    MESSAGE(STATUS "  ${_feature}_${_variable}: ${${_feature}_${_variable}}")
  ENDIF()

  IF(${_feature}_FOUND)
    #
    # Deduplicate entries:
    #
    FOREACH(_suffix ${DEAL_II_LIST_SUFFIXES})
      IF(_suffix MATCHES "INCLUDE_DIRS$")
        REMOVE_DUPLICATES(${_feature}_${_suffix})
      ELSE()
        REMOVE_DUPLICATES(${_feature}_${_suffix} REVERSE)
      ENDIF()
    ENDFOREACH()

    MESSAGE(STATUS "Found ${_feature}")

    MARK_AS_ADVANCED(${_feature}_DIR ${_feature}_ARCH)

  ELSE()

    MESSAGE(STATUS "Could NOT find ${_feature}")
  ENDIF()
ENDMACRO()
