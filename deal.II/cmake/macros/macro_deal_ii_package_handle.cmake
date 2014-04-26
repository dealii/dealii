## ---------------------------------------------------------------------
## $Id$
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
# DEAL_II_PACKAGE_HANDLE
#
# TODO: Documentation
#

MACRO(DEAL_II_PACKAGE_HANDLE _feature _var)

  IF(DEFINED ${_feature}_VERSION)
    MESSAGE(STATUS "  ${_feature}_VERSION: ${${_feature}_VERSION}")
  ENDIF()

  SET(${_feature}_FOUND TRUE)

  SET(_variable ${_var})
  SET(_cleanup ${_var})
  SET(${_feature}_${_variable} "")
  SET(_required TRUE)
  SET(_fine TRUE)

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
      LIST(APPEND _cleanup ${_var})
      SET(${_feature}_${_variable} "")
      SET(_required TRUE)
      SET(_fine TRUE)

    ELSEIF("${_arg}" STREQUAL "REQUIRED")
      SET(_required TRUE)
    ELSEIF("${_arg}" STREQUAL "OPTIONAL")
      SET(_required FALSE)
    ELSEIF( _arg MATCHES "^(optimized|debug|general)$"
            AND "${_variable}" STREQUAL "LIBRARIES")
      #
      # Keywords are special...
      #
      LIST(APPEND ${_feature}_${_variable} ${_arg})
    ELSE()
      MARK_AS_ADVANCED(${_arg})
      IF(NOT DEFINED ${_arg} OR ${_arg} MATCHES "-NOTFOUND")
        IF(_required AND _fine)
          IF(NOT DEFINED ${_arg})
            MESSAGE(STATUS
              "  ${_feature}_${_variable}: *** Required variable \"${_arg}\" undefined ***"
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
  ENDFOREACH()

  IF(_fine)
    IF(_variable MATCHES "^CXX_FLAGS(|_DEBUG|_RELEASE)"
       OR _variable MATCHES "^LINKER_FLAGS(|_DEBUG|_RELEASE)")
      TO_STRING(${_feature}_${_variable} ${${_feature}_${_variable}})
    ENDIF()
    MESSAGE(STATUS "  ${_feature}_${_variable}: ${${_feature}_${_variable}}")
  ENDIF()

  IF(${_feature}_FOUND)
    #
    # Deduplicate entries in *_INCLUDE_DIRS and *_LIBRARIES
    #
    FOREACH(_suffix INCLUDE_DIRS USER_INCLUDE_DIRS BUNDLED_INCLUDE_DIRS)
      REMOVE_DUPLICATES(${_feature}_${_suffix})
    ENDFOREACH()
    FOREACH(_suffix
        LIBRARIES LIBRARIES_RELEASE LIBRARIES_DEBUG
        USER_DEFINITIONS USER_DEFINITIONS_DEBUG USER_DEFINITIONS_RELEASE
        DEFINITIONS DEFINITIONS_DEBUG DEFINITIONS_RELEASE
        )
      REMOVE_DUPLICATES(${_feature}_${_suffix} REVERSE)
    ENDFOREACH()

    MESSAGE(STATUS "Found ${_feature}")
    MARK_AS_ADVANCED(${_feature}_DIR)

  ELSE()

    FOREACH(_v _cleanup)
      SET(${_feature}_${_v})
    ENDFOREACH()
    MESSAGE(STATUS "Could NOT find ${_feature}")
  ENDIF()
ENDMACRO()
