## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II
#
# Usage:
#     REGISTER_FEATURE(feature)
#
# This macro will
#
#   - add ${feature}_INCLUDE_DIRS and ${feature}_INCLUDE_PATH
#     to the list of (internal) include directories
#   - and if ${feature}_ADD_TO_USER_INCLUDE_DIRS is defined also to
#     DEAL_II_USER_INCLUDE_DIRS
#
#   - add ${feature}_LINKER_FLAGS and ${feature}_LINK_FLAGS to
#     DEAL_II_LINKER_FLAGS
#
#   - add ${feature}_CXX_FLAGS and ${feature}_COMPILE_FLAGS to
#     CMAKE_CXX_FLAGS
#
#   - add ${feature}_LIBRARIES to the list of deal.II libraries depending
#     on general, optimized or debug keyword
#


MACRO(REGISTER_FEATURE _feature)
  # variables for include directories:
  FOREACH(_var ${_feature}_INCLUDE_DIRS ${_feature}_INCLUDE_PATH)
    IF(DEFINED ${_var})
      INCLUDE_DIRECTORIES(${${_var}})
      IF(${_feature}_ADD_TO_USER_INCLUDE_DIRS)
        LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${${_var}})
      ENDIF()
    ENDIF()
  ENDFOREACH()

  # variables for linker flags:
  FOREACH(_var ${_feature}_LINKER_FLAGS ${_feature}_LINK_FLAGS)
    IF(DEFINED ${_var})
      ADD_FLAGS(DEAL_II_LINKER_FLAGS "${${_var}}")
    ENDIF()
  ENDFOREACH()

  # variables for compiler flags:
  FOREACH(_var ${_feature}_CXX_FLAGS ${_feature}_COMPILE_FLAGS)
    IF(DEFINED ${_var})
      ADD_FLAGS(CMAKE_CXX_FLAGS "${${_var}}")
    ENDIF()
  ENDFOREACH()

  IF(DEFINED ${_feature}_LIBRARIES)
    #
    # Add ${_feature}_LIBRARIES to
    #   DEAL_II_EXTERNAL_LIBRARIES
    #   DEAL_II_EXTERNAL_LIBRARIES_DEBUG
    #   DEAL_II_EXTERNAL_LIBRARIES_RELEASE
    # depending on the "optmized", "debug" or "general" keyword
    #
    SET(_toggle "general")
    FOREACH(_tmp ${${_feature}_LIBRARIES})
      IF( "${_tmp}" STREQUAL "debug" OR
          "${_tmp}" STREQUAL "optimized" OR
          "${_tmp}" STREQUAL "general" )
        SET(_toggle "${_tmp}")
      ELSE()
        IF("${_toggle}" STREQUAL "general")
          LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${_tmp})
        ELSEIF("${_toggle}" STREQUAL "debug")
          LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_DEBUG ${_tmp})
        ELSEIF("${_toggle}" STREQUAL "optimized")
          LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_RELEASE ${_tmp})
        ENDIF()
      ENDIF()
    ENDFOREACH()
  ENDIF()

ENDMACRO()
