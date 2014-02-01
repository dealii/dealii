## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 - 2014 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II. It adds
# individual FEATURE_* configuration variables to the corresponding
# DEAL_II_* variables
#
# Usage:
#     REGISTER_FEATURE(feature)
#
# This macro will add
#
#   <FEATURE>_LIBRARIES (respecting general, optimized, debug keyword)
#   <FEATURE>_LIBRARIES(_DEBUG|_RELEASE)
#   <FEATURE>_(|BUNDLED_|USER_)INCLUDE_DIRS
#   <FEATURE>_DEFINITIONS(|_DEBUG|_RELEASE)
#   <FEATURE>_USER_DEFINITIONS(|_DEBUG|_RELEASE)
#   <FEATURE>_CXX_FLAGS(|_DEBUG|_RELEASE)
#   <FEATURE>_LINKER_FLAGS(|_DEBUG|_RELEASE)
#
# to the corresponding DEAL_II_* variables
#

MACRO(REGISTER_FEATURE _feature)

  IF(DEFINED ${_feature}_LIBRARIES)
    #
    # Add ${_feature}_LIBRARIES to
    #   DEAL_II_LIBRARIES
    #   DEAL_II_LIBRARIES_DEBUG
    #   DEAL_II_LIBRARIES_RELEASE
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
          LIST(APPEND DEAL_II_LIBRARIES ${_tmp})
        ELSEIF("${_toggle}" STREQUAL "debug")
          LIST(APPEND DEAL_II_LIBRARIES_DEBUG ${_tmp})
        ELSEIF("${_toggle}" STREQUAL "optimized")
          LIST(APPEND DEAL_II_LIBRARIES_RELEASE ${_tmp})
        ENDIF()
      ENDIF()
    ENDFOREACH()
  ENDIF()

  FOREACH(_var
    LIBRARIES_DEBUG LIBRARIES_RELEASE
    INCLUDE_DIRS BUNDLED_INCLUDE_DIRS USER_INCLUDE_DIRS
    DEFINITIONS DEFINITIONS_DEBUG DEFINITIONS_RELEASE
    USER_DEFINITIONS USER_DEFINITIONS_DEBUG USER_DEFINITIONS_RELEASE
    )
    IF(DEFINED ${_feature}_${_var})
      LIST(APPEND DEAL_II_${_var} ${${_feature}_${_var}})
    ENDIF()
  ENDFOREACH()

  FOREACH(_var
    CXX_FLAGS CXX_FLAGS_DEBUG CXX_FLAGS_RELEASE
    LINKER_FLAGS LINKER_FLAGS_DEBUG LINKER_FLAGS_RELEASE
    )
    IF(DEFINED ${_feature}_${_var})
      ADD_FLAGS(DEAL_II_${_var} "${${_feature}_${_var}}")
    ENDIF()
  ENDFOREACH()

ENDMACRO()
