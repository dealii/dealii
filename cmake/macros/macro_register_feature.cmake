## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2015 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II. It adds
# individual FEATURE_* configuration variables to the corresponding
# DEAL_II_* variables
#
# Usage:
#     register_feature(feature)
#
# This macro will add
#
#   <FEATURE>_LIBRARIES (respecting general, optimized, debug keyword)
#
# and all other suffixes defined in DEAL_II_LIST_SUFFIXES and
# DEAL_II_STRING_SUFFIXES to the corresponding DEAL_II_* variables
#

macro(register_feature _feature)

  if(DEFINED ${_feature}_LIBRARIES)
    #
    # Add ${_feature}_LIBRARIES to
    #   DEAL_II_LIBRARIES
    #   DEAL_II_LIBRARIES_DEBUG
    #   DEAL_II_LIBRARIES_RELEASE
    # depending on the "optimized", "debug" or "general" keyword
    #
    set(_toggle "general")
    foreach(_tmp ${${_feature}_LIBRARIES})
      if( "${_tmp}" STREQUAL "debug" OR
          "${_tmp}" STREQUAL "optimized" OR
          "${_tmp}" STREQUAL "general" )
        set(_toggle "${_tmp}")
      else()
        if("${_toggle}" STREQUAL "general")
          list(APPEND DEAL_II_LIBRARIES ${_tmp})
        elseif("${_toggle}" STREQUAL "debug")
          list(APPEND DEAL_II_LIBRARIES_DEBUG ${_tmp})
        elseif("${_toggle}" STREQUAL "optimized")
          list(APPEND DEAL_II_LIBRARIES_RELEASE ${_tmp})
        endif()
      endif()
    endforeach()
  endif()

  foreach(_var ${DEAL_II_LIST_SUFFIXES})
    if(NOT "${_var}" STREQUAL "LIBRARIES" AND DEFINED ${_feature}_${_var})
      list(APPEND DEAL_II_${_var} ${${_feature}_${_var}})
    endif()
  endforeach()

  foreach(_var ${DEAL_II_STRING_SUFFIXES})
    if(DEFINED ${_feature}_${_var})
      add_flags(DEAL_II_${_var} "${${_feature}_${_var}}")
    endif()
  endforeach()

endmacro()
