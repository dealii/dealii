## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# This file provides an insource version of the DEAL_II_SETUP_TARGET macro.
#
# Usage:
#       DEAL_II_INSOURCE_SETUP_TARGET(target build)
#
# This appends necessary include directories, linker flags, compile
# definitions and the deal.II library link interface to the given target.
#
#

MACRO(DEAL_II_INSOURCE_SETUP_TARGET _target _build)

  SET_TARGET_PROPERTIES(${_target} PROPERTIES
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    LINKER_LANGUAGE "CXX"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_test_short}"
    )
  SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
    INCLUDE_DIRECTORIES
      "${CMAKE_BINARY_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include/deal.II/"
    )

GET_PROPERTY(_type TARGET ${_target} PROPERTY TYPE)
IF(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
  TARGET_LINK_LIBRARIES(${_target}
    ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
    )
ENDIF()

ENDMACRO()
