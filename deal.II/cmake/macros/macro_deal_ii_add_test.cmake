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
# Usage:
#     DEAL_II_ADD_TEST(category test_name [configurations])
#

MACRO(DEAL_II_ADD_TEST _category _test_target)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})

    ITEM_MATCHES(_match "${_build}" ${ARGN})
    IF(_match OR "${ARGN}" STREQUAL "")

      STRING(TOLOWER ${_build} _build_lowercase)
      ADD_EXECUTABLE(${_test_target}.${_build_lowercase}
        EXCLUDE_FROM_ALL
        ${_test_target}.cc
        )
      SET_TARGET_PROPERTIES(${_test_target}.${_build_lowercase} PROPERTIES
        LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
        COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
        COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${_build}}"
        LINKER_LANGUAGE "CXX"
        )
      SET_PROPERTY(TARGET ${_test_target}.${_build_lowercase} APPEND PROPERTY
        INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/include"
          "${CMAKE_SOURCE_DIR}/include"
        )
      TARGET_LINK_LIBRARIES(${_test_target}.${_build_lowercase}
        ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
        )

      ADD_TEST(${_test_target}.${_build_lowercase}.build
        ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${_test_target}.${_build_lowercase}
        )

    ENDIF()

  ENDFOREACH()

ENDMACRO()
