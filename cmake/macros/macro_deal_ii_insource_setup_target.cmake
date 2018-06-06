## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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
    LINKER_LANGUAGE "CXX"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

  TARGET_INCLUDE_DIRECTORIES(${_target}
    PRIVATE
      "${CMAKE_BINARY_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include"
      ${DEAL_II_BUNDLED_INCLUDE_DIRS}
    )
  TARGET_INCLUDE_DIRECTORIES(${_target} SYSTEM PRIVATE ${DEAL_II_INCLUDE_DIRS})

  IF(CMAKE_VERSION VERSION_LESS 3.9 OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY
      COMPILE_FLAGS " ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
      )
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
      COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS};${DEAL_II_USER_DEFINITIONS_${_build}}"
      )

  ELSE()

    SET(_flags "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}")
    SEPARATE_ARGUMENTS(_flags)
    TARGET_COMPILE_OPTIONS(${_target} PUBLIC
      $<$<COMPILE_LANGUAGE:CXX>:${_flags}>
      )

    TARGET_COMPILE_DEFINITIONS(${_target}
      PUBLIC ${DEAL_II_DEFINITIONS} ${DEAL_II_DEFINITIONS_${_build}}
      )

    IF(DEAL_II_WITH_CUDA)
      #
      # Add cxx compiler and cuda compilation flags to cuda source files:
      #

      SET(_cuda_flags "${DEAL_II_CUDA_FLAGS} ${DEAL_II_CUDA_FLAGS_${_build}}")
      SEPARATE_ARGUMENTS(_cuda_flags)

      #
      # Workaround: cuda will split every compiler option with a comma
      # (','), so remove all compiler flags that contain a comma:
      #
      STRING(REGEX REPLACE "[^ ]*,[^ ]*" "" _cxx_flags
        "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
        )

      TARGET_COMPILE_OPTIONS(${_target} PUBLIC
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler ${_cxx_flags}>
        $<$<COMPILE_LANGUAGE:CUDA>:${_cuda_flags}>
        )

      SET_TARGET_PROPERTIES(${_target} PROPERTIES
        CUDA_SEPARABLE_COMPILATION FALSE
        )
    ENDIF()
  ENDIF()

GET_PROPERTY(_type TARGET ${_target} PROPERTY TYPE)
IF(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
  TARGET_LINK_LIBRARIES(${_target}
    ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
    )
ENDIF()

ENDMACRO()
