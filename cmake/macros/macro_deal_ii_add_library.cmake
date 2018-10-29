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
# A small wrapper around ADD_LIBRARY that will define a target for each
# build type specified in DEAL_II_BUILD_TYPES. Only compatible with object
# targets (as used in the build system).
#
# It is assumed that the desired compilation configuration is set via
#   DEAL_II_CXX_FLAGS_${build}
#   DEAL_II_DEFINITIONS_${build}
#
# as well as the global (for all build types)
#   DEAL_II_CXX_FLAGS
#   DEAL_II_DEFINITIONS
#

MACRO(DEAL_II_ADD_LIBRARY _library)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    ADD_LIBRARY(${_library}_${_build_lowercase}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${_library}_${_build_lowercase} PROPERTIES
      LINKER_LANGUAGE "CXX"
      )

    IF(CMAKE_VERSION VERSION_LESS 3.9 OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

      SET_TARGET_PROPERTIES(${_library}_${_build_lowercase} PROPERTIES
        COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
        COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
        )

    ELSE()

      SET(_flags "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}")
      SEPARATE_ARGUMENTS(_flags)
      TARGET_COMPILE_OPTIONS(${_library}_${_build_lowercase} PUBLIC
        $<$<COMPILE_LANGUAGE:CXX>:${_flags}>
        )

      TARGET_COMPILE_DEFINITIONS(${_library}_${_build_lowercase}
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

        TARGET_COMPILE_OPTIONS(${_library}_${_build_lowercase} PUBLIC
          $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler ${_cxx_flags}>
          $<$<COMPILE_LANGUAGE:CUDA>:${_cuda_flags}>
          )

        SET_TARGET_PROPERTIES(${_library}_${_build_lowercase} PROPERTIES
          CUDA_SEPARABLE_COMPILATION FALSE
          )
      ENDIF()

    ENDIF()

    SET_PROPERTY(GLOBAL APPEND PROPERTY DEAL_II_OBJECTS_${_build}
      "$<TARGET_OBJECTS:${_library}_${_build_lowercase}>"
      )
  ENDFOREACH()

ENDMACRO()
