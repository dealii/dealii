## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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
# A small wrapper around ADD_LIBRARY that will define a target for each
# build type specified in DEAL_II_BUILD_TYPES
#
# It is assumed that the desired compilation configuration is set via
#   DEAL_II_LINKER_FLAGS_${build}
#   DEAL_II_CXX_FLAGS_${build}
#   DEAL_II_DEFINITIONS_${build}
#
# as well as the global (for all build types)
#   DEAL_II_LINKER_FLAGS
#   DEAL_II_CXX_FLAGS
#   DEAL_II_DEFINITIONS
#

MACRO(DEAL_II_ADD_LIBRARY _library)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    ADD_LIBRARY(${_library}.${_build_lowercase}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${_library}.${_build_lowercase} PROPERTIES
      LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
      COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
      LINKER_LANGUAGE "CXX"
      )

    SET_PROPERTY(GLOBAL APPEND PROPERTY DEAL_II_OBJECTS_${_build}
      "$<TARGET_OBJECTS:${_library}.${_build_lowercase}>"
      )

    #
    # Cuda specific target setup:
    #
    IF(DEAL_II_WITH_CUDA)
      CUDA_WRAP_SRCS(${_library}.${_build_lowercase}
        OBJ _generated_cuda_files ${ARGN} SHARED
        )

      ADD_CUSTOM_TARGET(${_library}.${_build_lowercase}_cuda
        DEPENDS
        ${_generated_cuda_files}
        )
      ADD_DEPENDENCIES(${_library}.${_build_lowercase}
        ${_library}.${_build_lowercase}_cuda
        )

      SET_PROPERTY(GLOBAL APPEND PROPERTY DEAL_II_OBJECTS_${_build}
        "${_generated_cuda_files}"
        )
    ENDIF()

  ENDFOREACH()

ENDMACRO()
