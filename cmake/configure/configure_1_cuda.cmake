## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2017 by the deal.II authors
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
# Configuration for cuda support:
#

#
# cuda support is experimental. Therefore, disable the feature per default:
#
SET(DEAL_II_WITH_CUDA FALSE CACHE BOOL "")

MACRO(FEATURE_CUDA_FIND_EXTERNAL var)

  FIND_PACKAGE(CUDA)

  IF(CUDA_FOUND)
    #
    # CUDA was found, check whether we can actually use it:
    #
    SET(${var} TRUE)

    #
    # CUDA support requires CMake version 3.9 or newer
    #
    IF(CMAKE_VERSION VERSION_LESS 3.9)
      SET(${var} FALSE)
      MESSAGE(STATUS "deal.II requires CMake version 3.9, or newer for CUDA support")
      SET(CUDA_ADDITIONAL_ERROR_STRING
        ${CUDA_ADDITIONAL_ERROR_STRING}
        "deal.II requires CMake version 3.9, or newer for CUDA support.\n"
        "Reconfigure with a sufficient cmake version."
        )
    ENDIF()

    IF("${DEAL_II_CUDA_FLAGS}" MATCHES "-arch[ ]*sm_([0-9]*)")
      SET(CUDA_COMPUTE_CAPABILITY "${CMAKE_MATCH_1}")
    ELSEIF("${DEAL_II_CUDA_FLAGS}" MATCHES "-arch=sm_([0-9]*)")
      SET(CUDA_COMPUTE_CAPABILITY "${CMAKE_MATCH_1}")
    ELSE()
      #
      # Assume a cuda compute capability of 35
      #
      SET(CUDA_COMPUTE_CAPABILITY "35")
      ADD_FLAGS(DEAL_II_CUDA_FLAGS "-arch=sm_35")
    ENDIF()

    IF("${CUDA_COMPUTE_CAPABILITY}" LESS "35")
      MESSAGE(STATUS "Too low CUDA Compute Capability specified -- deal.II requires at least 3.5 ")
      SET(CUDA_ADDITIONAL_ERROR_STRING
        ${CUDA_ADDITIONAL_ERROR_STRING}
        "Too low CUDA Compute Capability specified: ${CUDA_COMPUTE_CAPABILITY}\n"
        "deal.II requires at least Compute Capability 3.5\n"
        "which is used as default if nothing is specified."
        )
      SET(${var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()


MACRO(FEATURE_CUDA_CONFIGURE_EXTERNAL)

  ENABLE_LANGUAGE(CUDA)

  #
  # Set up cuda flags:
  #
  ADD_FLAGS(DEAL_II_CUDA_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

  #
  # Export definitions:
  #
  STRING(SUBSTRING "${CUDA_COMPUTE_CAPABILITY}" 0 1 CUDA_COMPUTE_CAPABILITY_MAJOR)
  STRING(SUBSTRING "${CUDA_COMPUTE_CAPABILITY}" 1 1 CUDA_COMPUTE_CAPABILITY_MINOR)
ENDMACRO()


MACRO(FEATURE_CUDA_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable cuda library!\n"
    ${CUDA_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that a cuda library is installed on your computer\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(CUDA)
