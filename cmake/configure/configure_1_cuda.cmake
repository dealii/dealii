## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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

#
# FindCUDA needs a compiler set up with C++11 support. Thus, only configure
# if deal.II was configured with C++11 support.
#
SET(FEATURE_CUDA_DEPENDS CXX11)


MACRO(FEATURE_CUDA_FIND_EXTERNAL var)

  #
  # FIXME Restructure this call into a ./modules/FindCUDA.cmake file and
  # use DEAL_II_PACKAGE_HANDLE
  #
  FIND_PACKAGE(CUDA)
  MARK_AS_ADVANCED(
    CUDA_HOST_COMPILER
    CUDA_SDK_ROOT_DIR
    CUDA_TOOLKIT_ROOT_DIR
    CUDA_USE_STATIC_CUDA_RUNTIME
    )

  IF(CUDA_FOUND)

    SET(${var} TRUE)

    IF(DEFINED CUDA_DIR)
      SET(CUDA_TOOLKIT_ROOT_DIR "${CUDA_DIR}")
    ENDIF()
    MESSAGE(STATUS "Configured to use CUDA installation at ${CUDA_TOOLKIT_ROOT_DIR}")

    IF("${CUDA_NVCC_FLAGS}" MATCHES "-arch")

      # Compute Capability specified explicitly.
      # Now parse:

      IF("${CUDA_NVCC_FLAGS}" MATCHES "-arch[ ]*sm_([0-9]*)")
        SET(CUDA_COMPUTE_CAPABILITY "${CMAKE_MATCH_1}")
      ELSEIF("${CUDA_NVCC_FLAGS}" MATCHES "-arch=sm_([0-9]*)")
        SET(CUDA_COMPUTE_CAPABILITY "${CMAKE_MATCH_1}")
      ELSE()
        STRING(REGEX MATCH "(-arch[ ]*[^ ]*)" match "${CUDA_NVCC_FLAGS}")
        MESSAGE(STATUS "Ill-formed Compute Capability specified.")
        SET(CUDA_ADDITIONAL_ERROR_STRING
          ${CUDA_ADDITIONAL_ERROR_STRING}
          "An ill-formed Compute Capability was passed in CUDA_NVCC_FLAGS: ${match}\n"
          "deal.II requires at least Compute Capability 3.5\n"
          "which is used as default is nothing is specified."
          )
        SET(${var} FALSE)
      ENDIF()


      IF("${CUDA_COMPUTE_CAPABILITY}" LESS "35")
        MESSAGE(STATUS "Too low CUDA Compute Capability specified -- deal.II requires at least 3.5 ")
        SET(CUDA_ADDITIONAL_ERROR_STRING
          ${CUDA_ADDITIONAL_ERROR_STRING}
          "Too low CUDA Compute Capability specified: ${CUDA_COMPUTE_CAPABILITY}\n"
          "deal.II requires at least Compute Capability 3.5\n"
          "which is used as default is nothing is specified."
          )
        SET(${var} FALSE)
      ENDIF()
    ENDIF()

    # Configuration was successful
    IF(${var})

      IF( NOT DEFINED CUDA_COMPUTE_CAPABILITY)
        # Set to use compute capability 3.5 by default
        SET(CUDA_COMPUTE_CAPABILITY "35")
        SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_35)
      ENDIF()

      # Export further definitions
      STRING(SUBSTRING "${CUDA_COMPUTE_CAPABILITY}" 0 1 CUDA_COMPUTE_CAPABILITY_MAJOR)
      STRING(SUBSTRING "${CUDA_COMPUTE_CAPABILITY}" 1 1 CUDA_COMPUTE_CAPABILITY_MINOR)
      SET(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE FALSE)

    ENDIF()
  ENDIF()

ENDMACRO()

MACRO(FEATURE_CUDA_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable cuda library!\n"
    ${CUDA_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that a cuda library is installed on your computer\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(CUDA)
