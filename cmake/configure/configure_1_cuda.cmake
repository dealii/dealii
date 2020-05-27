## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2019 by the deal.II authors
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
# Configuration for cuda support:
#

#
# cuda support is experimental. Therefore, disable the feature per default:
#
SET(DEAL_II_WITH_CUDA FALSE CACHE BOOL "")

MACRO(FEATURE_CUDA_FIND_EXTERNAL var)

  # We need to set CUDA_USE_STATIC_CUDA_RUNTIME before FIND_PACKAGE(CUDA) and to
  # force the value otherwise it is overwritten by FIND_PACKAGE(CUDA)
  IF(BUILD_SHARED_LIBS)
    SET(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "" FORCE)
  ENDIF()

  #
  # TODO: Ultimately, this find_package call is not needed any more. We
  # still use it because it is very convenient to (a) check that CUDA is
  # installed, (b) get compiler path and include directories / libraries.
  #
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

    #
    # disable CUDA support older than 9.0:
    #
    IF(CUDA_VERSION_MAJOR VERSION_LESS 9.0)
      MESSAGE(FATAL_ERROR "\n"
        "deal.II requires CUDA version 9 or newer."
      )
    ENDIF()

    #
    # CUDA Toolkit 9 and CUDA Toolkit 10 are incompatible with C++17.
    # Make sure that deal.II is configured appropriately
    #
    MACRO(_cuda_ensure_feature_off _version _cpp_version_bad _cpp_version_good)
      IF(${CUDA_VERSION_MAJOR} EQUAL ${_version})
        IF(${DEAL_II_HAVE_CXX${_cpp_version_bad}})
          SET(${var} FALSE)
          MESSAGE(STATUS "CUDA ${_version} requires ${_feature} to be set to off.")
          SET(CUDA_ADDITIONAL_ERROR_STRING
            ${CUDA_ADDITIONAL_ERROR_STRING}
            "CUDA ${_version} is not compatible with the C++${_cpp_version_bad} standard.\n"
            "Please explicitly set the standard version to C++${_cpp_version_good}, e.g. by reconfiguring with\n"
            "  cmake -DDEAL_II_CXX_FLAGS=\"-std=c++${_cpp_version_good}\" ."
            )
        ENDIF()
      ENDIF()
    ENDMACRO()

    _cuda_ensure_feature_off(9 17 14)
    _cuda_ensure_feature_off(10 17 14)

    IF("${DEAL_II_CUDA_FLAGS_SAVED}" MATCHES "-arch[ ]*sm_([0-9]*)")
      SET(CUDA_COMPUTE_CAPABILITY "${CMAKE_MATCH_1}")
    ELSEIF("${DEAL_II_CUDA_FLAGS_SAVED}" MATCHES "-arch=sm_([0-9]*)")
      SET(CUDA_COMPUTE_CAPABILITY "${CMAKE_MATCH_1}")
    ELSEIF(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)
      #
      # Try to autodetect the CUDA Compute Capability by asking the device
      #
      SET(_binary_test_dir ${CMAKE_CURRENT_BINARY_DIR}/cmake/configure/CUDAComputeCapabilityWorkdir)
      FILE(REMOVE_RECURSE ${_binary_test_dir})
      FILE(MAKE_DIRECTORY ${_binary_test_dir})

      EXECUTE_PROCESS(
        COMMAND ${CUDA_NVCC_EXECUTABLE}
          -ccbin=${CMAKE_CXX_COMPILER}
          ${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure/CUDAComputeCapability/cuda_compute_capability.cu
          -o cuda_compute_capability
        WORKING_DIRECTORY ${_binary_test_dir}
        OUTPUT_QUIET
        ERROR_QUIET
        )
      EXECUTE_PROCESS(COMMAND ${_binary_test_dir}/cuda_compute_capability
                      RESULT_VARIABLE _result
                      OUTPUT_VARIABLE CUDA_COMPUTE_CAPABILITY)
      IF(${_result} EQUAL 0)
        ADD_FLAGS(DEAL_II_CUDA_FLAGS "-arch=sm_${CUDA_COMPUTE_CAPABILITY}")
        MESSAGE(STATUS "Detected CUDA Compute Capability ${CUDA_COMPUTE_CAPABILITY}")
      ELSE()
        MESSAGE(STATUS "Couldn't detect CUDA Compute Capability! "
                       "The error message was: ${CUDA_COMPUTE_CAPABILITY}")
        SET(CUDA_ADDITIONAL_ERROR_STRING
          ${CUDA_ADDITIONAL_ERROR_STRING}
          "Couldn't detect CUDA Compute Capability! "
          "The error message was: ${CUDA_COMPUTE_CAPABILITY}\n"
          "Please check the return value of ${_binary_test_dir}/cuda_compute_capability.\n"
          "If you want to disable the autodetection, set the compute capability to be used manually."
          )
        SET(${var} FALSE)
      ENDIF()
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

    # cuSOLVER requires OpenMP
    FIND_PACKAGE(OpenMP REQUIRED)
    SET(DEAL_II_LINKER_FLAGS "${DEAL_II_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")

    ADD_FLAGS(DEAL_II_CUDA_FLAGS_DEBUG "-G")
  ENDIF()
ENDMACRO()


MACRO(FEATURE_CUDA_CONFIGURE_EXTERNAL)

  #
  # Ensure that we enable CMake-internal CUDA support with the right
  # compiler:
  #
  SET(CMAKE_CUDA_COMPILER "${CUDA_NVCC_EXECUTABLE}")
  SET(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
  ENABLE_LANGUAGE(CUDA)

  MARK_AS_ADVANCED(CMAKE_CUDA_HOST_COMPILER)

  #
  # Work around a cmake 3.10 bug, see https://gitlab.kitware.com/cmake/cmake/issues/17797
  # because make does not support rsp link commands
  #
  SET(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_INCLUDES 0)
  SET(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_LIBRARIES 0)
  SET(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_OBJECTS 0)

  #
  # Set up cuda flags:
  #
  ADD_FLAGS(DEAL_II_CUDA_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

  # We cannot use -pedantic as compiler flags. nvcc generates code that
  # produces a lot of warnings when pedantic is enabled. So filter out the
  # flag:
  #
  STRING(REPLACE "-pedantic" "" DEAL_II_CXX_FLAGS "${DEAL_II_CXX_FLAGS}")

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
