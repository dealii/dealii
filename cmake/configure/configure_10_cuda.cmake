## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2022 by the deal.II authors
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
    # disable CUDA support older than 10.2:
    #
    IF(CUDA_VERSION VERSION_LESS 10.2)
      MESSAGE(FATAL_ERROR "\n"
        "deal.II requires CUDA version 10.2 or newer."
      )
    ENDIF()

    #
    # CUDA Toolkit 10 is incompatible with C++17.
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

    _cuda_ensure_feature_off(10 17 14)

    # cuSOLVER requires OpenMP
    FIND_PACKAGE(OpenMP REQUIRED)
    SET(DEAL_II_LINKER_FLAGS "${DEAL_II_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDMACRO()


MACRO(FEATURE_CUDA_CONFIGURE_EXTERNAL)
  # We cannot use -pedantic as compiler flags. nvcc generates code that
  # produces a lot of warnings when pedantic is enabled. So filter out the
  # flag:
  #
  STRING(REPLACE "-pedantic" "" DEAL_II_CXX_FLAGS "${DEAL_II_CXX_FLAGS}")
ENDMACRO()


MACRO(FEATURE_CUDA_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable cuda library!\n"
    ${CUDA_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that a cuda library is installed on your computer\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(CUDA)
