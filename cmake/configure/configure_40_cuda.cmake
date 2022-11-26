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
set(DEAL_II_WITH_CUDA FALSE CACHE BOOL "")

macro(feature_cuda_find_external var)
  if(NOT Kokkos_ENABLE_CUDA)
    set(CUDA_ADDITIONAL_ERROR_STRING
      ${CUDA_ADDITIONAL_ERROR_STRING}
      "deal.II can only be compiled with Cuda support if Kokkos was built with Cuda support!"
      )
    set(${var} FALSE)
  else()
    # FIXME We need to also find and link with cuSolver and cuSparse even though
    # relying on Kokkos for linking with Cuda. That's why we keep the code below for now.

    # We need to set CUDA_USE_STATIC_CUDA_RUNTIME before find_package(CUDA) and to
    # force the value otherwise it is overwritten by find_package(CUDA)
    if(BUILD_SHARED_LIBS)
      set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "" FORCE)
    endif()

    #
    # TODO: Ultimately, this find_package call is not needed any more. We
    # still use it because it is very convenient to (a) check that CUDA is
    # installed, (b) get compiler path and include directories / libraries.
    #
    find_package(DEAL_II_CUDA)

    if(CUDA_FOUND)
      #
      # CUDA was found, check whether we can actually use it:
      #
      set(${var} TRUE)

      #
      # CUDA support requires CMake version 3.9 or newer
      #
      if(CMAKE_VERSION VERSION_LESS 3.9)
        set(${var} FALSE)
        message(STATUS "deal.II requires CMake version 3.9, or newer for CUDA support")
        set(CUDA_ADDITIONAL_ERROR_STRING
          ${CUDA_ADDITIONAL_ERROR_STRING}
          "deal.II requires CMake version 3.9, or newer for CUDA support.\n"
          "Reconfigure with a sufficient cmake version."
          )
      endif()

      #
      # disable CUDA support older than 10.2:
      #
      if(CUDA_VERSION VERSION_LESS 10.2)
        message(FATAL_ERROR "\n"
          "deal.II requires CUDA version 10.2 or newer."
        )
      endif()

      #
      # CUDA Toolkit 10 is incompatible with C++17.
      # Make sure that deal.II is configured appropriately
      #
      macro(_cuda_ensure_feature_off _version _cpp_version_bad _cpp_version_good)
        if(${CUDA_VERSION_MAJOR} EQUAL ${_version})
          if(${DEAL_II_HAVE_CXX${_cpp_version_bad}})
            set(${var} FALSE)
            message(STATUS "CUDA ${_version} requires ${_feature} to be set to off.")
            set(CUDA_ADDITIONAL_ERROR_STRING
              ${CUDA_ADDITIONAL_ERROR_STRING}
              "CUDA ${_version} is not compatible with the C++${_cpp_version_bad} standard.\n"
              "Please explicitly set the standard version to C++${_cpp_version_good}, e.g. by reconfiguring with\n"
              "  cmake -DDEAL_II_CXX_FLAGS=\"-std=c++${_cpp_version_good}\" ."
              )
          endif()
        endif()
      endmacro()

      _cuda_ensure_feature_off(10 17 14)

      # cuSOLVER requires OpenMP
      find_package(OpenMP REQUIRED)
      set(DEAL_II_LINKER_FLAGS "${DEAL_II_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
  endif()
endmacro()


macro(feature_cuda_configure_external)
  # We cannot use -pedantic as compiler flags. nvcc generates code that
  # produces a lot of warnings when pedantic is enabled. So filter out the
  # flag:
  #
  string(REPLACE "-pedantic" "" DEAL_II_CXX_FLAGS "${DEAL_II_CXX_FLAGS}")
endmacro()


macro(feature_cuda_error_message)
  message(FATAL_ERROR "\n"
    "Could not find any suitable cuda library!\n"
    ${CUDA_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that a cuda library is installed on your computer\n"
    )
endmacro()


configure_feature(CUDA)
