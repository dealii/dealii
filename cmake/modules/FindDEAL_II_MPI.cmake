## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2014 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Find MPI
#
# This module exports:
#   MPI_CXX_FLAGS
#   MPI_INCLUDE_DIRS
#   MPI_LIBRARIES
#   MPI_LINKER_FLAGS
#   MPI_VERSION
#   MPI_VERSION_MAJOR
#   MPI_VERSION_MINOR
#   MPI_WITH_DEVICE_SUPPORT
#   OMPI_VERSION
#

#
# Configuration for mpi support:
#
# We look for the C and Fortran libraries as well because they are needed
# by some external libraries for the link interface.
#

if(MPI_CXX_FOUND)
  set(MPI_FOUND TRUE)
endif()

#
# Call the system FindMPI.cmake module:
#

#
# Make sure we pick up the correct MPI implementation for the case that
# environment variables MPIEXEC_EXECUTABLE, or MPIEXEC are set. If
# MPIEXEC_EXECUTABLE is already set as a CMake variable simply ignore the
# environment variables.
#
if(NOT MPIEXEC_EXECUTABLE)
  if(DEFINED ENV{MPIEXEC_EXECUTABLE})
    find_program(MPIEXEC_EXECUTABLE $ENV{MPIEXEC_EXECUTABLE})
  elseif(DEFINED ENV{MPIEXEC})
    find_program(MPIEXEC_EXECUTABLE $ENV{MPIEXEC})
  endif()
  # For backwards compatibility with old cmake versions:
  set(MPIEXEC "${MPIEXEC_EXECUTABLE}")
endif()

find_package(MPI)

#
# Newer versions of FindMPI.cmake only populate MPI_CXX_* (and MPI_C_*,
# MPI_Fortran_*) variables. Let's rename these version names
#

if(NOT DEFINED MPI_VERSION AND DEFINED MPI_CXX_VERSION)
  set(MPI_VERSION ${MPI_CXX_VERSION})
  set(MPI_VERSION_MAJOR ${MPI_CXX_VERSION_MAJOR})
  set(MPI_VERSION_MINOR ${MPI_CXX_VERSION_MINOR})
endif()

#
# Really old versions of CMake do not export any version information. In
# this case, query the mpi.h header for the necessary information:
#

deal_ii_find_file(MPI_MPI_H
  NAMES mpi.h
  HINTS ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH}
  )
if(NOT MPI_MPI_H MATCHES "-NOTFOUND" AND NOT DEFINED MPI_VERSION)
  file(STRINGS "${MPI_MPI_H}" MPI_VERSION_MAJOR_STRING
    REGEX "#define.*MPI_VERSION")
  string(REGEX REPLACE "^.*MPI_VERSION[ ]+([0-9]+).*" "\\1"
    MPI_VERSION_MAJOR "${MPI_VERSION_MAJOR_STRING}"
    )
  file(STRINGS ${MPI_MPI_H} MPI_VERSION_MINOR_STRING
    REGEX "#define.*MPI_SUBVERSION")
  string(REGEX REPLACE "^.*MPI_SUBVERSION[ ]+([0-9]+).*" "\\1"
    MPI_VERSION_MINOR "${MPI_VERSION_MINOR_STRING}"
    )
  set(MPI_VERSION "${MPI_VERSION_MAJOR}.${MPI_VERSION_MINOR}")
endif()

#
# Except - this doesn't always work. Some distributions install a header
# stub mpi.h that includes the right mpi header depending on the
# architecture. In this case we are really out of luck. It is not
# straightforward to find the correct header file to query the version
# information from. Just set a very conservative default:
#
if(NOT DEFINED MPI_VERSION OR MPI_VERSION STREQUAL ".")
  set(MPI_VERSION "0.0")
  set(MPI_VERSION_MAJOR "0")
  set(MPI_VERSION_MINOR "0")
endif()

#
# Check whether the MPI implementation is able to work with buffers that
# reside in device memory directly, instead of requiring an explicit copy to
# host memory first (known as "device-aware" or "CUDA-aware" MPI).
#
# There is no portable way of querying this at configure time. But the two
# major MPI implementations both advertise the capability with preprocessor
# symbols exported by the (non standard) mpi-ext.h header.
#
# Should the heuristic fail for a given MPI implementation the check can
# always be overridden on the command line by configuring with
# -DDEAL_II_MPI_WITH_DEVICE_SUPPORT=ON.
#

clear_cmake_required()
set(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS_SAVED} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})

unset_if_changed(CHECK_MPI_DEVICE_SUPPORT_SAVED
  "${MPI_CXX_COMPILER}${CMAKE_REQUIRED_INCLUDES}"
  MPI_WITH_DEVICE_SUPPORT
  )

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <mpi.h>
  #include <mpi-ext.h>

  #if !defined(MPIX_CUDA_AWARE_SUPPORT) || MPIX_CUDA_AWARE_SUPPORT == 0
  #  if !defined(MPIX_ROCM_AWARE_SUPPORT) || MPIX_ROCM_AWARE_SUPPORT == 0
  #    error \"This MPI library does not advertise device support.\"
  #  endif
  #endif

  int main()
  {
  }
  "
  MPI_WITH_DEVICE_SUPPORT)

reset_cmake_required()

#
# Make sure that we do not run into underlinking on Debian/Ubuntu systems with
# lld / ld.gold and missing libopen-pal.so on the link line:
#

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  check_compiler_setup(
    "${DEAL_II_CXX_FLAGS_SAVED} ${DEAL_II_CXX_FLAGS}"
    "${DEAL_II_LINKER_FLAGS_SAVED} ${DEAL_II_LINKER_FLAGS}"
    MPI_UNDERLINKAGE_OK
    ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES}
    )

  if(NOT MPI_UNDERLINKAGE_OK AND NOT "${MPI_CXX_LIBRARIES}" STREQUAL "")
    # This check only works if MPI_CXX_LIBRARIES is non-empty, otherwise we will just give up
    # and hope for the best...
    message(STATUS "Trying to avoid underlinkage by expliclitly adding libopen-pal to link line")

    list(GET MPI_CXX_LIBRARIES 0 _lib)
    get_filename_component(_hint ${_lib} DIRECTORY)
    deal_ii_find_library(_mpi_libopen_pal_library
      NAMES open-pal
      HINTS ${_hint}
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH
      )

    #
    # Note: We don't need to check whether the find library call is
    # successful: If libopen-pal cannot be found then the
    # process_feature will drop the library automatically.
    #
    # In this case the sanity check in cmake/setup_finalize.cmake will fail
    # and we start dropping -fuse-ld=lld and -fuse-ld=ld.gold from the
    # command line.
    #
  endif()
endif()

process_feature(MPI
  LIBRARIES
    OPTIONAL MPI_CXX_LIBRARIES MPI_Fortran_LIBRARIES MPI_C_LIBRARIES _mpi_libopen_pal_library
  INCLUDE_DIRS
    OPTIONAL MPI_CXX_INCLUDE_PATH MPI_C_INCLUDE_PATH
  CXX_FLAGS OPTIONAL MPI_CXX_COMPILE_FLAGS
  LINKER_FLAGS OPTIONAL MPI_CXX_LINK_FLAGS
  CLEAR
    CHECK_MPI_DEVICE_SUPPORT_SAVED
    MPI_WITH_DEVICE_SUPPORT
    MPI_C_COMPILER
    MPI_CXX_COMPILER
    MPIEXEC
    MPI_EXTRA_LIBRARY
    MPI_Fortran_COMPILER
    MPI_HEADER_PATH
    MPI_LIB
    MPI_LIBRARY
    MPI_MPI_H
  )
