## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
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
# Setup default compiler flags: This file sets up sensible default compiler
# flags for the various platforms, compilers and build targets supported by
# the deal.II library.
#
#
# ####################
# #     FAT NOTE:    #
# ####################
#
# All configuration in setup_compiler_flags.cmake and
# setup_compiler_flags_<compiler>.cmake shall ONLY modify:
#
#   DEAL_II_CXX_FLAGS
#   DEAL_II_CXX_FLAGS_DEBUG
#   DEAL_II_CXX_FLAGS_RELEASE
#   DEAL_II_LINKER_FLAGS
#   DEAL_II_LINKER_FLAGS_DEBUG
#   DEAL_II_LINKER_FLAGS_RELEASE
#   DEAL_II_DEFINITIONS
#   DEAL_II_DEFINITIONS_DEBUG
#   DEAL_II_DEFINITIONS_RELEASE
#
# All modifications shall be guarded with the ENABLE_IF_SUPPORTED
# or ENABLE_IF_LINKS macro, e.g.
#
#   enable_if_supported(DEAL_II_CXX_FLAGS "-fpic")
#   enable_if_links(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")
#
# Checks for compiler features (such as C++14 support) and compiler
# specific bugs that
#   - usually set up further configuration (such as preprocessor
#     definitions)
#   - disable a specific flag for a specific compiler version.
#
# belong the corresponding file:
#
#   ./cmake/checks/check_01_cpu_features.cmake
#   ./cmake/checks/check_01_cxx_features.cmake
#   ./cmake/checks/check_02_compiler_features.cmake
#   ./cmake/checks/check_02_system_features.cmake
#   ./cmake/checks/check_03_compiler_bugs.cmake
#


########################################################################
#                                                                      #
#                            Sanity checks:                            #
#                                                                      #
########################################################################

#
# Check the user provided CXX flags:
#

foreach(build ${DEAL_II_BUILD_TYPES})
  check_compiler_setup(
    "${DEAL_II_CXX_FLAGS_SAVED} ${DEAL_II_CXX_FLAGS_${build}_SAVED}"
    "${DEAL_II_LINKER_FLAGS_SAVED} ${DEAL_II_LINKER_FLAGS_${build}_SAVED}"
    DEAL_II_HAVE_USABLE_USER_FLAGS_${build}
    ${DEAL_II_LIBRARIES} ${DEAL_II_LIBRARIES_${build}}
    )

  if(NOT DEAL_II_HAVE_USABLE_USER_FLAGS_${build})
    message(FATAL_ERROR "
  Configuration error: Cannot compile with the user supplied flags:
    CXX flags (${build}): ${DEAL_II_CXX_FLAGS_SAVED} ${DEAL_II_CXX_FLAGS_${build}_SAVED}
    LD flags  (${build}): ${DEAL_II_LINKER_FLAGS_SAVED} ${DEAL_II_LINKER_FLAGS_${build}_SAVED}
    LIBRARIES (${build}): ${DEAL_II_LIBRARIES};${DEAL_II_LIBRARIES_${build}}
  Please check the CMake variables
    DEAL_II_CXX_FLAGS, DEAL_II_CXX_FLAGS_${build},
    DEAL_II_LINKER_FLAGS, DEAL_II_CXX_FLAGS_${build}
  and the environment variables CXXFLAGS, LDFLAGS.\n\n"
      )
  endif()
endforeach()


########################################################################
#                                                                      #
#                           Compiler setup:                            #
#                                                                      #
########################################################################

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  #
  # General setup for GCC and compilers sufficiently close to GCC:
  #
  verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags_gnu.cmake)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  #
  # Setup for ICC compiler (version >= 10):
  #
  verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags_intel.cmake)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  #
  # Setup for MSVC compiler (version >= 2012):
  #
  verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags_msvc.cmake)

else()
  message(WARNING "\nUnknown compiler!\n"
    "Please populate the CMake variables DEAL_II_CXX_FLAGS(|DEBUG|RELEASE) "
    "and DEAL_II_LINKER_FLAGS(|DEBUG|RELEASE) as needed."
    )
endif()
