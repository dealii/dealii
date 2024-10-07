## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
# Checks for compiler features (such as C++17 support) and compiler
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
  set(CMAKE_TRY_COMPILE_CONFIGURATION ${build})
  check_compiler_setup(
    "${DEAL_II_CXX_FLAGS_SAVED} ${DEAL_II_CXX_FLAGS_${build}_SAVED}"
    "${DEAL_II_LINKER_FLAGS_SAVED} ${DEAL_II_LINKER_FLAGS_${build}_SAVED}"
    DEAL_II_HAVE_USABLE_USER_FLAGS_${build}
    ${DEAL_II_LIBRARIES} ${DEAL_II_LIBRARIES_${build}}
    )
  unset(CMAKE_TRY_COMPILE_CONFIGURATION)

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
    CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
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
