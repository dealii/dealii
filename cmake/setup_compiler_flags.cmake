## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
#
# All modifications shall be guarded with the ENABLE_IF_SUPPORTED
# or ENABLE_IF_LINKS macro, e.g.
#
#   ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-fpic")
#   ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")
#
# Checks for compiler features (such as C++11 support) and compiler
# specific bugs that
#   - usually set up further configuration (such as preprocessor
#     definitions)
#   - disable a specific flag for a specific compiler version.
#
# belong the corresponding file:
#
#   ./cmake/checks/check_01_compiler_features.cmake
#   ./cmake/checks/check_01_cpu_features.cmake
#   ./cmake/checks/check_01_cxx_features.cmake
#   ./cmake/checks/check_01_system_features.cmake
#   ./cmake/checks/check_02_compiler_bugs.cmake
#


########################################################################
#                                                                      #
#                            Sanity checks:                            #
#                                                                      #
########################################################################

#
# Check the user provided CXX flags:
#

IF(NOT "${DEAL_II_CXX_FLAGS_SAVED}" STREQUAL "${CACHED_DEAL_II_CXX_FLAGS_SAVED}"
   OR NOT "${DEAL_II_LINKER_FLAGS_SAVED}" STREQUAL "${CACHED_DEAL_II_LINKER_FLAGS_SAVED}")
  MESSAGE(STATUS "")
  # Rerun this test if cxx flags changed:
  UNSET(DEAL_II_HAVE_USABLE_CXX_FLAGS CACHE)
ELSE()
  SET(DEAL_II_HAVE_USABLE_CXX_FLAGS TRUE CACHE INTERNAL "")
ENDIF()
SET(CACHED_DEAL_II_CXX_FLAGS_SAVED "${DEAL_II_CXX_FLAGS_SAVED}" CACHE INTERNAL "" FORCE)
SET(CACHED_DEAL_II_LINKER_FLAGS_SAVED "${DEAL_II_LINKER_FLAGS_SAVED}" CACHE INTERNAL "" FORCE)

# Initialize all CMAKE_REQUIRED_* variables a this point:
RESET_CMAKE_REQUIRED()

CHECK_CXX_SOURCE_COMPILES(
  "int main(){ return 0; }"
  DEAL_II_HAVE_USABLE_CXX_FLAGS)

IF(NOT DEAL_II_HAVE_USABLE_CXX_FLAGS)
  UNSET(DEAL_II_HAVE_USABLE_CXX_FLAGS CACHE)
  MESSAGE(FATAL_ERROR "
Configuration error: Cannot compile with the user supplied flags:
CXX flags: ${DEAL_II_CXX_FLAGS_SAVED}
LD flags: ${DEAL_II_LINKER_FLAGS_SAVED}
Please check the CMake variables DEAL_II_CXX_FLAGS, DEAL_II_LINKER_FLAGS
and the environment variables CXXFLAGS, LDFLAGS.\n\n"
    )
ENDIF()


########################################################################
#                                                                      #
#                           Compiler setup:                            #
#                                                                      #
########################################################################

IF(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  #
  # *Hooray* We are allowed to set compiler flags :-]
  #
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Setting up default compiler flags.")

  #
  # General setup for GCC and compilers sufficiently close to GCC:
  #
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
      CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    INCLUDE(setup_compiler_flags_gnu)
    SET(DEAL_II_KNOWN_COMPILER TRUE)
  ENDIF()

  #
  # Setup for ICC compiler (version >= 10):
  #
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    INCLUDE(setup_compiler_flags_intel)
    SET(DEAL_II_KNOWN_COMPILER TRUE)
  ENDIF()

  #
  # Setup for MSVC compiler (version >= 2012):
  #
   IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    INCLUDE(setup_compiler_flags_msvc)
    SET(DEAL_II_KNOWN_COMPILER TRUE)
  ENDIF()

  IF(NOT DEAL_II_KNOWN_COMPILER)
    MESSAGE(FATAL_ERROR "\n"
      "Unknown compiler!\n"
      "If you're serious about it, set DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS=OFF "
      "and set the relevant compiler options by hand.\n\n"
      )
  ENDIF()

ELSE(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)

  MESSAGE(STATUS "")
  MESSAGE(STATUS
    "Skipped setup of default compiler flags "
    "(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS=OFF)"
    )
ENDIF(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
