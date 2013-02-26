#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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
#   CMAKE_CXX_FLAGS
#   DEAL_II_CXX_FLAGS_DEBUG
#   DEAL_II_CXX_FLAGS_RELEASE
#   CMAKE_SHARED_LINKER_FLAGS
#   DEAL_II_SHARED_LINKER_FLAGS_DEBUG
#   DEAL_II_SHARED_LINKER_FLAGS_RELEASE
#
#
# Checks for compiler features (such as C++11 support) and compiler
# specific bugs that
#   - usually set up further configuration (such as definitions)
#   - disable a specific flag for a specific compiler version.
#
# belong to
#
#   ./check/check_for_compiler_features.cmake
#   ./check/check_for_compiler_bugs.cmake
#   ./check/check_for_compiler_bugs_*.cmake
#   ./check/check_for_cxx_features.cmake
#


#
# CMAKE_C_COMPILER and CMAKE_CXX_COMPILER have to be of the same brand.
#
IF(NOT ( "${CMAKE_C_COMPILER_ID}" STREQUAL "${CMAKE_CXX_COMPILER_ID}" AND
         "${CMAKE_C_COMPILER_VERSION}" STREQUAL "${CMAKE_CXX_COMPILER_VERSION}" ) )
    MESSAGE(FATAL_ERROR "\n"
      "Configuration error: The specified C and CXX compiler have to be of the "
      "same family, but cmake found:\n"
      "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}\n"
      "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\n"
      )
ENDIF()


IF(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  #
  # *Hooray* We are allowed to set compiler flags :-]
  #
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Set up default compiler flags.")

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
