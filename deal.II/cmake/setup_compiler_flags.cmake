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
#   DEAL_II_LINKER_FLAGS
#   DEAL_II_LINKER_FLAGS_DEBUG
#   DEAL_II_LINKER_FLAGS_RELEASE
#
# All modifications shall be guarded with the ENABLE_IF_SUPPORTED
# or ENABLE_IF_LINKS macro, e.g.
#
#   ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-fpic")
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


###########################################################################
#                                                                         #
#                             Sanity checks:                              #
#                                                                         #
###########################################################################

#
# Check the user provided CXX flags:
#
#
# Omit to test for a sane C and Fortran toolchain for now..
#
IF(NOT "${CMAKE_CXX_FLAGS_SAVED}" STREQUAL "${DEAL_II_CXX_FLAGS_SAVED}")
  UNSET(DEAL_II_HAVE_USABLE_CXX_FLAGS CACHE)
ENDIF()
SET(DEAL_II_CXX_FLAGS_SAVED "${CMAKE_CXX_FLAGS_SAVED}" CACHE INTERNAL "" FORCE)

SET(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS_SAVED}")
CHECK_CXX_SOURCE_COMPILES(
  "int main(){ return 0; }"
  DEAL_II_HAVE_USABLE_CXX_FLAGS)
SET(CMAKE_REQUIRED_FLAGS "")

IF(NOT DEAL_II_HAVE_USABLE_CXX_FLAGS)
  UNSET(DEAL_II_HAVE_USABLE_CXX_FLAGS CACHE)
  MESSAGE(FATAL_ERROR "\n"
    "Configuration error: Cannot compile with the specified CXX flags: "
    "${CMAKE_CXX_FLAGS_SAVED}\n"
    )
ENDIF()


###########################################################################
#                                                                         #
#                            Compiler setup:                              #
#                                                                         #
###########################################################################

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
