#
# Setup default compiler flags:
#
# This file sets up sensible default compiler flags for the various
# platforms, compilers and build targets supported by the deal.II
# library.
#

#
#   ####################
#   #     FAT NOTE:    #
#   ####################
#
# As a rule of thumb: All configuration in setup_compiler_flags.cmake and
# setup_compiler_flags_<compiler>.cmake shall ONLY consist of CFLAGS,
# CXXFLAGS and LINKER_FLAGS being set.
#
# Checks for compiler features (such as C++11 support) and compiler
# specific bugs that set up further configuration (such as definitions)
# belong to
#
# ./check/check_for_compiler_features.cmake
# ./check/check_for_compiler_bugs.cmake
# ./check/check_for_cxx_features.cmake
#
# TODO: There is a bit of ambiguity. Clarify with Wolfgang.
#


#
# For the moment we assume that CC and CXX are the same compilers.
# (We only need CC for the compilation of the bundled umfpack library.)
# So, give a prominent error message in case CC and CXX differ:
#
IF(NOT ( ${CMAKE_C_COMPILER_ID} STREQUAL ${CMAKE_CXX_COMPILER_ID} AND
         ${CMAKE_C_COMPILER_VERSION} STREQUAL ${CMAKE_CXX_COMPILER_VERSION} ) )
    MESSAGE(SEND_ERROR "\n"
      "Configuration error: The specified C and CXX compiler have to be the "
      "same, but found:\n"
      "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}\n"
      "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\n"
      )
ENDIF()


#
# General setup for GCC and compilers sufficiently close to GCC:
#

IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  INCLUDE(setup_compiler_flags_gnu)
  SET(DEAL_II_KNOWN_COMPILER TRUE)
ENDIF()


IF(NOT DEAL_II_KNOWN_COMPILER)
  MESSAGE(WARNING "\n"
    "Unrecognized compiler!\n\n"
    "Please set the relevant compiler options by hand.\n")
ENDIF()

