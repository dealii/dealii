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
# All configuration in setup_compiler_flags.cmake and
# setup_compiler_flags_<compiler>.cmake shall ONLY consist of CFLAGS,
# CXXFLAGS and LINKER_FLAGS being set.
#
# Checks for compiler features (such as C++11 support) and compiler
# specific bugs that
#   - usually set up further configuration (such as definitions)
#   - disable a specific flag for a specific compiler version.
#
# belong to
#
# ./check/check_for_compiler_features.cmake
#
# ./check/check_for_compiler_bugs.cmake
#
# ./check/check_for_cxx_features.cmake
#


#
#   ######################
#   #     FAT NOTE 2:    #
#   ######################
#
# For the moment we assume that CC and CXX are the same compiler.
# We need CC for the compilation of the bundled umfpack library only...
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
# Setup CMAKE_CXX_FLAGS<_RELEASE|_DEBUG>:
#

IF(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  MESSAGE(STATUS "Set up default compiler flags.")

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

ELSE(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  MESSAGE(STATUS
    "Skipped setting up default compiler flags "
    "(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS=OFF)"
    )
ENDIF(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)



#
# Setup CMAKE_C_FLAGS<_RELEASE|_DEBUG>:
#

#
# For the moment we assume that CC and CXX are the same compiler and that
# we can set (almost) the same default flags for both.
#

SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})

#
# Strip flags not supported by a C target:
#
STRIP_FLAG(CMAKE_C_FLAGS "-Wsynth")
STRIP_FLAG(CMAKE_C_FLAGS_RELEASE "-felide-constructors")

#
# Disable some warnings:
#
STRIP_FLAG(CMAKE_C_FLAGS "-Wall") # There is no other way to disable -Wunknown-pragma atm...
STRIP_FLAG(CMAKE_C_FLAGS "-Wsign-compare")
STRIP_FLAG(CMAKE_C_FLAGS "-Wwrite-strings")
