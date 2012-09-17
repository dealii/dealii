#
# FAT NOTE:
#
# We set up the default compiler configuration, i.e. CMAKE_BUILD_TYPE,
# CMAKE_C_FLAGS, CMAKE_CXX_FLAGS as CACHED variables, so that the user can
# actually see and change them.
#
# Beware of the fact that compiler flags set later in the configuration via
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} <...>" are UNCACHED variables.
#
# cmake later concatenates cached and uncached variables together. So for a
# release build we will end up with (more or less) something like the
# following:
#
# $ <compiler>  <cached CMAKE_CXX_FLAGS> <uncached CMAKE_CXX_FLAGS> \
#      <cached CMAKE_CXX_FLAGS_RELEASE> <uncached CMAKE_CXX_FLAGS_RELEASE>\
#      <Definitions> <...>     <...>.c     -o <...>.o
#


#
# Setup CMAKE_BUILD_TYPE:
#

SET(CMAKE_BUILD_TYPE
  "Release"
  CACHE STRING
  "Choose the type of build, options are: Debug Release.")

IF( NOT CMAKE_BUILD_TYPE MATCHES "Release" AND
    NOT CMAKE_BUILD_TYPE MATCHES "Debug" )

  MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE does neither match Release, nor Debug.")

ENDIF()


#
# Setup default compiler flags:
#

# TODO: Copy the CFLAGS and CXXFLAGS logic from the former build system

SET(CMAKE_C_FLAGS
  "-Wfatal-errors -D_REENTRANT -fPIC -O2 -march=native"
  CACHE STRING
  "Flags used by the compiler during all build types.")

SET(CMAKE_CXX_FLAGS
  "-Wfatal-errors -D_REENTRANT -fPIC -O2 -march=native"
  CACHE STRING
  "Flags used by the compiler during all build types.")


SET(CMAKE_C_FLAGS_RELEASE
  "-O2"
  CACHE STRING
  "Flags used by the compiler during all release builds.")
SET(CMAKE_CXX_FLAGS_RELEASE
  "-O2"
  CACHE STRING
  "Flags used by the compiler during all release builds.")


SET(CMAKE_C_FLAGS_DEBUG
  "-O0 -ggdb -DDEBUG"
  CACHE STRING
  "Flags used by the compiler during all debug builds.")
SET(CMAKE_CXX_FLAGS_DEBUG
  "-O0 -ggdb -DDEBUG"
  CACHE STRING
  "Flags used by the compiler during all debug builds.")
