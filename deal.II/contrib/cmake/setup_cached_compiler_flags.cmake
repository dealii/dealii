#
# FAT NOTE:
#
# We have a problem: We would like to setup our choice of C_FLAGS and
# CXX_FLAGS but let the user overwrite it (if desired).
#
# Unfortunately this is not as easy as it sounds:
#
# We have to call PROJECT(deal.II) in order to set up the project and run
# the C- and CXX-compiler detection and configuration. (Otherwise we cannot
# compile anything and have no idea which compiler is selected.) In this
# process CMAKE_{C|CXX}_FLAGS are already set to stupid default values.
# (And _cannot_ be sanely set from this script afterwards...)
#
# To mitigate this problem, we do the following: We set cached
# CMAKE_{C|CXX}_FLAGS variables with empty strings prior to initializing
# the compiler, so that no default values are set.
#
# The compiler-flag setup in setup_compiler_flags.cmake later adds our
# (target and platform dependend) default choice of flags.
#
# We add the cached variables _to the end_ of this default choice to allow
# the user to overwrite our choice of compiler flags.
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
# Set cached compiler flags to an empty string:
#

SET(deal_ii_used_flags
  CMAKE_C_FLAGS
  CMAKE_CXX_FLAGS
  CMAKE_C_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_C_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_DEBUG
  )

FOREACH(flags ${deal_ii_used_flags})
  SET(${flags} "" CACHE STRING
   "The user supplied cache variable will be appended _at the end_ of the auto generated ${flags} variable"
   )

 #
 # Save the cached variable at this point and clear it.
 # ${flags}_SAVED will be appended to ${flags} in
 # setup_cached_compiler_flags_finalize.cmake (called at the end of the
 # main CMakeLists.txt file).
 #
 SET(${flags}_SAVED "${${flags}}")
 SET(${flags} "")
ENDFOREACH()
