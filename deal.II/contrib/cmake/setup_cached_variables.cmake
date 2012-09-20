#
# Setup cached variables prior to the PROJECT(deal.II) call
#
#
# We have a problem: We would like to setup our choice of C_FLAGS and
# CXX_FLAGS but let the user overwrite it (if desired).
#
# Unfortunately this is not as easy as it sounds:
#
# We have to call PROJECT(deal.II) in order to set up the project and run
# the C- and CXX-compiler detection and configuration. (Otherwise we cannot
# compile anything and, furthermore,  have no idea which compiler is
# selected.) In this process CMAKE_CXX_FLAGS[...] are already set to stupid
# default values.
# (And _cannot_ be sanely set from this script afterwards...)
#
# To mitigate this problem, we do the following:
#
#   - We initialize the cached CMAKE_CXX_FLAGS[...] variables with empty
#     strings prior to initializing the compiler, so that no default values
#     are set.
#
#   - We save the cached variables (possibly altered by the user via command
#     line or ccmake) in <variable>_SAVED and set <variable> to an empty
#     string.
#
#   - This way we can happily set our default flags in
#     setup_compiler_flags.cmake (if the user lets us, see
#     DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
#
#   - At the end of the configuration step we add <variables>_SAVED
#     _AT THE END_ of the respective <variable> allowing the user to
#     effectively overwrite our default settings.
#


#
# Setup CMAKE_BUILD_TYPE:
#
SET(CMAKE_BUILD_TYPE
  "Release"
  CACHE STRING
  "Choose the type of build, options are: Debug Release.")


#
# This is cruel, I know. But it is better to only have a known number of
# options for CMAKE_BUILD_TYPE...
#
IF( NOT CMAKE_BUILD_TYPE MATCHES "Release" AND
    NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE does neither match Release, nor Debug.")
ENDIF()


#
# Set BUILD_SHARED_LIBS to default to ON:
#
SET(BUILD_SHARED_LIBS ON CACHE OPTION "Build a shared library")


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

  # "CACHE" ensures that we only set the variable if it is not already set
  # as a  cached variable, effectively we're setting a default value:
  SET(${flags} "" CACHE STRING
   "The user supplied cache variable will be appended _at the end_ of the auto generated ${flags} variable"
   )

  #
  # Save the initial (cached) variable at this point and clear it.
  # ${flags}_SAVED will be appended to ${flags} in
  # setup_cached_compiler_flags_finalize.cmake (called at the end of the
  # main CMakeLists.txt file).
  #
  SET(${flags}_SAVED "${${flags}}")
  SET(${flags} "")

ENDFOREACH()

