#####
##
## Copyright (C) 2012 by the deal.II authors
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
# Setup cached variables prior to the PROJECT(deal.II) call
#


#
# Setup CMAKE_BUILD_TYPE:
#
SET(CMAKE_BUILD_TYPE
  "DebugRelease"
  CACHE STRING
  "Choose the type of build, options are: Debug, Release and DebugRelease.")


#
# This is cruel, I know. But it is better to only have a known number of
# options for CMAKE_BUILD_TYPE...
#
IF( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  MESSAGE(FATAL_ERROR
    "CMAKE_BUILD_TYPE does neither match Release, Debug, nor DebugRelease!"
    )
ENDIF()


#
# Set BUILD_SHARED_LIBS to default to ON and promote to cache so that the
# user can see the value.
#
SET(BUILD_SHARED_LIBS "ON" CACHE BOOL
  "Build a shared library"
  )
MARK_AS_ADVANCED(BUILD_SHARED_LIBS)


#
# Set CMAKE_INSTALL_RPATH_USE_LINK_PATH to default to ON and promote to
# cache so that the user can see the value.
#
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE BOOL
  "Set the rpath of the library to the external link pathes on installation"
  )
MARK_AS_ADVANCED(CMAKE_INSTALL_RPATH_USE_LINK_PATH)


#
# Tell the user very prominently, that we're doing things differently w.r.t
# CMAKE_(C|CXX)_FLAGS_(DEBUG|RELEASE)
#
SET(flags C_FLAGS_RELEASE CXX_FLAGS_RELEASE C_FLAGS_DEBUG CXX_FLAGS_DEBUG)
FOREACH(flag ${flags})
  IF(NOT "${CMAKE_${flag}}" STREQUAL "")
    MESSAGE(FATAL_ERROR
      "\nThe deal.II cmake build system does not use CMAKE_${flag}.\n"
      "Use DEAL_II_${flag}, instead!\n\n"
      )
  ENDIF()
ENDFOREACH()


#
# Hide all unused compiler flag variables:
#
SET(flags
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELWITHDEBINFO
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELWITHDEBINFO
  )
FOREACH(flag ${flags})
  #
  # Go away...
  #
  SET(${flag} "" CACHE INTERNAL "" FORCE)
ENDFOREACH()


#
# Read in CFLAGS, CXXFLAGS and LDFLAGS from environment
#
SET_IF_EMPTY(CMAKE_C_FLAGS "$ENV{CFLAGS}")
SET_IF_EMPTY(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS}")
SET_IF_EMPTY(CMAKE_SHARED_LINKER_FLAGS "$ENV{LDFLAGS}")


#
# Set cached compiler flags to an empty string:
#
SET(deal_ii_used_flags
  CMAKE_C_FLAGS
  CMAKE_CXX_FLAGS
  CMAKE_SHARED_LINKER_FLAGS
  DEAL_II_C_FLAGS_DEBUG
  DEAL_II_CXX_FLAGS_DEBUG
  DEAL_II_SHARED_LINKER_FLAGS_DEBUG
  DEAL_II_C_FLAGS_RELEASE
  DEAL_II_CXX_FLAGS_RELEASE
  DEAL_II_SHARED_LINKER_FLAGS_RELEASE
  )
FOREACH(flag ${deal_ii_used_flags})

  #
  # "CACHE" ensures that we only set the variable if it is not already set
  # as a  cached variable. Effectively we're setting a default value:
  #
  SET(${flag} "" CACHE STRING
   "The user supplied cache variable will be appended _at the end_ of the auto generated ${flag} variable"
   )

  #
  # Save the initial (cached) variable at this point and clear it.
  # ${flags}_SAVED will be appended to ${flags} in
  # setup_finalize.cmake (called at the end of the
  # main CMakeLists.txt file).
  #
  SET(${flag}_SAVED "${${flag}}")
  SET(${flag} "")

  #
  # Mark these flags as advanced.
  #
  MARK_AS_ADVANCED(${flag})
ENDFOREACH()

