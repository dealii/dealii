## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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
# Try to find the boost libraries
#
# This module exports:
#
#     BOOST_FOUND
#     BOOST_LIBRARIES
#     BOOST_INCLUDE_DIRS
#     BOOST_VERSION
#     BOOST_VERSION_MAJOR
#     BOOST_VERSION_MINOR
#     BOOST_VERSION_SUBMINOR
#

SET(BOOST_DIR "" CACHE PATH "An optional hint to a BOOST installation")
SET_IF_EMPTY(BOOST_DIR "$ENV{BOOST_DIR}")

IF(NOT "${BOOST_DIR}" STREQUAL "")
  SET(BOOST_ROOT "${BOOST_DIR}")
ENDIF()

#
# Prefer static libs if BUILD_SHARED_LIBS=OFF:
#
IF(NOT BUILD_SHARED_LIBS)
  SET(Boost_USE_STATIC_LIBS TRUE)
ENDIF()

LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
FIND_PACKAGE(Boost 1.48 COMPONENTS iostreams serialization system thread)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Fall back to dynamic libraries if no static libraries could be found:
#
IF(NOT Boost_FOUND AND Boost_USE_STATIC_LIBS)
  SET(Boost_USE_STATIC_LIBS FALSE)
  FIND_PACKAGE(Boost 1.48 COMPONENTS iostreams serialization system thread)
ENDIF()

IF(Boost_FOUND)
  #
  # Remove "pthread" from Boost_LIBRARIES. Threading, if necessary, is
  # already set up via configure_1_threads.cmake.
  #
  LIST(REMOVE_ITEM Boost_LIBRARIES "pthread")

  SET(BOOST_VERSION_MAJOR "${Boost_MAJOR_VERSION}")
  SET(BOOST_VERSION_MINOR "${Boost_MINOR_VERSION}")
  SET(BOOST_VERSION_SUBMINOR "${Boost_SUBMINOR_VERSION}")
  SET(BOOST_VERSION
    "${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_SUBMINOR}"
    )
ENDIF()

DEAL_II_PACKAGE_HANDLE(BOOST
  LIBRARIES REQUIRED Boost_LIBRARIES
  INCLUDE_DIRS REQUIRED Boost_INCLUDE_DIRS
  USER_INCLUDE_DIRS Boost_INCLUDE_DIRS
  CLEAR
    Boost_INCLUDE_DIR Boost_IOSTREAMS_LIBRARY_DEBUG
    Boost_IOSTREAMS_LIBRARY_RELEASE Boost_LIBRARY_DIR
    Boost_SERIALIZATION_LIBRARY_DEBUG Boost_SERIALIZATION_LIBRARY_RELEASE
    Boost_SYSTEM_LIBRARY_DEBUG Boost_SYSTEM_LIBRARY_RELEASE
    Boost_THREAD_LIBRARY_DEBUG Boost_THREAD_LIBRARY_RELEASE
    _Boost_COMPONENTS_SEARCHED _Boost_INCLUDE_DIR_LAST
    _Boost_LIBRARY_DIR_LAST _Boost_USE_MULTITHREADED_LAST
  )
