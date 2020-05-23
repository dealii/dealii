## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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

#
# We require at least boost 1.59.
# - Boost::container::small_vector was introduced in 1.58 and some
#   serialization bugs in 1.58 were not fixed until 1.59.
#
SET(BOOST_VERSION_REQUIRED 1.59)

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

# temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

# Work around a CMake compatibility issue with boost-1.70.0
# compare https://gitlab.kitware.com/cmake/cmake/issues/18865
# and https://lists.boost.org/Archives/boost/2019/02/245016.php
SET(Boost_NO_BOOST_CMAKE ON)

IF(DEAL_II_WITH_ZLIB)
  FIND_PACKAGE(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS
    iostreams serialization system thread
    )
ELSE()
  FIND_PACKAGE(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS
    serialization system thread
    )
ENDIF()
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Fall back to dynamic libraries if no static libraries could be found:
#
IF(NOT Boost_FOUND AND Boost_USE_STATIC_LIBS)
  SET(Boost_USE_STATIC_LIBS FALSE)

  # temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
  LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  IF(DEAL_II_WITH_ZLIB)
    FIND_PACKAGE(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS iostreams serialization system thread)
  ELSE()
    FIND_PACKAGE(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS serialization system thread)
  ENDIF()
  LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
ENDIF()

UNSET(Boost_NO_BOOST_CMAKE)

IF(Boost_FOUND)
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
    Boost_DIR Boost_INCLUDE_DIRS Boost_IOSTREAMS_LIBRARY_DEBUG
    Boost_IOSTREAMS_LIBRARY_RELEASE Boost_LIBRARY_DIR
    Boost_SERIALIZATION_LIBRARY_DEBUG Boost_SERIALIZATION_LIBRARY_RELEASE
    Boost_SYSTEM_LIBRARY_DEBUG Boost_SYSTEM_LIBRARY_RELEASE
    Boost_THREAD_LIBRARY_DEBUG Boost_THREAD_LIBRARY_RELEASE
    Boost_LIBRARY_DIR_DEBUG Boost_LIBRARY_DIR_RELEASE
    _Boost_COMPONENTS_SEARCHED _Boost_INCLUDE_DIR_LAST
    _Boost_LIBRARY_DIR_LAST _Boost_USE_MULTITHREADED_LAST
    BOOST_IOSTREAMS_USABLE # clean up check in configure_boost.cmake
    BOOST_SERIALIZATION_USABLE # clean up check in configure_boost.cmake
  )
