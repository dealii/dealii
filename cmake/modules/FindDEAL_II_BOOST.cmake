## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2020 by the deal.II authors
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
set(BOOST_VERSION_REQUIRED 1.59)

set(BOOST_DIR "" CACHE PATH "An optional hint to a BOOST installation")
set_if_empty(BOOST_DIR "$ENV{BOOST_DIR}")

if(NOT "${BOOST_DIR}" STREQUAL "")
  set(BOOST_ROOT "${BOOST_DIR}")
endif()

#
# Prefer static libs if BUILD_SHARED_LIBS=OFF:
#
if(NOT BUILD_SHARED_LIBS)
  set(Boost_USE_STATIC_LIBS TRUE)
endif()

# Work around a CMake compatibility issue with boost-1.70.0
# compare https://gitlab.kitware.com/cmake/cmake/issues/18865
# and https://lists.boost.org/Archives/boost/2019/02/245016.php
set(Boost_NO_BOOST_CMAKE ON)

set(Boost_NO_WARN_NEW_VERSIONS TRUE)
if(DEAL_II_WITH_ZLIB)
  find_package(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS
    iostreams serialization system thread
    )
else()
  find_package(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS
    serialization system thread
    )
endif()

#
# Fall back to dynamic libraries if no static libraries could be found:
#
if(NOT Boost_FOUND AND Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS FALSE)

  # temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  if(DEAL_II_WITH_ZLIB)
    find_package(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS iostreams serialization system thread)
  else()
    find_package(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS serialization system thread)
  endif()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
endif()

unset(Boost_NO_BOOST_CMAKE)

if(Boost_FOUND)
  set(BOOST_VERSION_MAJOR "${Boost_MAJOR_VERSION}")
  set(BOOST_VERSION_MINOR "${Boost_MINOR_VERSION}")
  set(BOOST_VERSION_SUBMINOR "${Boost_SUBMINOR_VERSION}")
  set(BOOST_VERSION
    "${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_SUBMINOR}"
    )
endif()

process_feature(BOOST
  LIBRARIES REQUIRED Boost_LIBRARIES
  INCLUDE_DIRS REQUIRED Boost_INCLUDE_DIRS
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
