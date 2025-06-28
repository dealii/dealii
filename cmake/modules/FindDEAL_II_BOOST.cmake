## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
# We require at least boost 1.74.
#
set(BOOST_VERSION_REQUIRED 1.74)

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

if(Boost_FOUND)
  set(BOOST_VERSION_MAJOR "${Boost_MAJOR_VERSION}")
  set(BOOST_VERSION_MINOR "${Boost_MINOR_VERSION}")
  set(BOOST_VERSION_SUBMINOR "${Boost_SUBMINOR_VERSION}")
  set(BOOST_VERSION
    "${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_SUBMINOR}"
    )
endif()

process_feature(BOOST
  TARGETS REQUIRED Boost_LIBRARIES
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
