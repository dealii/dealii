## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2016 by the deal.II authors
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
# We require at least boost 1.59 since boost::container::small_vector was
# introduced in 1.58 and some serialization bugs in 1.58 were not fixed until
# 1.59.

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
FIND_PACKAGE(Boost 1.59 COMPONENTS
  iostreams serialization system thread
  )
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Fall back to dynamic libraries if no static libraries could be found:
#
IF(NOT Boost_FOUND AND Boost_USE_STATIC_LIBS)
  SET(Boost_USE_STATIC_LIBS FALSE)

  # temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
  LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  FIND_PACKAGE(Boost 1.59 COMPONENTS iostreams serialization system thread)
  LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
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

  IF(DEAL_II_WITH_ZLIB)
    #
    # Test that Boost.Iostreams is usuable.
    #
    RESET_CMAKE_REQUIRED()
    PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
    PUSH_CMAKE_REQUIRED("-L${Boost_LIBRARY_DIRS}")
    SET(CMAKE_REQUIRED_LIBRARIES "-lboost_iostreams")
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <string>
      #include <boost/iostreams/device/back_inserter.hpp>
      #include <boost/iostreams/filter/gzip.hpp>
      #include <boost/iostreams/filtering_stream.hpp>

      int main()
      {
        std::string decompressed_buffer;
        char test[1] = {'c'};

        boost::iostreams::filtering_ostream decompressing_stream;
        decompressing_stream.push(boost::iostreams::gzip_decompressor());
        decompressing_stream.push(boost::iostreams::back_inserter(decompressed_buffer));
        decompressing_stream.write (test, 1);
      }
      "
      DEAL_II_BOOST_IOSTREAMS_USUABLE
      )
    MESSAGE(STATUS "${DEAL_II_BOOST_IOSTREAMS_USUABLE}")
    IF(${DEAL_II_BOOST_IOSTREAMS_USUABLE})
      MESSAGE(STATUS "Boost.Iostreams is usuable.")
    ELSE()
      MESSAGE(STATUS
              "DEAL_II_WITH_ZLIB=ON requires Boost.Iostreams to be compiled "
              "with zlib support but a simple test failed! "
              "Therefore, the bundled boost package is used."
             )
      SET(Boost_FOUND FALSE)
      SET(Boost_LIBRARIES "")
      SET(Boost_INCLUDE_DIRS "")
    ENDIF()
    RESET_CMAKE_REQUIRED()
  ENDIF()

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
  )
