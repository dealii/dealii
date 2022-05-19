## ---------------------------------------------------------------------
##
## Copyright (C) 2022 by the deal.II authors
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
# Try to find the CGAL libraries
#
# This module exports
#
#   CGAL_INCLUDE_DIRS
#

SET(CGAL_DIR "" CACHE PATH "An optional hint to a CGAL installation")
SET_IF_EMPTY(CGAL_DIR "$ENV{CGAL_DIR}")

IF(NOT "${CGAL_DIR}" STREQUAL "")
  SET(CGAL_DIR ${CGAL_DIR})
ENDIF()

#
# CGAL requires C++17 and an externally configured Boost, otherwise the
# call to FIND_PACKAGE(CGAL) will fail. Guard the call to FIND_PACKAGE to
# fail cracefully:
#
IF(DEAL_II_HAVE_CXX17 AND NOT FEATURE_BOOST_BUNDLED_CONFIGURED)
  # temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
  LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  SET(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
  FIND_PACKAGE(CGAL QUIET)
  # Check version manually. Older binary distros don't do this properly.
  IF(CGAL_MAJOR_VERSION LESS 5)
    SET(CGAL_FOUND FALSE)
    SET(CGAL_INCLUDE_DIRS "-NOTFOUND")
    SET(CGAL_LIBRARIES "-NOTFOUND")
    MESSAGE(STATUS "CGAL wrappers require CGAL version 5 and above.")
  ENDIF()
  LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  IF(CGAL_FOUND)
    GET_TARGET_PROPERTY(CGAL_LIBRARIES CGAL::CGAL INTERFACE_LINK_LIBRARIES)
    # Make sure we dont' pass Boost::Boost over to deal.II.
    LIST(FILTER CGAL_LIBRARIES EXCLUDE REGEX "::")
  ENDIF()
ELSE()
  SET(CGAL_FOUND FALSE)
  SET(CGAL_INCLUDE_DIRS "-NOTFOUND")
  SET(CGAL_LIBRARIES "-NOTFOUND")
  MESSAGE(STATUS "CGAL wrappers require C++17. Disabling CGAL Support.")
ENDIF()

DEAL_II_PACKAGE_HANDLE(CGAL
  INCLUDE_DIRS
    REQUIRED CGAL_INCLUDE_DIRS
  LIBRARIES
    OPTIONAL CGAL_LIBRARIES
  USER_INCLUDE_DIRS
    REQUIRED CGAL_INCLUDE_DIRS
  CLEAR 
    CGAL_INCLUDE_DIRS
    CGAL_LIBRARIES
  )
