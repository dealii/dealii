## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2022 - 2023 by the deal.II authors
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
# Try to find the CGAL libraries
#
# This module exports:
#
# CGAL_INCLUDE_DIRS
# CGAL_VERSION_MAJOR
# CGAL_VERSION_MINOR
# CGAL_VERSION_SUBMINOR
#

set(CGAL_DIR "" CACHE PATH "An optional hint to a CGAL installation")
set_if_empty(CGAL_DIR "$ENV{CGAL_DIR}")

if(NOT "${CGAL_DIR}" STREQUAL "")
  set(CGAL_DIR ${CGAL_DIR})
endif()

#
# CGAL requires C++17 and an externally configured Boost, otherwise the
# call to find_package(CGAL) will fail. Guard the call to FIND_PACKAGE to
# fail gracefully:
#
if(DEAL_II_HAVE_CXX17 AND NOT DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED)
  set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
  find_package(CGAL QUIET)

  # Check version manually. Older binary distros don't do this properly.
  if(CGAL_MAJOR_VERSION LESS 5)
    set(CGAL_FOUND FALSE)
    set(CGAL_INCLUDE_DIRS "-NOTFOUND")
    set(CGAL_LIBRARIES "-NOTFOUND")
    message(STATUS "CGAL wrappers require CGAL version 5 and above.")
  endif()

  if(CGAL_FOUND)
    get_target_property(CGAL_LIBRARIES CGAL::CGAL INTERFACE_LINK_LIBRARIES)

    if(DEFINED CGAL_VERSION)
      set(CGAL_VERSION "${CGAL_VERSION}")

      string(REGEX REPLACE
        "^([0-9]+).*$" "\\1"
        CGAL_VERSION_MAJOR "${CGAL_VERSION}")

      string(REGEX REPLACE
        "^[0-9]+\\.([0-9]+).*$" "\\1"
        CGAL_VERSION_MINOR "${CGAL_VERSION}")

      string(REGEX REPLACE
        "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
        CGAL_VERSION_SUBMINOR "${CGAL_VERSION}"
      )

      if("${CGAL_VERSION_SUBMINOR}" MATCHES "^(|${CGAL_VERSION})$")
        set(CGAL_VERSION_SUBMINOR "0")
      endif()
    endif()

    # Make sure we don't pass Boost::Boost over to deal.II.
    list(FILTER CGAL_LIBRARIES EXCLUDE REGEX "::")
  endif()
else()
  set(CGAL_FOUND FALSE)
  set(CGAL_INCLUDE_DIRS "-NOTFOUND")
  set(CGAL_LIBRARIES "-NOTFOUND")
  message(STATUS "CGAL wrappers require C++17. Disabling CGAL Support.")
endif()

process_feature(CGAL
  LIBRARIES OPTIONAL CGAL_LIBRARIES
  INCLUDE_DIRS REQUIRED CGAL_INCLUDE_DIRS
  CLEAR
    CGAL_INCLUDE_DIRS
    CGAL_LIBRARIES
)
