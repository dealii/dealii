## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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
# Configuration for the boost library:
#

SET(DEAL_II_WITH_BOOST ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for boost." FORCE
  )


MACRO(FEATURE_BOOST_CONFIGURE_COMMON)
  #
  # Boost version 1.62 - 1.63 checks for the availability of "emplace_hint"
  # incorrectly: It tests for the preprocessor define
  # BOOST_NO_CXX11_HDR_UNORDERED_MAP in .../boost/serialization/map.h
  # thinking that that this define is characteristic for the presence of
  # std::(multi)map::emplace_hint. This is generally correct, except for
  # GCC before 4.8, for which the preprocessor variable is defined, but the
  # function does not exist [1].
  #
  # Thus, simply define a BOOST_NO_CXX11_HDR_UNORDERED_MAP if the gcc
  # compiler version is less than 4.8.
  #
  # [1] https://svn.boost.org/trac/boost/ticket/12755
  #
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
      CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
    LIST(APPEND BOOST_DEFINITIONS "BOOST_NO_CXX11_HDR_UNORDERED_MAP")
    LIST(APPEND BOOST_USER_DEFINITIONS "BOOST_NO_CXX11_HDR_UNORDERED_MAP")
  ENDIF()

  ENABLE_IF_SUPPORTED(BOOST_CXX_FLAGS "-Wno-unused-local-typedefs")
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_BUNDLED)
  #
  # Add rt to the link interface as well, boost/chrono needs it.
  #
  IF(NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    FIND_SYSTEM_LIBRARY(rt_LIBRARY NAMES rt)
    MARK_AS_ADVANCED(rt_LIBRARY)
    IF(NOT rt_LIBRARY MATCHES "-NOTFOUND")
      SET(BOOST_LIBRARIES ${rt_LIBRARY})
    ENDIF()
  ENDIF()

  FEATURE_BOOST_CONFIGURE_COMMON()

  SET(BOOST_BUNDLED_INCLUDE_DIRS ${BOOST_FOLDER}/include)

  IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
    #
    # Bundled boost tries to (dl)open itself as a dynamic library on
    # Windows. Disable this undesired behavior by exporting
    # BOOST_ALL_NO_LIB on Windows platforms (for bundled boost).
    #
    LIST(APPEND BOOST_DEFINITIONS "BOOST_ALL_NO_LIB")
    LIST(APPEND BOOST_USER_DEFINITIONS "BOOST_ALL_NO_LIB")
  ENDIF()
ENDMACRO()

MACRO(FEATURE_BOOST_FIND_EXTERNAL var)
  FIND_PACKAGE(BOOST)

  IF(BOOST_FOUND)
    SET(${var} TRUE)

    #
    # Blacklist version 1.58 because we get serialization errors with it. At
    # least version 1.56 and 1.59 are known to work.
    #
    IF("${BOOST_VERSION_MAJOR}" STREQUAL "1" AND "${BOOST_VERSION_MINOR}" STREQUAL "58")
      MESSAGE(STATUS "Boost version 1.58 is not compatible with deal.II!")
      SET(${var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_EXTERNAL)
  FEATURE_BOOST_CONFIGURE_COMMON()
ENDMACRO()


CONFIGURE_FEATURE(BOOST)


#
# DEAL_II_WITH_BOOST is always required.
#
IF(NOT DEAL_II_WITH_BOOST)
  IF(DEAL_II_FEATURE_AUTODETECTION)
    FEATURE_ERROR_MESSAGE("BOOST")
  ELSE()
    MESSAGE(FATAL_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_BOOST required, but set to OFF!.\n\n"
      )
  ENDIF()
ENDIF()
