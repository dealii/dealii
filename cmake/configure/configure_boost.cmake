## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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

  ENABLE_IF_SUPPORTED(BOOST_CXX_FLAGS "-Wno-unused-local-typedefs")

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
  ENABLE_IF_SUPPORTED(BOOST_CXX_FLAGS "-Wno-unused-local-typedefs")
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
