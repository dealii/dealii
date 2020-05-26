## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
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
# Configuration for thread support in deal.II with the help of the tbb
# library:
#

MACRO(FEATURE_TBB_FIND_EXTERNAL var)
  FIND_PACKAGE(TBB)

  IF(TBB_FOUND)
    SET(${var} TRUE)
  ENDIF()

  #
  # TBB currently uses the version numbering scheme
  #
  #     YYYY.X
  #
  # (e.g., 2018.0) where YYYY is the year of the release and X is the yearly
  # release number. Older versions use
  #
  #     X.Y.Z
  #
  # (e.g., 4.2.1). Since we are compatible with all versions that use the new
  # numbering scheme we only check for very old versions here.
  #
  # TBB versions before 4.2 are missing some explicit calls to std::atomic::load
  # in ternary expressions; these cause compilation errors in some compilers
  # (such as GCC 8.1 and newer). To fix this we simply blacklist all older
  # versions:
  #
  IF(TBB_VERSION VERSION_LESS "4.2")
    # Clear the previously determined version numbers to avoid confusion
    SET(TBB_VERSION "bundled")
    SET(TBB_VERSION_MAJOR "")
    SET(TBB_VERSION_MINOR "")

    MESSAGE(STATUS
      "The externally provided TBB library is older than version 4.2.0, which "
      "cannot be used with deal.II."
      )
    SET(TBB_ADDITIONAL_ERROR_STRING
      "The externally provided TBB library is older than version\n"
      "4.2.0, which is the oldest version compatible with deal.II and its\n"
      "supported compilers."
      )
    SET(${var} FALSE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_TBB_CONFIGURE_BUNDLED)
  #
  # We have to disable a bunch of warnings:
  #
  ENABLE_IF_SUPPORTED(TBB_CXX_FLAGS "-Wno-parentheses")

  #
  # tbb uses dlopen/dlclose, so link against libdl.so as well:
  #
  LIST(APPEND TBB_LIBRARIES ${CMAKE_DL_LIBS})

  LIST(APPEND TBB_BUNDLED_INCLUDE_DIRS ${TBB_FOLDER}/include)
ENDMACRO()


CONFIGURE_FEATURE(TBB)
