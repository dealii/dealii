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
# Configuration for the CGAL library:
#

IF(NOT FEATURE_BOOST_PROCESSED)
  MESSAGE(FATAL_ERROR "\n"
    "Internal build system error: The configuration of "
    "DEAL_II_WITH_CGAL depends on "
    "DEAL_II_WITH_BOOST, but CONFIGURE_FEATURE(CGAL) "
    "was called before CONFIGURE_FEATURE(BOOST).\n\n"
    )
ENDIF()


MACRO(FEATURE_CGAL_FIND_EXTERNAL var)
  FIND_PACKAGE(CGAL)

  IF(CGAL_FOUND)
    SET(${var} TRUE)
  ENDIF()

  #
  # CGAL requires an full, externally installed Boost library. We can thus
  # not configure our internal boost and try to use CGAL at the same time.
  #
  IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
    MESSAGE(STATUS
      "Could not find a sufficient CGAL installation: "
      "CGAL links against external Boost but deal.II was configured "
      "with bundled Boost."
      )
    SET(CGAL_ADDITIONAL_ERROR_STRING
      ${CGAL_ADDITIONAL_ERROR_STRING}
      "Could not find a sufficient CGAL installation:\n"
      "CGAL links against external Boost but deal.II was configured "
      "with bundled Boost.\n\n"
      )
    SET(${var} FALSE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_CGAL_CONFIGURE_EXTERNAL)
  # Similarly to the DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS check run
  # in configure_20_boost.cmake we have to check whether cgal includes
  # a deprecated boost header. If yes, disable the boost deprecated header
  # warning as well.

  LIST(APPEND CMAKE_REQUIRED_INCLUDES
    ${BOOST_INCLUDE_DIRS} ${BOOST_BUNDLED_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS}
    )

  CHECK_CXX_COMPILER_BUG(
    "
    #define BOOST_CONFIG_HEADER_DEPRECATED_HPP_INCLUDED
    #define BOOST_HEADER_DEPRECATED(a) _Pragma(\"GCC error \\\"stop compilation\\\"\");
    #include <CGAL/make_mesh_3.h>
    int main() { return 0; }
    "
    DEAL_II_CGAL_HAS_DEPRECATED_BOOST_INCLUDES)

  RESET_CMAKE_REQUIRED()
ENDMACRO()

CONFIGURE_FEATURE(CGAL)
