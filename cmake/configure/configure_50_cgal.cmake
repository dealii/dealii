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
# Configuration for the CGAL library:
#

if(NOT FEATURE_BOOST_PROCESSED)
  message(FATAL_ERROR "\n"
    "Internal build system error: The configuration of "
    "DEAL_II_WITH_CGAL depends on "
    "DEAL_II_WITH_BOOST, but configure_feature(CGAL) "
    "was called before configure_feature(BOOST).\n\n"
    )
endif()


macro(feature_cgal_find_external var)
  find_package(DEAL_II_CGAL)

  if(CGAL_FOUND)
    set(${var} TRUE)
  endif()

  #
  # CGAL requires an full, externally installed Boost library. We can thus
  # not configure our internal boost and try to use CGAL at the same time.
  #
  if(DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED)
    message(STATUS
      "Could not find a sufficient CGAL installation: "
      "CGAL links against external Boost but deal.II was configured "
      "with bundled Boost."
      )
    set(CGAL_ADDITIONAL_ERROR_STRING
      ${CGAL_ADDITIONAL_ERROR_STRING}
      "Could not find a sufficient CGAL installation:\n"
      "CGAL links against external Boost but deal.II was configured "
      "with bundled Boost.\n\n"
      )
    set(${var} FALSE)
  endif()
endmacro()


macro(feature_cgal_configure_external)
  # Similarly to the DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS check run
  # in configure_20_boost.cmake we have to check whether cgal includes
  # a deprecated boost header. If yes, disable the boost deprecated header
  # warning as well.

  list(APPEND CMAKE_REQUIRED_INCLUDES
    ${BOOST_INCLUDE_DIRS} ${BOOST_BUNDLED_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS}
    )

  check_cxx_compiler_bug(
    "
    #define BOOST_CONFIG_HEADER_DEPRECATED_HPP_INCLUDED
    #define BOOST_HEADER_DEPRECATED(a) _Pragma(\"GCC error \\\"stop compilation\\\"\");
    #include <CGAL/make_mesh_3.h>
    int main() { return 0; }
    "
    DEAL_II_CGAL_HAS_DEPRECATED_BOOST_INCLUDES)

  reset_cmake_required()
endmacro()

configure_feature(CGAL)
