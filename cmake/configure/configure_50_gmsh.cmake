## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Configuration for the gmsh executable:
#

set(FEATURE_GMSH_AFTER BOOST)

macro(feature_gmsh_find_external var)
  feature_find_external(GMSH ${var})

  if(${var})
    reset_cmake_required()
    if(DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED)
      list(APPEND CMAKE_REQUIRED_INCLUDES ${BOOST_FOLDER}/include)
    else()
      list(APPEND CMAKE_REQUIRED_INCLUDES
        ${BOOST_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
        )
      list(APPEND CMAKE_REQUIRED_LIBRARIES
        ${BOOST_TARGETS} ${BOOST_LIBRARIES} ${Boost_LIBRARIES}
        )
    endif()

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <boost/process/io.hpp>
      #include <boost/process/system.hpp>

      int main()
      {
        namespace bp = boost::process;

        return bp::system(\"gmsh\",
                          \"--version\",
                          bp::std_out > \"gmsh.log\",
                          bp::std_err > \"gmsh_warn.log\");
      }
      "
      DEAL_II_GMSH_WITH_BOOST_PROCESS
      )

    if(NOT DEAL_II_GMSH_WITH_BOOST_PROCESS)
      message(STATUS
        "DEAL_II_WITH_GMSH requires a usable Boost.Process setup, "
        "but a simple compile test failed."
        )
      set(GMSH_ADDITIONAL_ERROR_STRING
        "DEAL_II_WITH_GMSH requires a usable Boost.Process setup, "
        "but a simple compile test failed.\n"
        )
      set(${var} FALSE)
    endif()
    reset_cmake_required()
  endif()
endmacro()

macro(feature_gmsh_configure_external)
  set(DEAL_II_GMSH_WITH_API ${GMSH_WITH_API})
endmacro()

configure_feature(GMSH)
