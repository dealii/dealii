## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2023 by the deal.II authors
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
# Configuration for the ADOL-C library:
#

set(FEATURE_ADOLC_AFTER BOOST TRILINOS)

macro(feature_adolc_find_external var)
  find_package(DEAL_II_ADOLC)

  if(ADOLC_FOUND)
    #
    # So, we have a library. Let's see whether we can use it:
    #
    set(${var} TRUE)

    #
    # If Adolc is configured to use the Boost allocator (of an external
    # boost library) we must not use a bundled Boost library for deal.II.
    #
    if(ADOLC_WITH_BOOST_ALLOCATOR AND DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED)
      message(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "ADOL-C links against external Boost but deal.II was configured "
        "with bundled Boost."
        )
      set(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "ADOL-C links against external Boost but deal.II was configured "
        "with bundled Boost.\n\n"
        )
      set(${var} FALSE)
    endif()

    #
    # We have to avoid a symbol clash with Trilinos' SEACASChaco library
    # (the libchaco.so shared object exports the global symbol 'divide' but
    # so does adolc itself).
    #
    item_matches(_module_found SEACASChaco ${Trilinos_PACKAGE_LIST})
    if(_module_found)
      message(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "Possible symbol clash between the ADOL-C library and Trilinos' SEACASChaco detected"
        )
      set(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "Possible symbol clash between the ADOL-C library and Trilinos' SEACASChaco detected."
        "If you want to use ADOL-C, please configure deal.II to use a "
        "Trilinos library with disabled SEACASChaco.\n\n"
        )
      set(${var} FALSE)
    endif()

    #
    # We have to avoid another symbol clash with the netcdf library that
    # might get transitively pulled in by Trilinos (the libnetcdf.so shared
    # object exports the global symbol 'function' but so does adolc
    # itself).
    #
    if("${Trilinos_TPL_LIBRARIES}" MATCHES "netcdf")
      message(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "Possible symbol clash between the ADOL-C library and netcdf "
        "(pulled in as optional external dependency of Trilinos) detected"
        )
      set(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "Possible symbol clash between the ADOL-C library and netcdf "
        "(pulled in as optional external dependency of Trilinos). "
        "If you want to use ADOL-C, please configure deal.II to use a "
        "Trilinos library with disabled netcdf bindings.\n\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Check whether we have a recent enough ADOL-C library that can return
    # values from constant objects.
    #

    list(APPEND CMAKE_REQUIRED_LIBRARIES ${ADOLC_LIBRARIES})
    list(APPEND CMAKE_REQUIRED_INCLUDES ${ADOLC_INCLUDE_DIRS})

    CHECK_CXX_SOURCE_COMPILES("
      #include <adolc/adouble.h>
      #include <iostream>
      void print_double(const adouble &val)
      {
        std::cout << \"val (non-const): \" << static_cast<double>(val) << std::endl;
      }
      int main (int argc, char *argv[])
      {
        const adouble val = 1.0;
        print_double(val);
      }"
      ADOLC_DOUBLE_CAST_CHECK)

    CHECK_CXX_SOURCE_COMPILES("
      #include <adolc/adouble.h>
      #include <adolc/adtl.h>
      #include <sstream>
      int main (int argc, char *argv[])
      {
        const adouble val_taped = 1.0;
        const adtl::adouble val_tapeless = 1.0;

        std::ostringstream ss;
        ss << val_taped;
        ss << val_tapeless;
      }"
      ADOLC_ADOUBLE_OSTREAM_CHECK)

    reset_cmake_required()

    if(NOT ADOLC_DOUBLE_CAST_CHECK)
      message(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "deal.II needs ADOL-C version 2.6.4 or newer."
        )
      set(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "ADOL-C cast check failed.\n"
        "deal.II needs ADOL-C version 2.6.4 or newer.\n\n"
        )
      set(${var} FALSE)
    endif()

    if(NOT ADOLC_ADOUBLE_OSTREAM_CHECK)
      message(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "deal.II needs ADOL-C version 2.6.4 or newer."
        )
      set(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "ADOL-C stream output check failed.\n"
        "deal.II needs ADOL-C version 2.6.4 or newer.\n\n"
        )
      set(${var} FALSE)
    endif()
  endif()
endmacro()


macro(feature_adolc_configure_external)
  set(DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING ${ADOLC_WITH_ADVANCED_BRANCHING})
  set(DEAL_II_ADOLC_WITH_ATRIG_ERF ${ADOLC_WITH_ATRIG_ERF})
  set(DEAL_II_ADOLC_WITH_TAPELESS_REFCOUNTING ${ADOLC_WITH_TAPELESS_REFCOUNTING})
  set(DEAL_II_ADOLC_WITH_BOOST_ALLOCATOR ${ADOLC_WITH_BOOST_ALLOCATOR})

  set(DEAL_II_EXPAND_ADOLC_TYPES "adouble; adtl::adouble")
endmacro()


configure_feature(ADOLC)
