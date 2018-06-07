## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2018 by the deal.II authors
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
# Configuration for the ADOL-C library:
#

SET(FEATURE_ADOLC_AFTER BOOST)

MACRO(FEATURE_ADOLC_FIND_EXTERNAL var)
  FIND_PACKAGE(ADOLC)

  IF(ADOLC_FOUND)
    #
    # So, we have a library. Let's see whether we can use it:
    #
    SET(${var} TRUE)

    #
    # If Adolc is configured to use the Boost allocator (of an external
    # boost library) we must not use a bundled Boost library for deal.II.
    #
    IF(ADOLC_WITH_BOOST_ALLOCATOR AND FEATURE_BOOST_BUNDLED_CONFIGURED)
      MESSAGE(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "ADOL-C links against external Boost but deal.II was configured "
        "with bundled Boost."
        )
      SET(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "ADOL-C links against external Boost but deal.II was configured "
        "with bundled Boost.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Check whether we have a recent enough ADOL-C library that can return
    # values from constant objects.
    #

    SET(CMAKE_REQUIRED_LIBRARIES ${ADOLC_LIBRARIES})
    SET(CMAKE_REQUIRED_INCLUDES ${ADOLC_INCLUDE_DIRS})
    SET(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS}")
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

    IF(NOT ADOLC_DOUBLE_CAST_CHECK)
      MESSAGE(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "deal.II needs ADOL-C version 2.6.4 or newer."
        )
      SET(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "ADOL-C cast check failed.\n"
        "deal.II needs ADOL-C version 2.6.4 or newer.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()

    IF(NOT ADOLC_ADOUBLE_OSTREAM_CHECK)
      MESSAGE(STATUS
        "Could not find a sufficient ADOL-C installation: "
        "deal.II needs ADOL-C version 2.6.4 or newer."
        )
      SET(ADOLC_ADDITIONAL_ERROR_STRING
        ${ADOLC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ADOL-C installation:\n"
        "ADOL-C stream output check failed.\n"
        "deal.II needs ADOL-C version 2.6.4 or newer.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()


MACRO(FEATURE_ADOLC_CONFIGURE_EXTERNAL)
  SET(DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING ${ADOLC_WITH_ADVANCED_BRANCHING})
  SET(DEAL_II_ADOLC_WITH_ATRIG_ERF ${ADOLC_WITH_ATRIG_ERF})
  SET(DEAL_II_ADOLC_WITH_BOOST_ALLOCATOR ${ADOLC_WITH_BOOST_ALLOCATOR})

  SET(DEAL_II_EXPAND_ADOLC_TYPES "adouble; adtl::adouble")
ENDMACRO()


CONFIGURE_FEATURE(ADOLC)
