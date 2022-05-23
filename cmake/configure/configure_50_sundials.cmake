## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2021 by the deal.II authors
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
# Configuration for the SUNDIALS library:
#

MACRO(FEATURE_SUNDIALS_FIND_EXTERNAL var)
  FIND_PACKAGE(SUNDIALS)

  IF(SUNDIALS_FOUND)
    SET(${var} TRUE)

    #
    # We require at least sundials 5.4.0
    #
    SET(_version_required 5.4.0)
    IF(SUNDIALS_VERSION VERSION_LESS ${_version_required})
      MESSAGE(STATUS "Could not find a sufficient Sundials installation: "
        "deal.II requires at least version ${_version_required}, "
        "but version ${SUNDIALS_VERSION} was found."
        )
      SET(SUNDIALS_ADDITIONAL_ERROR_STRING
        ${SUNDIALS_ADDITIONAL_ERROR_STRING}
        "The SUNDIALS installation (found at \"${SUNDIALS_DIR}\")\n"
        "with version ${SUNDIALS_VERSION} is too old.\n"
        "deal.II requires at least version ${_version_required}.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO(FEATURE_SUNDIALS_CONFIGURE_EXTERNAL)
  SET(DEAL_II_SUNDIALS_WITH_IDAS ${SUNDIALS_WITH_IDAS})
ENDMACRO()

CONFIGURE_FEATURE(SUNDIALS)
