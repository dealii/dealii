## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2019 by the deal.II authors
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
    # We don't support version 4.0.0 or later yet.
    #
    SET(_first_unsupported_sundials_version 4.0.0)
    IF(NOT SUNDIALS_VERSION VERSION_LESS ${_first_unsupported_sundials_version})
      MESSAGE(STATUS
              "Insufficient SUNDIALS installation found: "
              "version ${_first_unsupported_sundials_version} "
              "or later is not yet supported, "
              "but version ${SUNDIALS_VERSION} was found."
        )
      SET(SUNDIALS_ADDITIONAL_ERROR_STRING
          "Insufficient SUNDIALS installation found!\n"
          "Version ${_first_unsupported_sundials_version} "
          "or later is not yet supported, "
          "but version ${SUNDIALS_VERSION} was found.\n"
        )
      SET(${var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO(FEATURE_SUNDIALS_CONFIGURE_EXTERNAL)
  SET(DEAL_II_SUNDIALS_WITH_IDAS ${SUNDIALS_WITH_IDAS})
ENDMACRO()

CONFIGURE_FEATURE(SUNDIALS)
