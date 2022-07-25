## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2021 by the deal.II authors
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
# Configuration for the Ginkgo library:
#

MACRO(FEATURE_GINKGO_FIND_EXTERNAL var)
  FIND_PACKAGE(GINKGO)

  IF(GINKGO_FOUND)
    SET(${var} TRUE)

    #
    # We require at least version 1.4.0
    # - The interface requires in fact only GINKGO 1.3.0, however, below 1.4.0
    #   the LD_LIBRARY_PATH has to be set manually by:
    #   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GINKGO_DIR/lib
    #
    SET(_version_required 1.4.0)
    IF(GINKGO_VERSION VERSION_LESS ${_version_required})
      MESSAGE(STATUS "Insufficient ginkgo installation found: "
        "At least version ${_version_required} is required."
        )
      SET(GINKGO_ADDITIONAL_ERROR_STRING
        "Insufficient ginkgo installation found!\n"
        "At least version ${_version_required} is required.\n"
        )
      SET(${var} FALSE)

    ENDIF()
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(GINKGO)
