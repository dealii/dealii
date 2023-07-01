## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2022 by the deal.II authors
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

macro(feature_ginkgo_find_external var)
  find_package(DEAL_II_GINKGO)

  if(GINKGO_FOUND)
    set(${var} TRUE)

    #
    # We require at least version 1.4.0
    # - The interface requires in fact only GINKGO 1.3.0, however, below 1.4.0
    #   the LD_LIBRARY_PATH has to be set manually by:
    #   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GINKGO_DIR/lib
    #
    set(_version_required 1.4.0)
    if(GINKGO_VERSION VERSION_LESS ${_version_required})
      message(STATUS "Insufficient ginkgo installation found: "
        "At least version ${_version_required} is required."
        )
      set(GINKGO_ADDITIONAL_ERROR_STRING
        "Insufficient ginkgo installation found!\n"
        "At least version ${_version_required} is required.\n"
        )
      set(${var} FALSE)

    endif()
  endif()
endmacro()

configure_feature(GINKGO)
