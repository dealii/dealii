## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
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
# Configuration for the SUNDIALS library:
#

macro(feature_sundials_find_external var)
  find_package(DEAL_II_SUNDIALS)

  if(SUNDIALS_FOUND)
    set(${var} TRUE)

    #
    # We require at least sundials 5.4.0
    #
    set(_version_required 5.4.0)
    if(SUNDIALS_VERSION VERSION_LESS ${_version_required})
      message(STATUS "Could not find a sufficient Sundials installation: "
        "deal.II requires at least version ${_version_required}, "
        "but version ${SUNDIALS_VERSION} was found."
        )
      set(SUNDIALS_ADDITIONAL_ERROR_STRING
        ${SUNDIALS_ADDITIONAL_ERROR_STRING}
        "The SUNDIALS installation (found at \"${SUNDIALS_DIR}\")\n"
        "with version ${SUNDIALS_VERSION} is too old.\n"
        "deal.II requires at least version ${_version_required}.\n\n"
        )
      set(${var} FALSE)
    endif()
  endif()
endmacro()

macro(feature_sundials_configure_external)
  set(DEAL_II_SUNDIALS_WITH_IDAS ${SUNDIALS_WITH_IDAS})
endmacro()

configure_feature(SUNDIALS)
