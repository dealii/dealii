## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# A small wrapper around
# set_target_property(... PROPERTIES COMPILE_DEFINITIONS ...)
# to _add_ compile definitions to every target we have specified.
#

macro(deal_ii_add_definitions _name)

  foreach(_build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)

    set_property(TARGET ${_name}_${_build_lowercase}
      APPEND PROPERTY COMPILE_DEFINITIONS "${ARGN}"
      )
  endforeach()

endmacro()
