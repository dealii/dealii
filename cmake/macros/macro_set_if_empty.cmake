## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# If 'variable' is empty it will be set to 'value'
#
macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()
