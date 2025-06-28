## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
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
# A small macro used for (string-)appending a string "${flags}" to a
# string "${variable}"
#
# Usage:
#     add_flags(variable flags)
#

macro(add_flags _variable _flags)
  string(STRIP "${_flags}" _flags_stripped)
  if(NOT "${_flags_stripped}" STREQUAL "")
    set(${_variable} "${${_variable}} ${_flags}")
    string(STRIP "${${_variable}}" ${_variable})
  endif()
endmacro()
