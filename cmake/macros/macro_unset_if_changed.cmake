## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2018 - 2022 by the deal.II authors
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
# Usage:
#   unset_if_changed(<internal variable> "string"
#     [cached variable names]
#     )
#
# This macro caches the supplied "string" internally in ${<internal
# variable>} and unsets all supplied (cached) variables if this string
# changes.
#
macro(unset_if_changed _variable _string)
  if(DEFINED ${_variable})
    if(NOT "${${_variable}}" STREQUAL "${_string}")
      foreach(_arg ${ARGN})
        message(STATUS
          "Configuration changed. Unsetting cached variable \"${_arg}\" and rerunning checks.")
        unset(${_arg} CACHE)
      endforeach()
    endif()
  endif()
  set(${_variable} "${_string}" CACHE INTERNAL "" FORCE)
endmacro()
