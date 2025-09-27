## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
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
# Search for a system library. In contrast to normal libraries we do this
# purely via "-l<library name>" instead of selecting a full library path..
#
# USAGE:
#   find_system_library(variable NAMES [list of possible names])
#

macro(find_system_library)
  set(_argn ${ARGN})
  list(GET _argn 0 _variable)
  list(REMOVE_AT _argn 0 1)

  if(NOT DEFINED ${_variable})
    foreach(_arg ${_argn})
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${_arg}")
      CHECK_CXX_SOURCE_COMPILES("int main(){}" ${_variable})
      reset_cmake_required()

      if(${_variable})
        unset(${_variable} CACHE)
        set(${_variable} ${_arg} CACHE STRING "A system library.")
        set(${_variable} ${_arg})
        mark_as_advanced(${_variable})
        break()
      else()
        unset(${_variable} CACHE)
      endif()
    endforeach()

    if(NOT ${_variable})
      set(${_variable} "${_variable}-NOTFOUND")
    endif()
  endif()
endmacro()
