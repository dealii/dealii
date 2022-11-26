## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2020 by the deal.II authors
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
