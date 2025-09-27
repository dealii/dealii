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
# A small wrapper around FIND_FILE to be a bit more verbose
#

macro(deal_ii_find_file _file_name)
  # Save a string representation of the arguments before cmake's
  # FIND_FILE gets its hands on it.
  to_string(_str ${ARGN})

  find_file(${_file_name} ${ARGN})

  if(${_file_name} MATCHES "-NOTFOUND")
    message(STATUS "${_file_name} not found! The call was:")
    message(STATUS "    find_file(${_file_name} ${_str})")
  else()
    message(STATUS "Found ${_file_name}")
  endif()
endmacro()
