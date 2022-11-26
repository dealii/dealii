## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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
# A macro that evaluates an expression as supplied by a string. Suggestion
# from the Wiki http://cmake.org/Wiki/CMake/Language_Syntax
#
# USAGE:
#
# evaluate_expression("<expression>")
#

macro(evaluate_expression _the_expression)
  set(_tmp_name
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/evaluate_expression.tmp"
    )
  file(WRITE ${_tmp_name} "${_the_expression}")
  include("${_tmp_name}")
endmacro()
