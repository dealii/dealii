## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
