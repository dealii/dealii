## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2016 - 2022 by the deal.II authors
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
# Formats a string to be of length 75 and surrounded with '**'s, and, if
# necessary, centered with extra spaces. This allows for nicer looking output
# messages: for example, the middle two lines of
#
# ***************************************************************************
# **    Error: Could not build and install disabled component examples.    **
# **        Please reconfigure with -DDEAL_II_COMPONENT_EXAMPLES=ON        **
# ***************************************************************************
#
# could be built automatically with this macro. If the line is too long to fit
# in such a box then decorate it with just two stars at the beginning.
#
# Usage:
#     decorate_with_stars(message decorated_message)
#
#
macro(decorate_with_stars _message _decorated_message)
  string(LENGTH ${_message} _message_length)
  set(_line_length 75)
  math(EXPR _unpadded_line_length "${_line_length} - 6")

  if(${_message_length} LESS ${_unpadded_line_length})
    math(EXPR _left_pad_size "(${_unpadded_line_length} - ${_message_length} + 1)/2")
    math(EXPR _right_pad_size "(${_unpadded_line_length} - ${_message_length})/2")
    # Unfortunately, it looks like taking substrings is the only way to pad
    # strings with run time dependent length.
    set(_pad_strings_values "                                              ")
    string(SUBSTRING "${_pad_strings_values}" 0 "${_left_pad_size}" _left_pad)
    string(SUBSTRING "${_pad_strings_values}" 0 "${_right_pad_size}" _right_pad)
    set(${_decorated_message} "** ${_left_pad}${_message}${_right_pad} **")
  else()
    set(${_decorated_message} "** ${_message}")
  endif()
endmacro()
