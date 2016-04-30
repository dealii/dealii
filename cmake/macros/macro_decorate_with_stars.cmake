## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

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
#     DECORATE_WITH_STARS(message decorated_message)
#
#
MACRO(DECORATE_WITH_STARS _message _decorated_message)
  STRING(LENGTH ${_message} _message_length)
  SET(_line_length 75)
  MATH(EXPR _unpadded_line_length "${_line_length} - 6")

  IF(${_message_length} LESS ${_unpadded_line_length})
    MATH(EXPR _left_pad_size "(${_unpadded_line_length} - ${_message_length} + 1)/2")
    MATH(EXPR _right_pad_size "(${_unpadded_line_length} - ${_message_length})/2")
    # Unfortunately, it looks like taking substrings is the only way to pad
    # strings with run time dependent length.
    SET(_pad_strings_values "                                              ")
    STRING(SUBSTRING "${_pad_strings_values}" 0 "${_left_pad_size}" _left_pad)
    STRING(SUBSTRING "${_pad_strings_values}" 0 "${_right_pad_size}" _right_pad)
    SET(${_decorated_message} "** ${_left_pad}${_message}${_right_pad} **")
  ELSE()
    SET(${_decorated_message} "** ${_message}")
  ENDIF()
ENDMACRO()
