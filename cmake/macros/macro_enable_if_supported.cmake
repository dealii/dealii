## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Tests whether the cxx compiler understands a flag.
# If so, add it to 'variable'.
#
# Usage:
#     ENABLE_IF_SUPPORTED(variable flag)
#

MACRO(ENABLE_IF_SUPPORTED _variable _flag)
  #
  # Clang is too conservative when reporting unsupported compiler flags.
  # Therefore, we promote all warnings for an unsupported compiler flag to
  # actual errors with the -Werror switch:
  #
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(_werror_string "-Werror ")
  ELSE()
    SET(_werror_string "")
  ENDIF()

  STRING(STRIP "${_flag}" _flag_stripped)
  SET(_flag_stripped_orig "${_flag_stripped}")

  #
  # Gcc does not emit a warning if testing -Wno-... flags which leads to
  # false positive detection. Unfortunately it later warns that an unknown
  # warning option is used if another warning is emited in the same
  # compilation unit.
  # Therefore we invert the test for -Wno-... flags:
  #
  IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    STRING(REPLACE "-Wno-" "-W" _flag_stripped "${_flag_stripped}")
  ENDIF()

  IF(NOT "${_flag_stripped}" STREQUAL "")
    STRING(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    STRING(REPLACE "," "" _flag_name "${_flag_name}")
    STRING(REPLACE "-" "_" _flag_name "${_flag_name}")
    STRING(REPLACE "+" "_" _flag_name "${_flag_name}")
    CHECK_CXX_COMPILER_FLAG(
      "${_werror_string}${_flag_stripped}"
      DEAL_II_HAVE_FLAG_${_flag_name}
      )
    IF(DEAL_II_HAVE_FLAG_${_flag_name})
      SET(${_variable} "${${_variable}} ${_flag_stripped_orig}")
      STRING(STRIP "${${_variable}}" ${_variable})
    ENDIF()
  ENDIF()
ENDMACRO()

