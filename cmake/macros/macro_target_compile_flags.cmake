## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2023 - 2025 by the deal.II authors
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
# A small macro that takes a target, keyword and a string as argument and
# applies the given string as compile options to the target:
#
# Usage:
#   target_compile_flags(target INTERFACE|PUBLIC|PRIVATE "compile flags")
#
#   target_compile_flags(target INTERFACE|PUBLIC|PRIVATE
#      "$<generator expression>" "compile flags"
#      )
#
# In case of the second variant the compile options will be sandwiched
# between "$<$<generator expression>:" and ">".
#
function(target_compile_flags _target _keyword _string)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "\"${_target}\" is not a valid target")
  endif()
  if(NOT ${_keyword} MATCHES "(INTERFACE|PUBLIC|PRIVATE)")
    message(FATAL_ERROR
      "The supplied keyword has to be one of INTERFACE, PUBLIC, or PRIVATE"
      )
  endif()

  set(_guard_left "")
  set(_guard_right "")
  if("${_string}" MATCHES "^\\$<.*>$")
    set(_guard_left "$<${_string}:")
    set(_guard_right ">")
    set(_string "${ARGN}")
  endif()

  separate_arguments(_compile_options UNIX_COMMAND "${_string}")
  shell_escape_option_groups(_compile_options)
  target_compile_options(${_target} ${_keyword}
    ${_guard_left}${_compile_options}${_guard_right}
    )
endfunction()
