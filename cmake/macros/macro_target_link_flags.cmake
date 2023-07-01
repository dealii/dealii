## ---------------------------------------------------------------------
##
## Copyright (C) 2023 - 2023 by the deal.II authors
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
# A small macro that takes a target, keyword and a string as argument and
# applies the given string as link options to the target:
#
# Usage:
#   target_link_flags(target INTERFACE|PUBLIC|PRIVATE "link flags")
#
#   target_link_flags(target INTERFACE|PUBLIC|PRIVATE
#      "$<generator expression>" "link flags"
#      )
#
# In case of the second variant the compile options will be sandwiched
# between "$<$<generator expression>:" and ">".
#
function(target_link_flags _target _keyword _string)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "»${_target}« is not a valid target")
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

  separate_arguments(_link_options UNIX_COMMAND "${_string}")
  shell_escape_option_groups(_link_options)
  target_link_options(${_target} ${_keyword}
    ${_guard_left}${_link_options}${_guard_right}
    )
endfunction()
