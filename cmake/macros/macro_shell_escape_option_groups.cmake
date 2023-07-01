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
# Shell escape known option groups such as -Xlinker [...], or
# -Xcudafe [...] so that the toggle parameter (-Xlinker) does not get
# deduplicated by CMake.
#
# Prerequisite: variable must be the name of a _list_ variable containing
# compile or link options.
#
# Usage:
#     shell_escape_option_groups(variable)
#

macro(shell_escape_option_groups _variable)
  #
  # We have to capture simple option groups of the form
  #   "-Xlinker --as-needed"
  # but also multiple escape sequences such as
  #   "-Xlinker -rpath -Xlinker /path".
  # If we happen to escape the latter in the following form:
  # "SHELL:-Xlinker -rpath" "SHELL:-Xlinker /path" then the first statement
  # "SHELL:-Xlinker -rpath" might still get deduplicated...
  #
  # Let's play this game only over the two options groups above in order to
  # keep the regex sane.
  #

  # Matches: "-Xlinker;[option]" and replaces it with "SHELL:-Xlinker [option]"
  string(REGEX REPLACE
    "(-Xcudafe|-Xlinker);([^;]+)"
    "SHELL:\\1 \\2" ${_variable} "${${_variable}}"
    )
  # Matches: "SHELL:-Xlinker [option1];SHELL:-Xlinker [option2]" and replaces
  # it with "SHELL:-Xlinker [option1] -Xlinker [option2]"
  string(REGEX REPLACE
    "SHELL:(-Xcudafe [^;]+|-Xlinker [^;]+);SHELL:(-Xcudafe [^;]+|-Xlinker [^;]+)"
    "SHELL:\\1 \\2" ${_variable} "${${_variable}}"
    )

  #
  # In addition try to merge options of the form "-Wl,-flag -Wl,/path". We
  # do this by detecting all occurences of a flag ("-Wl,-[-]flag") followed
  # by an option that doesn't start with a dash ("-Wl,[option]"):
  #
  string(REGEX REPLACE
    "(-Wl,[-]+[^;]*);(-Wl,[^-][^;]+)"
    "SHELL:\\1 \\2" ${_variable} "${${_variable}}"
    )
endmacro()
