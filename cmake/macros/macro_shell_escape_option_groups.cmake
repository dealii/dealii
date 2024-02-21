## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2023 by the deal.II authors
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
  # do this by detecting all occurrences of a flag ("-Wl,-[-]flag") followed
  # by an option that doesn't start with a dash ("-Wl,[option]"):
  #
  string(REGEX REPLACE
    "(-Wl,[-]+[^;]*);(-Wl,[^-][^;]+)"
    "SHELL:\\1 \\2" ${_variable} "${${_variable}}"
    )
endmacro()
