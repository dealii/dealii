## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2023 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# strip_known_generator_expressions(<variable>)
#
# Strip an enclosing generator expression from the variable. This macro is
# primarily used in copy_target_properties
#

macro(strip_known_generator_expressions _variable)
  set(generator_expressions
    "\\$<LINK_ONLY:([^>]*)>"
    "\\$<\\$<LINK_LANGUAGE:CXX>:([^>]*)>"
    "\\$<\\$<COMPILE_LANGUAGE:CXX>:([^>]*)>"
    )

  foreach(expression ${generator_expressions})
    string(REGEX REPLACE ${expression} "\\1" ${_variable} "${${_variable}}")
  endforeach()
endmacro()
