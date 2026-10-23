## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Remove the linker selection flags "-fuse-ld=mold" and "-fuse-ld=lld" that
# were added automatically during configuration from the string variable.
# A flag is kept if the user supplied it, i.e., if it is also part of
# "${variable}_SAVED". The list of removed flags is stored in
# "${variable}_DETECTED".
#
# Usage:
#     strip_detected_linker_flags(variable)
#

macro(strip_detected_linker_flags _variable)
  set(${_variable}_DETECTED "")
  foreach(_flag "-fuse-ld=mold" "-fuse-ld=lld")
    if(" ${${_variable}} " MATCHES " ${_flag} " AND
       NOT " ${${_variable}_SAVED} " MATCHES " ${_flag} ")
      strip_flag(${_variable} "${_flag}")
      list(APPEND ${_variable}_DETECTED "${_flag}")
    endif()
  endforeach()
endmacro()
