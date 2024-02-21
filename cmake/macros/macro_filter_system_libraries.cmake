## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
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
# This macro replaces absolute paths to system libraries with the
# corresponding short name within the FEATURE_LIBRARIES(|_DEBUG|_RELEASE)
# variables
#
# Usage:
#     filter_system_libraries(feature)
#

macro(filter_system_libraries _feature)
  foreach(_variable
    ${_feature}_LIBRARIES
    ${_feature}_LIBRARIES_DEBUG
    ${_feature}_LIBRARIES_RELEASE
    )
    if(DEFINED ${_variable})
      set(_tmp_${_variable} ${${_variable}})
      set(${_variable} "")
      foreach(_lib ${_tmp_${_variable}})
        if(_lib MATCHES "lib(bfd|c|dl|gfortran|iberty|m|nsl|opcodes|pthread|quadmath|rt)\\.(a|so)$")
          string(REGEX REPLACE ".*lib([a-z]+).so$" "\\1" _lib ${_lib})
        endif()
        list(APPEND ${_variable} ${_lib})
      endforeach()
    endif()
  endforeach()
endmacro()
