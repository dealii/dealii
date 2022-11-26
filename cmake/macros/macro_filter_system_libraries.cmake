## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal.II authors
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
