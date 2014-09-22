## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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
# This macro replaces absolute paths to system libraries with the
# corresponding short name within the FEATURE_LIBRARIES(|_DEBUG|_RELEASE)
# variables
#
# Usage:
#     FILTER_SYSTEM_LIBRARIES(feature)
#

MACRO(FILTER_SYSTEM_LIBRARIES _feature)
  FOREACH(_variable
    ${_feature}_LIBRARIES
    ${_feature}_LIBRARIES_DEBUG
    ${_feature}_LIBRARIES_RELEASE
    )
    IF(DEFINED ${_variable})
      SET(_tmp_${_variable} ${${_variable}})
      SET(${_variable} "")
      FOREACH(_lib ${_tmp_${_variable}})
        IF(_lib MATCHES "lib(c|quadmath|gfortran|m|rt|nsl|dl|pthread)\\.(a|so)$")
          string(REGEX REPLACE ".*lib([a-z]+).so$" "\\1" _lib ${_lib})
        ENDIF()
        LIST(APPEND ${_variable} ${_lib})
      ENDFOREACH()
    ENDIF()
  ENDFOREACH()
ENDMACRO()
