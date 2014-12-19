## ---------------------------------------------------------------------
##
## Copyright (C) 2013 by the deal.II authors
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
# A small wrapper around FIND_PACKAGE.
# We guard the invocation of FIND_PACKAGE(package <...>) by
# ${package}_FOUND and ${package}_LIBRARIES to allow easy custom overrides
#

MACRO(FIND_PACKAGE _package_name)
  STRING(TOUPPER ${_package_name} _package_name_uppercase)

  IF( NOT DEFINED ${_package_name_uppercase}_FOUND AND
      NOT DEFINED ${_package_name_uppercase}_LIBRARIES )
    _FIND_PACKAGE (${_package_name} ${ARGN})
  ELSE()
    IF(NOT DEFINED ${_package_name_uppercase}_FOUND)
      SET(${_package_name_uppercase}_FOUND TRUE)
    ENDIF()
  ENDIF()
ENDMACRO()
