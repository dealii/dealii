#####
##
## Copyright (C) 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# A small wrapper around FIND_PACKAGE.
# We guard the invocation of FIND_PACKAGE(package <...>) by
# ${package}_FOUND to allow easy custom overrides
#

MACRO(FIND_PACKAGE _package_name)
  STRING(TOUPPER ${_package_name} _package_name_uppercase)

  IF(NOT DEFINED ${_package_name_uppercase}_FOUND)
    _FIND_PACKAGE (${_package_name} ${ARGN})
  ENDIF()
ENDMACRO()
