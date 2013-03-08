#####
##
## Copyright (C) 2012 by the deal.II authors
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
# If 'variable' is empty it will be set to 'value'
#
MACRO(SET_IF_EMPTY _variable _value)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${_value})
  ENDIF()
ENDMACRO()

