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
# If bool is "true" (in a cmake fashion...), set variable to "yes",
# otherwise to "no".
#
# Usage:
#     COND_SET_TO_YES(bool variable)
#

MACRO(COND_SET_TO_YES _bool _variable)
  IF(${_bool})
    SET(${_variable} "yes")
  ELSE()
    SET(${_variable} "no")
  ENDIF()
ENDMACRO()

