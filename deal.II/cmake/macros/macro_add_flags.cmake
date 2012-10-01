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
# A small macro used for (string-)appending a string "${flags}" to a
# string "${variable}"
#
# Usage:
#     ADD_FLAGS(variable flags)
#

MACRO(ADD_FLAGS variable flags)
  SET(${variable} "${${variable}} ${flags}")
ENDMACRO()

