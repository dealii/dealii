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
# Tests whether the cxx compiler understands a flag.
# If so, add it to 'variable'.
#
# Usage:
#     ENABLE_IF_SUPPORTED(variable flag)
#

MACRO(ENABLE_IF_SUPPORTED variable flag)
  CHECK_CXX_COMPILER_FLAG(
    "${flag}"
    DEAL_II_HAVE_FLAG_${flag}
    )
  IF(DEAL_II_HAVE_FLAG_${flag})
    SET(${variable} "${${variable}} ${flag}")
  ENDIF()
ENDMACRO()

