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
# This file implements the DEAL_II_SETUP_COMPILER macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_SETUP_COMPILER()
#

MACRO(DEAL_II_SETUP_TARGET target)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDE)
    MESSAGE(FATAL_ERROR
      "DEAL_II_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use."
      )
  ENDIF()

  SET(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER})

ENDMACRO()

