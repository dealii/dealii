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
# This file implements the DEAL_II_IMPORT_LIBRARY macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_IMPORT_LIBRARY()
#
# This sets some cached variables to the values used for compiling the
# deal.II library.
#
# This macro has to be called before PROJECT()!
#

MACRO(DEAL_II_IMPORT_LIBRARY)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use.\n\n"
      )
  ENDIF()

  INCLUDE(${DEAL_II_TARGET_CONFIG})

  #
  # Setup a ${DEAL_II_BASE_NAME} target and populate it with our
  #

ENDMACRO()

