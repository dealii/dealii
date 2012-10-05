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
      "DEAL_II_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use."
      )
  ENDIF()

  INCLUDE(${DEAL_II_TARGET_CONFIG})

  #
  # Fixup the CONFIGURATION types for debug and release targets:
  #
  FOREACH(build ${DEAL_II_BUILD_TYPES})
    SET_TARGET_PROPERTIES(${DEAL_II_TARGET_${build}}
      PROPERTIES
      MAP_IMPORTED_CONFIG ${build}
      )
  ENDFOREACH()

ENDMACRO()

