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
# A small wrapper around
# SET_TARGET_PROPERTY(... PROPERTIES COMPILE_DEFINITIONS ...)
# to _add_ compile definitions to every target we have specified.
#

MACRO(DEAL_II_ADD_DEFINITIONS name)

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${build} build_lowercase)

    GET_TARGET_PROPERTY(macro_definitions ${name}.${build_lowercase} COMPILE_DEFINITIONS)
    SET_TARGET_PROPERTIES(${name}.${build_lowercase} PROPERTIES
      COMPILE_DEFINITIONS "${ARGN};${macro_definitions}"
      )
  ENDFOREACH()

ENDMACRO()

