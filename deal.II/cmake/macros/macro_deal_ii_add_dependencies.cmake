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
# A small wrapper around ADD_DEPENDENCIES to add the specified dependencies
# to every ${target}_${build} target, where build runs through all build
# types specified in DEAL_II_BUILD_TYPES
#

MACRO(DEAL_II_ADD_DEPENDENCIES name target)

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${build} build_lowercase)
    ADD_DEPENDENCIES(${name}.${build_lowercase}
      ${target}.${build_lowercase}
      )
  ENDFOREACH()

ENDMACRO()
