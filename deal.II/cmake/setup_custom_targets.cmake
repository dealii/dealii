## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Setup some convenience custom targets for the build system:
#


########################################################################
#                                                                      #
#   Custom targets for library, documentation and compat_files comp.:  #
#                                                                      #
########################################################################

ADD_CUSTOM_TARGET(library)
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  ADD_DEPENDENCIES(library ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})
ENDFOREACH()

IF(DEAL_II_COMPONENT_DOCUMENTATION)
  ADD_CUSTOM_TARGET(documentation
    DEPENDS doxygen
    )
ENDIF()

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  ADD_CUSTOM_TARGET(compat_files
    DEPENDS
      expand_instantiations
      make_dependencies
      report_features
    )
ENDIF()
