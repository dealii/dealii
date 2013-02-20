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
# Setup some convenience custom targets for the build system:
#

#
# Dump log files, configuration and system information to a convenient
# place:
#
ADD_CUSTOM_TARGET(log
  COMMAND ${CMAKE_COMMAND} --system-information ${CMAKE_BINARY_DIR}/cmake.log
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Create cmake.log"
  )


###########################################################################
#                                                                         #
#  Custom targets for library, documentation and compat_files components: #
#                                                                         #
###########################################################################

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
  ADD_CUSTOM_TARGET(compat_files)
  ADD_DEPENDENCIES(compat_files
    expand_instantiations
    make_dependencies
    report_features
    )
ENDIF()

IF(DEAL_II_COMPONENT_CONTRIB)
  ADD_CUSTOM_TARGET(contrib
    DEPENDS
      mesh_conversion
      parameter_gui
    )
ENDIF()

