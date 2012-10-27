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


###########################################################################
#                                                                         #
#                  Custom targets for the autopilot mode:                 #
#                                                                         #
###########################################################################

IF( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" AND
    NOT DISABLE_AUTOPILOT)
  #
  # Setup the "distclean" target:
  #
  ADD_CUSTOM_TARGET(distclean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/distclean.cmake
    )

  #
  # Targets for changing the build type:
  #
  ADD_CUSTOM_TARGET(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
    )
  ADD_CUSTOM_TARGET(release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMENT "Switch CMAKE_BUILD_TYPE to Release"
    )
  ADD_CUSTOM_TARGET(debugrelease
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=DebugRelease ${CMAKE_SOURCE_DIR}
    COMMENT "Switch CMAKE_BUILD_TYPE to DebugRelease"
    )

  #
  # A cheesy trick to automatically install and print a summary
  #
  ADD_CUSTOM_TARGET("install-intree" ALL
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install/fast
    )
  ADD_DEPENDENCIES(install-intree compat_files)
  ADD_DEPENDENCIES(install-intree contrib)
  ADD_DEPENDENCIES(install-intree documentation)
  ADD_DEPENDENCIES(install-intree library)
ENDIF()

