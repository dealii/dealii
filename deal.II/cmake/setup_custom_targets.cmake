
#
# Setup some convenience custom targets for the build system, i.e.
#
#   $ make <custom_target>.
#
# We add custom targets for building all targets necessary to install a
# specific component (too bad, we have to do this by hand. There is no cmake
# internal way to do this, yet...):
#
#   library, documentation, compat_files, project_config
#
# And a release and debug target (depending on configuration)
#

ADD_CUSTOM_TARGET(library)

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_CUSTOM_TARGET(debug)

  ADD_DEPENDENCIES(library ${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX})
  ADD_DEPENDENCIES(debug ${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX})
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  ADD_CUSTOM_TARGET(release)

  ADD_DEPENDENCIES(library ${DEAL_II_BASE_NAME})
  ADD_DEPENDENCIES(release ${DEAL_II_BASE_NAME})
ENDIF()


IF(DEAL_II_COMPONENT_DOCUMENTATION)

  ADD_CUSTOM_TARGET(documentation)
  ADD_DEPENDENCIES(documentation doxygen)

ENDIF()


IF(DEAL_II_COMPONENT_COMPAT_FILES)

  ADD_CUSTOM_TARGET(compat_files)
  ADD_DEPENDENCIES(compat_files
    expand_instantiations
    make_dependencies
    report_features
    )

ENDIF()


IF(DEAL_II_COMPONENT_PROJECT_CONFIG)

  ADD_CUSTOM_TARGET(project_config)

ENDIF()

