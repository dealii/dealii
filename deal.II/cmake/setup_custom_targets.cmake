
#
# Setup some convenience custom targets for the build system:
#

#
# Custom targets for building all targets necessary to install a specific
# component. (Too bad, we have to do this by hand. There is no cmake
# internal way to do this, yet...)
#

ADD_CUSTOM_TARGET(library)
ADD_DEPENDENCIES(library ${DEAL_II_BASE_NAME})


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

