
#
# Setup some convenience custom targets for the build system:
#

#
# Custom targets for building all targets necessary to install a specific
# component. (Too bad, we have to do this by hand. There is no cmake
# internal way to do this, yet...)
#

ADD_CUSTOM_TARGET(library)
ADD_DEPENDENCIES(library deal_II)


IF(DEAL_II_WITH_DOXYGEN)

  ADD_CUSTOM_TARGET(documentation)
  ADD_DEPENDENCIES(documentation doxygen)

ENDIF()


IF(DEAL_II_INSTALL_COMPAT_FILES)

  ADD_CUSTOM_TARGET(compat_files)

ENDIF()


IF(DEAL_II_INSTALL_PROJECT_CONFIG)

  ADD_CUSTOM_TARGET(project_config)

ENDIF()

