
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

ADD_CUSTOM_TARGET(library)

FOREACH(build ${DEAL_II_BUILD_TYPES})
  ADD_DEPENDENCIES(library ${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX})
ENDFOREACH()

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

IF(DEAL_II_COMPONENT_CONTRIB) 
  ADD_CUSTOM_TARGET(contrib)
  ADD_DEPENDENCIES(contrib
    mesh_conversion
    )
ENDIF()

