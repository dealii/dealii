
ADD_EXECUTABLE(expand_instantiations expand_instantiations.cc)


IF(DEAL_II_COMPONENT_COMPAT_FILES)

  #
  # Build and install the old dependency_resolution and report_features
  # executables:
  #

  #INCLUDE_DIRECTORIES(
  #  ${CMAKE_SOURCE_DIR}/include/
  #  ${CMAKE_BINARY_DIR}/include/
  #  )

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)

  #ADD_EXECUTABLE(report_features report_features.cc)
  #TARGET_LINK_LIBRARIES(report_features ${DEAL_II_EXTERNAL_LIBRARIES})

  INSTALL(TARGETS
    make_dependencies
    #report_features
    expand_instantiations
    RUNTIME
    DESTINATION common/scripts
    COMPONENT compat_files
    )

ENDIF()

