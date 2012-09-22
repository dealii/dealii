CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )


IF(DEAL_II_COMPONENT_PROJECT_CONFIG)
  #
  # Configure and install the cmake project configuration:
  #

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_NAME}Config.cmake
    )

  GET_CMAKE_PROPERTY(res VARIABLES)
  FOREACH(var ${res})
    IF(var MATCHES "DEAL_II_WITH")
      FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_NAME}Config.cmake
        "SET(${var} ${${var}})\n"
        )
    ENDIF()
  ENDFOREACH()

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_NAME}ConfigVersion.cmake
    @ONLY
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
    COMPONENT project_config
    )
ENDIF()


IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Configure and install the old Makefile configuration:
  #

  #
  # Therefore, transform some cmake lists into a string that a Makefile
  # actually understands:
  #

  TO_STRING(MAKEFILE_LIBRARIES
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_LIBRARIES}
    )

  TO_STRING(MAKEFILE_DEFINITIONS
    ${DEAL_II_EXTERNAL_DEFINITIONS}
    )

  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS
    "${CMAKE_INCLUDE_FLAG_C}"
    ${DEAL_II_EXTERNAL_INCLUDE_DIRS}
    ${DEAL_II_INCLUDE_DIRS}
    )


  #
  # Add an rpath directive so that libraries outside of the default search
  # directories will be found by the runtime dynamic linker:
  #
  FOREACH(lib ${DEAL_II_EXTERNAL_LIBRARIES} ${DEAL_II_LIBRARIES})
    GET_FILENAME_COMPONENT(path ${lib} PATH)
    LIST(APPEND MAKEFILE_RPATHFLAGS "-Wl,-rpath  -Wl,${path}")
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES MAKEFILE_RPATHFLAGS)
  TO_STRING(MAKEFILE_RPATHFLAGS ${MAKEFILE_RPATHFLAGS})

  SET(MAKEFILE_LDFLAGS "${CMAKE_SHARED_LINKER_FLAGS}")

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options.${CMAKE_BUILD_TYPE}
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.stub
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    COPYONLY)

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Version.in
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options.${CMAKE_BUILD_TYPE}
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    DESTINATION common
    COMPONENT compat_files
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    COMPONENT compat_files
    )

ENDIF()

