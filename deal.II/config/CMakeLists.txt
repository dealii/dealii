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
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
    )

  GET_CMAKE_PROPERTY(res VARIABLES)
  FOREACH(var ${res})
    IF(var MATCHES "DEAL_II_WITH")
      FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
        "SET(${var} ${${var}})\n"
        )
    ENDIF()
  ENDFOREACH()

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
    @ONLY
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
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

  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS "${CMAKE_INCLUDE_FLAG_C}"
    ${DEAL_II_USER_INCLUDE_DIRS}
    ${DEAL_II_INCLUDE_DIRS}
    )

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND MAKEFILE_BUILD_TYPES "DEBUG")
  ENDIF()
  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    LIST(APPEND MAKEFILE_BUILD_TYPES "RELEASE")
  ENDIF()

  FOREACH(type ${MAKEFILE_BUILD_TYPES})
    TO_STRING(MAKEFILE_DEFINITIONS_${type}
      ${DEAL_II_USER_DEFINITIONS}
      ${DEAL_II_USER_DEFINITIONS_${type}}
      )
    #
    # Add an rpath directive in front of each library, so that libraries
    # outside of the default search directories will be found by the runtime
    # dynamic linker:
    #
    FOREACH(lib
        ${DEAL_II_EXTERNAL_LIBRARIES}
        ${DEAL_II_EXTERNAL_LIBRARIES_${type}}
        ${DEAL_II_LIBRARIES_${type}})
      GET_FILENAME_COMPONENT(path ${lib} PATH)
      LIST(APPEND MAKEFILE_LIBS_${type} "-Wl,-rpath  -Wl,${path}")
      LIST(APPEND MAKEFILE_LIBS_${type} ${lib})
    ENDFOREACH()
    TO_STRING(MAKEFILE_LIBS_${type} ${MAKEFILE_LIBS_${type}})
    #
    # Put our linker flags directly in front of this string:
    #
    SET(MAKEFILE_LIBS_${type}
      "${CMAKE_SHARED_LINKER_FLAGS} ${DEAL_II_SHARED_LINKER_FLAGS_${type}} ${MAKEFILE_LIBS_${type}}"
      )
  ENDFOREACH()

  #
  # Boilerplate: The Make.global_options expects variables to be set to
  # yes, as is common for Makefiles.
  #
  COND_SET_TO_YES(DEAL_II_WITH_TBB MAKEFILE_enablethreads)
  COND_SET_TO_YES(DEAL_II_WITH_FUNCTIONPARSER MAKEFILE_enableparser)
  COND_SET_TO_YES(BUILD_SHARED_LIBS MAKEFILE_enableshared)

  COND_SET_TO_YES(DEAL_II_WITH_PETSC MAKEFILE_PETSC)
  COND_SET_TO_YES(DEAL_II_USE_PETSC_DEV MAKEFILE_PETSC_DEV)
  COND_SET_TO_YES(DEAL_II_WITH_TRILINOS MAKEFILE_TRILINOS)
  COND_SET_TO_YES(DEAL_II_WITH_BLAS MAKEFILE_BLAS)
  COND_SET_TO_YES(DEAL_II_WITH_LAPACK MAKEFILE_LAPACK)
  COND_SET_TO_YES(DEAL_II_WITH_ARPACK MAKEFILE_ARPACK)
  COND_SET_TO_YES(DEAL_II_WITH_METIS MAKEFILE_METIS)
  COND_SET_TO_YES(DEAL_II_WITH_UMFPACK MAKEFILE_UMFPACK)
  COND_SET_TO_YES(DEAL_II_WITH_P4EST MAKEFILE_P4EST)
  COND_SET_TO_YES(DEAL_II_WITH_MPI MAKEFILE_MPI)

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Version.in
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
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

