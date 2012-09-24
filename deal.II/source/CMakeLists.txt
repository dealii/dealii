#
# Compile the deal.II library
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(distributed)
ADD_SUBDIRECTORY(dofs)
ADD_SUBDIRECTORY(fe)
ADD_SUBDIRECTORY(grid)
ADD_SUBDIRECTORY(hp)
ADD_SUBDIRECTORY(lac)
ADD_SUBDIRECTORY(multigrid)
ADD_SUBDIRECTORY(numerics)

IF(CMAKE_BUILD_TYPE MATCHES "Debug")

  #
  # Combine all debug OBJECT targets to a debug library:
  #
  ADD_LIBRARY(${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX}
    ${deal_ii_objects.g}
    )

  SET_TARGET_PROPERTIES(${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX}
    PROPERTIES
    VERSION ${VERSION}
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION ${VERSION}
    LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_DEBUG}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_DEBUG}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_DEBUG}"
    )

  TARGET_LINK_LIBRARIES(${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX}
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_EXTERNAL_LIBRARIES_DEBUG}
    )

  INSTALL(TARGETS ${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX}
    LIBRARY
    DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    COMPONENT library
    )

ENDIF()


IF(CMAKE_BUILD_TYPE MATCHES "Release")

  #
  # Combine all release OBJECT targets to a release library:
  #
  ADD_LIBRARY(${DEAL_II_BASE_NAME}
    ${deal_ii_objects}
    )

  SET_TARGET_PROPERTIES(${DEAL_II_BASE_NAME}
    PROPERTIES
    VERSION ${VERSION}
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION ${VERSION}
    LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_RELEASE}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_RELEASE}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_RELEASE}"
    )

  TARGET_LINK_LIBRARIES(${DEAL_II_BASE_NAME}
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_EXTERNAL_LIBRARIES_RELEASE}
    )

  INSTALL(TARGETS ${DEAL_II_BASE_NAME}
    LIBRARY
    DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    COMPONENT library
    )

ENDIF()

