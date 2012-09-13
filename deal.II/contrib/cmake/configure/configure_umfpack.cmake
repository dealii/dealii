FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(BLAS REQUIRED)

# TODO: A deal.II specific error message if blas or lapack is not found

IF(NOT DEAL_II_FORCE_CONTRIB_UMFPACK)

  FIND_PACKAGE(UMFPACK)
  FIND_PACKAGE(AMD)

  IF(NOT DEAL_II_ALLOW_CONTRIB)
    IF(NOT UMFPACK_FOUND)
      macro_message_not_found("umfpack" "UMFPACK")
    ENDIF()
    IF(NOT AMD_FOUND)
      macro_message_not_found("amd" "AMD")
    ENDIF()
  ENDIF()

ENDIF()

IF(UMFPACK_FOUND AND AMD_FOUND) # TODO
  SET(UMFPACKAMD_FOUND TRUE)
ELSE()
  SET(UMFPACKAMD_FOUND FALSE)
ENDIF()

IF(DEAL_II_FORCE_CONTRIB_UMFPACK OR NOT UMFPACKAMD_FOUND)

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/umfpack/UMFPACK/Include
    ${CMAKE_SOURCE_DIR}/contrib/umfpack/AMD/Include
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/umfpack/UMFPACK/Source)
  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/umfpack/AMD/Source)

  #
  # Add umfpack and amd directly to the object files of deal.II
  #
  LIST(APPEND deal_ii_additional_object_files
    ${obj_umfpack_object_files}
    $<TARGET_OBJECTS:obj_amd_int>
    $<TARGET_OBJECTS:obj_amd_long>
    $<TARGET_OBJECTS:obj_amd_global>
    )

ELSE()

  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

  LIST(APPEND deal_ii_required_linker_flags
    ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}
    )

  LIST(APPEND deal_ii_external_libraries
    ${UMFPACK_LIBRARY} ${AMD_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    )

ENDIF()


SET(HAVE_LIBBLAS TRUE)
SET(HAVE_LIBLAPACK TRUE)
SET(HAVE_LIBUMFPACK TRUE)
