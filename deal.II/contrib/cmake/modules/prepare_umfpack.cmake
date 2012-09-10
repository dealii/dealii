FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(BLAS REQUIRED)

IF(DEAL_II_USE_CONTRIB)
  #
  # Add umfpack and amd directly to the object files of deal.II
  #

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/umfpack/UMFPACK/Include
    ${CMAKE_SOURCE_DIR}/contrib/umfpack/AMD/Include
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/umfpack/UMFPACK/Source)
  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/umfpack/AMD/Source)

  LIST(APPEND deal_ii_additional_object_files
    ${obj_umfpack_object_files}
    $<TARGET_OBJECTS:obj_amd_int>
    $<TARGET_OBJECTS:obj_amd_long>
    $<TARGET_OBJECTS:obj_amd_global>
    )
ELSE()
  FIND_PACKAGE(Umfpack REQUIRED)
  FIND_PACKAGE(AMD REQUIRED)

  INCLUDE_DIRECTORIES(${Umfpack_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

  #
  # We skip *_INCLUDE_DIR because it is not needed for the use of the
  # deal.II library
  #

  LIST(APPEND deal_ii_external_libraries
    ${Umfpack_LIBRARY} ${AMD_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    )

  LIST(APPEND deal_ii_external_debug_libraries
    ${Umfpack_LIBRARY} ${AMD_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    )
ENDIF()


SET(HAVE_LIBBLAS TRUE)
SET(HAVE_LIBLAPACK TRUE)
SET(HAVE_LIBUMFPACK TRUE)
