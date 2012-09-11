FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(BLAS REQUIRED)

IF(DEAL_II_ALLOW_CONTRIB)
  FIND_PACKAGE(Umfpack)
  FIND_PACKAGE(AMD)
ELSE()
  FIND_PACKAGE(Umfpack REQUIRED)
  FIND_PACKAGE(AMD REQUIRED)
ENDIF()

IF(UMFPACK_FOUND AND AMD_FOUND)
  SET(UmfpackAMD_FOUND TRUE)
ELSE()
  SET(UmfpackAMD_FOUND FALSE)
ENDIF()


IF(DEAL_II_FORCE_CONTRIB_UMFPACK OR NOT UmfpackAMD_FOUND)

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

  INCLUDE_DIRECTORIES(${Umfpack_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

  # TODO: Set necessary linker flags...

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
