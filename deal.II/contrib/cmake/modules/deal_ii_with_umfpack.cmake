FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(BLAS REQUIRED)

IF(DEAL_II_USE_CONTRIB)
  # Compiles and links libumfpack, exports the Umfpack_* variables as well:
  ADD_SUBDIRECTORY(contrib/umfpack)
ELSE()
  FIND_PACKAGE(Umfpack REQUIRED)
  FIND_PACKAGE(AMD REQUIRED)
ENDIF()

INCLUDE_DIRECTORIES(${Umfpack_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

SET(deal_ii_external_libraries
  ${deal_ii_external_libraries}
  ${Umfpack_LIBRARY} ${AMD_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
  )

SET(deal_ii_external_debug_libraries
  ${deal_ii_external_debug_libraries}
  ${Umfpack_LIBRARY} ${AMD_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
  )
