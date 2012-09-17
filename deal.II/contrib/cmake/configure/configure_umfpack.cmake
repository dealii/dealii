#
# Configuration for the umfpack and amd libraries:
#


SET(FEATURE_UMFPACK_DEPENDS
  DEAL_II_WITH_BLAS
  DEAL_II_WITH_LAPACK
  )


MACRO(FEATURE_UMFPACK_FIND_EXTERNAL var)

  FIND_PACKAGE(UMFPACK)
  FIND_PACKAGE(AMD)

  IF(UMFPACK_FOUND AND AMD_FOUND)
    #SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(FEATURE_UMFPACK_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

  LIST(APPEND deal_ii_external_libraries
    ${UMFPACK_LIBRARY} ${AMD_LIBRARY}
    )

  SET(HAVE_LIBUMFPACK TRUE)

  SET(${var} TRUE)

ENDMACRO()


SET(FEATURE_UMFPACK_HAVE_CONTRIB TRUE)


MACRO(FEATURE_UMFPACK_CONFIGURE_CONTRIB var)

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

  SET(HAVE_LIBUMFPACK TRUE)

  SET(${var} TRUE)

ENDMACRO()


SET(FEATURE_UMFPACK_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_UMFPACK_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "
Could not find the umfpack and amd libraries!

Please ensure that the libraries are installed on your computer.
If the libraries are not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

Alternatively you may choose to compile the bundled contrib libraries
by setting DEAL_II_ALLOW_CONTRIB=on or DEAL_II_FORCE_CONTRIB_UMFPACK=on.

")
ENDMACRO()


CONFIGURE_FEATURE(UMFPACK)
