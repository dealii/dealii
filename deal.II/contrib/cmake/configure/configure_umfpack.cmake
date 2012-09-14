#
# Configuration for the umfpack and amd libraries:
#


MACRO(FIND_FEATURE_UMFPACK_EXTERNAL var)

  FIND_PACKAGE(UMFPACK)
  FIND_PACKAGE(AMD)

  IF(UMFPACK_FOUND AND AMD_FOUND)
    SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_UMFPACK_EXTERNAL var)

  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

  LIST(APPEND deal_ii_external_libraries
    ${UMFPACK_LIBRARY} ${AMD_LIBRARY}
    )

  SET(HAVE_LIBUMFPACK TRUE)

  SET(${var} TRUE)

ENDMACRO()


SET(HAVE_CONTRIB_FEATURE_UMFPACK TRUE)


MACRO(CONFIGURE_FEATURE_UMFPACK_CONTRIB var)

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


MACRO(CONFIGURE_FEATURE_BOOST_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "
Could not find the umfpack and amd libraries!

Please ensure that the libraries are installed on your computer.
If the libraries are not at a default location, either provide some hints
via environment variables:
UMFPACK_LIBRARY_DIR UMFPACK_INCLUDE_DIR
AMD_LIBRARY_DIR AMD_INCLUDE_DIR
Or set the relevant variables by hand in ccmake.

Alternatively you may choose to compile the bundled contrib libraries
by setting DEAL_II_ALLOW_CONTRIB=on or DEAL_II_FORCE_CONTRIB_UMFPACK=on.

")
ENDMACRO()

#
# UMFPACK needs BLAS and LAPACK to be configured:
#
IF(DEAL_II_WITH_BLAS AND DEAL_II_WITH_LAPACK)
  CONFIGURE_FEATURE(UMFPACK)
ELSE()
  IF(DEAL_II_FEATURE_AUTODETECT)
    MESSAGE(STATUS
      "DEAL_II_WITH_UMFPACK has unmet configuration requirements: Both, DEAL_II_WITH_BLAS and DEAL_II_WITH_LAPACK have to be set."
      )
    SET_CACHED_OPTION(DEAL_II_WITH_UMFPACK OFF)
  ELSE()
    MESSAGE(SEND_ERROR
      "DEAL_II_WITH_UMFPACK has unmet configuration requirements: Both, DEAL_II_WITH_BLAS and DEAL_II_WITH_LAPACK have to be set."
      )
  ENDIF()
ENDIF()
