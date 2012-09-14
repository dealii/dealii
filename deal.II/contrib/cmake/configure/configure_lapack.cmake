#
# Configuration for the lapack library:
#

# TODO:
# - include dir?
# - unit checks and definitions.

MACRO(FIND_FEATURE_LAPACK_EXTERNAL var)

  FIND_PACKAGE(LAPACK)

  IF(LAPACK_FOUND)
    SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_LAPACK_EXTERNAL var)

  SET(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}"
    )

  LIST(APPEND deal_ii_external_libraries
    ${LAPACK_LIBRARIES}
    )

  SET(HAVE_LIBLAPACK TRUE)

  SET(${var} TRUE)

ENDMACRO()


MACRO(CONFIGURE_FEATURE_BLAS_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "
Could not find the lapack library!

Please ensure that the lapack library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

")
ENDMACRO()


CONFIGURE_FEATURE(LAPACK)
