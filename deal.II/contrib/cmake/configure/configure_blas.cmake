#
# Configuration for the blas library:
#

# TODO:
# - include dir?
# - unit checks and definitions.

MACRO(FIND_FEATURE_BLAS_EXTERNAL var)

  FIND_PACKAGE(BLAS)

  IF(BLAS_FOUND)
    SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_BLAS_EXTERNAL var)

  SET(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    )

  LIST(APPEND deal_ii_external_libraries
    ${BLAS_LIBRARIES}
    )

  SET(HAVE_LIBBLAS TRUE)

  SET(${var} TRUE)
ENDMACRO()


MACRO(CONFIGURE_FEATURE_BLAS_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "
Could not find the blas library!

Please ensure that the blas library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

")
ENDMACRO()


CONFIGURE_FEATURE(BLAS)
