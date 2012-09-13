# Try to find Umfpack

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(Umfpack_INCLUDE_DIR umfpack.h
  HINTS
  $ENV{Umfpack_INCLUDE_DIR}
  ${Umfpack_INCLUDE_DIR}
)

FIND_LIBRARY(Umfpack_LIBRARY
  NAMES libumfpack.so
  PATHS
  $ENV{Umfpack_LIBRARY_DIR}
  ${Umfpack_LIBRARY_DIR}
  PATH_SUFFIXES lib64 lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Umfpack DEFAULT_MSG Umfpack_LIBRARY Umfpack_INCLUDE_DIR)

IF(Umfpack_FOUND)
  MARK_AS_ADVANCED(
    Umfpack_LIBRARY
    Umfpack_INCLUDE_DIR
  )
ENDIF()
