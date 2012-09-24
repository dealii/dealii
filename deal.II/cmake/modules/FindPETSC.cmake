#
# Try to find the petsc library
#

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(PETSC_INCLUDE_DIR petscconf.h
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}/${PETSC_ARCH}/include
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include
  PATH_SUFFIXES petsc
)

FIND_LIBRARY(PETSC_LIBRARY
  NAMES petsc
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}/${PETSC_ARCH}
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSC DEFAULT_MSG PETSC_LIBRARY PETSC_INCLUDE_DIR)

IF(PETSC_FOUND)
  MARK_AS_ADVANCED(
    PETSC_LIBRARY
    PETSC_INCLUDE_DIR
  )
ENDIF()

