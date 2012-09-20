#
# Try to find the Trilinos library
#

INCLUDE(FindPackageHandleStandardArgs)

#
# Include the trilinos package configuration:
#
find_package(TRILINOS
  CONFIG
  NAMES Trilinos TRILINOS
  HINTS ${TRILINOS_DIR}
  PATH_PREFIXES lib64/cmake/Trilinos lib/cmake/Trilinos
  )

#
# Extract the major and minor version numbers:
#
STRING(REGEX REPLACE
  "^([0-9]+).*$" "\\1"
  TRILINOS_MAJOR "${Trilinos_VERSION}")

STRING(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*$" "\\1"
  TRILINOS_MINOR "${Trilinos_VERSION}")

STRING(REGEX REPLACE
  "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
  TRILINOS_SUBMINOR "${Trilinos_VERSION}")


FIND_PATH(TRILINOS_INCLUDE_DIR Trilinos_version.h
  HINTS ${Trilinos_INCLUDE_DIRS}
  )

#
# We'd like to have the full library names but the Trilinos package only
# exports a list with short names...
# So we check again for every lib and store the full path:
#
FOREACH(macro_library ${Trilinos_LIBRARIES})
  FIND_LIBRARY(TRILINOS_LIBRARY_${macro_library}
    NAMES ${macro_library}
    HINTS ${Trilinos_LIBRARY_DIRS}
    )
  MARK_AS_ADVANCED(TRILINOS_LIBRARY_${macro_library})

  LIST(APPEND TRILINOS_LIBRARIES ${TRILINOS_LIBRARY_${macro_library}})
ENDFOREACH()


FIND_PACKAGE_HANDLE_STANDARD_ARGS(TRILINOS DEFAULT_MSG
  TRILINOS_DIR
  TRILINOS_INCLUDE_DIR
  TRILINOS_LIBRARIES
  )

IF(TRILINOS_FOUND)
  MARK_AS_ADVANCED(
    TRILINOS_DIR
    TRILINOS_INCLUDE_DIR
    TRILINOS_LIBRARIES
    )
ENDIF()

