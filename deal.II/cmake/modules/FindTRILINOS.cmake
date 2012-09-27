#
# Try to find the Trilinos library
#
# This module exports:
#
#   TRILINOS_DIR
#   TRILINOS_INCLUDE_DIRS
#   TRILINOS_LIBRARY_*
#   TRILINOS_LIBRARIES
#   TRILINOS_VERSION_MAJOR
#   TRILINOS_VERSION_MINOR
#   TRILINOS_VERSION_SUBMINOR
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(TRILINOS_DIR "$ENV{TRILINOS_DIR}")

#
# Include the trilinos package configuration:
#
FIND_PACKAGE(TRILINOS
  QUIET CONFIG
  NAMES Trilinos TRILINOS
  HINTS
    ${TRILINOS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX}/cmake/Trilinos
    lib64/cmake/Trilinos
    lib/cmake/Trilinos
  )

#
# Extract the major and minor version numbers:
#
STRING(REGEX REPLACE
  "^([0-9]+).*$" "\\1"
  TRILINOS_VERSION_MAJOR "${Trilinos_VERSION}")

STRING(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*$" "\\1"
  TRILINOS_VERSION_MINOR "${Trilinos_VERSION}")

STRING(REGEX REPLACE
  "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
  TRILINOS_VERSION_SUBMINOR "${Trilinos_VERSION}")


SET(TRILINOS_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS})

#
# We'd like to have the full library names but the Trilinos package only
# exports a list with short names...
# So we check again for every lib and store the full path:
#
FOREACH(library ${Trilinos_LIBRARIES})
  FIND_LIBRARY(TRILINOS_LIBRARY_${library}
    NAMES ${library}
    HINTS ${Trilinos_LIBRARY_DIRS}
    )

  LIST(APPEND TRILINOS_LIBRARIES ${TRILINOS_LIBRARY_${library}})

  #
  # Remove the variables from the cache, so that updating TRILINOS_DIR will
  # find the new libraries..
  #
  UNSET(TRILINOS_LIBRARY_${library} CACHE)
ENDFOREACH()


FIND_PACKAGE_HANDLE_STANDARD_ARGS(TRILINOS DEFAULT_MSG
  TRILINOS_DIR
  TRILINOS_INCLUDE_DIRS
  TRILINOS_LIBRARIES
  )

IF(TRILINOS_FOUND)
  MARK_AS_ADVANCED(
    TRILINOS_DIR
    TRILINOS_INCLUDE_DIRS
    TRILINOS_LIBRARIES
    )
ELSE()
  IF(TRILINOS_DIR MATCHES "-NOTFOUND")
    SET(TRILINOS_DIR "" CACHE STRING
      "A hint to a Trilinos installation"
      FORCE
      )
  ENDIF()
ENDIF()

