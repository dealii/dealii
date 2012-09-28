#
# Try to find the petsc library
#
# This module exports:
#
#     PETSC_FOUND
#     PETSC_LIBRARIES
#     PETSC_INCLUDE_DIRS
#     PETSC_VERSION
#     PETSC_VERSION_MAJOR
#     PETSC_VERSION_MINOR
#     PETSC_VERSION_SUBMINOR
#     PETSC_VERSION_PATCH
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(PETSC_DIR "$ENV{PETSC_DIR}")
SET_IF_EMPTY(PETSC_ARCH "$ENV{PETSC_ARCH}")

#
# So, well, yes. I'd like to include the PETScConfig.cmake file via
# FIND_PACKAGE(), but it is broken beyond belief:
#
# - In source, i.e. PETSC_DIR/PETSC_ARCH, it sets BUILD_SHARED_LIBS.
# - It does not contain its very own version number
# - It does not contain its very own library location(s) or name(s)
# - It does not contain necessary includes
#
# - It writes a lot of FIND_LIBRARY(..) statements. Seriously. What the
#   heck? If its not the same library you're linking against, you cannot
#   assume to be API compatible, so why not just give a list of libraries?
#
# - It is not even considered to be installed in the gentoo petsc package
#   :-]
#

#
# TODO: We'll have to guess which external libraries we'll have to link
# against someday to avoid underlinkage
#

FIND_LIBRARY(PETSC_LIBRARIES
  NAMES petsc
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)


#
# So, up to this point it was easy. Now, the tricky part:
#


#
# Search for the first part of the includes:
#
FIND_PATH(PETSC_INCLUDE_DIR_ARCH petscconf.h
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}/include
    ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES include petsc
)

#
# Sometimes, this is not enough...
# If petsc is not installed but in source tree layout, there will be
#   ${PETSC_DIR}/${PETSC_ARCH}/include - which we should have found by now.
#   ${PETSC_DIR}/include               - which we still have to find.
#
# Or it is installed in a non standard layout in the system (e.g. in
# Gentoo), where there will be
#   ${PETSC_DIR}/${PETSC_ARCH}/include
#   /usr/include/petsc ...
#
# Either way, petscversion.h should lie around:
#
FIND_PATH(PETSC_INCLUDE_DIR_COMMON petscversion.h
  HINTS
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}/include
    ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc
)

#
# And finally set PETSC_INCLUDE_DIRS depending on the outcome of our crude
# guess:
#
IF( PETSC_INCLUDE_DIR_ARCH MATCHES "-NOTFOUND" OR
    PETSC_INCLUDE_DIR_COMMON MATCHES "-NOTFOUND" )
  SET(PETSC_INCLUDE_DIRS "PETSC_INCLUDE_DIRS-NOTFOUND"
    CACHE STRING "Include paths for petsc"
    FORCE
    )
  UNSET(PETSC_INCLUDE_DIR_ARCH CACHE)
  UNSET(PETSC_INCLUDE_DIR_COMMON CACHE)
ELSE()
  UNSET(PETSC_INCLUDE_DIRS CACHE)
  SET(PETSC_INCLUDE_DIRS
    ${PETSC_INCLUDE_DIR_ARCH}
    ${PETSC_INCLUDE_DIR_COMMON}
    )
ENDIF()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSC DEFAULT_MSG
  PETSC_LIBRARIES
  PETSC_INCLUDE_DIRS
  )

IF(PETSC_FOUND)

  #
  # Is petsc compiled with support for MPIUNI?
  #
  FILE(STRINGS "${PETSC_PETSCCONF_H}" PETSC_MPIUNI_STRING
    REGEX "#define.*PETSC_HAVE_MPIUNI 1")
  IF("${PETSC_MPIUNI_STRING}" STREQUAL "")
    SET(PETSC_WITH_MPIUNI FALSE)
  ELSE()
    SET(PETSC_WITH_MPIUNI TRUE)
  ENDIF()

  IF(PETSC_WITH_MPIUNI)
    #
    # TODO: Still needed? We have a way bigger problem with underlinkage so
    #       far...
    #
    # If yes, add libmpiuni.so/a (if available)
    # We need to link with it on some systems where PETSc is built without
    # a real MPI and we need to handle trivial (one process) MPI
    # functionality.
    #
    FIND_LIBRARY(PETSC_LIBMPIUNI
      NAMES mpiuni
      HINTS
        ${PETSC_DIR}/${PETSC_ARCH}
      PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
      )
    IF(NOT PETSC_LIBMPIUNI MATCHES "-NOTFOUND")
      LIST(APPEND PETSC_LIBRARIES "${PETSC_LIBMPIUNI}")
    ELSE()
      SET(PETSC_LIBMPIUNI "")
    ENDIF()
    MARK_AS_ADVANCED(PETSC_LIBMPIUNI)
  ENDIF()


  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MAJOR_STRING
    REGEX "#define.*PETSC_VERSION_MAJOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MAJOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_MAJOR "${PETSC_VERSION_MAJOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MINOR_STRING
    REGEX "#define.*PETSC_VERSION_MINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MINOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_MINOR "${PETSC_VERSION_MINOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_SUBMINOR_STRING
    REGEX "#define.*PETSC_VERSION_SUBMINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_SUBMINOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_SUBMINOR "${PETSC_VERSION_SUBMINOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_PATCH_STRING
    REGEX "#define.*PETSC_VERSION_PATCH")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_PATCH.*([0-9]+).*" "\\1"
    PETSC_VERSION_PATCH "${PETSC_VERSION_PATCH_STRING}"
    )

  SET(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}")

  MARK_AS_ADVANCED(
    PETSC_LIBRARIES
    PETSC_INCLUDE_DIRS
    PETSC_DIR
    PETSC_ARCH
  )
ELSE()
  SET(PETSC_DIR "" CACHE STRING
    "An optional hint to a PETSc directory"
    )
  SET(PETSC_ARCH "" CACHE STRING
    "An optional hint to a PETSc arch"
    )
ENDIF()

