## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

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
#     PETSC_WITH_64BIT_INDICES
#     PETSC_WITH_COMPLEX
#     PETSC_WITH_HYPRE
#     PETSC_WITH_KOKKOS
#     PETSC_WITH_MPIUNI
#     PETSC_WITH_MUMPS
#

set(PETSC_DIR "" CACHE PATH "An optional hint to a PETSc directory")
set(PETSC_ARCH "" CACHE STRING "An optional hint to a PETSc arch")
set_if_empty(PETSC_DIR "$ENV{PETSC_DIR}")
set_if_empty(PETSC_ARCH "$ENV{PETSC_ARCH}")

deal_ii_find_library(PETSC_LIBRARY
  NAMES petsc libpetsc
  HINTS ${PETSC_DIR} ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Search for the first part of the includes:
#

deal_ii_find_path(PETSC_INCLUDE_DIR_ARCH petscconf.h
  HINTS ${PETSC_DIR} ${PETSC_DIR}/${PETSC_ARCH} ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc include include/petsc
)

set(PETSC_PETSCCONF_H "${PETSC_INCLUDE_DIR_ARCH}/petscconf.h")

macro(_petsc_feature_check _var _regex)
  file(STRINGS "${PETSC_PETSCCONF_H}" PETSC_${_var}_STRING
    REGEX "${_regex}")
  if("${PETSC_${_var}_STRING}" STREQUAL "")
    set(PETSC_WITH_${_var} FALSE)
  else()
    set(PETSC_WITH_${_var} TRUE)
  endif()
endmacro()

if(EXISTS ${PETSC_PETSCCONF_H})
  _petsc_feature_check(64BIT_INDICES "#define.*PETSC_USE_64BIT_INDICES 1")
  _petsc_feature_check(COMPLEX "#define.*PETSC_USE_COMPLEX 1")
  _petsc_feature_check(HYPRE "#define.*PETSC_HAVE_HYPRE 1")
  _petsc_feature_check(KOKKOS "#define.*PETSC_HAVE_KOKKOS 1")
  _petsc_feature_check(MPIUNI "#define.*PETSC_HAVE_MPIUNI 1")
  _petsc_feature_check(MUMPS "#define.*PETSC_HAVE_MUMPS 1")
endif()

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
# Either way, we must be able to find petscversion.h:
#

deal_ii_find_path(PETSC_INCLUDE_DIR_COMMON petscversion.h
  HINTS ${PETSC_DIR} ${PETSC_DIR}/${PETSC_ARCH} ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc include include/petsc
)

set(PETSC_PETSCVERSION_H "${PETSC_INCLUDE_DIR_COMMON}/petscversion.h")
if(EXISTS ${PETSC_PETSCVERSION_H})
  file(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MAJOR_STRING
    REGEX "^#[ \t]*define[ \t]+PETSC_VERSION_MAJOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PETSC_VERSION_MAJOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    PETSC_VERSION_MAJOR "${PETSC_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MINOR_STRING
    REGEX "^#[ \t]*define[ \t]+PETSC_VERSION_MINOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PETSC_VERSION_MINOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    PETSC_VERSION_MINOR "${PETSC_VERSION_MINOR_STRING}"
    )
  file(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_SUBMINOR_STRING
    REGEX "^#[ \t]*define[ \t]+PETSC_VERSION_SUBMINOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PETSC_VERSION_SUBMINOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    PETSC_VERSION_SUBMINOR "${PETSC_VERSION_SUBMINOR_STRING}"
    )
  file(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_PATCH_STRING
    REGEX "^#[ \t]*define[ \t]+PETSC_VERSION_PATCH[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PETSC_VERSION_PATCH[ \t]+([0-9]+)[ \t]*$" "\\1"
    PETSC_VERSION_PATCH "${PETSC_VERSION_PATCH_STRING}"
    )
  set(PETSC_VERSION
    "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.${PETSC_VERSION_PATCH}"
    )
endif()

#
# So, up to this point it was easy. Now, the tricky part. Search for
# petscvariables and determine the includes and the link interface from
# that file:
#

deal_ii_find_file(PETSC_PETSCVARIABLES
  NAMES petscvariables
  HINTS ${PETSC_DIR}/${PETSC_ARCH} ${PETSC_DIR}
  PATH_SUFFIXES conf lib/petsc/conf
  )

if(NOT PETSC_PETSCVARIABLES MATCHES "-NOTFOUND")
  #
  # Includes:
  #

  file(STRINGS "${PETSC_PETSCVARIABLES}" _external_includes
    REGEX "^PETSC_CC_INCLUDES =.*")
  separate_arguments(_external_includes)

  set(_petsc_includes)
  foreach(_token ${_external_includes})
    #
    # workaround: Do not pull in scotch include directory. It clashes with
    # our use of the metis headers...
    #
    if(_token MATCHES "^-I" AND NOT _token MATCHES "scotch$")
      string(REGEX REPLACE "^-I" "" _token "${_token}")
      list(APPEND _petsc_includes ${_token})
    endif()
  endforeach()

  #
  # Link line:
  #

  file(STRINGS "${PETSC_PETSCVARIABLES}" PETSC_EXTERNAL_LINK_LINE
    REGEX "^PETSC_WITH_EXTERNAL_LIB =.*")

  separate_arguments(PETSC_EXTERNAL_LINK_LINE)

  set(_hints)
  set(_petsc_libraries)
  set(_cleanup_variables)
  foreach(_token ${PETSC_EXTERNAL_LINK_LINE})
    if(_token MATCHES "^-L")
      # Build up hints with the help of all tokens passed with -L:
      string(REGEX REPLACE "^-L" "" _token "${_token}")
      list(APPEND _hints ${_token})
    elseif(_token MATCHES "^-l")
      # Search for every library that was specified with -l:
      string(REGEX REPLACE "^-l" "" _token "${_token}")

      if(NOT _token MATCHES "(petsc|stdc\\+\\+|gcc_s|clang_rt)")
        list(APPEND _cleanup_variables PETSC_LIBRARY_${_token})

        if(_token MATCHES "^(c|quadmath|gfortran|m|rt|nsl|dl|pthread)$")
          find_system_library(PETSC_LIBRARY_${_token} NAMES ${_token})
        else()
          deal_ii_find_library(PETSC_LIBRARY_${_token}
            NAMES ${_token}
            HINTS ${_hints}
            )
        endif()
        if(NOT PETSC_LIBRARY_${_token} MATCHES "-NOTFOUND")
          list(APPEND _petsc_libraries ${PETSC_LIBRARY_${_token}})
        endif()

      endif()

    endif()
  endforeach()
endif()

if(PETSC_WITH_MPIUNI)
  #
  # Workaround: Some distributions happen to not install petscvariables and
  # we consequently might miss some essential include directories. Let's
  # try at least to find the mpiuni include directory.
  #
  deal_ii_find_path(PETSC_INCLUDE_DIR_MPIUNI mpiuni/mpi.h
    HINTS ${PETSC_INCLUDE_DIR_COMMON} ${PETSC_INCLUDE_DIR_ARCH} ${_petsc_includes}
    PATH_SUFFIXES petsc
    )
  set(PETSC_INCLUDE_DIR_MPIUNI "${PETSC_INCLUDE_DIR_MPIUNI}/mpiuni")
endif()

process_feature(PETSC
  LIBRARIES
    REQUIRED PETSC_LIBRARY
    OPTIONAL _petsc_libraries
  INCLUDE_DIRS
    REQUIRED PETSC_INCLUDE_DIR_COMMON PETSC_INCLUDE_DIR_ARCH
    OPTIONAL PETSC_INCLUDE_DIR_MPIUNI _petsc_includes
  CLEAR
    PETSC_LIBRARY PETSC_INCLUDE_DIR_COMMON PETSC_INCLUDE_DIR_ARCH
    PETSC_PETSCVARIABLES ${_cleanup_variables}
  )

if(PETSC_FOUND)
  mark_as_advanced(PETSC_ARCH)
else()
  mark_as_advanced(CLEAR PETSC_ARCH)
endif()
