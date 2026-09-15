## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Try to find the GMSH library
#
# This module exports
#
#   GMSH_EXECUTABLE
#   GMSH_LIBRARY
#   GMSH_INCLUDE_DIR
#   GMSH_WITH_API
#   GMSH_WITH_API_VERSION_MAJOR
#   GMSH_WITH_API_VERSION_MINOR
#   GMSH_WITH_API_VERSION_PATCH
#   GMSH_WITH_API_VERSION
#

set(GMSH_DIR "" CACHE PATH "An optional hint to a Gmsh installation containing the gmsh executable")
set_if_empty(GMSH_DIR "$ENV{GMSH_DIR}")

set(GMSH_LIBRARY_DIR "" CACHE PATH "An optional hint to a Gmsh SDK installation")
set_if_empty(GMSH_LIBRARY_DIR "${GMSH_DIR}")

deal_ii_find_program(GMSH_EXE gmsh${CMAKE_EXECUTABLE_SUFFIX}
  HINTS ${GMSH_DIR}
  PATH_SUFFIXES bin
  )

deal_ii_find_library(GMSH_LIBRARY
  NAMES gmsh
  HINTS ${GMSH_LIBRARY_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_path(GMSH_INCLUDE_DIR gmsh.h
  HINTS ${GMSH_LIBRARY_DIR}
  PATH_SUFFIXES include
  )

if(GMSH_LIBRARY MATCHES "-NOTFOUND" OR 
   GMSH_INCLUDE_DIR MATCHES "-NOTFOUND")
  set(GMSH_WITH_API FALSE)
else()
  set(GMSH_WITH_API TRUE)
endif()

if(GMSH_WITH_API)
  set(GMSH_WITH_API_CONFIG_H "${GMSH_INCLUDE_DIR}/gmsh.h")
  # version 4.4 and newer of gmsh have GMSH_API_VERSION_MAJOR etc., but older
  # versions just define a version string: for compatibility just unpack the
  # version string
  file(STRINGS "${GMSH_WITH_API_CONFIG_H}" GMSH_WITH_API_VERSION_STRING
    REGEX "#define.*GMSH_API_VERSION ")

  # The first two group matches must be the major and minor versions. Older
  # versions of gmsh don't define a PATCH version so that match is optional:
  string(REGEX MATCH "GMSH_API_VERSION.*\"([0-9]+)\\.([0-9]+)\\.?([0-9]+)?\""
    _ "${GMSH_WITH_API_VERSION_STRING}")
  set(GMSH_WITH_API_VERSION_MAJOR "${CMAKE_MATCH_1}")
  set(GMSH_WITH_API_VERSION_MINOR "${CMAKE_MATCH_2}")
  if("${CMAKE_MATCH_3}" STREQUAL "")
    set(GMSH_WITH_API_VERSION_PATCH "0")
  else()
    set(GMSH_WITH_API_VERSION_PATCH "${CMAKE_MATCH_3}")
  endif()
endif()


process_feature(GMSH
  EXECUTABLE REQUIRED GMSH_EXE
  LIBRARIES OPTIONAL GMSH_LIBRARY
  INCLUDE_DIRS OPTIONAL GMSH_INCLUDE_DIR
  CLEAR
    GMSH_EXE GMSH_LIBRARY GMSH_INCLUDE_DIR GMSH_WITH_API
    DEAL_II_GMSH_WITH_BOOST_PROCESS
  )
