## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Try to find the GMSH library
#
# This module exports
#
#   GMSH_EXECUTABLE
#   GMSH_LIBRARY
#   GMSH_INCLUDE_DIR
#   GMSH_WITH_API
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

process_feature(GMSH
  EXECUTABLE REQUIRED GMSH_EXE
  LIBRARIES OPTIONAL GMSH_LIBRARY
  INCLUDE_DIRS OPTIONAL GMSH_INCLUDE_DIR
  CLEAR
    GMSH_EXE GMSH_LIBRARY GMSH_INCLUDE_DIR GMSH_WITH_API
  )
