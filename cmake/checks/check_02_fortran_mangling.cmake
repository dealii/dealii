## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2021 - 2022 by the deal.II authors
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
# This file sets up
#
#   DEAL_II_FORTRAN_MANGLE
#   DEAL_II_FORTRAN_MANGLE_UNDERSCORE
#

#
# Determine the mangling procedure used to generate Fortran symbols that can be
# called from C. For example - defining a Fortran subroutine DYSEVR may result
# in a symbol named dysevr, dysevr_, DYSEVR, or DYSEVR_ depending on the
# compiler. Different strings may be prepended or appended to the symbol if it
# contains underscores.
#
# Fortunately, the rules for checking symbols are something provided by CMake.
# If we cannot inspect the Fortran compiler then default to using GNU-style
# mangling.
#

if(${DEAL_II_Fortran_COMPILER_WORKS})
  message(STATUS "Found a valid Fortran compiler - using it to set up mangling")
  include(FortranCInterface)
  if("${FortranCInterface_GLOBAL_CASE}" STREQUAL "LOWER")
    set(_name "name")
  else()
    set(_name "NAME")
  endif()

  # string(JOIN ...) is in CMake 3.12 and newer, which we don't yet require
  set(DEAL_II_FORTRAN_MANGLE "${_name}")
  if(NOT "${FortranCInterface_GLOBAL_PREFIX}" STREQUAL "")
    set(DEAL_II_FORTRAN_MANGLE
      "${FortranCInterface_GLOBAL_PREFIX} ## ${DEAL_II_FORTRAN_MANGLE}")
  endif()
  if(NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
    set(DEAL_II_FORTRAN_MANGLE
      "${DEAL_II_FORTRAN_MANGLE} ## ${FortranCInterface_GLOBAL_SUFFIX}")
  endif()

  # Same issue for identifiers containing underscores
  if("${FortranCInterface_GLOBAL__CASE}" STREQUAL "LOWER")
    set(_name "name")
  else()
    set(_name "NAME")
  endif()

  set(DEAL_II_FORTRAN_MANGLE_UNDERSCORE "${_name}")
  if(NOT "${FortranCInterface_GLOBAL__PREFIX}" STREQUAL "")
    set(DEAL_II_FORTRAN_MANGLE_UNDERSCORE
      "${FortranCInterface_GLOBAL__PREFIX} ## ${DEAL_II_FORTRAN_MANGLE_UNDERSCORE}")
  endif()
  if(NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")
    set(DEAL_II_FORTRAN_MANGLE_UNDERSCORE
      "${DEAL_II_FORTRAN_MANGLE_UNDERSCORE} ## ${FortranCInterface_GLOBAL__SUFFIX}")
  endif()
else()
  # Otherwise just use GNU mangling. This is almost always the right choice
  # anyway - until 2021 we did not even support other mangling options
  message(STATUS "Could NOT find a valid Fortran compiler - using default mangling")

  set(DEAL_II_FORTRAN_MANGLE "name ## _")
  set(DEAL_II_FORTRAN_MANGLE_UNDERSCORE "name ## _")
endif()
