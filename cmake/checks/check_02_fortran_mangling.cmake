## ---------------------------------------------------------------------
##
## Copyright (C) 2021 by the deal.II authors
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

IF(${DEAL_II_Fortran_COMPILER_WORKS})
  MESSAGE(STATUS "Found a valid Fortran compiler - using it to set up mangling")
  INCLUDE(FortranCInterface)
  IF("${FortranCInterface_GLOBAL_CASE}" STREQUAL "LOWER")
    SET(_name "name")
  ELSE()
    SET(_name "NAME")
  ENDIF()

  # STRING(JOIN ...) is in CMake 3.12 and newer, which we don't yet require
  SET(DEAL_II_FORTRAN_MANGLE "${_name}")
  IF(NOT "${FortranCInterface_GLOBAL_PREFIX}" STREQUAL "")
    SET(DEAL_II_FORTRAN_MANGLE
      "${FortranCInterface_GLOBAL_PREFIX} ## ${DEAL_II_FORTRAN_MANGLE}")
  ENDIF()
  IF(NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
    SET(DEAL_II_FORTRAN_MANGLE
      "${DEAL_II_FORTRAN_MANGLE} ## ${FortranCInterface_GLOBAL_SUFFIX}")
  ENDIF()

  # Same issue for identifiers containing underscores
  IF("${FortranCInterface_GLOBAL__CASE}" STREQUAL "LOWER")
    SET(_name "name")
  ELSE()
    SET(_name "NAME")
  ENDIF()

  SET(DEAL_II_FORTRAN_MANGLE_UNDERSCORE "${_name}")
  IF(NOT "${FortranCInterface_GLOBAL__PREFIX}" STREQUAL "")
    SET(DEAL_II_FORTRAN_MANGLE_UNDERSCORE
      "${FortranCInterface_GLOBAL__PREFIX} ## ${DEAL_II_FORTRAN_MANGLE_UNDERSCORE}")
  ENDIF()
  IF(NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")
    SET(DEAL_II_FORTRAN_MANGLE_UNDERSCORE
      "${DEAL_II_FORTRAN_MANGLE_UNDERSCORE} ## ${FortranCInterface_GLOBAL__SUFFIX}")
  ENDIF()
ELSE()
  # Otherwise just use GNU mangling. This is almost always the right choice
  # anyway - until 2021 we did not even support other mangling options
  MESSAGE(STATUS "Could NOT find a valid Fortran compiler - using default mangling")

  SET(DEAL_II_FORTRAN_MANGLE "name ## _")
  SET(DEAL_II_FORTRAN_MANGLE_UNDERSCORE "name ## _")
ENDIF()
