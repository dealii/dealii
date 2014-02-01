## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Configuration for the lapack library:
#

MACRO(FEATURE_LAPACK_FIND_EXTERNAL var)
  FIND_PACKAGE(DEALII_LAPACK)

  IF(LAPACK_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


#
# We do a fine grained check for availability of every single LAPACK
# function we use. We have to ensure that this check is repeated every time
# the lapack library or DEAL_II_WITH_LAPACK changes.
#

SET(DEAL_II_LAPACK_FUNCTIONS
  daxpy_ dgeev_ dgeevx_ dgelsd_ dgemm_ dgemv_ dgeqrf_ dgesdd_ dgesvd_
  dgetrf_ dgetri_ dgetrs_ dorgqr_ dormqr_ dstev_ dsyevx_  dsygv_ dsygvx_
  dtrtrs_ saxpy_ sgeev_ sgeevx_ sgelsd_ sgemm_ sgemv_ sgeqrf_ sgesdd_
  sgesvd_ sgetrf_ sgetri_ sgetrs_ sorgqr_ sormqr_ sstev_ ssyevx_ ssygv_
  ssygvx_ strtrs_
  )

MACRO(CHECK_FOR_LAPACK_FUNCTIONS)
  #
  # Clear the test flags because the following test will use a C compiler
  #
  CLEAR_CMAKE_REQUIRED()
  SET(CMAKE_REQUIRED_FLAGS "${LAPACK_LINKER_FLAGS}")
  SET(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  #
  # Push -pthread as well:
  #
  ENABLE_IF_SUPPORTED(CMAKE_REQUIRED_FLAGS "-pthread")

  IF(CMAKE_C_COMPILER_WORKS)
    FOREACH(_func ${DEAL_II_LAPACK_FUNCTIONS})
      STRING(TOUPPER ${_func} _func_uppercase)
      CHECK_FUNCTION_EXISTS(${_func} HAVE_${_func_uppercase})
    ENDFOREACH()
  ELSE()
    MESSAGE(STATUS
      "No suitable C compiler was found! Skipping LAPACK symbol check."
      )
    FOREACH(_func ${DEAL_II_LAPACK_FUNCTIONS})
      SET_IF_EMPTY(HAVE_${_func_uppercase} TRUE)
    ENDFOREACH()
  ENDIF()

  RESET_CMAKE_REQUIRED()
ENDMACRO()


MACRO(RESET_LAPACK_FUNCTIONS_CACHE)
  FOREACH(_func ${DEAL_II_LAPACK_FUNCTIONS})
    STRING(TOUPPER ${_func} _func_uppercase)
    UNSET(HAVE_${_func_uppercase} CACHE)
  ENDFOREACH()
ENDMACRO()


MACRO(FEATURE_LAPACK_CONFIGURE_EXTERNAL)
  CHECK_FOR_LAPACK_FUNCTIONS()
ENDMACRO()


CONFIGURE_FEATURE(LAPACK)

#
# Call RESET_LAPACK_FUNCTIONS_CHECK if DEAL_II_WITH_LAPACK is unset to
# clean the configuration
#
IF(NOT DEAL_II_WITH_LAPACK)
  RESET_LAPACK_FUNCTIONS_CACHE()
ENDIF()
