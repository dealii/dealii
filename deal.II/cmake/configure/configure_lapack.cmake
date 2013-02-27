#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Configuration for the lapack library:
#


MACRO(FEATURE_LAPACK_FIND_EXTERNAL var)
  FIND_PACKAGE(LAPACK)

  IF(LAPACK_FOUND)
    MARK_AS_ADVANCED(
      atlas_LIBRARY
      blas_LIBRARY
      lapack_LIBRARY
      refblas_LIBRARY
      reflapack_LIBRARY
      )
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
  SET(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${LAPACK_LINKER_FLAGS}")

  FOREACH(_func ${DEAL_II_LAPACK_FUNCTIONS})
    STRING(TOUPPER ${_func} _func_uppercase)
    CHECK_FUNCTION_EXISTS(${_func} HAVE_${_func_uppercase})
  ENDFOREACH()

  SET(CMAKE_REQUIRED_LIBRARIES)
  STRIP_FLAG(CMAKE_REQUIRED_FLAGS "${LAPACK_LINKER_FLAGS}")
ENDMACRO()


MACRO(RESET_LAPACK_FUNCTIONS_CACHE)
  FOREACH(_func ${DEAL_II_LAPACK_FUNCTIONS})
    STRING(TOUPPER ${_func} _func_uppercase)
    UNSET(HAVE_${_func_uppercase} CACHE)
  ENDFOREACH()
ENDMACRO()


MACRO(FEATURE_LAPACK_CONFIGURE_EXTERNAL var)
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${LAPACK_LINKER_FLAGS}")

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    ${LAPACK_LIBRARIES}
    )

  CHECK_FOR_LAPACK_FUNCTIONS()

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(LAPACK)

#
# Call RESET_LAPACK_FUNCTIONS_CHECK if DEAL_II_WITH_LAPACK is unset to
# clean the configuration
#
IF(NOT DEAL_II_WITH_LAPACK)
  RESET_LAPACK_FUNCTIONS_CACHE()
ENDIF()
