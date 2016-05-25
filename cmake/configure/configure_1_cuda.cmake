## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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
# Configuration for cuda support:
#

#
# FindCUDA needs a compiler set up with C++11 support. Thus, only configure
# if deal.II was configured with C++11 support.
#
SET(FEATURE_CUDA_DEPENDS CXX11)

MACRO(FEATURE_CUDA_FIND_EXTERNAL var)

  FIND_PACKAGE(CUDA)

  IF(CUDA_FOUND)

    SET(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE FALSE)

    SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11)

    # FIXME: CUDA compiler NVCC doesn't support C++14.

    SET(${var} TRUE)
  ENDIF()

ENDMACRO()

MACRO(FEATURE_CUDA_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable cuda library!\n"
    ${CUDA_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that a cuda library is installed on your computer\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(CUDA)
