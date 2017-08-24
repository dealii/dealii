## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2014 by the deal.II authors
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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  matrix_free.cc
  evaluation_selector.cc
  )

SET(_inst
  matrix_free.inst.in
  evaluation_selector.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/matrix_free/*.h
  )

IF(DEAL_II_WITH_CUDA)
  SET(_src
    cuda_matrix_free.cu
    ${_src}
    )
  SET(_inst
    cuda_matrix_free.inst.in
    ${_inst}
    )
  FILE(GLOB _header
    ${CMAKE_SOURCE_DIR}/include/deal.II/matrix_free/*.cuh
    ${_header}
    )
ENDIF()

DEAL_II_ADD_LIBRARY(obj_matrix_free OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_matrix_free "${_inst}")
