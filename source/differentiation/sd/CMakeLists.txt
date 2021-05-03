## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2020 by the deal.II authors
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
  symengine_math.cc
  symengine_number_types.cc
  symengine_number_visitor_internal.cc
  symengine_optimizer.cc
  symengine_scalar_operations.cc
  symengine_tensor_operations.cc
  symengine_types.cc
  symengine_utilities.cc
  )

SET(_inst
  symengine_number_visitor_internal.inst.in
  symengine_optimizer.inst.in
  symengine_tensor_operations.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/differentiation/sd/*.h
  )

DEAL_II_ADD_LIBRARY(obj_differentiation_sd OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_differentiation_sd "${_inst}")
