## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2019 by the deal.II authors
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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  ad_drivers.cc
  ad_helpers.cc
  adolc_number_types.cc
  sacado_number_types.cc
  )

SET(_inst
  ad_drivers.inst.in
  ad_drivers.inst1.in
  ad_drivers.inst2.in
  ad_helpers.inst.in
  ad_helpers.inst1.in
  ad_helpers.inst2.in
  adolc_number_types.inst.in
  sacado_number_types.inst1.in
  sacado_number_types.inst2.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/differentiation/ad/*.h
  )

DEAL_II_ADD_LIBRARY(obj_differentiation_ad OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_differentiation_ad "${_inst}")
