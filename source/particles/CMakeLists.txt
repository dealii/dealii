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
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------


include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

set(_src
  data_out.cc
  particle.cc
  particle_handler.cc
  generators.cc
  property_pool.cc
  utilities.cc
  )

set(_inst
  data_out.inst.in
  particle.inst.in
  particle_handler.inst.in
  generators.inst.in
  utilities.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/particles/*.h
  )

deal_ii_add_library(obj_particle OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(obj_particle "${_inst}")
