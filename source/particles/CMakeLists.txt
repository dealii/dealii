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


INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  data_out.cc
  particle.cc
  particle_accessor.cc
  particle_iterator.cc
  particle_handler.cc
  generators.cc
  property_pool.cc
  utilities.cc
  )

SET(_inst
  data_out.inst.in
  particle.inst.in
  particle_accessor.inst.in
  particle_iterator.inst.in
  particle_handler.inst.in
  generators.inst.in
  utilities.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/particles/*.h
  )

DEAL_II_ADD_LIBRARY(obj_particle OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_particle "${_inst}")
